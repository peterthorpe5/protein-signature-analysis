#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  # Create campaign.yaml, validate it, then stop for scientific review.
  ./start_from_inputs.sh \
    --work-dir /persistent/path/campaign_001 \
    --campaign-id campaign_001 \
    --sequences-fasta /path/proteins.faa \
    --label-assignments /path/label_assignments.tsv \
    [--profile e3] [--domains /path/domains.tsv] \
    [--features /path/features.tsv] [--structures /path/structures.tsv] \
    [--structure-comparisons /path/comparisons.tsv] \
    [--structural-alignment-resource /path/completed_structural_result] \
    [--alphafold-accessions /path/alphafold_accessions.tsv --enable-alphafold] \
    [--orthofinder-resource /path/published_orthofinder_resource] \
    [--orthofinder-results /path/raw_orthofinder_results] \
    [--enable-foldseek] [--foldseek-maximum-hits 1000] \
    [--threads 8] [--initialise-only]

  # Validate and run an existing, reviewed campaign.yaml.
  ./start_from_inputs.sh \
    --work-dir /persistent/path/campaign_001 \
    --resume [--threads 8]

The script creates or synchronises a conda environment, writes and validates a
new campaign configuration, and normally runs the complete workflow. Use
--initialise-only to stop after validation so campaign.yaml can be reviewed and
edited before analysis. When campaign.yaml already exists, it is the sole
configuration authority: --resume is required and config-defining command-line
options are rejected. The script consumes OrthoFinder output but never runs
OrthoFinder itself.
EOF
}

require_option_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

WORK_DIR=""
CAMPAIGN_ID=""
PROFILE="e3"
SEQUENCES_FASTA=""
LABEL_ASSIGNMENTS=""
CONDA_ENVIRONMENT="protein_signature_analysis"
THREADS="1"
LOG_LEVEL="INFO"
RESUME="false"
INITIALISE_ONLY="false"
ENABLE_ALPHAFOLD="false"
ENABLE_FOLDSEEK="false"
FOLDSEEK_MAXIMUM_HITS="1000"
OPTIONAL_ARGUMENTS=()
CONFIG_DEFINING_OPTIONS=()
OPTIONAL_ARGUMENTS_PRESENT="false"
CONFIG_DEFINING_OPTIONS_PRESENT="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-dir)
            require_option_value "$@"
            WORK_DIR="${2:-}"
            shift 2
            ;;
        --campaign-id)
            require_option_value "$@"
            CAMPAIGN_ID="${2:-}"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift 2
            ;;
        --profile)
            require_option_value "$@"
            PROFILE="${2:-}"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift 2
            ;;
        --sequences-fasta)
            require_option_value "$@"
            SEQUENCES_FASTA="${2:-}"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift 2
            ;;
        --label-assignments)
            require_option_value "$@"
            LABEL_ASSIGNMENTS="${2:-}"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift 2
            ;;
        --conda-environment)
            require_option_value "$@"
            CONDA_ENVIRONMENT="${2:-}"
            shift 2
            ;;
        --threads)
            require_option_value "$@"
            THREADS="${2:-}"
            shift 2
            ;;
        --log-level)
            require_option_value "$@"
            LOG_LEVEL="${2:-}"
            shift 2
            ;;
        --features|--domains|--redundancy-clusters|--structures|--structure-comparisons|\
        --structural-alignment-resource|--alphafold-accessions|\
        --orthofinder-resource|--orthofinder-results|--orthofinder-group-type|\
        --orthofinder-hierarchy-node|--orthofinder-run-id)
            require_option_value "$@"
            CONFIG_DEFINING_OPTIONS+=("$1")
            OPTIONAL_ARGUMENTS+=("$1" "${2:-}")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            OPTIONAL_ARGUMENTS_PRESENT="true"
            shift 2
            ;;
        --enable-alphafold)
            ENABLE_ALPHAFOLD="true"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift
            ;;
        --enable-foldseek)
            ENABLE_FOLDSEEK="true"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift
            ;;
        --foldseek-maximum-hits)
            require_option_value "$@"
            FOLDSEEK_MAXIMUM_HITS="${2:-}"
            CONFIG_DEFINING_OPTIONS+=("$1")
            CONFIG_DEFINING_OPTIONS_PRESENT="true"
            shift 2
            ;;
        --initialise-only) INITIALISE_ONLY="true"; shift ;;
        --resume) RESUME="true"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "${WORK_DIR}" ]]; then
    echo "--work-dir is required." >&2
    usage >&2
    exit 2
fi
if [[ "${WORK_DIR}" == -* ]]; then
    echo "--work-dir must not begin with '-': ${WORK_DIR}" >&2
    exit 2
fi
if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "--threads must be a positive integer: ${THREADS}" >&2
    exit 2
fi
if [[ ! "${FOLDSEEK_MAXIMUM_HITS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "--foldseek-maximum-hits must be a positive integer: ${FOLDSEEK_MAXIMUM_HITS}" >&2
    exit 2
fi
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--conda-environment must be a valid Conda environment name." >&2
    exit 2
fi
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${WORK_DIR}" == "/" || "${WORK_DIR}" == "${HOME}" || \
        "${WORK_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing unsafe --work-dir destination: ${WORK_DIR}" >&2
    exit 2
fi
readonly CONFIG_PATH="${WORK_DIR}/campaign.yaml"
readonly RESULT_DIR="${WORK_DIR}/result"

if [[ -f "${CONFIG_PATH}" ]]; then
    if [[ "${RESUME}" != "true" ]]; then
        echo "Configuration already exists: ${CONFIG_PATH}; pass --resume to use it." >&2
        exit 2
    fi
    if [[ "${CONFIG_DEFINING_OPTIONS_PRESENT}" == "true" ]]; then
        echo "Existing ${CONFIG_PATH} is the sole configuration authority." >&2
        echo "Do not pass config-defining options with --resume: ${CONFIG_DEFINING_OPTIONS[*]}" >&2
        exit 2
    fi
    echo "Using existing campaign YAML as the sole configuration authority: ${CONFIG_PATH}"
elif [[ -z "${CAMPAIGN_ID}" || -z "${SEQUENCES_FASTA}" || \
        -z "${LABEL_ASSIGNMENTS}" ]]; then
    echo "A new campaign requires --campaign-id, --sequences-fasta and --label-assignments." >&2
    usage >&2
    exit 2
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found on PATH." >&2
    exit 2
fi

mkdir -p -- "${WORK_DIR}"

if conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
    echo "Synchronising existing conda environment from environment.yml: ${CONDA_ENVIRONMENT}"
    conda env update --file "${SCRIPT_DIR}/environment.yml" \
        --name "${CONDA_ENVIRONMENT}" --prune
else
    conda env create --file "${SCRIPT_DIR}/environment.yml" --name "${CONDA_ENVIRONMENT}"
fi
conda run --name "${CONDA_ENVIRONMENT}" \
    python -m pip install --no-deps --force-reinstall --editable "${SCRIPT_DIR}"

if [[ ! -f "${CONFIG_PATH}" ]]; then
    INITIALISE_COMMAND=(
        conda run --no-capture-output --name "${CONDA_ENVIRONMENT}"
        protein-signatures initialise
        --config "${CONFIG_PATH}"
        --campaign-id "${CAMPAIGN_ID}"
        --profile "${PROFILE}"
        --sequences-fasta "${SEQUENCES_FASTA}"
        --label-assignments "${LABEL_ASSIGNMENTS}"
        --foldseek-maximum-hits "${FOLDSEEK_MAXIMUM_HITS}"
        --log-level "${LOG_LEVEL}"
    )
    if [[ "${OPTIONAL_ARGUMENTS_PRESENT}" == "true" ]]; then
        INITIALISE_COMMAND+=("${OPTIONAL_ARGUMENTS[@]}")
    fi
    if [[ "${ENABLE_ALPHAFOLD}" == "true" ]]; then
        INITIALISE_COMMAND+=(--enable-alphafold)
    fi
    if [[ "${ENABLE_FOLDSEEK}" == "true" ]]; then
        INITIALISE_COMMAND+=(--enable-foldseek)
    fi
    "${INITIALISE_COMMAND[@]}"
fi

conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    protein-signatures validate --config "${CONFIG_PATH}" --log-level "${LOG_LEVEL}"

if [[ "${INITIALISE_ONLY}" == "true" ]]; then
    echo "Initialisation and validation complete; analysis was not started: ${CONFIG_PATH}"
    exit 0
fi

RUN_COMMAND=(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}"
    protein-signatures run-all
    --config "${CONFIG_PATH}"
    --output-dir "${RESULT_DIR}"
    --threads "${THREADS}"
    --log-level "${LOG_LEVEL}"
)
if [[ "${RESUME}" == "true" ]]; then
    RUN_COMMAND+=(--resume)
fi
"${RUN_COMMAND[@]}"
