#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  # Phase 1: verify the completed predecessor and prepare review-gated inputs.
  ./run_completed_e3_workflow.sh \
    --phase prepare \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/signature_campaign \
    [--minimum-mean-plddt 50]

  # Phase 2: after copying and curating the generated label template.
  ./run_completed_e3_workflow.sh \
    --phase initialise \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/signature_campaign \
    --campaign-id e3_all1972_signatures_20260914 \
    --label-assignments /persistent/path/reviewed_label_assignments.tsv \
    [--threads 24]

  # Phase 3: run the reviewed campaign, normally inside a scheduler allocation.
  ./run_completed_e3_workflow.sh \
    --phase run \
    --work-dir /persistent/path/signature_campaign \
    [--threads 24] [--resume]

  # Phase 4: independently verify every published result checksum.
  ./run_completed_e3_workflow.sh \
    --phase verify \
    --work-dir /persistent/path/signature_campaign

The predecessor is never modified and OrthoFinder is never run. Phase prepare
creates proteins.faa, an explicit Pfam assessment ledger, a checksum-verified
structure inventory and a label curation worksheet. Its generated assignments
are all UNMAPPED. Phase initialise therefore requires a separate reviewed label
authority, imports Stage 09b US-align/TM-align and pocket evidence, reuses Stage
09 AlphaFold coordinate assets, enables campaign-wide Foldseek and consumes the
completed raw OrthoFinder 2.5.5/3 results from Stage 04.
EOF
}

require_option_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

PHASE=""
RUN_ROOT=""
WORK_DIR=""
CAMPAIGN_ID=""
LABEL_ASSIGNMENTS=""
MINIMUM_MEAN_PLDDT="50"
CONDA_ENVIRONMENT="protein_signature_analysis"
THREADS="1"
LOG_LEVEL="INFO"
RESUME="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phase)
            require_option_value "$@"
            PHASE="${2:-}"
            shift 2
            ;;
        --run-root)
            require_option_value "$@"
            RUN_ROOT="${2:-}"
            shift 2
            ;;
        --work-dir)
            require_option_value "$@"
            WORK_DIR="${2:-}"
            shift 2
            ;;
        --campaign-id)
            require_option_value "$@"
            CAMPAIGN_ID="${2:-}"
            shift 2
            ;;
        --label-assignments)
            require_option_value "$@"
            LABEL_ASSIGNMENTS="${2:-}"
            shift 2
            ;;
        --minimum-mean-plddt)
            require_option_value "$@"
            MINIMUM_MEAN_PLDDT="${2:-}"
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
        --resume)
            RESUME="true"
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "${PHASE}" || -z "${WORK_DIR}" ]]; then
    echo "--phase and --work-dir are required." >&2
    usage >&2
    exit 2
fi
if [[ "${PHASE}" != "prepare" && "${PHASE}" != "initialise" && \
        "${PHASE}" != "run" && "${PHASE}" != "verify" ]]; then
    echo "--phase must be prepare, initialise, run or verify: ${PHASE}" >&2
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
if [[ ! "${MINIMUM_MEAN_PLDDT}" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]] || \
        ! awk -v value="${MINIMUM_MEAN_PLDDT}" \
            'BEGIN { exit !(value >= 0 && value <= 100) }'; then
    echo "--minimum-mean-plddt must be a number from 0 to 100." >&2
    exit 2
fi
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--conda-environment must be a valid Conda environment name." >&2
    exit 2
fi

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${WORK_DIR}" == "/" || "${WORK_DIR}" == "${HOME:-}" || \
        "${WORK_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing unsafe --work-dir destination: ${WORK_DIR}" >&2
    exit 2
fi
readonly PREPARED_DIR="${WORK_DIR}/prepared_inputs"
readonly CONFIG_PATH="${WORK_DIR}/campaign.yaml"
readonly RESULT_DIR="${WORK_DIR}/result"

ensure_environment() {
    if ! command -v conda >/dev/null 2>&1; then
        echo "conda is required but was not found on PATH." >&2
        exit 2
    fi
    if ! conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
        conda env create \
            --file "${SCRIPT_DIR}/environment.yml" \
            --name "${CONDA_ENVIRONMENT}"
    fi
    conda run --name "${CONDA_ENVIRONMENT}" \
        python -m pip install --no-deps --force-reinstall --editable "${SCRIPT_DIR}"
}

if [[ "${PHASE}" == "prepare" ]]; then
    if [[ -z "${RUN_ROOT}" || ! -d "${RUN_ROOT}" ]]; then
        echo "prepare requires an existing --run-root directory: ${RUN_ROOT}" >&2
        exit 2
    fi
    mkdir -p -- "${WORK_DIR}"
    ensure_environment
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
        protein-signatures prepare-e3-workflow \
        --run-root "${RUN_ROOT}" \
        --output-dir "${PREPARED_DIR}" \
        --minimum-mean-plddt "${MINIMUM_MEAN_PLDDT}" \
        --log-level "${LOG_LEVEL}"
    echo "Prepared inputs require label and control curation before initialisation:"
    echo "  ${PREPARED_DIR}/e3_label_curation_review.tsv"
    echo "  ${PREPARED_DIR}/label_assignments.REVIEW_REQUIRED.tsv"
    exit 0
fi

if [[ "${PHASE}" == "initialise" ]]; then
    if [[ -z "${RUN_ROOT}" || ! -d "${RUN_ROOT}" ]]; then
        echo "initialise requires an existing --run-root directory: ${RUN_ROOT}" >&2
        exit 2
    fi
    if [[ -z "${CAMPAIGN_ID}" || -z "${LABEL_ASSIGNMENTS}" ]]; then
        echo "initialise requires --campaign-id and --label-assignments." >&2
        exit 2
    fi
    if [[ ! -s "${PREPARED_DIR}/PREPARED.json" ]]; then
        echo "Prepared input marker is missing; complete the prepare phase first." >&2
        exit 2
    fi
    if [[ ! -s "${LABEL_ASSIGNMENTS}" ]]; then
        echo "Reviewed label authority is missing or empty: ${LABEL_ASSIGNMENTS}" >&2
        exit 2
    fi
    if [[ "$(basename "${LABEL_ASSIGNMENTS}")" == \
            "label_assignments.REVIEW_REQUIRED.tsv" ]]; then
        echo "Refusing the generated UNMAPPED template; supply a reviewed copy." >&2
        exit 2
    fi
    readonly STRUCTURE_COUNT="$(
        awk -F '\t' '
            NR == 1 {
                for (column = 1; column <= NF; column++) {
                    if ($column == "analysis_eligibility_status") {
                        eligibility_column = column
                    }
                }
                if (!eligibility_column) {
                    exit 2
                }
                next
            }
            $eligibility_column == "ELIGIBLE" { count++ }
            END {
                if (!eligibility_column) {
                    exit 2
                }
                print count + 0
            }
        ' "${PREPARED_DIR}/structures.tsv"
    )"
    if [[ ! "${STRUCTURE_COUNT}" =~ ^[0-9]+$ ]] || (( STRUCTURE_COUNT < 2 )); then
        echo "At least two Foldseek-eligible prepared structures are required: ${STRUCTURE_COUNT}" >&2
        exit 2
    fi
    readonly ORTHOFINDER_RESULTS="${RUN_ROOT}/04_orthofinder/Results"
    "${SCRIPT_DIR}/start_from_inputs.sh" \
        --work-dir "${WORK_DIR}" \
        --campaign-id "${CAMPAIGN_ID}" \
        --profile e3 \
        --sequences-fasta "${PREPARED_DIR}/proteins.faa" \
        --label-assignments "${LABEL_ASSIGNMENTS}" \
        --domains "${PREPARED_DIR}/domains.tsv" \
        --structures "${PREPARED_DIR}/structures.tsv" \
        --structural-alignment-resource "${RUN_ROOT}" \
        --orthofinder-results "${ORTHOFINDER_RESULTS}" \
        --orthofinder-group-type HOG \
        --orthofinder-hierarchy-node N0 \
        --enable-foldseek \
        --foldseek-maximum-hits "${STRUCTURE_COUNT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --threads "${THREADS}" \
        --log-level "${LOG_LEVEL}" \
        --initialise-only
    echo "Review and freeze ${CONFIG_PATH}; the expensive analysis has not started."
    exit 0
fi

if [[ ! -s "${CONFIG_PATH}" ]]; then
    echo "Campaign configuration is missing; complete initialise first: ${CONFIG_PATH}" >&2
    exit 2
fi

if [[ "${PHASE}" == "run" ]]; then
    RUN_ARGUMENTS=(
        --config "${CONFIG_PATH}"
        --output-dir "${RESULT_DIR}"
        --conda-environment "${CONDA_ENVIRONMENT}"
        --threads "${THREADS}"
        --log-level "${LOG_LEVEL}"
    )
    if [[ "${RESUME}" == "true" ]]; then
        RUN_ARGUMENTS+=(--resume)
    fi
    "${SCRIPT_DIR}/run_protein_signature_analysis.sh" "${RUN_ARGUMENTS[@]}"
    exit 0
fi

ensure_environment
conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    protein-signatures verify --resource "${RESULT_DIR}" --log-level "${LOG_LEVEL}"
