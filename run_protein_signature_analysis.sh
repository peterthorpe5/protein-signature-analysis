#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./run_protein_signature_analysis.sh \
    --config /absolute/path/campaign.yaml \
    --output-dir /absolute/path/new_result \
    [--workflow-state-dir /persistent/path/workflow_state] \
    [--profile local|slurm|/absolute/custom/profile] \
    [--conda-environment protein_signature_analysis] \
    [--threads 8] [--memory-mb 64000] [--runtime-minutes 240] \
    [--max-jobs 10] [--account barton] [--partition general] \
    [--resume] [--dry-run] [--unlock] [--log-level INFO]

This is the generic input-driven Snakemake launcher. It consumes campaign.yaml
created from any supported FASTA, reviewed labels, optional domains, structures,
AlphaFold requests and OrthoFinder 2.5.5/3 results. It never runs OrthoFinder and
does not depend on the completed-E3 workflow adapter.

The local profile executes the complete transactional analysis in the current
process or an enclosing batch allocation. The slurm profile uses Snakemake's
Slurm executor and must remain under a durable controller; use
submit_protein_signature_workflow_slurm.sh for logout-safe execution.
EOF
}

require_option_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" || "${2:-}" == --* ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

CONFIG=""
OUTPUT_DIR=""
WORKFLOW_STATE_DIR=""
PROFILE="local"
CONDA_ENVIRONMENT="protein_signature_analysis"
THREADS="1"
MEMORY_MB="64000"
RUNTIME_MINUTES="240"
MAX_JOBS="10"
ACCOUNT="barton"
PARTITION="general"
LOG_LEVEL="INFO"
RESUME="false"
DRY_RUN="false"
UNLOCK="false"
ALLOW_INSIDE_SLURM="false"
SKIP_ENVIRONMENT_SYNC="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            require_option_value "$@"
            CONFIG="${2:-}"
            shift 2
            ;;
        --output-dir)
            require_option_value "$@"
            OUTPUT_DIR="${2:-}"
            shift 2
            ;;
        --workflow-state-dir)
            require_option_value "$@"
            WORKFLOW_STATE_DIR="${2:-}"
            shift 2
            ;;
        --profile)
            require_option_value "$@"
            PROFILE="${2:-}"
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
        --memory-mb)
            require_option_value "$@"
            MEMORY_MB="${2:-}"
            shift 2
            ;;
        --runtime-minutes)
            require_option_value "$@"
            RUNTIME_MINUTES="${2:-}"
            shift 2
            ;;
        --max-jobs)
            require_option_value "$@"
            MAX_JOBS="${2:-}"
            shift 2
            ;;
        --account)
            require_option_value "$@"
            ACCOUNT="${2:-}"
            shift 2
            ;;
        --partition)
            require_option_value "$@"
            PARTITION="${2:-}"
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
        --dry-run)
            DRY_RUN="true"
            shift
            ;;
        --unlock)
            UNLOCK="true"
            shift
            ;;
        --allow-inside-slurm)
            ALLOW_INSIDE_SLURM="true"
            shift
            ;;
        --skip-environment-sync)
            SKIP_ENVIRONMENT_SYNC="true"
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

if [[ -z "${CONFIG}" || -z "${OUTPUT_DIR}" ]]; then
    echo "--config and --output-dir are required." >&2
    usage >&2
    exit 2
fi
if [[ ! -f "${CONFIG}" ]]; then
    echo "Configuration file does not exist: ${CONFIG}" >&2
    exit 2
fi
CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd -P)/$(basename "${CONFIG}")"
for option_value in "${THREADS}" "${MEMORY_MB}" "${RUNTIME_MINUTES}" "${MAX_JOBS}"; do
    if [[ ! "${option_value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Threads, memory, runtime and max-jobs must be positive integers." >&2
        exit 2
    fi
done
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--conda-environment must be a valid Conda environment name." >&2
    exit 2
fi
if [[ ! "${LOG_LEVEL}" =~ ^(DEBUG|INFO|WARNING|ERROR|CRITICAL)$ ]]; then
    echo "--log-level must be DEBUG, INFO, WARNING, ERROR or CRITICAL." >&2
    exit 2
fi
if [[ ! "${ACCOUNT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ || \
        ! "${PARTITION}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--account and --partition must contain only safe scheduler characters." >&2
    exit 2
fi

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${OUTPUT_DIR}" == -* ]]; then
    echo "--output-dir must not begin with '-': ${OUTPUT_DIR}" >&2
    exit 2
fi
if [[ "${OUTPUT_DIR}" == "/" || \
        "${OUTPUT_DIR}" == "${HOME:-}" || "${OUTPUT_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing unsafe --output-dir destination: ${OUTPUT_DIR}" >&2
    exit 2
fi
if [[ -z "${WORKFLOW_STATE_DIR}" ]]; then
    WORKFLOW_STATE_DIR="$(dirname "${OUTPUT_DIR}")/workflow_state"
fi
if [[ "${WORKFLOW_STATE_DIR}" == -* || "${WORKFLOW_STATE_DIR}" == "/" || \
        "${WORKFLOW_STATE_DIR}" == "${HOME:-}" || \
        "${WORKFLOW_STATE_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing unsafe --workflow-state-dir destination: ${WORKFLOW_STATE_DIR}" >&2
    exit 2
fi
mkdir -p -- "$(dirname "${OUTPUT_DIR}")" "${WORKFLOW_STATE_DIR}"
OUTPUT_DIR="$(cd "$(dirname "${OUTPUT_DIR}")" && pwd -P)/$(basename "${OUTPUT_DIR}")"
WORKFLOW_STATE_DIR="$(cd "${WORKFLOW_STATE_DIR}" && pwd -P)"
if [[ "${PROFILE}" == "local" || "${PROFILE}" == "slurm" ]]; then
    PROFILE="${SCRIPT_DIR}/profiles/${PROFILE}"
elif [[ "${PROFILE}" != /* ]]; then
    echo "A custom --profile must be an absolute directory." >&2
    exit 2
fi
if [[ ! -d "${PROFILE}" ]]; then
    echo "Snakemake profile does not exist: ${PROFILE}" >&2
    exit 2
fi
if [[ "${PROFILE}" == "${SCRIPT_DIR}/profiles/slurm" && \
        -n "${SLURM_JOB_ID:-}" && "${ALLOW_INSIDE_SLURM}" != "true" ]]; then
    echo "Refusing to start the Slurm executor inside job ${SLURM_JOB_ID}." >&2
    echo "Use submit_protein_signature_workflow_slurm.sh from a login node." >&2
    exit 2
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found on PATH." >&2
    exit 2
fi

if [[ "${SKIP_ENVIRONMENT_SYNC}" != "true" ]]; then
    if conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
        conda env update --file "${SCRIPT_DIR}/environment.yml" \
            --name "${CONDA_ENVIRONMENT}" --prune
    else
        conda env create --file "${SCRIPT_DIR}/environment.yml" \
            --name "${CONDA_ENVIRONMENT}"
    fi
    conda run --name "${CONDA_ENVIRONMENT}" \
        python -m pip install --no-deps --force-reinstall --editable "${SCRIPT_DIR}"
fi
if ! conda run --name "${CONDA_ENVIRONMENT}" snakemake --version >/dev/null 2>&1; then
    echo "Snakemake is unavailable in Conda environment ${CONDA_ENVIRONMENT}." >&2
    exit 2
fi

SNAKEMAKE_COMMAND=(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}"
    snakemake
    --snakefile "${SCRIPT_DIR}/workflow/Snakefile"
    --configfile "${CONFIG}"
    --directory "${WORKFLOW_STATE_DIR}"
    --profile "${PROFILE}"
    --cores "${THREADS}"
    --jobs "${MAX_JOBS}"
    --rerun-triggers mtime
    --config
    "workflow_state_dir=${WORKFLOW_STATE_DIR}"
    "workflow_result_dir=${OUTPUT_DIR}"
    "workflow_threads=${THREADS}"
    "workflow_memory_mb=${MEMORY_MB}"
    "workflow_runtime_minutes=${RUNTIME_MINUTES}"
    "workflow_resume=${RESUME}"
    "workflow_log_level=${LOG_LEVEL}"
)
if [[ "${PROFILE}" == "${SCRIPT_DIR}/profiles/slurm" ]]; then
    SNAKEMAKE_COMMAND+=(
        --default-resources
        "slurm_account=${ACCOUNT}"
        "slurm_partition=${PARTITION}"
        "mem_mb=8000"
        "runtime=60"
    )
fi
if [[ "${UNLOCK}" == "true" ]]; then
    SNAKEMAKE_COMMAND+=(--unlock)
elif [[ "${DRY_RUN}" == "true" ]]; then
    SNAKEMAKE_COMMAND+=(
        --dry-run
        --printshellcmds
        --forcerun
        validate_campaign
        run_campaign
        verify_campaign
    )
else
    SNAKEMAKE_COMMAND+=(
        --forcerun
        validate_campaign
        run_campaign
        verify_campaign
    )
fi
"${SNAKEMAKE_COMMAND[@]}"
