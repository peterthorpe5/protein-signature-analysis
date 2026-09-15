#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./submit_protein_signature_workflow_slurm.sh \
    --config /absolute/path/campaign.yaml \
    --work-dir /persistent/path/signature_campaign \
    [--threads 24] [--memory-mb 64000] [--runtime-minutes 1440] \
    [--max-jobs 10] [--account barton] [--partition barton] \
    [--controller-memory 4G] [--controller-time 3-00:00:00] \
    [--conda-environment protein_signature_analysis] [--resume] [--dry-run]

This login-node launcher submits a small durable controller. The controller runs
the generic Snakemake workflow with its Slurm executor; scientific rules receive
the requested threads, memory and runtime. OrthoFinder inputs are consumed, never
generated. Scheduler and workflow logs remain under WORK_DIR/workflow_state.
EOF
}

require_option_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" || "${2:-}" == --* ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

validate_slurm_time() {
    local value="$1"
    if [[ "${value}" =~ ^([0-9]+)-([0-9]{2}):([0-9]{2}):([0-9]{2})$ ]]; then
        (( 10#${BASH_REMATCH[2]} <= 23 )) || return 1
        (( 10#${BASH_REMATCH[3]} <= 59 )) || return 1
        (( 10#${BASH_REMATCH[4]} <= 59 )) || return 1
        return 0
    fi
    if [[ "${value}" =~ ^([0-9]+):([0-9]{2}):([0-9]{2})$ ]]; then
        (( 10#${BASH_REMATCH[2]} <= 59 )) || return 1
        (( 10#${BASH_REMATCH[3]} <= 59 )) || return 1
        return 0
    fi
    return 1
}

CONFIG=""
WORK_DIR=""
OUTPUT_DIR=""
WORKFLOW_STATE_DIR=""
THREADS="24"
MEMORY_MB="64000"
RUNTIME_MINUTES="1440"
MAX_JOBS="10"
ACCOUNT="barton"
PARTITION="barton"
CONTROLLER_MEMORY="4G"
CONTROLLER_TIME="3-00:00:00"
CONDA_ENVIRONMENT="protein_signature_analysis"
LOG_LEVEL="INFO"
RESUME="false"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            require_option_value "$@"
            CONFIG="${2:-}"
            shift 2
            ;;
        --work-dir)
            require_option_value "$@"
            WORK_DIR="${2:-}"
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
        --controller-memory)
            require_option_value "$@"
            CONTROLLER_MEMORY="${2:-}"
            shift 2
            ;;
        --controller-time)
            require_option_value "$@"
            CONTROLLER_TIME="${2:-}"
            shift 2
            ;;
        --conda-environment)
            require_option_value "$@"
            CONDA_ENVIRONMENT="${2:-}"
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

if [[ -z "${CONFIG}" || -z "${WORK_DIR}" ]]; then
    echo "--config and --work-dir are required." >&2
    usage >&2
    exit 2
fi
if [[ "${CONFIG}" != /* || ! -f "${CONFIG}" ]]; then
    echo "--config must be an existing absolute file: ${CONFIG}" >&2
    exit 2
fi
if [[ "${WORK_DIR}" != /* || "${WORK_DIR}" == "/" || \
        "${WORK_DIR}" == "${HOME:-}" ]]; then
    echo "--work-dir must be a safe absolute persistent directory." >&2
    exit 2
fi
for value in "${THREADS}" "${MEMORY_MB}" "${RUNTIME_MINUTES}" "${MAX_JOBS}"; do
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Threads, memory, runtime and max-jobs must be positive integers." >&2
        exit 2
    fi
done
if [[ ! "${ACCOUNT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ || \
        ! "${PARTITION}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--account and --partition contain unsafe characters." >&2
    exit 2
fi
if [[ ! "${CONTROLLER_MEMORY}" =~ ^[1-9][0-9]*[KMGT]?$ ]]; then
    echo "--controller-memory must be a positive Slurm size such as 4G." >&2
    exit 2
fi
if ! validate_slurm_time "${CONTROLLER_TIME}"; then
    echo "--controller-time must use valid HH:MM:SS or D-HH:MM:SS syntax." >&2
    exit 2
fi
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--conda-environment must be a valid Conda environment name." >&2
    exit 2
fi
if [[ ! "${LOG_LEVEL}" =~ ^(DEBUG|INFO|WARNING|ERROR|CRITICAL)$ ]]; then
    echo "--log-level must be DEBUG, INFO, WARNING, ERROR or CRITICAL." >&2
    exit 2
fi
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "Submit the workflow controller from a login node, not job ${SLURM_JOB_ID}." >&2
    exit 2
fi

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly CONTROLLER_WORKER="${SCRIPT_DIR}/slurm/protein_signature_workflow_controller.sbatch"
readonly RUNNER="${SCRIPT_DIR}/run_protein_signature_analysis.sh"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_DIR}/result}"
WORKFLOW_STATE_DIR="${WORKFLOW_STATE_DIR:-${WORK_DIR}/workflow_state}"
if [[ "${OUTPUT_DIR}" != /* || "${WORKFLOW_STATE_DIR}" != /* ]]; then
    echo "Output and workflow-state directories must be absolute." >&2
    exit 2
fi
if [[ "${WORK_DIR}" == "${SCRIPT_DIR}" || "${OUTPUT_DIR}" == "/" || \
        "${OUTPUT_DIR}" == "${HOME:-}" || "${WORKFLOW_STATE_DIR}" == "/" || \
        "${WORKFLOW_STATE_DIR}" == "${HOME:-}" || \
        "${WORKFLOW_STATE_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing an unsafe work, output or workflow-state destination." >&2
    exit 2
fi
if [[ ! -x "${CONTROLLER_WORKER}" || ! -x "${RUNNER}" ]]; then
    echo "Workflow runner or Slurm controller worker is missing or not executable." >&2
    exit 2
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found on PATH." >&2
    exit 2
fi
if [[ "${DRY_RUN}" != "true" ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is unavailable; submit from a Slurm login node." >&2
    exit 2
fi

RUNNER_ARGUMENTS=(
    --config "${CONFIG}"
    --output-dir "${OUTPUT_DIR}"
    --workflow-state-dir "${WORKFLOW_STATE_DIR}"
    --profile slurm
    --conda-environment "${CONDA_ENVIRONMENT}"
    --threads "${THREADS}"
    --memory-mb "${MEMORY_MB}"
    --runtime-minutes "${RUNTIME_MINUTES}"
    --max-jobs "${MAX_JOBS}"
    --account "${ACCOUNT}"
    --partition "${PARTITION}"
    --log-level "${LOG_LEVEL}"
    --allow-inside-slurm
    --skip-environment-sync
)
if [[ "${RESUME}" == "true" ]]; then
    RUNNER_ARGUMENTS+=(--resume)
fi
SBATCH_ARGUMENTS=(
    --parsable
    --job-name=protein_signature_controller
    "--account=${ACCOUNT}"
    "--partition=${PARTITION}"
    "--mem=${CONTROLLER_MEMORY}"
    "--time=${CONTROLLER_TIME}"
    --cpus-per-task=1
    "--chdir=${SCRIPT_DIR}"
    "--output=${WORKFLOW_STATE_DIR}/slurm/controller_%j.out"
    "--error=${WORKFLOW_STATE_DIR}/slurm/controller_%j.err"
)

if [[ "${DRY_RUN}" == "true" ]]; then
    printf 'Validated controller command; nothing was submitted:\n  env -u SLURM_CPUS_PER_TASK sbatch'
    printf ' %q' "${SBATCH_ARGUMENTS[@]}" "${CONTROLLER_WORKER}" \
        --state-dir "${WORKFLOW_STATE_DIR}" -- "${RUNNER}" "${RUNNER_ARGUMENTS[@]}"
    printf '\n'
    exit 0
fi

if conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
    conda env update --file "${SCRIPT_DIR}/environment.yml" \
        --name "${CONDA_ENVIRONMENT}" --prune
else
    conda env create --file "${SCRIPT_DIR}/environment.yml" \
        --name "${CONDA_ENVIRONMENT}"
fi
conda run --name "${CONDA_ENVIRONMENT}" \
    python -m pip install --no-deps --force-reinstall --editable "${SCRIPT_DIR}"
conda run --name "${CONDA_ENVIRONMENT}" snakemake --version >/dev/null

mkdir -p -- "${WORKFLOW_STATE_DIR}/slurm"
SUBMISSION_RESULT="$(
    env -u SLURM_CPUS_PER_TASK sbatch \
        "${SBATCH_ARGUMENTS[@]}" \
        "${CONTROLLER_WORKER}" \
        --state-dir "${WORKFLOW_STATE_DIR}" \
        -- "${RUNNER}" "${RUNNER_ARGUMENTS[@]}"
)"
JOB_ID="${SUBMISSION_RESULT%%;*}"
if [[ ! "${JOB_ID}" =~ ^[0-9]+$ ]]; then
    echo "sbatch returned an unexpected job identifier: ${SUBMISSION_RESULT}" >&2
    exit 2
fi
echo "Submitted protein-signature Snakemake controller job ${JOB_ID}."
echo "Controller stdout: ${WORKFLOW_STATE_DIR}/slurm/controller_${JOB_ID}.out"
echo "Controller stderr: ${WORKFLOW_STATE_DIR}/slurm/controller_${JOB_ID}.err"
echo "Rule logs: ${WORKFLOW_STATE_DIR}/logs"
echo "Monitor: squeue --job ${JOB_ID}"
