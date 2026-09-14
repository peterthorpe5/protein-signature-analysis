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

  # Submit the memory-intensive prepare phase from a Slurm login node.
  ./run_completed_e3_workflow.sh \
    --phase prepare \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/signature_campaign \
    --submit-slurm \
    --slurm-account ACCOUNT \
    --slurm-partition PARTITION \
    --slurm-memory 64G \
    --slurm-time 04:00:00 \
    --threads 4

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

Slurm options:
  --submit-slurm          Submit prepare or run instead of executing locally.
  --slurm-account NAME   Account; omitted by default so the site default applies.
  --slurm-partition NAME Partition; omitted by default so the site default applies.
  --slurm-memory SIZE    Memory request (default: 64G).
  --slurm-time TIME      Wall time as HH:MM:SS or D-HH:MM:SS (default: 04:00:00).
  --slurm-job-name NAME  Job name (default: protein_signature_PHASE).
  --slurm-log-dir PATH   Absolute log directory (default: WORK_DIR/slurm_logs).
  --slurm-dry-run        Print the exact validated sbatch command without submitting.
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
SUBMIT_SLURM="false"
SLURM_ACCOUNT=""
SLURM_PARTITION=""
SLURM_MEMORY="64G"
SLURM_TIME="04:00:00"
SLURM_JOB_NAME=""
SLURM_LOG_DIR=""
SLURM_DRY_RUN="false"
SLURM_OPTION_SEEN="false"

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
        --submit-slurm)
            SUBMIT_SLURM="true"
            shift
            ;;
        --slurm-account)
            require_option_value "$@"
            SLURM_ACCOUNT="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-partition)
            require_option_value "$@"
            SLURM_PARTITION="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-memory)
            require_option_value "$@"
            SLURM_MEMORY="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-time|--slurm-walltime)
            require_option_value "$@"
            SLURM_TIME="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-job-name)
            require_option_value "$@"
            SLURM_JOB_NAME="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-log-dir)
            require_option_value "$@"
            SLURM_LOG_DIR="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-dry-run)
            SLURM_DRY_RUN="true"
            SLURM_OPTION_SEEN="true"
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
if [[ "${SLURM_OPTION_SEEN}" == "true" && "${SUBMIT_SLURM}" != "true" ]]; then
    echo "Slurm options require --submit-slurm." >&2
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
readonly SLURM_WORKER="${SCRIPT_DIR}/slurm/run_completed_e3_workflow.sbatch"

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

submit_slurm_phase() {
    local job_name="${SLURM_JOB_NAME:-protein_signature_${PHASE}}"
    local log_dir="${SLURM_LOG_DIR:-${WORK_DIR}/slurm_logs}"
    local submission_result=""
    local job_id=""
    local -a worker_arguments=(
        --phase "${PHASE}"
        --work-dir "${WORK_DIR}"
        --conda-environment "${CONDA_ENVIRONMENT}"
        --threads "${THREADS}"
        --log-level "${LOG_LEVEL}"
    )
    local -a sbatch_arguments=(
        --parsable
        "--job-name=${job_name}"
        "--mem=${SLURM_MEMORY}"
        "--time=${SLURM_TIME}"
        "--cpus-per-task=${THREADS}"
        "--export=ALL,PROTEIN_SIGNATURE_REQUESTED_CPUS=${THREADS}"
        "--chdir=${SCRIPT_DIR}"
        "--output=${log_dir}/%x_%j.out"
        "--error=${log_dir}/%x_%j.err"
    )

    if [[ "${PHASE}" == "prepare" ]]; then
        worker_arguments+=(
            --run-root "${RUN_ROOT}"
            --minimum-mean-plddt "${MINIMUM_MEAN_PLDDT}"
        )
    elif [[ "${RESUME}" == "true" ]]; then
        worker_arguments+=(--resume)
    fi
    if [[ -n "${SLURM_ACCOUNT}" ]]; then
        sbatch_arguments+=("--account=${SLURM_ACCOUNT}")
    fi
    if [[ -n "${SLURM_PARTITION}" ]]; then
        sbatch_arguments+=("--partition=${SLURM_PARTITION}")
    fi

    if [[ "${SLURM_DRY_RUN}" == "true" ]]; then
        printf 'Validated Slurm command; nothing was submitted:\n  env -u SLURM_CPUS_PER_TASK sbatch'
        printf ' %q' "${sbatch_arguments[@]}" "${SLURM_WORKER}" \
            "${SCRIPT_DIR}/run_completed_e3_workflow.sh" "${worker_arguments[@]}"
        printf '\n'
        return 0
    fi

    command -v sbatch >/dev/null 2>&1 || {
        echo "sbatch is unavailable; submit from a Slurm login node." >&2
        exit 2
    }
    if [[ ! -x "${SLURM_WORKER}" ]]; then
        echo "Slurm worker is missing or not executable: ${SLURM_WORKER}" >&2
        exit 2
    fi
    mkdir -p -- "${log_dir}"
    submission_result="$(
        env -u SLURM_CPUS_PER_TASK sbatch \
            "${sbatch_arguments[@]}" \
            "${SLURM_WORKER}" \
            "${SCRIPT_DIR}/run_completed_e3_workflow.sh" \
            "${worker_arguments[@]}"
    )"
    job_id="${submission_result%%;*}"
    if [[ ! "${job_id}" =~ ^[0-9]+$ ]]; then
        echo "sbatch returned an unexpected job identifier: ${submission_result}" >&2
        exit 2
    fi
    echo "Submitted Slurm job ${job_id} for phase ${PHASE}."
    echo "Slurm stdout: ${log_dir}/${job_name}_${job_id}.out"
    echo "Slurm stderr: ${log_dir}/${job_name}_${job_id}.err"
    echo "Monitor: squeue --job ${job_id}"
}

if [[ "${SUBMIT_SLURM}" == "true" ]]; then
    if [[ "${PHASE}" != "prepare" && "${PHASE}" != "run" ]]; then
        echo "--submit-slurm currently supports prepare and run phases only." >&2
        exit 2
    fi
    if [[ "${WORK_DIR}" != /* ]]; then
        echo "Slurm submission requires an absolute --work-dir." >&2
        exit 2
    fi
    if [[ "${PHASE}" == "prepare" && ("${RUN_ROOT}" != /* || ! -d "${RUN_ROOT}") ]]; then
        echo "Slurm prepare requires an existing absolute --run-root: ${RUN_ROOT}" >&2
        exit 2
    fi
    if [[ "${PHASE}" == "run" && ! -s "${CONFIG_PATH}" ]]; then
        echo "Campaign configuration is missing; complete initialise first: ${CONFIG_PATH}" >&2
        exit 2
    fi
    if [[ ! "${SLURM_MEMORY}" =~ ^[1-9][0-9]*[KMGT]?$ ]]; then
        echo "--slurm-memory must be a positive Slurm size such as 64G." >&2
        exit 2
    fi
    if ! validate_slurm_time "${SLURM_TIME}"; then
        echo "--slurm-time must use valid HH:MM:SS or D-HH:MM:SS syntax." >&2
        exit 2
    fi
    if [[ -n "${SLURM_ACCOUNT}" && \
            ! "${SLURM_ACCOUNT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "--slurm-account contains unsafe characters." >&2
        exit 2
    fi
    if [[ -n "${SLURM_PARTITION}" && \
            ! "${SLURM_PARTITION}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "--slurm-partition contains unsafe characters." >&2
        exit 2
    fi
    if [[ -n "${SLURM_JOB_NAME}" && \
            ! "${SLURM_JOB_NAME}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ ]]; then
        echo "--slurm-job-name must be a safe name of at most 64 characters." >&2
        exit 2
    fi
    if [[ -n "${SLURM_LOG_DIR}" && "${SLURM_LOG_DIR}" != /* ]]; then
        echo "--slurm-log-dir must be absolute." >&2
        exit 2
    fi
    submit_slurm_phase
    exit 0
fi

ensure_environment() {
    if ! command -v conda >/dev/null 2>&1; then
        echo "conda is required but was not found on PATH." >&2
        exit 2
    fi
    if conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
        conda env update \
            --file "${SCRIPT_DIR}/environment.yml" \
            --name "${CONDA_ENVIRONMENT}" \
            --prune
    else
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
