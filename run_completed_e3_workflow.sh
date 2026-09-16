#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  # Phase 1: Snakemake verifies/prepares inputs and safely stages the review file.
  ./run_completed_e3_workflow.sh \
    --phase prepare \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/signature_campaign \
    --campaign-id e3_all1972_signatures_20260914 \
    [--minimum-mean-plddt 50]

  # After editing reviewed_label_assignments.tsv, issue checksum-bound approval.
  ./run_completed_e3_workflow.sh \
    --phase approve \
    --work-dir /persistent/path/signature_campaign \
    --curator "Curator name" \
    [--review-note "Reviewed against named evidence authorities"]

  # Start-to-finish Snakemake execution after approval, submitted to Slurm.
  ./run_completed_e3_workflow.sh \
    --phase all \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/signature_campaign \
    --campaign-id e3_all1972_signatures_20260914 \
    --submit-slurm \
    --slurm-account barton \
    --slurm-partition barton \
    --slurm-memory 128G \
    --slurm-time 2-00:00:00 \
    --threads 24

  # Evidence-led PROVISIONAL run with automated labels and matched controls.
  ./run_completed_e3_workflow.sh \
    --phase all \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/e3_evidence_signatures \
    --campaign-id e3_evidence_signatures \
    --evidence-led-labels \
    --accept-provisional-evidence-labels \
    [--evidence-rules e3] \
    [--seed-assignments /absolute/path/reviewed_seed_labels.tsv] \
    [--seed-catalogue /absolute/path/e3_seed_catalogue.tsv] \
    [--external-annotations /absolute/path/external_annotations.tsv] \
    --submit-slurm \
    --slurm-account barton \
    --slurm-partition barton \
    --slurm-memory 128G \
    --slurm-time 2-00:00:00 \
    --threads 24

  # Fully automated SOFTWARE SMOKE TEST with synthetic labels and no human review.
  # Use a separate work directory and campaign ID containing "smoke" or "test".
  ./run_completed_e3_workflow.sh \
    --phase all \
    --run-root /absolute/path/completed_e3_end_to_end_run \
    --work-dir /persistent/path/e3_automated_smoke_test \
    --campaign-id e3_automated_smoke_test \
    --automated-test-labels \
    --test-target-label ALL \
    --test-samples-per-class 20 \
    --submit-slurm \
    --slurm-account barton \
    --slurm-partition barton \
    --slurm-memory 128G \
    --slurm-time 2-00:00:00 \
    --threads 24

  # Optional checkpoints and backward-compatible operations.
  ./run_completed_e3_workflow.sh \
    --phase initialise|run|verify \
    --work-dir /persistent/path/signature_campaign [...]

The predecessor is never modified and OrthoFinder is never run. The E3
Snakemake DAG owns preparation, review staging, approval verification,
initialisation, atomic analysis and independent final verification. It never
overwrites reviewed_label_assignments.tsv. Generated assignments begin as
UNMAPPED, so the one intentional pause is human curation plus explicit approval.
After approval, --phase all is unattended and safely resumable.

--automated-test-labels is a separate end-to-end software-test route. It creates
deterministic synthetic target/control cohorts, automated test-only approval and
conspicuous provenance. Its results MUST NOT be interpreted scientifically.
  --test-target-label ID  Profile target/alias, or ALL (default: ALL).
  --test-samples-per-class N  Synthetic samples per direct cohort (default: 20;
                             minimum: 20).

--evidence-led-labels uses versioned annotation/Pfam rules, optional curated
direct labels, one-generation unanimous HOG propagation, structure eligibility
and prespecified outcome-blind matching. Without
--accept-provisional-evidence-labels it stops after the evidence bundle is ready.
With that explicit flag it runs unattended, but results remain provisional
hypothesis-generation evidence until human scientific review.

Slurm options:
  --submit-slurm          Submit prepare, initialise, all or run.
  --slurm-account NAME   Account (default: barton).
  --slurm-partition NAME Partition (default: barton).
  --slurm-memory SIZE    Memory request (default: 128G).
  --slurm-time TIME      Wall time as HH:MM:SS or D-HH:MM:SS (default: 2-00:00:00).
  --slurm-job-name NAME  Job name (default: protein_signature_PHASE).
  --slurm-log-dir PATH   Absolute log directory (default: WORK_DIR/slurm_logs).
  --slurm-scratch-base PATH
                         Optional preferred scratch base. The worker otherwise
                         tries SLURM_TMPDIR, TMPDIR, node /tmp and finally a
                         persistent WORK_DIR fallback.
  --slurm-min-scratch-free-gib N
                         Minimum free space accepted for scratch (default: 10 GiB).
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
REVIEW_APPROVAL=""
CURATOR=""
REVIEW_NOTE=""
PROFILE="e3"
AUTOMATED_TEST_LABELS="false"
AUTOMATED_TEST_OPTION_SEEN="false"
EVIDENCE_LED_LABELS="false"
ACCEPT_PROVISIONAL_EVIDENCE="false"
EVIDENCE_OPTION_SEEN="false"
EVIDENCE_RULES="e3"
SEED_ASSIGNMENTS=""
SEED_CATALOGUE=""
EXTERNAL_ANNOTATIONS=""
TEST_TARGET_LABEL="ALL"
TEST_SAMPLES_PER_CLASS="20"
MINIMUM_MEAN_PLDDT="50"
CONDA_ENVIRONMENT="protein_signature_analysis"
THREADS="1"
LOG_LEVEL="INFO"
RESUME="false"
SUBMIT_SLURM="false"
SLURM_ACCOUNT="barton"
SLURM_PARTITION="barton"
SLURM_MEMORY="128G"
SLURM_TIME="2-00:00:00"
SLURM_JOB_NAME=""
SLURM_LOG_DIR=""
SLURM_SCRATCH_BASE="${PROTEIN_SIGNATURE_SCRATCH_BASE:-}"
SLURM_MIN_SCRATCH_FREE_GIB="${PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB:-10}"
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
        --review-approval)
            require_option_value "$@"
            REVIEW_APPROVAL="${2:-}"
            shift 2
            ;;
        --curator)
            require_option_value "$@"
            CURATOR="${2:-}"
            shift 2
            ;;
        --review-note)
            require_option_value "$@"
            REVIEW_NOTE="${2:-}"
            shift 2
            ;;
        --profile)
            require_option_value "$@"
            PROFILE="${2:-}"
            shift 2
            ;;
        --automated-test-labels)
            AUTOMATED_TEST_LABELS="true"
            shift
            ;;
        --evidence-led-labels)
            EVIDENCE_LED_LABELS="true"
            shift
            ;;
        --accept-provisional-evidence-labels)
            ACCEPT_PROVISIONAL_EVIDENCE="true"
            EVIDENCE_OPTION_SEEN="true"
            shift
            ;;
        --evidence-rules)
            require_option_value "$@"
            EVIDENCE_RULES="${2:-}"
            EVIDENCE_OPTION_SEEN="true"
            shift 2
            ;;
        --seed-assignments)
            require_option_value "$@"
            SEED_ASSIGNMENTS="${2:-}"
            EVIDENCE_OPTION_SEEN="true"
            shift 2
            ;;
        --seed-catalogue)
            require_option_value "$@"
            SEED_CATALOGUE="${2:-}"
            EVIDENCE_OPTION_SEEN="true"
            shift 2
            ;;
        --external-annotations)
            require_option_value "$@"
            EXTERNAL_ANNOTATIONS="${2:-}"
            EVIDENCE_OPTION_SEEN="true"
            shift 2
            ;;
        --test-target-label)
            require_option_value "$@"
            TEST_TARGET_LABEL="${2:-}"
            AUTOMATED_TEST_OPTION_SEEN="true"
            shift 2
            ;;
        --test-samples-per-class)
            require_option_value "$@"
            TEST_SAMPLES_PER_CLASS="${2:-}"
            AUTOMATED_TEST_OPTION_SEEN="true"
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
        --slurm-scratch-base)
            require_option_value "$@"
            SLURM_SCRATCH_BASE="${2:-}"
            SLURM_OPTION_SEEN="true"
            shift 2
            ;;
        --slurm-min-scratch-free-gib)
            require_option_value "$@"
            SLURM_MIN_SCRATCH_FREE_GIB="${2:-}"
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
if [[ "${PHASE}" != "prepare" && "${PHASE}" != "approve" && \
        "${PHASE}" != "initialise" && "${PHASE}" != "all" && \
        "${PHASE}" != "run" && "${PHASE}" != "verify" ]]; then
    echo "--phase must be prepare, approve, initialise, all, run or verify: ${PHASE}" >&2
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
if [[ ! "${TEST_SAMPLES_PER_CLASS}" =~ ^[1-9][0-9]*$ ]] || \
        (( TEST_SAMPLES_PER_CLASS < 20 )); then
    echo "--test-samples-per-class must be an integer of at least 20." >&2
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
if [[ "${AUTOMATED_TEST_OPTION_SEEN}" == "true" && \
        "${AUTOMATED_TEST_LABELS}" != "true" ]]; then
    echo "Automated test label options require --automated-test-labels." >&2
    exit 2
fi
if [[ "${AUTOMATED_TEST_LABELS}" == "true" && \
        "${EVIDENCE_LED_LABELS}" == "true" ]]; then
    echo "--automated-test-labels and --evidence-led-labels are mutually exclusive." >&2
    exit 2
fi
if [[ "${EVIDENCE_OPTION_SEEN}" == "true" && \
        "${EVIDENCE_LED_LABELS}" != "true" ]]; then
    echo "Evidence-label options require --evidence-led-labels." >&2
    exit 2
fi
if [[ "${EVIDENCE_LED_LABELS}" == "true" && "${PHASE}" != "all" ]]; then
    echo "--evidence-led-labels is supported only with --phase all." >&2
    exit 2
fi
if [[ -z "${EVIDENCE_RULES}" ]]; then
    echo "--evidence-rules must not be empty." >&2
    exit 2
fi
if [[ "${EVIDENCE_RULES}" == */* ]]; then
    if [[ "${EVIDENCE_RULES}" != /* || ! -s "${EVIDENCE_RULES}" ]]; then
        echo "Custom --evidence-rules must be an existing absolute file." >&2
        exit 2
    fi
elif [[ ! "${EVIDENCE_RULES}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "Built-in --evidence-rules contains unsafe characters." >&2
    exit 2
fi
for evidence_file in "${SEED_ASSIGNMENTS}" "${SEED_CATALOGUE}" \
        "${EXTERNAL_ANNOTATIONS}"; do
    if [[ -n "${evidence_file}" && \
            ("${evidence_file}" != /* || ! -s "${evidence_file}") ]]; then
        echo "Evidence authorities must be existing non-empty absolute files: ${evidence_file}" >&2
        exit 2
    fi
done

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${WORK_DIR}" == "/" || "${WORK_DIR}" == "${HOME:-}" || \
        "${WORK_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing unsafe --work-dir destination: ${WORK_DIR}" >&2
    exit 2
fi
readonly PREPARED_DIR="${WORK_DIR}/prepared_inputs"
readonly CONFIG_PATH="${WORK_DIR}/campaign.yaml"
readonly RESULT_DIR="${WORK_DIR}/result"
readonly E3_STATE_DIR="${WORK_DIR}/workflow_state/e3"
readonly PREPARATION_MARKER="${E3_STATE_DIR}/01_preparation/PREPARED_VERIFIED.json"
readonly REVIEW_MARKER="${E3_STATE_DIR}/02_label_review/REVIEW_READY.json"
CAMPAIGN_ID="${CAMPAIGN_ID:-$(basename "${WORK_DIR}")}"
readonly AUTOMATED_LABELS_PATH="${WORK_DIR}/AUTOMATED_TEST_ONLY.label_assignments.tsv"
readonly AUTOMATED_LABELS_MARKER="${E3_STATE_DIR}/02_label_review/AUTOMATED_TEST_ONLY.LABELS.json"
readonly AUTOMATED_APPROVAL_PATH="${E3_STATE_DIR}/02_label_review/AUTOMATED_TEST_ONLY.REVIEW_APPROVED.json"
readonly EVIDENCE_BUNDLE_DIR="${WORK_DIR}/evidence_label_bundle"
readonly EVIDENCE_LABELS_PATH="${EVIDENCE_BUNDLE_DIR}/label_assignments.tsv"
readonly EVIDENCE_APPROVAL_PATH="${E3_STATE_DIR}/02_label_review/EVIDENCE_PROVISIONAL.REVIEW_APPROVED.json"
if [[ "${AUTOMATED_TEST_LABELS}" == "true" ]]; then
    if [[ "${PHASE}" != "all" ]]; then
        echo "--automated-test-labels is supported only with --phase all." >&2
        exit 2
    fi
    if [[ ! "${CAMPAIGN_ID}" =~ ([Ss][Mm][Oo][Kk][Ee]|[Tt][Ee][Ss][Tt]) ]]; then
        echo "Automated labels require 'smoke' or 'test' in --campaign-id." >&2
        exit 2
    fi
    if [[ ! "$(basename "${WORK_DIR}")" =~ ([Ss][Mm][Oo][Kk][Ee]|[Tt][Ee][Ss][Tt]) ]]; then
        echo "Automated labels require a separate work directory containing 'smoke' or 'test'." >&2
        exit 2
    fi
    if [[ -n "${LABEL_ASSIGNMENTS}" && "${LABEL_ASSIGNMENTS}" != "${AUTOMATED_LABELS_PATH}" ]]; then
        echo "Automated mode reserves its fixed AUTOMATED_TEST_ONLY label path." >&2
        exit 2
    fi
    if [[ -n "${REVIEW_APPROVAL}" && "${REVIEW_APPROVAL}" != "${AUTOMATED_APPROVAL_PATH}" ]]; then
        echo "Automated mode reserves its fixed AUTOMATED_TEST_ONLY approval path." >&2
        exit 2
    fi
    LABEL_ASSIGNMENTS="${AUTOMATED_LABELS_PATH}"
    REVIEW_APPROVAL="${AUTOMATED_APPROVAL_PATH}"
elif [[ "${EVIDENCE_LED_LABELS}" == "true" ]]; then
    if [[ -n "${LABEL_ASSIGNMENTS}" && \
            "${LABEL_ASSIGNMENTS}" != "${EVIDENCE_LABELS_PATH}" ]]; then
        echo "Evidence-led mode reserves its fixed evidence-bundle label path." >&2
        exit 2
    fi
    if [[ -n "${REVIEW_APPROVAL}" && \
            "${REVIEW_APPROVAL}" != "${EVIDENCE_APPROVAL_PATH}" ]]; then
        echo "Evidence-led mode reserves its fixed provisional approval path." >&2
        exit 2
    fi
    LABEL_ASSIGNMENTS="${EVIDENCE_LABELS_PATH}"
    REVIEW_APPROVAL="${EVIDENCE_APPROVAL_PATH}"
else
    LABEL_ASSIGNMENTS="${LABEL_ASSIGNMENTS:-${WORK_DIR}/reviewed_label_assignments.tsv}"
    REVIEW_APPROVAL="${REVIEW_APPROVAL:-${E3_STATE_DIR}/02_label_review/REVIEW_APPROVED.json}"
fi
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
        --campaign-id "${CAMPAIGN_ID}"
        --label-assignments "${LABEL_ASSIGNMENTS}"
        --review-approval "${REVIEW_APPROVAL}"
        --profile "${PROFILE}"
    )
    local -a sbatch_arguments=(
        --parsable
        "--job-name=${job_name}"
        "--mem=${SLURM_MEMORY}"
        "--time=${SLURM_TIME}"
        "--cpus-per-task=${THREADS}"
        "--export=ALL"
        "--chdir=${SCRIPT_DIR}"
        "--output=${log_dir}/%x_%j.out"
        "--error=${log_dir}/%x_%j.err"
    )

    if [[ "${PHASE}" == "prepare" || "${PHASE}" == "initialise" || \
            "${PHASE}" == "all" ]]; then
        worker_arguments+=(
            --run-root "${RUN_ROOT}"
            --minimum-mean-plddt "${MINIMUM_MEAN_PLDDT}"
        )
    elif [[ "${RESUME}" == "true" ]]; then
        worker_arguments+=(--resume)
    fi
    if [[ "${AUTOMATED_TEST_LABELS}" == "true" ]]; then
        worker_arguments+=(
            --automated-test-labels
            --test-target-label "${TEST_TARGET_LABEL}"
            --test-samples-per-class "${TEST_SAMPLES_PER_CLASS}"
        )
    fi
    if [[ "${EVIDENCE_LED_LABELS}" == "true" ]]; then
        worker_arguments+=(
            --evidence-led-labels
            --evidence-rules "${EVIDENCE_RULES}"
        )
        if [[ "${ACCEPT_PROVISIONAL_EVIDENCE}" == "true" ]]; then
            worker_arguments+=(--accept-provisional-evidence-labels)
        fi
        if [[ -n "${SEED_ASSIGNMENTS}" ]]; then
            worker_arguments+=(--seed-assignments "${SEED_ASSIGNMENTS}")
        fi
        if [[ -n "${SEED_CATALOGUE}" ]]; then
            worker_arguments+=(--seed-catalogue "${SEED_CATALOGUE}")
        fi
        if [[ -n "${EXTERNAL_ANNOTATIONS}" ]]; then
            worker_arguments+=(--external-annotations "${EXTERNAL_ANNOTATIONS}")
        fi
    fi
    if [[ -n "${SLURM_ACCOUNT}" ]]; then
        sbatch_arguments+=("--account=${SLURM_ACCOUNT}")
    fi
    if [[ -n "${SLURM_PARTITION}" ]]; then
        sbatch_arguments+=("--partition=${SLURM_PARTITION}")
    fi

    if [[ "${SLURM_DRY_RUN}" == "true" ]]; then
        printf 'Validated Slurm command; nothing was submitted:\n  env -u SLURM_CPUS_PER_TASK'
        printf ' PROTEIN_SIGNATURE_REQUESTED_CPUS=%q' "${THREADS}"
        printf ' PROTEIN_SIGNATURE_WORK_DIR=%q' "${WORK_DIR}"
        printf ' PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB=%q' \
            "${SLURM_MIN_SCRATCH_FREE_GIB}"
        printf ' PROTEIN_SIGNATURE_SCRATCH_BASE=%q sbatch' "${SLURM_SCRATCH_BASE}"
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
        env -u SLURM_CPUS_PER_TASK \
            PROTEIN_SIGNATURE_REQUESTED_CPUS="${THREADS}" \
            PROTEIN_SIGNATURE_WORK_DIR="${WORK_DIR}" \
            PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB="${SLURM_MIN_SCRATCH_FREE_GIB}" \
            PROTEIN_SIGNATURE_SCRATCH_BASE="${SLURM_SCRATCH_BASE}" \
            sbatch \
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
    if [[ "${PHASE}" != "prepare" && "${PHASE}" != "initialise" && \
            "${PHASE}" != "all" && "${PHASE}" != "run" ]]; then
        echo "--submit-slurm supports prepare, initialise, all and run phases." >&2
        exit 2
    fi
    if [[ "${WORK_DIR}" != /* ]]; then
        echo "Slurm submission requires an absolute --work-dir." >&2
        exit 2
    fi
    if [[ "${PHASE}" != "run" && ("${RUN_ROOT}" != /* || ! -d "${RUN_ROOT}") ]]; then
        echo "Slurm ${PHASE} requires an existing absolute --run-root: ${RUN_ROOT}" >&2
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
    if [[ -n "${SLURM_SCRATCH_BASE}" && "${SLURM_SCRATCH_BASE}" != /* ]]; then
        echo "--slurm-scratch-base must be absolute." >&2
        exit 2
    fi
    if [[ ! "${SLURM_MIN_SCRATCH_FREE_GIB}" =~ ^[1-9][0-9]*$ ]]; then
        echo "--slurm-min-scratch-free-gib must be a positive integer." >&2
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

run_e3_snakemake() {
    local target="$1"
    local -a command=(
        conda run --no-capture-output --name "${CONDA_ENVIRONMENT}"
        snakemake
        --snakefile "${SCRIPT_DIR}/workflow/E3Snakefile"
        --directory "${E3_STATE_DIR}"
        --profile "${SCRIPT_DIR}/profiles/local"
        --cores "${THREADS}"
        --jobs 1
        --rerun-incomplete
        --printshellcmds
        "${target}"
        --config
        "e3_run_root=${RUN_ROOT}"
        "e3_work_dir=${WORK_DIR}"
        "e3_campaign_id=${CAMPAIGN_ID}"
        "e3_profile=${PROFILE}"
        "e3_minimum_mean_plddt=${MINIMUM_MEAN_PLDDT}"
        "e3_reviewed_labels=${LABEL_ASSIGNMENTS}"
        "e3_review_approval=${REVIEW_APPROVAL}"
        "e3_automated_test_labels=${AUTOMATED_TEST_LABELS}"
        "e3_evidence_led_labels=${EVIDENCE_LED_LABELS}"
        "e3_accept_provisional_evidence_labels=${ACCEPT_PROVISIONAL_EVIDENCE}"
        "e3_evidence_rules=${EVIDENCE_RULES}"
        "e3_seed_assignments=${SEED_ASSIGNMENTS}"
        "e3_seed_catalogue=${SEED_CATALOGUE}"
        "e3_external_annotations=${EXTERNAL_ANNOTATIONS}"
        "e3_test_target_label=${TEST_TARGET_LABEL}"
        "e3_test_samples_per_class=${TEST_SAMPLES_PER_CLASS}"
        "e3_threads=${THREADS}"
        "e3_memory_mb=128000"
        "e3_runtime_minutes=2880"
        "e3_log_level=${LOG_LEVEL}"
    )

    if [[ -z "${RUN_ROOT}" || ! -d "${RUN_ROOT}" ]]; then
        echo "${PHASE} requires an existing --run-root directory: ${RUN_ROOT}" >&2
        exit 2
    fi
    mkdir -p -- "${WORK_DIR}" "${E3_STATE_DIR}"
    ensure_environment
    if ! conda run --name "${CONDA_ENVIRONMENT}" snakemake --version >/dev/null 2>&1; then
        echo "Snakemake is unavailable in Conda environment ${CONDA_ENVIRONMENT}." >&2
        exit 2
    fi
    "${command[@]}"
}

if [[ "${PHASE}" == "approve" ]]; then
    if [[ -z "${CURATOR}" ]]; then
        echo "approve requires --curator with the responsible person's name." >&2
        exit 2
    fi
    if [[ ! -s "${PREPARATION_MARKER}" || ! -s "${REVIEW_MARKER}" ]]; then
        echo "Review staging markers are missing; complete --phase prepare first." >&2
        exit 2
    fi
    if [[ ! -s "${LABEL_ASSIGNMENTS}" ]]; then
        echo "Reviewed label authority is missing or empty: ${LABEL_ASSIGNMENTS}" >&2
        exit 2
    fi
    ensure_environment
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
        protein-signatures approve-e3-review \
        --preparation-marker "${PREPARATION_MARKER}" \
        --review-marker "${REVIEW_MARKER}" \
        --reviewed-labels "${LABEL_ASSIGNMENTS}" \
        --approval-marker "${REVIEW_APPROVAL}" \
        --curator "${CURATOR}" \
        --note "${REVIEW_NOTE}" \
        --profile "${PROFILE}" \
        --log-level "${LOG_LEVEL}"
    echo "Approved current reviewed-label checksum: ${REVIEW_APPROVAL}"
    echo "The campaign can now continue with --phase all."
    exit 0
fi

if [[ "${PHASE}" == "prepare" ]]; then
    run_e3_snakemake review_ready
    echo "Snakemake prepared and staged the editable review authority:"
    echo "  ${LABEL_ASSIGNMENTS}"
    echo "Curate that file, then run --phase approve with --curator."
    exit 0
fi

if [[ "${PHASE}" == "initialise" ]]; then
    if [[ ! -s "${REVIEW_APPROVAL}" ]]; then
        echo "Checksum-bound label approval is missing: ${REVIEW_APPROVAL}" >&2
        echo "Complete --phase approve before initialisation." >&2
        exit 2
    fi
    run_e3_snakemake campaign_ready
    echo "Snakemake created or adopted and validated ${CONFIG_PATH}."
    echo "The expensive analysis has not started."
    exit 0
fi

if [[ "${PHASE}" == "all" ]]; then
    if [[ "${AUTOMATED_TEST_LABELS}" == "true" ]]; then
        run_e3_snakemake all
        echo "Completed AUTOMATED TEST ONLY campaign: ${RESULT_DIR}"
        echo "These synthetic-label results test software execution and MUST NOT be " \
            "interpreted scientifically."
        exit 0
    fi
    if [[ "${EVIDENCE_LED_LABELS}" == "true" ]]; then
        if [[ "${ACCEPT_PROVISIONAL_EVIDENCE}" != "true" ]]; then
            run_e3_snakemake review_ready
            echo "Evidence-led label and matched-control bundle is ready:"
            echo "  ${EVIDENCE_BUNDLE_DIR}"
            echo "No analysis was started. Inspect its TSV/XLSX audits, then rerun the same"
            echo "command with --accept-provisional-evidence-labels to continue unattended."
            exit 0
        fi
        run_e3_snakemake all
        echo "Completed provisional evidence-led E3 signature result: ${RESULT_DIR}"
        echo "Interpret as hypothesis-generation evidence until human scientific review."
        exit 0
    fi
    if [[ ! -s "${REVIEW_APPROVAL}" ]]; then
        run_e3_snakemake review_ready
        echo "Workflow paused cleanly at the required human-curation checkpoint."
        echo "Edit ${LABEL_ASSIGNMENTS}, then run --phase approve with --curator."
        echo "After approval, resubmit this same --phase all command."
        exit 0
    fi
    if [[ ! -s "${LABEL_ASSIGNMENTS}" ]]; then
        echo "Approved label authority is missing or empty: ${LABEL_ASSIGNMENTS}" >&2
        exit 2
    fi
    run_e3_snakemake all
    echo "Completed and independently verified E3 signature result: ${RESULT_DIR}"
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
