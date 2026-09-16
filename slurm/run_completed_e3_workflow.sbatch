#!/usr/bin/env bash
set -Eeuo pipefail

JOB_SCRATCH=""
SCRATCH_BASE=""
SCRATCH_SOURCE=""
SCRATCH_AVAILABLE_KIB=""

warn() {
    printf 'WARNING: %s\n' "$1" >&2
}

cleanup_scratch() {
    local exit_code=$?

    trap - EXIT
    if [[ -z "${JOB_SCRATCH}" ]]; then
        exit "${exit_code}"
    fi
    if (( exit_code == 0 )); then
        if [[ "${JOB_SCRATCH}" == "${SCRATCH_BASE%/}/protein_signature_${SLURM_JOB_ID}" && \
                "${JOB_SCRATCH}" != "/" ]]; then
            if rm -rf -- "${JOB_SCRATCH}"; then
                echo "Removed completed job scratch: ${JOB_SCRATCH}"
            else
                warn "Could not remove completed job scratch: ${JOB_SCRATCH}"
            fi
        else
            warn "Refusing cleanup because the scratch path failed its safety check: ${JOB_SCRATCH}"
        fi
    else
        warn "Job scratch retained after exit ${exit_code}: ${JOB_SCRATCH}"
    fi
    exit "${exit_code}"
}

try_scratch_candidate() {
    local source="$1"
    local candidate="$2"
    local required_kib="$3"
    local available_kib=""
    local proposed_job_scratch=""
    local write_probe=""

    if [[ -z "${candidate}" || "${candidate}" != /* || "${candidate}" == "/" ]]; then
        warn "Skipping invalid ${source} scratch base: ${candidate:-<empty>}"
        return 1
    fi
    if ! mkdir -p -- "${candidate}" 2>/dev/null || \
            [[ ! -d "${candidate}" || ! -w "${candidate}" || ! -x "${candidate}" ]]; then
        warn "Skipping unavailable ${source} scratch base: ${candidate}"
        return 1
    fi
    candidate="$(cd -- "${candidate}" && pwd -P)"
    if [[ "${candidate}" == "/" ]]; then
        warn "Skipping ${source} scratch base resolved to the filesystem root."
        return 1
    fi
    available_kib="$(df -Pk "${candidate}" 2>/dev/null | awk 'END {print $4}')"
    if [[ ! "${available_kib}" =~ ^[0-9]+$ ]]; then
        warn "Skipping ${source}; free space could not be read: ${candidate}"
        return 1
    fi
    if (( available_kib < required_kib )); then
        warn "Skipping ${source}; scratch has $((available_kib / 1024 / 1024)) GiB "\
"free but $((required_kib / 1024 / 1024)) GiB is required: ${candidate}"
        return 1
    fi

    proposed_job_scratch="${candidate%/}/protein_signature_${SLURM_JOB_ID}"
    if [[ -L "${proposed_job_scratch}" || \
            (-e "${proposed_job_scratch}" && ! -d "${proposed_job_scratch}") ]]; then
        warn "Skipping ${source}; per-job scratch is not an ordinary directory: "\
"${proposed_job_scratch}"
        return 1
    fi
    if ! mkdir -p -- "${proposed_job_scratch}/generic_tmp"; then
        warn "Skipping unwritable ${source} per-job scratch: ${proposed_job_scratch}"
        return 1
    fi
    write_probe="${proposed_job_scratch}/.write_probe_${SLURM_JOB_ID}"
    if ! (umask 077 && : > "${write_probe}"); then
        warn "Skipping ${source}; write probe failed: ${proposed_job_scratch}"
        return 1
    fi
    rm -f -- "${write_probe}"

    SCRATCH_SOURCE="${source}"
    SCRATCH_BASE="${candidate%/}"
    SCRATCH_AVAILABLE_KIB="${available_kib}"
    JOB_SCRATCH="${proposed_job_scratch}"
    return 0
}

select_scratch() {
    local required_gib="${PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB:-10}"
    local required_kib=""
    local user_name="${USER:-unknown_user}"
    local persistent_fallback=""

    if [[ ! "${required_gib}" =~ ^[1-9][0-9]*$ ]]; then
        echo "PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB must be a positive integer." >&2
        exit 2
    fi
    required_kib=$((required_gib * 1024 * 1024))
    if [[ -n "${PROTEIN_SIGNATURE_WORK_DIR:-}" ]]; then
        persistent_fallback="${PROTEIN_SIGNATURE_WORK_DIR%/}/workflow_state/tmp"
    fi

    if [[ -n "${PROTEIN_SIGNATURE_SCRATCH_BASE:-}" ]] && \
            try_scratch_candidate \
                "PROTEIN_SIGNATURE_SCRATCH_BASE" \
                "${PROTEIN_SIGNATURE_SCRATCH_BASE}" \
                "${required_kib}"; then
        return 0
    fi
    if [[ -n "${SLURM_TMPDIR:-}" ]] && \
            try_scratch_candidate "SLURM_TMPDIR" "${SLURM_TMPDIR}" "${required_kib}"; then
        return 0
    fi
    if [[ -n "${TMPDIR:-}" ]] && \
            try_scratch_candidate "TMPDIR" "${TMPDIR}" "${required_kib}"; then
        return 0
    fi
    if try_scratch_candidate "NODE_TMP" "/tmp/${user_name}" "${required_kib}"; then
        return 0
    fi
    if [[ -n "${persistent_fallback}" ]] && \
            try_scratch_candidate \
                "WORK_DIR_FALLBACK" "${persistent_fallback}" "${required_kib}"; then
        return 0
    fi
    echo "No writable scratch location has at least ${required_gib} GiB free." >&2
    exit 2
}

write_scratch_provenance() {
    local provenance_dir=""
    local provenance_path=""
    local staging_path=""

    if [[ -z "${PROTEIN_SIGNATURE_WORK_DIR:-}" ]]; then
        return 0
    fi
    provenance_dir="${PROTEIN_SIGNATURE_WORK_DIR%/}/slurm_logs"
    provenance_path="${provenance_dir}/protein_signature_scratch_${SLURM_JOB_ID}.tsv"
    staging_path="${provenance_path}.tmp.${SLURM_JOB_ID}"
    mkdir -p -- "${provenance_dir}"
    {
        printf 'field\tvalue\n'
        printf 'job_id\t%s\n' "${SLURM_JOB_ID}"
        printf 'source\t%s\n' "${SCRATCH_SOURCE}"
        printf 'scratch_base\t%s\n' "${SCRATCH_BASE}"
        printf 'job_scratch\t%s\n' "${JOB_SCRATCH}"
        printf 'tmpdir\t%s\n' "${TMPDIR}"
        printf 'available_gib\t%s\n' "$((SCRATCH_AVAILABLE_KIB / 1024 / 1024))"
        printf 'minimum_free_gib\t%s\n' \
            "${PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB:-10}"
    } > "${staging_path}"
    mv -- "${staging_path}" "${provenance_path}"
}

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This worker must be launched by sbatch." >&2
    exit 2
fi
if (( $# < 1 )); then
    echo "The Slurm worker requires a launcher command." >&2
    exit 2
fi

readonly REQUESTED_CPUS="${PROTEIN_SIGNATURE_REQUESTED_CPUS:-}"
if [[ ! "${REQUESTED_CPUS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "PROTEIN_SIGNATURE_REQUESTED_CPUS is missing or invalid." >&2
    exit 2
fi
if [[ ! "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]]; then
    echo "SLURM_CPUS_PER_TASK is missing or invalid." >&2
    exit 2
fi
if (( SLURM_CPUS_PER_TASK < REQUESTED_CPUS )); then
    echo "Allocated CPUs (${SLURM_CPUS_PER_TASK}) are fewer than requested CPUs " \
        "(${REQUESTED_CPUS})." >&2
    exit 2
fi
if [[ ! -x "$1" ]]; then
    echo "Submitted launcher is missing or not executable: $1" >&2
    exit 2
fi

select_scratch
export TMPDIR="${JOB_SCRATCH}/generic_tmp"
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"
write_scratch_provenance
trap cleanup_scratch EXIT
trap 'exit 130' INT TERM

echo "Starting protein-signature Slurm worker."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "Launcher: $1"
echo "Scratch source: ${SCRATCH_SOURCE}"
echo "Job scratch: ${JOB_SCRATCH}"
echo "Scratch free: $((SCRATCH_AVAILABLE_KIB / 1024 / 1024)) GiB"
unset PROTEIN_SIGNATURE_REQUESTED_CPUS
"$@"
