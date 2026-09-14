#!/usr/bin/env bash
set -Eeuo pipefail

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
if [[ "${SLURM_CPUS_PER_TASK}" != "${REQUESTED_CPUS}" ]]; then
    echo "Allocated CPUs (${SLURM_CPUS_PER_TASK}) differ from requested CPUs " \
        "(${REQUESTED_CPUS})." >&2
    exit 2
fi
if [[ ! -x "$1" ]]; then
    echo "Submitted launcher is missing or not executable: $1" >&2
    exit 2
fi

echo "Starting protein-signature Slurm worker."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "Launcher: $1"
unset PROTEIN_SIGNATURE_REQUESTED_CPUS
exec "$@"
