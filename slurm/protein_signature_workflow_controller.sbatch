#!/usr/bin/env bash
set -Eeuo pipefail

STATE_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --state-dir)
            if (( $# < 2 )) || [[ -z "${2:-}" ]]; then
                echo "--state-dir requires a value." >&2
                exit 2
            fi
            STATE_DIR="${2:-}"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unknown controller-worker option: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This controller worker must be launched by sbatch." >&2
    exit 2
fi
if [[ -z "${STATE_DIR}" || "${STATE_DIR}" != /* || "${STATE_DIR}" == "/" ]]; then
    echo "The controller requires a safe absolute workflow state directory." >&2
    exit 2
fi
if (( $# < 1 )) || [[ ! -x "$1" ]]; then
    echo "The controller requires an executable workflow launcher." >&2
    exit 2
fi
if ! command -v flock >/dev/null 2>&1; then
    echo "flock is required for exclusive controller ownership." >&2
    exit 2
fi

mkdir -p -- "${STATE_DIR}/control" "${STATE_DIR}/logs"
exec 9>"${STATE_DIR}/control/controller.lock"
if ! flock --nonblock 9; then
    echo "Another controller already owns workflow state ${STATE_DIR}." >&2
    exit 75
fi

echo "Starting protein-signature Snakemake controller."
echo "Controller job ID: ${SLURM_JOB_ID}"
echo "Controller host: $(hostname)"
echo "Workflow state: ${STATE_DIR}"
printf 'Controller command:'
printf ' %q' "$@"
printf '\n'
exec "$@"
