#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./run_protein_signature_analysis.sh \
    --config /absolute/path/campaign.yaml \
    --output-dir /absolute/path/new_result \
    [--conda-environment protein_signature_analysis] \
    [--threads 8] [--resume] [--log-level INFO]

This launcher creates the declared conda environment when necessary, installs
the local package and analyses supplied inputs. It imports published resources
or raw layouts tested with OrthoFinder 2.5.5 and 3 when configured; it never
runs OrthoFinder. Raw layouts must contain OrthoFinder's exact completion marker.
EOF
}

require_option_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

CONFIG=""
OUTPUT_DIR=""
CONDA_ENVIRONMENT="protein_signature_analysis"
THREADS="1"
LOG_LEVEL="INFO"
RESUME="false"

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

if [[ -z "${CONFIG}" || -z "${OUTPUT_DIR}" ]]; then
    echo "--config and --output-dir are required." >&2
    usage >&2
    exit 2
fi
if [[ "${OUTPUT_DIR}" == -* ]]; then
    echo "--output-dir must not begin with '-': ${OUTPUT_DIR}" >&2
    exit 2
fi
if [[ ! -f "${CONFIG}" ]]; then
    echo "Configuration file does not exist: ${CONFIG}" >&2
    exit 2
fi
if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "--threads must be a positive integer: ${THREADS}" >&2
    exit 2
fi
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--conda-environment must be a valid Conda environment name." >&2
    exit 2
fi
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${OUTPUT_DIR}" == "/" || "${OUTPUT_DIR}" == "${HOME}" || \
        "${OUTPUT_DIR}" == "${SCRIPT_DIR}" ]]; then
    echo "Refusing unsafe --output-dir destination: ${OUTPUT_DIR}" >&2
    exit 2
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found on PATH." >&2
    exit 2
fi
if ! conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
    conda env create --file "${SCRIPT_DIR}/environment.yml" --name "${CONDA_ENVIRONMENT}"
fi

conda run --name "${CONDA_ENVIRONMENT}" \
    python -m pip install --no-deps --force-reinstall --editable "${SCRIPT_DIR}"

COMMAND=(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}"
    protein-signatures run-all
    --config "${CONFIG}"
    --output-dir "${OUTPUT_DIR}"
    --threads "${THREADS}"
    --log-level "${LOG_LEVEL}"
)
if [[ "${RESUME}" == "true" ]]; then
    COMMAND+=(--resume)
fi
"${COMMAND[@]}"
