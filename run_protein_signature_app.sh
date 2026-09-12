#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./run_protein_signature_app.sh \
    --resource /absolute/path/completed_result \
    [--conda-environment protein_signature_analysis] \
    [--port 8501] \
    [--address localhost]
EOF
}

require_option_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

RESOURCE=""
CONDA_ENVIRONMENT="protein_signature_analysis"
PORT="8501"
ADDRESS="localhost"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --resource)
            require_option_value "$@"
            RESOURCE="${2:-}"
            shift 2
            ;;
        --conda-environment)
            require_option_value "$@"
            CONDA_ENVIRONMENT="${2:-}"
            shift 2
            ;;
        --port)
            require_option_value "$@"
            PORT="${2:-}"
            shift 2
            ;;
        --address)
            require_option_value "$@"
            ADDRESS="${2:-}"
            shift 2
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

if [[ -z "${RESOURCE}" ]]; then
    echo "--resource is required." >&2
    usage >&2
    exit 2
fi
if [[ ! -d "${RESOURCE}" ]]; then
    echo "Resource directory does not exist: ${RESOURCE}" >&2
    exit 2
fi
if [[ ! "${PORT}" =~ ^[1-9][0-9]*$ ]] || (( PORT > 65535 )); then
    echo "--port must be an integer from 1 to 65535: ${PORT}" >&2
    exit 2
fi
if [[ -z "${ADDRESS}" ]] || [[ "${ADDRESS}" =~ [[:space:]] ]]; then
    echo "--address must be non-empty and contain no whitespace: ${ADDRESS}" >&2
    exit 2
fi
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "--conda-environment must be a valid Conda environment name." >&2
    exit 2
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found on PATH." >&2
    exit 2
fi

conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    protein-signature-app --resource "${RESOURCE}" --port "${PORT}" --address "${ADDRESS}"
