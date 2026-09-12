#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

python -m compileall -q src tests
python -m pycodestyle src tests --max-line-length=100
python -m pydocstyle --convention=google --add-ignore=D105,D107,D202 \
    src/protein_signatures src/protein_signature_app
python -m ruff format --check src tests
python -m ruff check src tests
python -m coverage erase
python -m coverage run --branch -m pytest -q
python -m coverage report --fail-under=95
