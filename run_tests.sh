#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

bash -n \
    run_tests.sh \
    run_protein_signature_analysis.sh \
    run_protein_signature_app.sh \
    start_from_inputs.sh \
    run_completed_e3_workflow.sh \
    submit_protein_signature_workflow_slurm.sh \
    slurm/protein_signature_workflow_controller.sbatch \
    slurm/run_completed_e3_workflow.sbatch
python -m compileall -q src tests
python -m pycodestyle src tests --max-line-length=100
python -m pydocstyle --convention=google --add-ignore=D105,D107,D202 \
    src/protein_signatures src/protein_signature_app
python -m ruff format --check src tests
python -m ruff check src tests
python -m coverage erase
python -m coverage run --branch -m pytest -q
python -m coverage report --fail-under=95
