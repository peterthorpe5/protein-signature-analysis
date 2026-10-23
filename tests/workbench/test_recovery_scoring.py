"""Interrupted-worker recovery and real frozen-signature integration checks."""

import os
import socket
import subprocess
from dataclasses import replace
from types import SimpleNamespace

import pandas as pd
import pytest

from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic
from protein_signatures.workbench.analysis import terminate_process
from protein_signatures.workbench.jobs import (
    artifact_path,
    job_directory,
    read_job,
    recover_interrupted_job,
    run_next,
    submit_job,
)


def test_stopped_worker_recovery_retains_files_and_never_retries(workbench_config):
    """A running record blocks further claims until an explicit safe recovery action."""
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    path = job_directory(config=workbench_config, job_id=job)
    state = read_job(config=workbench_config, job_id=job)
    state.update(status="RUNNING", worker_pid=2_000_000_000, worker_host=socket.gethostname())
    write_json_atomic(path=path / "state.json", value=state)
    assert run_next(config=workbench_config) is None
    recover_interrupted_job(config=workbench_config, job_id=job)
    assert read_job(config=workbench_config, job_id=job)["status"] == "FAILED"
    assert (path / "request.json").is_file()
    assert run_next(config=workbench_config) is None
    with pytest.raises(InputValidationError):
        recover_interrupted_job(config=workbench_config, job_id=job)


def test_recovery_refuses_live_worker_or_tool_group(workbench_config, monkeypatch):
    """Neither GUI nor CLI can mark a job failed while its recorded processes are live."""
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    path = job_directory(config=workbench_config, job_id=job)
    state = read_job(config=workbench_config, job_id=job)
    state.update(status="RUNNING", worker_pid=os.getpid(), worker_host=socket.gethostname())
    write_json_atomic(path=path / "state.json", value=state)
    with pytest.raises(InputValidationError, match="still live"):
        recover_interrupted_job(config=workbench_config, job_id=job)
    state["worker_pid"] = 2_000_000_000
    write_json_atomic(path=path / "state.json", value=state)
    (path / "outputs").mkdir()
    write_json_atomic(
        path=path / "outputs" / "active_process.json", value={"process_group": os.getpgrp()}
    )
    monkeypatch.setattr(os, "killpg", lambda *_args: None)
    with pytest.raises(InputValidationError, match="tool group is still live"):
        recover_interrupted_job(config=workbench_config, job_id=job)


def test_process_group_shutdown_escalation_and_exit_race(monkeypatch):
    """An unresponsive tool gets a final kill; an already exited child is still reaped."""
    observed = []
    waits = iter([subprocess.TimeoutExpired(cmd="tool", timeout=5), None])

    def fake_wait(*, timeout=None):
        """Simulate one grace-period timeout and a successful final reap."""
        value = next(waits)
        if isinstance(value, Exception):
            raise value
        return value

    monkeypatch.setattr(os, "killpg", lambda *args: observed.append(args))
    terminate_process(process=SimpleNamespace(pid=1234, wait=fake_wait))
    assert len(observed) == 2

    def exited_group(*_args):
        """Simulate the child exiting immediately before termination."""
        raise ProcessLookupError

    monkeypatch.setattr(os, "killpg", exited_group)
    reaped = []
    terminate_process(process=SimpleNamespace(pid=1234, wait=lambda: reaped.append(True)))
    assert reaped == [True]


def test_real_reference_signature_scoring_keeps_query_output_separate(
    completed_result, workbench_config, seed_bytes
):
    """The existing frozen-signature engine accepts new candidate FASTA through the queue."""
    config = replace(workbench_config, resource=str(completed_result))
    before = (completed_result / "manifest.json").read_bytes()
    job = submit_job(config=config, operation="score", uploads={"family": seed_bytes})
    run_next(config=config)
    state = read_job(config=config, job_id=job)
    assert state["status"] == "COMPLETE", state.get("error")
    path = artifact_path(config=config, job_id=job, name="signature_evidence/candidate_summary.tsv")
    summary = pd.read_csv(path, sep="\t")
    assert not summary.empty and set(summary.unit_id) == {"local_seed"}
    assert (summary.unknown_feature_count >= 0).all()
    assert (completed_result / "manifest.json").read_bytes() == before
    assert "signature_evidence/candidate_feature_hits.tsv" in state["artifacts"] or any(
        "hits" in name for name in state["artifacts"]
    )
