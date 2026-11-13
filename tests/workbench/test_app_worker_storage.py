"""Exercise exact job selection, private EBI contacts and conservative scratch retention."""

import fcntl
import os
import stat
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.io_utils import write_json_atomic
from protein_signatures.workbench import operations, providers
from protein_signatures.workbench.analysis import run_command
from protein_signatures.workbench.jobs import (
    artifact_path,
    job_directory,
    read_job,
    run_next,
    submit_job,
)
from protein_signatures.workbench.local_worker import start_local_job
from protein_signatures.workbench.storage import clear_scratch, preview_cleanup, temporary_root


def test_exact_worker_and_app_child(workbench_config, config_path):
    """A real app-host child starts only the selected request, bypassing the normal queue."""
    queued = submit_job(
        config=workbench_config,
        operation="resolve",
        parameters={"identifiers": ["B"]},
    )
    app = submit_job(
        config=workbench_config,
        operation="resolve",
        parameters={"identifiers": ["C"]},
        execution="app",
    )
    assert start_local_job(config=workbench_config, config_path=config_path, job_id=app) > 0
    deadline = time.monotonic() + 8
    state = read_job(config=workbench_config, job_id=app)
    while state["status"] in {"QUEUED", "RUNNING"} and time.monotonic() < deadline:
        time.sleep(0.05)
        state = read_job(config=workbench_config, job_id=app)
    assert state["status"] == "COMPLETE", state.get("error")
    assert read_job(config=workbench_config, job_id=queued)["status"] == "QUEUED"
    assert run_next(config=workbench_config) == queued
    with pytest.raises(InputValidationError, match="queued"):
        start_local_job(config=workbench_config, config_path=config_path, job_id=app)


def test_normal_worker_skips_reserved_request(workbench_config):
    """A normal worker does not claim an app reservation, which can be started explicitly."""
    reserved = submit_job(
        config=workbench_config,
        operation="resolve",
        execution="app",
        parameters={"identifiers": ["B"]},
    )
    assert run_next(config=workbench_config) is None
    assert run_next(config=workbench_config, job_id=reserved) == reserved
    assert read_job(config=workbench_config, job_id=reserved)["status"] == "COMPLETE"


def test_private_interpro_contact_is_not_published(workbench_config, monkeypatch):
    """Pass the contact to the opted-in provider and keep it outside result files."""
    observed = []

    def provider(*, config, proteins, output):
        """Capture the request without making a network call."""
        observed.append((config.interpro_email, tuple(proteins.values())))
        return []

    monkeypatch.setattr(providers, "remote_interpro", provider)
    email = b"contact@example.org"
    job = submit_job(
        config=workbench_config,
        operation="annotate_proteins",
        uploads={"family": b">one\nACDEFGHIK\n", "ebi_email": email},
        parameters={"engine": "interpro_remote", "remote_allowed": True},
    )
    assert run_next(config=workbench_config) == job
    state = read_job(config=workbench_config, job_id=job)
    assert state["status"] == "COMPLETE", state.get("error")
    assert observed == [(email.decode(), ("ACDEFGHIK",))]
    root = job_directory(config=workbench_config, job_id=job)
    assert stat.S_IMODE((root / "ebi_email.private").stat().st_mode) == 0o600
    assert email not in (root / "request.json").read_bytes()
    for name in state["artifacts"]:
        assert (
            email not in artifact_path(config=workbench_config, job_id=job, name=name).read_bytes()
        )


def test_scratch_preview_rechecks_changes_and_preserves_results(
    workbench_config,
    tmp_path,
    monkeypatch,
):
    """Respect retention, reject links and changed previews, and preserve checked artifacts."""

    def sample(*, config, request, output):
        """Write private scratch and an independent published result."""
        (temporary_root(config=config) / "disposable.idx").write_bytes(b"fixture scratch")
        (output / "retained.tsv").write_bytes(b"id\tlength\none\t9\n")

    monkeypatch.setattr(operations, "execute_job", sample)
    identifier = submit_job(config=workbench_config, operation="resolve")
    assert run_next(config=workbench_config) == identifier
    assert not preview_cleanup(config=workbench_config, retention_days=1)
    root = job_directory(config=workbench_config, job_id=identifier)
    scratch = root / "scratch"
    assert scratch.is_dir()
    path = scratch / "disposable.idx"
    old = datetime.now(tz=timezone.utc) - timedelta(days=40)
    write_json_atomic(
        path=root / "state.json",
        value={
            **read_job(config=workbench_config, job_id=identifier),
            "finished_at": old.isoformat(),
        },
    )
    for item in (scratch, *scratch.iterdir()):
        os.utime(item, (old.timestamp(), old.timestamp()))
    first = preview_cleanup(config=workbench_config)
    assert len(first) == 1 and first[0].bytes == path.stat().st_size
    with Path(workbench_config.workspace, "worker.lock").open(mode="a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with pytest.raises(InputValidationError, match="active worker"):
            clear_scratch(config=workbench_config, entries=first)
    outside = tmp_path / "must_keep"
    outside.write_bytes(b"original data")
    link = scratch / "unsafe"
    link.symlink_to(outside)
    os.utime(scratch, (old.timestamp(), old.timestamp()))
    assert not preview_cleanup(config=workbench_config)
    link.unlink()
    os.utime(scratch, (old.timestamp(), old.timestamp()))
    path.write_bytes(b"changed bytes")
    os.utime(path, (old.timestamp(), old.timestamp()))
    with pytest.raises(InputValidationError, match="changed"):
        clear_scratch(config=workbench_config, entries=first)
    new = preview_cleanup(config=workbench_config)
    assert clear_scratch(config=workbench_config, entries=new) == len(b"changed bytes")
    assert not scratch.exists()
    assert outside.read_bytes() == b"original data"
    kept = artifact_path(config=workbench_config, job_id=identifier, name="retained.tsv")
    assert kept.read_bytes() == b"id\tlength\none\t9\n"


def test_bounded_command_response(workbench_config, tmp_path):
    """Stream a finite response and terminate a command that exceeds its byte limit."""
    source = tmp_path / "input"
    source.write_bytes(b"0123456789")
    output = tmp_path / "output"
    output.mkdir()
    reader = [sys.executable, "-c", "import sys; sys.stdout.buffer.write(sys.stdin.buffer.read())"]
    response = tmp_path / "response"
    run_command(
        config=workbench_config,
        arguments=reader,
        output=output,
        stdin_path=source,
        stdout_path=response,
        stdout_limit=10,
    )
    assert response.read_bytes() == source.read_bytes()
    with pytest.raises(ExternalToolError, match="size limit"):
        run_command(
            config=workbench_config,
            arguments=reader,
            output=output,
            stdin_path=source,
            stdout_path=tmp_path / "limited",
            stdout_limit=9,
        )
    assert not (output / "active_process.json").exists()
