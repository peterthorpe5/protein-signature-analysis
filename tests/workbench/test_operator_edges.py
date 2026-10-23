"""Operator-facing recovery, stale-request and CLI lifecycle contracts."""

import fcntl
import os
import runpy
import signal
import socket
import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from protein_signature_app import workbench_viewers
from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic
from protein_signatures.workbench import cli, jobs, operations, seeds
from protein_signatures.workbench.config import configure_denbi
from protein_signatures.workbench.jobs import (
    artifact_path,
    job_directory,
    read_job,
    recover_interrupted_job,
    run_next,
    submit_job,
)


def test_cli_import_and_continuous_worker_interrupt(config_path, tmp_path, monkeypatch):
    """A long-lived worker exits cleanly on interruption and restores signal handlers."""
    destination = tmp_path / "new.toml"
    assert (
        cli.main(
            argv=[
                "--config",
                str(destination),
                "--action",
                "configure-denbi",
                "--source-config",
                str(config_path),
                "--workspace",
                str(tmp_path / "independent"),
            ]
        )
        == 0
    )
    old = signal.getsignal(signal.SIGTERM)
    monkeypatch.setattr(
        cli.time, "sleep", lambda *_args: (_ for _ in ()).throw(KeyboardInterrupt())
    )
    assert cli.main(argv=["--config", str(config_path), "--action", "worker"]) == 130
    assert signal.getsignal(signal.SIGTERM) == old
    with pytest.raises(InputValidationError):
        configure_denbi(
            source=tmp_path / "missing",
            destination=tmp_path / "other.toml",
            workspace=tmp_path / "other",
        )


def test_cli_source_job_features_and_upload_bounds(config_path, workbench_config, tmp_path):
    """CLI references retain source checksums and oversized uploads are rejected before reading."""
    job = submit_job(
        config=workbench_config, operation="retrieve", parameters={"cluster_ids": ["UNIPROT_P1"]}
    )
    run_next(config=workbench_config)
    assert (
        cli.main(
            argv=[
                "--config",
                str(config_path),
                "--action",
                "submit",
                "--operation",
                "score",
                "--family-job",
                job,
                "--family-artifact",
                "family.fasta",
                "--seed-job",
                job,
                "--seed-artifact",
                "family.fasta",
                "--features-job",
                job,
                "--features-artifact",
                "members.tsv",
            ]
        )
        == 0
    )
    path = tmp_path / "large.fa"
    with path.open(mode="wb") as handle:
        handle.truncate(21 * 1024 * 1024)
    assert (
        cli.main(
            argv=[
                "--config",
                str(config_path),
                "--action",
                "submit",
                "--operation",
                "align",
                "--family-fasta",
                str(path),
            ]
        )
        == 2
    )


def test_module_entry_point_uses_named_actions(config_path, monkeypatch):
    """The documented python -m route behaves like the installed entry point."""
    monkeypatch.setattr(
        sys, "argv", ["workbench", "--config", str(config_path), "--action", "doctor"]
    )
    with pytest.raises(SystemExit) as captured:
        runpy.run_module(mod_name="protein_signatures.workbench.cli", run_name="__main__")
    assert captured.value.code == 0


def test_corrupt_job_and_changed_request_or_package_fail_closed(workbench_config, monkeypatch):
    """Changing requests or running a new package cannot silently alter queued analyses."""
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    path = job_directory(config=workbench_config, job_id=job)
    original = (path / "request.json").read_text()
    (path / "request.json").write_text(data=original + "\n", encoding="utf-8")
    run_next(config=workbench_config)
    assert "request changed" in read_job(config=workbench_config, job_id=job)["error"]
    other = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    monkeypatch.setattr(jobs, "__version__", "999.0")
    run_next(config=workbench_config)
    assert "Package version changed" in read_job(config=workbench_config, job_id=other)["error"]
    (path / "state.json").write_text(data='{"job_id":"wrong"}', encoding="utf-8")
    with pytest.raises(InputValidationError, match="unavailable"):
        read_job(config=workbench_config, job_id=job)


def test_sources_changed_during_execution_and_malformed_operation(workbench_config, monkeypatch):
    """Mid-job source changes and undeclared controls prevent completed publication."""
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"unrecognised": True}
    )
    run_next(config=workbench_config)
    assert "Unsupported" in read_job(config=workbench_config, job_id=job)["error"]
    other = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    monkeypatch.setattr(
        operations,
        "execute_job",
        lambda **_kwargs: Path(workbench_config.sequence_map).write_bytes(data=b"changed"),
    )
    run_next(config=workbench_config)
    assert "while the job was running" in read_job(config=workbench_config, job_id=other)["error"]


def test_checksum_failure_with_unchanged_size_and_symlink_escape(workbench_config, tmp_path):
    """File length alone cannot authenticate a reused analysis artifact."""
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    run_next(config=workbench_config)
    path = artifact_path(config=workbench_config, job_id=job, name="mapping.tsv")
    payload = path.read_bytes()
    path.write_bytes(data=b"X" + payload[1:])
    with pytest.raises(InputValidationError, match="checksum"):
        artifact_path(config=workbench_config, job_id=job, name="mapping.tsv")
    root = Path(workbench_config.workspace) / "jobs"
    invalid = "job_" + "a" * 32
    (root / invalid).symlink_to(target=tmp_path, target_is_directory=True)
    with pytest.raises(InputValidationError, match="escapes"):
        job_directory(config=workbench_config, job_id=invalid)


def test_recovery_permission_and_host_uncertainty_are_not_guessed(workbench_config, monkeypatch):
    """Recovery requires positive evidence that the recorded worker and tools have stopped."""
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    path = job_directory(config=workbench_config, job_id=job)
    state = read_job(config=workbench_config, job_id=job)
    state.update(status="RUNNING", worker_pid=2_000_000_000, worker_host=socket.gethostname())
    write_json_atomic(path=path / "state.json", value=state)
    with (Path(workbench_config.workspace) / "worker.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(InputValidationError, match="active"):
            recover_interrupted_job(config=workbench_config, job_id=job)
    monkeypatch.setattr(os, "kill", lambda *_args: (_ for _ in ()).throw(PermissionError()))
    with pytest.raises(InputValidationError, match="Cannot establish"):
        recover_interrupted_job(config=workbench_config, job_id=job)
    monkeypatch.setattr(os, "kill", lambda *_args: (_ for _ in ()).throw(ProcessLookupError()))
    (path / "outputs").mkdir()
    write_json_atomic(
        path=path / "outputs" / "active_process.json", value={"process_group": 2_000_000_000}
    )
    monkeypatch.setattr(os, "killpg", lambda *_args: (_ for _ in ()).throw(PermissionError()))
    with pytest.raises(InputValidationError, match="tool group has stopped"):
        recover_interrupted_job(config=workbench_config, job_id=job)
    monkeypatch.setattr(os, "killpg", lambda *_args: (_ for _ in ()).throw(ProcessLookupError()))
    recover_interrupted_job(config=workbench_config, job_id=job)
    assert read_job(config=workbench_config, job_id=job)["status"] == "FAILED"


def test_reference_seed_payload_limit_and_bad_tree_parser(tmp_path, monkeypatch):
    """The browser's payload bound and invalid tree syntax remain explicit failures."""
    calls = iter(
        [
            pd.DataFrame(data={"protein_id": ["p1"]}),
            pd.DataFrame(data={"protein_id": ["p1"], "sequence": ["ACD"]}),
        ]
    )
    monkeypatch.setattr(seeds, "query_dataframe", lambda **_kwargs: next(calls))
    monkeypatch.setattr(seeds, "MAX_UPLOAD_BYTES", 1)
    with pytest.raises(InputValidationError, match="20 MiB"):
        seeds.export_seeds(database=tmp_path / "db", selection="protein", values=["p1"])
    with pytest.raises(InputValidationError, match="parse"):
        workbench_viewers.tree_figure(newick="((a,b,c);", labels={"a": "A", "b": "B", "c": "C"})


def test_build_and_annotation_operations_are_queued(
    workbench_config, seed_bytes, fake_tools, tmp_path
):
    """Optional local analyses use the same durable dispatch and publication contracts."""
    centroid = tmp_path / "centroids.fasta"
    centroid.write_bytes(data=seed_bytes)
    hmm = tmp_path / "Pfam-A.hmm"
    hmm.write_bytes(data=b"fixture")
    config = replace(workbench_config, centroid_fasta=str(centroid), pfam_hmm=str(hmm))
    job = submit_job(config=config, operation="build_search")
    run_next(config=config)
    assert read_job(config=config, job_id=job)["status"] == "COMPLETE"
    annotation = submit_job(
        config=config, operation="annotate", uploads={"family": b">B\nACDEFGHIK\n"}
    )
    run_next(config=config)
    assert read_job(config=config, job_id=annotation)["status"] == "COMPLETE"
    assert artifact_path(config=config, job_id=annotation, name="query_features.tsv").is_file()
