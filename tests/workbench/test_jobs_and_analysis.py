"""Durable queue, real bounded SQL and mocked external-tool integration tests."""

import fcntl
import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.workbench.analysis import (
    align_family,
    index_fasta,
    iter_fasta,
    run_command,
    thresholds,
)
from protein_signatures.workbench.catalogue import open_catalogue
from protein_signatures.workbench.jobs import (
    artifact_path,
    cancel_job,
    job_directory,
    list_jobs,
    read_job,
    run_next,
    submit_job,
    timestamp,
)
from protein_signatures.workbench.operations import annotate_family, build_search_database


def test_job_chain_retrieval_screening_preparation(workbench_config, seed_bytes, fake_tools):
    "Retrieve complete families, recheck every member and export provenance without relabelling."
    retrieved = submit_job(
        config=workbench_config, operation="retrieve", parameters={"cluster_ids": ["UNIPROT_P1"]}
    )
    assert read_job(config=workbench_config, job_id=retrieved)["status"] == "QUEUED"
    assert run_next(config=workbench_config) == retrieved
    assert read_job(config=workbench_config, job_id=retrieved)["status"] == "COMPLETE"
    screened = submit_job(
        config=workbench_config,
        operation="screen",
        uploads={"seeds": seed_bytes},
        sources={"family": {"job_id": retrieved, "name": "family.fasta"}},
        parameters={"deduplicate": True, "taxon_ids": [3702]},
    )
    run_next(config=workbench_config)
    state = read_job(config=workbench_config, job_id=screened)
    assert state["status"] == "COMPLETE", state.get("error")
    audit_path = artifact_path(config=workbench_config, job_id=screened, name="candidate_audit.tsv")
    audit = pd.read_csv(audit_path, sep="\t").set_index("sequence_id")
    assert audit.loc["B", "selection_status"] == "SELECTED"
    assert audit.loc["UNIPROT_P1", "selection_status"] == "EXACT_SEQUENCE_DUPLICATE"
    assert audit.loc["C", "selection_status"] == "TAXONOMY_FILTERED_OR_UNKNOWN"
    prepared = submit_job(
        config=workbench_config,
        operation="prepare",
        uploads={"seeds": seed_bytes},
        sources={"family": {"job_id": screened, "name": "candidates.fasta"}},
    )
    run_next(config=workbench_config)
    state = read_job(config=workbench_config, job_id=prepared)
    assert state["status"] == "COMPLETE", state.get("error")
    fasta = artifact_path(
        config=workbench_config, job_id=prepared, name="expanded.fasta"
    ).read_text()
    assert fasta == ">seed__local_seed\nACDEFGHIK\n"
    labels = pd.read_csv(
        artifact_path(
            config=workbench_config, job_id=prepared, name="label_assignments_to_review.tsv"
        ),
        sep="\t",
    )
    assert set(labels.curation_status) == {"PROPOSED"}
    assert labels.label_id.isna().all()
    assert "--query-cover" in fake_tools[1] and "--subject-cover" in fake_tools[1]


def test_discovery_search_job_keeps_search_evidence_separate(
    workbench_config, seed_bytes, fake_tools, tmp_path
):
    """A representative hit is explicitly a similarity hit and can remain unindexed."""
    database = tmp_path / "representatives.dmnd"
    database.write_bytes(data=b"fixture")
    config = replace(workbench_config, centroid_database=str(database))
    job = submit_job(config=config, operation="search", uploads={"seeds": seed_bytes})
    run_next(config=config)
    state = read_job(config=config, job_id=job)
    assert state["status"] == "COMPLETE", state.get("error")
    hits = pd.read_csv(
        artifact_path(config=config, job_id=job, name="representative_hits.tsv"), sep="\t"
    )
    assert hits.loc[0, "mapping_evidence"] == "SIMILARITY_TO_CENTROID"
    assert hits.loc[0, "member_count"] == 3
    assert "--very-sensitive" in fake_tools[0]
    assert "--memory-limit" not in fake_tools[0]


def test_missing_search_database_is_a_visible_failed_job(workbench_config, seed_bytes):
    """Unavailable optional dependencies must not masquerade as an empty biological result."""
    job = submit_job(config=workbench_config, operation="search", uploads={"seeds": seed_bytes})
    run_next(config=workbench_config)
    state = read_job(config=workbench_config, job_id=job)
    assert state["status"] == "FAILED" and "centroid_database" in state["error"]
    assert (job_directory(config=workbench_config, job_id=job) / "failure.log").is_file()
    with pytest.raises(InputValidationError):
        artifact_path(config=workbench_config, job_id=job, name="anything")


def test_queue_bounds_cancellation_and_exclusive_worker(workbench_config, tmp_path):
    """Duplicate workers cannot claim work and cancellation cannot rewrite a started task."""
    assert timestamp().endswith("+00:00")
    assert run_next(config=workbench_config) is None
    job = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    with (Path(workbench_config.workspace) / "worker.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert run_next(config=workbench_config) is None
        with pytest.raises(InputValidationError, match="worker is active"):
            cancel_job(config=workbench_config, job_id=job)
    cancel_job(config=workbench_config, job_id=job)
    assert read_job(config=workbench_config, job_id=job)["status"] == "CANCELLED"
    assert run_next(config=workbench_config) is None
    with pytest.raises(InputValidationError):
        cancel_job(config=workbench_config, job_id=job)
    with pytest.raises(InputValidationError):
        list_jobs(config=workbench_config, maximum=0)
    with pytest.raises(InputValidationError):
        job_directory(config=workbench_config, job_id="../../elsewhere")
    with pytest.raises(InputValidationError):
        submit_job(config=workbench_config, operation="shell")
    with pytest.raises(InputValidationError):
        submit_job(config=workbench_config, operation="align", uploads={"family": b""})
    with pytest.raises(InputValidationError):
        submit_job(config=workbench_config, operation="align", uploads={"other": b"a"})


def test_modified_inputs_profiles_and_artifacts_are_rejected(workbench_config, seed_bytes):
    """Checksums and configuration identity bind every queued task to its intended inputs."""
    job = submit_job(config=workbench_config, operation="align", uploads={"family": seed_bytes})
    (job_directory(config=workbench_config, job_id=job) / "family.fasta").write_bytes(
        data=b"changed"
    )
    run_next(config=workbench_config)
    assert "changed" in read_job(config=workbench_config, job_id=job)["error"]
    mapping = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    run_next(config=replace(workbench_config, threads=2))
    assert "changed" in read_job(config=workbench_config, job_id=mapping)["error"]
    mapping = submit_job(
        config=workbench_config, operation="resolve", parameters={"identifiers": ["B"]}
    )
    run_next(config=workbench_config)
    path = artifact_path(config=workbench_config, job_id=mapping, name="mapping.tsv")
    with pytest.raises(InputValidationError):
        artifact_path(config=workbench_config, job_id=mapping, name="../request.json")
    path.write_text(data="changed", encoding="utf-8")
    with pytest.raises(InputValidationError):
        artifact_path(config=workbench_config, job_id=mapping, name="mapping.tsv")


def test_alignment_preserves_ids_and_emits_conservation_tree(
    workbench_config, tmp_path, fake_tools
):
    """Tool-safe tip IDs must round-trip and conservation must describe the exact selection."""
    source = tmp_path / "family.fa"
    source.write_text(
        data=">a|one\nACDEFGHIK\n>b:two\nACDEFGHIL\n>c\nACDEFGHIM\n", encoding="utf-8"
    )
    output = tmp_path / "alignment"
    output.mkdir()
    align_family(config=workbench_config, source=source, parameters={"tree": True}, output=output)
    assert ">b:two" in (output / "alignment.fasta").read_text()
    conservation = pd.read_csv(output / "conservation.tsv", sep="\t")
    assert conservation.loc[0, "consensus_fraction_non_gap"] == 1
    assert conservation.loc[8, "consensus_fraction_non_gap"] == pytest.approx(1 / 3)
    assert (output / "tree.nwk").read_text().endswith(";\n")
    assert (output / "tip_labels.tsv").is_file()


@pytest.mark.parametrize(
    "parameters", [{"selected_ids": ["missing"]}, {"selected_ids": ["a"]}, {"tree": "yes"}]
)
def test_alignment_rejects_invalid_or_insufficient_selection(
    workbench_config, tmp_path, fake_tools, parameters
):
    """Selections may not silently drop requested proteins or coerce scientific controls."""
    source = tmp_path / "input.fa"
    source.write_text(data=">a\nACD\n>b\nACE\n", encoding="utf-8")
    output = tmp_path / "out"
    output.mkdir()
    with pytest.raises(InputValidationError):
        align_family(config=workbench_config, source=source, parameters=parameters, output=output)


def test_pfam_annotation_and_database_build(workbench_config, tmp_path, fake_tools):
    """Local annotation records domain coordinates and explicit no-hit assessments."""
    source = tmp_path / "input.fa"
    source.write_text(data=">B\nACDEFGHIK\n>C\nACDEFGHIL\n", encoding="utf-8")
    hmm = tmp_path / "Pfam-A.hmm"
    hmm.write_text(data="fixture", encoding="utf-8")
    output = tmp_path / "annotation"
    output.mkdir()
    config = replace(workbench_config, pfam_hmm=str(hmm), centroid_fasta=str(source))
    annotate_family(config=config, source=source, output=output)
    domains = pd.read_csv(output / "domains.tsv", sep="\t").set_index("protein_id")
    assert domains.loc["B", "domain_id"] == "PF00646"
    assert domains.loc["C", "assessment_status"] == "ASSESSED_NO_HIT"
    database_output = tmp_path / "database"
    database_output.mkdir()
    build_search_database(config=config, output=database_output)
    assert (database_output / "representatives.dmnd").is_file()
    for function in (build_search_database,):
        with pytest.raises(InputValidationError):
            function(config=workbench_config, output=database_output)


def test_streaming_input_validation_and_thresholds(workbench_config, tmp_path):
    """Malformed input, duplicates and non-finite controls are rejected before analysis."""
    source = tmp_path / "input.fa"
    for value in ("ACD\n", ">\nACD\n", "", ">a\nAC-\n"):
        source.write_text(data=value, encoding="utf-8")
        with pytest.raises(InputValidationError):
            list(iter_fasta(path=source))
    source.write_text(data=">a\nACD\n>a\nACE\n", encoding="utf-8")
    with open_catalogue(config=workbench_config, with_index=False) as connection:
        with pytest.raises(InputValidationError, match="duplicate"):
            index_fasta(connection=connection, path=source, name="family", maximum=10)
    with pytest.raises(InputValidationError):
        thresholds(parameters={"identity": float("nan")})
    assert thresholds(parameters={}) == (30, 50, 1e-5)


def test_command_errors_are_audited_and_timeout_is_bounded(workbench_config, tmp_path):
    """External failures retain an argument audit and cannot hang the worker indefinitely."""
    Path(workbench_config.workspace, "temporary").mkdir(parents=True)
    with pytest.raises(ExternalToolError, match="unavailable"):
        run_command(config=workbench_config, arguments=["/does/not/exist"], output=tmp_path)
    with pytest.raises(ExternalToolError, match="status 3"):
        run_command(
            config=workbench_config,
            arguments=[sys.executable, "-c", "raise SystemExit(3)"],
            output=tmp_path,
        )
    with pytest.raises(ExternalToolError, match="timed out"):
        run_command(
            config=replace(workbench_config, command_timeout=0.01),
            arguments=[sys.executable, "-c", "import time; time.sleep(10)"],
            output=tmp_path,
        )
    assert (tmp_path / "commands.jsonl").is_file()
