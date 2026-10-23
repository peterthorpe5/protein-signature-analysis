"""Operator profiles, named CLI actions and canonical seed selection contracts."""

import json
from dataclasses import replace
from pathlib import Path

import duckdb
import pytest

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.cli import build_parser, main
from protein_signatures.workbench.config import (
    bounded_number,
    config_identity,
    configure_denbi,
    diagnostics,
    load_config,
)
from protein_signatures.workbench.jobs import list_jobs, read_job, run_next, submit_job
from protein_signatures.workbench.seeds import export_seeds


def test_load_configuration_and_capability_report(config_path, workbench_config):
    """Installed resource paths remain explicit and diagnostics do not execute analyses."""
    config = load_config(path=config_path)
    assert config == workbench_config
    assert config_identity(config=config) == config_identity(config=workbench_config)
    assert config_identity(config=replace(config, threads=2)) != config_identity(config=config)
    report = {row["capability"]: row for row in diagnostics(config=config)}
    assert report["index"]["status"] == "READY"
    assert report["centroid_database"]["status"] == "UNAVAILABLE"


@pytest.mark.parametrize(
    "content",
    [
        "",
        'workspace="x"\nunknown="x"',
        'workspace="x"\nthreads=0',
        'workspace="x"\nthreads=true',
        'workspace="x"\nidentifier_prefixes=[1]',
        'workspace="x"\nsequences=1',
        'workspace="x"\nresource="x/result"',
        'workspace="x"\nindex="x/source.parquet"',
        "invalid TOML",
    ],
)
def test_invalid_profiles_fail_before_mutation(tmp_path, content):
    """Malformed limits and overlapping immutable/writable paths must fail explicitly."""
    path = tmp_path / "invalid.toml"
    path.write_text(data=content, encoding="utf-8")
    with pytest.raises(InputValidationError):
        load_config(path=path)
    assert not (tmp_path / "x").exists()


def test_configuration_import_preserves_e3_and_rejects_existing_destination(config_path, tmp_path):
    """Import shares read-only sources but creates independent signature job storage."""
    original = config_path.read_bytes()
    destination = tmp_path / "new" / "workbench.toml"
    result = configure_denbi(
        source=config_path, destination=destination, workspace=tmp_path / "separate"
    )
    config = load_config(path=result)
    assert config.workspace == str(tmp_path / "separate")
    assert config_path.read_bytes() == original
    with pytest.raises(InputValidationError, match="already exists"):
        configure_denbi(source=config_path, destination=destination, workspace=tmp_path / "other")
    with pytest.raises(InputValidationError, match="separate workspace"):
        configure_denbi(
            source=config_path,
            destination=tmp_path / "duplicate.toml",
            workspace=Path(load_config(path=config_path).workspace),
        )


@pytest.mark.parametrize("value", [None, "bad", True, float("inf"), -1, 1.5])
def test_numeric_defences(value):
    """All worker controls reject non-finite, Boolean and non-integral integer values."""
    with pytest.raises(InputValidationError):
        bounded_number(value=value, name="limit", minimum=1, maximum=10, integer=True)
    assert bounded_number(value="3", name="limit", minimum=1, maximum=10, integer=True) == 3


def test_named_cli_can_submit_execute_and_report_mapping(config_path, tmp_path, capsys):
    """The production command-line route exercises the same queue as the GUI."""
    assert (
        build_parser().parse_args(args=["--config", str(config_path), "--action", "doctor"]).action
        == "doctor"
    )
    assert main(argv=["--config", str(config_path), "--action", "doctor"]) == 0
    assert "capability\tstatus\tdetail" in capsys.readouterr().out
    parameters = tmp_path / "parameters.json"
    parameters.write_text(data=json.dumps(obj={"identifiers": ["B", "missing"]}), encoding="utf-8")
    base = ["--config", str(config_path)]
    assert (
        main(
            argv=[
                *base,
                "--action",
                "submit",
                "--operation",
                "resolve",
                "--parameters",
                str(parameters),
            ]
        )
        == 0
    )
    identifier = capsys.readouterr().out.strip()
    assert main(argv=[*base, "--action", "worker", "--once"]) == 0
    assert read_job(config=load_config(path=config_path), job_id=identifier)["status"] == "COMPLETE"
    assert main(argv=[*base, "--action", "status"]) == 0
    assert "COMPLETE" in capsys.readouterr().out
    assert main(argv=[*base, "--action", "worker", "--once"]) == 0
    assert main(argv=["--config", str(tmp_path / "missing"), "--action", "doctor"]) == 2


def test_cli_upload_source_and_cancellation_validation(config_path, seed_bytes, tmp_path):
    """Named FASTA input is copied, missing input partners fail and waiting jobs cancel."""
    fasta = tmp_path / "seed.fa"
    fasta.write_bytes(data=seed_bytes)
    base = ["--config", str(config_path)]
    assert (
        main(
            argv=[*base, "--action", "submit", "--operation", "align", "--family-fasta", str(fasta)]
        )
        == 0
    )
    config = load_config(path=config_path)
    identifier = list_jobs(config=config)[0]["job_id"]
    assert main(argv=[*base, "--action", "cancel", "--job-id", identifier]) == 0
    assert (
        main(argv=[*base, "--action", "submit", "--operation", "align", "--family-job", identifier])
        == 2
    )
    assert main(argv=[*base, "--action", "cancel"]) == 2


def test_canonical_seed_selectors_are_complete_bounded_and_exact(tmp_path):
    """HOGs, labels and features export actual protein sequences, never global guesses."""
    database = tmp_path / "reference.duckdb"
    with duckdb.connect(database=str(database)) as connection:
        connection.execute(
            (
                "CREATE TABLE proteins AS SELECT 'p1' AS protein_id,'ACD' AS sequence "
                "UNION ALL SELECT 'p2','ACE'"
            )
        )
        connection.execute(
            (
                "CREATE TABLE orthofinder_memberships AS SELECT 'HOG1' AS "
                "group_id,protein_id FROM proteins"
            )
        )
        connection.execute(
            (
                "CREATE TABLE label_memberships AS SELECT 'class:a' AS "
                "label_id,protein_id FROM proteins"
            )
        )
        connection.execute(
            "CREATE TABLE features AS SELECT 'k3:ACD' AS feature_id,'p1' AS protein_id"
        )
    assert export_seeds(database=database, selection="protein", values=["p1"]) == b">p1\nACD\n"
    assert export_seeds(database=database, selection="group", values=["HOG1"]).count(b">") == 2
    assert export_seeds(database=database, selection="label", values=["class:a"]).count(b">") == 2
    assert export_seeds(database=database, selection="feature", values=["k3:ACD"]).count(b">") == 1
    with pytest.raises(InputValidationError, match="limit"):
        export_seeds(database=database, selection="group", values=["HOG1"], maximum=1)
    with pytest.raises(InputValidationError, match="not found"):
        export_seeds(database=database, selection="protein", values=["missing"])
    with pytest.raises(InputValidationError):
        export_seeds(database=database, selection="invalid", values=["p1"])


def test_score_requires_a_reference_and_preserves_unknown_evidence(workbench_config, seed_bytes):
    """A missing reference result is an explicit failed job, never an invented score."""
    identifier = submit_job(
        config=workbench_config, operation="score", uploads={"family": seed_bytes}
    )
    run_next(config=workbench_config)
    state = read_job(config=workbench_config, job_id=identifier)
    assert state["status"] == "FAILED" and "reference result" in state["error"]
