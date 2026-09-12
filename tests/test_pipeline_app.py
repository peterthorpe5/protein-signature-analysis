"""End-to-end pipeline, publication, CLI and application-backend tests."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

import duckdb
import pytest

import protein_signatures.cli as cli_module
from protein_signature_app import backend, launcher
from protein_signatures.cli import build_parser, main
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.pipeline import run_campaign, validate_campaign
from protein_signatures.publication import (
    verify_completed_result,
    verify_input_authorities,
)
from protein_signatures.schemas import schema_for, table_schemas


def test_offline_example_publishes_all_authorities(completed_result: Path) -> None:
    """The packaged campaign should publish portable TSV, Parquet and DuckDB evidence."""

    verify_completed_result(result_dir=completed_result)
    verify_input_authorities(result_dir=completed_result)
    metadata = json.loads((completed_result / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["campaign"]["campaign"]["campaign_id"] == "minimal_e3_fbox_demo"
    assert metadata["profile"]["require_structural_evidence"] is True
    assert (
        "mechanism- or component-role-specific"
        in metadata["profile"]["default_comparison_description"]
    )
    assert metadata["evidence_availability"]["pfam"] == "COMPLETE"
    assert metadata["counts"]["proteins"] == 24
    assert metadata["counts"]["signatures"] > 0
    assert metadata["counts"]["orthofinder_group_context"] == 0
    assert metadata["profile_structural_evidence"]["status"] == "COMPLETE"
    assert metadata["profile_structural_evidence"]["eligible_evidence_item_count"] == 15
    assert metadata["human_reports"]["status"] == "COMPLETE"
    assert metadata["human_reports"]["excel_workbook_count"] >= len(table_schemas())
    manifest = json.loads((completed_result / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "COMPLETE"
    expected_tables = set(table_schemas())
    with duckdb.connect(
        str(completed_result / "protein_signatures.duckdb"), read_only=True
    ) as connection:
        relations = {row[0] for row in connection.execute("SHOW TABLES").fetchall()}
        assert expected_tables <= relations
        assert "signature_evidence" in relations
        assert connection.execute("SELECT count(*) FROM proteins").fetchone()[0] == 24
        assert (
            connection.execute(
                "SELECT count(*) FROM signatures "
                "WHERE evidence_class = 'DECISION_CANDIDATE__VALIDATED_STUDY_WIDE'"
            ).fetchone()[0]
            > 0
        )
        assert connection.execute(
            "SELECT default_analysis, default_background_label_id "
            "FROM profile_labels WHERE label_id = ?",
            ("e3:ubiquitin:crl:crl1_scf:f_box",),
        ).fetchone() == (True, "control:matched_substrate_receptor_reference")
    for table_name in expected_tables:
        assert (completed_result / "tables" / f"{table_name}.tsv").is_file()
        assert (completed_result / "tables" / f"{table_name}.parquet").is_file()
    assert (completed_result / "analysis/99_final_results/tables/signatures.xlsx").is_file()
    assert (
        completed_result
        / "analysis/05_association_statistics/figures/00_signature_evidence_classes.pdf"
    ).is_file()
    assert (completed_result / "analysis/00_run_information/report_inventory.tsv").is_file()


def test_validate_and_resume_are_deterministic(tmp_path: Path, example_dir: Path) -> None:
    """Validation should be read-only and resume should verify the same completed run."""

    summary = validate_campaign(config_path=example_dir / "campaign.yaml")
    assert summary == {
        "status": "VALID",
        "campaign_id": "minimal_e3_fbox_demo",
        "profile_id": "e3_default",
        "protein_count": 24,
        "assignment_count": 24,
        "external_feature_record_count": 0,
        "domain_hit_count": 12,
        "domain_assessment_count": 24,
        "redundancy_membership_count": 0,
        "structure_count": 8,
        "structure_comparison_count": 7,
        "alphafold_request_count": 0,
        "comparison_count": 1,
        "alphafold_enabled": False,
        "foldseek_preflight": {
            "status": "NOT_SELECTED",
            "tool_version": "",
            "candidate_model_count": 0,
            "candidate_protein_count": 0,
            "maximum_hits": 1000,
        },
        "explainable_ml_enabled": True,
        "orthofinder_version": "",
        "orthofinder_source_mode": "",
        "structural_alignment_resource": "NOT_SELECTED",
        "profile_structural_evidence": {
            "required": True,
            "status": "COMPLETE",
            "eligible_coordinate_model_count": 0,
            "eligible_fold_assignment_count": 8,
            "complete_structure_comparison_count": 7,
            "imported_structural_feature_count": 0,
            "pending_alphafold_request_count": 0,
            "pending_structural_search_count": 0,
            "completed_structural_search_count": 0,
            "eligible_evidence_item_count": 15,
        },
    }
    destination = tmp_path / "campaign"
    first = run_campaign(
        config_path=example_dir / "campaign.yaml",
        output_dir=destination,
        threads=2,
    )
    assert (
        run_campaign(
            config_path=example_dir / "campaign.yaml",
            output_dir=destination,
            threads=99,
            resume=True,
        )
        == first
    )
    with pytest.raises(PublicationError):
        run_campaign(
            config_path=example_dir / "campaign.yaml",
            output_dir=destination,
        )
    with pytest.raises(InputValidationError):
        run_campaign(
            config_path=example_dir / "campaign.yaml",
            output_dir=tmp_path / "bad_threads",
            threads=0,
        )


def test_resume_and_verification_reject_tampering(
    completed_result: Path, example_dir: Path
) -> None:
    """A changed output or input must invalidate verification and computational resume."""

    protein_table = completed_result / "tables" / "proteins.tsv"
    protein_table.write_text(
        protein_table.read_text(encoding="utf-8") + "tamper\n", encoding="utf-8"
    )
    with pytest.raises(PublicationError, match="size mismatch"):
        verify_completed_result(result_dir=completed_result)
    with pytest.raises(PublicationError):
        run_campaign(
            config_path=example_dir / "campaign.yaml",
            output_dir=completed_result,
            resume=True,
        )


def test_resume_rejects_changed_run_identity(tmp_path: Path, example_dir: Path) -> None:
    """A valid result cannot be reused for a changed configuration identity."""

    destination = run_campaign(
        config_path=example_dir / "campaign.yaml",
        output_dir=tmp_path / "result",
    )
    metadata_path = destination / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["run_identity_sha256"] = "0" * 64
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    _refresh_output_checksum(result_dir=destination, relative="run_metadata.json")
    with pytest.raises(InputValidationError, match="different configuration"):
        run_campaign(
            config_path=example_dir / "campaign.yaml",
            output_dir=destination,
            resume=True,
        )


def test_app_backend_reads_only_canonical_tables(completed_result: Path) -> None:
    """The app backend should validate resources and parameterise read-only queries."""

    database = backend.resolve_database(resource=completed_result)
    assert backend.resolve_database(resource=database) == database
    assert backend.table_count(database=database, table_name="proteins") == 24
    assert backend.distinct_values(
        database=database, table_name="signatures", column_name="comparison_id"
    ) == ("fbox_vs_matched_receptor_reference",)
    metadata = backend.load_metadata(database=database)
    assert metadata["package"] == "protein-signature-analysis"
    frame = backend.query_dataframe(
        database=database,
        sql="SELECT protein_id FROM proteins WHERE protein_id = ?",
        parameters=("fbox01",),
    )
    assert frame.iloc[0]["protein_id"] == "fbox01"
    with pytest.raises(InputValidationError, match="read-only"):
        backend.query_dataframe(database=database, sql="DELETE FROM proteins")
    with pytest.raises(InputValidationError, match="comment-free"):
        backend.query_dataframe(
            database=database,
            sql="SELECT protein_id FROM proteins; SELECT protein_id FROM proteins",
        )
    with pytest.raises(InputValidationError, match="non-canonical"):
        backend.query_dataframe(database=database, sql="SELECT * FROM duckdb_tables()")
    with pytest.raises(InputValidationError, match="canonical result relation"):
        backend.query_dataframe(database=database, sql="SELECT read_blob(?)", parameters=("x",))
    with pytest.raises(InputValidationError, match="Unknown canonical"):
        backend.table_count(database=database, table_name="missing")
    with pytest.raises(InputValidationError, match="Unknown canonical"):
        backend.distinct_values(database=database, table_name="proteins", column_name="missing")


def test_app_backend_exposes_every_complete_canonical_download(
    completed_result: Path,
) -> None:
    """Every canonical table should have full TSV/XLSX assets and a bounded preview."""

    database = backend.resolve_database(resource=completed_result)
    names = backend.canonical_table_names()
    assert names == tuple(table_schemas())
    assert len(names) == 26
    for table_name in names:
        count = backend.table_count(database=database, table_name=table_name)
        preview = backend.canonical_table_preview(
            database=database,
            table_name=table_name,
            limit=5,
            offset=0,
        )
        assert tuple(preview.columns) == tuple(table_schemas()[table_name].names)
        assert len(preview) == min(5, count)
        assets = backend.load_canonical_table_assets(
            database=database,
            table_name=table_name,
        )
        assert tuple(asset.file_format for asset in assets) == ("TSV", "XLSX")
        assert all(asset.row_count == count for asset in assets)
        assert assets[0].payload.startswith(table_schemas()[table_name].names[0].encode())
        assert assets[1].payload.startswith(b"PK")
    inventory_assets = backend.load_report_inventory_assets(database=database)
    assert tuple(asset.file_format for asset in inventory_assets) == ("TSV", "XLSX")
    assert {asset.row_count for asset in inventory_assets} == {181}


@pytest.mark.parametrize(
    ("limit", "offset", "message"),
    [
        (0, 0, "limit"),
        (5_001, 0, "limit"),
        (True, 0, "limit"),
        (5, -1, "offset"),
        (5, True, "offset"),
    ],
)
def test_canonical_preview_rejects_invalid_bounds(
    completed_result: Path,
    limit: int,
    offset: int,
    message: str,
) -> None:
    """Canonical previews should reject booleans and out-of-contract pagination."""

    database = completed_result / "protein_signatures.duckdb"
    with pytest.raises(InputValidationError, match=message):
        backend.canonical_table_preview(
            database=database,
            table_name="proteins",
            limit=limit,
            offset=offset,
        )


def test_app_backend_rejects_invalid_resources(tmp_path: Path) -> None:
    """The app should fail closed on missing, unsupported and malformed resources."""

    with pytest.raises(InputValidationError, match="not a valid completed"):
        backend.resolve_database(resource=tmp_path / "missing")
    wrong = tmp_path / "wrong.duckdb"
    wrong.write_bytes(b"x")
    with pytest.raises(InputValidationError):
        backend.resolve_database(resource=wrong)
    metadata = tmp_path / "run_metadata.json"
    metadata.write_text("[]", encoding="utf-8")
    with pytest.raises(InputValidationError, match="JSON object"):
        backend.load_metadata(database=tmp_path / "database.duckdb")
    metadata.write_text("{", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Could not load"):
        backend.load_metadata(database=tmp_path / "database.duckdb")
    with pytest.raises(InputValidationError, match="query failed"):
        backend.query_dataframe(
            database=tmp_path / "missing.duckdb",
            sql="SELECT * FROM proteins",
        )


def test_app_backend_cannot_execute_external_file_views(tmp_path: Path) -> None:
    """A canonical-named DuckDB view must not read an external local file."""

    secret = tmp_path / "outside.txt"
    secret.write_text("must-not-be-read\n", encoding="utf-8")
    database = tmp_path / "malicious.duckdb"
    escaped = str(secret).replace("'", "''")
    with duckdb.connect(str(database)) as connection:
        connection.execute(
            f"CREATE VIEW proteins AS SELECT content AS protein_id FROM read_text('{escaped}')"
        )
    with pytest.raises(InputValidationError, match="DuckDB query failed"):
        backend.query_dataframe(database=database, sql="SELECT * FROM proteins")


def test_app_backend_rejects_missing_database_and_wrong_result_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A verified root still requires its one canonical database path."""

    monkeypatch.setattr(backend, "verify_completed_result", lambda **_kwargs: None)
    with pytest.raises(InputValidationError, match="lacks its DuckDB"):
        backend.resolve_database(resource=tmp_path)
    wrong = tmp_path / "other.duckdb"
    wrong.write_bytes(b"not canonical")
    canonical = tmp_path / "protein_signatures.duckdb"
    canonical.write_bytes(b"canonical")
    with pytest.raises(InputValidationError, match="Unsupported resource file"):
        backend.resolve_database(resource=wrong)


def test_app_launcher_builds_safe_argv_and_runs(
    completed_result: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The launcher should use an argument vector and propagate the child exit code."""

    database = completed_result / "protein_signatures.duckdb"
    command = launcher.build_streamlit_command(database=database, port=8123, address="127.0.0.1")
    assert command[-2:] == ("--resource", str(completed_result))
    assert "8123" in command
    observed: dict[str, object] = {}

    def fake_run(command_value: tuple[str, ...], *, check: bool) -> SimpleNamespace:
        observed.update(command=command_value, check=check)
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(launcher.subprocess, "run", fake_run)
    monkeypatch.setattr(launcher, "check_pdf_export_runtime", lambda: None)
    assert launcher.main(["--resource", str(completed_result), "--port", "8123"]) == 7
    assert observed["check"] is False
    with pytest.raises(ValueError):
        launcher.build_streamlit_command(database=database, port=0, address="localhost")
    with pytest.raises(ValueError):
        launcher.build_streamlit_command(database=database, port=80, address="bad address")
    assert launcher.main(["--resource", str(completed_result), "--port", "0"]) == 2
    monkeypatch.setattr(
        launcher,
        "check_pdf_export_runtime",
        lambda: (_ for _ in ()).throw(PublicationError("missing Chrome")),
    )
    assert launcher.main(["--resource", str(completed_result)]) == 2


def test_app_launcher_preflights_pdf_download_runtime(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The normal launcher should prove PDF export works before exposing the app."""

    observed: list[object] = []
    monkeypatch.setattr(
        launcher,
        "plotly_figure_to_pdf_bytes",
        lambda *, figure: observed.append(figure) or b"%PDF-1.7\n%%EOF\n",
    )
    launcher.check_pdf_export_runtime()
    assert len(observed) == 1


def test_app_launcher_reports_a_missing_plotly_dependency(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The PDF pre-flight should contextualise an unavailable Plotly installation."""

    import builtins

    original_import = builtins.__import__

    def blocked_import(name: str, *args: object, **kwargs: object) -> object:
        """Reject only the lazy Plotly import used by the pre-flight."""

        if name == "plotly.graph_objects":
            raise ImportError("blocked")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked_import)
    with pytest.raises(PublicationError, match="Plotly dependency"):
        launcher.check_pdf_export_runtime()


def test_shell_app_launcher_preserves_argument_boundaries(
    tmp_path: Path,
) -> None:
    """The shell convenience launcher should pass resource, port and address safely."""

    root = Path(__file__).parents[1]
    resource = tmp_path / "result with spaces"
    resource.mkdir()
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_conda = fake_bin / "conda"
    fake_conda.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$@" > "${APP_ARGV_LOG}"\n',
        encoding="utf-8",
    )
    fake_conda.chmod(0o755)
    log = tmp_path / "argv.log"
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["APP_ARGV_LOG"] = str(log)
    result = subprocess.run(
        (
            "bash",
            str(root / "run_protein_signature_app.sh"),
            "--resource",
            str(resource),
            "--conda-environment",
            "app_env",
            "--port",
            "8123",
            "--address",
            "127.0.0.1",
        ),
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert result.returncode == 0, result.stderr
    assert log.read_text(encoding="utf-8").splitlines() == [
        "run",
        "--no-capture-output",
        "--name",
        "app_env",
        "protein-signature-app",
        "--resource",
        str(resource),
        "--port",
        "8123",
        "--address",
        "127.0.0.1",
    ]


def test_cli_commands_and_schema_errors(
    tmp_path: Path, example_dir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """CLI commands should return stable codes and JSON or contextual errors."""

    parser = build_parser()
    assert parser.prog == "protein-signatures"
    assert main(["validate", "--config", str(example_dir / "campaign.yaml")]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "VALID"
    assert main(["describe-profile", "--profile", "e3"]) == 0
    profile_description = json.loads(capsys.readouterr().out)
    assert profile_description["profile_id"] == "e3_default"
    assert profile_description["require_structural_evidence"] is True
    assert profile_description["default_comparison"]["background_label_id"] == (
        "control:prespecified_non_e3_reference"
    )
    f_box = next(
        label
        for label in profile_description["labels"]
        if label["label_id"] == "e3:ubiquitin:crl:crl1_scf:f_box"
    )
    assert f_box["default_analysis"] is True
    assert f_box["default_background_label_id"] == ("control:matched_substrate_receptor_reference")
    result = tmp_path / "cli_result"
    assert (
        main(
            [
                "run-all",
                "--config",
                str(example_dir / "campaign.yaml"),
                "--output-dir",
                str(result),
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)["status"] == "COMPLETE"
    assert main(["verify", "--resource", str(result)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "VALID"
    assert main(["validate", "--config", str(tmp_path / "missing.yaml")]) == 2
    assert "ERROR:" in capsys.readouterr().err
    with pytest.raises(PublicationError):
        schema_for(table_name="not_a_table")


def test_cli_routes_catalogue_initialisation_and_interrupts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Thin CLI routes should bind named arguments and preserve interruption status."""

    starter = tmp_path / "starter"
    config = tmp_path / "campaign.yaml"
    monkeypatch.setattr(cli_module, "prepare_catalogue", lambda **_kwargs: starter)
    assert (
        cli_module.main(
            [
                "prepare-catalogue",
                "--catalogue",
                str(tmp_path / "catalogue.tsv"),
                "--output-dir",
                str(starter),
                "--id-column",
                "id",
                "--sequence-column",
                "sequence",
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)["starter_dir"] == str(starter)
    monkeypatch.setattr(cli_module, "initialise_campaign", lambda **_kwargs: config)
    assert (
        cli_module.main(
            [
                "initialise",
                "--config",
                str(config),
                "--campaign-id",
                "campaign",
                "--sequences-fasta",
                str(tmp_path / "proteins.faa"),
                "--label-assignments",
                str(tmp_path / "labels.tsv"),
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)["config"] == str(config)
    monkeypatch.setattr(
        cli_module,
        "configure_logging",
        lambda **_kwargs: (_ for _ in ()).throw(KeyboardInterrupt()),
    )
    assert cli_module.main(["describe-profile", "--profile", "e3"]) == 130


def _refresh_output_checksum(*, result_dir: Path, relative: str) -> None:
    """Update a test manifest after controlled metadata mutation.

    Args:
        result_dir: Mutable test copy of a published result.
        relative: Output path whose checksum must be refreshed.
    """

    from protein_signatures.checksums import sha256_file

    manifest_path = result_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    changed = result_dir / relative
    for record in manifest["outputs"]:
        if record["relative_path"] == relative:
            record["size_bytes"] = changed.stat().st_size
            record["sha256"] = sha256_file(path=changed)
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")
    marker_path = result_dir / "COMPLETED.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    marker["manifest_sha256"] = sha256_file(path=manifest_path)
    marker_path.write_text(json.dumps(marker), encoding="utf-8")
