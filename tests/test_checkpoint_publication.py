"""Focused tests for bounded analytical checkpoints and large-table reports."""

from __future__ import annotations

import gzip
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest

import protein_signatures.checkpoint as checkpoint_module
import protein_signatures.disk_reporting as disk_reporting_module
import protein_signatures.runtime_resources as runtime_resources_module
from protein_signatures.checkpoint import (
    RecordSequence,
    create_analysis_checkpoint,
    verify_analysis_checkpoint,
)
from protein_signatures.disk_reporting import build_checkpoint_human_reports
from protein_signatures.errors import PublicationError
from protein_signatures.runtime_resources import configure_duckdb_runtime
from protein_signatures.schemas import table_schemas


def test_checkpoint_streams_large_table_and_preserves_record_view_identity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Large tables should use TSV.GZ and bounded Parquet without duplicate models."""

    monkeypatch.setattr(checkpoint_module, "_BATCH_ROWS", 2)
    monkeypatch.setattr(checkpoint_module, "_PARQUET_PART_ROWS", 2)
    monkeypatch.setattr(checkpoint_module, "_LARGE_TSV_ROWS", 2)
    values = tuple(_Record(index=index) for index in range(3))
    view = RecordSequence(values=values)
    assert len(view) == 3
    assert view._values is values  # noqa: SLF001 - explicit memory-regression contract
    tables = _empty_tables()
    tables["features"] = view
    authority = tmp_path / "authority.tsv"
    authority.write_text("authority\n", encoding="utf-8")

    checkpoint = create_analysis_checkpoint(
        checkpoint_dir=tmp_path / "checkpoint",
        run_identity_sha256="a" * 64,
        tables=tables,
        state={"status": "ANALYSIS_COMPLETE"},
        input_paths=(authority,),
    )

    features = checkpoint.table(name="features")
    assert features.tsv_format == "TSV.GZ"
    assert features.tsv_path.name == "features.tsv.gz"
    assert len(features.parquet_parts) == 2
    assert sum(pq.ParquetFile(path).metadata.num_rows for path in features.parquet_parts) == 3
    assert features.tsv_path.read_bytes()[4:8] == b"\x00\x00\x00\x00"
    with gzip.open(features.tsv_path, mode="rt", encoding="utf-8") as handle:
        assert len(handle.readlines()) == 4
    assert (
        verify_analysis_checkpoint(
            checkpoint_dir=checkpoint.root,
            run_identity_sha256="a" * 64,
            verify_inputs=True,
        ).counts["features"]
        == 3
    )


def test_checkpoint_rejects_changed_inputs_and_table_files(tmp_path: Path) -> None:
    """Resume must fail closed when an authority or checkpoint file changes."""

    authority = tmp_path / "authority.tsv"
    authority.write_text("authority\n", encoding="utf-8")
    checkpoint = create_analysis_checkpoint(
        checkpoint_dir=tmp_path / "checkpoint",
        run_identity_sha256="b" * 64,
        tables=_empty_tables(),
        state={"status": "ANALYSIS_COMPLETE"},
        input_paths=(authority,),
    )
    authority.write_text("changed\n", encoding="utf-8")
    with pytest.raises(PublicationError, match="input .*changed"):
        verify_analysis_checkpoint(
            checkpoint_dir=checkpoint.root,
            run_identity_sha256="b" * 64,
            verify_inputs=True,
        )
    authority.write_text("authority\n", encoding="utf-8")
    proteins = checkpoint.table(name="proteins").parquet_parts[0]
    proteins.write_bytes(proteins.read_bytes() + b"tamper")
    with pytest.raises(PublicationError, match="size differs"):
        verify_analysis_checkpoint(
            checkpoint_dir=checkpoint.root,
            run_identity_sha256="b" * 64,
            verify_inputs=True,
        )


def test_large_feature_report_uses_summary_excel_and_self_contained_html(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A large feature authority should never become a monolithic workbook."""

    monkeypatch.setattr(checkpoint_module, "_LARGE_TSV_ROWS", 2)
    monkeypatch.setattr(disk_reporting_module, "_FULL_EXCEL_MAX_ROWS", 2)
    tables = _empty_tables()
    tables["features"] = tuple(_feature_record(index=index) for index in range(3))
    authority = tmp_path / "authority.tsv"
    authority.write_text("authority\n", encoding="utf-8")
    checkpoint = create_analysis_checkpoint(
        checkpoint_dir=tmp_path / "checkpoint",
        run_identity_sha256="c" * 64,
        tables=tables,
        state={"status": "ANALYSIS_COMPLETE"},
        input_paths=(authority,),
    )

    reports = build_checkpoint_human_reports(
        checkpoint=checkpoint,
        cache_dir=tmp_path / "reports",
        fdr_threshold=0.05,
        summary_context={"campaign_id": "large_demo", "profile_display_name": "Demo"},
    )
    assets = dict(reports.assets)

    summary = "analysis/03_sequence_and_domains/tables/features_summary.xlsx"
    assert summary in assets
    assert assets[summary].read_bytes().startswith(b"PK")
    assert "analysis/03_sequence_and_domains/tables/features.xlsx" not in assets
    html_path = assets["analysis/00_run_information/results_summary.html"]
    html_text = html_path.read_text(encoding="utf-8")
    assert "large_demo" in html_text
    assert "Feature memberships" in html_text
    assert "119" not in html_text
    assert any(
        row["asset_kind"] == "CANONICAL_TABLE"
        and row["content_id"] == "features"
        and row["file_format"] == "TSV.GZ"
        for row in reports.inventory
    )


def test_duckdb_runtime_uses_bounded_memory_and_validated_spill(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """DuckDB should retain allocation headroom and require a real spill directory."""

    monkeypatch.setattr(runtime_resources_module, "_available_memory_mb", lambda: 2_048)
    with duckdb.connect(":memory:") as connection:
        assert (
            configure_duckdb_runtime(
                connection=connection,
                temporary_directory=tmp_path,
            )
            == 1_433
        )
    with duckdb.connect(":memory:") as connection:
        with pytest.raises(OSError, match="temporary directory is missing"):
            configure_duckdb_runtime(
                connection=connection,
                temporary_directory=tmp_path / "missing",
            )


class _Record:
    """Small model-like value used to exercise lazy record conversion."""

    def __init__(self, *, index: int) -> None:
        """Store the feature index."""

        self.index = index

    def to_record(self) -> dict[str, object]:
        """Return one schema-complete feature record."""

        return _feature_record(index=self.index)


def _empty_tables() -> dict[str, tuple[dict[str, object], ...] | RecordSequence]:
    """Return every canonical relation with no rows."""

    return {name: () for name in table_schemas()}


def _feature_record(*, index: int) -> dict[str, object]:
    """Return one valid canonical amino-acid feature record."""

    return {
        "protein_id": f"p{index}",
        "feature_type": "AMINO_ACID_KMER",
        "feature_id": f"k3:A{index}A",
        "feature_name": f"k3 A{index}A",
        "start": 1,
        "end": 3,
        "evidence_status": "PRESENT",
        "evidence_source": "test",
        "evidence_reference": "test",
        "derivation_scope": "DISCOVERY_DERIVED",
        "feature_definition_sha256": f"{index + 1:064x}",
        "derivation_cohort_sha256": "f" * 64,
    }
