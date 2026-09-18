"""Bounded-memory reports built from verified analytical checkpoints."""

from __future__ import annotations

import html
import logging
import os
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

from .checkpoint import AnalysisCheckpoint, CheckpointTable
from .errors import InputValidationError, PublicationError
from .exports import dataframe_to_tsv_bytes, dataframe_to_xlsx_bytes
from .io_utils import write_text_atomic
from .reporting import (
    _FIGURE_FORMATS,
    _FINAL_TABLES,
    _TABLE_SECTIONS,
    HumanReportResult,
    _build_static_figures,
    _final_shap_relative_path,
    _inventory_row,
    _validate_fdr_threshold,
    _write_binary_atomic,
)
from .runtime_resources import configure_duckdb_runtime
from .schemas import table_schemas

LOGGER = logging.getLogger(__name__)

_FULL_EXCEL_MAX_ROWS = 500_000
_FIGURE_TABLES = frozenset(
    {
        "associations",
        "class_labelling_summary",
        "comparisons",
        "control_matching_audit",
        "domain_assessments",
        "label_assignments",
        "ml_feature_importance",
        "ml_models",
        "ml_predictions",
        "orthofinder_group_context",
        "partitions",
        "signatures",
        "structures",
    }
)

_REPORT_GUIDE = """# Human-readable analysis reports

The numbered directories contain bounded, human-readable summaries and formatted Excel
workbooks. Complete canonical data remain under `tables/` as typed Parquet and either TSV
or compressed TSV.GZ. The DuckDB database is the query and application authority.

Very large canonical tables are deliberately not converted into monolithic Excel files.
Their summary workbook explains the complete row count and available formats. Use the
application's filtered feature exporter to create practical TSV or Excel subsets.

Open `results_summary.html` for a self-contained overview of this completed campaign.
Verify `COMPLETED.json` and `manifest.json` before interpreting copied results.
"""


class _DuckDbFrameMapping(Mapping[str, pd.DataFrame]):
    """Load one bounded figure-support table only when it is requested."""

    def __init__(
        self,
        *,
        connection: duckdb.DuckDBPyConnection,
        table_names: Sequence[str],
    ) -> None:
        """Store a trusted connection and canonical allow-list."""

        self._connection = connection
        self._table_names = tuple(table_names)

    def __getitem__(self, key: str) -> pd.DataFrame:
        """Read one canonical table into one short-lived data frame."""

        if key not in self._table_names:
            raise KeyError(key)
        LOGGER.info("Loading bounded figure-support table table=%s", key)
        return self._connection.execute(
            f'SELECT * FROM "{key}"'  # noqa: S608 - canonical allow-list above
        ).fetchdf()

    def __iter__(self) -> Iterator[str]:
        """Iterate the canonical allow-list."""

        return iter(self._table_names)

    def __len__(self) -> int:
        """Return the allow-list size."""

        return len(self._table_names)


def build_checkpoint_human_reports(
    *,
    checkpoint: AnalysisCheckpoint,
    cache_dir: Path,
    fdr_threshold: float,
    summary_context: Mapping[str, Any],
    existing_plot_assets: Sequence[tuple[str, Path]] = (),
) -> HumanReportResult:
    """Build bounded reports from checkpoint files rather than in-memory tables.

    Args:
        checkpoint: Verified complete analytical checkpoint.
        cache_dir: Run-specific report cache.
        fdr_threshold: Configured false-discovery-rate threshold.
        summary_context: Campaign and profile metadata for the HTML overview.
        existing_plot_assets: Native SHAP plots already present on disk.

    Returns:
        Portable report assets and inventory.

    Raises:
        InputValidationError: If the checkpoint table inventory is incomplete.
        PublicationError: If a report cannot be written.
    """

    threshold = _validate_fdr_threshold(fdr_threshold=fdr_threshold)
    tables = {table.name: table for table in checkpoint.tables}
    if set(tables) != set(table_schemas()) or set(tables) != set(_TABLE_SECTIONS):
        raise InputValidationError("Checkpoint/report table routing inventory differs.")
    root = Path(cache_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    LOGGER.info(
        "Building bounded human reports for %d checkpoint tables in %s",
        len(tables),
        root,
    )
    assets: dict[str, Path] = {}
    inventory: list[dict[str, Any]] = []
    scratch_parent = Path(os.environ.get("TMPDIR", "")).expanduser()
    if not scratch_parent.is_dir():
        scratch_parent = root
    with tempfile.TemporaryDirectory(
        prefix="protein_signature_reports_",
        dir=scratch_parent,
    ) as temporary_directory:
        with duckdb.connect(":memory:") as connection:
            configure_duckdb_runtime(
                connection=connection,
                temporary_directory=Path(temporary_directory),
            )
            _attach_checkpoint_views(connection=connection, tables=tables)
            feature_counts = connection.execute(
                "SELECT feature_type, count(*) AS membership_row_count, "
                "count(DISTINCT feature_id) AS feature_count, "
                "count(DISTINCT protein_id) AS protein_count "
                "FROM features GROUP BY feature_type ORDER BY membership_row_count DESC, "
                "feature_type"
            ).fetchdf()
            for table_index, table_name in enumerate(sorted(tables), start=1):
                table = tables[table_name]
                LOGGER.info(
                    "Building human table report %d/%d table=%s rows=%d",
                    table_index,
                    len(tables),
                    table_name,
                    table.row_count,
                )
                _register_checkpoint_table(
                    connection=connection,
                    table=table,
                    section=_TABLE_SECTIONS[table_name],
                    root=root,
                    assets=assets,
                    inventory=inventory,
                    feature_counts=feature_counts,
                )
                if table_name in _FINAL_TABLES:
                    _register_final_table_copy(
                        table=table,
                        root=root,
                        assets=assets,
                        inventory=inventory,
                    )
            frames = _DuckDbFrameMapping(
                connection=connection,
                table_names=sorted(_FIGURE_TABLES),
            )
            figure_count = _build_static_figures(
                frames=frames,
                feature_counts=feature_counts,
                fdr_threshold=threshold,
                root=root,
                assets=assets,
                inventory=inventory,
            )
            del frames
            summary_path = root / "analysis/00_run_information/results_summary.html"
            summary_html = _build_html_summary(
                connection=connection,
                checkpoint=checkpoint,
                summary_context=summary_context,
                feature_counts=feature_counts,
                fdr_threshold=threshold,
            )
            write_text_atomic(path=summary_path, text=summary_html)
    summary_relative = "analysis/00_run_information/results_summary.html"
    assets[summary_relative] = summary_path
    inventory.append(
        _inventory_row(
            section="00_run_information",
            asset_kind="DOCUMENTATION",
            content_id="results_summary",
            file_format="HTML",
            relative_path=summary_relative,
            row_count=0,
            description="Self-contained campaign results and interpretation overview.",
        )
    )
    _register_existing_plots(
        existing_plot_assets=existing_plot_assets,
        assets=assets,
        inventory=inventory,
    )
    guide_relative = "analysis/00_run_information/README.md"
    guide_path = root / guide_relative
    write_text_atomic(path=guide_path, text=_REPORT_GUIDE)
    assets[guide_relative] = guide_path
    inventory.append(
        _inventory_row(
            section="00_run_information",
            asset_kind="DOCUMENTATION",
            content_id="report_guide",
            file_format="MARKDOWN",
            relative_path=guide_relative,
            row_count=0,
            description="Reading order, large-table policy and result authorities.",
        )
    )
    ordered_inventory = tuple(_ordered_inventory(rows=inventory))
    inventory_frame = pd.DataFrame.from_records(ordered_inventory)
    inventory_root = "analysis/00_run_information/report_inventory"
    inventory_tsv = root / f"{inventory_root}.tsv"
    inventory_xlsx = root / f"{inventory_root}.xlsx"
    write_text_atomic(
        path=inventory_tsv,
        text=dataframe_to_tsv_bytes(frame=inventory_frame).decode("utf-8"),
    )
    _write_binary_atomic(
        path=inventory_xlsx,
        payload=dataframe_to_xlsx_bytes(
            frame=inventory_frame,
            title="Human-readable report inventory",
        ),
    )
    assets[f"{inventory_root}.tsv"] = inventory_tsv
    assets[f"{inventory_root}.xlsx"] = inventory_xlsx
    result = HumanReportResult(
        assets=tuple(sorted(assets.items())),
        inventory=ordered_inventory,
        figure_count=figure_count
        + len(
            {
                Path(relative_path).with_suffix("").as_posix()
                for relative_path, _source in existing_plot_assets
            }
        ),
        workbook_count=sum(relative.endswith(".xlsx") for relative in assets),
    )
    LOGGER.info(
        "Built bounded human reports files=%d workbooks=%d figures=%d html=1",
        len(result.assets),
        result.workbook_count,
        result.figure_count,
    )
    return result


def _attach_checkpoint_views(
    *, connection: duckdb.DuckDBPyConnection, tables: Mapping[str, CheckpointTable]
) -> None:
    """Expose trusted checkpoint Parquet tables to an internal DuckDB connection."""

    connection.execute("SET preserve_insertion_order = false")
    for table_name, table in sorted(tables.items()):
        escaped_path = str(table.parquet_path / "part-*.parquet").replace("'", "''")
        connection.execute(
            f"CREATE VIEW \"{table_name}\" AS SELECT * FROM read_parquet('{escaped_path}')"
        )


def _register_checkpoint_table(
    *,
    connection: duckdb.DuckDBPyConnection,
    table: CheckpointTable,
    section: str,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
    feature_counts: pd.DataFrame,
) -> None:
    """Register one full bounded table or a large-table summary workbook."""

    if table.row_count <= _FULL_EXCEL_MAX_ROWS:
        frame = connection.execute(f'SELECT * FROM "{table.name}"').fetchdf()
        _register_full_table_assets(
            frame=frame,
            table=table,
            section=section,
            root=root,
            assets=assets,
            inventory=inventory,
        )
        del frame
        return
    source_relative = f"tables/{table.tsv_path.name}"
    inventory.append(
        _inventory_row(
            section=section,
            asset_kind="CANONICAL_TABLE",
            content_id=table.name,
            file_format=table.tsv_format,
            relative_path=source_relative,
            row_count=table.row_count,
            description=(
                "Complete canonical large table; query through DuckDB or export a filtered "
                "subset in the application."
            ),
        )
    )
    if table.name == "features":
        summary = feature_counts.copy(deep=True)
    else:
        summary = pd.DataFrame.from_records(
            (
                {
                    "table_name": table.name,
                    "query_relation": table.name,
                },
            )
        )
    summary.insert(0, "complete_table_row_count", table.row_count)
    summary.insert(1, "canonical_tsv", source_relative)
    summary.insert(2, "canonical_parquet", f"tables/{table.parquet_path.name}")
    relative_root = Path("analysis") / section / "tables" / f"{table.name}_summary"
    _register_summary_assets(
        frame=summary,
        table=table,
        relative_root=relative_root,
        root=root,
        assets=assets,
        inventory=inventory,
    )


def _register_full_table_assets(
    *,
    frame: pd.DataFrame,
    table: CheckpointTable,
    section: str,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
) -> None:
    """Register checkpoint TSV and one formatted workbook for a manageable table."""

    relative_root = Path("analysis") / section / "tables" / table.name
    tsv_relative = relative_root.with_suffix(
        ".tsv.gz" if table.tsv_format == "TSV.GZ" else ".tsv"
    ).as_posix()
    xlsx_relative = relative_root.with_suffix(".xlsx").as_posix()
    if tsv_relative in assets or xlsx_relative in assets:
        raise PublicationError(f"Human table report path is duplicated: {relative_root}")
    xlsx_path = root / xlsx_relative
    _write_binary_atomic(
        path=xlsx_path,
        payload=dataframe_to_xlsx_bytes(
            frame=frame,
            title=table.name.replace("_", " ").title(),
        ),
    )
    assets[tsv_relative] = table.tsv_path
    assets[xlsx_relative] = xlsx_path
    for file_format, relative_path in (
        (table.tsv_format, tsv_relative),
        ("XLSX", xlsx_relative),
    ):
        inventory.append(
            _inventory_row(
                section=section,
                asset_kind="TABLE",
                content_id=table.name,
                file_format=file_format,
                relative_path=relative_path,
                row_count=table.row_count,
                description="Complete human-facing copy of the canonical result table.",
            )
        )


def _register_summary_assets(
    *,
    frame: pd.DataFrame,
    table: CheckpointTable,
    relative_root: Path,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
) -> None:
    """Write compact TSV and Excel navigation assets for a large table."""

    tsv_relative = relative_root.with_suffix(".tsv").as_posix()
    xlsx_relative = relative_root.with_suffix(".xlsx").as_posix()
    tsv_path = root / tsv_relative
    xlsx_path = root / xlsx_relative
    write_text_atomic(
        path=tsv_path,
        text=dataframe_to_tsv_bytes(frame=frame).decode("utf-8"),
    )
    _write_binary_atomic(
        path=xlsx_path,
        payload=dataframe_to_xlsx_bytes(
            frame=frame,
            title=f"{table.name.replace('_', ' ').title()} Summary",
        ),
    )
    assets[tsv_relative] = tsv_path
    assets[xlsx_relative] = xlsx_path
    for file_format, relative_path in (("TSV", tsv_relative), ("XLSX", xlsx_relative)):
        inventory.append(
            _inventory_row(
                section=str(relative_root.parts[1]),
                asset_kind="TABLE_INDEX",
                content_id=table.name,
                file_format=file_format,
                relative_path=relative_path,
                row_count=len(frame),
                description=(
                    f"Compact summary for {table.row_count:,} canonical rows; this is not "
                    "a replacement for the complete Parquet/TSV.GZ table."
                ),
            )
        )


def _register_final_table_copy(
    *,
    table: CheckpointTable,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
) -> None:
    """Register a decision-facing copy without regenerating workbook bytes."""

    source_section = _TABLE_SECTIONS[table.name]
    source_root = Path("analysis") / source_section / "tables"
    if table.row_count <= _FULL_EXCEL_MAX_ROWS:
        source_paths = (
            source_root / f"{table.name}{'.tsv.gz' if table.tsv_format == 'TSV.GZ' else '.tsv'}",
            source_root / f"{table.name}.xlsx",
        )
    else:
        source_paths = (
            source_root / f"{table.name}_summary.tsv",
            source_root / f"{table.name}_summary.xlsx",
        )
    for source_relative in source_paths:
        source_key = source_relative.as_posix()
        source = assets[source_key]
        final_relative = (
            Path("analysis") / "99_final_results" / "tables" / source_relative.name
        ).as_posix()
        if final_relative in assets:
            raise PublicationError(f"Final table report path is duplicated: {final_relative!r}")
        assets[final_relative] = source
        inventory.append(
            _inventory_row(
                section="99_final_results",
                asset_kind=("TABLE" if table.row_count <= _FULL_EXCEL_MAX_ROWS else "TABLE_INDEX"),
                content_id=table.name,
                file_format=("XLSX" if final_relative.endswith(".xlsx") else "TSV"),
                relative_path=final_relative,
                row_count=(table.row_count if table.row_count <= _FULL_EXCEL_MAX_ROWS else 0),
                description="Decision-facing copy of the corresponding numbered report.",
            )
        )


def _register_existing_plots(
    *,
    existing_plot_assets: Sequence[tuple[str, Path]],
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
) -> None:
    """Register native and decision-facing SHAP plots."""

    for relative_path, source_path in existing_plot_assets:
        if Path(relative_path).suffix.lstrip(".").casefold() not in _FIGURE_FORMATS:
            raise PublicationError(f"Unsupported existing plot format: {relative_path!r}")
        if relative_path in assets:
            raise PublicationError(f"Human-report asset path is duplicated: {relative_path!r}")
        assets[relative_path] = source_path
        inventory.append(
            _inventory_row(
                section="06_explainable_models",
                asset_kind="FIGURE",
                content_id=Path(relative_path).stem,
                file_format=Path(relative_path).suffix.lstrip(".").upper(),
                relative_path=relative_path,
                row_count=0,
                description="Native SHAP model explanation.",
            )
        )
        final_relative = _final_shap_relative_path(relative_path=relative_path)
        assets[final_relative] = source_path
        inventory.append(
            _inventory_row(
                section="99_final_results",
                asset_kind="FIGURE",
                content_id=f"shap_{Path(relative_path).stem}",
                file_format=Path(relative_path).suffix.lstrip(".").upper(),
                relative_path=final_relative,
                row_count=0,
                description="Decision-facing copy of a native SHAP explanation.",
            )
        )


def _build_html_summary(
    *,
    connection: duckdb.DuckDBPyConnection,
    checkpoint: AnalysisCheckpoint,
    summary_context: Mapping[str, Any],
    feature_counts: pd.DataFrame,
    fdr_threshold: float,
) -> str:
    """Return a portable HTML overview containing only bounded aggregates."""

    campaign_id = html.escape(str(summary_context.get("campaign_id", "unknown")))
    profile_name = html.escape(str(summary_context.get("profile_display_name", "unknown")))
    evidence = connection.execute(
        "SELECT evidence_class, count(*) AS signature_count FROM signatures "
        "GROUP BY evidence_class ORDER BY signature_count DESC, evidence_class"
    ).fetchdf()
    comparison_status = connection.execute(
        "SELECT comparison_id, count(*) AS association_rows, "
        "sum(CASE WHEN status = 'COMPLETE' THEN 1 ELSE 0 END) AS complete_rows, "
        "min(study_q_value) AS minimum_study_q_value "
        "FROM associations GROUP BY comparison_id ORDER BY comparison_id"
    ).fetchdf()
    top_signatures = connection.execute(
        "SELECT comparison_id, feature_type, feature_id, feature_name, "
        "CASE WHEN discovery_prevalence_difference > 0 THEN 'ENRICHED_IN_TARGET' "
        "WHEN discovery_prevalence_difference < 0 THEN 'DEPLETED_IN_TARGET' "
        "ELSE 'NO_DIFFERENCE' END AS direction, discovery_q_value, "
        "validation_q_value, discovery_study_q_value, validation_study_q_value, "
        "evidence_class FROM signatures ORDER BY validation_study_q_value NULLS LAST, "
        "discovery_study_q_value NULLS LAST, discovery_q_value NULLS LAST, "
        "comparison_id, feature_type, feature_id LIMIT 100"
    ).fetchdf()
    count_rows = pd.DataFrame.from_records(
        (
            {"canonical_table": name, "row_count": count}
            for name, count in sorted(checkpoint.counts.items())
        )
    )
    css = """
body{font-family:Arial,Helvetica,sans-serif;color:#17324d;margin:0;background:#f3f7fa}
main{max-width:1280px;margin:0 auto;padding:28px}h1,h2{color:#17324d}
.hero{background:linear-gradient(135deg,#17324d,#2f6f95);color:white;padding:28px;
border-radius:14px}.hero h1{color:white;margin-top:0}.grid{display:grid;
grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:18px 0}
.card{background:white;border:1px solid #d7e2ea;border-radius:10px;padding:16px}
.metric{font-size:1.8rem;font-weight:700;color:#1f4e78}.note{background:#fff8dc;
border-left:5px solid #d4a017;padding:12px;margin:16px 0}table{border-collapse:collapse;
width:100%;background:white;font-size:.88rem;margin-bottom:24px}th{background:#1f4e78;
color:white;text-align:left}th,td{border:1px solid #d9e2f3;padding:7px;vertical-align:top}
tr:nth-child(even){background:#eef5fa}a{color:#145f8c}code{background:#edf2f5;padding:2px 4px}
"""
    key_counts = {
        "Proteins": checkpoint.counts.get("proteins", 0),
        "Feature memberships": checkpoint.counts.get("features", 0),
        "Associations": checkpoint.counts.get("associations", 0),
        "Signatures": checkpoint.counts.get("signatures", 0),
        "Comparisons": checkpoint.counts.get("comparisons", 0),
    }
    cards = "".join(
        f'<div class="card"><div>{html.escape(label)}</div>'
        f'<div class="metric">{value:,}</div></div>'
        for label, value in key_counts.items()
    )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{campaign_id} protein-signature results</title><style>{css}</style></head>"
        "<body><main>"
        f'<section class="hero"><h1>Protein-signature analysis</h1><p><strong>Campaign:</strong> '
        f"{campaign_id}<br><strong>Profile:</strong> {profile_name}<br>"
        f"<strong>FDR threshold:</strong> {fdr_threshold:g}</p></section>"
        f'<section class="grid">{cards}</section>'
        '<div class="note"><strong>Interpretation:</strong> signatures are prioritisation '
        "evidence, not proof of biochemical activity. Discovery, held-out validation, "
        "assessment universes and independence blocks remain explicit.</div>"
        "<h2>Feature evidence coverage</h2>"
        f"{feature_counts.to_html(index=False, border=0, escape=True)}"
        "<h2>Signature evidence classes</h2>"
        f"{evidence.to_html(index=False, border=0, escape=True)}"
        "<h2>Comparison completion</h2>"
        f"{comparison_status.to_html(index=False, border=0, escape=True)}"
        "<h2>Top prioritised signatures</h2><p>The first 100 rows are ordered by study-wide "
        "and discovery false-discovery rates.</p>"
        f"{top_signatures.to_html(index=False, border=0, escape=True)}"
        "<h2>Canonical table inventory</h2>"
        f"{count_rows.to_html(index=False, border=0, escape=True)}"
        "<h2>Data access</h2><p>Use <code>protein_signatures.duckdb</code> for queries, "
        "the application for filtered TSV/Excel exports, and <code>tables/</code> for complete "
        "Parquet and TSV/TSV.GZ authorities.</p>"
        "</main></body></html>"
    )


def _ordered_inventory(*, rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Return deterministic report inventory records."""

    return sorted(
        (dict(row) for row in rows),
        key=lambda row: (
            str(row["section"]),
            str(row["asset_kind"]),
            str(row["content_id"]),
            str(row["file_format"]),
            str(row["relative_path"]),
        ),
    )
