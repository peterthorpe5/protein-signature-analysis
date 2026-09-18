"""Numbered, human-readable table and figure reports for completed campaigns."""

from __future__ import annotations

import hashlib
import logging
import math
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from .errors import InputValidationError, PublicationError
from .exports import dataframe_to_tsv_bytes, dataframe_to_xlsx_bytes
from .io_utils import write_text_atomic
from .schemas import schema_for

matplotlib.use("Agg", force=True)

LOGGER = logging.getLogger(__name__)

_TABLE_SECTIONS = {
    "proteins": "01_proteins_and_curation",
    "profile_labels": "01_proteins_and_curation",
    "label_assignments": "01_proteins_and_curation",
    "label_evidence_audit": "01_proteins_and_curation",
    "control_matching_audit": "01_proteins_and_curation",
    "label_definition_features": "01_proteins_and_curation",
    "class_labelling_summary": "01_proteins_and_curation",
    "unresolved_assignments": "01_proteins_and_curation",
    "label_memberships": "01_proteins_and_curation",
    "comparisons": "01_proteins_and_curation",
    "redundancy_clusters": "02_homology_and_partitions",
    "orthofinder_memberships": "02_homology_and_partitions",
    "orthofinder_group_context": "02_homology_and_partitions",
    "partitions": "02_homology_and_partitions",
    "features": "03_sequence_and_domains",
    "feature_assessments": "03_sequence_and_domains",
    "domain_hits": "03_sequence_and_domains",
    "domain_assessments": "03_sequence_and_domains",
    "domain_sequences": "03_sequence_and_domains",
    "structures": "04_structures_and_folds",
    "alphafold_acquisitions": "04_structures_and_folds",
    "structure_comparisons": "04_structures_and_folds",
    "structure_clusters": "04_structures_and_folds",
    "imported_structural_group_summaries": "04_structures_and_folds",
    "associations": "05_association_statistics",
    "signatures": "05_association_statistics",
    "ml_models": "06_explainable_models",
    "ml_feature_importance": "06_explainable_models",
    "ml_predictions": "06_explainable_models",
    "ml_explanations": "06_explainable_models",
    "ml_plot_inventory": "06_explainable_models",
}
_FINAL_TABLES = (
    "comparisons",
    "class_labelling_summary",
    "label_evidence_audit",
    "control_matching_audit",
    "label_definition_features",
    "unresolved_assignments",
    "signatures",
    "associations",
    "structures",
    "structure_clusters",
    "ml_models",
    "ml_feature_importance",
    "ml_predictions",
    "ml_explanations",
    "ml_plot_inventory",
)
_FIGURE_FORMATS = ("png", "svg", "pdf")
_REPORT_GUIDE = """# Human-readable analysis reports

Read the numbered directories in order. Every canonical table is copied here as TSV and a
formatted, filterable Excel workbook. The `tables/` directory at result root remains the
machine authority in TSV and Parquet, and `protein_signatures.duckdb` remains the app/query
authority.

- `01_proteins_and_curation`: FASTA inventory, labels, automated evidence decisions,
  matched-control diagnostics, circularity exclusions, abstentions and comparisons.
- `02_homology_and_partitions`: exact/near redundancy, OrthoFinder and frozen splits.
- `03_sequence_and_domains`: unified amino-acid/Pfam evidence plus the complete imported
  feature-assessment and derivation-provenance ledger.
- `04_structures_and_folds`: models, AlphaFold outcomes, alignments and structural clusters.
- `05_association_statistics`: complete tests and discovery/validation signature summaries.
- `06_explainable_models`: model status, held-out predictions, SHAP values and figures.
- `99_final_results`: decision-facing copies of principal tables and figures.

`report_inventory.tsv` lists the generated human-facing files. Missing scientific evidence
is represented by table status values; absence of a figure is not proof of biological
absence. Verify `COMPLETED.json` and `manifest.json` before interpreting a copied result.
"""


@dataclass(frozen=True)
class HumanReportResult:
    """Portable report assets and their human-facing inventory."""

    assets: tuple[tuple[str, Path], ...]
    inventory: tuple[dict[str, Any], ...]
    figure_count: int
    workbook_count: int


def _validate_fdr_threshold(*, fdr_threshold: float) -> float:
    """Return a finite false-discovery-rate threshold in the unit interval.

    Args:
        fdr_threshold: Candidate false-discovery-rate significance threshold.

    Returns:
        Validated floating-point threshold.

    Raises:
        InputValidationError: If the threshold is not numeric or lies outside
            the interval ``0 < threshold <= 1``.
    """

    if isinstance(fdr_threshold, bool):
        raise InputValidationError("FDR threshold must be numeric and satisfy 0 < threshold <= 1.")
    try:
        threshold = float(fdr_threshold)
    except (TypeError, ValueError) as exc:
        raise InputValidationError(
            "FDR threshold must be numeric and satisfy 0 < threshold <= 1."
        ) from exc
    if not math.isfinite(threshold) or threshold <= 0.0 or threshold > 1.0:
        raise InputValidationError("FDR threshold must be numeric and satisfy 0 < threshold <= 1.")
    return threshold


def build_human_reports(
    *,
    tables: Mapping[str, Sequence[Mapping[str, Any]]],
    cache_dir: Path,
    fdr_threshold: float,
    existing_plot_assets: Sequence[tuple[str, Path]] = (),
) -> HumanReportResult:
    """Build numbered TSV/Excel reports and static figures before publication.

    Args:
        tables: Complete canonical table records, before report publication.
        cache_dir: Run-specific local report cache.
        fdr_threshold: Configured false-discovery-rate significance threshold.
        existing_plot_assets: Already-rendered SHAP result paths and sources.

    Returns:
        Result-relative asset mappings, inventory rows and report counts.

    Raises:
        InputValidationError: If the FDR threshold is invalid or a canonical table
            lacks a report section.
        PublicationError: If a table or figure cannot be generated atomically.
    """

    validated_fdr_threshold = _validate_fdr_threshold(fdr_threshold=fdr_threshold)
    unknown = sorted(set(tables) - set(_TABLE_SECTIONS))
    missing = sorted(set(_TABLE_SECTIONS) - set(tables))
    if unknown or missing:
        raise InputValidationError(
            f"Human-report table routing mismatch; unknown={unknown}, missing={missing}."
        )
    root = Path(cache_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    LOGGER.info(
        "Building numbered human reports for %d canonical tables in %s",
        len(tables),
        root,
    )
    assets: dict[str, Path] = {}
    inventory: list[dict[str, Any]] = []
    frames: dict[str, pd.DataFrame] = {}
    for table_name in sorted(tables):
        frame = _table_frame(table_name=table_name, records=tables[table_name])
        frames[table_name] = frame
        _register_table_exports(
            frame=frame,
            table_name=table_name,
            section=_TABLE_SECTIONS[table_name],
            root=root,
            assets=assets,
            inventory=inventory,
        )
    for table_name in _FINAL_TABLES:
        _register_table_exports(
            frame=frames[table_name],
            table_name=table_name,
            section="99_final_results",
            root=root,
            assets=assets,
            inventory=inventory,
        )
    figure_count = _build_static_figures(
        frames=frames,
        fdr_threshold=validated_fdr_threshold,
        root=root,
        assets=assets,
        inventory=inventory,
    )
    for relative_path, source_path in existing_plot_assets:
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
        if final_relative in assets:
            raise PublicationError(f"Final SHAP asset path is duplicated: {final_relative!r}")
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
    guide_relative = "analysis/00_run_information/README.md"
    guide_source = root / guide_relative
    write_text_atomic(path=guide_source, text=_REPORT_GUIDE)
    assets[guide_relative] = guide_source
    inventory.append(
        _inventory_row(
            section="00_run_information",
            asset_kind="DOCUMENTATION",
            content_id="report_guide",
            file_format="MARKDOWN",
            relative_path=guide_relative,
            row_count=0,
            description="Reading order and authority guide.",
        )
    )
    ordered_inventory = tuple(
        sorted(
            inventory,
            key=lambda row: (
                str(row["section"]),
                str(row["asset_kind"]),
                str(row["content_id"]),
                str(row["file_format"]),
                str(row["relative_path"]),
            ),
        )
    )
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
                for relative_path, _source_path in existing_plot_assets
            }
        ),
        workbook_count=sum(path.endswith(".xlsx") for path in assets),
    )
    LOGGER.info(
        "Built %d human-report files, %d Excel workbooks and %d logical figures",
        len(result.assets),
        result.workbook_count,
        result.figure_count,
    )
    return result


def _table_frame(*, table_name: str, records: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    """Create a schema-ordered frame for one canonical table.

    Args:
        table_name: Canonical table identifier.
        records: Table records.

    Returns:
        Data frame whose columns match the Arrow schema exactly.
    """

    columns = schema_for(table_name=table_name).names
    return pd.DataFrame.from_records((dict(row) for row in records), columns=columns)


def _register_table_exports(
    *,
    frame: pd.DataFrame,
    table_name: str,
    section: str,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
) -> None:
    """Write and register one human-facing TSV/XLSX table pair.

    Args:
        frame: Schema-ordered table.
        table_name: Canonical content identifier.
        section: Numbered report section.
        root: Report cache root.
        assets: Mutable result-asset mapping.
        inventory: Mutable report inventory.

    Raises:
        PublicationError: If a destination would be repeated.
    """

    relative_root = Path("analysis") / section / "tables"
    for file_format in ("TSV", "XLSX"):
        suffix = file_format.casefold()
        relative = (relative_root / f"{table_name}.{suffix}").as_posix()
        if relative in assets:
            raise PublicationError(f"Human table report path is duplicated: {relative!r}")
        source = root / relative
        if file_format == "TSV":
            write_text_atomic(
                path=source,
                text=dataframe_to_tsv_bytes(frame=frame).decode("utf-8"),
            )
        else:
            _write_binary_atomic(
                path=source,
                payload=dataframe_to_xlsx_bytes(
                    frame=frame,
                    title=table_name.replace("_", " ").title(),
                ),
            )
        assets[relative] = source
        LOGGER.debug("Built %s table report %s", file_format, relative)
        inventory.append(
            _inventory_row(
                section=section,
                asset_kind="TABLE",
                content_id=table_name,
                file_format=file_format,
                relative_path=relative,
                row_count=len(frame),
                description="Human-facing copy of the canonical result table.",
            )
        )


def _build_static_figures(
    *,
    frames: Mapping[str, pd.DataFrame],
    fdr_threshold: float,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
    feature_counts: pd.DataFrame | None = None,
) -> int:
    """Build global and comparison-specific static scientific figures.

    Args:
        frames: Canonical tables as data frames.
        fdr_threshold: Configured false-discovery-rate significance threshold.
        root: Report cache root.
        assets: Mutable result-asset mapping.
        inventory: Mutable report inventory.
        feature_counts: Optional pre-aggregated feature coverage. Supplying this
            avoids loading a very large protein-feature membership table.

    Returns:
        Number of logical figures, excluding format variants.

    Raises:
        InputValidationError: If the FDR threshold is invalid.
    """

    validated_fdr_threshold = _validate_fdr_threshold(fdr_threshold=fdr_threshold)
    matplotlib.rcParams["svg.hashsalt"] = "protein-signature-analysis-reports"
    figure_count = 0

    def emit(specification: tuple[Any, str, str, str, str, bool]) -> None:
        """Write and close one figure before constructing the next."""

        nonlocal figure_count
        figure, stem, section, description, content_id, copy_to_final = specification
        try:
            _register_figure_set(
                figure=figure,
                stem=stem,
                section=section,
                description=description,
                content_id=content_id,
                root=root,
                assets=assets,
                inventory=inventory,
                copy_to_final=copy_to_final,
            )
        finally:
            plt.close(figure)
        figure_count += 1

    assignments = frames["label_assignments"]
    assignment_counts = (
        assignments.groupby("label_id", sort=True)
        .size()
        .sort_values(ascending=False, kind="stable")
        .head(40)
        .sort_values(kind="stable")
    )
    figure, axis = plt.subplots(
        figsize=(10.5, max(5.5, 0.25 * max(1, len(assignment_counts)) + 2.0))
    )
    if assignment_counts.empty:
        _draw_no_data(axis=axis, message="No direct label assignments were published")
    else:
        axis.barh(assignment_counts.index.astype(str), assignment_counts.values, color="#3274A1")
        axis.set_xlabel("Directly assigned proteins")
    axis.set_title("Most populated direct protein labels")
    emit(
        (
            figure,
            "00_direct_label_counts",
            "01_proteins_and_curation",
            "Counts for the forty most populated direct labels.",
            "direct_label_counts",
            True,
        )
    )
    del assignments, assignment_counts
    class_summary = frames["class_labelling_summary"]
    figure, axis = plt.subplots(figsize=(10.5, 6.0))
    target_summary = class_summary[class_summary["label_type"].astype(str) == "TARGET"].copy()
    if target_summary.empty:
        _draw_no_data(
            axis=axis,
            message="Automated evidence-led labelling was not selected",
        )
    else:
        target_summary["direct_positive_protein_count"] = pd.to_numeric(
            target_summary["direct_positive_protein_count"], errors="coerce"
        )
        target_summary = target_summary.nlargest(
            30,
            "direct_positive_protein_count",
        ).sort_values("direct_positive_protein_count")
        axis.barh(
            target_summary["label_id"].astype(str),
            target_summary["direct_positive_protein_count"],
            color="#2A9D8F",
        )
        axis.set_xlabel("Evidence-supported proteins")
    axis.set_title("Automated provisional target-label coverage")
    emit(
        (
            figure,
            "01_evidence_supported_label_coverage",
            "01_proteins_and_curation",
            "Evidence-supported protein counts for populated provisional target classes.",
            "evidence_supported_label_coverage",
            True,
        )
    )
    del class_summary, target_summary
    control_matches = frames["control_matching_audit"]
    figure, axis = plt.subplots(figsize=(9.0, 5.8))
    match_scores = pd.to_numeric(
        control_matches.loc[
            control_matches["status"].astype(str) == "MATCHED",
            "match_score",
        ],
        errors="coerce",
    ).dropna()
    if match_scores.empty:
        _draw_no_data(axis=axis, message="No matched-control audit was available")
    else:
        axis.hist(
            match_scores,
            bins=min(30, max(5, int(math.sqrt(len(match_scores))))),
            color="#E9C46A",
            edgecolor="#8C6D1F",
        )
        axis.set_xlabel("Prespecified matching distance (lower is closer)")
        axis.set_ylabel("Target-control matches")
    axis.set_title("Outcome-blind matched-control distance")
    emit(
        (
            figure,
            "02_matched_control_distance",
            "01_proteins_and_curation",
            "Distance distribution for controls passing every prespecified caliper.",
            "matched_control_distance",
            True,
        )
    )
    del control_matches, match_scores
    partitions = frames["partitions"]
    figure, axis = plt.subplots(figsize=(8.5, 5.5))
    if partitions.empty:
        _draw_no_data(axis=axis, message="No analysis partitions were published")
    else:
        partition_counts = partitions.groupby("partition", sort=True).size()
        axis.bar(partition_counts.index.astype(str), partition_counts.values, color="#4C956C")
        axis.set_ylabel("Proteins")
    axis.set_title("Frozen discovery/validation partition sizes")
    emit(
        (
            figure,
            "00_partition_sizes",
            "02_homology_and_partitions",
            "Protein counts in frozen analysis partitions.",
            "partition_sizes",
            True,
        )
    )
    del partitions
    orthogroups = frames["orthofinder_group_context"]
    figure, axis = plt.subplots(figsize=(9.5, 5.8))
    group_sizes = pd.to_numeric(orthogroups.get("member_count"), errors="coerce").dropna()
    if group_sizes.empty:
        _draw_no_data(axis=axis, message="No OrthoFinder group context was available")
    else:
        axis.hist(
            group_sizes, bins=min(30, max(5, int(math.sqrt(len(group_sizes))))), color="#8172B3"
        )
        axis.set_xlabel("Proteins per selected group")
        axis.set_ylabel("Groups")
    axis.set_title("Selected OrthoFinder group-size distribution")
    emit(
        (
            figure,
            "01_orthofinder_group_sizes",
            "02_homology_and_partitions",
            "Member-count distribution for selected OrthoFinder groups.",
            "orthofinder_group_sizes",
            True,
        )
    )
    del orthogroups, group_sizes
    if feature_counts is None:
        features = frames["features"]
        feature_counts = (
            features.groupby("feature_type", sort=True)
            .agg(
                feature_count=("feature_id", "nunique"),
                protein_count=("protein_id", "nunique"),
            )
            .reset_index()
        )
        del features
    figure, axis = plt.subplots(figsize=(10, 5.8))
    if feature_counts.empty:
        _draw_no_data(axis=axis, message="No positive feature evidence was available")
    else:
        axis.bar(feature_counts["feature_type"], feature_counts["protein_count"], color="#3274A1")
        axis.tick_params(axis="x", rotation=35)
        axis.set_ylabel("Proteins")
    axis.set_title("Evidence coverage by feature type")
    emit(
        (
            figure,
            "00_evidence_coverage",
            "03_sequence_and_domains",
            "Evidence coverage by feature type.",
            "evidence_coverage",
            True,
        )
    )
    signatures = frames["signatures"]
    figure, axis = plt.subplots(figsize=(10, 5.8))
    if signatures.empty:
        _draw_no_data(axis=axis, message="No signature summaries were published")
    else:
        signature_counts = signatures.groupby("evidence_class", sort=True).size().sort_values()
        axis.barh(signature_counts.index.astype(str), signature_counts.values, color="#C44E52")
        axis.set_xlabel("Signatures")
    axis.set_title("Signature evidence classes across all comparisons")
    emit(
        (
            figure,
            "00_signature_evidence_classes",
            "05_association_statistics",
            "Counts of discovery and validation evidence classes.",
            "signature_evidence_classes",
            True,
        )
    )
    del signatures
    models = frames["ml_models"]
    figure, axis = plt.subplots(figsize=(10, 5.8))
    if models.empty:
        _draw_no_data(axis=axis, message="No explainable-model statuses were published")
    else:
        model_counts = models.groupby("status", sort=True).size().sort_values()
        axis.barh(model_counts.index.astype(str), model_counts.values, color="#937860")
        axis.set_xlabel("Comparisons")
    axis.set_title("Explainable-model completion status")
    emit(
        (
            figure,
            "00_model_statuses",
            "06_explainable_models",
            "Completion and non-fitted status counts for mandatory models.",
            "model_statuses",
            True,
        )
    )
    domains = frames["domain_assessments"]
    figure, axis = plt.subplots(figsize=(10, 5.8))
    if domains.empty:
        _draw_no_data(axis=axis, message="No domain assessments were published")
    else:
        domain_counts = (
            domains.groupby(["domain_authority", "assessment_status"], sort=True)
            .size()
            .unstack(fill_value=0)
        )
        domain_counts.plot(kind="bar", stacked=True, ax=axis, colormap="Set2")
        axis.set_xlabel("Domain authority")
        axis.set_ylabel("Proteins")
        axis.legend(title="Assessment", bbox_to_anchor=(1.02, 1), loc="upper left")
    axis.set_title("Domain assessment coverage")
    emit(
        (
            figure,
            "01_domain_assessment_coverage",
            "03_sequence_and_domains",
            "Hit, no-hit, not-assessed and failed domain coverage.",
            "domain_assessment_coverage",
            True,
        )
    )
    del domains
    structures = frames["structures"]
    figure, axis = plt.subplots(figsize=(10, 5.8))
    if structures.empty:
        _draw_no_data(axis=axis, message="No structure inventory was published")
    else:
        structure_counts = (
            structures.groupby(["structure_source", "availability_status"], sort=True)
            .size()
            .unstack(fill_value=0)
        )
        structure_counts.plot(kind="bar", stacked=True, ax=axis, colormap="Paired")
        axis.set_xlabel("Structure source")
        axis.set_ylabel("Models")
        axis.legend(title="Status", bbox_to_anchor=(1.02, 1), loc="upper left")
    axis.set_title("Structure evidence coverage")
    emit(
        (
            figure,
            "00_structure_evidence_coverage",
            "04_structures_and_folds",
            "Structure source and availability coverage.",
            "structure_evidence_coverage",
            True,
        )
    )
    del structures
    comparison_rows = frames["comparisons"][["comparison_id", "display_name"]]
    association = frames["associations"]
    importance = frames["ml_feature_importance"]
    predictions = frames["ml_predictions"]
    for comparison_row in comparison_rows.itertuples(index=False):
        comparison_id = str(comparison_row.comparison_id)
        raw_display_name = comparison_row.display_name
        display_name = "" if pd.isna(raw_display_name) else str(raw_display_name).strip()
        if not display_name:
            display_name = comparison_id
        subset = association[
            (association["comparison_id"].astype(str) == comparison_id)
            & (association["partition"].astype(str) == "DISCOVERY")
            & (association["status"].astype(str) == "COMPLETE")
        ].copy()
        if not subset.empty:
            figure = _association_figure(
                frame=subset,
                comparison_display_name=display_name,
                fdr_threshold=validated_fdr_threshold,
            )
            emit(
                (
                    figure,
                    f"{_safe_report_token(value=comparison_id)}_association_landscape",
                    "05_association_statistics",
                    "Discovery effect size against local false-discovery rate.",
                    f"{comparison_id}:association_landscape",
                    True,
                )
            )
            figure = _prevalence_figure(
                frame=subset,
                comparison_display_name=display_name,
            )
            emit(
                (
                    figure,
                    f"{_safe_report_token(value=comparison_id)}_top_prevalence",
                    "05_association_statistics",
                    "Target and background prevalence for prioritised discovery features.",
                    f"{comparison_id}:top_prevalence",
                    True,
                )
            )
        model_rows = models[models["comparison_id"].astype(str) == comparison_id]
        complete = not model_rows.empty and str(model_rows.iloc[0]["status"]).startswith("COMPLETE")
        selected_importance = importance[importance["comparison_id"].astype(str) == comparison_id]
        if complete and not selected_importance.empty:
            figure = _model_importance_figure(
                frame=selected_importance,
                comparison_display_name=display_name,
            )
            emit(
                (
                    figure,
                    f"{_safe_report_token(value=comparison_id)}_model_coefficients",
                    "06_explainable_models",
                    "Largest absolute elastic-net coefficients.",
                    f"{comparison_id}:model_coefficients",
                    True,
                )
            )
        selected_predictions = predictions[
            predictions["comparison_id"].astype(str) == comparison_id
        ]
        if complete and not selected_predictions.empty:
            figure = _prediction_figure(
                frame=selected_predictions,
                comparison_display_name=display_name,
            )
            emit(
                (
                    figure,
                    f"{_safe_report_token(value=comparison_id)}_prediction_distributions",
                    "06_explainable_models",
                    "Discovery and validation predicted-probability distributions.",
                    f"{comparison_id}:prediction_distributions",
                    True,
                )
            )
    del association, comparison_rows, importance, models, predictions
    return figure_count


def _association_figure(
    *,
    frame: pd.DataFrame,
    comparison_display_name: str,
    fdr_threshold: float,
) -> Any:
    """Create an effect-versus-FDR discovery landscape.

    Args:
        frame: Complete discovery association rows.
        comparison_display_name: Human-readable comparison name for the title.
        fdr_threshold: Configured false-discovery-rate significance threshold.

    Returns:
        Matplotlib figure.

    Raises:
        InputValidationError: If the FDR threshold is invalid.
    """

    validated_fdr_threshold = _validate_fdr_threshold(fdr_threshold=fdr_threshold)
    figure, axis = plt.subplots(figsize=(10.5, 6.5))
    prepared = frame.dropna(subset=["prevalence_difference", "q_value"]).copy()
    if prepared.empty:
        _draw_no_data(axis=axis, message="No finite discovery association values")
    else:
        prepared["negative_log10_q"] = prepared["q_value"].map(
            lambda value: -math.log10(max(float(value), 1e-300))
        )
        for feature_type in sorted(prepared["feature_type"].astype(str).unique()):
            rows = prepared[prepared["feature_type"].astype(str) == feature_type]
            axis.scatter(
                rows["prevalence_difference"],
                rows["negative_log10_q"],
                label=feature_type,
                alpha=0.72,
                s=26,
            )
        axis.axvline(0.0, color="#555555", linewidth=0.8)
        axis.axhline(
            -math.log10(validated_fdr_threshold),
            color="#B22222",
            linestyle="--",
            linewidth=1,
        )
        axis.set_xlabel("Target - background prevalence")
        axis.set_ylabel("-log10(local discovery q-value)")
        if prepared["feature_type"].nunique() <= 12:
            axis.legend(title="Feature type", bbox_to_anchor=(1.02, 1), loc="upper left")
    axis.set_title(f"Discovery association landscape\n{comparison_display_name}")
    return figure


def _prevalence_figure(*, frame: pd.DataFrame, comparison_display_name: str) -> Any:
    """Create a paired prevalence plot for the strongest discovery features.

    Args:
        frame: Complete discovery association rows.
        comparison_display_name: Human-readable comparison name for the title.

    Returns:
        Matplotlib figure.
    """

    prepared = frame.dropna(subset=["target_prevalence", "background_prevalence"]).copy()
    prepared["ranking_q"] = prepared["study_q_value"].fillna(prepared["q_value"]).fillna(1.0)
    prepared = prepared.sort_values(
        ["ranking_q", "q_value", "feature_type", "feature_id"],
        kind="stable",
    ).head(25)
    height = max(5.5, 0.30 * max(1, len(prepared)) + 2.0)
    figure, axis = plt.subplots(figsize=(11.5, height))
    if prepared.empty:
        _draw_no_data(axis=axis, message="No finite prevalence estimates")
    else:
        labels = [
            f"{row.feature_type} | {str(row.feature_name)[:55]}"
            for row in prepared.itertuples(index=False)
        ]
        positions = list(range(len(prepared)))
        axis.barh(
            [position - 0.18 for position in positions],
            prepared["target_prevalence"],
            height=0.34,
            label="Target",
            color="#3274A1",
        )
        axis.barh(
            [position + 0.18 for position in positions],
            prepared["background_prevalence"],
            height=0.34,
            label="Background",
            color="#E1812C",
        )
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.set_xlim(0, 1)
        axis.set_xlabel("Fraction of independent blocks")
        axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    axis.set_title(f"Prioritised discovery feature prevalence\n{comparison_display_name}")
    return figure


def _model_importance_figure(*, frame: pd.DataFrame, comparison_display_name: str) -> Any:
    """Create a signed coefficient plot for one fitted model.

    Args:
        frame: Model feature-importance rows.
        comparison_display_name: Human-readable comparison name for the title.

    Returns:
        Matplotlib figure.
    """

    prepared = frame.copy()
    prepared["absolute"] = prepared["coefficient_log_odds"].abs()
    prepared = prepared.sort_values(
        ["absolute", "feature_type", "feature_id"],
        ascending=[False, True, True],
        kind="stable",
    ).head(30)
    prepared = prepared.sort_values("coefficient_log_odds", kind="stable")
    height = max(5.5, 0.27 * max(1, len(prepared)) + 2.0)
    figure, axis = plt.subplots(figsize=(11.5, height))
    labels = [
        f"{row.feature_type} | {str(row.feature_name)[:55]}"
        for row in prepared.itertuples(index=False)
    ]
    colours = ["#3274A1" if value >= 0 else "#C44E52" for value in prepared["coefficient_log_odds"]]
    axis.barh(labels, prepared["coefficient_log_odds"], color=colours)
    axis.axvline(0.0, color="#555555", linewidth=0.8)
    axis.set_xlabel("Elastic-net coefficient (log odds)")
    axis.set_title(f"Largest model coefficients\n{comparison_display_name}")
    return figure


def _prediction_figure(*, frame: pd.DataFrame, comparison_display_name: str) -> Any:
    """Create class- and partition-aware probability histograms.

    Args:
        frame: Model prediction rows.
        comparison_display_name: Human-readable comparison name for the title.

    Returns:
        Matplotlib figure.
    """

    figure, axis = plt.subplots(figsize=(10.5, 6.2))
    colours = {
        ("DISCOVERY", "TARGET"): "#3274A1",
        ("DISCOVERY", "BACKGROUND"): "#E1812C",
        ("VALIDATION", "TARGET"): "#5BA4CF",
        ("VALIDATION", "BACKGROUND"): "#F0A35E",
    }
    for partition in ("DISCOVERY", "VALIDATION"):
        for true_class in ("TARGET", "BACKGROUND"):
            rows = frame[
                (frame["partition"].astype(str) == partition)
                & (frame["true_class"].astype(str) == true_class)
            ]
            if rows.empty:
                continue
            axis.hist(
                rows["predicted_probability"].astype(float),
                bins=20,
                range=(0, 1),
                alpha=0.52,
                label=f"{partition.title()} {true_class.title()}",
                color=colours[(partition, true_class)],
            )
    axis.axvline(0.5, color="#555555", linestyle="--", linewidth=1)
    axis.set_xlim(0, 1)
    axis.set_xlabel("Predicted target probability")
    axis.set_ylabel("Proteins")
    handles, labels = axis.get_legend_handles_labels()
    if handles:
        axis.legend(handles, labels)
    axis.set_title(f"Model probability distributions\n{comparison_display_name}")
    return figure


def _register_figure_set(
    *,
    figure: Any,
    stem: str,
    section: str,
    description: str,
    content_id: str,
    root: Path,
    assets: dict[str, Path],
    inventory: list[dict[str, Any]],
    copy_to_final: bool,
) -> None:
    """Save PNG/SVG/PDF forms and optionally register final-result copies.

    Args:
        figure: Matplotlib figure.
        stem: Safe base filename.
        section: Owning numbered report section.
        description: Human-readable figure meaning.
        content_id: Stable inventory identity.
        root: Report cache root.
        assets: Mutable result-asset mapping.
        inventory: Mutable report inventory.
        copy_to_final: Also expose the same source in `99_final_results`.
    """

    figure.tight_layout()
    for file_format in _FIGURE_FORMATS:
        relative = f"analysis/{section}/figures/{stem}.{file_format}"
        if relative in assets:
            raise PublicationError(f"Static figure path is duplicated: {relative!r}")
        source = root / relative
        _save_figure_atomic(figure=figure, path=source, file_format=file_format)
        assets[relative] = source
        LOGGER.debug("Built %s static figure %s", file_format.upper(), relative)
        inventory.append(
            _inventory_row(
                section=section,
                asset_kind="FIGURE",
                content_id=content_id,
                file_format=file_format.upper(),
                relative_path=relative,
                row_count=0,
                description=description,
            )
        )
        if copy_to_final:
            final_relative = f"analysis/99_final_results/figures/{stem}.{file_format}"
            if final_relative in assets:
                raise PublicationError(f"Final figure path is duplicated: {final_relative!r}")
            assets[final_relative] = source
            inventory.append(
                _inventory_row(
                    section="99_final_results",
                    asset_kind="FIGURE",
                    content_id=content_id,
                    file_format=file_format.upper(),
                    relative_path=final_relative,
                    row_count=0,
                    description=f"Decision-facing copy. {description}",
                )
            )


def _save_figure_atomic(*, figure: Any, path: Path, file_format: str) -> None:
    """Write one deterministic Matplotlib figure through an atomic rename.

    Args:
        figure: Matplotlib figure.
        path: Final cache destination.
        file_format: One of PNG, SVG or PDF in lower case.

    Raises:
        PublicationError: If the format or write is invalid.
    """

    if file_format not in _FIGURE_FORMATS:
        raise PublicationError(f"Unsupported static report format: {file_format!r}")
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=f".{file_format}",
    )
    os.close(descriptor)
    if file_format == "png":
        metadata = {"Software": "protein-signature-analysis"}
    elif file_format == "svg":
        metadata = {"Creator": "protein-signature-analysis", "Date": None}
    else:
        metadata = {
            "Creator": "protein-signature-analysis",
            "CreationDate": None,
            "ModDate": None,
        }
    try:
        figure.savefig(
            temporary_name,
            format=file_format,
            dpi=180,
            bbox_inches="tight",
            metadata=metadata,
        )
        os.replace(temporary_name, destination)
    except (OSError, TypeError, ValueError) as error:
        Path(temporary_name).unlink(missing_ok=True)
        raise PublicationError(
            f"Could not write static report figure {destination}: {error}"
        ) from error


def _write_binary_atomic(*, path: Path, payload: bytes) -> None:
    """Write non-empty binary report content through an atomic rename.

    Args:
        path: Final cache destination.
        payload: Non-empty binary content.

    Raises:
        PublicationError: If content is empty or cannot be written.
    """

    if not payload:
        raise PublicationError(f"Refusing to write an empty report file: {path}")
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    except OSError as error:
        Path(temporary_name).unlink(missing_ok=True)
        raise PublicationError(f"Could not write report file {destination}: {error}") from error


def _draw_no_data(*, axis: Any, message: str) -> None:
    """Draw a clear unavailable-evidence state on an otherwise empty axis.

    Args:
        axis: Matplotlib axis.
        message: Human-readable reason no marks are displayed.
    """

    axis.text(0.5, 0.5, message, ha="center", va="center", transform=axis.transAxes)
    axis.set_xticks([])
    axis.set_yticks([])


def _safe_report_token(*, value: str) -> str:
    """Create a bounded readable filename token with a collision suffix.

    Args:
        value: Stable scientific identifier.

    Returns:
        Filesystem-safe token.

    Raises:
        InputValidationError: If the identifier is blank.
    """

    text = str(value).strip()
    if not text:
        raise InputValidationError("A report identifier must not be blank.")
    readable = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._-") or "item"
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()[:10]
    return f"{readable[:80]}_{digest}"


def _final_shap_relative_path(*, relative_path: str) -> str:
    """Map one native SHAP result path into the final-results subtree.

    Args:
        relative_path: Native result-relative SHAP asset path.

    Returns:
        Collision-resistant result-relative final-results path.

    Raises:
        InputValidationError: If the path is not within the declared SHAP root.
    """

    source = Path(relative_path)
    prefix = Path("analysis") / "06_explainable_models" / "figures" / "shap"
    try:
        suffix = source.relative_to(prefix)
    except ValueError as error:
        raise InputValidationError(
            f"SHAP asset is outside the numbered report hierarchy: {relative_path!r}."
        ) from error
    if not suffix.parts or any(part in {"", ".", ".."} for part in suffix.parts):
        raise InputValidationError(f"SHAP asset path is invalid: {relative_path!r}.")
    return (Path("analysis") / "99_final_results" / "figures" / "shap" / suffix).as_posix()


def _inventory_row(
    *,
    section: str,
    asset_kind: str,
    content_id: str,
    file_format: str,
    relative_path: str,
    row_count: int,
    description: str,
) -> dict[str, Any]:
    """Build one stable human-report inventory record.

    Args:
        section: Numbered report section.
        asset_kind: TABLE, FIGURE or DOCUMENTATION.
        content_id: Stable content identity.
        file_format: Upper-case file format.
        relative_path: Portable result-relative path.
        row_count: Table rows, or zero for non-tabular assets.
        description: Human-readable meaning.

    Returns:
        Inventory record.
    """

    return {
        "section": section,
        "asset_kind": asset_kind,
        "content_id": content_id,
        "file_format": file_format,
        "relative_path": relative_path,
        "row_count": row_count,
        "status": "COMPLETE",
        "description": description,
    }
