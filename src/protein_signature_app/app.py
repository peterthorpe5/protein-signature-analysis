"""Streamlit user interface for sequence, domain, fold and structure evidence."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from protein_signature_app.backend import (
    canonical_table_names,
    canonical_table_preview,
    canonical_table_storage,
    distinct_values,
    filtered_feature_exports,
    load_canonical_table_assets,
    load_metadata,
    load_report_inventory_assets,
    query_dataframe,
    resolve_database_for_app,
    result_inventory_identity,
    table_count,
)
from protein_signature_app.structure_viewer import (
    ANNOTATION_COLUMNS,
    ModelTrace,
    align_sequences,
    alignment_fasta,
    canonical_accession,
    enrichment_track,
    external_links,
    feature_intervals,
    fetch_alphafold_model,
    pair_links,
    parse_annotation_tsv,
    parse_mmcif_trace,
    parse_pdb_trace,
    project_kmer_intervals,
    read_published_model,
    significant_intervals,
)
from protein_signature_app.viewer_help import GLOSSARY, PAGE_HELP
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.exports import (
    dataframe_to_tsv_bytes,
    dataframe_to_xlsx_bytes,
    normalise_dataframe,
    plotly_figure_to_pdf_bytes,
    safe_download_stem,
)

_MAX_APP_ASSET_BYTES = 100 * 1024 * 1024
_ASSET_FORMATS = {
    "PNG": (".png", "image/png"),
    "SVG": (".svg", "image/svg+xml"),
    "PDF": (".pdf", "application/pdf"),
}


@st.cache_resource(show_spinner="Verifying the completed result")
def _verified_app_database(*, resource: str, inventory_identity: str) -> Path:
    """Verify once per immutable result snapshot across interactive reruns.

    Args:
        resource: Completed result path supplied to the application.
        inventory_identity: File metadata and inode digest that invalidates
            the cache whenever the result inventory changes.

    Returns:
        Fully verified read-only DuckDB path.
    """

    if not inventory_identity:
        raise InputValidationError("Could not identify the result inventory.")
    return resolve_database_for_app(resource=Path(resource), inventory_identity=inventory_identity)


def main() -> None:
    """Render the complete read-only result interrogation application."""

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--resource", required=True, type=Path)
    arguments, _ = parser.parse_known_args()
    st.set_page_config(
        page_title="Protein Signature Analysis",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.markdown(
        """
        <style>
        .block-container {padding-top: 1.6rem; padding-bottom: 3rem;}
        [data-testid="stMetric"] {border: 1px solid #dbe4ea; border-radius: 12px;
          padding: 0.8rem; background: linear-gradient(145deg, #f8fbfc, #ffffff);}
        </style>
        """,
        unsafe_allow_html=True,
    )
    try:
        inventory_identity = result_inventory_identity(resource=arguments.resource)
        database = _verified_app_database(
            resource=str(arguments.resource.expanduser().resolve()),
            inventory_identity=inventory_identity,
        )
        metadata = load_metadata(database=database)
    except Exception as error:
        st.error(f"Could not open the result: {error}")
        st.stop()
    campaign = metadata.get("campaign", {}).get("campaign", {})
    st.sidebar.title("Protein signatures")
    st.sidebar.caption(f"Campaign: {campaign.get('campaign_id', 'unknown')}")
    st.sidebar.caption(f"Package: {metadata.get('package_version', 'unknown')}")
    page = st.sidebar.radio(
        "Explore",
        (
            "Overview",
            "Signature explorer",
            "Explainable prediction",
            "Protein & Pfam",
            "Classes & roles",
            "Structures & folds",
            "Model & alignment explorer",
            "Orthology & partitions",
            "Canonical data & downloads",
            "Data quality & provenance",
            "Glossary & help",
        ),
    )
    with st.sidebar.expander("How to use this page"):
        st.markdown(PAGE_HELP[page])
    if page == "Overview":
        _render_overview(
            database=database, metadata=metadata, inventory_identity=inventory_identity
        )
    elif page == "Signature explorer":
        _render_signatures(database=database)
    elif page == "Explainable prediction":
        _render_explainable_ml(database=database)
    elif page == "Protein & Pfam":
        _render_proteins(database=database)
    elif page == "Classes & roles":
        _render_classes(database=database)
    elif page == "Structures & folds":
        _render_structures(database=database)
    elif page == "Model & alignment explorer":
        _render_model_explorer(database=database)
    elif page == "Orthology & partitions":
        _render_orthology(database=database)
    elif page == "Canonical data & downloads":
        _render_canonical_data(database=database)
    elif page == "Data quality & provenance":
        _render_quality(database=database, metadata=metadata)
    else:
        _render_glossary()


def _render_canonical_data(*, database: Path) -> None:
    """Render every canonical dataset with complete published downloads.

    Args:
        database: Verified result database.
    """

    st.title("Canonical data & downloads")
    names = canonical_table_names()
    st.caption(
        "Browse a bounded preview and inspect the complete checksum-verified storage files. "
        "Manageable tables have full formatted Excel workbooks; very large tables have a "
        f"compact workbook and filtered exporter. All {len(names)} datasets are available."
    )
    table_name = st.selectbox("Canonical dataset", names)
    row_count = table_count(database=database, table_name=table_name)
    st.subheader(str(table_name).replace("_", " ").title())
    st.caption(f"{row_count:,} rows; previewing at most 500 rows below.")
    preview = canonical_table_preview(
        database=database,
        table_name=table_name,
        limit=500,
        offset=0,
    )
    _render_downloadable_table(frame=preview, download_name=f"{table_name}_visible_preview")
    try:
        storage = canonical_table_storage(database=database, table_name=table_name)
    except InputValidationError as error:
        st.error(f"Could not inspect canonical storage: {error}")
        return
    st.caption("Complete canonical storage (paths are relative to the result directory).")
    _render_downloadable_table(frame=storage, download_name=f"{table_name}_storage_inventory")
    try:
        assets = load_canonical_table_assets(database=database, table_name=table_name)
    except InputValidationError as error:
        st.error(f"Could not prepare canonical downloads: {error}")
        return
    for asset in assets:
        st.download_button(
            label=(
                f"Download complete {asset.file_format}"
                if asset.complete
                else f"Download compact {asset.file_format} summary"
            ),
            data=asset.payload,
            file_name=Path(asset.relative_path).name,
            mime=asset.mime_type,
            key=f"canonical_{table_name}_{asset.file_format}",
        )
    if table_name == "features":
        _render_filtered_feature_export(database=database)
    st.subheader("Complete report inventory")
    st.caption("The inventory indexes every numbered table, figure and documentation asset.")
    try:
        inventory_assets = load_report_inventory_assets(database=database)
    except InputValidationError as error:
        st.error(f"Could not prepare report-inventory downloads: {error}")
        return
    for asset in inventory_assets:
        st.download_button(
            label=f"Download report inventory as {asset.file_format}",
            data=asset.payload,
            file_name=Path(asset.relative_path).name,
            mime=asset.mime_type,
            key=f"report_inventory_{asset.file_format}",
        )


def _render_filtered_feature_export(*, database: Path) -> None:
    """Render bounded filters for the complete feature-membership relation.

    Args:
        database: Verified result database.
    """

    st.subheader("Filtered feature export")
    st.caption(
        "Select at least one filter. Exports are generated from DuckDB and capped before "
        "conversion to TSV and formatted Excel."
    )
    available_types = distinct_values(
        database=database,
        table_name="features",
        column_name="feature_type",
    )
    selected_types = tuple(st.multiselect("Feature types to export", available_types))
    protein_id = st.text_input("Exact protein identifier (optional)")
    feature_id_contains = st.text_input("Feature identifier contains (optional)")
    maximum_rows = st.selectbox(
        "Maximum exported rows",
        (10_000, 50_000, 100_000, 250_000),
        index=2,
    )
    if not st.button("Prepare filtered feature downloads"):
        return
    try:
        exports = filtered_feature_exports(
            database=database,
            feature_types=selected_types,
            protein_id=protein_id,
            feature_id_contains=feature_id_contains,
            maximum_rows=int(maximum_rows),
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    st.success(f"Prepared {exports[0].row_count:,} matching feature rows.")
    for asset in exports:
        st.download_button(
            label=f"Download filtered features as {asset.file_format}",
            data=asset.payload,
            file_name=Path(asset.relative_path).name,
            mime=asset.mime_type,
            key=f"filtered_features_{asset.file_format}",
        )


def _overview_count(*, database: Path, metadata: dict[str, object], table_name: str) -> int:
    """Use verified published counts without rescanning a large physical table.

    Args:
        database: Verified result database.
        metadata: Published run metadata.
        table_name: Canonical table whose count is requested.

    Returns:
        Published non-negative row count, or the queried count when unavailable.
    """

    counts = metadata.get("counts")
    if isinstance(counts, dict):
        value = counts.get(table_name)
        if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
            return value
    return table_count(database=database, table_name=table_name)


def _render_overview(
    *, database: Path, metadata: dict[str, object], inventory_identity: str = ""
) -> None:
    """Render campaign metrics and evidence composition.

    Args:
        database: Verified result database.
        metadata: Run metadata.
        inventory_identity: Current verified inventory for optional full coverage.
    """

    st.title("Protein signature analysis")
    st.caption("Sequence · domains · folds · pairwise structural evidence")
    st.info(
        "Signatures are prioritisation evidence, not proof of biochemical activity. "
        "Discovery and held-out validation are kept separate."
    )
    columns = st.columns(4)
    columns[0].metric(
        "Proteins",
        f"{_overview_count(database=database, metadata=metadata, table_name='proteins'):,}",
    )
    columns[1].metric(
        "Candidate signatures",
        f"{_overview_count(database=database, metadata=metadata, table_name='signatures'):,}",
    )
    columns[2].metric(
        "Structure models",
        f"{_overview_count(database=database, metadata=metadata, table_name='structures'):,}",
    )
    columns[3].metric(
        "Domain hits",
        f"{_overview_count(database=database, metadata=metadata, table_name='domain_hits'):,}",
    )
    signature_counts = query_dataframe(
        database=database,
        sql=(
            "SELECT feature_type, count(*) AS signature_count FROM signatures "
            "WHERE status = 'COMPLETE' GROUP BY feature_type "
            "ORDER BY signature_count DESC, feature_type"
        ),
    )
    left, right = st.columns((3, 2))
    with left:
        st.subheader("Completed signatures by feature type")
        if signature_counts.empty:
            st.info("No completed signature evidence was available.")
        else:
            figure = px.bar(
                signature_counts,
                x="feature_type",
                y="signature_count",
                color="feature_type",
                labels={
                    "signature_count": "Completed signatures",
                    "feature_type": "Feature type",
                },
            )
            _render_plotly_figure(figure=figure, download_name="overview_signature_types")
    with right:
        st.subheader("Availability states")
        availability = metadata.get("evidence_availability", {})
        if isinstance(availability, dict):
            _render_downloadable_table(
                frame={"evidence": list(availability), "status": list(availability.values())},
                download_name="overview_availability_states",
            )
    with st.expander("Exact feature coverage (large calculation)", expanded=False):
        st.caption(
            "Calculating distinct features and proteins for each feature type scans "
            "the complete feature table and may take several minutes."
        )
        cache_key = f"overview-feature-coverage-{database}-{inventory_identity}"
        if st.button("Calculate exact feature coverage"):
            try:
                st.session_state[cache_key] = query_dataframe(
                    database=database,
                    sql=(
                        "SELECT feature_type, count(DISTINCT feature_id) AS feature_count, "
                        "count(DISTINCT protein_id) AS protein_count FROM features "
                        "GROUP BY feature_type ORDER BY protein_count DESC, feature_type"
                    ),
                )
            except InputValidationError as error:
                st.warning(f"Could not calculate full feature coverage: {error}")
        coverage = st.session_state.get(cache_key)
        if coverage is not None:
            _render_downloadable_table(
                frame=coverage,
                download_name="overview_exact_feature_coverage",
            )
            if not coverage.empty:
                figure = px.bar(
                    coverage,
                    x="feature_type",
                    y="protein_count",
                    color="feature_count",
                    labels={"protein_count": "Proteins", "feature_type": "Feature type"},
                )
                _render_plotly_figure(
                    figure=figure, download_name="overview_exact_feature_coverage_chart"
                )


def _render_signatures(*, database: Path) -> None:
    """Render filterable discovery and validation signature evidence.

    Args:
        database: Verified result database.
    """

    st.title("Signature explorer")
    comparisons = distinct_values(
        database=database, table_name="signatures", column_name="comparison_id"
    )
    feature_types = distinct_values(
        database=database, table_name="signatures", column_name="feature_type"
    )
    if not comparisons:
        st.info("No comparisons are present in this result.")
        return
    completed = query_dataframe(
        database=database,
        sql=(
            "SELECT comparison_id, count(*) AS enriched_count FROM signatures "
            "WHERE status = 'COMPLETE' AND discovery_prevalence_difference > 0 "
            "AND discovery_q_value <= 0.05 GROUP BY comparison_id "
            "ORDER BY enriched_count DESC, comparison_id"
        ),
    )
    counts = (
        {
            str(row.comparison_id): int(row.enriched_count)
            for row in completed.itertuples(index=False)
        }
        if {"comparison_id", "enriched_count"}.issubset(completed.columns)
        else {}
    )
    comparisons = tuple(counts) + tuple(value for value in comparisons if value not in counts)
    st.caption(
        "Target prevalence describes how common a feature is. The q-value controls FDR "
        "within this comparison and evidence family; study q-value also corrects across "
        "all configured comparisons in the same evidence family. "
        f"{len(counts)} of {len(comparisons)} comparisons have completed enriched signatures."
    )
    selected_comparison = st.selectbox(
        "Comparison",
        comparisons,
        format_func=lambda value: (
            f"{value} · {counts[value]:,} enriched" if value in counts else value
        ),
    )
    selected_types = st.multiselect("Feature types", feature_types, default=feature_types)
    if not selected_types:
        st.info("Select at least one feature type.")
        return
    placeholders = ",".join("?" for _ in selected_types)
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM signatures WHERE comparison_id = ? "
            f"AND feature_type IN ({placeholders}) "
            "ORDER BY discovery_q_value NULLS LAST, feature_type, feature_id"
        ),
        parameters=(selected_comparison, *selected_types),
    )
    _render_downloadable_table(
        frame=frame,
        download_name=f"{selected_comparison}_signatures",
        height=480,
    )
    if not counts.get(selected_comparison):
        st.info("No completed positive enriched signature is available for this comparison.")
    chart_data = frame.dropna(subset=["discovery_prevalence_difference"])
    if not chart_data.empty:
        figure = px.scatter(
            chart_data,
            x="discovery_prevalence_difference",
            y="discovery_q_value",
            color="feature_type",
            symbol="evidence_class",
            hover_data=["feature_id", "feature_name"],
            labels={
                "discovery_prevalence_difference": "Target − background prevalence",
                "discovery_q_value": "Discovery q-value",
            },
        )
        figure.update_yaxes(autorange="reversed")
        _render_plotly_figure(
            figure=figure,
            download_name=f"{selected_comparison}_signature_scatter",
        )
    with st.expander("Association counts, intervals and multiplicity", expanded=False):
        association_frame = query_dataframe(
            database=database,
            sql=(
                "SELECT a.* "
                "FROM associations a JOIN signatures s USING "
                "(comparison_id, feature_type, feature_id) "
                "WHERE a.comparison_id = ? "
                f"AND a.feature_type IN ({placeholders}) "
                "ORDER BY a.partition, a.study_q_value NULLS LAST, a.feature_type, "
                "a.feature_id LIMIT 10000"
            ),
            parameters=(selected_comparison, *selected_types),
        )
        _render_downloadable_table(
            frame=association_frame,
            download_name=f"{selected_comparison}_association_counts",
            height=440,
        )


def _render_explainable_ml(*, database: Path) -> None:
    """Render held-out prediction, global importance and local explanations.

    Args:
        database: Verified result database.
    """

    st.title("Explainable prediction")
    st.info(
        "Prediction is a separate corroborating layer, not a replacement for association "
        "testing and not evidence of biochemical causation. Validation groups never enter "
        "model fitting or feature selection."
    )
    comparisons = distinct_values(
        database=database,
        table_name="ml_models",
        column_name="comparison_id",
    )
    if not comparisons:
        st.info("No configured comparison is present in the modelling results.")
        return
    comparison_id = st.selectbox("Model comparison", comparisons)
    model = query_dataframe(
        database=database,
        sql="SELECT * FROM ml_models WHERE comparison_id = ?",
        parameters=(comparison_id,),
    )
    _render_downloadable_table(
        frame=model,
        download_name=f"{comparison_id}_model_summary",
    )
    if model.empty or not str(model.iloc[0]["status"]).startswith("COMPLETE"):
        st.warning("This comparison has an explicit non-fitted model status.")
        return
    columns = st.columns(4)
    columns[0].metric("CV ROC AUC", _metric_text(model.iloc[0]["cv_roc_auc"]))
    columns[1].metric(
        "Validation ROC AUC",
        _metric_text(model.iloc[0]["validation_roc_auc"]),
    )
    columns[2].metric(
        "Validation PR AUC",
        _metric_text(model.iloc[0]["validation_average_precision"]),
    )
    columns[3].metric(
        "Validation MCC",
        _metric_text(model.iloc[0]["validation_matthews_correlation"]),
    )
    global_plots = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM ml_plot_inventory WHERE comparison_id = ? "
            "AND plot_type IN ('SHAP_BEESWARM', 'SHAP_GLOBAL_BAR') "
            "ORDER BY plot_type, file_format"
        ),
        parameters=(comparison_id,),
    )
    _render_shap_assets(
        database=database,
        inventory=global_plots,
        heading="SHAP graphical explanations",
    )
    importance = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM ml_feature_importance WHERE comparison_id = ? "
            "ORDER BY importance_rank LIMIT 200"
        ),
        parameters=(comparison_id,),
    )
    st.subheader("Global model evidence")
    st.caption(
        "Coefficients are signed log-odds effects. Held-out permutation importance is "
        "model-specific and is only calculated when both validation classes are present."
    )
    _render_downloadable_table(
        frame=importance,
        download_name=f"{comparison_id}_model_feature_importance",
        height=440,
    )
    if not importance.empty:
        chart = importance.head(30).sort_values("coefficient_log_odds")
        figure = px.bar(
            chart,
            x="coefficient_log_odds",
            y="feature_name",
            color="feature_type",
            orientation="h",
            hover_data=["feature_id", "validation_permutation_importance_mean"],
        )
        _render_plotly_figure(
            figure=figure,
            download_name=f"{comparison_id}_model_coefficients",
        )
    explained_partition = str(model.iloc[0].get("shap_explained_partition") or "VALIDATION")
    predictions = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM ml_predictions WHERE comparison_id = ? AND partition = ? "
            "ORDER BY predicted_probability DESC, protein_id"
        ),
        parameters=(comparison_id, explained_partition),
    )
    st.subheader(f"{explained_partition.title()} predictions")
    _render_downloadable_table(
        frame=predictions,
        download_name=f"{comparison_id}_{explained_partition}_predictions",
        height=360,
    )
    if predictions.empty:
        return
    figure = px.histogram(
        predictions,
        x="predicted_probability",
        color="true_class",
        barmode="overlay",
        nbins=20,
        labels={"predicted_probability": "Predicted target probability"},
    )
    _render_plotly_figure(
        figure=figure,
        download_name=f"{comparison_id}_{explained_partition}_probabilities",
    )
    protein_id = st.selectbox(
        "Explain validation protein",
        tuple(str(value) for value in predictions["protein_id"].tolist()),
    )
    local = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM ml_explanations WHERE comparison_id = ? AND protein_id = ? "
            "ORDER BY absolute_rank"
        ),
        parameters=(comparison_id, protein_id),
    )
    st.subheader("Local additive explanation")
    st.caption(
        "Linear SHAP contributions use the discovery background and are additive in "
        "model log-odds space; correlated features must be interpreted together."
    )
    _render_downloadable_table(
        frame=local,
        download_name=f"{comparison_id}_{protein_id}_shap_values",
    )
    local_plots = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM ml_plot_inventory WHERE comparison_id = ? AND protein_id = ? "
            "AND plot_type = 'SHAP_WATERFALL' ORDER BY file_format"
        ),
        parameters=(comparison_id, protein_id),
    )
    _render_shap_assets(
        database=database,
        inventory=local_plots,
        heading="SHAP waterfall",
    )


def _render_shap_assets(*, database: Path, inventory: object, heading: str) -> None:
    """Render and expose downloads for one set of published SHAP graphics.

    Args:
        database: Canonical DuckDB path inside the completed result.
        inventory: Pandas-like plot inventory frame.
        heading: User-facing section heading.
    """

    if inventory.empty:
        st.info(f"{heading}: no graphic is available for this selection.")
        return
    st.subheader(heading)
    valid_assets: list[tuple[object, Path, bytes, str]] = []
    for _, row in inventory.iterrows():
        try:
            asset_path, payload, mime = _read_result_asset(
                database=database,
                relative_path=str(row["asset_path"]),
                file_format=str(row["file_format"]),
            )
        except InputValidationError as error:
            st.warning(str(error))
            continue
        valid_assets.append((row, asset_path, payload, mime))
    pdf_identities = {
        (str(row["comparison_id"]), str(row["plot_type"]), str(row.get("protein_id") or ""))
        for row, _, _, _ in valid_assets
        if str(row["file_format"]).upper() == "PDF"
    }
    for row, _, payload, _ in valid_assets:
        if str(row["file_format"]).upper() != "PNG":
            continue
        protein = str(row.get("protein_id") or "")
        identity = (str(row["comparison_id"]), str(row["plot_type"]), protein)
        if identity not in pdf_identities:
            st.warning(
                f"{row['plot_type']} {protein or 'global'} was not displayed because its "
                "required PDF companion is absent."
            )
            continue
        caption = str(row["plot_type"]).replace("_", " ").title()
        if protein:
            caption += f" · {protein}"
        st.image(payload, caption=caption, width="stretch")
    for row, asset_path, payload, mime in valid_assets:
        file_format = str(row["file_format"]).upper()
        protein = str(row.get("protein_id") or "global")
        st.download_button(
            label=f"Download {row['plot_type']} {protein} ({file_format})",
            data=payload,
            file_name=asset_path.name,
            mime=mime,
            key=(f"shap-{row['comparison_id']}-{row['plot_type']}-{protein}-{file_format}"),
        )


def _result_asset_path(*, database: Path, relative_path: str) -> Path:
    """Resolve a verified result-relative asset without allowing path escape.

    Args:
        database: Canonical database path within a completed result.
        relative_path: Portable path recorded in a canonical inventory table.

    Returns:
        Existing resolved asset path.

    Raises:
        InputValidationError: If the path is absolute, escapes the result or is absent.
    """

    result_root = Path(database).expanduser().resolve().parent
    relative = Path(str(relative_path).strip())
    if not str(relative_path).strip() or relative.is_absolute():
        raise InputValidationError("A SHAP asset path must be non-empty and result-relative.")
    candidate = (result_root / relative).resolve()
    if result_root not in candidate.parents or not candidate.is_file():
        raise InputValidationError(
            f"Published SHAP asset is absent or outside the completed result: {relative_path!r}."
        )
    return candidate


def _read_result_asset(
    *, database: Path, relative_path: str, file_format: str
) -> tuple[Path, bytes, str]:
    """Read one bounded, type-checked asset inside a completed result.

    Args:
        database: Canonical database path within a completed result.
        relative_path: Portable path recorded in the plot inventory.
        file_format: Declared PNG, SVG or PDF format.

    Returns:
        Resolved path, validated bytes and MIME type.

    Raises:
        InputValidationError: If the format, extension, size or content is unsafe.
    """

    normalised_format = str(file_format).strip().upper()
    if normalised_format not in _ASSET_FORMATS:
        raise InputValidationError(f"Unsupported published plot format: {file_format!r}.")
    suffix, mime = _ASSET_FORMATS[normalised_format]
    path = _result_asset_path(database=database, relative_path=relative_path)
    if path.suffix.casefold() != suffix:
        raise InputValidationError(
            f"Published {normalised_format} asset has an inconsistent extension: {relative_path!r}."
        )
    size = path.stat().st_size
    if size <= 0 or size > _MAX_APP_ASSET_BYTES:
        raise InputValidationError(
            f"Published plot asset size is outside the allowed range: {relative_path!r}."
        )
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise InputValidationError(
            f"Could not read published plot asset {relative_path!r}."
        ) from error
    if normalised_format == "PNG" and not payload.startswith(b"\x89PNG\r\n\x1a\n"):
        raise InputValidationError(
            f"Published PNG asset has an invalid signature: {relative_path!r}."
        )
    if normalised_format == "PDF" and not payload.startswith(b"%PDF-"):
        raise InputValidationError(
            f"Published PDF asset has an invalid signature: {relative_path!r}."
        )
    if normalised_format == "SVG":
        prefix = payload[:4096].decode("utf-8", errors="ignore")
        if not re.search(r"<svg(?:\s|>)", prefix, flags=re.IGNORECASE):
            raise InputValidationError(
                f"Published SVG asset is not recognisable: {relative_path!r}."
            )
        if re.search(r"<script\b|\bon\w+\s*=|javascript:", prefix, flags=re.IGNORECASE):
            raise InputValidationError(
                f"Published SVG asset contains active content: {relative_path!r}."
            )
    return path, payload, mime


def _render_downloadable_table(
    *, frame: object, download_name: str, height: int | None = None
) -> None:
    """Render one table with TSV and on-demand formatted Excel downloads.

    Args:
        frame: Pandas-like tabular value.
        download_name: Stable human-readable export identity.
        height: Optional Streamlit table height in pixels.
    """

    normalised = normalise_dataframe(value=frame)
    display_options: dict[str, object] = {
        "hide_index": True,
        "width": "stretch",
    }
    if height is not None:
        display_options["height"] = height
    stem = safe_download_stem(value=download_name)
    try:
        tsv = dataframe_to_tsv_bytes(frame=normalised)
    except (InputValidationError, PublicationError) as error:
        st.warning(f"Could not prepare table downloads: {error}")
        return
    fingerprint = hashlib.sha256(tsv).hexdigest()
    workbook_key = f"table-xlsx-payload-{stem}"
    if st.session_state.get(workbook_key, (None, None))[0] != fingerprint:
        st.session_state.pop(workbook_key, None)
    st.dataframe(normalised, **display_options)
    st.download_button(
        label="Download table as TSV",
        data=tsv,
        file_name=f"{stem}.tsv",
        mime="text/tab-separated-values",
        key=f"table-tsv-{stem}",
    )
    if st.button("Prepare formatted Excel workbook", key=f"table-xlsx-prepare-{stem}"):
        try:
            workbook = dataframe_to_xlsx_bytes(
                frame=normalised,
                title=str(download_name).replace("_", " ").title(),
            )
            st.session_state[workbook_key] = (fingerprint, workbook)
        except (InputValidationError, PublicationError) as error:
            st.warning(f"Could not prepare formatted Excel workbook: {error}")
    cached = st.session_state.get(workbook_key)
    if cached is not None and cached[0] == fingerprint:
        st.download_button(
            label="Download formatted Excel workbook",
            data=cached[1],
            file_name=f"{stem}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key=f"table-xlsx-{stem}",
        )


def _render_plotly_figure(*, figure: object, download_name: str) -> None:
    """Render an interactive figure with on-demand PDF, PNG and HTML exports.

    Args:
        figure: Plotly-compatible figure.
        download_name: Stable human-readable export identity.
    """

    stem = safe_download_stem(value=download_name)
    signature = (
        hashlib.sha256(figure.to_json().encode("utf-8")).hexdigest()
        if hasattr(figure, "to_json")
        else repr(figure)
    )
    if st.session_state.get(f"plot-signature-{stem}") != signature:
        st.session_state.pop(f"plot-png-{stem}", None)
        st.session_state.pop(f"plot-html-{stem}", None)
        st.session_state[f"plot-signature-{stem}"] = signature
    st.plotly_chart(figure, width="stretch")
    if st.button(
        "Prepare plot image (PNG)",
        key=f"plot-png-prepare-{stem}",
        help="Render an image of the exact visible plot on demand.",
    ):
        try:
            image = bytes(figure.to_image(format="png"))
            if not image.startswith(b"\x89PNG\r\n\x1a\n"):
                raise PublicationError("The plot renderer returned an invalid PNG image.")
            st.session_state[f"plot-png-{stem}"] = image
        except Exception as error:
            st.warning(f"Could not render this image; check the Kaleido/Chrome runtime: {error}")
    if st.session_state.get(f"plot-png-{stem}") is not None:
        st.download_button(
            label="Download plot image (PNG)",
            data=st.session_state[f"plot-png-{stem}"],
            file_name=f"{stem}.png",
            mime="image/png",
            key=f"plot-png-download-{stem}",
        )
    if st.button(
        "Prepare interactive plot (HTML)",
        key=f"plot-html-prepare-{stem}",
        help="Bundle Plotly into a standalone interactive HTML file.",
    ):
        try:
            st.session_state[f"plot-html-{stem}"] = figure.to_html(
                full_html=True,
                include_plotlyjs=True,
            ).encode("utf-8")
        except (AttributeError, RuntimeError, ValueError) as error:
            st.warning(f"Could not bundle this interactive plot: {error}")
    if st.session_state.get(f"plot-html-{stem}") is not None:
        st.download_button(
            label="Download interactive plot (HTML)",
            data=st.session_state[f"plot-html-{stem}"],
            file_name=f"{stem}.html",
            mime="text/html",
            key=f"plot-html-download-{stem}",
        )
    if not st.button(
        "Prepare plot PDF download",
        key=f"plot-pdf-prepare-{stem}",
        help="Generate the vector PDF only when needed.",
    ):
        return
    try:
        pdf = plotly_figure_to_pdf_bytes(figure=figure)
    except PublicationError as error:
        st.warning(str(error))
        return
    st.download_button(
        label="Download plot as PDF",
        data=pdf,
        file_name=f"{stem}.pdf",
        mime="application/pdf",
        key=f"plot-pdf-{stem}",
        on_click="ignore",
    )


def _metric_text(value: object) -> str:
    """Format an optional metric for compact display.

    Args:
        value: Numeric metric or missing value.

    Returns:
        Three-decimal text or ``NA``.
    """

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    if numeric != numeric:
        return "NA"
    return f"{numeric:.3f}"


def _render_proteins(*, database: Path) -> None:
    """Render protein-level labels, features and Pfam assessment states.

    Args:
        database: Verified result database.
    """

    st.title("Protein & Pfam explorer")
    protein_ids = distinct_values(
        database=database, table_name="proteins", column_name="protein_id"
    )
    if not protein_ids:
        st.info("No proteins are present.")
        return
    protein_id = st.selectbox("Protein", protein_ids)
    inventory = query_dataframe(
        database=database,
        sql=(
            "SELECT protein_id, description, sequence_length, sequence_sha256 "
            "FROM proteins WHERE protein_id = ?"
        ),
        parameters=(protein_id,),
    )
    _render_downloadable_table(
        frame=inventory,
        download_name=f"{protein_id}_protein_inventory",
    )
    labels_tab, domain_tab, feature_tab, assessment_tab = st.tabs(
        ("Labels", "Domains", "Positive features", "Feature assessments")
    )
    with labels_tab:
        _render_downloadable_table(
            frame=query_dataframe(
                database=database,
                sql=(
                    "SELECT m.label_id, p.display_name, m.membership_source, m.direct_label_id "
                    "FROM label_memberships m LEFT JOIN profile_labels p USING (label_id) "
                    "WHERE m.protein_id = ? ORDER BY p.level, m.label_id"
                ),
                parameters=(protein_id,),
            ),
            download_name=f"{protein_id}_labels",
        )
    with domain_tab:
        st.caption("No-hit, not-assessed and failed scans are deliberately different states.")
        _render_downloadable_table(
            frame=query_dataframe(
                database=database,
                sql=(
                    "SELECT * FROM domain_assessments WHERE protein_id = ? "
                    "ORDER BY domain_authority"
                ),
                parameters=(protein_id,),
            ),
            download_name=f"{protein_id}_domain_assessments",
        )
        _render_downloadable_table(
            frame=query_dataframe(
                database=database,
                sql=(
                    'SELECT domain_authority, domain_id, domain_name, "start", "end", '
                    "score, e_value FROM domain_hits WHERE protein_id = ? "
                    'ORDER BY "start", "end"'
                ),
                parameters=(protein_id,),
            ),
            download_name=f"{protein_id}_domain_hits",
        )
    with feature_tab:
        _render_downloadable_table(
            frame=query_dataframe(
                database=database,
                sql=(
                    "SELECT * FROM features WHERE protein_id = ? "
                    'ORDER BY feature_type, "start", feature_id'
                ),
                parameters=(protein_id,),
            ),
            download_name=f"{protein_id}_features",
            height=450,
        )
    with assessment_tab:
        st.caption(
            "Assessed absence, unassessed, failed and excluded evidence remain distinct "
            "from positive feature evidence."
        )
        _render_downloadable_table(
            frame=query_dataframe(
                database=database,
                sql=(
                    "SELECT * FROM feature_assessments WHERE protein_id = ? "
                    'ORDER BY feature_type, feature_id, "start"'
                ),
                parameters=(protein_id,),
            ),
            download_name=f"{protein_id}_feature_assessments",
            height=450,
        )


def _render_classes(*, database: Path) -> None:
    """Render the configured hierarchy and observed class/role coverage.

    Args:
        database: Verified result database.
    """

    st.title("Classes & component roles")
    st.caption(
        "Mechanistic class, system class and component role are independent fields; "
        "a substrate receptor is not silently relabelled as a catalytic protein."
    )
    provisional = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM class_labelling_summary "
            "ORDER BY label_type, direct_positive_protein_count DESC, label_id"
        ),
    )
    if not provisional.empty:
        st.warning(
            "These automated labels are evidence-supported proposals for hypothesis "
            "generation. They are not equivalent to human-reviewed biochemical truth."
        )
        st.subheader("Automated evidence-label coverage")
        _render_downloadable_table(
            frame=provisional,
            download_name="automated_evidence_label_coverage",
            height=420,
        )
    coverage = query_dataframe(
        database=database,
        sql=(
            "SELECT p.label_id, p.display_name, p.parent_label_id, p.level, "
            "p.system_class, p.mechanistic_class, p.component_role, p.family, "
            "p.active_site_expected, p.active_site_residue, "
            "count(DISTINCT m.protein_id) AS reviewed_proteins "
            "FROM profile_labels p LEFT JOIN label_memberships m USING (label_id) "
            "GROUP BY ALL ORDER BY p.label_id"
        ),
    )
    populated = coverage[coverage["reviewed_proteins"] > 0]
    left, right = st.columns((3, 2))
    with left:
        st.subheader("Observed hierarchy")
        if populated.empty:
            st.info("No reviewed-positive profile memberships are present.")
        else:
            figure = px.sunburst(
                populated,
                ids="label_id",
                names="display_name",
                parents="parent_label_id",
                values="reviewed_proteins",
                color="reviewed_proteins",
                hover_data=["mechanistic_class", "component_role"],
            )
            _render_plotly_figure(
                figure=figure,
                download_name="class_hierarchy_sunburst",
            )
    with right:
        st.subheader("Role coverage")
        roles = query_dataframe(
            database=database,
            sql=(
                "SELECT p.component_role, count(DISTINCT m.protein_id) AS proteins "
                "FROM label_memberships m JOIN profile_labels p USING (label_id) "
                "WHERE p.component_role <> '' GROUP BY p.component_role "
                "ORDER BY proteins DESC, p.component_role"
            ),
        )
        _render_downloadable_table(
            frame=roles,
            download_name="component_role_coverage",
            height=360,
        )
    with st.expander("Full controlled vocabulary", expanded=False):
        _render_downloadable_table(
            frame=coverage,
            download_name="controlled_profile_vocabulary",
            height=520,
        )


def _render_structures(*, database: Path) -> None:
    """Render structure availability, folds, clusters and alignment evidence.

    Args:
        database: Verified result database.
    """

    st.title("Structures & folds")
    source_counts = query_dataframe(
        database=database,
        sql=(
            "SELECT structure_source, availability_status, analysis_eligibility_status, "
            "fold_evidence_status, count(*) AS models, "
            "round(avg(mean_confidence), 2) AS mean_confidence "
            "FROM structures GROUP BY structure_source, availability_status, "
            "analysis_eligibility_status, fold_evidence_status "
            "ORDER BY models DESC"
        ),
    )
    _render_downloadable_table(
        frame=source_counts,
        download_name="structure_source_coverage",
    )
    if not source_counts.empty and {
        "analysis_eligibility_status",
        "models",
        "fold_evidence_status",
        "structure_source",
        "availability_status",
        "mean_confidence",
    }.issubset(source_counts.columns):
        figure = px.bar(
            source_counts,
            x="analysis_eligibility_status",
            y="models",
            color="fold_evidence_status",
            hover_data=["structure_source", "availability_status", "mean_confidence"],
            labels={"analysis_eligibility_status": "Model eligibility", "models": "Models"},
        )
        figure.update_layout(title="Published model coverage and eligibility", height=400)
        _render_plotly_figure(figure=figure, download_name="structure_eligibility_summary")
    with st.expander("Structure model records and comparison universes", expanded=False):
        _render_downloadable_table(
            frame=query_dataframe(
                database=database,
                sql=(
                    "SELECT * FROM structures ORDER BY analysis_eligibility_status, "
                    "structure_source, protein_id, structure_id LIMIT 5000"
                ),
            ),
            download_name="structure_model_records",
            height=440,
        )
    fold_counts = query_dataframe(
        database=database,
        sql=(
            "SELECT fold_authority, fold_id, fold_name, count(DISTINCT protein_id) AS proteins "
            "FROM structures WHERE fold_id <> '' GROUP BY ALL ORDER BY proteins DESC, fold_id"
        ),
    )
    left, right = st.columns(2)
    with left:
        st.subheader("Fold assignments")
        if fold_counts.empty:
            st.info(
                "No named fold assignments were published. The alignment-derived "
                "clusters alongside are whole-model similarity evidence, not fold names."
            )
        else:
            _render_downloadable_table(
                frame=fold_counts,
                download_name="fold_assignments",
                height=380,
            )
    with right:
        st.subheader("Alignment-derived clusters")
        clusters = query_dataframe(
            database=database,
            sql=(
                "SELECT cluster_id, max(reference_partition) AS reference_partition, "
                "max(reference_member_count) AS discovery_members, "
                "count_if(membership_method = 'VALIDATION_PROJECTION') AS projected_members, "
                "max(member_count) AS total_members, max(edge_count) AS discovery_edges, "
                "max(tm_score_threshold) AS tm_threshold, max(minimum_coverage) AS coverage "
                "FROM structure_clusters GROUP BY cluster_id "
                "ORDER BY total_members DESC, cluster_id"
            ),
        )
        _render_downloadable_table(
            frame=clusters,
            download_name="structural_clusters",
            height=380,
        )
    st.subheader("Pairwise structural alignments")
    comparisons = query_dataframe(
        database=database,
        sql=("SELECT * FROM structure_comparisons ORDER BY tm_score DESC NULLS LAST LIMIT 5000"),
    )
    if not comparisons.empty:
        plot_data = comparisons.dropna(subset=["tm_score", "coverage_a", "coverage_b"]).copy()
        if not plot_data.empty:
            plot_data["minimum_coverage"] = plot_data[["coverage_a", "coverage_b"]].min(axis=1)
            figure = px.scatter(
                plot_data,
                x="minimum_coverage",
                y="tm_score",
                color="comparison_tool",
                hover_data=["protein_a_id", "protein_b_id", "rmsd_angstrom"],
                labels={"minimum_coverage": "Minimum bilateral coverage"},
            )
            _render_plotly_figure(
                figure=figure,
                download_name="structural_alignment_score_coverage",
            )
    _render_downloadable_table(
        frame=comparisons,
        download_name="pairwise_structural_alignments",
        height=480,
    )
    imported = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM imported_structural_group_summaries "
            "ORDER BY group_support_fraction DESC NULLS LAST, cluster_id"
        ),
    )
    if not imported.empty:
        st.subheader("Imported within-group pocket conservation")
        st.caption(
            "Imported US-align/TM-align evidence complements, rather than replaces, "
            "cross-protein Foldseek fold discovery."
        )
        _render_downloadable_table(
            frame=imported,
            download_name="imported_pocket_conservation",
            height=420,
        )


def _render_model_explorer(*, database: Path) -> None:
    """Explore enrichment on exact sequence positions, models and aligned pairs."""

    st.title("Model & alignment explorer")
    st.caption(
        "White: no significant mapped enrichment; blue: weaker significant enrichment; "
        "red: stronger enrichment. A white residue may simply have no localisable feature. "
        "Colour represents target enrichment for one selected comparison, not pLDDT."
    )
    with st.expander("What can be mapped onto a model?", expanded=False):
        st.markdown(
            "Recorded 1-based feature intervals, exact occurrences of enriched published "
            "k-mers, and optional residue annotations can be coloured. A full-length model "
            "must exactly match the published protein sequence and residue numbering. "
            "Whole-model folds or Foldseek clusters do not specify motif positions. "
            "The pair view is a **sequence** alignment; the "
            "published structural comparisons contain scores and coverage, but no "
            "residue-to-residue superposition."
        )
    overview = query_dataframe(
        database=database,
        sql=(
            "SELECT c.comparison_id, c.display_name, c.target_label_ids, "
            "coalesce(count_if(s.status = 'COMPLETE' AND "
            "s.discovery_prevalence_difference > 0 "
            "AND s.discovery_q_value <= 0.05), 0) AS enriched_count, "
            "coalesce(count_if(s.status = 'COMPLETE'), 0) AS complete_count, "
            "coalesce(count_if(s.status = 'INSUFFICIENT_SAMPLE_SIZE'), 0) "
            "AS insufficient_count, "
            "coalesce(count_if(s.status = 'NO_SIGNIFICANT_SIGNATURE'), 0) "
            "AS no_signature_count "
            "FROM comparisons c LEFT JOIN signatures s USING (comparison_id) "
            "GROUP BY c.comparison_id, c.display_name, c.target_label_ids "
            "ORDER BY enriched_count DESC, c.display_name"
        ),
    )
    if overview.empty:
        st.info("This result has no comparisons to inspect.")
        return
    descriptions = {
        str(row.comparison_id): f"{row.display_name} · {int(row.enriched_count):,} enriched"
        for row in overview.itertuples(index=False)
    }
    enriched_comparisons = tuple(
        str(row.comparison_id)
        for row in overview.itertuples(index=False)
        if int(row.enriched_count) > 0
    )
    st.caption(
        f"{len(enriched_comparisons):,} of {len(overview):,} configured comparisons "
        "have significant positive discovery results in this completed run. "
        "All comparison outcomes are available below."
    )
    with st.expander("All comparison outcomes", expanded=False):
        outcomes = overview[
            [
                "comparison_id",
                "display_name",
                "enriched_count",
                "complete_count",
                "insufficient_count",
                "no_signature_count",
            ]
        ].copy()
        outcomes["analysis_status"] = outcomes.apply(
            lambda row: (
                "SIGNIFICANT_POSITIVE"
                if row["enriched_count"]
                else "COMPLETE_NO_POSITIVE_ENRICHMENT"
                if row["complete_count"]
                else "INSUFFICIENT_SAMPLE_SIZE"
                if row["insufficient_count"]
                else "NO_SIGNIFICANT_SIGNATURE"
                if row["no_signature_count"]
                else "NO_PUBLISHED_SIGNATURE"
            ),
            axis=1,
        )
        _render_downloadable_table(
            frame=outcomes,
            download_name="all_comparison_outcomes",
            height=350,
        )
    comparison_choices = _comparison_choices(descriptions=descriptions)
    comparison_id = comparison_choices[
        st.selectbox("Comparison for protein and model mapping", tuple(comparison_choices))
    ]
    selected = overview.loc[overview["comparison_id"] == comparison_id].iloc[0]
    if int(selected["enriched_count"]) == 0:
        st.info("This comparison has no completed positive significant signatures.")
    target_labels = selected["target_label_ids"]
    proteins = _target_model_proteins(
        database=database,
        target_label_ids=str(target_labels) if pd.notna(target_labels) else "",
    )
    if not proteins:
        proteins = distinct_values(
            database=database, table_name="proteins", column_name="protein_id"
        )
        st.info("No target proteins have a packaged model; search any published protein ID.")
    if not proteins:
        st.info("This result has no protein sequences to inspect.")
        return
    st.caption(
        f"{len(proteins):,} suggested proteins are available. Select one below "
        "or enter an exact published protein ID."
    )
    suggested_id = st.selectbox("Protein to inspect", proteins)
    typed_id = st.text_input("Exact protein ID (optional; overrides selection)").strip()
    protein_id = typed_id or suggested_id
    _render_ranked_associations(
        database=database,
        enriched_comparisons=enriched_comparisons,
        descriptions=descriptions,
    )
    protein = query_dataframe(
        database=database,
        sql="SELECT sequence FROM proteins WHERE protein_id = ?",
        parameters=(protein_id,),
    )
    if protein.empty:
        st.warning(f"Protein {protein_id!r} is not in this completed result.")
        return
    sequence = str(protein.iloc[0]["sequence"])
    upload = st.file_uploader(
        "Optional residue annotations (TSV; pockets or other regions)",
        type=["tsv"],
        help=(
            "The file is held in this browser session. It must contain comparison_id, "
            "protein_id, region_type, region_id, start, end, q_value, "
            "prevalence_difference and evidence_source. Positions are 1-based inclusive."
        ),
    )
    st.download_button(
        "Download annotation template",
        data=("\t".join(ANNOTATION_COLUMNS) + "\n").encode("utf-8"),
        file_name="residue_annotations_template.tsv",
        mime="text/tab-separated-values",
    )
    annotations = None
    if upload is not None:
        try:
            annotations = parse_annotation_tsv(payload=upload.getvalue())
        except InputValidationError as error:
            st.warning(str(error))
    mode_options = (
        ("Published signatures", "Uploaded residue annotations")
        if annotations is not None
        else ("Published signatures",)
    )
    mode = st.radio("Enrichment colour source", mode_options, horizontal=True)
    evidence_tier = (
        st.selectbox(
            "Published evidence tier",
            (
                "Discovery within comparison",
                "Validated within comparison",
                "Validated study-wide",
            ),
            help="Stronger tiers require positive held-out enrichment and the selected q-value.",
        )
        if mode == "Published signatures"
        else "Uploaded annotations"
    )
    maximum_features = (
        st.selectbox(
            "Most significant feature memberships to colour",
            (10, 25, 50, 100, 250, 500, 1_000),
            index=2,
            help=(
                "Limits published membership rows before mapping exact occurrences; "
                "avoids covering the whole protein with many weaker overlapping k-mers."
            ),
        )
        if mode == "Published signatures"
        else 50
    )
    rows = _enriched_intervals(
        database=database,
        protein_id=protein_id,
        comparison_id=comparison_id,
        sequence=sequence,
        annotations=annotations if mode == "Uploaded residue annotations" else None,
        evidence_tier=evidence_tier,
        maximum_features=maximum_features,
    )
    limit_caption = (
        f"At most {maximum_features:,} published feature memberships are mapped. "
        if mode == "Published signatures"
        else "Uploaded annotations are mapped separately. "
    )
    st.caption(
        f"Showing positive target enrichment at {evidence_tier.lower()} q ≤ 0.05. "
        f"{limit_caption}"
        "Overlaps use the smallest q-value; colours saturate at q ≤ 10⁻⁸. "
        "Uploaded annotations form a separate source and are never combined "
        "statistically with published q-values."
    )
    mapped = feature_intervals(rows=rows, sequence_length=len(sequence))
    if len(mapped) != len(rows):
        st.warning("Some interval positions are outside this protein sequence and were omitted.")
    valid = significant_intervals(rows=mapped, sequence_length=len(sequence))
    scores, descriptions = enrichment_track(rows=valid, sequence_length=len(sequence))
    st.subheader("Most enriched mapped regions")
    if valid:
        _render_downloadable_table(
            frame=pd.DataFrame(valid),
            download_name=f"{comparison_id}_{protein_id}_mapped_regions_{mode}",
            height=300,
        )
    else:
        st.info("No localisable enriched regions are available for this selection.")
    _render_enrichment_track(
        protein_id=protein_id,
        scores=scores,
        descriptions=descriptions,
        download_name=f"{comparison_id}_{protein_id}_enrichment_map",
    )
    st.subheader("3D Cα model and external views")
    _render_model_for_protein(
        database=database,
        protein_id=protein_id,
        sequence=sequence,
        scores=scores,
        descriptions=descriptions,
        comparison_id=comparison_id,
    )
    st.subheader("Paired sequence alignment")
    _render_pair_alignment(
        database=database,
        protein_id=protein_id,
        comparison_id=comparison_id,
        reference_sequence=sequence,
        reference_scores=scores,
        annotations=annotations if mode == "Uploaded residue annotations" else None,
        evidence_tier=evidence_tier,
        maximum_features=maximum_features,
    )


def _comparison_choices(*, descriptions: dict[str, str]) -> dict[str, str]:
    """Map distinct visible comparison labels back to their published IDs.

    Args:
        descriptions: Display labels keyed by comparison ID.

    Returns:
        Ordered labels mapped to their comparison IDs.

    Raises:
        InputValidationError: If the generated labels are not unique.
    """

    labels = tuple(descriptions.values())
    if len(set(labels)) == len(labels):
        return {label: comparison_id for comparison_id, label in descriptions.items()}
    choices = {
        f"{label} [{comparison_id}]": comparison_id for comparison_id, label in descriptions.items()
    }
    if len(choices) != len(descriptions):
        raise InputValidationError("Comparison display labels must be unique.")
    return choices


def _render_ranked_associations(
    *,
    database: Path,
    enriched_comparisons: tuple[str, ...],
    descriptions: dict[str, str],
) -> None:
    """Render balanced positive associations across selected target classes.

    Args:
        database: Verified result database.
        enriched_comparisons: Comparison IDs with significant positive discovery rows.
        descriptions: Human-readable comparison labels keyed by ID.
    """

    st.subheader("Most significant association results")
    st.caption(
        "This table can compare classes independently of the protein and model selection. "
        "It includes positive discovery enrichment at q ≤ 0.05, balanced by comparison; "
        "validation and study-wide q-values remain visible. The Signature explorer "
        "shows the full set for any one comparison."
    )
    choices = _comparison_choices(descriptions=descriptions)
    available_labels = tuple(
        label for label, comparison_id in choices.items() if comparison_id in enriched_comparisons
    )
    selected_labels = st.multiselect(
        "Comparisons in ranked association table",
        available_labels,
        default=available_labels[:10],
    )
    selected_results = tuple(choices[label] for label in selected_labels)
    if selected_results:
        rows_per_comparison = st.selectbox("Maximum results per comparison", (20, 100, 500, 1_000))
        placeholders = ",".join("?" for _ in selected_results)
        ranked = query_dataframe(
            database=database,
            sql=(
                "SELECT s.*, row_number() OVER (PARTITION BY s.comparison_id "
                "ORDER BY s.discovery_q_value ASC NULLS LAST, "
                "s.discovery_prevalence_difference DESC, s.feature_type, s.feature_id) "
                "AS within_comparison_rank FROM signatures s "
                f"WHERE s.comparison_id IN ({placeholders}) "
                "AND s.status = 'COMPLETE' "
                "AND s.discovery_prevalence_difference > 0 "
                "AND s.discovery_q_value <= 0.05 "
                "QUALIFY within_comparison_rank <= ? "
                "ORDER BY within_comparison_rank, s.comparison_id LIMIT 5000"
            ),
            parameters=(*selected_results, rows_per_comparison),
        )
        _render_downloadable_table(
            frame=ranked,
            download_name="selected_comparisons_top_enriched_signatures",
            height=300,
        )
        st.caption(
            f"Up to {rows_per_comparison:,} rows per comparison and 5,000 total "
            "are shown here. The selected protein's mapped regions below still use "
            "only the protein/model comparison."
        )
        strongest = (
            ranked.dropna(subset=["discovery_q_value"])
            .groupby("comparison_id", sort=False)
            .head(5)
            .copy()
        )
    else:
        st.info("Select one or more comparisons to review association results.")
        strongest = pd.DataFrame()
    if not strongest.empty:
        strongest["comparison"] = strongest["comparison_id"].map(
            lambda value: descriptions[str(value)].split(" · ")[0]
        )
        strongest["label"] = (
            strongest["comparison"].astype(str).str.slice(0, 30)
            + " · "
            + strongest["feature_type"].astype(str)
            + ": "
            + strongest["feature_id"].astype(str).str.slice(0, 45)
        )
        figure = px.bar(
            strongest,
            x="discovery_prevalence_difference",
            y="label",
            color="comparison",
            orientation="h",
            hover_data=[
                "feature_name",
                "feature_type",
                "discovery_q_value",
                "evidence_class",
            ],
            labels={
                "discovery_prevalence_difference": "Target − background prevalence",
                "label": "Enriched feature",
                "comparison": "Comparison",
            },
        )
        figure.update_yaxes(autorange="reversed")
        figure.update_layout(
            title="Top ranked enriched features across selected comparisons",
            height=max(450, 30 * len(strongest)),
        )
        _render_plotly_figure(
            figure=figure,
            download_name="selected_comparisons_top_enriched_features",
        )


def _target_model_proteins(*, database: Path, target_label_ids: str) -> tuple[str, ...]:
    """Suggest model-bearing members of the selected comparison's target class.

    Args:
        database: Verified result database.
        target_label_ids: Pipe-delimited target labels from the comparison table.

    Returns:
        At most 5,000 protein IDs with published coordinates, prioritised by
        recorded mean model confidence.
    """

    labels = tuple(label for label in target_label_ids.split("|") if label)
    if not labels:
        return ()
    placeholders = ",".join("?" for _ in labels)
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT s.protein_id, max(s.mean_confidence) AS model_confidence "
            "FROM structures s JOIN label_memberships l USING (protein_id) "
            f"WHERE l.label_id IN ({placeholders}) AND s.coordinate_path <> '' "
            "GROUP BY s.protein_id ORDER BY model_confidence DESC NULLS LAST, "
            "s.protein_id LIMIT 5000"
        ),
        parameters=labels,
    )
    return tuple(str(value) for value in frame["protein_id"])


def _enriched_intervals(
    *,
    database: Path,
    protein_id: str,
    comparison_id: str,
    sequence: str | None = None,
    annotations: object = None,
    evidence_tier: str = "Discovery within comparison",
    maximum_features: int = 50,
) -> list[dict[str, object]]:
    """Find position-bearing published signatures or exact uploaded rows."""

    if annotations is not None:
        selected = annotations[
            (annotations["protein_id"] == protein_id)
            & (annotations["comparison_id"] == comparison_id)
        ].copy()
        selected = selected.rename(columns={"region_id": "feature_name"})
        return selected.to_dict(orient="records")
    if maximum_features < 1 or maximum_features > 10_000:
        raise InputValidationError("Mapped feature limit must be between 1 and 10,000.")
    columns = {
        "Discovery within comparison": ("discovery_q_value", "discovery_prevalence_difference"),
        "Validated within comparison": ("validation_q_value", "validation_prevalence_difference"),
        "Validated study-wide": ("validation_study_q_value", "validation_prevalence_difference"),
    }
    if evidence_tier not in columns:
        raise InputValidationError("Unknown enrichment evidence tier.")
    q_column, difference_column = columns[evidence_tier]
    accepted_classes = {
        "Discovery within comparison": "",
        "Validated within comparison": (
            "AND s.evidence_class IN ('DECISION_CANDIDATE__VALIDATED_STUDY_WIDE', "
            "'DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON') "
        ),
        "Validated study-wide": (
            "AND s.evidence_class = 'DECISION_CANDIDATE__VALIDATED_STUDY_WIDE' "
        ),
    }
    frame = query_dataframe(
        database=database,
        sql=(
            'SELECT f.feature_type, f.feature_id, f.feature_name, f."start", f."end", '
            "f.evidence_source, s.discovery_q_value, s.discovery_prevalence_difference, "
            "s.validation_q_value, s.validation_study_q_value, "
            "s.validation_prevalence_difference, s.evidence_class "
            "FROM features f JOIN signatures s ON f.feature_type = s.feature_type "
            "AND f.feature_id = s.feature_id "
            "WHERE f.protein_id = ? AND s.comparison_id = ? AND s.status = 'COMPLETE' "
            "AND (f.feature_type = 'AMINO_ACID_KMER' OR "
            '(f."start" IS NOT NULL AND f."end" IS NOT NULL)) '
            f"AND s.{difference_column} > 0 AND s.{q_column} <= 0.05 "
            f"{accepted_classes[evidence_tier]}"
            f"ORDER BY s.{q_column}, f.feature_id LIMIT ?"
        ),
        parameters=(protein_id, comparison_id, maximum_features),
    )
    if evidence_tier != "Discovery within comparison":
        accepted = "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE"
        if evidence_tier == "Validated within comparison":
            accepted = (accepted, "DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON")
        frame = frame[
            frame["evidence_class"].isin(accepted if isinstance(accepted, tuple) else (accepted,))
        ].copy()
    frame["q_value"] = frame[q_column]
    frame["prevalence_difference"] = frame[difference_column]
    rows = frame.to_dict(orient="records")
    return project_kmer_intervals(rows=rows, sequence=sequence) if sequence else rows


def _render_enrichment_track(
    *, protein_id: str, scores: list[float], descriptions: list[str], download_name: str
) -> None:
    """Render a downloadable positional heat map with exact hover evidence."""

    figure = go.Figure(
        go.Heatmap(
            z=[scores],
            x=list(range(1, len(scores) + 1)),
            y=[protein_id],
            text=[descriptions],
            hovertemplate="Position %{x}<br>%{text}<extra></extra>",
            zmin=0,
            zmax=1,
            colorscale=[
                [0, "#ffffff"],
                [0.099, "#ffffff"],
                [0.1, "#2463c5"],
                [1, "#b41821"],
            ],
            colorbar={
                "title": "Enrichment",
                "tickvals": [0, 0.1, 1],
                "ticktext": ["No mapped hit", "q ≤ 0.05", "q ≤ 10⁻⁸"],
            },
        )
    )
    figure.update_layout(
        title="2D enrichment along the protein sequence",
        xaxis_title="1-based residue position",
        height=270,
    )
    _render_plotly_figure(figure=figure, download_name=download_name)


def _render_model_for_protein(
    *,
    database: Path,
    protein_id: str,
    sequence: str,
    scores: list[float],
    descriptions: list[str],
    comparison_id: str,
) -> None:
    """Show a verified model, an opt-in AFDB fetch and its exact external links."""

    models = query_dataframe(
        database=database,
        sql=(
            "SELECT structure_id, structure_source, coordinate_path, coordinate_sha256, "
            "mean_confidence, analysis_eligibility_status FROM structures "
            "WHERE protein_id = ? ORDER BY (coordinate_path <> '') DESC, "
            "(analysis_eligibility_status = 'ELIGIBLE') DESC, "
            "mean_confidence DESC NULLS LAST, structure_id"
        ),
        parameters=(protein_id,),
    )
    acquisitions = query_dataframe(
        database=database,
        sql=(
            "SELECT uniprot_accession FROM alphafold_acquisitions "
            "WHERE protein_id = ? AND sequence_match = true ORDER BY uniprot_accession"
        ),
        parameters=(protein_id,),
    )
    accession = next(
        (
            value
            for raw in acquisitions["uniprot_accession"]
            if (value := canonical_accession(value=raw)) is not None
        ),
        canonical_accession(value=protein_id),
    )
    for label, url in external_links(accession=accession).items():
        st.link_button(f"Open {label}: {accession}", url)
    if accession is None:
        st.caption("External accession-specific links require a canonical UniProt ID.")
    payload = None
    filename = ""
    if not models.empty:
        model_options = tuple(
            f"{row.structure_id} · {row.structure_source} · record {index + 1}"
            for index, row in enumerate(models.itertuples(index=False))
        )
        selected = model_options.index(st.selectbox("Published structure record", model_options))
        row = models.iloc[selected]
        relative = str(row["coordinate_path"] or "")
        if relative:
            try:
                payload = read_published_model(
                    database=database,
                    relative_path=relative,
                    sha256=str(row["coordinate_sha256"]),
                )
                filename = Path(relative).name.removesuffix(".gz")
                if relative.endswith(".gz"):
                    st.caption("Published compressed model was decompressed for download.")
            except InputValidationError as error:
                st.warning(str(error))
        else:
            st.info("This published structure record has no packaged coordinates.")
    session_model_key = (
        f"afdb-{accession}-{hashlib.sha256(sequence.encode('ascii')).hexdigest()}"
        if accession is not None
        else None
    )
    if payload is None and accession is not None and session_model_key is not None:
        st.caption(
            "Fetch a current AlphaFold DB PDB only if you need it; the bundle stays unchanged."
        )
        if st.button("Fetch exact AlphaFold model", key=f"fetch-{protein_id}"):
            try:
                st.session_state[session_model_key] = fetch_alphafold_model(
                    accession=accession,
                    sequence=sequence,
                )
            except InputValidationError as error:
                st.warning(str(error))
        payload = st.session_state.get(session_model_key)
        if payload is not None:
            filename = f"AF-{accession}-verified-session.pdb"
    if payload is None:
        st.info("No local coordinates are available for an in-app 3D view.")
        return
    st.download_button(
        "Download selected coordinate model",
        data=payload,
        file_name=filename,
        mime="chemical/x-pdb" if filename.endswith(".pdb") else "chemical/x-cif",
    )
    st.link_button(
        "Search selected model with Foldseek",
        "https://search.foldseek.com/search",
    )
    try:
        trace = (
            parse_pdb_trace(payload=payload)
            if filename.endswith(".pdb")
            else parse_mmcif_trace(payload=payload)
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    exact = trace.matches_sequence(sequence=sequence)
    if not exact:
        st.warning(
            "This model does not match the full published sequence and residue numbering. "
            "The grey 3D trace is shown without projected enrichment."
        )
    _render_trace_figure(
        trace=trace,
        scores=scores if exact else None,
        descriptions=descriptions if exact else None,
        download_name=f"{comparison_id}_{protein_id}_three_dimensional_model",
    )


def _render_trace_figure(
    *,
    trace: ModelTrace,
    scores: list[float] | None,
    descriptions: list[str] | None,
    download_name: str,
) -> None:
    """Render a rotatable Cα backbone coloured only for exact mapping."""

    residues = trace.residues
    x, y, z = ([getattr(item, axis) for item in residues] for axis in ("x", "y", "z"))
    figure = go.Figure()
    figure.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            line={"color": "#8293a1", "width": 4},
            name="Cα backbone",
            hoverinfo="skip",
        )
    )
    text = (
        [
            f"{item.amino_acid}{item.position} · {descriptions[index]}"
            for index, item in enumerate(residues)
        ]
        if descriptions is not None
        else [f"{item.amino_acid}{item.position}" for item in residues]
    )
    figure.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            name="Residue enrichment",
            text=text,
            hovertemplate="%{text}<extra></extra>",
            marker={
                "size": 4,
                "color": scores if scores is not None else "#8c9aa6",
                "cmin": 0,
                "cmax": 1,
                "colorscale": [
                    [0, "#ffffff"],
                    [0.099, "#ffffff"],
                    [0.1, "#2463c5"],
                    [1, "#b41821"],
                ],
                "line": {"color": "#587086", "width": 0.4},
                "showscale": scores is not None,
            },
        )
    )
    figure.update_layout(
        title=f"Cα trace · chain {trace.chain} · drag to rotate",
        height=700,
        scene={"aspectmode": "data"},
    )
    _render_plotly_figure(figure=figure, download_name=download_name)


def _render_pair_alignment(
    *,
    database: Path,
    protein_id: str,
    comparison_id: str,
    reference_sequence: str,
    reference_scores: list[float],
    annotations: object,
    evidence_tier: str,
    maximum_features: int,
) -> None:
    """Show structurally compared partners with a distinct sequence alignment."""

    candidates = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM structure_comparisons WHERE protein_a_id = ? OR protein_b_id = ? "
            "ORDER BY tm_score DESC NULLS LAST LIMIT 500"
        ),
        parameters=(protein_id, protein_id),
    )
    if candidates.empty:
        st.info("No published pairwise structural comparison includes this protein.")
        return
    pair_options = tuple(
        f"{row.protein_a_id} ↔ {row.protein_b_id} · "
        f"TM {_metric_text(row.tm_score)} · record {index + 1}"
        for index, row in enumerate(candidates.itertuples(index=False))
    )
    selected = pair_options.index(st.selectbox("Published structural comparison", pair_options))
    selected_pair = candidates.iloc[selected]
    other_id = str(
        selected_pair["protein_b_id"]
        if selected_pair["protein_a_id"] == protein_id
        else selected_pair["protein_a_id"]
    )
    _render_downloadable_table(
        frame=candidates.iloc[[selected]],
        download_name=f"{comparison_id}_{protein_id}_{other_id}_structural_comparison",
    )
    st.caption(
        "The structural row gives aggregate TM-score and coverage. The alignment below "
        "is a fresh global amino-acid alignment (+2 match, −1 substitution, −2 gap). "
        "Its columns are not Foldseek superposed-residue coordinates."
    )
    other = query_dataframe(
        database=database,
        sql="SELECT sequence FROM proteins WHERE protein_id = ?",
        parameters=(other_id,),
    )
    if other.empty:
        st.warning("The partner sequence is missing from the published result.")
        return
    other_sequence = str(other.iloc[0]["sequence"])
    try:
        aligned = align_sequences(reference=reference_sequence, comparison=other_sequence)
    except InputValidationError as error:
        st.info(str(error))
        return
    other_rows = _enriched_intervals(
        database=database,
        protein_id=other_id,
        comparison_id=comparison_id,
        sequence=other_sequence,
        annotations=annotations,
        evidence_tier=evidence_tier,
        maximum_features=maximum_features,
    )
    other_scores, _ = enrichment_track(rows=other_rows, sequence_length=len(other_sequence))
    alignment_rows = []
    for index, column in enumerate(aligned, start=1):
        a, b = column.reference_position, column.comparison_position
        alignment_rows.append(
            {
                "alignment_column": index,
                "reference_position": a,
                "reference_residue": column.reference,
                "reference_enrichment": reference_scores[a - 1] if a else 0.0,
                "comparison_position": b,
                "comparison_residue": column.comparison,
                "comparison_enrichment": other_scores[b - 1] if b else 0.0,
                "identity": a is not None
                and b is not None
                and column.reference == column.comparison,
            }
        )
    st.download_button(
        "Download aligned pair FASTA",
        data=alignment_fasta(
            reference_id=protein_id,
            comparison_id=other_id,
            columns=aligned,
        ),
        file_name=f"{safe_download_stem(value=protein_id + '_' + other_id)}_aligned.fasta",
        mime="text/plain",
    )
    window_start = (
        st.slider(
            "Alignment window start",
            min_value=1,
            max_value=len(aligned) - 99,
            value=1,
            help="Scroll through the sequence alignment in windows of 100 columns.",
        )
        if len(aligned) > 100
        else 1
    )
    window = alignment_rows[window_start - 1 : window_start + 99]
    figure = go.Figure(
        go.Heatmap(
            z=[
                [row["reference_enrichment"] for row in window],
                [row["comparison_enrichment"] for row in window],
            ],
            x=[row["alignment_column"] for row in window],
            y=[protein_id, other_id],
            zmin=0,
            zmax=1,
            colorscale=[[0, "#ffffff"], [0.099, "#ffffff"], [0.1, "#2463c5"], [1, "#b41821"]],
            customdata=[
                [
                    f"{row['reference_residue']} {row['reference_position'] or 'gap'}"
                    for row in window
                ],
                [
                    f"{row['comparison_residue']} {row['comparison_position'] or 'gap'}"
                    for row in window
                ],
            ],
            hovertemplate="Column %{x}<br>%{y}: %{customdata}<extra></extra>",
        )
    )
    figure.update_layout(title="Enrichment on aligned sequence columns", height=300)
    _render_plotly_figure(
        figure=figure,
        download_name=f"{comparison_id}_{protein_id}_{other_id}_alignment_window_{window_start}",
    )
    _render_downloadable_table(
        frame=pd.DataFrame(alignment_rows),
        download_name=f"{comparison_id}_{protein_id}_{other_id}_full_sequence_alignment",
        height=330,
    )
    accessions = query_dataframe(
        database=database,
        sql=(
            "SELECT protein_id, uniprot_accession FROM alphafold_acquisitions "
            "WHERE protein_id IN (?, ?) AND sequence_match = true"
        ),
        parameters=(protein_id, other_id),
    )
    mapping = dict(zip(accessions["protein_id"], accessions["uniprot_accession"]))
    for label, url in pair_links(
        reference=mapping.get(protein_id, protein_id),
        comparison=mapping.get(other_id, other_id),
    ).items():
        st.link_button(label, url)


def _render_glossary() -> None:
    """Show a searchable, downloadable dictionary for any protein profile."""

    st.title("Glossary & help")
    term = st.text_input("Search definitions", help="Search terms, categories and explanations.")
    frame = pd.DataFrame(GLOSSARY, columns=["category", "term", "definition"])
    if term.strip():
        mask = frame.apply(
            lambda column: column.astype(str).str.contains(term.strip(), case=False, regex=False)
        ).any(axis=1)
        frame = frame.loc[mask].reset_index(drop=True)
    _render_downloadable_table(frame=frame, download_name="protein_signature_glossary")


def _render_orthology(*, database: Path) -> None:
    """Render OrthoFinder context and leakage-safe data partitions.

    Args:
        database: Verified result database.
    """

    st.title("Orthology & partitions")
    st.caption(
        "The composite OrthoFinder authority is run ID + group type + hierarchy node + "
        "group ID. Whole connected homology/redundancy blocks stay in one partition."
    )
    columns = st.columns(4)
    columns[0].metric(
        "Memberships",
        f"{table_count(database=database, table_name='orthofinder_memberships'):,}",
    )
    groups = query_dataframe(
        database=database,
        sql=(
            "SELECT count(DISTINCT (run_id, group_type, hierarchy_node, group_id)) AS n "
            "FROM orthofinder_memberships"
        ),
    )
    columns[1].metric("Groups", f"{int(groups.iloc[0]['n']):,}")
    blocks = query_dataframe(
        database=database,
        sql="SELECT count(DISTINCT partition_key) AS n FROM partitions",
    )
    columns[2].metric("Partition blocks", f"{int(blocks.iloc[0]['n']):,}")
    near = query_dataframe(
        database=database,
        sql=(
            "SELECT count(DISTINCT cluster_id) AS n FROM redundancy_clusters "
            "WHERE cluster_type = 'NEAR_REDUNDANCY'"
        ),
    )
    columns[3].metric("Near-redundancy clusters", f"{int(near.iloc[0]['n']):,}")
    partitions = query_dataframe(
        database=database,
        sql=(
            "SELECT partition, partition_unit, count(*) AS proteins, "
            "count(DISTINCT partition_key) AS blocks FROM partitions "
            "GROUP BY ALL ORDER BY partition, partition_unit"
        ),
    )
    st.subheader("Discovery/validation allocation")
    _render_downloadable_table(
        frame=partitions,
        download_name="discovery_validation_allocation",
    )
    context = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM orthofinder_group_context "
            "ORDER BY member_count DESC, run_id, group_type, hierarchy_node, group_id "
            "LIMIT 5000"
        ),
    )
    st.subheader("Published OrthoFinder group context")
    if context.empty:
        st.info(
            "No OrthoFinder group context was supplied. Exact-sequence blocks still "
            "prevent duplicate-sequence leakage."
        )
    else:
        _render_downloadable_table(
            frame=context,
            download_name="orthofinder_group_context",
            height=500,
        )


def _render_quality(*, database: Path, metadata: dict[str, object]) -> None:
    """Render provenance, assessment coverage and controlled failure states.

    Args:
        database: Verified result database.
        metadata: Run metadata.
    """

    st.title("Data quality & provenance")
    st.success("Completion marker and all checksums verified when this resource was opened.")
    label_evidence = metadata.get("automated_label_evidence", {})
    if isinstance(label_evidence, dict) and label_evidence.get("status") != "NOT_SELECTED":
        st.subheader("Automated label evidence and circularity safeguards")
        warning = str(label_evidence.get("warning") or "").strip()
        if warning:
            st.warning(warning)
        evidence_tabs = st.tabs(
            ("Decisions", "Matched controls", "Excluded label features", "Abstentions")
        )
        evidence_queries = (
            (
                "SELECT * FROM label_evidence_audit ORDER BY protein_id, label_id, rule_id",
                "label_evidence_audit",
            ),
            (
                "SELECT * FROM control_matching_audit "
                "ORDER BY background_label_id, target_unit_id, control_unit_id",
                "control_matching_audit",
            ),
            (
                "SELECT * FROM label_definition_features "
                "ORDER BY label_id, feature_type, feature_id",
                "label_definition_features",
            ),
            (
                "SELECT * FROM unresolved_assignments ORDER BY curation_status, protein_id",
                "unresolved_assignments",
            ),
        )
        for tab, (sql, download_name) in zip(
            evidence_tabs,
            evidence_queries,
            strict=True,
        ):
            with tab:
                _render_downloadable_table(
                    frame=query_dataframe(database=database, sql=sql),
                    download_name=download_name,
                    height=420,
                )
    st.subheader("Feature assessment coverage")
    st.caption(
        "A missing positive row is not treated as absence: explicit assessment state and "
        "derivation scope determine each feature's tested universe."
    )
    _render_downloadable_table(
        frame=query_dataframe(
            database=database,
            sql=(
                "SELECT feature_type, evidence_status, derivation_scope, "
                "count(DISTINCT feature_id) AS features, "
                "count(DISTINCT protein_id) AS proteins FROM feature_assessments "
                "GROUP BY ALL ORDER BY feature_type, evidence_status, derivation_scope"
            ),
        ),
        download_name="feature_assessment_coverage",
    )
    st.subheader("Domain assessment coverage")
    _render_downloadable_table(
        frame=query_dataframe(
            database=database,
            sql=(
                "SELECT domain_authority, assessment_status, count(*) AS proteins "
                "FROM domain_assessments GROUP BY ALL ORDER BY domain_authority, assessment_status"
            ),
        ),
        download_name="domain_assessment_coverage",
    )
    st.subheader("AlphaFold acquisition outcomes")
    _render_downloadable_table(
        frame=query_dataframe(
            database=database,
            sql=(
                "SELECT acquisition_status, count(*) AS proteins, "
                "round(avg(mean_plddt), 2) AS mean_plddt "
                "FROM alphafold_acquisitions GROUP BY acquisition_status "
                "ORDER BY proteins DESC"
            ),
        ),
        download_name="alphafold_acquisition_outcomes",
    )
    st.subheader("Curation states")
    _render_downloadable_table(
        frame=query_dataframe(
            database=database,
            sql=(
                "SELECT curation_status, count(DISTINCT protein_id) AS proteins, "
                "count(*) AS assignments FROM label_assignments "
                "GROUP BY curation_status ORDER BY proteins DESC, curation_status"
            ),
        ),
        download_name="curation_state_coverage",
    )
    st.subheader("Inferential independence blocks")
    st.caption(
        "Blocks containing both a target and background member are excluded from that "
        "comparison rather than counted in both classes."
    )
    _render_downloadable_table(
        frame=query_dataframe(
            database=database,
            sql=(
                "SELECT comparison_id, partition, max(target_unit_count) AS target_blocks, "
                "max(background_unit_count) AS background_blocks, "
                "max(excluded_mixed_unit_count) AS excluded_mixed_blocks "
                "FROM associations GROUP BY comparison_id, partition "
                "ORDER BY comparison_id, partition"
            ),
        ),
        download_name="inferential_independence_blocks",
        height=400,
    )
    with st.expander("Run metadata", expanded=False):
        st.json(metadata)


if __name__ == "__main__":  # pragma: no cover
    main()
