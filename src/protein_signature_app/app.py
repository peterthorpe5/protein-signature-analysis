"""Streamlit user interface for sequence, domain, fold and structure evidence."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import plotly.express as px
import streamlit as st

from protein_signature_app.backend import (
    canonical_table_names,
    canonical_table_preview,
    distinct_values,
    load_canonical_table_assets,
    load_metadata,
    load_report_inventory_assets,
    query_dataframe,
    resolve_database,
    table_count,
)
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
        database = resolve_database(resource=arguments.resource)
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
            "Orthology & partitions",
            "Canonical data & downloads",
            "Data quality & provenance",
        ),
    )
    if page == "Overview":
        _render_overview(database=database, metadata=metadata)
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
    elif page == "Orthology & partitions":
        _render_orthology(database=database)
    elif page == "Canonical data & downloads":
        _render_canonical_data(database=database)
    else:
        _render_quality(database=database, metadata=metadata)


def _render_canonical_data(*, database: Path) -> None:
    """Render every canonical dataset with complete published downloads.

    Args:
        database: Verified result database.
    """

    st.title("Canonical data & downloads")
    names = canonical_table_names()
    st.caption(
        "Browse a bounded preview, then download the complete checksum-verified "
        f"TSV or formatted Excel workbook. All {len(names)} canonical datasets are available."
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
    st.dataframe(preview, use_container_width=True, hide_index=True)
    try:
        assets = load_canonical_table_assets(database=database, table_name=table_name)
    except InputValidationError as error:
        st.error(f"Could not prepare canonical downloads: {error}")
        return
    for asset in assets:
        st.download_button(
            label=f"Download complete {asset.file_format}",
            data=asset.payload,
            file_name=Path(asset.relative_path).name,
            mime=asset.mime_type,
            key=f"canonical_{table_name}_{asset.file_format}",
        )
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


def _render_overview(*, database: Path, metadata: dict[str, object]) -> None:
    """Render campaign metrics and evidence composition.

    Args:
        database: Verified result database.
        metadata: Run metadata.
    """

    st.title("Protein signature analysis")
    st.caption("Sequence · domains · folds · pairwise structural evidence")
    st.info(
        "Signatures are prioritisation evidence, not proof of biochemical activity. "
        "Discovery and held-out validation are kept separate."
    )
    columns = st.columns(4)
    columns[0].metric("Proteins", f"{table_count(database=database, table_name='proteins'):,}")
    columns[1].metric(
        "Candidate signatures", f"{table_count(database=database, table_name='signatures'):,}"
    )
    columns[2].metric(
        "Structure models", f"{table_count(database=database, table_name='structures'):,}"
    )
    columns[3].metric(
        "Domain hits", f"{table_count(database=database, table_name='domain_hits'):,}"
    )
    feature_counts = query_dataframe(
        database=database,
        sql=(
            "SELECT feature_type, count(DISTINCT feature_id) AS feature_count, "
            "count(DISTINCT protein_id) AS protein_count FROM features "
            "GROUP BY feature_type ORDER BY protein_count DESC, feature_type"
        ),
    )
    left, right = st.columns((3, 2))
    with left:
        st.subheader("Evidence coverage")
        if feature_counts.empty:
            st.info("No positive feature evidence was available.")
        else:
            figure = px.bar(
                feature_counts,
                x="feature_type",
                y="protein_count",
                color="feature_count",
                labels={"protein_count": "Proteins", "feature_type": "Feature type"},
            )
            _render_plotly_figure(figure=figure, download_name="overview_evidence_coverage")
    with right:
        st.subheader("Availability states")
        availability = metadata.get("evidence_availability", {})
        if isinstance(availability, dict):
            _render_downloadable_table(
                frame={"evidence": list(availability), "status": list(availability.values())},
                download_name="overview_availability_states",
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
    st.caption(
        "Target prevalence describes how common a feature is. The q-value controls FDR "
        "within this comparison and evidence family; study q-value also corrects across "
        "all configured comparisons in the same evidence family."
    )
    selected_comparison = st.selectbox("Comparison", comparisons)
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
    """Render one table with matching TSV and formatted Excel downloads.

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
        workbook = dataframe_to_xlsx_bytes(
            frame=normalised,
            title=str(download_name).replace("_", " ").title(),
        )
    except (InputValidationError, PublicationError) as error:
        st.warning(f"Could not prepare table downloads: {error}")
        return
    st.dataframe(normalised, **display_options)
    st.download_button(
        label="Download table as TSV",
        data=tsv,
        file_name=f"{stem}.tsv",
        mime="text/tab-separated-values",
        key=f"table-tsv-{stem}",
    )
    st.download_button(
        label="Download formatted Excel workbook",
        data=workbook,
        file_name=f"{stem}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key=f"table-xlsx-{stem}",
    )


def _render_plotly_figure(*, figure: object, download_name: str) -> None:
    """Render an interactive Plotly figure with an on-demand PDF download.

    Args:
        figure: Plotly-compatible figure.
        download_name: Stable human-readable export identity.
    """

    stem = safe_download_stem(value=download_name)
    st.plotly_chart(figure, width="stretch")
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
