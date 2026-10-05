"""Streamlit user interface for sequence, domain, fold and structure evidence."""

from __future__ import annotations

import argparse
import hashlib
import html
import logging
import re
from pathlib import Path

import numpy as np
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
from protein_signature_app.help_content import (
    METRIC_HELP,
    PAGE_METHODS,
    PAGE_TERMS,
    campaign_limit_rows,
    feature_explanation,
    glossary_rows,
    graph_explanation,
    metric_reading,
    status_explanation,
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
from protein_signatures.result_help import column_definition

LOGGER = logging.getLogger(__name__)

_MAX_APP_ASSET_BYTES = 100 * 1024 * 1024
_SIGNATURE_Q_PLOT_FLOOR = 1e-300
_ASSET_FORMATS = {
    "PNG": (".png", "image/png"),
    "SVG": (".svg", "image/svg+xml"),
    "PDF": (".pdf", "application/pdf"),
}
_AMINO_ACIDS = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
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
    st.session_state["help-result-metadata"] = metadata
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
        st.caption("Open the question-mark panels on the page and hover over table headings.")
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
        _render_classes(database=database, metadata=metadata)
    elif page == "Structures & folds":
        _render_structures(database=database)
    elif page == "Model & alignment explorer":
        _render_model_explorer(database=database)
    elif page == "Orthology & partitions":
        _render_orthology(database=database, metadata=metadata)
    elif page == "Canonical data & downloads":
        _render_canonical_data(database=database)
    elif page == "Data quality & provenance":
        _render_quality(database=database, metadata=metadata)
    else:
        _render_glossary()


def _feature_explanation(*, feature_type: str, feature_id: str, feature_name: str) -> str:
    """Explain a canonical feature without replacing its stable identifier.

    Args:
        feature_type: Published evidence family.
        feature_id: Machine-stable feature key.
        feature_name: Published short name, if available.

    Returns:
        Readable biological meaning and a key scope limitation.
    """

    if feature_type == "AMINO_ACID_KMER":
        match = re.fullmatch(r"k(\d+):([A-Z]+)", feature_id)
        if match and len(match.group(2)) == int(match.group(1)):
            sequence = match.group(2)
            expanded = "–".join(_AMINO_ACIDS.get(letter, letter) for letter in sequence)
            return (
                f"Exact {len(sequence)}-residue sequence {sequence} ({expanded}) "
                "somewhere in a protein; this alone does not establish a functional motif."
            )
        return "Exact short amino-acid sequence present somewhere in a protein."
    if feature_type == "STRUCTURE_CLUSTER":
        return (
            "Whole-model structural-similarity cluster built from aligned model "
            "relationships; the SC_ identifier is a stable fingerprint, not a named "
            "fold, local motif or residue interval."
        )
    if feature_type in {"PFAM_DOMAIN", "DOMAIN"}:
        return (
            f"Recognised sequence domain ({feature_name or feature_id}); inspect domain "
            "hits for coordinates and assessment status."
        )
    if feature_type in {"PFAM_ARCHITECTURE", "DOMAIN_ARCHITECTURE"}:
        return (
            "Ordered list of detected domains along one protein; identifiers after the "
            "colon are domain accessions, not amino-acid positions."
        )
    if feature_type in {"FOLD", "STRUCTURE_FOLD"}:
        return f"Assigned whole-model fold ({feature_name or feature_id}); no motif interval."
    return f"{feature_name or feature_id}: consult the feature family and source record."


def _feature_plot_label(*, feature_type: str, feature_id: str, feature_name: str) -> str:
    """Give a chart feature a short human-readable label.

    Args:
        feature_type: Published evidence family.
        feature_id: Reproducible feature identity.
        feature_name: Published short name.

    Returns:
        Compact label retaining the source identity when useful.
    """

    if feature_type == "AMINO_ACID_KMER":
        match = re.fullmatch(r"k(\d+):([A-Z]+)", feature_id)
        if match:
            expanded = "–".join(_AMINO_ACIDS.get(letter, letter) for letter in match.group(2))
            return f"{match.group(2)} · {expanded} (k={match.group(1)})"
    if feature_type == "STRUCTURE_CLUSTER":
        return f"Whole-model cluster · {feature_id}"
    return f"{feature_name or feature_id} · {feature_id}"[:75]


def _explain_feature_rows(*, frame: pd.DataFrame) -> pd.DataFrame:
    """Add readable meaning to a bounded signature table.

    Args:
        frame: Signature rows with canonical type, ID and optional short name.

    Returns:
        A copy with an explanatory column; canonical columns are preserved.
    """

    explained = frame.copy()
    if {"feature_type", "feature_id"}.issubset(explained.columns):
        explained.insert(
            min(3, len(explained.columns)),
            "feature_explanation",
            [
                _feature_explanation(
                    feature_type=str(row.feature_type),
                    feature_id=str(row.feature_id),
                    feature_name=(
                        str(row.feature_name)
                        if hasattr(row, "feature_name") and pd.notna(row.feature_name)
                        else ""
                    ),
                )
                for row in explained.itertuples(index=False)
            ],
        )
    return explained


def _feature_key() -> None:
    """Explain opaque IDs beside the tables where readers encounter them."""

    with st.expander("❔ What do these feature IDs mean?"):
        st.markdown(
            "**`k3:LPD`** means the three consecutive amino acids L–P–D "
            "(leucine–proline–aspartate) occur somewhere in a protein. `k5:` uses "
            "five amino acids. Prevalence is measured per protein, not per occurrence. "
            "A frequent short word is not automatically a functional motif.\n\n"
            "**`SC_...`** is a stable fingerprint for a whole-model similarity "
            "cluster derived from structural comparisons. It is not a named fold and "
            "does not identify a local 3D motif or pocket. Inspect its members under "
            "Structures & folds.\n\n"
            "**`Pfam:PF...`** names a detected domain; an architecture lists domains "
            "in sequence order with `>`. The raw ID remains in exports so a result "
            "can be traced to the original analysis."
        )


def _render_page_help(*, page: str) -> None:
    """Show contextual methods, terms and actual campaign decision limits.

    Args:
        page: Exact application navigation label.

    Raises:
        InputValidationError: If the page has no registered help content.
    """
    if page not in PAGE_METHODS or page not in PAGE_TERMS:
        raise InputValidationError(f"No contextual help is registered for page {page!r}.")
    with st.expander(label="❓ What can I answer on this page?", expanded=False):
        st.markdown(body=PAGE_HELP[page])
    with st.expander(label="❓ Methods, evidence and limitations", expanded=False):
        st.markdown(body=PAGE_METHODS[page])
    definitions = {term: definition for _, term, definition in glossary_rows(base_rows=GLOSSARY)}
    with st.expander(label="❓ Terms used on this page", expanded=False):
        for term in PAGE_TERMS[page]:
            st.markdown(body=f"**{term}** — {definitions[term]}")
        st.caption(body="Every displayed field is also defined in Glossary & help.")
    with st.expander(label="❓ Sample-size limits and thresholds for this result", expanded=False):
        for row in campaign_limit_rows(metadata=st.session_state.get("help-result-metadata", {})):
            st.markdown(body=f"**{row['meaning']}: {row['limit']}** — {row['unit']}.")
            st.caption(body=row["setting"])
        st.markdown(
            body=(
                "Association limits count **independent blocks**, despite the configuration "
                "names ending in proteins. Model fitting and held-out evaluation require "
                "both the **protein** and **group** minimum in each class. These are "
                "execution limits, not a statistical power calculation. Historical missing "
                "settings are never replaced with current defaults. Positional colouring "
                "and positive-discovery suggestions use the viewer's stated q ≤ 0.05 rule."
            )
        )


def _table_column_configuration(*, frame: pd.DataFrame) -> dict[str, object]:
    """Build header tooltips and preserve small probabilities and q-values.

    Args:
        frame: Normalised visible table.

    Returns:
        Streamlit column specifications with readable labels and exact-name help.
    """
    configuration = {}
    for name in frame.columns:
        label = str(name).replace("_", " ")
        help_text = f"{name}: {column_definition(column_name=name)}"
        if pd.api.types.is_numeric_dtype(frame[name]) and not pd.api.types.is_bool_dtype(
            frame[name]
        ):
            scientific = name.endswith(("q_value", "p_value")) or name == "e_value"
            configuration[name] = st.column_config.NumberColumn(
                label=label,
                help=help_text,
                format="%.3g"
                if scientific
                else "%.6f"
                if name == "predicted_probability"
                else None,
            )
        else:
            configuration[name] = st.column_config.Column(label=label, help=help_text)
    return configuration


def _render_table_help(*, frame: pd.DataFrame) -> None:
    """Explain visible fields and outcome codes next to their table.

    Args:
        frame: Bounded visible table, never the full unqueried result.
    """
    with st.expander(label="❓ Table fields, units and status codes", expanded=False):
        for name in frame.columns:
            st.markdown(body=f"**{name}** — {column_definition(column_name=name)}")
        context = (
            "model" if {"model_type", "cv_roc_auc"}.intersection(frame.columns) else "association"
        )
        for name in (
            "status",
            "analysis_status",
            "assessment_status",
            "evidence_status",
            "curation_status",
            "analysis_eligibility_status",
            "fold_evidence_status",
            "availability_status",
            "acquisition_status",
            "decision",
            "derivation_scope",
            "membership_source",
            "membership_method",
            "evidence_class",
        ):
            if name not in frame.columns:
                continue
            for code in sorted(str(value) for value in frame[name].dropna().unique() if str(value)):
                text = status_explanation(
                    status=code,
                    metadata=st.session_state.get("help-result-metadata", {}),
                    context=context,
                )
                st.markdown(body=f"**{code}** — {text}")
        st.caption(body="Header labels are readable; downloads retain the exact field names above.")


def _render_graph_help(*, graph_name: str, explanation: str = "") -> None:
    """Place interpretation help beside an interactive or published graphic.

    Args:
        graph_name: Stable chart identity or published plot type.
        explanation: Optional plot-specific interpretation from the source renderer.
    """
    with st.expander(
        label="❓ What does this graph show and how should I interpret it?", expanded=False
    ):
        st.markdown(body=explanation or graph_explanation(graph_name=graph_name))


def _render_canonical_data(*, database: Path) -> None:
    """Render every canonical dataset with complete published downloads.

    Args:
        database: Verified result database.
    """

    st.title("Canonical data & downloads")
    _render_page_help(page="Canonical data & downloads")
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


def _comparison_outcomes(*, database: Path) -> pd.DataFrame:
    """Summarise every configured comparison, including those without a discovery.

    Args:
        database: Verified result database.

    Returns:
        One row per comparison with separate discovery and validation counts.
    """

    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT c.comparison_id, c.display_name, c.target_label_ids, "
            "coalesce(count_if(s.status = 'COMPLETE' AND "
            "s.discovery_prevalence_difference > 0 "
            "AND s.discovery_q_value <= 0.05 AND "
            "starts_with(s.evidence_class, 'DECISION_CANDIDATE__')), 0) "
            "AS enriched_count, "
            "coalesce(count_if(s.status = 'COMPLETE' AND "
            "s.discovery_prevalence_difference > 0 "
            "AND s.discovery_q_value <= 0.05 AND "
            "s.evidence_class IN ('DECISION_CANDIDATE__VALIDATED_STUDY_WIDE', "
            "'DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON')), 0) "
            "AS validated_within_count, "
            "coalesce(count_if(s.status = 'COMPLETE' AND "
            "s.discovery_prevalence_difference > 0 "
            "AND s.discovery_q_value <= 0.05 AND "
            "s.evidence_class = 'DECISION_CANDIDATE__VALIDATED_STUDY_WIDE'), 0) "
            "AS validated_study_count, "
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
    if frame.empty:
        return frame
    frame["analysis_status"] = frame.apply(
        lambda row: (
            "DISCOVERY_ENRICHED"
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
    return frame


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
    _render_page_help(page="Overview")
    st.caption("Sequence · domains · folds · pairwise structural evidence")
    st.info(
        "Signatures are prioritisation evidence, not proof of biochemical activity. "
        "Discovery and held-out validation are kept separate."
    )
    columns = st.columns(4)
    columns[0].metric(
        "Proteins",
        f"{_overview_count(database=database, metadata=metadata, table_name='proteins'):,}",
        help="Distinct published protein records; statistical independence is counted by blocks.",
    )
    columns[1].metric(
        "Signature outcome rows",
        f"{_overview_count(database=database, metadata=metadata, table_name='signatures'):,}",
        help="All signature rows, including explicit insufficient/no-signature outcomes. "
        "Completed evidence is shown separately below; this is not a unique-feature count.",
    )
    columns[2].metric(
        "Structure models",
        f"{_overview_count(database=database, metadata=metadata, table_name='structures'):,}",
        help=(
            "Coordinate records; inspect eligibility and confidence before interpreting similarity."
        ),
    )
    columns[3].metric(
        "Domain hits",
        f"{_overview_count(database=database, metadata=metadata, table_name='domain_hits'):,}",
        help="Retained domain-hit intervals; several hits can belong to the same protein.",
    )
    st.subheader("Where can this run support a discovery?")
    outcomes = _comparison_outcomes(database=database)
    if outcomes.empty:
        st.info("This result has no configured comparisons.")
    else:
        status_counts = outcomes["analysis_status"].value_counts()
        summary = st.columns(3)
        summary[0].metric("Comparisons", f"{len(outcomes):,}")
        summary[1].metric(
            "Discovery enriched",
            f"{int(status_counts.get('DISCOVERY_ENRICHED', 0)):,}",
        )
        summary[2].metric(
            "Insufficient sample size",
            f"{int(status_counts.get('INSUFFICIENT_SAMPLE_SIZE', 0)):,}",
        )
        st.caption(
            "An insufficient comparison cannot establish absence of a class signature. "
            "Discovery counts include sequence and whole-model features; held-out "
            "validation is counted separately below."
        )
        enriched = outcomes[outcomes["enriched_count"] > 0].head(12).copy()
        if not enriched.empty:
            st.markdown("**Classes with positive discovery results**")
            visible = enriched[
                [
                    "display_name",
                    "enriched_count",
                    "validated_within_count",
                    "validated_study_count",
                ]
            ]
            _render_downloadable_table(frame=visible, download_name="overview_discovery_shortlist")
            figure = px.bar(
                enriched.sort_values("enriched_count"),
                x="enriched_count",
                y="display_name",
                orientation="h",
                hover_data=["validated_within_count", "validated_study_count"],
                labels={
                    "enriched_count": "Positive discovery signatures",
                    "display_name": "Target versus background",
                },
            )
            figure.update_layout(height=max(300, 65 * len(enriched)))
            _render_plotly_figure(
                figure=figure,
                download_name="overview_discovery_shortlist_chart",
                explanation=(
                    "Each bar counts positive discovery signatures for one target versus "
                    "background comparison. Hover for held-out validation counts. "
                    "A larger bar can contain many correlated short sequence words; "
                    "it is not a count of independent motifs."
                ),
            )
        with st.expander("All comparison outcomes and validation counts"):
            _render_downloadable_table(
                frame=outcomes.drop(columns=["target_label_ids"]),
                download_name="overview_all_comparison_outcomes",
                height=370,
            )
    label_evidence = metadata.get("automated_label_evidence", {})
    if isinstance(label_evidence, dict) and label_evidence.get("human_review_completed") is False:
        st.warning(
            "Target labels in this campaign have not completed human scientific review. "
            "Treat these enrichments as provisional hypotheses."
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
            _render_plotly_figure(
                figure=figure,
                download_name="overview_signature_types",
                explanation=(
                    "Bars count published completed signature rows in each feature "
                    "family across comparisons. Sequence words, domains and "
                    "whole-model groups are different biological scales, so compare "
                    "family composition rather than treating rows as unique motifs."
                ),
            )
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
                    figure=figure,
                    download_name="overview_exact_feature_coverage_chart",
                    explanation=(
                        "The vertical axis counts distinct proteins with a positive "
                        "feature membership by family; colour indicates the number "
                        "of distinct features. This is input coverage, not evidence "
                        "of significant class enrichment."
                    ),
                )


def _signature_evidence(*, frame: pd.DataFrame, tier: str) -> pd.DataFrame:
    """Select positive signatures meeting an explicit evidence threshold.

    Args:
        frame: Published signatures for one comparison and selected feature types.
        tier: Discovery, within-comparison validation or study-wide validation.

    Returns:
        Rows ordered by the tier-specific q-value and prevalence difference.

    Raises:
        InputValidationError: If the requested tier is unsupported.
    """

    columns = {
        "Discovery": ("discovery_q_value", "discovery_prevalence_difference"),
        "Validated within comparison": (
            "validation_q_value",
            "validation_prevalence_difference",
        ),
        "Validated study-wide": (
            "validation_study_q_value",
            "validation_prevalence_difference",
        ),
    }
    if tier not in columns:
        raise InputValidationError(f"Unknown signature evidence tier: {tier}")
    q_column, difference_column = columns[tier]
    if not {q_column, difference_column}.issubset(frame.columns):
        return frame.iloc[:0].copy()
    q_values = pd.to_numeric(frame[q_column], errors="coerce")
    differences = pd.to_numeric(frame[difference_column], errors="coerce")
    accepted = (q_values <= 0.05) & (differences > 0)
    if "status" in frame:
        accepted &= frame["status"].eq("COMPLETE")
    if "evidence_class" in frame:
        accepted &= frame["evidence_class"].astype(str).str.startswith("DECISION_CANDIDATE__")
    if "evidence_class" in frame and tier != "Discovery":
        classes = ("DECISION_CANDIDATE__VALIDATED_STUDY_WIDE",)
        if tier == "Validated within comparison":
            classes += ("DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON",)
        accepted &= frame["evidence_class"].isin(classes)
    return frame.loc[accepted].sort_values(
        [q_column, difference_column, "feature_type", "feature_id"],
        ascending=[True, False, True, True],
    )


def _signature_plot_data(*, frame: pd.DataFrame) -> pd.DataFrame:
    """Prepare finite signature effects and explicitly capped q-value coordinates.

    Args:
        frame: Published signature rows with discovery q-values and effects.

    Returns:
        Detached plottable rows, retaining recorded values and adding display
        coordinates. Recorded zeros and q-values below 1e-300 are flagged.

    Raises:
        InputValidationError: If the input is not a dataframe or lacks fields.
    """
    required = {"discovery_q_value", "discovery_prevalence_difference"}
    if not isinstance(frame, pd.DataFrame) or not required.issubset(frame.columns):
        raise InputValidationError("Signature plotting requires recorded q-values and effects.")
    q_values = pd.to_numeric(arg=frame["discovery_q_value"], errors="coerce").astype("float64")
    effects = pd.to_numeric(arg=frame["discovery_prevalence_difference"], errors="coerce").astype(
        "float64"
    )
    valid = (
        np.isfinite(q_values)
        & q_values.between(left=0, right=1)
        & np.isfinite(effects)
        & effects.between(left=-1, right=1)
    )
    omitted = len(frame) - int(valid.sum())
    if omitted:
        LOGGER.debug(
            "Omitted %d unavailable or invalid signature coordinates from the plot.", omitted
        )
    result = frame.loc[valid].copy()
    selected_q = q_values.loc[valid]
    result["plot_q_value_clipped"] = selected_q < _SIGNATURE_Q_PLOT_FLOOR
    result["plot_negative_log10_q_value"] = -np.log10(
        selected_q.clip(lower=_SIGNATURE_Q_PLOT_FLOOR)
    )
    return result


def _render_signatures(*, database: Path) -> None:
    """Render filterable discovery and validation signature evidence.

    Args:
        database: Verified result database.
    """

    st.title("Signature explorer")
    _render_page_help(page="Signature explorer")
    _feature_key()
    outcomes = _comparison_outcomes(database=database)
    comparisons = tuple(outcomes["comparison_id"].astype(str)) if not outcomes.empty else ()
    feature_types = distinct_values(
        database=database, table_name="signatures", column_name="feature_type"
    )
    if not comparisons:
        st.info("No comparisons are present in this result.")
        return
    counts = dict(zip(outcomes["comparison_id"], outcomes["enriched_count"], strict=True))
    display_names = dict(zip(outcomes["comparison_id"], outcomes["display_name"], strict=True))
    st.caption(
        "Target prevalence describes how common a feature is. The q-value controls FDR "
        "within this comparison and evidence family; study q-value also corrects across "
        "all configured comparisons in the same evidence family. "
        f"{sum(value > 0 for value in counts.values())} of {len(comparisons)} "
        "comparisons have positive decision-candidate discovery signatures."
    )
    comparison_choices = _comparison_choices(
        descriptions={
            comparison: f"{display_names[comparison]} · {int(counts[comparison]):,} enriched"
            for comparison in comparisons
        }
    )
    selected_comparison = comparison_choices[st.selectbox("Comparison", tuple(comparison_choices))]
    selected_outcome = outcomes.loc[outcomes["comparison_id"] == selected_comparison].iloc[0]
    st.caption(
        f"Comparison outcome: {selected_outcome['analysis_status'].replace('_', ' ').lower()}."
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
    evidence_tier = st.selectbox(
        "Evidence to prioritise",
        ("Discovery", "Validated within comparison", "Validated study-wide"),
        help=(
            "Discovery proposes candidates. Validation requires positive enrichment "
            "in held-out proteins; study-wide validation also corrects across comparisons."
        ),
    )
    ranked = _signature_evidence(frame=frame, tier=evidence_tier)
    st.caption(
        f"{len(ranked):,} positive {evidence_tier.lower()} decision candidates among the selected "
        "feature types. Sequence k-mers and whole-model clusters represent different "
        "scales of evidence."
    )
    if ranked.empty:
        st.info(
            "No positive signature meets this tier for the selected comparison and "
            "feature types. An insufficient comparison is not evidence of absence."
        )
    else:
        shortlist = ranked.groupby("feature_type", sort=False).head(5).copy()
        shortlist["feature_label"] = [
            _feature_plot_label(
                feature_type=str(row.feature_type),
                feature_id=str(row.feature_id),
                feature_name=str(row.feature_name) if pd.notna(row.feature_name) else "",
            )
            for row in shortlist.itertuples(index=False)
        ]
        st.subheader("Shortlist across evidence families")
        st.caption(
            "Up to five features from each family, ranked first by the selected q-value "
            "and then by target − background prevalence. Inspect counts and controls "
            "before interpreting a candidate."
        )
        columns = [
            "feature_type",
            "feature_name",
            "feature_id",
            "discovery_prevalence_difference",
            "discovery_q_value",
            "validation_prevalence_difference",
            "validation_q_value",
            "validation_study_q_value",
            "evidence_class",
        ]
        _render_downloadable_table(
            frame=_explain_feature_rows(
                frame=shortlist[[column for column in columns if column in shortlist]]
            ),
            download_name=f"{selected_comparison}_{evidence_tier}_shortlist",
            height=360,
        )
        difference_column = (
            "discovery_prevalence_difference"
            if evidence_tier == "Discovery"
            else "validation_prevalence_difference"
        )
        figure = px.bar(
            shortlist.sort_values(difference_column),
            x=difference_column,
            y="feature_label",
            color="feature_type",
            orientation="h",
            hover_data=["feature_id", "discovery_q_value", "evidence_class"],
            labels={
                difference_column: "Target − background prevalence",
                "feature_label": "Feature",
                "feature_type": "Evidence family",
            },
        )
        figure.update_layout(height=max(360, 33 * len(shortlist)))
        _render_plotly_figure(
            figure=figure,
            download_name=f"{selected_comparison}_{evidence_tier}_shortlist_chart",
            explanation=(
                "Each bar is the target-minus-background prevalence for a shortlisted "
                "feature at the selected evidence tier. More positive means the "
                "feature is more common in the target; colour distinguishes evidence "
                "families. Hover to see its ID and discovery q-value. Examine "
                "held-out results and cohort counts before making a biological claim."
            ),
        )
    with st.expander("Complete signature ledger and downloads", expanded=False):
        _render_downloadable_table(
            frame=_explain_feature_rows(frame=frame),
            download_name=f"{selected_comparison}_signatures",
            height=480,
        )
    chart_data = _signature_plot_data(frame=frame)
    if not chart_data.empty:
        axis = st.selectbox(
            label="q-value axis",
            options=("−log10 q-value", "Recorded q-value"),
            help=(
                "The logarithmic view separates very small q-values. Recorded zeros and "
                "values below 1e-300 are capped at 300 for display and flagged in hover. "
                "The table and downloads retain the original q-values."
            ),
        )
        logarithmic = axis == "−log10 q-value"
        figure = px.scatter(
            data_frame=chart_data,
            x="discovery_prevalence_difference",
            y="plot_negative_log10_q_value" if logarithmic else "discovery_q_value",
            color="feature_type",
            symbol="evidence_class",
            hover_data={
                "feature_id": True,
                "feature_name": True,
                "discovery_q_value": ":.3g",
                "plot_q_value_clipped": True,
            },
            labels={
                "discovery_prevalence_difference": "Target − background prevalence",
                "discovery_q_value": "Discovery q-value",
                "plot_negative_log10_q_value": "−log10 discovery q-value (display cap: 300)",
                "plot_q_value_clipped": "q-value below display floor",
            },
        )
        if logarithmic:
            figure.add_hline(
                y=-np.log10(0.05),
                line_dash="dash",
                annotation_text="q = 0.05 (viewer reference)",
            )
            clipped = int(chart_data["plot_q_value_clipped"].sum())
            if clipped:
                st.caption(
                    body=f"{clipped:,} recorded zero or below-floor q-values appear at "
                    "the display cap of 300. Their exact significance is unresolved at "
                    "that scale; recorded values remain in hover and downloads."
                )
        else:
            figure.update_yaxes(autorange="reversed")
        _render_plotly_figure(figure=figure, download_name="signature_scatter")
    if len(chart_data) < len(frame):
        st.caption(
            body=f"{len(frame) - len(chart_data):,} rows without finite, in-range q-values "
            "and prevalence differences are omitted from this chart and retained in the table."
        )
    st.caption("Association counts may take longer to load for large campaigns.")
    association_key = (
        f"signature-associations-{database}-{selected_comparison}-{'-'.join(selected_types)}"
    )
    if st.button("Load association counts, intervals and multiplicity"):
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
        st.session_state[association_key] = association_frame
    if association_key in st.session_state:
        _render_downloadable_table(
            frame=st.session_state[association_key],
            download_name=f"{selected_comparison}_association_counts",
            height=440,
        )


def _render_explainable_ml(*, database: Path) -> None:
    """Render held-out prediction, global importance and local explanations.

    Args:
        database: Verified result database.
    """

    st.title("Explainable prediction")
    _render_page_help(page="Explainable prediction")
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
    catalogue = query_dataframe(
        database=database,
        sql=(
            "SELECT m.comparison_id, c.display_name, m.status, "
            "m.validation_target_count FROM ml_models m "
            "LEFT JOIN comparisons c USING (comparison_id) "
            "ORDER BY m.comparison_id"
        ),
    )
    if catalogue.empty:
        st.info("No model outcomes were published for these comparisons.")
        return
    catalogue["has_validation"] = catalogue["status"].eq("COMPLETE") & pd.to_numeric(
        catalogue["validation_target_count"], errors="coerce"
    ).gt(0)
    catalogue = catalogue.sort_values(
        ["has_validation", "display_name", "comparison_id"],
        ascending=[False, True, True],
    )
    st.caption(
        f"{int(catalogue['has_validation'].sum()):,} of {len(catalogue):,} "
        "comparisons have a complete model with held-out target proteins. "
        "Other comparisons retain their explicit reason for not being fitted."
    )
    model_choices = _comparison_choices(
        descriptions={
            str(row.comparison_id): (
                f"{row.display_name} · {str(row.status).replace('_', ' ').lower()}"
            )
            for row in catalogue.itertuples(index=False)
        }
    )
    comparison_id = model_choices[st.selectbox("Model comparison", tuple(model_choices))]
    LOGGER.debug("Rendering published model comparison %s", comparison_id)
    model = query_dataframe(
        database=database,
        sql="SELECT * FROM ml_models WHERE comparison_id = ?",
        parameters=(comparison_id,),
    )
    with st.expander("Complete model record and download"):
        _render_downloadable_table(
            frame=model,
            download_name=f"{comparison_id}_model_summary",
        )
    if model.empty or not str(model.iloc[0]["status"]).startswith("COMPLETE"):
        reason = str(model.iloc[0].get("status_message") or "") if not model.empty else ""
        st.warning(
            "No fitted model is available for this comparison. "
            + (reason or "Inspect the published model status and sample counts above.")
        )
        return
    record = model.iloc[0]
    st.caption(body=str(record.get("status_message") or ""))
    if str(record["status"]) != "COMPLETE":
        st.warning(
            body=status_explanation(
                status=str(record["status"]),
                metadata=st.session_state.get("help-result-metadata", {}),
                context="model",
            )
        )
    columns = st.columns(4)
    columns[0].metric(
        label="CV ROC AUC",
        value=_metric_text(record["cv_roc_auc"]),
        help=METRIC_HELP["cv_roc_auc"],
    )
    columns[1].metric(
        label="Validation ROC AUC",
        value=_metric_text(record["validation_roc_auc"]),
        help=METRIC_HELP["validation_roc_auc"],
    )
    columns[2].metric(
        label="Validation average precision (AP)",
        value=_metric_text(record["validation_average_precision"]),
        help=METRIC_HELP["validation_average_precision"],
    )
    columns[3].metric(
        label="Validation MCC",
        value=_metric_text(record["validation_matthews_correlation"]),
        help=METRIC_HELP["validation_matthews_correlation"],
    )
    _render_metric_guide(database=database, comparison_id=comparison_id, record=record.to_dict())
    _render_individual_shap_summary(database=database, comparison_id=comparison_id)
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
        "balanced-accuracy loss on pure blocks, calculated only when both per-class "
        "sample/group limits pass. This preview contains up to 200 features ranked "
        "by absolute coefficient; these ranks are different from the SHAP bar ordering."
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
            labels={
                "coefficient_log_odds": "Fitted coefficient (log-odds)",
                "feature_name": "Feature",
            },
        )
        _render_plotly_figure(
            figure=figure,
            download_name=f"{comparison_id}_model_coefficients",
            explanation=(
                "Signed model coefficients show which features raise or lower the "
                "fitted target log-odds, conditional on the other model features. "
                "They are model contributions, not causal effects or enrichment tests."
            ),
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
        explanation=(
            "The bars show the distribution of predicted target probabilities for "
            "proteins in the selected partition, split by their known class. "
            "Overlapping groups indicate uncertain discrimination; consult the "
            "held-out metrics and sample counts."
        ),
    )
    protein_id = st.selectbox(
        f"Explain {explained_partition.lower()} protein",
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


def _render_individual_shap_summary(*, database: Path, comparison_id: str) -> None:
    """Show the largest published individual validation SHAP contributions.

    Args:
        database: Verified result database.
        comparison_id: Exact selected model comparison.
    """
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT feature_type, feature_id, feature_name, "
            "mean_absolute_validation_contribution FROM ml_feature_importance "
            "WHERE comparison_id = ? AND mean_absolute_validation_contribution IS NOT NULL "
            "AND isfinite(mean_absolute_validation_contribution) "
            "AND mean_absolute_validation_contribution >= 0 "
            "ORDER BY mean_absolute_validation_contribution DESC, feature_type, feature_id LIMIT 30"
        ),
        parameters=(comparison_id,),
    )
    st.subheader("Individual feature contributions in validation")
    required = {
        "feature_type",
        "feature_id",
        "feature_name",
        "mean_absolute_validation_contribution",
    }
    if frame.empty or not required.issubset(frame.columns):
        st.info("No individual validation SHAP summary was published for this comparison.")
        return
    chart = frame.copy()
    chart["plot_feature_label"] = (
        chart["feature_name"].astype(str).str.slice(stop=45)
        + " · "
        + chart["feature_id"].astype(str)
    )
    st.caption(
        body="Up to 30 individual features ranked by mean absolute SHAP contribution across "
        "validation proteins, in model log-odds. Each bar is one feature; the aggregate "
        "other-features remainder is excluded. Magnitude does not show direction or causation. "
        "The original published SHAP figures remain available below."
    )
    figure = px.bar(
        data_frame=chart,
        x="mean_absolute_validation_contribution",
        y="plot_feature_label",
        color="feature_type",
        orientation="h",
        hover_data=["feature_id", "feature_name"],
        labels={
            "mean_absolute_validation_contribution": "Mean absolute SHAP contribution (log-odds)",
            "plot_feature_label": "Individual feature",
            "feature_type": "Evidence family",
        },
    )
    figure.update_yaxes(autorange="reversed", automargin=True)
    figure.update_layout(height=max(400, 24 * len(chart)))
    _render_plotly_figure(
        figure=figure, download_name=f"{comparison_id}_validation_individual_shap"
    )


def _render_metric_guide(*, database: Path, comparison_id: str, record: dict[str, object]) -> None:
    """Interpret every model metric beside actual pure-block class counts.

    Args:
        database: Verified result database.
        comparison_id: Exact selected model comparison.
        record: Published selected model record.
    """
    counts = query_dataframe(
        database=database,
        sql=(
            "SELECT coalesce(count_if(true_class = 'TARGET'), 0) AS target_blocks, "
            "coalesce(count_if(true_class = 'BACKGROUND'), 0) AS background_blocks FROM ("
            "SELECT partition_key, min(true_class) AS true_class FROM ml_predictions "
            "WHERE comparison_id = ? AND partition = 'VALIDATION' GROUP BY partition_key "
            "HAVING count(DISTINCT true_class) = 1)"
        ),
        parameters=(comparison_id,),
    )
    baseline = None
    with st.expander(label="❓ ROC, AUC, AP and MCC — interpret these values", expanded=False):
        if not counts.empty and {"target_blocks", "background_blocks"}.issubset(counts.columns):
            target = int(counts.iloc[0]["target_blocks"])
            background = int(counts.iloc[0]["background_blocks"])
            if target + background:
                baseline = target / (target + background)
                st.markdown(
                    body=f"Held-out pure blocks: **{target} target**, **{background} background**. "
                    f"Target fraction for the AP baseline: **{baseline:.3f}**."
                )
        for name, definition in METRIC_HELP.items():
            st.markdown(body=f"**{name.replace('_', ' ')}** — {definition}")
            st.caption(
                body=metric_reading(metric_name=name, value=record.get(name), baseline=baseline)
            )
        st.markdown(
            body=(
                "Metrics use aggregated pure independence blocks; prediction/SHAP tables show "
                "individual proteins. Their apparent separation can differ. A high value is "
                "performance against supplied labels in this sampled cohort, not a probability "
                "of biological correctness. No confidence interval or external-proteome "
                "evaluation is supplied by these metric cards."
            )
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
        _render_graph_help(graph_name=str(row["plot_type"]))
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


def _with_feature_explanations(*, frame: pd.DataFrame) -> pd.DataFrame:
    """Add readable feature semantics to a detached visible table.

    Args:
        frame: Normalised visible table with unique field names.

    Returns:
        Copy preserving original values and any authoritative explanation column.
    """
    result = frame.copy(deep=True)
    if not {"feature_type", "feature_id"}.issubset(result.columns):
        return result
    if "feature_explanation" in result.columns:
        return result
    descriptions = {}
    values = []
    for kind, identifier in zip(result["feature_type"], result["feature_id"], strict=True):
        if not isinstance(kind, str) or not isinstance(identifier, str) or not identifier:
            values.append("No feature identifier is recorded on this outcome row.")
            continue
        key = (kind, identifier)
        if key not in descriptions:
            descriptions[key] = feature_explanation(feature_type=kind, feature_id=identifier)
        values.append(descriptions[key])
    result["feature_explanation"] = values
    return result


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
    normalised = _with_feature_explanations(frame=normalised)
    display_options: dict[str, object] = {
        "hide_index": True,
        "width": "stretch",
        "column_config": _table_column_configuration(frame=normalised),
    }
    display_options["height"] = min(height or 460, 38 + 35 * max(1, min(len(normalised), 12)))
    priority = (
        "protein_id",
        "feature_name",
        "feature_explanation",
        "display_name",
        "status",
        "evidence_class",
        "true_class",
        "predicted_probability",
        "predicted_class",
        "discovery_prevalence_difference",
        "discovery_q_value",
        "validation_q_value",
        "validation_study_q_value",
        "prevalence_difference",
        "q_value",
        "study_q_value",
    )
    front = [name for name in priority if name in normalised.columns]
    display_options["column_order"] = front + [
        name for name in normalised.columns if name not in front
    ]
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
    display = normalised.copy()
    for column in display:
        if column.endswith("q_value") or column == "q_value":
            display[column] = pd.to_numeric(display[column], errors="coerce").map(
                lambda value: "—" if pd.isna(value) else f"{value:.2e}"
            )
        elif column.endswith("prevalence_difference"):
            display[column] = pd.to_numeric(display[column], errors="coerce").map(
                lambda value: "—" if pd.isna(value) else f"{value:+.3f}"
            )
    display_options["column_config"] = _table_column_configuration(frame=display)
    st.dataframe(display, **display_options)
    _render_table_help(frame=normalised)
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


def _render_plotly_figure(*, figure: object, download_name: str, explanation: str = "") -> None:
    """Render an interactive figure with on-demand PDF, PNG and HTML exports.

    Args:
        figure: Plotly-compatible figure.
        download_name: Stable human-readable export identity.
        explanation: What this plot encodes and how to interpret it.
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
    _render_graph_help(graph_name=download_name, explanation=explanation)
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
    _render_page_help(page="Protein & Pfam")
    _feature_key()
    suggestions = query_dataframe(
        database=database,
        sql=(
            "SELECT m.protein_id, "
            "max(CASE WHEN s.coordinate_path <> '' THEN 1 ELSE 0 END) AS has_model, "
            "count(DISTINCT m.label_id) AS label_count "
            "FROM label_memberships m LEFT JOIN structures s USING (protein_id) "
            "GROUP BY m.protein_id ORDER BY has_model DESC, label_count DESC, "
            "m.protein_id LIMIT 250"
        ),
    )
    protein_ids = tuple(suggestions["protein_id"].astype(str)) if not suggestions.empty else ()
    st.caption(
        "Suggested proteins have published class memberships, with available models "
        "first. Enter an exact protein ID to inspect any other published sequence."
    )
    selected_id = st.selectbox("Suggested labelled protein", protein_ids) if protein_ids else ""
    exact_id = st.text_input("Exact published protein ID (optional)").strip()
    protein_id = exact_id or selected_id
    if not protein_id:
        st.info("Enter a published protein ID to inspect its evidence.")
        return
    inventory = query_dataframe(
        database=database,
        sql=(
            "SELECT protein_id, description, sequence_length, sequence_sha256 "
            "FROM proteins WHERE protein_id = ?"
        ),
        parameters=(protein_id,),
    )
    if inventory.empty:
        st.warning(f"Protein {protein_id!r} is not present in this completed result.")
        return
    _render_downloadable_table(
        frame=inventory,
        download_name=f"{protein_id}_protein_inventory",
    )
    section = st.radio(
        "Protein evidence",
        ("Labels", "Domains", "Positive features", "Feature assessments"),
        horizontal=True,
        help="Load one evidence source at a time; feature membership can be very large.",
    )
    if section == "Labels":
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
    elif section == "Domains":
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
    elif section == "Positive features":
        _render_downloadable_table(
            frame=_explain_feature_rows(
                frame=query_dataframe(
                    database=database,
                    sql=(
                        "SELECT * FROM features WHERE protein_id = ? "
                        'ORDER BY feature_type, "start", feature_id'
                    ),
                    parameters=(protein_id,),
                )
            ),
            download_name=f"{protein_id}_features",
            height=450,
        )
    else:
        st.caption(
            "Assessed absence, unassessed, failed and excluded evidence remain distinct "
            "from positive feature evidence."
        )
        _render_downloadable_table(
            frame=_explain_feature_rows(
                frame=query_dataframe(
                    database=database,
                    sql=(
                        "SELECT * FROM feature_assessments WHERE protein_id = ? "
                        'ORDER BY feature_type, feature_id, "start"'
                    ),
                    parameters=(protein_id,),
                )
            ),
            download_name=f"{protein_id}_feature_assessments",
            height=450,
        )


def _matched_control_coverage(*, database: Path) -> pd.DataFrame:
    """Count target units with a match within each pooled background stratum.

    Args:
        database: Verified result database.

    Returns:
        Per-background coverage with its explicit denominator.
    """

    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT background_label_id, "
            "count(DISTINCT target_unit_id) AS target_units, "
            "count(DISTINCT CASE WHEN status = 'MATCHED' "
            "THEN target_unit_id END) AS covered_target_units, "
            "count(DISTINCT CASE WHEN status = 'MATCHED' "
            "THEN control_unit_id END) AS matched_control_units "
            "FROM control_matching_audit GROUP BY background_label_id "
            "ORDER BY background_label_id"
        ),
    )
    if not frame.empty:
        frame["target_coverage_fraction"] = frame["covered_target_units"] / frame[
            "target_units"
        ].replace(0, float("nan"))
    return frame


def _matched_comparison_coverage(
    *, database: Path, metadata: dict[str, object] | None
) -> pd.DataFrame:
    """Read the exact cohorts used for each audited comparison.

    Args:
        database: Verified result database.
        metadata: Published campaign metadata with matched cohort counts.

    Returns:
        One row per audited comparison, with its target-unit denominator.
    """

    evidence = (metadata or {}).get("automated_label_evidence")
    cohorts = evidence.get("matched_comparison_cohorts") if isinstance(evidence, dict) else None
    if not isinstance(cohorts, dict) or not cohorts:
        return pd.DataFrame()
    comparisons = query_dataframe(
        database=database,
        sql="SELECT comparison_id, display_name FROM comparisons ORDER BY comparison_id",
    )
    names = (
        dict(zip(comparisons["comparison_id"], comparisons["display_name"], strict=True))
        if not comparisons.empty
        else {}
    )
    records = []
    fields = (
        "matched_target_unit_count",
        "excluded_unmatched_target_unit_count",
        "target_protein_count",
        "excluded_unmatched_target_protein_count",
        "control_unit_count",
        "control_protein_count",
    )
    for comparison_id, cohort in sorted(cohorts.items()):
        if not isinstance(cohort, dict) or any(
            not isinstance(cohort.get(field), int)
            or isinstance(cohort[field], bool)
            or cohort[field] < 0
            for field in fields
        ):
            continue
        total = cohort["matched_target_unit_count"] + cohort["excluded_unmatched_target_unit_count"]
        records.append(
            {
                "comparison_id": comparison_id,
                "display_name": names.get(comparison_id, comparison_id),
                "target_units": total,
                "matched_target_units": cohort["matched_target_unit_count"],
                "excluded_unmatched_target_units": cohort["excluded_unmatched_target_unit_count"],
                "target_coverage_fraction": (
                    cohort["matched_target_unit_count"] / total if total else float("nan")
                ),
                "analysed_target_proteins": cohort["target_protein_count"],
                "excluded_unmatched_target_proteins": cohort[
                    "excluded_unmatched_target_protein_count"
                ],
                "matched_control_units": cohort["control_unit_count"],
                "matched_control_proteins": cohort["control_protein_count"],
            }
        )
    return pd.DataFrame.from_records(records)


def _render_classes(*, database: Path, metadata: dict[str, object] | None = None) -> None:
    """Render the configured hierarchy and observed class/role coverage.

    Args:
        database: Verified result database.
        metadata: Published comparison-specific cohort counts, when available.
    """

    st.title("Classes & component roles")
    _render_page_help(page="Classes & roles")
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
        comparison_coverage = _matched_comparison_coverage(database=database, metadata=metadata)
        if not comparison_coverage.empty:
            st.subheader("Matched cohort used for each comparison")
            st.caption(
                "The denominator includes matched and excluded target units. Only "
                "matched targets and their allocated controls enter that comparison's "
                "association and prediction analysis."
            )
            _render_downloadable_table(
                frame=comparison_coverage,
                download_name="comparison_specific_matched_cohorts",
                height=420,
            )
        else:
            st.info(
                "No comparison-specific matched cohort summary was published for "
                "this resource. Background control totals cannot substitute for it."
            )
        st.subheader("Provisional target labels and background pools")
        st.caption(
            "These direct label counts and pooled background controls are descriptive. "
            "Several comparisons may share one background. Where published, use the "
            "comparison-specific matched cohort for the actual analysis denominator."
        )
        target_rows = provisional.loc[provisional["label_type"] == "TARGET"].copy()
        _render_downloadable_table(
            frame=target_rows,
            download_name="provisional_target_label_and_control_coverage",
            height=420,
        )
        with st.expander("Background labels and full assignment summary"):
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
            "count(DISTINCT m.protein_id) AS observed_proteins "
            "FROM profile_labels p LEFT JOIN label_memberships m USING (label_id) "
            "GROUP BY ALL ORDER BY p.label_id"
        ),
    )
    populated = coverage[coverage["observed_proteins"] > 0]
    left, right = st.columns((3, 2))
    with left:
        st.subheader("Populated class labels")
        st.caption(
            "Parent and child memberships overlap. Bars are not additive parts "
            "of one total; labels can be reviewed or provisional under the recorded authority."
        )
        if populated.empty:
            st.info("No positive profile memberships are present.")
        else:
            chart_rows = populated.loc[populated["parent_label_id"].ne("")].nlargest(
                15, "observed_proteins"
            )
            if not chart_rows.empty:
                figure = px.bar(
                    chart_rows.sort_values("observed_proteins"),
                    x="observed_proteins",
                    y="display_name",
                    orientation="h",
                    hover_data=["label_id", "parent_label_id", "component_role"],
                    labels={
                        "observed_proteins": "Proteins with this label",
                        "display_name": "Class label",
                    },
                )
                figure.update_layout(height=max(420, 34 * len(chart_rows)))
                _render_plotly_figure(
                    figure=figure,
                    download_name="most_populated_class_labels",
                    explanation=(
                        "Bars show up to 15 populated child labels. The same protein "
                        "can appear in a parent and several related labels, so bars "
                        "must not be added to obtain a unique protein total. The "
                        "chart shows membership, not class-specific enrichment."
                    ),
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
    _render_page_help(page="Structures & folds")
    _feature_key()
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
        _render_plotly_figure(
            figure=figure,
            download_name="structure_eligibility_summary",
            explanation=(
                "Bars count model records by analysis eligibility and fold-evidence "
                "status. Hover to identify sources and average model confidence. "
                "A coordinate model can be eligible without a named fold or "
                "local structural motif."
            ),
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
        st.caption(
            "These are connected whole-model groups. At permissive similarity and "
            "coverage thresholds, single-linkage chains can join diverse proteins; "
            "a large group is not a shared local motif."
        )
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
        if not clusters.empty:
            clusters.insert(
                1,
                "cluster_explanation",
                [
                    _feature_explanation(
                        feature_type="STRUCTURE_CLUSTER",
                        feature_id=str(cluster_id),
                        feature_name="",
                    )
                    for cluster_id in clusters["cluster_id"]
                ],
            )
        _render_downloadable_table(
            frame=clusters,
            download_name="structural_clusters",
            height=380,
        )
    st.subheader("Pairwise structural alignments")
    density = query_dataframe(
        database=database,
        sql=(
            "SELECT comparison_tool, "
            "round(floor(least(coverage_a, coverage_b) * 20) / 20, 2) "
            "AS minimum_coverage_bin, "
            "round(floor(tm_score * 20) / 20, 2) AS tm_score_bin, "
            "count(*) AS comparison_count FROM structure_comparisons "
            "WHERE tm_score IS NOT NULL AND coverage_a IS NOT NULL "
            "AND coverage_b IS NOT NULL GROUP BY ALL "
            "ORDER BY minimum_coverage_bin, tm_score_bin, comparison_tool"
        ),
    )
    if not density.empty:
        st.caption(
            f"All {int(density['comparison_count'].sum()):,} recorded pairs with scores "
            "and bilateral coverage contribute to this binned view."
        )
        figure = px.scatter(
            density,
            x="minimum_coverage_bin",
            y="tm_score_bin",
            size="comparison_count",
            color="comparison_tool",
            size_max=32,
            hover_data=["comparison_count"],
            labels={
                "minimum_coverage_bin": "Minimum bilateral coverage (0.05 bins)",
                "tm_score_bin": "TM-score (0.05 bins)",
                "comparison_count": "Pairs in bin",
            },
        )
        _render_plotly_figure(
            figure=figure,
            download_name="all_structural_alignment_score_coverage_bins",
            explanation=(
                "Every recorded pair with both coverage values contributes to one "
                "bin. Larger circles contain more pairs. TM-score indicates model "
                "similarity, and the horizontal position is the smaller aligned "
                "fraction. The previous high-score-only sample was unsuitable for "
                "interpreting the full distribution."
            ),
        )
        _render_downloadable_table(
            frame=density,
            download_name="all_structural_alignment_score_coverage_bins",
        )
    if st.checkbox("Inspect the 5,000 highest-scoring structural pairs"):
        st.caption(
            "This is a score-ranked subset for inspecting individual pairs. "
            "It is deliberately not used to estimate the score distribution."
        )
        comparisons = query_dataframe(
            database=database,
            sql=(
                "SELECT * FROM structure_comparisons ORDER BY tm_score DESC NULLS LAST LIMIT 5000"
            ),
        )
        _render_downloadable_table(
            frame=comparisons,
            download_name="highest_scoring_pairwise_structural_alignments",
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
        st.subheader("Imported exploratory structural summaries")
        st.caption(
            "These Stage 09b group and pocket-assessment summaries were imported for "
            "context only. They were not enrichment-tested as residue-level pockets "
            "in this campaign and provide no motif coordinates to colour."
        )
        _render_downloadable_table(
            frame=imported,
            download_name="imported_pocket_conservation",
            height=420,
        )


def _render_model_explorer(*, database: Path) -> None:
    """Explore enrichment on exact sequence positions, models and aligned pairs."""

    st.title("Model & alignment explorer")
    _render_page_help(page="Model & alignment explorer")
    _feature_key()
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
    overview = _comparison_outcomes(database=database)
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
        outcomes = overview.drop(columns=["target_label_ids"])
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
            index=(
                2
                if int(selected["validated_study_count"]) > 0
                else 1
                if int(selected["validated_within_count"]) > 0
                else 0
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
    st.caption(
        f"{len(valid):,} localisable feature occurrences cover "
        f"{sum(value > 0 for value in scores):,} of {len(sequence):,} residues. "
        "Whole-model structural clusters have no residue boundaries and cannot colour "
        "a specific motif."
    )
    if valid:
        _render_downloadable_table(
            frame=_explain_feature_rows(frame=pd.DataFrame(valid)),
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
    if st.checkbox("Inspect a paired sequence alignment"):
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
    if st.checkbox("Compare ranked signatures across classes"):
        _render_ranked_associations(
            database=database,
            enriched_comparisons=enriched_comparisons,
            descriptions=descriptions,
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

    st.subheader("Discovery candidates across classes")
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
                "AND starts_with(s.evidence_class, 'DECISION_CANDIDATE__') "
                "QUALIFY within_comparison_rank <= ? "
                "ORDER BY within_comparison_rank, s.comparison_id LIMIT 5000"
            ),
            parameters=(*selected_results, rows_per_comparison),
        )
        _render_downloadable_table(
            frame=_explain_feature_rows(frame=ranked),
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
        strongest["label"] = [
            f"{str(row.comparison)[:24]} · "
            + _feature_plot_label(
                feature_type=str(row.feature_type),
                feature_id=str(row.feature_id),
                feature_name=str(row.feature_name) if pd.notna(row.feature_name) else "",
            )
            for row in strongest.itertuples(index=False)
        ]
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
            explanation=(
                "Bars show target-minus-background prevalence for the top selected "
                "discovery features, coloured by comparison. Rows are balanced "
                "between selected classes. The shortlist is ranked by discovery "
                "q-value, so also compare effect size and validation."
            ),
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
        "Discovery within comparison": (
            "AND starts_with(s.evidence_class, 'DECISION_CANDIDATE__') "
        ),
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
    _render_plotly_figure(
        figure=figure,
        download_name=download_name,
        explanation=(
            "The horizontal axis is the selected protein's 1-based sequence position. "
            "White has no mapped significant feature; blue to red encodes smaller "
            "adjusted q-values of overlapping enriched features. A coloured residue "
            "inherits a protein-level feature result; it does not have its own test."
        ),
    )


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
        height=560,
        margin={"l": 0, "r": 0, "t": 48, "b": 0},
        scene={
            "aspectmode": "data",
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
        },
    )
    plot_column, explanation_column = st.columns((3, 2))
    with plot_column:
        _render_plotly_figure(
            figure=figure,
            download_name=download_name,
            explanation=(
                "This is a rotatable Cα backbone trace of one exact-sequence model. "
                "When the model matches, coloured markers inherit the mapped "
                "feature significance from the 2D track. White is unmapped; "
                "colour does not show pLDDT, binding pockets or a structural "
                "superposition."
            ),
        )
    with explanation_column:
        st.markdown("**How to read this model**")
        st.caption(f"{len(residues):,} Cα coordinates; drag the model to rotate it.")
        if scores is not None:
            st.caption(
                f"{sum(value > 0 for value in scores):,} residues have localised "
                "significant enrichment. White marks no mapped enrichment. Colour "
                "strength reflects adjusted significance, not pLDDT or effect size."
            )
        else:
            st.caption("No enrichment was projected because the model did not match exactly.")


def _alignment_residue_strip(
    *, rows: list[dict[str, object]], reference_id: str, comparison_id: str
) -> str:
    """Build an escaped, horizontally scrollable amino-acid alignment window.

    Args:
        rows: Consecutive sequence-alignment columns with enrichment scores.
        reference_id: Published identifier for the first sequence.
        comparison_id: Published identifier for the second sequence.

    Returns:
        Static HTML with coloured residues, gaps and identity markers.
    """

    def cells(*, prefix: str) -> str:
        """Render one safe sequence row from the current alignment window.

        Args:
            prefix: ``reference`` or ``comparison``.

        Returns:
            Consecutive, tooltip-labelled residue cells.
        """

        values = []
        for row in rows:
            amino_acid = html.escape(str(row[f"{prefix}_residue"]), quote=True)
            position = row[f"{prefix}_position"]
            score = max(0.0, min(1.0, float(row[f"{prefix}_enrichment"])))
            hue = round(215 * (1 - max(0.0, score - 0.1) / 0.9))
            background = f"hsl({hue}, 68%, 78%)" if score > 0 else "#ffffff"
            foreground = "#192b40" if score > 0 else "#586573"
            label = html.escape(
                f"{prefix} position {position or 'gap'}; enrichment colour {score:.2f}",
                quote=True,
            )
            values.append(
                f'<span title="{label}" style="display:inline-block;width:1.35em;'
                "box-sizing:border-box;"
                f"text-align:center;background:{background};color:{foreground};"
                f'border:1px solid #dae1e7">{amino_acid}</span>'
            )
        return "".join(values)

    identities = "".join(
        '<span style="display:inline-block;width:1.35em;box-sizing:border-box;'
        'text-align:center">' + ("|" if row["identity"] else "&nbsp;") + "</span>"
        for row in rows
    )
    reference = html.escape(reference_id, quote=True)
    comparison = html.escape(comparison_id, quote=True)
    return (
        '<div style="max-width:100%;overflow-x:auto;border:1px solid #dae1e7;'
        'border-radius:8px;padding:0.75rem;font:14px monospace">'
        f'<div style="white-space:nowrap"><strong title="{reference}" '
        f'style="display:inline-block;width:12rem;overflow:hidden">{reference}</strong>'
        f"{cells(prefix='reference')}</div>"
        '<div style="white-space:nowrap"><strong style="display:inline-block;'
        f'width:12rem">Identity</strong>{identities}</div>'
        f'<div style="white-space:nowrap"><strong title="{comparison}" '
        'style="display:inline-block;width:12rem;overflow:hidden">'
        f"{comparison}</strong>"
        f"{cells(prefix='comparison')}</div></div>"
    )


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
    st.markdown("**Residues in this alignment window**")
    st.caption(
        "Each letter is an aligned amino acid; dashes are gaps and vertical marks are "
        "exact matches. Hover over a letter for its sequence position. These are "
        "sequence columns, not structurally superposed coordinates."
    )
    st.markdown(
        _alignment_residue_strip(rows=window, reference_id=protein_id, comparison_id=other_id),
        unsafe_allow_html=True,
    )
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
        explanation=(
            "The two rows show a window of the newly computed amino-acid alignment; "
            "columns line up by sequence alignment. Colour comes independently from "
            "significant features in each protein. Read the letter strip above for "
            "residues and gaps. These are not Foldseek superposed-residue columns."
        ),
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
    _render_page_help(page="Glossary & help")
    term = st.text_input("Search definitions", help="Search terms, categories and explanations.")
    frame = pd.DataFrame(
        glossary_rows(base_rows=GLOSSARY), columns=["category", "term", "definition"]
    )
    if term.strip():
        mask = frame.apply(
            lambda column: column.astype(str).str.contains(term.strip(), case=False, regex=False)
        ).any(axis=1)
        frame = frame.loc[mask].reset_index(drop=True)
    _render_downloadable_table(frame=frame, download_name="protein_signature_glossary")


def _render_orthology(*, database: Path, metadata: dict[str, object] | None = None) -> None:
    """Render OrthoFinder context and leakage-safe data partitions.

    Args:
        database: Verified result database.
        metadata: Optional published evidence-availability states.
    """

    st.title("Orthology & partitions")
    _render_page_help(page="Orthology & partitions")
    st.caption(
        "The composite OrthoFinder authority is run ID + group type + hierarchy node + "
        "group ID. Whole connected homology/redundancy blocks stay in one partition."
    )
    columns = st.columns(4)
    columns[0].metric(
        "Memberships",
        f"{table_count(database=database, table_name='orthofinder_memberships'):,}",
        help=(
            "Published protein-to-group records; one protein can belong to several hierarchy nodes."
        ),
    )
    groups = query_dataframe(
        database=database,
        sql=(
            "SELECT count(DISTINCT (run_id, group_type, hierarchy_node, group_id)) AS n "
            "FROM orthofinder_memberships"
        ),
    )
    columns[1].metric(
        "Groups",
        f"{int(groups.iloc[0]['n']):,}",
        help="Distinct composite keys: run ID, group type, hierarchy node and group ID.",
    )
    blocks = query_dataframe(
        database=database,
        sql="SELECT count(DISTINCT partition_key) AS n FROM partitions",
    )
    columns[2].metric(
        "Partition blocks",
        f"{int(blocks.iloc[0]['n']):,}",
        help=column_definition(column_name="partition_key"),
    )
    near = query_dataframe(
        database=database,
        sql=(
            "SELECT count(DISTINCT cluster_id) AS n FROM redundancy_clusters "
            "WHERE cluster_type = 'NEAR_REDUNDANCY'"
        ),
    )
    near_count = int(near.iloc[0]["n"])
    metadata = (
        metadata if metadata is not None else st.session_state.get("help-result-metadata", {})
    )
    availability = metadata.get("evidence_availability", {}) if isinstance(metadata, dict) else {}
    near_status = availability.get("near_redundancy") if isinstance(availability, dict) else None
    near_value = (
        "Not assessed"
        if near_count == 0 and near_status in {"INPUT_UNAVAILABLE", "NOT_SELECTED", "NOT_ASSESSED"}
        else f"{near_count:,}"
    )
    columns[3].metric(
        "Near-redundancy clusters",
        near_value,
        help="Supplied NEAR_REDUNDANCY groups. Zero means none recorded, "
        "not proof that the proteome contains no closely related sequences.",
    )
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
    if not partitions.empty and {"partition", "partition_unit", "blocks"}.issubset(
        partitions.columns
    ):
        figure = px.bar(
            data_frame=partitions,
            x="partition",
            y="blocks",
            color="partition_unit",
            barmode="group",
            hover_data=["proteins"],
            labels={"blocks": "Independent partition blocks", "partition": "Partition"},
        )
        _render_plotly_figure(
            figure=figure,
            download_name="discovery_validation_block_allocation",
            explanation=(
                "Bars count independent allocation blocks by partition and unit "
                "type. Hover for the number of proteins they contain. Whole linked "
                "groups are kept together to reduce discovery/validation leakage; "
                "a block count is not a count of distinct functional orthogroups."
            ),
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
            "No optional OrthoFinder group-context summary was supplied. The membership "
            "and partition counts above remain authoritative. OrthoFinder memberships, "
            "when supplied, still contribute to the connected partition blocks."
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
    _render_page_help(page="Data quality & provenance")
    st.success("Completion marker and all checksums verified when this resource was opened.")
    label_evidence = metadata.get("automated_label_evidence", {})
    if isinstance(label_evidence, dict) and label_evidence.get("status") != "NOT_SELECTED":
        st.subheader("Automated label evidence and circularity safeguards")
        warning = str(label_evidence.get("warning") or "").strip()
        if warning:
            st.warning(warning)
        comparison_coverage = _matched_comparison_coverage(database=database, metadata=metadata)
        if not comparison_coverage.empty:
            st.markdown("**Comparison-specific matched cohorts**")
            _render_downloadable_table(
                frame=comparison_coverage,
                download_name="quality_comparison_matched_cohorts",
            )
        _render_control_matching_coverage(database=database)
        _render_label_evidence_audit(database=database)
    st.subheader("Feature assessment coverage")
    st.caption(
        "A missing positive row is not treated as absence: explicit assessment state and "
        "derivation scope determine each feature's tested universe."
    )
    assessment_coverage = query_dataframe(
        database=database,
        sql=(
            "SELECT feature_type, evidence_status, derivation_scope, "
            "count(DISTINCT feature_id) AS features, "
            "count(DISTINCT protein_id) AS proteins FROM feature_assessments "
            "GROUP BY ALL ORDER BY feature_type, evidence_status, derivation_scope"
        ),
    )
    if "STRUCTURAL_POCKET" in set(assessment_coverage.get("feature_type", ())):
        st.info(
            "STRUCTURAL_POCKET rows here are imported exploratory assessment "
            "records from an upstream structural resource. They are not "
            "enrichment-tested, residue-mapped pocket discoveries in this campaign."
        )
    _render_downloadable_table(
        frame=assessment_coverage,
        download_name="feature_assessment_coverage",
    )
    st.subheader("Domain assessment coverage")
    domain_coverage = query_dataframe(
        database=database,
        sql=(
            "SELECT domain_authority, assessment_status, count(*) AS proteins "
            "FROM domain_assessments GROUP BY ALL ORDER BY domain_authority, assessment_status"
        ),
    )
    if not domain_coverage.empty:
        unassessed = int(
            domain_coverage.loc[
                domain_coverage["assessment_status"] == "NOT_ASSESSED", "proteins"
            ].sum()
        )
        if unassessed:
            total = int(domain_coverage["proteins"].sum())
            st.warning(
                f"{unassessed:,} of {total:,} recorded domain assessments are "
                "NOT_ASSESSED. Missing Pfam evidence for these proteins is unknown, "
                "not a verified domain absence."
            )
        figure = px.bar(
            domain_coverage,
            x="assessment_status",
            y="proteins",
            color="domain_authority",
            labels={"assessment_status": "Assessment state", "proteins": "Proteins"},
        )
        _render_plotly_figure(
            figure=figure,
            download_name="domain_assessment_state_chart",
            explanation=(
                "Bars count proteins by explicit domain-assessment state. A "
                "completed scan with no hit differs from NOT_ASSESSED; the latter "
                "must not be used as negative evidence for a Pfam feature."
            ),
        )
    _render_downloadable_table(
        frame=domain_coverage,
        download_name="domain_assessment_coverage",
    )
    st.subheader("AlphaFold acquisition outcomes")
    acquisition_outcomes = query_dataframe(
        database=database,
        sql=(
            "SELECT acquisition_status, count(*) AS proteins, "
            "round(avg(mean_plddt), 2) AS mean_plddt "
            "FROM alphafold_acquisitions GROUP BY acquisition_status "
            "ORDER BY proteins DESC"
        ),
    )
    if acquisition_outcomes.empty:
        st.info(
            "No new AlphaFold acquisition was selected for this campaign. Packaged "
            "coordinate models from the upstream workflow can still be available."
        )
    _render_downloadable_table(
        frame=acquisition_outcomes,
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
    st.caption("Calculating this audit scans the full association ledger.")
    audit_key = f"inferential-blocks-{database}"
    if st.button("Calculate independence-block summary"):
        st.session_state[audit_key] = query_dataframe(
            database=database,
            sql=(
                "SELECT comparison_id, partition, max(target_unit_count) AS target_blocks, "
                "max(background_unit_count) AS background_blocks, "
                "max(excluded_mixed_unit_count) AS excluded_mixed_blocks "
                "FROM associations GROUP BY comparison_id, partition "
                "ORDER BY comparison_id, partition"
            ),
        )
    if audit_key in st.session_state:
        _render_downloadable_table(
            frame=st.session_state[audit_key],
            download_name="inferential_independence_blocks",
            height=400,
        )
    with st.expander("Run metadata", expanded=False):
        st.json(metadata)


def _render_control_matching_coverage(*, database: Path) -> None:
    """Show pooled matching coverage without treating it as a test denominator.

    Args:
        database: Verified result database.
    """
    coverage = _matched_control_coverage(database=database)
    st.subheader(body="Pooled matched-control coverage")
    st.caption(
        body="Fraction of requested target blocks with at least one matched control. "
        "Unmatched targets are excluded from matched comparisons; review selection bias. "
        "Pools can serve several comparisons, so these totals are not test denominators."
    )
    _render_downloadable_table(frame=coverage, download_name="pooled_control_matching_coverage")
    if not coverage.empty:
        figure = px.bar(
            data_frame=coverage,
            x="background_label_id",
            y="target_coverage_fraction",
            labels={
                "background_label_id": "Control pool",
                "target_coverage_fraction": "Matched target fraction",
            },
            range_y=[0, 1],
            hover_data=["target_units", "covered_target_units", "matched_control_units"],
        )
        _render_plotly_figure(figure=figure, download_name="control_matching_coverage")


def _render_label_evidence_audit(*, database: Path) -> None:
    """Load only the selected, explicitly bounded label-audit preview.

    Args:
        database: Verified result database.
    """
    queries = {
        "Decisions": ("label_evidence_audit", "protein_id, label_id, rule_id"),
        "Matched controls": (
            "control_matching_audit",
            "background_label_id, target_unit_id, control_unit_id",
        ),
        "Excluded label features": (
            "label_definition_features",
            "label_id, feature_type, feature_id",
        ),
        "Abstentions": ("unresolved_assignments", "curation_status, protein_id"),
    }
    selected = st.selectbox(
        label="Evidence audit to preview",
        options=tuple(queries),
        help=(
            "Only this audit is queried. The preview and its downloads contain at most 5,000 rows."
        ),
    )
    table_name, ordering = queries[selected]
    st.caption(
        body="Preview: at most 5,000 rows. For the complete audit, use its published "
        "TSV/Parquet files in Canonical data & downloads."
    )
    audit_key = f"bounded-evidence-audit-{database}-{table_name}"
    if st.button(label="Load selected detailed audit"):
        st.session_state[audit_key] = query_dataframe(
            database=database, sql=f"SELECT * FROM {table_name} ORDER BY {ordering} LIMIT 5000"
        )
    if audit_key in st.session_state:
        _render_downloadable_table(
            frame=st.session_state[audit_key],
            download_name=f"{table_name}_preview",
            height=420,
        )


if __name__ == "__main__":  # pragma: no cover
    main()
