"""Run and explore audited DIAMOND partitions in the protein workbench."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from protein_signatures.workbench.config import WorkbenchConfig
from protein_signatures.workbench.context import CONTEXT_FILE
from protein_signatures.workbench.jobs import artifact_path, list_jobs
from protein_signatures.workbench.subcluster_review import (
    compare_partition,
    reviewed_taxon_sets,
    taxon_overlap,
)

from .subcluster_network import MAX_NODES, load_network_members, network_figure, network_graphml
from .workspace_ui import execution_choice, job_label, open_research_tool, queue_job, show_progress


def partition_summary(*, config: WorkbenchConfig, job: dict) -> pd.DataFrame:
    """Read the verified partition inventory with exact original representatives."""
    path = artifact_path(config=config, job_id=job["job_id"], name="subclusters.tsv")
    return pd.read_csv(path, sep="\t", dtype={"subcluster_id": str, "representative_id": str})


def partition_choice(
    *,
    config: WorkbenchConfig,
    job: dict,
    optional: bool,
) -> str:
    """Choose a partition without exposing an unbounded browser widget."""
    summary = partition_summary(config=config, job=job)
    key = f"wb-partition-{job['job_id']}" if optional else "wb-subcluster-review-id"
    query = st.text_input("Filter subclusters by ID or centroid", key=f"{key}-search").strip()
    if query:
        keep = summary.subcluster_id.str.contains(
            query, case=False, regex=False
        ) | summary.representative_id.str.contains(query, case=False, regex=False)
        visible = summary.loc[keep].head(200)
    else:
        visible = summary.head(200)
    previous = st.session_state.get(key)
    if previous and not query and previous not in set(visible.subcluster_id):
        visible = pd.concat([summary.loc[summary.subcluster_id.eq(previous)], visible])
    options = {
        row.subcluster_id: f"{row.representative_id} · {row.member_count:,} proteins · "
        f"{row.subcluster_id}"
        for row in visible.itertuples(index=False)
    }
    if optional:
        options = {"": "All saved input proteins", **options}
    if not options:
        st.info("No subcluster matches the search.")
        return ""
    return st.selectbox(
        "Audited subcluster",
        list(options),
        key=key,
        format_func=lambda value: options[value],
    )


def subcluster_review(*, config: WorkbenchConfig, job: dict) -> None:
    """Show complete-family comparisons and a bounded centroid network."""
    selected = partition_choice(config=config, job=job, optional=False)
    if not selected:
        return
    membership = artifact_path(config=config, job_id=job["job_id"], name="membership.parquet")
    context = (
        artifact_path(config=config, job_id=job["job_id"], name=CONTEXT_FILE)
        if CONTEXT_FILE in job["artifacts"]
        else None
    )
    comparisons = compare_partition(
        config=config,
        membership=membership,
        context=context,
        subcluster_id=selected,
    )
    st.caption("Comparisons include every saved input protein and every member of the selection.")
    for label, frame in comparisons.items():
        st.markdown(f"**{label}**")
        if frame.empty:
            st.caption("No recorded taxonomy context is available for this request.")
            continue
        st.dataframe(frame.head(1000), hide_index=True, width="stretch")
        name = f"subcluster_{label.lower().replace(' ', '_')}.tsv"
        st.download_button(
            f"Download {label.lower()}",
            frame.to_csv(sep="\t", index=False).encode(),
            name,
            on_click="ignore",
        )
    selection = {
        "wb-family-source": job_label(job=job, artifact="sequences.fasta"),
        f"wb-partition-{job['job_id']}": selected,
    }
    for tool in ("Start MAFFT alignment", "Predict annotations", "ESMFold & models"):
        st.button(
            f"Use this subcluster in {tool}",
            on_click=open_research_tool,
            kwargs={"tool": tool, "selection": selection},
            key=f"wb-partition-open-{tool}",
        )
    if not st.checkbox("Open interactive centroid network", key="wb-subcluster-network"):
        return
    summary = partition_summary(config=config, job=job)
    labels = dict(zip(summary.subcluster_id, summary.representative_id, strict=True))
    display = list(dict.fromkeys([selected, *summary.subcluster_id.head(200)]))
    groups = st.multiselect(
        "Subclusters to display",
        display,
        default=[selected],
        max_selections=12,
        format_func=lambda value: f"{labels[value]} · {value}",
        key="wb-network-groups",
    )
    if not groups:
        return
    members = load_network_members(path=membership, subclusters=groups)
    st.plotly_chart(network_figure(members=members), width="stretch")
    st.caption(
        f"Up to {MAX_NODES} displayed nodes, including singleton groups. Each line is a passing "
        "centroid-to-member alignment. There is no all-pairs edge evidence."
    )
    st.download_button(
        "Download displayed network (GraphML)",
        network_graphml(members=members),
        "audited_subcluster_network.graphml",
        on_click="ignore",
    )
    if context is not None:
        render_taxon_comparison(
            config=config, membership=membership, context=context, groups=groups
        )


def render_taxon_comparison(
    *,
    config: WorkbenchConfig,
    membership: Path,
    context: Path,
    groups: list[str],
) -> None:
    """Compare exact reviewed taxon sets for the chosen complete audited groups."""
    taxa = reviewed_taxon_sets(config=config, membership=membership, context=context, groups=groups)
    if not any(taxa.values()):
        st.caption("No reviewed NCBI taxon IDs are recorded in these subclusters.")
        return
    pairs, regions = taxon_overlap(groups=taxa)
    frame = pd.DataFrame(pairs)
    matrix = frame.pivot(index="subcluster_a", columns="subcluster_b", values="jaccard")
    st.plotly_chart(
        go.Figure(
            go.Heatmap(
                x=list(matrix.columns),
                y=list(matrix.index),
                z=matrix.to_numpy(),
                zmin=0,
                zmax=1,
                colorbar={"title": "Jaccard"},
                colorscale="Viridis",
            )
        ),
        width="stretch",
    )
    st.caption(
        "Overlap uses reviewed taxon IDs from complete groups; unreviewed names are excluded."
    )
    st.download_button(
        "Download exact taxon overlap counts",
        frame.to_csv(sep="\t", index=False).encode(),
        "reviewed_taxon_jaccard.tsv",
        on_click="ignore",
    )
    if regions:
        frame = pd.DataFrame(regions)
        st.dataframe(frame, hide_index=True)
        st.download_button(
            "Download exclusive overlap regions",
            frame.to_csv(sep="\t", index=False).encode(),
            "reviewed_taxon_overlap_regions.tsv",
            on_click="ignore",
        )


def render_subclustering(
    *,
    config: WorkbenchConfig,
    config_path: Path | None,
    uploads: dict,
    sources: dict,
) -> None:
    """Submit a stricter DeepClust run or inspect a completed audited partition."""
    st.caption(
        "DIAMOND DeepClust proposes centroid groups. Exact realignment checks identity, "
        "both sequence coverages and significance. Failed or missing links become singletons."
    )
    with st.expander("Start a new audited subclustering job", expanded=True):
        identity = st.number_input("Minimum exact identity (%)", 1.0, 100.0, 70.0, step=1.0)
        coverage = st.number_input(
            "Minimum coverage of both proteins (%)", 1.0, 100.0, 80.0, step=1.0
        )
        evalue = st.text_input("Maximum E-value", "1e-10")
        parameters = {
            "identity": identity,
            "mutual_coverage": coverage,
            "evalue": evalue,
            "transport": "local",
        }
        if config.subcluster_ssh_host and st.checkbox("Use the configured SSH analysis host"):
            parameters["transport"] = "ssh"
            parameters["remote_allowed"] = st.checkbox(
                f"Send the selected protein sequences to {config.subcluster_ssh_host}",
            )
        execution = execution_choice(config_path=config_path, key="wb-subcluster-execution")
        permitted = "family" in (uploads | sources) and (
            parameters["transport"] == "local" or parameters.get("remote_allowed") is True
        )
        if st.button("Start audited subclustering", disabled=not permitted):
            st.session_state["wb-subcluster-last"] = queue_job(
                config=config,
                config_path=config_path,
                operation="subcluster",
                execution=execution,
                parameters=parameters,
                uploads=uploads,
                sources=sources,
            )
        show_progress(config=config, key="wb-subcluster-last")
    jobs = {
        job_label(job=row): row
        for row in list_jobs(config=config)
        if row["operation"] == "subcluster" and row["status"] == "COMPLETE"
    }
    if jobs:
        chosen = st.selectbox(
            "Completed audited subclustering", list(jobs), key="wb-subcluster-job"
        )
        subcluster_review(config=config, job=jobs[chosen])
    else:
        st.info("No completed audited subclustering jobs yet.")
