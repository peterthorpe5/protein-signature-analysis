"""HOG alignment submission and exploration of immutable saved alignments."""

from __future__ import annotations

import csv
import hashlib
import io
import math
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.alignment_study import PAIR_CELLS, PAIR_ROWS, pairwise_identity
from protein_signatures.workbench.config import WorkbenchConfig
from protein_signatures.workbench.context import CONTEXT_FILE
from protein_signatures.workbench.jobs import list_jobs
from protein_signatures.workbench.species_tree import species_display_tree

from .hog_ui import released_inputs
from .workbench_viewers import tree_figure
from .workspace_ui import (
    execution_choice,
    job_bytes,
    job_label,
    open_research_tool,
    queue_job,
    show_progress,
)


def render_alignment_submit(
    *,
    config: WorkbenchConfig,
    config_path: Path | None,
    uploads: dict,
    sources: dict,
    key: str,
    allow_id_subset: bool = True,
    button_label: str = "Start MAFFT alignment",
) -> None:
    """Choose an exact reference and submit MAFFT, with an optional FastTree step."""
    st.caption(
        f"MAFFT limits: {config.max_alignment_sequences:,} proteins and "
        f"{config.max_alignment_cells:,} alignment cells."
    )
    selected = []
    if allow_id_subset:
        text = st.text_area(
            "Optional exact sequence IDs to align, one per line",
            key=f"{key}-ids",
            help="Leave blank for the complete input within the configured limits.",
        )
        selected = [line.strip() for line in text.splitlines() if line.strip()]
    reference = st.text_input(
        "Optional exact reference protein ID",
        key=f"{key}-ref",
        help="Controls the reference statistics. Leave blank to use the first ID in sorted order.",
    ).strip()
    tree = st.checkbox("Infer a FastTree WAG+gamma tree", key=f"{key}-tree")
    execution = execution_choice(config_path=config_path, key=f"{key}-execution")
    if st.button(button_label, key=f"{key}-submit", disabled="family" not in (uploads | sources)):
        st.session_state[f"{key}-last"] = queue_job(
            config=config,
            config_path=config_path,
            operation="align",
            execution=execution,
            parameters={"selected_ids": selected, "reference_id": reference, "tree": tree},
            uploads=uploads,
            sources=sources,
        )
    show_progress(config=config, key=f"{key}-last", tool="JBrowse MSA")


def context_for_alignment(*, data: bytes, records: dict[str, str]) -> dict[str, dict]:
    """Read species evidence only when the ungapped alignment row matches its saved hash."""
    result = {}
    rows = csv.DictReader(io.StringIO(data.decode()), delimiter="\t")
    for row in rows:
        identifier = row["sequence_id"]
        sequence = records.get(identifier)
        if sequence is None:
            continue
        digest = hashlib.sha256(sequence.replace("-", "").encode()).hexdigest()
        if row["sequence_sha256"] == digest:
            if identifier in result:
                raise InputValidationError(
                    "Alignment context contains a duplicate matching protein."
                )
            result[identifier] = row
    return result


def identity_figure(*, selected: list[str], pairs: list[dict]) -> go.Figure:
    """Plot exact pair identities with comparable counts and undefined cells."""
    indexes = {identifier: i for i, identifier in enumerate(selected)}
    shape = (len(selected), len(selected))
    values, overlap = np.full(shape, np.nan), np.zeros(shape, dtype=int)
    for pair in pairs:
        i, j = indexes[pair["sequence_a"]], indexes[pair["sequence_b"]]
        overlap[i, j] = overlap[j, i] = pair["comparable_residue_count"]
        if pair["identity"] is not None:
            values[i, j] = values[j, i] = pair["identity"]
    labels = [escape(value) for value in selected]
    return go.Figure(
        go.Heatmap(
            x=labels,
            y=labels,
            z=values,
            customdata=overlap,
            zmin=0,
            zmax=1,
            colorscale="Viridis",
            colorbar={"title": "Identity"},
            hovertemplate=(
                "%{y} × %{x}<br>Identity: %{z:.1%}<br>Comparable: %{customdata}<extra></extra>"
            ),
        )
    ).update_yaxes(autorange="reversed")


def render_alignment_statistics(
    *,
    config: WorkbenchConfig,
    job: dict,
    records: dict[str, str],
) -> None:
    """Display saved reference metrics and a bounded selected pairwise identity matrix."""
    if not st.checkbox("Inspect alignment statistics and identity", key="wb-jb-stats"):
        return
    st.caption(
        f"{len(records):,} aligned proteins · {len(next(iter(records.values()))):,} columns. "
        "Identity uses only positions containing a canonical amino acid in both rows. "
        "The raw BLOSUM62 sum has no gap penalty."
    )
    for artifact, label in (
        ("alignment_members.tsv", "Reference-to-member statistics"),
        ("alignment_columns.tsv", "Alignment column statistics"),
    ):
        if artifact in job["artifacts"]:
            payload = job_bytes(config=config, job=job, name=artifact)
            st.download_button(label, payload, artifact, on_click="ignore", key=f"wb-jb-{artifact}")
            st.dataframe(
                pd.read_csv(io.BytesIO(payload), sep="\t", nrows=1000),
                width="stretch",
                hide_index=True,
            )
    width = len(next(iter(records.values())))
    maximum = min(PAIR_ROWS, len(records), max(1, math.isqrt(2 * PAIR_CELLS // width) - 1))
    selected = st.multiselect(
        "Aligned proteins for the identity heatmap",
        options=list(records),
        default=list(records)[: min(maximum, 40)],
        max_selections=maximum,
        key=f"wb-jb-identity-{job['job_id']}",
    )
    if not selected:
        return
    pairs = pairwise_identity(records={identifier: records[identifier] for identifier in selected})
    st.plotly_chart(identity_figure(selected=selected, pairs=pairs), width="stretch")
    st.caption("A blank cell has no comparable residues. The TSV gives the exact denominators.")
    st.download_button(
        "Download displayed pairwise counts",
        pd.DataFrame(pairs).to_csv(sep="\t", index=False).encode(),
        "pairwise_identity.tsv",
        on_click="ignore",
    )
    if "pairwise_identity.tsv" in job["artifacts"] and len(selected) != len(records):
        st.download_button(
            "Download all saved pairwise counts",
            job_bytes(config=config, job=job, name="pairwise_identity.tsv"),
            "complete_pairwise_identity.tsv",
            on_click="ignore",
        )


def saved_tree(
    *,
    config: WorkbenchConfig,
    config_path: Path | None,
    database: Path | None,
    job: dict,
    records: dict[str, str],
) -> tuple[str, dict[str, str] | None]:
    """Select a matching saved FastTree result or a reviewed display backbone."""
    sources = {"alignment": {"job_id": job["job_id"], "name": "alignment.fasta"}}
    context_upload = {}
    evidence = {}
    if CONTEXT_FILE in job["artifacts"]:
        sources["context"] = {"job_id": job["job_id"], "name": CONTEXT_FILE}
        evidence = context_for_alignment(
            data=job_bytes(config=config, job=job, name=CONTEXT_FILE),
            records=records,
        )
    elif database is not None and st.checkbox("Recover HOG context for this saved alignment"):
        supplied, _ = released_inputs(database=database, config=config, key="wb-legacy-alignment")
        if "context" in supplied:
            evidence = context_for_alignment(data=supplied["context"], records=records)
            if evidence:
                context_upload["context"] = supplied["context"]
            st.caption(f"{len(evidence):,} of {len(records):,} exact IDs and sequences matched.")
    if evidence:
        with st.expander("Aligned member species and reviewed taxon IDs"):
            frame = pd.DataFrame(evidence.values())
            st.dataframe(frame.drop(columns=["sequence_sha256"]).head(1000), hide_index=True)
    related = [
        row
        for row in list_jobs(config=config)
        if (
            row["status"] == "COMPLETE"
            and row["operation"] == "tree"
            and row.get("input_jobs", {}).get("alignment") == job["job_id"]
        )
    ]
    if "tree.nwk" in job["artifacts"]:
        related.insert(0, job)
    options = {job_label(job=row): row for row in related}
    selected = st.selectbox(
        "Tree for this alignment",
        [*options, "No tree", "Reviewed species backbone"],
        key=f"wb-jb-tree-choice-{job['job_id']}",
    )
    if len(records) >= 3:
        with st.expander("Infer FastTree from this saved alignment"):
            st.caption("Uses the exact saved matrix and leaves MAFFT results as they are.")
            execution = execution_choice(config_path=config_path, key="wb-tree-execution")
            if st.button("Start FastTree from saved alignment"):
                st.session_state["wb-tree-last"] = queue_job(
                    config=config,
                    config_path=config_path,
                    operation="tree",
                    execution=execution,
                    sources=sources,
                    uploads=context_upload,
                )
            show_progress(config=config, key="wb-tree-last")
    if selected in options:
        source = options[selected]
        newick = job_bytes(config=config, job=source, name="tree.nwk").decode()
        rows = csv.DictReader(
            io.StringIO(job_bytes(config=config, job=source, name="tip_labels.tsv").decode()),
            delimiter="\t",
        )
        labels = {}
        for row in rows:
            if row["tree_tip"] in labels:
                raise InputValidationError("Duplicate FastTree tip mapping.")
            labels[row["tree_tip"]] = row["sequence_id"]
        if len(labels) != len(records) or set(labels.values()) != records.keys():
            raise InputValidationError(
                "The selected tree does not match this alignment's proteins."
            )
        if "auspice.json" in source["artifacts"]:
            st.download_button(
                "Download Auspice divergence JSON",
                job_bytes(config=config, job=source, name="auspice.json"),
                "auspice.json",
                on_click="ignore",
            )
        name = "protein_tree.nwk"
    elif selected == "Reviewed species backbone":
        backbone = st.file_uploader(
            "Reviewed species Newick (exact names or reviewed taxon IDs)",
            type=["nwk", "tree", "txt"],
            key="wb-jb-backbone",
        )
        if backbone is None:
            return "", None
        if backbone.size > 2_000_000:
            raise InputValidationError("The reviewed Newick must be at most 2 MB.")
        newick, audit = species_display_tree(
            newick=backbone.getvalue().decode(),
            sequence_ids=list(records),
            metadata=evidence,
        )
        labels, name = None, "species_backbone_display.nwk"
        st.download_button(
            "Download exact species mapping",
            pd.DataFrame(audit).to_csv(sep="\t", index=False).encode(),
            "species_backbone_mapping.tsv",
            on_click="ignore",
        )
        st.caption("Zero-length links group same-species proteins for display.")
    else:
        return "", None
    st.download_button("Download Newick", newick, name, on_click="ignore")
    if st.checkbox("Preview interactive tree", key="wb-jb-tree-preview"):
        st.plotly_chart(
            tree_figure(newick=newick, labels=labels or {name: name for name in records}),
            width="stretch",
        )
    return newick, labels


def alignment_handoffs(*, config: WorkbenchConfig, job: dict) -> None:
    """Offer the exact saved unaligned inputs for each downstream family tool."""
    name = "sequences.fasta"
    if name not in job["artifacts"]:
        return
    selected = {"wb-family-source": job_label(job=job, artifact=name)}
    for tool in ("Audited subclustering", "Predict annotations", "ESMFold & models"):
        st.button(
            f"Use unaligned input in {tool}",
            on_click=open_research_tool,
            kwargs={"tool": tool, "selection": selected},
            key=f"wb-jb-handoff-{tool}",
        )
    st.download_button(
        "Download exact unaligned inputs",
        job_bytes(config=config, job=job, name=name),
        "sequences.fasta",
        on_click="ignore",
    )
