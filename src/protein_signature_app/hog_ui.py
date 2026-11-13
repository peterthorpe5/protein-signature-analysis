"""Published HOG search, exact member selection and downstream snapshots."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.config import WorkbenchConfig
from protein_signatures.workbench.hogs import HogScope, hog_members, hog_snapshot, search_hogs

from .workspace_ui import open_research_tool


def choose_hog(*, database: Path, key: str) -> HogScope | None:
    """Choose a HOG from a literal search, with its run and hierarchy shown."""
    search = st.text_input("Search published HOG IDs", key=f"{key}-search", max_chars=200)
    scopes = search_hogs(database=database, text=search)
    if not scopes:
        st.info("No published HOG matches that search.")
        return None
    return st.selectbox(
        "Published HOG",
        scopes,
        index=None,
        key=f"{key}-hog",
        format_func=lambda scope: f"{scope.group_id} · {scope.run_id} · {scope.hierarchy_node}",
        placeholder="Choose a released HOG",
        help="Shows up to 50 HOGs. Search by ID to find a specific family.",
    )


def choose_hog_members(
    *,
    database: Path,
    config: WorkbenchConfig,
    key: str,
) -> tuple[HogScope | None, pd.DataFrame | None, bool]:
    """Choose the available released family or enter an explicit set of exact member IDs."""
    scope = choose_hog(database=database, key=key)
    if scope is None:
        return None, None, False
    subset = (
        st.radio(
            "Released members",
            ("All available released sequences", "Exact member subset"),
            horizontal=True,
            key=f"{key}-scope",
        )
        == "Exact member subset"
    )
    values = None
    if subset:
        text = st.text_area("Exact member IDs, one per line", key=f"{key}-ids")
        values = [line.strip() for line in text.splitlines() if line.strip()]
        if not values:
            return scope, None, subset
    frame = hog_members(
        database=database,
        scope=scope,
        maximum=config.max_seed_sequences,
        selected_ids=values,
    )
    count = int(frame.sequence.ne("").sum())
    st.caption(f"{len(frame):,} selected published members · {count:,} with an available sequence")
    if st.checkbox("Show selected member IDs and species", key=f"{key}-members"):
        visible = frame[["protein_id", "species_label"]].assign(
            sequence_available=frame.sequence.ne(""),
        )
        st.dataframe(visible.head(1000), width="stretch", hide_index=True)
    if subset and count != len(frame):
        raise InputValidationError("A selected member has no released sequence. Adjust the subset.")
    return scope, frame, subset


def released_inputs(
    *,
    database: Path,
    config: WorkbenchConfig,
    key: str,
) -> tuple[dict, dict]:
    """Return the exact released FASTA and sequence context for a new request."""
    scope, frame, _ = choose_hog_members(database=database, config=config, key=key)
    if scope is None or frame is None or not frame.sequence.ne("").any():
        return {}, {}
    fasta, context = hog_snapshot(config=config, scope=scope, members=frame)
    return {"family": fasta, "context": context}, {}


def send_hog_to_tool(*, tool: str, fasta: bytes, context: bytes) -> None:
    """Keep the selected released snapshot separate from the user's original search seeds."""
    st.session_state["wb-hog-snapshot"] = (fasta, context)
    open_research_tool(tool=tool, selection={"wb-family-source": "Current released HOG snapshot"})


def open_published_hog(*, scope: HogScope, tool: str) -> None:
    """Open the exact published scope from the canonical OrthoFinder page."""
    st.session_state["signature-app-page"] = "Protein research workbench"
    if tool == "JBrowse MSA":
        selection = {
            "wb-jb-source": "Align a published HOG",
            "wb-jb-hog-search": scope.group_id,
            "wb-jb-hog-hog": scope,
        }
    else:
        selection = {
            "wb-feature-source": "Published HOG features",
            "wb-feature-hog-search": scope.group_id,
            "wb-feature-hog-hog": scope,
        }
    open_research_tool(tool=tool, selection=selection)


def use_hog_as_seeds(*, fasta: bytes, context: bytes) -> None:
    """Use the explicitly chosen HOG as seeds for installed DeepClust discovery."""
    st.session_state.update(
        {
            "wb-seeds": fasta,
            "wb-seeds-context": context,
            "wb-stage": "2. Find families",
        }
    )
    st.session_state.pop("wb-map-ids", None)
