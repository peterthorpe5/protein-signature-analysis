"""Guided protein-family discovery, server jobs and bounded result inspection."""

from __future__ import annotations

import hashlib
import json
import re
import shlex
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from protein_signatures.errors import ProteinSignatureError
from protein_signatures.workbench.config import WorkbenchConfig, diagnostics, load_config
from protein_signatures.workbench.context import CONTEXT_FILE
from protein_signatures.workbench.jobs import (
    MAX_UPLOAD_BYTES,
    artifact_path,
    cancel_job,
    job_directory,
    list_jobs,
    recover_interrupted_job,
)
from protein_signatures.workbench.local_worker import start_local_job
from protein_signatures.workbench.seeds import export_seeds

from .alignment_ui import render_alignment_submit
from .hog_ui import released_inputs
from .research_tools import (
    portable_model_request,
    render_alignment_viewer,
    render_deepclust_import,
    render_feature_viewer,
    render_imported_model,
    render_predicted_models,
    render_structure_jobs,
)
from .subcluster_ui import partition_choice, render_subclustering
from .workbench_viewers import alignment_window, tree_figure
from .workspace_ui import (
    execution_choice,
    job_label,
    open_research_tool,
    queue_job,
    render_storage,
    show_progress,
)

WORKBENCH_HELP = """Choose seed proteins, map exact DeepClust membership or search global
representatives, retrieve complete families, and check each member against your seeds.
Then inspect conservation, infer a tree, scan local Pfam domains or apply frozen signatures.
Jobs run on the analysis server and survive closing the browser. New relatives remain
unreviewed candidates; global cluster membership does not establish orthology or function."""

GUIDE = """1. **Choose seeds:** upload/paste protein FASTA, or export exact proteins, local
   HOG members, profile labels or feature members from a completed result.
2. **Find families:** exact mapping finds existing membership. Representative search uses
   amino-acid similarity and also works when a seed accession is absent from DeepClust.
3. **Retrieve sequences:** select exact global centroid IDs from a completed mapping/search
   job. Retrieval checks complete member counts against the published row-group index.
4. **Analyse family:** screen every member against the original seeds. Select a manageable
   set for alignment/tree inference, local Pfam annotation or frozen-signature inspection.
5. **Prepare a new campaign:** export seeds plus selected relatives, deduplicated by exact
   sequence, with provenance and a curation template. Review labels, select controls and
   recompute homology blocks before a new statistical analysis.

**Keep the evidence separate.** A local HOG is not a global DeepClust cluster ID. An exact
member mapping is not a similarity search. A homologous protein is not automatically an
orthologue, a functional positive or an independent validation observation. Unknown domain,
taxonomy or structure evidence stays unknown. Additional relatives need not add independent
sample size. Search hit limits and thresholds are recorded; this is not an exhaustive claim
to recover every related sequence across the tree of life.
"""


def parse_ids(*, text: str) -> list[str]:
    """Read one exact ID per line, preserving source prefixes and accession versions.

    Args:
        text: Newline-delimited identifier text.

    Returns:
        Non-empty stripped lines.
    """
    return [value.strip() for value in text.splitlines() if value.strip()]


def render_workbench(*, config_path: Path | None, database: Path | None = None) -> None:
    """Render a guided workbench independently of an optional reference campaign.

    Args:
        config_path: Launch-time operator configuration.
        database: Optional verified reference database for seed selection.
    """
    st.title("Protein research workbench")
    st.write(
        "Expand a protein family using DeepClust, then investigate its sequences and signatures."
    )
    if config_path is None:
        st.info(
            (
                "Launch with --workbench-config pointing to the server profile to "
                "enable discovery and analysis jobs. See docs/DENBI_WORKBENCH.md for "
                "setup and the Mac launcher."
            )
        )
        with st.expander("? Workflow, inputs and scientific scope", expanded=True):
            st.markdown(GUIDE)
        return
    try:
        config = load_config(path=config_path)
        with st.expander("? Start here: what goes in, what comes out", expanded=False):
            st.markdown(GUIDE)
        stage = st.radio(
            label="Workflow step",
            options=(
                "Start here",
                "1. Choose seeds",
                "2. Find families",
                "3. Retrieve sequences",
                "4. Analyse family",
                "5. Research tools",
                "Jobs & results",
            ),
            horizontal=True,
            key="wb-stage",
        )
        if stage == "Start here":
            st.markdown(GUIDE)
            st.dataframe(
                data=pd.DataFrame(data=diagnostics(config=config)), hide_index=True, width="stretch"
            )
            st.caption(
                (
                    "Unavailable optional tools affect only the tasks that need them. "
                    "Mapping and retrieval do not require a representative search database."
                )
            )
            with st.expander("Operator details and representative database setup"):
                st.code(
                    body=(
                        f"Profile: {config_path}\nWorkspace: {config.workspace}\n"
                        f"Release: {config.release}"
                    )
                )
                st.write(
                    (
                        "If a representative .dmnd database is already installed, set "
                        "centroid_database in the profile. Otherwise set centroid_fasta to the "
                        "release's representative FASTA and queue its database build below. "
                        "This can use substantial disk space and time."
                    )
                )
                if st.button(
                    label="Queue representative database build",
                    disabled=not bool(config.centroid_fasta),
                ):
                    queue_job(config=config, operation="build_search")
                st.caption(
                    (
                        "After completion, set centroid_database to the resulting "
                        "representatives.dmnd path, then submit new searches."
                    )
                )
        elif stage == "1. Choose seeds":
            render_seeds(config=config, database=database)
        elif stage == "2. Find families":
            render_discovery(config=config)
        elif stage == "3. Retrieve sequences":
            render_retrieval(config=config)
        elif stage == "4. Analyse family":
            render_analysis(config=config, config_path=config_path, database=database)
        elif stage == "5. Research tools":
            render_research_tools(config=config, config_path=config_path, database=database)
        else:
            render_jobs(config=config, config_path=config_path)
    except (ProteinSignatureError, OSError, ValueError, KeyError) as error:
        st.error(str(error))


def render_seeds(*, config: WorkbenchConfig, database: Path | None) -> None:
    """Collect a seed FASTA snapshot from uploads or exact canonical selectors.

    Args:
        config: Submission bounds.
        database: Optional verified reference result.
    """
    st.subheader("1 · Choose the proteins you want to expand")
    st.caption(
        (
            "Start with one small family. Seeds are your query sequences, not "
            "automatic functional labels."
        )
    )
    route = st.selectbox(
        label="Seed source",
        options=(
            "Upload or paste protein FASTA",
            "Select from the completed result",
            *(("Search published HOGs",) if database is not None else ()),
        ),
        key="wb-seed-route",
    )
    if route == "Upload or paste protein FASTA":
        upload = st.file_uploader(
            label="Unaligned protein FASTA (maximum 20 MiB)",
            type=["fa", "faa", "fasta"],
            key="wb-seed-upload",
            help=(
                "Headers start with > followed by unique protein IDs. Use amino-acid "
                "sequences, with no alignment gaps. Example: >protein_A then MAAKLVV on"
                " the next line."
            ),
        )
        pasted = st.text_area(
            label="Or paste protein FASTA", key="wb-seed-paste", placeholder=">protein_A\nMAAKLVV"
        )
        if st.button(label="Use these seed proteins"):
            payload = upload.getvalue() if upload is not None else pasted.encode()
            if not payload or len(payload) > MAX_UPLOAD_BYTES:
                raise ValueError("Supply a protein FASTA of at most 20 MiB.")
            st.session_state["wb-seeds"] = payload
            st.session_state.pop("wb-seeds-context", None)
            st.success("Seed snapshot selected. Continue to Find families.")
    elif route == "Search published HOGs":
        uploads, _ = released_inputs(database=database, config=config, key="wb-seed-hog")
        if st.button("Use these released HOG members as seeds", disabled="family" not in uploads):
            st.session_state["wb-seeds"] = uploads["family"]
            st.session_state["wb-seeds-context"] = uploads["context"]
            st.success("Released HOG seed snapshot selected. Continue to Find families.")
    elif database is None:
        st.info(
            (
                "Open a completed reference result with --resource to use its proteins,"
                " HOGs, labels or features. You can upload seed FASTA without a "
                "reference result."
            )
        )
    else:
        selection = st.selectbox(
            label="Select proteins by",
            options=("protein", "group", "label", "feature"),
            format_func=lambda value: {
                "protein": "Exact protein IDs",
                "group": "Local HOG / orthogroup IDs",
                "label": "Profile label IDs",
                "feature": "Signature feature IDs",
            }[value],
            key="wb-selector",
            help=(
                "Uses the current completed result. HOGs select their member proteins; "
                "they are mapped to global DeepClust IDs in the next step."
            ),
        )
        values = st.text_area(label="Exact selector IDs, one per line", key="wb-selector-ids")
        maximum = st.number_input(
            label="Maximum complete seed count",
            min_value=1,
            max_value=min(config.max_seed_sequences, 10000),
            value=min(1000, config.max_seed_sequences),
            step=1,
            help="An oversized selection is rejected rather than silently sampled.",
        )
        if st.button(label="Prepare seed FASTA from this selection"):
            st.session_state["wb-seeds"] = export_seeds(
                database=database,
                selection=selection,
                values=parse_ids(text=values),
                maximum=int(maximum),
            )
            st.session_state.pop("wb-seeds-context", None)
            st.success("Complete seed set prepared. Continue to Find families.")
    if payload := st.session_state.get("wb-seeds"):
        count = sum(line.startswith(b">") for line in payload.splitlines())
        st.info(
            f"Current seed snapshot: {count} FASTA headers; {len(payload):,} bytes. "
            "The worker validates IDs and protein sequences before analysis."
        )
        st.download_button(
            label="Download selected seeds",
            data=payload,
            file_name="seed_proteins.fasta",
            mime="text/plain",
        )
        st.button(
            label="Open ESMFold for these seeds",
            on_click=open_research_tool,
            kwargs={
                "tool": "ESMFold & models",
                "selection": {"wb-family-source": "Current seed snapshot"},
            },
        )


def similarity_controls(*, key: str) -> dict:
    """Render identity, bilateral coverage and E-value controls with biological help.

    Args:
        key: Widget namespace.

    Returns:
        Named numerical thresholds.
    """
    identity = st.number_input(
        label="Minimum identity (%)",
        min_value=0.0,
        max_value=100.0,
        value=30.0,
        key=f"{key}-identity",
        help=(
            "Percent identical residues within the alignment. A starting point, not"
            " a universal boundary for shared function."
        ),
    )
    coverage = st.number_input(
        label="Minimum coverage of both proteins (%)",
        min_value=0.0,
        max_value=100.0,
        value=50.0,
        key=f"{key}-coverage",
        help=(
            "Both the query and matched protein must meet this coverage. Shared "
            "short domains can fail full-length coverage; reducing it broadens the "
            "search."
        ),
    )
    evalue = st.text_input(
        label="Maximum E-value",
        value="1e-5",
        key=f"{key}-evalue",
        help=(
            "Smaller is stricter. E-values depend on database size and do not establish orthology."
        ),
    )
    return {"identity": float(identity), "coverage": float(coverage), "evalue": float(evalue)}


def render_discovery(*, config: WorkbenchConfig) -> None:
    """Submit exact mapping or sequence-similarity search against representatives.

    Args:
        config: Operator sources and task bounds.
    """
    st.subheader("2 · Find the global families related to your seeds")
    route = st.radio(
        label="Discovery route",
        options=("Map exact protein identifiers", "Search representative sequences"),
        key="wb-discovery-route",
    )
    if route == "Map exact protein identifiers":
        st.write(
            (
                "Maps known accessions to all matching global clusters. Unmatched IDs "
                "and small/unindexed clusters remain visible. This does not compare "
                "sequences."
            )
        )
        current = st.session_state.get("wb-seeds", b"").decode(encoding="utf-8", errors="replace")
        defaults = "\n".join(
            line[1:].split()[0]
            for line in current.splitlines()
            if line.startswith(">") and line[1:].split()
        )
        text = st.text_area(
            label="Protein IDs, one per line",
            value=defaults,
            key="wb-map-ids",
            help=(
                "Use exact source IDs or accessions, for example NR_UGC79169.1. "
                "Versions and prefixes are preserved; configured namespace aliases are "
                "recorded."
            ),
        )
        if st.button(label="Queue identifier mapping"):
            queue_job(
                config=config, operation="resolve", parameters={"identifiers": parse_ids(text=text)}
            )
    else:
        st.write(
            (
                "DIAMOND compares the selected seed amino-acid sequences with the "
                "configured global representatives. Inspect retained hits before "
                "selecting families for retrieval."
            )
        )
        parameters = similarity_controls(key="wb-search")
        maximum = st.number_input(
            label="Maximum retained representatives per seed",
            min_value=1,
            max_value=1000,
            value=50,
            help=(
                "This bounds output and search scope. Increasing it can recover more "
                "families; no-hit is not proof of absence."
            ),
        )
        if st.button(
            label="Queue representative search", disabled=not bool(config.centroid_database)
        ):
            seeds = st.session_state.get("wb-seeds")
            if not seeds:
                raise ValueError("Choose seed FASTA in step 1 first.")
            uploads = {"seeds": seeds}
            if context := st.session_state.get("wb-seeds-context"):
                uploads["context"] = context
            queue_job(
                config=config,
                operation="search",
                parameters={**parameters, "maximum_hits": int(maximum)},
                uploads=uploads,
            )
        if not config.centroid_database:
            st.info(
                (
                    "Configure the global representative .dmnd database in Start here to "
                    "enable homology searches. Exact membership mapping is available with "
                    "the mapping Parquet."
                )
            )


def render_retrieval(*, config: WorkbenchConfig) -> None:
    """Let users explicitly select complete global families from prior results or IDs.

    Args:
        config: Catalogue and job workspace.
    """
    st.subheader("3 · Retrieve complete families from the server database")
    st.write(
        (
            "Choose exact global centroid IDs. The worker reads the indexed "
            "physical Parquet row groups and verifies the complete published member"
            " count. No complete cluster is silently sampled."
        )
    )
    jobs = [
        job
        for job in list_jobs(config=config)
        if job["status"] == "COMPLETE" and job["operation"] in {"resolve", "search"}
    ]
    selected = []
    if jobs:
        options = {job_label(job=job): job for job in jobs}
        label = st.selectbox(
            label="Completed discovery job", options=list(options), key="wb-retrieve-job"
        )
        job = options[label]
        filename = "mapping.tsv" if job["operation"] == "resolve" else "representative_hits.tsv"
        path = artifact_path(config=config, job_id=job["job_id"], name=filename)
        preview = pd.read_csv(filepath_or_buffer=path, sep="\t", nrows=1000, dtype=str)
        st.dataframe(data=preview, hide_index=True, width="stretch")
        values = sorted(
            preview.loc[preview["member_count"].notna(), "cluster_id"].dropna().unique()
        )
        selected = st.multiselect(
            label="Global clusters to retrieve",
            options=values,
            key="wb-clusters",
            help=(
                "Only indexed clusters can be retrieved. Choose one family to begin. If"
                " more than 1000 hit rows exist, use the complete result file in Jobs "
                "to select additional exact IDs."
            ),
        )
    text = st.text_area(
        label="Or enter exact global centroid IDs, one per line",
        key="wb-exact-centroids",
        help=(
            "A local HOG ID cannot be used here. Global examples include "
            "NR_UGC79169.1. The official index omits clusters smaller than three "
            "members."
        ),
    )
    if st.button(label="Queue complete cluster retrieval"):
        queue_job(
            config=config,
            operation="retrieve",
            parameters={"cluster_ids": list(dict.fromkeys([*selected, *parse_ids(text=text)]))},
        )


def family_input(*, config: WorkbenchConfig, database: Path | None = None) -> tuple[dict, dict]:
    """Select an exact unaligned artifact, audited partition, released HOG or FASTA upload."""
    st.caption("Choose unaligned proteins. Aligned matrices cannot be used as family FASTA.")
    candidates = {}
    for job in list_jobs(config=config):
        if job["status"] == "COMPLETE":
            for name in job.get("artifacts", {}):
                if (
                    name
                    in {
                        "family.fasta",
                        "candidates.fasta",
                        "sequences.fasta",
                        "expanded.fasta",
                        "seeds.fasta",
                        "model_inputs.fasta",
                        "annotation_input.fasta",
                    }
                    and job["artifacts"][name]["bytes"]
                ):
                    candidates[job_label(job=job, artifact=name)] = (job, name)
    seed_snapshot = st.session_state.get("wb-seeds")
    hog_snapshot = st.session_state.get("wb-hog-snapshot")
    selected = st.selectbox(
        label="Family FASTA from a completed job",
        options=[
            "Upload instead",
            *(["Current seed snapshot"] if seed_snapshot else []),
            *(["Current released HOG snapshot"] if hog_snapshot else []),
            *(["Select a published HOG"] if database is not None else []),
            *candidates,
        ],
        key="wb-family-source",
    )
    if selected == "Current seed snapshot":
        uploads = {"family": seed_snapshot}
        if context := st.session_state.get("wb-seeds-context"):
            uploads["context"] = context
        return uploads, {}
    if selected == "Current released HOG snapshot":
        return {"family": hog_snapshot[0], "context": hog_snapshot[1]}, {}
    if selected == "Select a published HOG":
        return released_inputs(database=database, config=config, key="wb-family-hog")
    if selected in candidates:
        job, name = candidates[selected]
        identifier = job["job_id"]
        uploads, sources = {}, {"family": {"job_id": identifier, "name": name}}
        if CONTEXT_FILE in job["artifacts"]:
            sources["context"] = {"job_id": identifier, "name": CONTEXT_FILE}
        elif database is not None and st.checkbox(
            "Match this saved input to a released HOG for species context",
            key="wb-legacy-context",
        ):
            matched, _ = released_inputs(database=database, config=config, key="wb-legacy-family")
            if "context" in matched:
                uploads["context"] = matched["context"]
                st.caption("The worker retains metadata only for exact ID and sequence matches.")
        if job["operation"] == "subcluster" and name == "sequences.fasta":
            partition = partition_choice(config=config, job=job, optional=True)
            if partition:
                sources["subcluster_membership"] = {
                    "job_id": identifier,
                    "name": "membership.parquet",
                    "subcluster_id": partition,
                }
        return uploads, sources
    upload = st.file_uploader(
        label="Family protein FASTA (maximum 20 MiB)",
        type=["fa", "faa", "fasta"],
        key="wb-family-upload",
    )
    return ({"family": upload.getvalue()} if upload is not None else {}), {}


def render_analysis(
    *,
    config: WorkbenchConfig,
    config_path: Path | None = None,
    database: Path | None = None,
) -> None:
    """Render family screening, alignment, annotation, signature and expansion tasks.

    Args:
        config: Operator profile.
        config_path: Trusted launch profile for app-host jobs.
        database: Optional read-only result database for HOG selection.
    """
    st.subheader("4 · Analyse and prepare the expanded family")
    task = st.selectbox(
        label="Task",
        options=(
            "Screen members against seeds",
            "Alignment and tree",
            "Local Pfam domains",
            "Frozen signature evidence",
            "Prepare expanded campaign inputs",
        ),
        key="wb-analysis-task",
    )
    uploads, sources = family_input(config=config, database=database)
    parameters, operation = (
        {},
        {
            "Screen members against seeds": "screen",
            "Alignment and tree": "align",
            "Local Pfam domains": "annotate",
            "Frozen signature evidence": "score",
            "Prepare expanded campaign inputs": "prepare",
        }[task],
    )
    if operation == "screen":
        st.write(
            (
                "Recheck individual family members against your original seeds; review "
                "the full selection audit before interpreting the retained candidates."
            )
        )
        parameters = similarity_controls(key="wb-screen")
        parameters["deduplicate"] = st.checkbox(
            label="Retain one protein per identical amino-acid sequence",
            value=True,
            help=(
                "Keeps the lexicographically first eligible ID for each exact sequence."
                " The audit preserves excluded IDs."
            ),
        )
        text = st.text_area(
            label="Optional exact NCBI taxon IDs, one per line",
            disabled=not bool(config.taxonomy),
            help=(
                "Requires configured reviewed sequence taxonomy. Exact taxa only; "
                "descendants are not inferred. Unknown/unreviewed assignments are "
                "excluded when a taxon filter is active."
            ),
        )
        parameters["taxon_ids"] = parse_ids(text=text)
    elif operation == "align":
        render_alignment_submit(
            config=config,
            config_path=config_path,
            uploads=uploads,
            sources=sources,
            key="wb-family-mafft",
            button_label="Queue alignment and tree",
        )
        return
    elif operation == "annotate":
        st.write(
            (
                "Run local HMMER hmmscan against the configured pressed Pfam database "
                "using curated gathering thresholds. Produces domain intervals and "
                "positive query-feature evidence; no-hit records are explicit."
            )
        )
        if not config.pfam_hmm:
            st.info(
                "Set pfam_hmm in the server profile and prepare the database with hmmpress first."
            )
    elif operation == "score":
        st.write(
            (
                "Apply the completed campaign's frozen validated signatures to the "
                "selected proteins. Exact k-mers are assessed from sequence; "
                "unavailable structural/domain evidence remains unknown. This is "
                "candidate prioritisation, not external validation or automatic "
                "function prediction."
            )
        )
        text = st.text_area(
            label="Optional exact comparison IDs, one per line",
            help=(
                "Blank uses all reference comparisons. Scoring requires resource in the"
                " workbench profile."
            ),
        )
        parameters = {
            "comparison_ids": parse_ids(text=text),
            "include_within_comparison": st.checkbox(
                label="Also include within-comparison validated signatures", value=False
            ),
        }
        annotation_jobs = [
            job
            for job in list_jobs(config=config)
            if job["status"] == "COMPLETE" and job["operation"] == "annotate"
        ]
        options = {job["job_id"]: job for job in annotation_jobs}
        chosen = st.selectbox(
            label="Optional Pfam evidence from the same protein input",
            options=["No extra evidence", *options],
        )
        if chosen in options:
            sources["features"] = {"job_id": chosen, "name": "query_features.tsv"}
        if not config.resource:
            st.info(
                (
                    "Set resource in the profile to the immutable completed reference "
                    "result before scoring."
                )
            )
    else:
        st.write(
            (
                "Combine original seeds with selected candidates, deduplicate exact "
                "sequences and preserve every source ID. Produces expanded.fasta, a "
                "provenance ledger and an unreviewed label template. Review labels and "
                "controls, then recompute homology blocks across the expanded set "
                "before starting a new campaign."
            )
        )
    if st.button(label=f"Queue {task.lower()}"):
        if "family" not in (uploads | sources):
            raise ValueError("Select a family FASTA first.")
        if operation in {"screen", "prepare"}:
            seeds = st.session_state.get("wb-seeds")
            if not seeds:
                raise ValueError("Choose the original seed FASTA in step 1 first.")
            uploads["seeds"] = seeds
        queue_job(
            config=config,
            operation=operation,
            parameters=parameters,
            uploads=uploads,
            sources=sources,
        )


def render_research_tools(
    *,
    config: WorkbenchConfig,
    config_path: Path | None = None,
    database: Path | None = None,
) -> None:
    """Connect portable inputs and new evidence providers to the durable job queue."""
    tool = st.selectbox(
        "Research tool",
        options=(
            "JBrowse MSA",
            "Start MAFFT alignment",
            "Audited subclustering",
            "ESMFold & models",
            "Nightingale features",
            "Predict annotations",
            "LoL-align structures",
            "Import DeepClust results",
        ),
        key="wb-research-tool",
    )
    if tool == "JBrowse MSA":
        render_alignment_viewer(config=config, config_path=config_path, database=database)
        return
    if tool == "Nightingale features":
        render_feature_viewer(config=config, database=database)
        return
    if tool == "LoL-align structures":
        render_structure_jobs(config=config)
        return
    if tool == "Import DeepClust results":
        render_deepclust_import(config=config)
        return
    if tool in {"Start MAFFT alignment", "Audited subclustering"}:
        st.subheader(tool)
        uploads, sources = family_input(config=config, database=database)
        if tool == "Start MAFFT alignment":
            render_alignment_submit(
                config=config,
                config_path=config_path,
                uploads=uploads,
                sources=sources,
                key="wb-research-mafft",
            )
        else:
            render_subclustering(
                config=config,
                config_path=config_path,
                uploads=uploads,
                sources=sources,
            )
        return
    model = tool == "ESMFold & models"
    st.subheader(tool)
    if model:
        st.caption(
            f"Predict up to 20 models, one per exact sequence of at most "
            f"{config.esm_max_length} residues."
        )
        render_predicted_models(config=config)
        render_imported_model(config=config)
        if not config.esm_endpoint:
            st.info(
                "Configure esm_endpoint for predictions on this worker. "
                "You can still prepare an exact FASTA or portable request below."
            )
        engine = "esmfold"
    else:
        st.caption(
            "Predict features for proteins with little existing annotation. Local Pfam remains "
            "available in Analyse family; only that route supplies frozen-signature Pfam evidence."
        )
        engines = {
            "Local InterProScan": "interpro_local",
            "EBI InterProScan service": "interpro_remote",
            "Local eggNOG-mapper": "eggnog",
        }
        engine = engines[
            st.selectbox("Annotation provider", list(engines), key="wb-annotation-engine")
        ]
    uploads, sources = family_input(config=config, database=database)
    parameters = {
        "selected_ids": parse_ids(
            text=st.text_area(
                "Exact protein IDs, one per line (blank selects the whole bounded input)",
                key=f"wb-selected-{engine}",
            )
        ),
    }
    can_submit = "family" in (uploads | sources)
    if model and not config.esm_endpoint:
        can_submit = False
    elif model or engine == "interpro_remote":
        service = config.esm_endpoint if model else "the EBI InterProScan service"
        parameters["remote_allowed"] = st.checkbox(
            f"Send the selected amino-acid sequences to {service}",
            key=f"wb-provider-consent-{engine}",
        )
        can_submit = can_submit and parameters["remote_allowed"]
    if not model:
        parameters["engine"] = engine
    if engine == "interpro_remote":
        email = st.text_input(
            "Private EBI contact email",
            type="password",
            max_chars=254,
            help="Used only for this EBI request. The contact is kept out of result files.",
        ).strip()
        if email:
            valid = bool(re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", email)) and email.isascii()
            if valid:
                uploads["ebi_email"] = email.encode()
            else:
                st.info("Enter a valid email address.")
                can_submit = False
        elif not config.interpro_email:
            st.info("Enter an EBI contact for this request.")
            can_submit = False
    if engine == "eggnog":
        parameters["tax_scope"] = st.text_input("eggNOG taxonomy scope", value="auto")
        parameters["target_taxon"] = st.text_input("Optional exact target taxon ID")
        if not config.eggnog_data:
            st.info("Configure the local eggNOG database directory to enable this route.")
            can_submit = False
    execution = execution_choice(config_path=config_path, key=f"wb-execution-{engine}")
    if st.button("Queue selected proteins", disabled=not can_submit, key=f"wb-queue-{engine}"):
        st.session_state[f"wb-provider-last-{engine}"] = queue_job(
            config_path=config_path,
            execution=execution,
            config=config,
            operation="predict_models" if model else "annotate_proteins",
            parameters=parameters,
            uploads=uploads,
            sources=sources,
        )
    show_progress(config=config, key=f"wb-provider-last-{engine}")
    if model and "family" in (uploads | sources):
        fingerprint = hashlib.sha256(
            uploads.get("family", b"")
            + uploads.get("context", b"")
            + json.dumps(
                [sources, parameters["selected_ids"], config.esm_max_length], sort_keys=True
            ).encode()
        ).hexdigest()
        if st.button("Prepare exact FASTA and portable request", key="wb-model-prepare"):
            prepared = portable_model_request(
                config=config,
                selected=parameters["selected_ids"],
                uploads=uploads,
                sources=sources,
            )
            st.session_state["wb-model-request"] = (fingerprint, *prepared)
        ready = st.session_state.get("wb-model-request")
        if ready and ready[0] == fingerprint:
            st.download_button(
                "Download exact FASTA", ready[1], "model_inputs.fasta", on_click="ignore"
            )
            st.download_button(
                "Download portable ESMFold request",
                ready[2],
                "protein_signature_esmfold_request.zip",
                on_click="ignore",
            )


def render_jobs(*, config: WorkbenchConfig, config_path: Path | None = None) -> None:
    """Show durable progress, bounded previews and verified on-demand downloads.

    Args:
        config: Job workspace.
        config_path: Trusted launch profile for starting an exact waiting job.
    """
    st.subheader("Jobs & results")
    render_storage(config=config)
    st.button(label="Refresh job status")
    st.caption(
        (
            "QUEUED waits for the server worker; RUNNING is executing; COMPLETE has"
            " verified output files. FAILED retains diagnostics. Closing this page "
            "does not cancel work. An interrupted RUNNING job is never silently "
            "retried; inspect its PID and logs before submitting a replacement."
        )
    )
    jobs = list_jobs(config=config)
    if not jobs:
        st.info("No workbench jobs yet. Start by choosing seed proteins.")
        return
    st.dataframe(
        data=pd.DataFrame(
            data=[
                {
                    name: row.get(name, "")
                    for name in ("job_id", "operation", "status", "created_at", "error")
                }
                for row in jobs
            ]
        ),
        hide_index=True,
        width="stretch",
    )
    options = {job_label(job=job, with_status=True): job for job in jobs}
    selected = st.selectbox(label="Inspect job", options=list(options), key="wb-inspect-job")
    job = options[selected]
    if job["status"] == "QUEUED" and config_path is not None:
        if st.button("Start this queued job on the app host"):
            start_local_job(config=config, config_path=config_path, job_id=job["job_id"])
            st.success("App-host start requested. Refresh to see its state.")
        if job.get("execution") == "app":
            st.caption("This request is reserved for an exact app-host start.")
        with st.expander("Start this request from a server terminal"):
            st.code(
                f"protein-signature-workbench --config {shlex.quote(str(config_path))} "
                f"--action worker --once --job-id {job['job_id']}",
                language="bash",
            )
    if job["status"] == "QUEUED" and st.button(label="Cancel this waiting job"):
        cancel_job(config=config, job_id=job["job_id"])
        st.success("Waiting job cancelled.")
    if job["status"] == "FAILED":
        st.error(job.get("error", "Job failed."))
        path = job_directory(config=config, job_id=job["job_id"]) / "failure.log"
        if path.is_file():
            st.code(body=path.read_text(encoding="utf-8")[-12000:])
        return
    if job["status"] != "COMPLETE":
        if job["status"] == "RUNNING":
            st.caption(
                f"Worker: {job.get('worker_host')} · PID {job.get('worker_pid')} · "
                f"started {job.get('started_at')}"
            )
            if st.button(label="Mark a stopped, interrupted job as failed"):
                recover_interrupted_job(config=config, job_id=job["job_id"])
                st.success(
                    (
                        "Interrupted job marked failed. Its files are retained; submit a new "
                        "request to retry."
                    )
                )
        return
    route = {
        "align": ("JBrowse MSA", "wb-jb-source", "Completed MAFFT job", "wb-jb-job", "wb-jb-open"),
        "annotate": (
            "Nightingale features",
            "wb-feature-source",
            "Completed annotation job",
            "wb-feature-job",
            "wb-feature-open",
        ),
        "annotate_proteins": (
            "Nightingale features",
            "wb-feature-source",
            "Completed annotation job",
            "wb-feature-job",
            "wb-feature-open",
        ),
    }.get(job["operation"])
    if route:
        tool, source_key, source, job_key, open_key = route
        selected_label = job_label(job=job)
        st.button(
            label=f"Open in {tool}",
            on_click=open_research_tool,
            kwargs={
                "tool": tool,
                "selection": {source_key: source, job_key: selected_label, open_key: True},
            },
        )
    if job["operation"] == "subcluster":
        st.button(
            "Inspect audited subclusters",
            on_click=open_research_tool,
            kwargs={
                "tool": "Audited subclustering",
                "selection": {"wb-subcluster-job": job_label(job=job)},
            },
        )
    family = next(
        (
            name
            for name in ("candidates.fasta", "family.fasta", "sequences.fasta", "expanded.fasta")
            if name in job["artifacts"] and job["artifacts"][name]["bytes"]
        ),
        "",
    )
    if family:
        selection = {"wb-family-source": job_label(job=job, artifact=family)}
        st.button(
            label="Use this FASTA for audited subclustering",
            on_click=open_research_tool,
            kwargs={"tool": "Audited subclustering", "selection": selection},
        )
        st.button(
            label="Use this FASTA in ESMFold",
            on_click=open_research_tool,
            kwargs={"tool": "ESMFold & models", "selection": selection},
        )
        st.button(
            label="Use this FASTA for annotation",
            on_click=open_research_tool,
            kwargs={"tool": "Predict annotations", "selection": selection},
        )
    if job["operation"] == "align":
        with st.expander("Inspect the alignment and tree", expanded=False):
            alignment_path = artifact_path(
                config=config, job_id=job["job_id"], name="alignment.fasta"
            )
            alignment_text = alignment_path.read_text(encoding="utf-8")
            first = int(
                st.number_input(label="First alignment position to display", min_value=1, value=1)
            )
            st.code(body=alignment_window(text=alignment_text, start=first))
            st.caption(
                (
                    "Display: up to 50 proteins and 80 columns from the chosen position. "
                    "Complete alignment FASTA remains downloadable below; gaps are "
                    "alignment positions, not protein residue numbers."
                )
            )
            if "tree.nwk" in job.get("artifacts", {}) and st.checkbox(
                label="Show interactive phylogram"
            ):
                tree_path = artifact_path(config=config, job_id=job["job_id"], name="tree.nwk")
                map_path = artifact_path(config=config, job_id=job["job_id"], name="tip_labels.tsv")
                mapping = pd.read_csv(filepath_or_buffer=map_path, sep="\t", dtype=str)
                figure = tree_figure(
                    newick=tree_path.read_text(encoding="utf-8"),
                    labels=dict(zip(mapping.tree_tip, mapping.sequence_id, strict=True)),
                )
                st.plotly_chart(figure_or_data=figure, width="stretch")
                st.caption(
                    (
                        "FastTree's supplied display root is retained; no biological rooting or"
                        " time calibration is inferred. Related sequences do not establish "
                        "shared function."
                    )
                )
    names = sorted(job.get("artifacts", {}))
    if not names:
        return
    name = st.selectbox(label="Result file", options=names, key=f"wb-artifact-{job['job_id']}")
    record = job["artifacts"][name]
    st.caption(
        f"Complete file: {record['bytes']:,} bytes. "
        "Downloads are checked against their saved SHA-256."
    )
    if name.endswith(".tsv") and st.checkbox(
        label="Show table preview (first 1000 rows)", key=f"wb-preview-{job['job_id']}-{name}"
    ):
        path = artifact_path(config=config, job_id=job["job_id"], name=name)
        frame = pd.read_csv(filepath_or_buffer=path, sep="\t", nrows=1000)
        st.dataframe(data=frame, hide_index=True, width="stretch")
        if name == "conservation.tsv":
            figure = px.line(
                data_frame=frame,
                x="alignment_position",
                y=["consensus_fraction_non_gap", "gap_fraction"],
                labels={
                    "value": "Fraction",
                    "alignment_position": "Alignment position",
                    "variable": "Measure",
                },
                color_discrete_sequence=["#1768ac", "#d47800"],
            )
            st.plotly_chart(figure_or_data=figure, width="stretch")
    if name.endswith((".nwk", ".txt", ".json")) and record["bytes"] <= 100_000:
        path = artifact_path(config=config, job_id=job["job_id"], name=name)
        st.code(body=path.read_text(encoding="utf-8"))
    if record["bytes"] <= 100 * 1024 * 1024:
        if st.button(label="Prepare verified download", key=f"wb-prepare-{job['job_id']}-{name}"):
            path = artifact_path(config=config, job_id=job["job_id"], name=name)
            st.session_state["wb-download"] = (job["job_id"], name, path.read_bytes())
        cached = st.session_state.get("wb-download")
        if cached and cached[:2] == (job["job_id"], name):
            st.download_button(
                label="Download complete file", data=cached[2], file_name=Path(name).name
            )
    else:
        remote_path = job_directory(config=config, job_id=job["job_id"]) / "outputs" / name
        st.code(body=f"rsync -av --progress {shlex.quote(f'SSH_ALIAS:{remote_path}')} ./")
        st.caption("Replace SSH_ALIAS with the Mac SSH alias for this host (for example stefan).")
        st.caption(
            "This file exceeds the browser download limit. Transfer it directly from the server."
        )
