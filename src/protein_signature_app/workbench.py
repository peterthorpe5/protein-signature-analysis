"""Guided protein-family discovery, server jobs and bounded result inspection."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from protein_signatures.errors import ProteinSignatureError
from protein_signatures.workbench.config import WorkbenchConfig, diagnostics, load_config
from protein_signatures.workbench.jobs import (
    MAX_UPLOAD_BYTES,
    artifact_path,
    cancel_job,
    job_directory,
    list_jobs,
    recover_interrupted_job,
    submit_job,
)
from protein_signatures.workbench.seeds import export_seeds

from .workbench_viewers import alignment_window, tree_figure

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


def queue_job(
    *,
    config: WorkbenchConfig,
    operation: str,
    parameters: dict | None = None,
    uploads: dict | None = None,
    sources: dict | None = None,
) -> None:
    """Submit one explicit user action and report its durable identifier.

    Args:
        config: Operator profile.
        operation: Task name.
        parameters: Named scientific controls.
        uploads: Private FASTA snapshots.
        sources: Completed job artifact references.
    """
    identifier = submit_job(
        config=config, operation=operation, parameters=parameters, uploads=uploads, sources=sources
    )
    st.success(f"Queued {operation}: {identifier}. Open Jobs & results to follow progress.")


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
            render_analysis(config=config)
        else:
            render_jobs(config=config)
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
        options=("Upload or paste protein FASTA", "Select from the completed result"),
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
            st.success("Seed snapshot selected. Continue to Find families.")
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
            queue_job(
                config=config,
                operation="search",
                parameters={**parameters, "maximum_hits": int(maximum)},
                uploads={"seeds": seeds},
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
        options = {
            f"{job['operation']} · {job['created_at']} · {job['job_id']}": job for job in jobs
        }
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


def family_input(*, config: WorkbenchConfig) -> tuple[dict, dict]:
    """Select an unaligned completed artifact or a new bounded upload.

    Args:
        config: Job workspace.

    Returns:
        Upload and source-reference dictionaries.
    """
    st.caption(
        (
            "Use an unaligned FASTA from a completed job, or upload your own. "
            "Alignment matrices are not accepted as protein input."
        )
    )
    candidates = {}
    for job in list_jobs(config=config):
        if job["status"] == "COMPLETE":
            for name in job.get("artifacts", {}):
                if name in {
                    "family.fasta",
                    "candidates.fasta",
                    "sequences.fasta",
                    "expanded.fasta",
                    "seeds.fasta",
                }:
                    candidates[
                        f"{job['operation']} · {name} · {job['created_at']} · {job['job_id']}"
                    ] = {"job_id": job["job_id"], "name": name}
    selected = st.selectbox(
        label="Family FASTA from a completed job",
        options=["Upload instead", *candidates],
        key="wb-family-source",
    )
    if selected != "Upload instead":
        return {}, {"family": candidates[selected]}
    upload = st.file_uploader(
        label="Family protein FASTA (maximum 20 MiB)",
        type=["fa", "faa", "fasta"],
        key="wb-family-upload",
    )
    return ({"family": upload.getvalue()} if upload is not None else {}), {}


def render_analysis(*, config: WorkbenchConfig) -> None:
    """Render family screening, alignment, annotation, signature and expansion tasks.

    Args:
        config: Operator profile.
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
    uploads, sources = family_input(config=config)
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
        st.write(
            (
                "Create a MAFFT alignment and residue conservation table. Large "
                "families require a smaller explicit selection; the app never samples "
                "silently."
            )
        )
        text = st.text_area(
            label="Optional exact sequence IDs to align, one per line",
            help=(
                "Blank selects the complete input, subject to configured limits. Copy "
                "IDs from members.tsv or candidate_audit.tsv."
            ),
        )
        parameters = {
            "selected_ids": parse_ids(text=text),
            "tree": st.checkbox(
                label="Infer a FastTree WAG+gamma tree",
                value=False,
                help=(
                    "At least three sequences. Branch lengths represent inferred "
                    "divergence, not time. tip_labels.tsv links safe tree labels to the "
                    "exact protein IDs."
                ),
            ),
        }
        st.caption(
            f"Limits: {config.max_alignment_sequences:,} sequences and "
            f"{config.max_alignment_cells:,} alignment cells."
        )
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
        if not uploads and "family" not in sources:
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


def render_jobs(*, config: WorkbenchConfig) -> None:
    """Show durable progress, bounded previews and verified on-demand downloads.

    Args:
        config: Job workspace.
    """
    st.subheader("Jobs & results")
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
    options = {
        f"{job['operation']} · {job['status']} · {job['created_at']} · {job['job_id']}": job
        for job in jobs
    }
    selected = st.selectbox(label="Inspect job", options=list(options), key="wb-inspect-job")
    job = options[selected]
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
        st.code(body=f"rsync -av --progress 'denbi:{remote_path}' ./")
        st.caption(
            "This file exceeds the browser download limit. Transfer it directly from the server."
        )
