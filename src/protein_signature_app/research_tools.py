"""Workflows for offline viewers, portable family input and coordinate comparisons."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pyarrow as pa
import streamlit as st
import streamlit.components.v1 as components

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.catalogue import normalise_sequence
from protein_signatures.workbench.config import WorkbenchConfig
from protein_signatures.workbench.context import CONTEXT_FILE, context_bytes
from protein_signatures.workbench.exchange import result_reader
from protein_signatures.workbench.hogs import hog_members, recorded_domains
from protein_signatures.workbench.jobs import MAX_UPLOAD_BYTES, artifact_path, list_jobs, submit_job
from protein_signatures.workbench.providers import folding_fasta, select_proteins
from protein_signatures.workbench.storage import temporary_root
from protein_signatures.workbench.structures import model_residues
from protein_signatures.workbench.subcluster_review import (
    filter_subcluster,
    validate_partition_link,
)

from .alignment_ui import (
    alignment_handoffs,
    context_for_alignment,
    render_alignment_statistics,
    render_alignment_submit,
    saved_tree,
)
from .hog_ui import choose_hog, released_inputs, send_hog_to_tool, use_hog_as_seeds
from .research_viewers import (
    EXAMPLE_ALIGNMENT,
    EXAMPLE_TREE,
    alignment_payload,
    fasta_records,
    feature_payload,
    viewer_document,
)
from .workspace_ui import bounded_bytes, job_bytes, job_label


def completed_job(*, config: WorkbenchConfig, operations: set[str], key: str) -> dict | None:
    """Choose a completed job of the required kind without opening its files."""
    jobs = {
        job_label(job=row): row
        for row in list_jobs(config=config)
        if row["status"] == "COMPLETE" and row["operation"] in operations
    }
    if not jobs:
        st.info("No matching completed jobs yet.")
        return None
    selected = st.selectbox("Completed job", options=list(jobs), key=key)
    return jobs[selected]


def show_viewer(*, viewer: str, payload: dict, key: str) -> None:
    """Display and offer the exact same offline viewer document for download."""
    document = viewer_document(viewer=viewer, payload=payload)
    height = (
        payload["height"] + 44
        if viewer == "jbrowse"
        else min(620, 160 + 43 * len(payload.get("features", [])))
    )
    components.html(document, height=height, scrolling=True)
    st.download_button(
        "Download offline viewer",
        data=document,
        file_name=f"{viewer}_protein_viewer.html",
        mime="text/html",
        key=key,
        on_click="ignore",
    )


def render_alignment_viewer(
    *,
    config: WorkbenchConfig,
    config_path: Path | None = None,
    database: Path | None = None,
) -> None:
    """Align a released HOG or open a completed or uploaded alignment in JBrowse."""
    st.subheader("JBrowse MSA")
    st.caption("Open an aligned FASTA and, optionally, a Newick tree with the same protein tips.")
    options = ["Synthetic example", "Completed MAFFT job", "Upload aligned FASTA"]
    if database is not None:
        options.append("Align a published HOG")
    mode = st.radio("Alignment source", options, key="wb-jb-source")
    mapping = None
    tree = ""
    if mode == "Align a published HOG":
        uploads, sources = released_inputs(database=database, config=config, key="wb-jb-hog")
        if "family" in uploads:
            st.download_button(
                "Download released HOG FASTA",
                uploads["family"],
                "released_hog.fasta",
                on_click="ignore",
            )
            st.download_button(
                "Download exact member context",
                uploads["context"],
                CONTEXT_FILE,
                on_click="ignore",
            )
            st.button(
                "Use these HOG proteins as DeepClust search seeds",
                on_click=use_hog_as_seeds,
                kwargs={"fasta": uploads["family"], "context": uploads["context"]},
            )
            for tool in ("Audited subclustering", "Predict annotations", "ESMFold & models"):
                st.button(
                    f"Use selected HOG in {tool}",
                    on_click=send_hog_to_tool,
                    kwargs={
                        "tool": tool,
                        "fasta": uploads["family"],
                        "context": uploads["context"],
                    },
                    key=f"wb-jb-hog-{tool}",
                )
        render_alignment_submit(
            config=config,
            config_path=config_path,
            uploads=uploads,
            sources=sources,
            key="wb-jb-hog-mafft",
            allow_id_subset=False,
        )
        return
    if mode == "Synthetic example":
        st.info("Synthetic example: six proteins and 56 alignment columns; demonstration only.")
        alignment, tree = EXAMPLE_ALIGNMENT, EXAMPLE_TREE
    elif mode == "Completed MAFFT job":
        job = completed_job(config=config, operations={"align"}, key="wb-jb-job")
        if job is None:
            return
        alignment = job_bytes(config=config, job=job, name="alignment.fasta").decode()
        records = fasta_records(text=alignment, aligned=True)
        render_alignment_statistics(config=config, job=job, records=records)
        tree, mapping = saved_tree(
            config=config,
            config_path=config_path,
            database=database,
            job=job,
            records=records,
        )
        alignment_handoffs(config=config, job=job)
        if not st.checkbox("Open JBrowse viewer", key="wb-jb-open"):
            return
    else:
        upload = st.file_uploader(
            "Aligned protein FASTA (up to 20 MiB)", type=["fasta", "faa", "fa"], key="wb-jb-msa"
        )
        optional = st.file_uploader(
            "Optional Newick tree", type=["nwk", "tree", "txt"], key="wb-jb-tree"
        )
        if upload is None or not st.checkbox("Open JBrowse viewer", key="wb-jb-upload-open"):
            return
        if upload.size > MAX_UPLOAD_BYTES or (optional is not None and optional.size > 1_000_000):
            raise InputValidationError("Upload at most 20 MiB of FASTA and a 1 MB Newick tree.")
        alignment = upload.getvalue().decode()
        tree = optional.getvalue().decode() if optional is not None else ""
    show_viewer(
        viewer="jbrowse",
        payload=alignment_payload(fasta=alignment, newick=tree, tip_labels=mapping),
        key="wb-jb-download",
    )
    st.caption(
        "Tree branches describe the supplied divergence estimate; no biological root is inferred."
    )


def render_released_features(*, config: WorkbenchConfig, database: Path) -> None:
    """Inspect one released member, including proteins without recorded domains or sequence."""
    scope = choose_hog(database=database, key="wb-feature-hog")
    if scope is None:
        return
    frame = hog_members(database=database, scope=scope, maximum=config.max_seed_sequences)
    records = {row.protein_id: row for row in frame.itertuples(index=False)}
    identifier = st.selectbox(
        "Published HOG member",
        list(records),
        key="wb-feature-hog-protein",
        format_func=lambda name: f"{name} · {records[name].species_label or 'Species unrecorded'}",
    )
    protein = records[identifier]
    rows = recorded_domains(database=database, identifier=identifier)
    st.caption(f"{identifier} · {protein.species_label or 'Species unrecorded'}")
    if rows:
        st.dataframe(rows, hide_index=True, width="stretch")
    else:
        st.info("No domain intervals were recorded for this released member.")
    if not protein.sequence:
        st.info("The released resource has no sequence for this member.")
        return
    sequence = normalise_sequence(value=protein.sequence)
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    if protein.sequence_sha256 and protein.sequence_sha256 != digest:
        raise InputValidationError("Published sequence checksum mismatch.")
    features = [
        {
            "sequence_id": identifier,
            "sequence_sha256": digest,
            "accession": row["domain_id"],
            "description": row["domain_name"],
            "start": row["start"],
            "end": row["end"],
        }
        for row in rows
    ]
    show_viewer(
        viewer="nightingale",
        payload=feature_payload(
            identifier=identifier,
            sequence=sequence,
            rows=features,
            evidence_label="published domain evidence",
        ),
        key="wb-feature-hog-download",
    )
    st.caption("The displayed intervals are those recorded in the published result.")


def render_feature_viewer(*, config: WorkbenchConfig, database: Path | None = None) -> None:
    """Show exact-sequence domain positions from verified jobs or a portable TSV."""
    st.subheader("Nightingale protein features")
    options = ["Completed annotation job", "Upload FASTA and annotation TSV"]
    if database is not None:
        options.append("Published HOG features")
    mode = st.radio("Feature source", options, key="wb-feature-source")
    if mode == "Published HOG features":
        render_released_features(config=config, database=database)
        return
    pfam = False
    job = None
    if mode == "Completed annotation job":
        job = completed_job(
            config=config, operations={"annotate", "annotate_proteins"}, key="wb-feature-job"
        )
        if job is None or not st.checkbox("Open feature viewer", key="wb-feature-open"):
            return
        pfam = job["operation"] == "annotate"
        fasta = job_bytes(
            config=config, job=job, name=("sequences.fasta" if pfam else "annotation_input.fasta")
        )
        tsv = job_bytes(config=config, job=job, name="domains.tsv" if pfam else "annotations.tsv")
    else:
        source = st.file_uploader(
            "Exact unaligned FASTA", type=["fa", "faa", "fasta"], key="wb-feature-fasta"
        )
        annotation = st.file_uploader(
            "annotations.tsv with sequence SHA-256", type=["tsv"], key="wb-feature-tsv"
        )
        if source is None or annotation is None:
            return
        if max(source.size, annotation.size) > MAX_UPLOAD_BYTES:
            raise InputValidationError("Each viewer upload must be at most 20 MiB.")
        fasta, tsv = source.getvalue(), annotation.getvalue()
    proteins = fasta_records(text=fasta.decode())
    identifier = st.selectbox(
        "Protein to inspect", options=list(proteins), key="wb-feature-protein"
    )
    if job and CONTEXT_FILE in job["artifacts"]:
        evidence = context_for_alignment(
            data=job_bytes(config=config, job=job, name=CONTEXT_FILE),
            records={identifier: proteins[identifier]},
        ).get(identifier, {})
        if evidence:
            st.caption(
                f"Released species: {evidence['released_species'] or 'unrecorded'} · "
                f"Reviewed taxon: {evidence['reviewed_taxon_id'] or 'unavailable'}"
            )
    reader = csv.DictReader(io.StringIO(tsv.decode()), delimiter="\t")
    required = (
        {"protein_id", "domain_id", "start", "end"}
        if pfam
        else ({"sequence_id", "sequence_sha256", "start", "end"})
    )
    if not required.issubset(reader.fieldnames or []):
        raise InputValidationError(
            "Annotation TSV lacks its required identity or coordinate columns."
        )
    rows = []
    for row in reader:
        if pfam and row["protein_id"] == identifier and row["domain_id"]:
            rows.append(
                {
                    "sequence_id": identifier,
                    "sequence_sha256": hashlib.sha256(proteins[identifier].encode()).hexdigest(),
                    "accession": row["domain_id"],
                    "description": row.get("domain_name", ""),
                    "start": row["start"],
                    "end": row["end"],
                }
            )
        elif not pfam and row["sequence_id"] == identifier:
            rows.append(row)
    show_viewer(
        viewer="nightingale",
        payload=feature_payload(identifier=identifier, sequence=proteins[identifier], rows=rows),
        key="wb-feature-download",
    )
    if rows:
        st.dataframe(rows[:1000], width="stretch", hide_index=True)
    else:
        st.info("No features were recorded for this exact protein in the selected annotation.")
    st.caption(
        "Feature positions are one-based and inclusive. Orthology transfers remain table evidence."
    )


def trace_figure(
    *,
    models: dict[str, list[dict]],
    transform: dict | None = None,
    prediction_confidence: bool = False,
) -> go.Figure:
    """Plot C-alpha traces, optionally applying a verified target-to-query transform.

    Args:
        models: Named coordinate residue records.
        transform: LoL-align transform for the model named target.
        prediction_confidence: The model provider explicitly uses B-factors for confidence.

    Returns:
        Interactive 3D figure; a single model is coloured by its reported B-factor.
    """
    figure = go.Figure()
    colours = ("#126987", "#bc5638")
    for index, (name, records) in enumerate(models.items()):
        xyz = np.array([[record[axis] for axis in ("x", "y", "z")] for record in records])
        if name == "target" and transform:
            xyz = xyz @ np.asarray(transform["rotation"]).T + np.asarray(transform["translation"])
        marker = {"size": 3, "color": colours[index % len(colours)]}
        if len(models) == 1 and prediction_confidence:
            marker.update(
                color=[record["bfactor"] for record in records],
                colorscale="Viridis",
                cmin=0,
                cmax=100,
                showscale=True,
                colorbar={"title": "Reported confidence"},
            )
        figure.add_trace(
            go.Scatter3d(
                x=xyz[:, 0],
                y=xyz[:, 1],
                z=xyz[:, 2],
                mode="lines+markers",
                name=name,
                marker=marker,
                line={"color": colours[index % len(colours)], "width": 3},
                text=[
                    f"{row['residue']}{row['model_residue']} · chain {row['chain']}"
                    for row in records
                ],
                hovertemplate="%{text}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        scene={"aspectmode": "data"}, height=525, margin={"l": 0, "r": 0, "t": 0, "b": 0}
    )
    return figure


def render_predicted_models(*, config: WorkbenchConfig) -> None:
    """Inspect a verified ESMFold model and its confidence without another API call."""
    job = completed_job(config=config, operations={"predict_models"}, key="wb-model-result")
    if job is None or not st.checkbox("Inspect a predicted model", key="wb-model-inspect"):
        return
    rows = list(
        csv.DictReader(
            io.StringIO(job_bytes(config=config, job=job, name="models.tsv").decode()),
            delimiter="\t",
        )
    )
    choices = {row["sequence_id"]: row["model_file"] for row in rows}
    selected = st.selectbox("Predicted protein", list(choices), key="wb-model-protein")
    content = job_bytes(config=config, job=job, name=choices[selected])
    _, _, residues = model_residues(content=content, maximum=config.max_structure_residues)
    st.plotly_chart(
        trace_figure(models={selected: residues}, prediction_confidence=True), width="stretch"
    )
    st.caption("The colour shows provider-reported confidence. The model is a prediction.")
    st.download_button(
        "Download model PDB",
        data=content,
        file_name=Path(choices[selected]).name,
        on_click="ignore",
    )


def render_imported_model(*, config: WorkbenchConfig) -> None:
    """Inspect a portable model only after matching it to the supplied exact sequence."""
    with st.expander("Inspect a model from another machine"):
        fasta = st.file_uploader(
            "Original one-protein FASTA (20 MiB)",
            type=["fasta", "faa", "fa"],
            key="wb-import-model-fasta",
        )
        model = st.file_uploader(
            "Monomer PDB or mmCIF (20 MiB)",
            type=["pdb", "cif", "mmcif"],
            key="wb-import-model-file",
        )
        if fasta is None or model is None:
            return
        if max(fasta.size, model.size) > MAX_UPLOAD_BYTES:
            raise InputValidationError("Each model upload must be at most 20 MiB.")
        proteins = fasta_records(text=fasta.getvalue().decode())
        if len(proteins) != 1:
            raise InputValidationError(
                "Supply the exact one-protein FASTA used to build the model."
            )
        identifier, expected = next(iter(proteins.items()))
        _, observed, residues = model_residues(
            content=model.getvalue(), maximum=config.max_structure_residues
        )
        if observed != expected:
            raise InputValidationError(
                "The coordinate model does not match the exact FASTA sequence."
            )
        st.plotly_chart(trace_figure(models={identifier: residues}), width="stretch")
        st.caption(
            "The supplied model matches the FASTA sequence. Its B-factors are not interpreted."
        )


def portable_model_request(
    *,
    config: WorkbenchConfig,
    selected: list[str],
    uploads: dict,
    sources: dict,
) -> tuple[bytes, bytes]:
    """Prepare the exact inputs and a portable worker request without calling a provider."""
    if "family" in sources:
        reference = sources["family"]
        path = artifact_path(config=config, job_id=reference["job_id"], name=reference["name"])
        if "subcluster_membership" in sources:
            partition = validate_partition_link(sources=sources)
            membership = artifact_path(
                config=config,
                job_id=reference["job_id"],
                name="membership.parquet",
            )
            with TemporaryDirectory(dir=temporary_root(config=config)) as directory:
                subset = Path(directory) / "selected.fasta"
                filter_subcluster(
                    config=config,
                    fasta=path,
                    membership=membership,
                    subcluster_id=partition,
                    destination=subset,
                )
                proteins = select_proteins(
                    source=subset,
                    selected_ids=selected,
                    maximum=20,
                    scan_limit=config.max_sequences,
                )
        else:
            proteins = select_proteins(
                source=path,
                selected_ids=selected,
                maximum=20,
                scan_limit=config.max_sequences,
            )
    else:
        proteins = fasta_records(text=uploads["family"].decode())
        if selected:
            if set(selected) - proteins.keys():
                raise InputValidationError("An exact selected protein is absent from the FASTA.")
            proteins = {name: proteins[name] for name in dict.fromkeys(selected)}
    fasta = folding_fasta(proteins=proteins, maximum_length=config.esm_max_length)
    context = uploads.get("context")
    if "context" in sources:
        reference = sources["context"]
        context = bounded_bytes(
            path=artifact_path(
                config=config,
                job_id=reference["job_id"],
                name=reference["name"],
            )
        )
    evidence = context_for_alignment(data=context, records=proteins) if context else {}
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, mode="w", compression=zipfile.ZIP_DEFLATED) as bundle:
        bundle.writestr("model_inputs.fasta", fasta)
        bundle.writestr("parameters.json", '{"remote_allowed": false}\n')
        if evidence:
            bundle.writestr(CONTEXT_FILE, context_bytes(rows=evidence.values()))
        bundle.writestr(
            "README.txt",
            (
                "Run on the worker machine after configuring its ESMFold endpoint.\n"
                "Review model_inputs.fasta. To approve sending the sequences to the configured "
                "provider, explicitly set remote_allowed to true in parameters.json.\n\n"
                'protein-signature-workbench --config "$PS_PROFILE" --action submit \\\n'
                "  --operation predict_models --family-fasta model_inputs.fasta \\\n"
                "  --parameters parameters.json\n\n"
                + (
                    "Add --context-tsv family_context.tsv to retain member context.\n"
                    if evidence
                    else ""
                )
                + "A running worker executes the job. Follow it in Jobs & results.\n"
            ),
        )
    return fasta, archive.getvalue()


def render_comparison_report() -> None:
    """Filter a portable rerun by the prior method and the observed LoL status."""
    with st.expander("Inspect a structural rerun report"):
        upload = st.file_uploader(
            "structural_method_comparison.tsv (20 MiB)", type=["tsv"], key="wb-comparison-report"
        )
        if upload is None:
            return
        if upload.size > MAX_UPLOAD_BYTES:
            raise InputValidationError("Use a comparison report of at most 20 MiB.")
        frame = pd.read_csv(
            io.BytesIO(upload.getvalue()), sep="\t", dtype=str, keep_default_na=False
        )
        required = {"previous_alignment_tool", "previous_cluster_id", "lol_status"}
        if not required.issubset(frame.columns):
            raise InputValidationError("Select a structural_method_comparison.tsv from a rerun.")
        methods = sorted(frame["previous_alignment_tool"].unique())
        method = st.selectbox("Previous method", ["All", *methods], key="wb-report-method")
        status = st.selectbox(
            "LoL-align status",
            ["All", *sorted(frame["lol_status"].unique())],
            key="wb-report-status",
        )
        if method != "All":
            frame = frame.loc[frame["previous_alignment_tool"] == method]
        if status != "All":
            frame = frame.loc[frame["lol_status"] == status]
        st.caption(
            f"{len(frame):,} matching rows. Up to 1,000 shown; the original TSV is unchanged."
        )
        st.dataframe(frame.head(1000), width="stretch", hide_index=True)


def render_structure_jobs(*, config: WorkbenchConfig) -> None:
    """Queue a pair of explicit monomer models and inspect the local comparison."""
    st.subheader("LoL-align structural comparison")
    st.caption("Compare two single-chain PDB or mmCIF models with local 3D correspondence.")
    models = {
        f"{job['job_id']} · {name}": {"job_id": job["job_id"], "name": name}
        for job in list_jobs(config=config)
        if job["status"] == "COMPLETE" and job["operation"] == "predict_models"
        for name in job["artifacts"]
        if name.endswith(".pdb")
    }
    uploads, sources = {}, {}
    for role in ("query", "target"):
        key = f"{role}_model"
        label = "Reference query" if role == "query" else "Mobile target"
        selected = st.selectbox(label, ["Upload a model", *models], key=f"wb-{role}-model-choice")
        if selected != "Upload a model":
            sources[key] = models[selected]
        else:
            upload = st.file_uploader(
                f"{label} coordinate file",
                type=["pdb", "cif", "mmcif"],
                key=f"wb-{role}-model-upload",
            )
            if upload is not None:
                uploads[key] = upload.getvalue()
    single = st.checkbox("Use single-domain LoL scoring", key="wb-lol-single")
    if st.button("Queue LoL-align pair", disabled=len(uploads) + len(sources) != 2):
        identifier = submit_job(
            config=config,
            operation="structure_pair",
            uploads=uploads,
            sources=sources,
            parameters={"multidomain": not single},
        )
        st.success(f"Queued {identifier}. Follow it in Jobs & results.")
    render_comparison_report()
    job = completed_job(config=config, operations={"structure_pair"}, key="wb-lol-result")
    if job is None or not st.checkbox("Inspect completed structural comparison", key="wb-lol-view"):
        return
    records = list(
        csv.DictReader(
            io.StringIO(job_bytes(config=config, job=job, name="lolalign.tsv").decode()),
            delimiter="\t",
        )
    )
    st.dataframe(records, hide_index=True, width="stretch")
    if "target_to_query.json" in job["artifacts"]:
        transform = json.loads(job_bytes(config=config, job=job, name="target_to_query.json"))
        traces = {}
        for role in ("query", "target"):
            name = next(
                candidate
                for candidate in (f"{role}.pdb", f"{role}.cif")
                if candidate in job["artifacts"]
            )
            _, _, traces[role] = model_residues(
                content=job_bytes(config=config, job=job, name=name),
                maximum=config.max_structure_residues,
            )
        st.plotly_chart(trace_figure(models=traces, transform=transform), width="stretch")
        st.caption(
            "Mobile C-alpha coordinates placed onto the reference using the local correspondence."
        )
    else:
        st.info("Foldseek returned no local correspondence for these exact models.")


def render_deepclust_import(*, config: WorkbenchConfig) -> None:
    """Import a portable result subset without touching the mounted global database."""
    st.subheader("Import a DeepClust results export")
    st.caption(
        "Import a selected Parquet result and reuse its complete FASTA in other workbench jobs."
    )
    upload = st.file_uploader(
        "DeepClust result Parquet (up to 20 MiB)", type=["parquet"], key="wb-deepclust-export"
    )
    if upload is None:
        return
    content = upload.getvalue()
    with result_reader(source=content) as reader:
        names = reader.schema_arrow.names
        st.dataframe(
            pa.Table.from_batches([next(reader.iter_batches(batch_size=20))]),
            hide_index=True,
            width="stretch",
        )
    identifier = st.selectbox("Protein identifier field", names, key="wb-import-id")
    sequence = st.selectbox(
        "Protein sequence field", names, index=min(1, len(names) - 1), key="wb-import-sequence"
    )
    cluster = st.selectbox("Optional exact cluster field", ["", *names], key="wb-import-cluster")
    value = st.text_input("Exact cluster ID", disabled=not cluster, key="wb-import-cluster-id")
    if st.button(
        "Queue portable family import",
        disabled=identifier == sequence or bool(cluster) != bool(value),
    ):
        job = submit_job(
            config=config,
            operation="import_deepclust",
            uploads={"deepclust_export": content},
            parameters={
                "id_column": identifier,
                "sequence_column": sequence,
                "cluster_column": cluster,
                "cluster_id": value,
            },
        )
        st.success(f"Queued {job}. Use its family.fasta when the job completes.")
