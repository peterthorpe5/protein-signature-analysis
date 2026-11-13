"""Shared workbench controls for job handoffs, execution and private scratch."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import streamlit as st

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.config import WorkbenchConfig
from protein_signatures.workbench.jobs import MAX_UPLOAD_BYTES, artifact_path, read_job, submit_job
from protein_signatures.workbench.local_worker import start_local_job
from protein_signatures.workbench.storage import clear_scratch, preview_cleanup


def bounded_bytes(*, path: Path) -> bytes:
    """Read at most one viewer-sized artifact, rejecting oversized input."""
    with path.open(mode="rb") as handle:
        content = handle.read(MAX_UPLOAD_BYTES + 1)
    if len(content) > MAX_UPLOAD_BYTES:
        raise InputValidationError("Select a smaller result for the 20 MiB viewer limit.")
    return content


def job_bytes(*, config: WorkbenchConfig, job: dict, name: str) -> bytes:
    """Return a bounded, manifest-verified completed artifact."""
    return bounded_bytes(path=artifact_path(config=config, job_id=job["job_id"], name=name))


def job_label(*, job: dict, artifact: str = "", with_status: bool = False) -> str:
    """Label completed inputs with their recorded HOG or DeepClust provenance."""
    groups = job.get("source_groups", [])
    prefix = ", ".join(groups[:3]) + (" …" if len(groups) > 3 else "") if groups else ""
    return " · ".join(
        filter(
            None,
            [
                prefix,
                job["operation"],
                artifact,
                job["status"] if with_status else "",
                job["created_at"],
                job["job_id"],
            ],
        )
    )


def open_research_tool(*, tool: str, selection: dict) -> None:
    """Switch to a research tool and carry its explicit input selection."""
    st.session_state.update(
        {
            "wb-stage": "5. Research tools",
            "wb-research-tool": tool,
            **selection,
        }
    )


def execution_choice(*, config_path: Path | None, key: str) -> str:
    """Choose the normal worker queue or start just this request on the app host."""
    if config_path is None:
        return "worker"
    return st.radio(
        "Execution",
        options=("worker", "app"),
        horizontal=True,
        key=key,
        format_func=lambda value: {"worker": "Queue for worker", "app": "Run on app host"}[value],
        help="Both options use the same workspace lock. App-host work survives browser closure.",
    )


def queue_job(
    *,
    config: WorkbenchConfig,
    operation: str,
    parameters: dict | None = None,
    uploads: dict | None = None,
    sources: dict | None = None,
    config_path: Path | None = None,
    execution: str = "worker",
) -> str:
    """Save one request and optionally start that exact request on the app host."""
    identifier = submit_job(
        config=config,
        operation=operation,
        parameters=parameters,
        uploads=uploads,
        sources=sources,
        execution=execution,
    )
    if execution == "app" and config_path is not None:
        st.success(f"Saved {operation}: {identifier}. Starting it on this app host.")
        start_local_job(config=config, config_path=config_path, job_id=identifier)
        st.caption(
            "If the workspace is busy, start this waiting request again from Jobs & results."
        )
    else:
        st.success(f"Queued {operation}: {identifier}. Open Jobs & results to follow progress.")
    return identifier


def show_progress(*, config: WorkbenchConfig, key: str, tool: str | None = None) -> None:
    """Report the most recent explicitly submitted job in a tool."""
    identifier = st.session_state.get(key)
    if not identifier:
        return
    job = read_job(config=config, job_id=identifier)
    st.caption(f"Latest request: {job['status']} · {identifier}")
    if job["status"] == "FAILED":
        st.error(job.get("error", "Request failed. Inspect it in Jobs & results."))
    if job["status"] == "COMPLETE" and tool == "JBrowse MSA":
        st.button(
            "Open the completed alignment",
            on_click=open_research_tool,
            kwargs={
                "tool": tool,
                "selection": {
                    "wb-jb-source": "Completed MAFFT job",
                    "wb-jb-job": job_label(job=job),
                    "wb-jb-open": True,
                },
            },
            key=f"{key}-open",
        )
    st.button("Refresh request status", key=f"{key}-refresh")


def render_storage(*, config: WorkbenchConfig) -> None:
    """Show current workspace paths and require an unchanged scratch preview before deletion."""
    with st.expander("Workspace storage and scratch retention"):
        st.code(
            f"Workspace: {config.workspace}\n"
            f"Saved jobs: {Path(config.workspace, 'jobs')}\n"
            f"Per-job scratch: {Path(config.workspace, 'jobs', '<job-id>', 'scratch')}"
        )
        st.caption(
            "Clear only owned scratch from finished jobs. Saved results, inputs, logs, source "
            "databases and configured search indexes stay available."
        )
        days = st.number_input("Keep scratch for at least this many days", 1, 3650, 30, step=1)
        key = "wb-storage-preview"
        if st.button("Preview disposable scratch"):
            st.session_state[key] = (
                config.workspace,
                days,
                preview_cleanup(config=config, retention_days=days),
            )
        preview = st.session_state.get(key)
        if preview is None or preview[:2] != (config.workspace, days):
            return
        if not preview[2]:
            st.info("No finished job scratch meets this retention threshold.")
            return
        st.dataframe(
            [
                {k: v for k, v in asdict(entry).items() if k != "fingerprint"}
                for entry in preview[2]
            ],
            hide_index=True,
            width="stretch",
        )
        if st.button("Clear the previewed scratch"):
            removed = clear_scratch(config=config, entries=preview[2], retention_days=days)
            st.session_state.pop(key, None)
            st.success(f"Cleared {removed:,} bytes of disposable scratch.")
