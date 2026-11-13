"""Durable filesystem queue with exclusive workers and checksummed job inputs."""

from __future__ import annotations

import fcntl
import json
import logging
import os
import re
import socket
import traceback
import uuid
from collections.abc import Mapping
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from protein_signatures import __version__
from protein_signatures.checksums import sha256_file
from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic

from .config import WorkbenchConfig, config_identity
from .storage import job_scratch

LOGGER = logging.getLogger(__name__)
OPERATIONS = (
    "resolve",
    "retrieve",
    "search",
    "screen",
    "align",
    "tree",
    "subcluster",
    "annotate",
    "score",
    "prepare",
    "build_search",
    "structure_pair",
    "predict_models",
    "annotate_proteins",
    "import_deepclust",
)
JOB_PATTERN = re.compile(pattern=r"job_[0-9a-f]{32}")
MAX_UPLOAD_BYTES = 20 * 1024 * 1024
INPUT_SUFFIXES = {
    "family": ".fasta",
    "alignment": ".fasta",
    "context": ".tsv",
    "ebi_email": ".private",
    "seeds": ".fasta",
    "query_model": ".model",
    "target_model": ".model",
    "deepclust_export": ".parquet",
}


def timestamp() -> str:
    """Return the current UTC timestamp for an audit record."""
    return datetime.now(tz=timezone.utc).isoformat()


def job_directory(*, config: WorkbenchConfig, job_id: str) -> Path:
    """Resolve an exact job directory, rejecting traversal and symlink escapes.

    Args:
        config: Trusted workspace.
        job_id: Previously returned job identifier.

    Returns:
        Contained job directory.

    Raises:
        InputValidationError: If the ID or directory is unsafe.
    """
    if not isinstance(job_id, str) or JOB_PATTERN.fullmatch(string=job_id) is None:
        raise InputValidationError("Invalid workbench job identifier.")
    root = (Path(config.workspace) / "jobs").resolve()
    candidate = root / job_id
    if candidate.is_symlink() or not candidate.resolve().is_relative_to(root):
        raise InputValidationError("Job directory escapes its workspace.")
    return candidate


def read_job(*, config: WorkbenchConfig, job_id: str) -> dict:
    """Read a durable state record for one known job.

    Args:
        config: Trusted workspace.
        job_id: Exact job ID.

    Returns:
        Current job record.

    Raises:
        InputValidationError: If the job is unavailable or corrupt.
    """
    try:
        document = json.loads(
            s=(job_directory(config=config, job_id=job_id) / "state.json").read_text(
                encoding="utf-8"
            )
        )
        if document["job_id"] != job_id:
            raise ValueError("Job identity changed")
        return document
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise InputValidationError(f"Job state is unavailable: {job_id}") from error


def list_jobs(*, config: WorkbenchConfig, maximum: int = 200) -> list[dict]:
    """List recent jobs; incomplete submissions are not exposed.

    Args:
        config: Trusted workspace.
        maximum: Maximum displayed jobs, between 1 and 1000.

    Returns:
        Most recent job records first.
    """
    if isinstance(maximum, bool) or not isinstance(maximum, int) or not 1 <= maximum <= 1000:
        raise InputValidationError("Job display limit must be between 1 and 1000.")
    root = Path(config.workspace) / "jobs"
    candidates = sorted(
        root.glob(pattern="job_*/state.json"), key=lambda p: p.stat().st_mtime_ns, reverse=True
    )
    return [read_job(config=config, job_id=p.parent.name) for p in candidates[:maximum]]


def artifact_path(*, config: WorkbenchConfig, job_id: str, name: str) -> Path:
    """Verify a completed output against its manifest before reuse or download.

    Args:
        config: Trusted workspace.
        job_id: Completed source job.
        name: Exact relative manifest entry.

    Returns:
        Contained, checksum-verified local path.

    Raises:
        InputValidationError: If the job, path or checksum is invalid.
    """
    state = read_job(config=config, job_id=job_id)
    root = job_directory(config=config, job_id=job_id) / "outputs"
    path = (root / name).resolve()
    if (
        state["status"] != "COMPLETE"
        or not path.is_relative_to(root.resolve())
        or path == root.resolve()
    ):
        raise InputValidationError("Select a completed job artifact within its output directory.")
    record = state.get("artifacts", {}).get(name)
    if not record or not path.is_file() or path.stat().st_size != record["bytes"]:
        raise InputValidationError("The selected job artifact is unavailable or changed.")
    if sha256_file(path=path) != record["sha256"]:
        raise InputValidationError("The selected job artifact failed its checksum.")
    return path


def submit_job(
    *,
    config: WorkbenchConfig,
    operation: str,
    parameters: Mapping[str, object] | None = None,
    uploads: Mapping[str, bytes] | None = None,
    sources: Mapping[str, Mapping[str, str]] | None = None,
    execution: str = "worker",
) -> str:
    """Submit named operations with private uploaded snapshots or verified source artifacts.

    Args:
        config: Trusted workspace and source profile.
        operation: Supported task name.
        parameters: JSON-serialisable scientific controls, validated by the worker.
        uploads: Named, bounded FASTA, structure or DeepClust result snapshots.
        sources: Named references containing source job_id and artifact name.
        execution: Reserve the request for an exact app-host start or the normal worker.

    Returns:
        Persistent job identifier.

    Raises:
        InputValidationError: If task, inputs or upload bounds are invalid.
    """
    if operation not in OPERATIONS:
        raise InputValidationError("Unsupported workbench operation.")
    if execution not in {"worker", "app"}:
        raise InputValidationError("Choose worker or app-host execution.")
    uploaded, linked = dict(uploads or {}), dict(sources or {})
    if "ebi_email" in linked or ("ebi_email" in uploaded and operation != "annotate_proteins"):
        raise InputValidationError("Supply EBI contact only as a private annotation input.")
    if set(uploaded) & set(linked) or (set(uploaded) | set(linked)) - (
        INPUT_SUFFIXES.keys() | {"features", "subcluster_membership"}
    ):
        raise InputValidationError("Supply each supported named input only once.")
    for name, reference in linked.items():
        expected = (
            {"job_id", "name", "subcluster_id"}
            if name == "subcluster_membership"
            else {
                "job_id",
                "name",
            }
        )
        if (
            not isinstance(reference, Mapping)
            or set(reference) != expected
            or not all(isinstance(value, str) for value in reference.values())
        ):
            raise InputValidationError("A source reference must name an exact completed artifact.")
    parent = linked.get("alignment" if operation == "tree" else "family")
    if parent and "context" not in (uploaded | linked):
        state = read_job(config=config, job_id=parent["job_id"])
        if "family_context.tsv" in state.get("artifacts", {}):
            linked["context"] = {"job_id": parent["job_id"], "name": "family_context.tsv"}
    if "subcluster_membership" in linked:
        from .subcluster_review import validate_partition_link

        validate_partition_link(sources=linked)
        state = read_job(config=config, job_id=linked["subcluster_membership"]["job_id"])
        if state["operation"] != "subcluster":
            raise InputValidationError("Use membership from a completed subcluster job.")
    for name, payload in uploaded.items():
        if (
            name not in INPUT_SUFFIXES
            or not isinstance(payload, bytes)
            or not 0 < len(payload) <= MAX_UPLOAD_BYTES
        ):
            raise InputValidationError(
                "Each upload must contain 1 byte to 20 MiB of supported input data."
            )
    if "ebi_email" in uploaded and len(uploaded["ebi_email"]) > 254:
        raise InputValidationError("EBI contact exceeds the maximum address length.")
    # Verify source integrity before making a queue record visible.
    references = {}
    for name, reference in linked.items():
        path = artifact_path(config=config, job_id=reference["job_id"], name=reference["name"])
        references[name] = {**reference, "sha256": sha256_file(path=path)}
    identifier = "job_" + uuid.uuid4().hex
    root = job_directory(config=config, job_id=identifier)
    root.mkdir(parents=True, mode=0o700)
    inputs = {}
    for name, payload in uploaded.items():
        path = root / f"{name}{INPUT_SUFFIXES[name]}"
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, mode="wb") as handle:
            handle.write(payload)
        inputs[name] = {"file": path.name, "sha256": sha256_file(path=path)}
    request = {
        "job_id": identifier,
        "operation": operation,
        "execution": execution,
        "parameters": dict(parameters or {}),
        "uploads": inputs,
        "sources": references,
        "config_identity": config_identity(config=config),
        "config": {
            name: value for name, value in asdict(config).items() if name != "interpro_email"
        },
        "package_version": __version__,
        "created_at": timestamp(),
    }
    write_json_atomic(path=root / "request.json", value=request)
    write_json_atomic(
        path=root / "state.json",
        value={
            "job_id": identifier,
            "operation": operation,
            "execution": execution,
            "status": "QUEUED",
            "created_at": request["created_at"],
            "request_sha256": sha256_file(path=root / "request.json"),
            "input_jobs": {name: ref["job_id"] for name, ref in references.items()},
        },
    )
    LOGGER.info("Queued %s (%s)", identifier, operation)
    return identifier


def resolve_input(*, config: WorkbenchConfig, request: dict, name: str) -> Path:
    """Verify a saved named input immediately before execution.

    Args:
        config: Trusted workspace.
        request: Immutable job request.
        name: Required named input.

    Returns:
        Validated input path.

    Raises:
        InputValidationError: If the input was modified or is missing.
    """
    if name in request["uploads"]:
        record = request["uploads"][name]
        path = job_directory(config=config, job_id=request["job_id"]) / (
            f"{name}{INPUT_SUFFIXES[name]}"
        )
    elif name in request["sources"]:
        record = request["sources"][name]
        path = artifact_path(config=config, job_id=record["job_id"], name=record["name"])
    else:
        raise InputValidationError(f"This task requires the {name} input.")
    if sha256_file(path=path) != record["sha256"]:
        raise InputValidationError("Job input changed after submission.")
    return path


def cancel_job(*, config: WorkbenchConfig, job_id: str) -> None:
    """Cancel only a waiting job under the same exclusive lock used by workers.

    Args:
        config: Workspace profile.
        job_id: Waiting job.

    Raises:
        InputValidationError: If a worker is active or the job has already started.
    """
    root = Path(config.workspace)
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (root / "worker.lock").open(mode="a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise InputValidationError(
                "A worker is active; retry cancellation after its current task."
            ) from error
        state = read_job(config=config, job_id=job_id)
        if state["status"] != "QUEUED":
            raise InputValidationError("Only queued jobs can be cancelled.")
        state.update(status="CANCELLED", finished_at=timestamp())
        write_json_atomic(
            path=job_directory(config=config, job_id=job_id) / "state.json", value=state
        )


def run_next(*, config: WorkbenchConfig, job_id: str | None = None) -> str | None:
    """Execute one waiting job with one exclusive worker per workspace.

    Args:
        config: Worker configuration, matching the submitted snapshot.
        job_id: Exact waiting job to start; otherwise take the oldest worker request.

    Returns:
        Claimed job ID or None when idle or another worker holds the lock.
    """
    from .operations import execute_job

    if job_id is not None:
        job_directory(config=config, job_id=job_id)

    root = Path(config.workspace)
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (root / "worker.lock").open(mode="a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return None
        queued = []
        # Execution checks the complete queue, even when the UI shows only recent jobs.
        for path in (root / "jobs").glob(pattern="job_*/state.json"):
            state = read_job(config=config, job_id=path.parent.name)
            if state["status"] == "RUNNING":
                LOGGER.warning(
                    "Inspect interrupted/running job %s before starting more work", state["job_id"]
                )
                return None
            if state["status"] == "QUEUED" and (
                state["job_id"] == job_id
                if job_id is not None
                else state.get("execution", "worker") == "worker"
            ):
                queued.append(state)
        if not queued:
            return None
        state = min(queued, key=lambda row: (row["created_at"], row["job_id"]))
        identifier = state["job_id"]
        directory = job_directory(config=config, job_id=identifier)
        state.update(
            status="RUNNING",
            started_at=timestamp(),
            worker_pid=os.getpid(),
            worker_host=socket.gethostname(),
        )
        write_json_atomic(path=directory / "state.json", value=state)
        try:
            request = json.loads(s=(directory / "request.json").read_text(encoding="utf-8"))
            if sha256_file(path=directory / "request.json") != state["request_sha256"]:
                raise InputValidationError("Job request changed after submission.")
            if request["package_version"] != __version__:
                raise InputValidationError(
                    "Package version changed after submission; submit a new job."
                )
            if (
                request["job_id"] != identifier
                or request["operation"] != state["operation"]
                or request.get("execution", "worker") != state.get("execution", "worker")
            ):
                raise InputValidationError("Job request and queue state disagree.")
            if config_identity(config=config) != request["config_identity"]:
                raise InputValidationError(
                    "Profile or source metadata changed after submission; submit a new job."
                )
            output = directory / "outputs"
            output.mkdir(mode=0o700)
            with job_scratch(directory=directory):
                execute_job(config=config, request=request, output=output)
            if config_identity(config=config) != request["config_identity"]:
                raise InputValidationError("Source metadata changed while the job was running.")
            artifacts = {}
            for path in sorted(output.rglob(pattern="*")):
                if path.is_file():
                    artifacts[str(path.relative_to(output))] = {
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path=path),
                    }
            write_json_atomic(
                path=directory / "manifest.json",
                value={
                    "request_sha256": sha256_file(path=directory / "request.json"),
                    "artifacts": artifacts,
                    "source_identity_policy": (
                        "Profile plus source size/mtime; global source bytes not hashed"
                    ),
                    "package_version": __version__,
                    "completed_at": timestamp(),
                },
            )
            state.update(status="COMPLETE", artifacts=artifacts)
            if "family_context.json" in artifacts:
                state["source_groups"] = json.loads(
                    s=(output / "family_context.json").read_text(encoding="utf-8")
                )["source_groups"]
        except Exception as error:
            LOGGER.exception("Workbench job %s failed", identifier)
            (directory / "failure.log").write_text(data=traceback.format_exc(), encoding="utf-8")
            state.update(status="FAILED", error=str(error))
        state["finished_at"] = timestamp()
        write_json_atomic(path=directory / "state.json", value=state)
        return identifier


def recover_interrupted_job(*, config: WorkbenchConfig, job_id: str) -> None:
    """Mark a demonstrably stopped job failed, without restarting or deleting its work.

    Args:
        config: Workspace on the recorded worker host.
        job_id: Interrupted RUNNING job selected explicitly by the operator.

    Raises:
        InputValidationError: If a worker/tool is live, identity is uncertain or state differs.
    """
    root = Path(config.workspace)
    with (root / "worker.lock").open(mode="a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise InputValidationError(
                "The workspace worker is active; recovery is unavailable."
            ) from error
        state = read_job(config=config, job_id=job_id)
        if state["status"] != "RUNNING" or state.get("worker_host") != socket.gethostname():
            raise InputValidationError("Recovery requires an interrupted RUNNING job on this host.")
        try:
            os.kill(state["worker_pid"], 0)
        except ProcessLookupError:
            pass
        except PermissionError as error:
            raise InputValidationError(
                "Cannot establish that the recorded worker has stopped."
            ) from error
        else:
            raise InputValidationError("The recorded worker process is still live.")
        marker = job_directory(config=config, job_id=job_id) / "outputs" / "active_process.json"
        if marker.is_file():
            process = json.loads(s=marker.read_text(encoding="utf-8"))
            try:
                os.killpg(process["process_group"], 0)
            except ProcessLookupError:
                pass
            except PermissionError as error:
                raise InputValidationError(
                    "Cannot establish that the recorded tool group has stopped."
                ) from error
            else:
                raise InputValidationError(
                    "The recorded external tool group is still live; inspect it first."
                )
        state.update(
            status="FAILED",
            finished_at=timestamp(),
            error=("Interrupted worker; explicitly recovered. Submit a new request to retry."),
        )
        write_json_atomic(
            path=job_directory(config=config, job_id=job_id) / "state.json", value=state
        )
