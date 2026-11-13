"""Private job scratch and previewed retention of disposable tool files."""

from __future__ import annotations

import fcntl
import hashlib
import logging
import re
import shutil
import stat
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

from protein_signatures.errors import InputValidationError

from .config import WorkbenchConfig, bounded_number

LOGGER = logging.getLogger(__name__)
_SCRATCH: ContextVar[Path | None] = ContextVar("workbench_scratch", default=None)
_MARKER = ".workbench-scratch-v1"


@dataclass(frozen=True, kw_only=True)
class ScratchEntry:
    """Fingerprint and size of one reviewed job scratch directory."""

    job_id: str
    bytes: int
    fingerprint: str


@contextmanager
def job_scratch(*, directory: Path):
    """Route tool temporary files to an owned job directory for this worker call."""
    path = directory / "scratch"
    path.mkdir(mode=0o700)
    (path / _MARKER).write_text(data=directory.name, encoding="utf-8")
    token = _SCRATCH.set(path)
    try:
        yield path
    finally:
        _SCRATCH.reset(token)


def temporary_root(*, config: WorkbenchConfig) -> Path:
    """Return the active job scratch or a legacy spill area for standalone commands."""
    path = _SCRATCH.get() or Path(config.workspace, "temporary")
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    return path


def _has_mount(*, root: Path) -> bool:
    """Detect Linux bind mounts that can share the parent filesystem device."""
    mountinfo = Path("/proc/self/mountinfo")
    if not mountinfo.is_file():
        return root.is_mount()
    for row in mountinfo.read_text(encoding="utf-8", errors="surrogateescape").splitlines():
        location = row.split()[4]
        location = re.sub(r"\\([0-7]{3})", lambda match: chr(int(match[1], 8)), location)
        if Path(location).is_relative_to(root):
            return True
    return False


def _entry(*, config: WorkbenchConfig, job: dict, cutoff: datetime) -> ScratchEntry | None:
    """Inspect a terminal job without following links or crossing filesystem boundaries."""
    from .jobs import job_directory

    if job["status"] not in {"COMPLETE", "FAILED"} or not job.get("finished_at"):
        return None
    if datetime.fromisoformat(job["finished_at"]) > cutoff:
        return None
    root = job_directory(config=config, job_id=job["job_id"]) / "scratch"
    marker = root / _MARKER
    if (
        root.is_symlink()
        or not root.is_dir()
        or _has_mount(root=root)
        or marker.is_symlink()
        or not marker.is_file()
        or marker.stat().st_size != len(job["job_id"])
        or marker.read_text(encoding="utf-8") != job["job_id"]
    ):
        return None
    for value in asdict(config).values():
        if isinstance(value, str) and value and Path(value).is_absolute():
            if Path(value).resolve().is_relative_to(root):
                return None
    device = root.stat().st_dev
    digest, size = hashlib.sha256(), 0
    for path in [root, *sorted(root.rglob(pattern="*"))]:
        info = path.lstat()
        if info.st_dev != device or not (stat.S_ISDIR(info.st_mode) or stat.S_ISREG(info.st_mode)):
            return None
        if datetime.fromtimestamp(info.st_mtime, tz=timezone.utc) > cutoff:
            return None
        digest.update(
            repr(
                (str(path.relative_to(root)), info.st_ino, info.st_size, info.st_mtime_ns)
            ).encode()
        )
        if stat.S_ISREG(info.st_mode) and path != marker:
            size += info.st_size
    if not size:
        return None
    return ScratchEntry(job_id=job["job_id"], bytes=size, fingerprint=digest.hexdigest())


def preview_cleanup(*, config: WorkbenchConfig, retention_days: int = 30) -> list[ScratchEntry]:
    """List old owned scratch eligible for removal; keep results and sources intact."""
    from .jobs import list_jobs

    days = bounded_number(
        value=retention_days, name="Retention days", minimum=1, maximum=3650, integer=True
    )
    cutoff = datetime.now(tz=timezone.utc) - timedelta(days=days)
    return [
        entry
        for job in list_jobs(config=config, maximum=1000)
        if (entry := _entry(config=config, job=job, cutoff=cutoff)) is not None
    ]


def clear_scratch(
    *, config: WorkbenchConfig, entries: list[ScratchEntry], retention_days: int = 30
) -> int:
    """Remove reviewed unchanged scratch under the worker lock and return bytes removed."""
    from .jobs import job_directory

    if not entries:
        return 0
    root = Path(config.workspace)
    with (root / "worker.lock").open(mode="a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise InputValidationError(
                "Wait for the active worker before clearing scratch."
            ) from error
        current = set(preview_cleanup(config=config, retention_days=retention_days))
        if not set(entries) <= current:
            raise InputValidationError("Scratch changed since preview. Review the updated list.")
        removed = 0
        for entry in dict.fromkeys(entries):
            path = job_directory(config=config, job_id=entry.job_id) / "scratch"
            shutil.rmtree(path=path)
            LOGGER.info("Cleared %d scratch bytes for %s", entry.bytes, entry.job_id)
            removed += entry.bytes
        return removed
