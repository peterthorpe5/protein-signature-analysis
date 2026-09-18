"""Runtime resource bounds for analytical database operations."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import duckdb

LOGGER = logging.getLogger(__name__)

_MINIMUM_DUCKDB_MEMORY_MB = 512
_DUCKDB_MEMORY_FRACTION = 0.70
_MAX_PLAUSIBLE_MEMORY_BYTES = 1024**5


def configure_duckdb_runtime(
    *, connection: duckdb.DuckDBPyConnection, temporary_directory: Path
) -> int | None:
    """Bound DuckDB memory to the scheduler/cgroup allocation and enable spill.

    Args:
        connection: Active trusted internal DuckDB connection.
        temporary_directory: Existing writable directory for external operations.

    Returns:
        Configured memory limit in MiB, or ``None`` when no reliable allocation
        limit was discoverable.

    Raises:
        OSError: If the temporary directory is unavailable.
        duckdb.Error: If DuckDB rejects a validated setting.
    """

    temporary = Path(temporary_directory).expanduser().resolve()
    if not temporary.is_dir():
        raise OSError(f"DuckDB temporary directory is missing: {temporary}")
    escaped = str(temporary).replace("'", "''")
    connection.execute(f"SET temp_directory = '{escaped}'")
    connection.execute("SET preserve_insertion_order = false")
    allocation_mb = _available_memory_mb()
    if allocation_mb is None:
        LOGGER.info("DuckDB memory limit was not overridden; no scheduler/cgroup limit found")
        return None
    limit_mb = max(
        _MINIMUM_DUCKDB_MEMORY_MB,
        int(allocation_mb * _DUCKDB_MEMORY_FRACTION),
    )
    connection.execute(f"SET memory_limit = '{limit_mb}MB'")
    LOGGER.info(
        "Configured DuckDB memory_limit=%dMB allocation=%dMB spill=%s",
        limit_mb,
        allocation_mb,
        temporary,
    )
    return limit_mb


def _available_memory_mb() -> int | None:
    """Return the smallest reliable scheduler or cgroup limit in MiB."""

    candidates: list[int] = []
    per_node = _positive_integer_environment(name="SLURM_MEM_PER_NODE")
    if per_node is not None:
        candidates.append(per_node)
    per_cpu = _positive_integer_environment(name="SLURM_MEM_PER_CPU")
    cpus = _positive_integer_environment(name="SLURM_CPUS_PER_TASK")
    if per_cpu is not None and cpus is not None:
        candidates.append(per_cpu * cpus)
    for path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            text = path.read_text(encoding="utf-8").strip()
        except OSError:
            continue
        if not text.isdigit():
            continue
        value = int(text)
        if 0 < value < _MAX_PLAUSIBLE_MEMORY_BYTES:
            candidates.append(value // (1024 * 1024))
    positive = [value for value in candidates if value >= _MINIMUM_DUCKDB_MEMORY_MB]
    return min(positive) if positive else None


def _positive_integer_environment(*, name: str) -> int | None:
    """Read one positive integer environment value conservatively."""

    text = os.environ.get(name, "").strip()
    if not text.isdigit():
        return None
    value = int(text)
    return value if value > 0 else None
