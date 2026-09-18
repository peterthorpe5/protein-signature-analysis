"""Checksum-bound, bounded-memory analytical table checkpoints."""

from __future__ import annotations

import csv
import gzip
import io
import logging
import os
import shutil
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .checksums import sha256_file
from .errors import InputValidationError, PublicationError
from .io_utils import read_json, write_json_atomic
from .schemas import schema_for, table_schemas

LOGGER = logging.getLogger(__name__)

_BATCH_ROWS = 50_000
_PARQUET_PART_ROWS = 1_000_000
_LARGE_TSV_ROWS = 1_000_000
_PROGRESS_ROWS = 1_000_000
_MARKER_NAME = "ANALYSIS_CHECKPOINT.json"
_MANIFEST_NAME = "checkpoint_manifest.json"
_STATE_NAME = "analysis_state.json"


@dataclass(frozen=True)
class CheckpointTable:
    """One verified table stored in bounded analytical formats."""

    name: str
    row_count: int
    parquet_path: Path
    parquet_parts: tuple[Path, ...]
    tsv_path: Path
    tsv_format: str


@dataclass(frozen=True)
class AnalysisCheckpoint:
    """A verified analytical checkpoint that is safe to resume."""

    root: Path
    run_identity_sha256: str
    tables: tuple[CheckpointTable, ...]
    state: Mapping[str, Any]
    input_manifest: tuple[Mapping[str, Any], ...]

    def table(self, *, name: str) -> CheckpointTable:
        """Return one checkpoint table by canonical name.

        Args:
            name: Canonical table name.

        Returns:
            Matching checkpoint table.

        Raises:
            InputValidationError: If the name is not present.
        """

        for table in self.tables:
            if table.name == name:
                return table
        raise InputValidationError(f"Checkpoint lacks canonical table {name!r}.")

    @property
    def counts(self) -> dict[str, int]:
        """Return canonical row counts keyed by table name."""

        return {table.name: table.row_count for table in self.tables}


class RecordSequence(Sequence[Mapping[str, Any]]):
    """Expose model objects as records without materialising a second tuple."""

    def __init__(self, *, values: Sequence[Any]) -> None:
        """Store a model sequence that provides ``to_record``.

        Args:
            values: Immutable or stable source model sequence.
        """

        self._values = values

    def __len__(self) -> int:
        """Return the source row count."""

        return len(self._values)

    def __getitem__(self, index: int | slice) -> Mapping[str, Any] | tuple[Mapping[str, Any], ...]:
        """Convert one source item, or a requested slice, to records."""

        value = self._values[index]
        if isinstance(index, slice):
            return tuple(item.to_record() for item in value)
        return value.to_record()

    def __iter__(self) -> Iterator[Mapping[str, Any]]:
        """Yield converted records one at a time."""

        for value in self._values:
            yield value.to_record()


def checkpoint_directory(*, cache_root: Path, run_identity_sha256: str) -> Path:
    """Return the stable location for one analytical checkpoint.

    Args:
        cache_root: Campaign cache root.
        run_identity_sha256: Configuration/profile/package identity.

    Returns:
        Absolute checkpoint directory.
    """

    identity = _validated_digest(value=run_identity_sha256, field="run identity")
    return Path(cache_root).expanduser().resolve() / "analysis_checkpoints" / identity


def create_analysis_checkpoint(
    *,
    checkpoint_dir: Path,
    run_identity_sha256: str,
    tables: Mapping[str, Sequence[Mapping[str, Any]]],
    state: Mapping[str, Any],
    input_paths: Sequence[Path],
) -> AnalysisCheckpoint:
    """Write all completed analytical tables through a bounded atomic checkpoint.

    Args:
        checkpoint_dir: Final checkpoint directory.
        run_identity_sha256: Configuration/profile/package identity.
        tables: Complete canonical result rows.
        state: JSON-compatible non-tabular analysis state.
        input_paths: Input authorities that must remain unchanged.

    Returns:
        Verified checkpoint.

    Raises:
        InputValidationError: If the table inventory is not canonical.
        PublicationError: If checkpoint publication or verification fails.
    """

    destination = Path(checkpoint_dir).expanduser().resolve()
    identity = _validated_digest(value=run_identity_sha256, field="run identity")
    expected_tables = set(table_schemas())
    if set(tables) != expected_tables:
        raise InputValidationError(
            "Analysis checkpoint table inventory differs from canonical schemas: "
            f"undeclared={sorted(set(tables) - expected_tables)}, "
            f"missing={sorted(expected_tables - set(tables))}."
        )
    if destination.exists():
        return verify_analysis_checkpoint(
            checkpoint_dir=destination,
            run_identity_sha256=identity,
            verify_inputs=True,
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        table_dir = staging / "tables"
        table_dir.mkdir()
        table_records: list[dict[str, Any]] = []
        for table_index, table_name in enumerate(sorted(tables), start=1):
            row_count = len(tables[table_name])
            LOGGER.info(
                "Checkpointing canonical table %d/%d table=%s rows=%d",
                table_index,
                len(tables),
                table_name,
                row_count,
            )
            record = _write_checkpoint_table(
                table_name=table_name,
                records=tables[table_name],
                table_dir=table_dir,
                expected_rows=row_count,
            )
            table_records.append(record)
            LOGGER.info(
                "Checkpointed canonical table %d/%d table=%s rows=%d",
                table_index,
                len(tables),
                table_name,
                row_count,
            )
        write_json_atomic(path=staging / _STATE_NAME, value=dict(state))
        input_manifest = _build_input_manifest(input_paths=input_paths)
        generated_assets = _build_generated_asset_manifest(state=state)
        manifest = {
            "schema_version": 1,
            "status": "COMPLETE",
            "run_identity_sha256": identity,
            "inputs": input_manifest,
            "generated_assets": generated_assets,
            "state": {
                "relative_path": _STATE_NAME,
                "size_bytes": (staging / _STATE_NAME).stat().st_size,
                "sha256": sha256_file(path=staging / _STATE_NAME),
            },
            "tables": table_records,
        }
        write_json_atomic(path=staging / _MANIFEST_NAME, value=manifest)
        write_json_atomic(
            path=staging / _MARKER_NAME,
            value={
                "status": "COMPLETE",
                "manifest_sha256": sha256_file(path=staging / _MANIFEST_NAME),
            },
        )
        os.replace(staging, destination)
    except Exception as error:
        shutil.rmtree(staging, ignore_errors=True)
        if isinstance(error, (InputValidationError, PublicationError)):
            raise
        raise PublicationError(
            f"Could not create analytical checkpoint {destination}: {error}"
        ) from error
    LOGGER.info("Published checksum-bound analytical checkpoint %s", destination)
    return verify_analysis_checkpoint(
        checkpoint_dir=destination,
        run_identity_sha256=identity,
        verify_inputs=True,
    )


def verify_analysis_checkpoint(
    *, checkpoint_dir: Path, run_identity_sha256: str, verify_inputs: bool
) -> AnalysisCheckpoint:
    """Validate a checkpoint marker, manifest, inputs and every table file.

    Args:
        checkpoint_dir: Candidate checkpoint directory.
        run_identity_sha256: Required configuration/profile/package identity.
        verify_inputs: Whether original input authorities must still match.

    Returns:
        Verified checkpoint descriptor.

    Raises:
        PublicationError: If any checkpoint invariant fails.
    """

    root = Path(checkpoint_dir).expanduser().resolve()
    identity = _validated_digest(value=run_identity_sha256, field="run identity")
    marker_path = root / _MARKER_NAME
    manifest_path = root / _MANIFEST_NAME
    if not marker_path.is_file() or not manifest_path.is_file():
        raise PublicationError(f"Analysis checkpoint is incomplete: {root}")
    marker = read_json(path=marker_path)
    manifest = read_json(path=manifest_path)
    if not isinstance(marker, dict) or marker.get("status") != "COMPLETE":
        raise PublicationError(f"Invalid analysis checkpoint marker: {marker_path}")
    if marker.get("manifest_sha256") != sha256_file(path=manifest_path):
        raise PublicationError(f"Analysis checkpoint manifest checksum mismatch: {manifest_path}")
    if not isinstance(manifest, dict) or manifest.get("status") != "COMPLETE":
        raise PublicationError(f"Invalid analysis checkpoint manifest: {manifest_path}")
    if manifest.get("run_identity_sha256") != identity:
        raise PublicationError("Analysis checkpoint was created for a different run identity.")
    input_manifest = manifest.get("inputs")
    if not isinstance(input_manifest, list):
        raise PublicationError("Analysis checkpoint input manifest is malformed.")
    if verify_inputs:
        _verify_input_manifest(records=input_manifest)
    generated_assets = manifest.get("generated_assets")
    if not isinstance(generated_assets, list):
        raise PublicationError("Analysis checkpoint generated-asset manifest is malformed.")
    _verify_generated_asset_manifest(records=generated_assets)
    state_record = manifest.get("state")
    if not isinstance(state_record, dict):
        raise PublicationError("Analysis checkpoint state record is malformed.")
    state_path = root / str(state_record.get("relative_path", ""))
    _verify_file_record(root=root, path=state_path, record=state_record)
    state = read_json(path=state_path)
    if not isinstance(state, dict):
        raise PublicationError("Analysis checkpoint state must be a JSON object.")
    raw_tables = manifest.get("tables")
    if not isinstance(raw_tables, list):
        raise PublicationError("Analysis checkpoint table manifest is malformed.")
    tables: list[CheckpointTable] = []
    names: set[str] = set()
    declared_paths = {_MARKER_NAME, _MANIFEST_NAME, _STATE_NAME}
    for record in raw_tables:
        if not isinstance(record, dict) or not isinstance(record.get("name"), str):
            raise PublicationError("Analysis checkpoint contains a malformed table record.")
        name = record["name"]
        if name in names:
            raise PublicationError(f"Duplicate checkpoint table record: {name!r}")
        names.add(name)
        schema_for(table_name=name)
        parquet_record = record.get("parquet")
        tsv_record = record.get("tsv")
        if not isinstance(parquet_record, dict) or not isinstance(tsv_record, dict):
            raise PublicationError(f"Checkpoint table {name!r} lacks bounded formats.")
        parquet_path = root / str(parquet_record.get("relative_path", ""))
        tsv_path = root / str(tsv_record.get("relative_path", ""))
        parquet_parts = _verify_parquet_dataset(
            root=root,
            dataset_path=parquet_path,
            record=parquet_record,
        )
        _verify_file_record(root=root, path=tsv_path, record=tsv_record)
        declared_paths.add(str(tsv_path.relative_to(root)))
        declared_paths.update(str(path.relative_to(root)) for path in parquet_parts)
        try:
            row_count = int(record["row_count"])
        except (KeyError, TypeError, ValueError) as error:
            raise PublicationError(
                f"Checkpoint table {name!r} has an invalid row count."
            ) from error
        part_row_limit = parquet_record.get("part_row_limit")
        if (
            isinstance(part_row_limit, bool)
            or not isinstance(part_row_limit, int)
            or part_row_limit < 1
        ):
            raise PublicationError(f"Checkpoint table {name!r} has an invalid part row limit.")
        parquet_part_counts = tuple(
            pq.ParquetFile(path).metadata.num_rows for path in parquet_parts
        )
        if any(count > part_row_limit for count in parquet_part_counts) or any(
            count != part_row_limit for count in parquet_part_counts[:-1]
        ):
            raise PublicationError(f"Checkpoint table {name!r} Parquet partitioning differs.")
        parquet_row_count = sum(parquet_part_counts)
        if row_count < 0 or parquet_row_count != row_count:
            raise PublicationError(f"Checkpoint table {name!r} Parquet row count differs.")
        tsv_format = str(record.get("tsv_format", ""))
        expected_suffix = ".tsv.gz" if tsv_format == "TSV.GZ" else ".tsv"
        if tsv_format not in {"TSV", "TSV.GZ"} or not tsv_path.name.endswith(expected_suffix):
            raise PublicationError(f"Checkpoint table {name!r} TSV format differs.")
        tables.append(
            CheckpointTable(
                name=name,
                row_count=row_count,
                parquet_path=parquet_path,
                parquet_parts=parquet_parts,
                tsv_path=tsv_path,
                tsv_format=tsv_format,
            )
        )
    if names != set(table_schemas()):
        raise PublicationError("Analysis checkpoint does not contain every canonical table.")
    actual_paths = {str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}
    if actual_paths != declared_paths:
        raise PublicationError(
            "Analysis checkpoint file inventory differs from its manifest: "
            f"undeclared={sorted(actual_paths - declared_paths)[:10]}, "
            f"missing={sorted(declared_paths - actual_paths)[:10]}."
        )
    LOGGER.info(
        "Reused checksum-bound analytical checkpoint tables=%d rows=%d path=%s",
        len(tables),
        sum(table.row_count for table in tables),
        root,
    )
    return AnalysisCheckpoint(
        root=root,
        run_identity_sha256=identity,
        tables=tuple(sorted(tables, key=lambda table: table.name)),
        state=state,
        input_manifest=tuple(dict(record) for record in input_manifest),
    )


def _write_checkpoint_table(
    *,
    table_name: str,
    records: Sequence[Mapping[str, Any]],
    table_dir: Path,
    expected_rows: int,
) -> dict[str, Any]:
    """Stream one canonical table to Parquet and TSV without whole-table copies."""

    schema = schema_for(table_name=table_name)
    compressed = expected_rows > _LARGE_TSV_ROWS
    tsv_suffix = ".tsv.gz" if compressed else ".tsv"
    parquet_path = table_dir / f"{table_name}.parquet"
    tsv_path = table_dir / f"{table_name}{tsv_suffix}"
    parquet_temp = table_dir / f".{table_name}.parquet.tmp"
    tsv_temp = table_dir / f".{table_name}{tsv_suffix}.tmp"
    writer: pq.ParquetWriter | None = None
    part_index = 0
    part_rows = 0
    count = 0
    batch: list[dict[str, Any]] = []
    fields = tuple(schema.names)
    required_fields = frozenset(fields)
    try:
        with tsv_temp.open(mode="wb") as binary_handle:
            compressed_handle = (
                gzip.GzipFile(filename="", fileobj=binary_handle, mode="wb", mtime=0)
                if compressed
                else binary_handle
            )
            text_handle = io.TextIOWrapper(compressed_handle, encoding="utf-8", newline="")
            try:
                tsv_writer = csv.writer(
                    text_handle,
                    delimiter="\t",
                    lineterminator="\n",
                )
                tsv_writer.writerow(fields)
                for raw_record in records:
                    record = dict(raw_record)
                    if record.keys() != required_fields:
                        missing = sorted(required_fields - record.keys())
                        extra = sorted(record.keys() - required_fields)
                        raise PublicationError(
                            f"Canonical table {table_name!r} row shape differs: "
                            f"missing={missing}, extra={extra}."
                        )
                    batch.append(record)
                    count += 1
                    if len(batch) >= _BATCH_ROWS:
                        _write_tsv_rows(writer=tsv_writer, fields=fields, rows=batch)
                        writer, part_index, part_rows = _write_parquet_batch(
                            dataset_dir=parquet_temp,
                            schema=schema,
                            rows=batch,
                            writer=writer,
                            part_index=part_index,
                            part_rows=part_rows,
                        )
                        batch.clear()
                    if count % _PROGRESS_ROWS == 0:
                        LOGGER.info(
                            "Checkpoint table progress table=%s rows=%d/%d",
                            table_name,
                            count,
                            expected_rows,
                        )
                if batch:
                    _write_tsv_rows(writer=tsv_writer, fields=fields, rows=batch)
                    writer, part_index, part_rows = _write_parquet_batch(
                        dataset_dir=parquet_temp,
                        schema=schema,
                        rows=batch,
                        writer=writer,
                        part_index=part_index,
                        part_rows=part_rows,
                    )
                    batch.clear()
                if writer is None and part_index == 0:
                    writer = _open_parquet_part(
                        dataset_dir=parquet_temp,
                        schema=schema,
                        part_index=part_index,
                    )
                text_handle.flush()
            finally:
                if compressed:
                    text_handle.close()
                else:
                    text_handle.detach()
            binary_handle.flush()
            os.fsync(binary_handle.fileno())
        if writer is not None:
            writer.close()
            writer = None
        if count != expected_rows:
            raise PublicationError(
                f"Canonical table {table_name!r} changed while checkpointing: "
                f"expected={expected_rows}, observed={count}."
            )
        os.replace(parquet_temp, parquet_path)
        os.replace(tsv_temp, tsv_path)
    except (OSError, TypeError, ValueError, pa.ArrowException, PublicationError) as error:
        if writer is not None:
            writer.close()
        shutil.rmtree(parquet_temp, ignore_errors=True)
        tsv_temp.unlink(missing_ok=True)
        if isinstance(error, PublicationError):
            raise
        raise PublicationError(
            f"Could not checkpoint canonical table {table_name!r}: {error}"
        ) from error
    return {
        "name": table_name,
        "row_count": count,
        "tsv_format": "TSV.GZ" if compressed else "TSV",
        "parquet": {
            "relative_path": str(parquet_path.relative_to(table_dir.parent)),
            "format": "PARTITIONED_PARQUET",
            "part_row_limit": _PARQUET_PART_ROWS,
            "parts": [
                _file_record(root=table_dir.parent, path=path)
                for path in sorted(parquet_path.glob("part-*.parquet"))
            ],
        },
        "tsv": _file_record(root=table_dir.parent, path=tsv_path),
    }


def _write_parquet_batch(
    *,
    dataset_dir: Path,
    schema: pa.Schema,
    rows: Sequence[Mapping[str, Any]],
    writer: pq.ParquetWriter | None,
    part_index: int,
    part_rows: int,
) -> tuple[pq.ParquetWriter | None, int, int]:
    """Append one bounded batch across deterministic Parquet part files."""

    if _PARQUET_PART_ROWS < 1:
        raise PublicationError("Parquet part row limit must be positive.")
    offset = 0
    active = writer
    current_index = part_index
    current_rows = part_rows
    while offset < len(rows):
        if active is None:
            active = _open_parquet_part(
                dataset_dir=dataset_dir,
                schema=schema,
                part_index=current_index,
            )
        take = min(_PARQUET_PART_ROWS - current_rows, len(rows) - offset)
        chunk = rows[offset : offset + take]
        table = pa.Table.from_pylist(chunk, schema=schema)
        active.write_table(table, row_group_size=min(_BATCH_ROWS, take))
        current_rows += take
        offset += take
        if current_rows == _PARQUET_PART_ROWS:
            active.close()
            active = None
            current_index += 1
            current_rows = 0
    return active, current_index, current_rows


def _write_tsv_rows(
    *,
    writer: Any,
    fields: Sequence[str],
    rows: Sequence[Mapping[str, Any]],
) -> None:
    """Write one already shape-validated TSV batch without dictionary rescans."""

    writer.writerows(tuple(row[field] for field in fields) for row in rows)


def _open_parquet_part(
    *,
    dataset_dir: Path,
    schema: pa.Schema,
    part_index: int,
) -> pq.ParquetWriter:
    """Open one deterministic Zstandard-compressed Parquet dataset part."""

    dataset_dir.mkdir(parents=True, exist_ok=True)
    return pq.ParquetWriter(
        dataset_dir / f"part-{part_index:05d}.parquet",
        schema=schema,
        compression="zstd",
        write_statistics=True,
    )


def _file_record(*, root: Path, path: Path) -> dict[str, Any]:
    """Return a checksum record for one checkpoint file."""

    return {
        "relative_path": str(path.relative_to(root)),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path=path),
    }


def _verify_parquet_dataset(
    *,
    root: Path,
    dataset_path: Path,
    record: Mapping[str, Any],
) -> tuple[Path, ...]:
    """Verify one partitioned Parquet dataset and return its ordered parts."""

    candidate = dataset_path.expanduser().resolve()
    if root not in candidate.parents or not candidate.is_dir() or candidate.is_symlink():
        raise PublicationError(f"Checkpoint Parquet dataset is missing or unsafe: {candidate}")
    if record.get("format") != "PARTITIONED_PARQUET":
        raise PublicationError(f"Checkpoint Parquet dataset format is invalid: {candidate}")
    raw_parts = record.get("parts")
    if not isinstance(raw_parts, list) or not raw_parts:
        raise PublicationError(f"Checkpoint Parquet dataset has no declared parts: {candidate}")
    parts: list[Path] = []
    names: set[str] = set()
    for part_index, part_record in enumerate(raw_parts):
        if not isinstance(part_record, Mapping):
            raise PublicationError(f"Checkpoint Parquet part record is malformed: {candidate}")
        part = root / str(part_record.get("relative_path", ""))
        resolved_part = part.expanduser().resolve()
        expected_name = f"part-{part_index:05d}.parquet"
        if resolved_part.parent != candidate or resolved_part.name != expected_name:
            raise PublicationError(f"Checkpoint Parquet part is outside its dataset: {part}")
        if resolved_part.name in names:
            raise PublicationError(f"Duplicate checkpoint Parquet part: {resolved_part.name}")
        names.add(resolved_part.name)
        _verify_file_record(root=root, path=resolved_part, record=part_record)
        parts.append(resolved_part)
    ordered = tuple(sorted(parts))
    actual = tuple(sorted(candidate.glob("part-*.parquet")))
    if ordered != actual:
        raise PublicationError(f"Checkpoint Parquet part inventory differs: {candidate}")
    return ordered


def _verify_file_record(*, root: Path, path: Path, record: Mapping[str, Any]) -> None:
    """Verify one checkpoint file record without permitting path escape."""

    candidate = path.expanduser().resolve()
    if root not in candidate.parents or not candidate.is_file():
        raise PublicationError(f"Checkpoint file is missing or unsafe: {candidate}")
    if candidate.stat().st_size != record.get("size_bytes"):
        raise PublicationError(f"Checkpoint file size differs: {candidate}")
    if sha256_file(path=candidate) != record.get("sha256"):
        raise PublicationError(f"Checkpoint file checksum differs: {candidate}")


def _build_input_manifest(*, input_paths: Sequence[Path]) -> list[dict[str, Any]]:
    """Build stable input-authority records for checkpoint invalidation."""

    records: list[dict[str, Any]] = []
    sources = sorted({Path(path).expanduser().resolve() for path in input_paths}, key=str)
    for source in sources:
        if not source.is_file():
            raise PublicationError(f"Checkpoint input is missing or not a file: {source}")
        before = source.stat()
        digest = sha256_file(path=source)
        after = source.stat()
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise PublicationError(f"Checkpoint input changed while checksumming: {source}")
        records.append({"path": str(source), "size_bytes": after.st_size, "sha256": digest})
    return records


def _verify_input_manifest(*, records: Sequence[Mapping[str, Any]]) -> None:
    """Verify the original inputs recorded by an analytical checkpoint."""

    for record in records:
        if not isinstance(record, Mapping) or not isinstance(record.get("path"), str):
            raise PublicationError("Checkpoint contains an invalid input authority record.")
        source = Path(str(record["path"])).expanduser().resolve()
        if not source.is_file():
            raise PublicationError(f"Checkpoint input is no longer available: {source}")
        if source.stat().st_size != record.get("size_bytes"):
            raise PublicationError(f"Checkpoint input size changed: {source}")
        if sha256_file(path=source) != record.get("sha256"):
            raise PublicationError(f"Checkpoint input checksum changed: {source}")


def _build_generated_asset_manifest(*, state: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Checksum report and scientific assets required after analysis."""

    records: list[dict[str, Any]] = []
    for source in _state_asset_paths(state=state):
        if not source.is_file() or source.stat().st_size == 0:
            raise PublicationError(f"Checkpoint generated asset is missing or empty: {source}")
        records.append(
            {
                "path": str(source),
                "size_bytes": source.stat().st_size,
                "sha256": sha256_file(path=source),
            }
        )
    return records


def _verify_generated_asset_manifest(*, records: Sequence[Mapping[str, Any]]) -> None:
    """Verify generated assets needed to resume reporting and publication."""

    for record in records:
        if not isinstance(record, Mapping) or not isinstance(record.get("path"), str):
            raise PublicationError("Checkpoint contains an invalid generated-asset record.")
        source = Path(str(record["path"])).expanduser().resolve()
        if not source.is_file() or source.stat().st_size == 0:
            raise PublicationError(f"Checkpoint generated asset is unavailable: {source}")
        if source.stat().st_size != record.get("size_bytes"):
            raise PublicationError(f"Checkpoint generated asset size changed: {source}")
        if sha256_file(path=source) != record.get("sha256"):
            raise PublicationError(f"Checkpoint generated asset checksum changed: {source}")


def _state_asset_paths(*, state: Mapping[str, Any]) -> tuple[Path, ...]:
    """Extract unique absolute asset sources from checkpoint state."""

    raw_sources: list[object] = []
    base = state.get("base_asset_sources", {})
    if isinstance(base, Mapping):
        raw_sources.extend(base.values())
    plots = state.get("existing_plot_assets", [])
    if isinstance(plots, Sequence) and not isinstance(plots, (str, bytes)):
        for item in plots:
            if isinstance(item, Sequence) and not isinstance(item, (str, bytes)) and len(item) == 2:
                raw_sources.append(item[1])
    return tuple(
        sorted(
            {Path(str(source)).expanduser().resolve() for source in raw_sources},
            key=str,
        )
    )


def _validated_digest(*, value: str, field: str) -> str:
    """Validate one lower-case SHA-256 identifier."""

    text = str(value).strip()
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise InputValidationError(f"{field.capitalize()} must be a lower-case SHA-256 digest.")
    return text
