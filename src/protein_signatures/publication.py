"""Atomic, checksummed TSV, Parquet and DuckDB result publication."""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .checksums import sha256_file
from .errors import PublicationError
from .io_utils import read_json, write_json_atomic, write_tsv_atomic
from .schemas import schema_for

LOGGER = logging.getLogger(__name__)


def publish_result(
    *,
    output_dir: Path,
    tables: Mapping[str, Sequence[Mapping[str, Any]]],
    metadata: Mapping[str, Any],
    input_paths: Sequence[Path],
    asset_sources: Mapping[str, Path],
    resume: bool,
) -> Path:
    """Publish a complete result directory through an atomic rename.

    Args:
        output_dir: New immutable result destination.
        tables: Records keyed by canonical table name.
        metadata: JSON-compatible campaign metadata.
        input_paths: Input authorities to checksum.
        asset_sources: Portable result-relative asset paths mapped to source files.
        resume: Reuse an existing result only after full checksum validation.

    Returns:
        Completed result directory.

    Raises:
        PublicationError: If output exists, publication fails or resume validation fails.
    """

    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        if resume:
            verify_completed_result(result_dir=destination)
            LOGGER.info("Reusing checksum-verified completed result %s", destination)
            return destination
        raise PublicationError(
            f"Output already exists and will not be overwritten: {destination}. "
            "Use --resume only for a completed checksum-valid result."
        )
    initial_input_manifest = _build_input_manifest(input_paths=input_paths)
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        _copy_assets(staging=staging, asset_sources=asset_sources)
        table_dir = staging / "tables"
        table_dir.mkdir()
        for table_name in sorted(tables):
            _write_table(
                table_name=table_name,
                records=tables[table_name],
                table_dir=table_dir,
            )
        _create_duckdb(path=staging / "protein_signatures.duckdb", table_dir=table_dir)
        write_json_atomic(path=staging / "run_metadata.json", value=dict(metadata))
        input_manifest = _build_input_manifest(input_paths=input_paths)
        if input_manifest != initial_input_manifest:
            raise PublicationError("An input authority changed during result publication.")
        output_manifest = _manifest_files(root=staging)
        manifest = {
            "schema_version": 1,
            "status": "COMPLETE",
            "inputs": input_manifest,
            "outputs": output_manifest,
        }
        write_json_atomic(path=staging / "manifest.json", value=manifest)
        write_json_atomic(
            path=staging / "COMPLETED.json",
            value={
                "status": "COMPLETE",
                "manifest_sha256": sha256_file(path=staging / "manifest.json"),
            },
        )
        os.replace(staging, destination)
    except Exception as error:
        shutil.rmtree(staging, ignore_errors=True)
        if isinstance(error, PublicationError):
            raise
        raise PublicationError(f"Could not publish result {destination}: {error}") from error
    LOGGER.info("Published immutable checksummed result to %s", destination)
    return destination


def _build_input_manifest(*, input_paths: Sequence[Path]) -> tuple[dict[str, Any], ...]:
    """Checksum every declared input authority without silently dropping one.

    Args:
        input_paths: Paths declared as authorities for the result.

    Returns:
        Deterministically ordered input-manifest records.

    Raises:
        PublicationError: If an authority is absent, not a file or changes while read.
    """

    records: list[dict[str, Any]] = []
    sources = sorted({Path(item).expanduser().resolve() for item in input_paths}, key=str)
    for source in sources:
        if not source.is_file():
            raise PublicationError(f"Declared input authority is missing or not a file: {source}")
        try:
            before = source.stat()
            digest = sha256_file(path=source)
            after = source.stat()
        except OSError as error:
            raise PublicationError(
                f"Could not checksum input authority {source}: {error}"
            ) from error
        before_identity = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        )
        after_identity = (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        )
        if before_identity != after_identity:
            raise PublicationError(f"Input authority changed while checksumming: {source}")
        records.append(
            {
                "path": str(source),
                "size_bytes": after.st_size,
                "sha256": digest,
            }
        )
    return tuple(records)


def verify_completed_result(*, result_dir: Path) -> None:
    """Verify completion marker, manifest checksum and every declared output.

    Args:
        result_dir: Candidate completed result directory.

    Raises:
        PublicationError: If any completion invariant fails.
    """

    root = Path(result_dir).expanduser().resolve()
    marker_path = root / "COMPLETED.json"
    manifest_path = root / "manifest.json"
    if not marker_path.is_file() or not manifest_path.is_file():
        raise PublicationError(f"Result lacks completion marker or manifest: {root}")
    marker = read_json(path=marker_path)
    manifest = read_json(path=manifest_path)
    if not isinstance(marker, dict) or marker.get("status") != "COMPLETE":
        raise PublicationError(f"Invalid completion marker: {marker_path}")
    if marker.get("manifest_sha256") != sha256_file(path=manifest_path):
        raise PublicationError(f"Manifest checksum mismatch: {manifest_path}")
    if not isinstance(manifest, dict) or not isinstance(manifest.get("outputs"), list):
        raise PublicationError(f"Invalid output manifest structure: {manifest_path}")
    declared_paths: set[str] = set()
    for record in manifest["outputs"]:
        if not isinstance(record, dict) or not isinstance(record.get("relative_path"), str):
            raise PublicationError(f"Invalid output record in {manifest_path}")
        relative_path = record["relative_path"]
        if relative_path in declared_paths:
            raise PublicationError(f"Duplicate output path in manifest: {relative_path!r}")
        declared_paths.add(relative_path)
        candidate = (root / relative_path).resolve()
        if root not in candidate.parents or not candidate.is_file():
            raise PublicationError(f"Manifest output is missing or unsafe: {candidate}")
        if candidate.stat().st_size != record.get("size_bytes"):
            raise PublicationError(f"Output size mismatch: {candidate}")
        if sha256_file(path=candidate) != record.get("sha256"):
            raise PublicationError(f"Output checksum mismatch: {candidate}")
    expected_paths = declared_paths | {"manifest.json", "COMPLETED.json"}
    actual_paths = {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() or path.is_symlink()
    }
    if actual_paths != expected_paths:
        undeclared = sorted(actual_paths - expected_paths)
        missing = sorted(expected_paths - actual_paths)
        raise PublicationError(
            "Result file inventory differs from the manifest; "
            f"undeclared={undeclared[:10]}, missing={missing[:10]}."
        )


def verify_input_authorities(*, result_dir: Path) -> None:
    """Verify that original local inputs still match a completed manifest.

    This check is used for computational resume, not for portable read-only app
    viewing, where the original inputs are not required to remain present.

    Args:
        result_dir: Completed result directory.

    Raises:
        PublicationError: If an input is missing or its checksum changed.
    """

    root = Path(result_dir).expanduser().resolve()
    manifest = read_json(path=root / "manifest.json")
    if not isinstance(manifest, dict) or not isinstance(manifest.get("inputs"), list):
        raise PublicationError(f"Invalid input manifest structure: {root / 'manifest.json'}")
    for record in manifest["inputs"]:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise PublicationError("Invalid input authority record in manifest.")
        source = Path(record["path"]).expanduser().resolve()
        if not source.is_file():
            raise PublicationError(f"Input authority is no longer available: {source}")
        if source.stat().st_size != record.get("size_bytes"):
            raise PublicationError(f"Input authority size changed: {source}")
        if sha256_file(path=source) != record.get("sha256"):
            raise PublicationError(f"Input authority checksum changed: {source}")


def _write_table(*, table_name: str, records: Sequence[Mapping[str, Any]], table_dir: Path) -> None:
    """Write one canonical result table as TSV and typed Parquet.

    Args:
        table_name: Canonical table name.
        records: Table rows.
        table_dir: Result table directory.
    """

    schema = schema_for(table_name=table_name)
    normalised = tuple(dict(record) for record in records)
    write_tsv_atomic(
        path=table_dir / f"{table_name}.tsv",
        fieldnames=schema.names,
        records=normalised,
    )
    try:
        table = pa.Table.from_pylist(list(normalised), schema=schema)
    except (pa.ArrowException, TypeError, ValueError) as error:
        raise PublicationError(f"Rows violate schema for table {table_name!r}: {error}") from error
    destination = table_dir / f"{table_name}.parquet"
    descriptor, temporary_name = tempfile.mkstemp(
        dir=table_dir, prefix=f".{destination.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        pq.write_table(table, temporary_name, compression="zstd", write_statistics=True)
        os.replace(temporary_name, destination)
    except (OSError, pa.ArrowException) as error:
        Path(temporary_name).unlink(missing_ok=True)
        raise PublicationError(f"Could not publish Parquet table {destination}: {error}") from error


def _copy_assets(*, staging: Path, asset_sources: Mapping[str, Path]) -> None:
    """Copy portable coordinates and human-report assets into staging.

    Args:
        staging: Exact staging result directory.
        asset_sources: Relative destinations mapped to source files.

    Raises:
        PublicationError: If an asset path escapes staging or copying fails.
    """

    for relative_name, source_value in sorted(asset_sources.items()):
        relative = Path(relative_name)
        if relative.is_absolute() or ".." in relative.parts:
            raise PublicationError(f"Unsafe result asset path: {relative_name!r}")
        source = Path(source_value).expanduser().resolve()
        if not source.is_file() or source.stat().st_size == 0:
            raise PublicationError(f"Missing or empty result asset source: {source}")
        destination = (staging / relative).resolve()
        if staging not in destination.parents:
            raise PublicationError(f"Result asset escapes staging directory: {relative_name!r}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            shutil.copyfile(source, destination)
        except OSError as error:
            raise PublicationError(
                f"Could not copy result asset {source} to {destination}: {error}"
            ) from error


def _create_duckdb(*, path: Path, table_dir: Path) -> None:
    """Create a physical DuckDB database from canonical Parquet files.

    Args:
        path: DuckDB destination.
        table_dir: Directory containing canonical Parquet tables.
    """

    try:
        with duckdb.connect(str(path)) as connection:
            for parquet_path in sorted(table_dir.glob("*.parquet")):
                table_name = parquet_path.stem
                schema_for(table_name=table_name)
                escaped_path = str(parquet_path).replace("'", "''")
                connection.execute(
                    f"CREATE TABLE \"{table_name}\" AS SELECT * FROM read_parquet('{escaped_path}')"
                )
            connection.execute(
                "CREATE VIEW signature_evidence AS "
                "SELECT s.*, c.display_name AS comparison_name "
                "FROM signatures s LEFT JOIN comparisons c USING (comparison_id)"
            )
            connection.execute("CHECKPOINT")
    except (duckdb.Error, OSError) as error:
        raise PublicationError(f"Could not create DuckDB database {path}: {error}") from error


def _manifest_files(*, root: Path) -> list[dict[str, Any]]:
    """Describe every current file beneath a staging result.

    Args:
        root: Staging result directory.

    Returns:
        Deterministically ordered checksum records.
    """

    return [
        {
            "relative_path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in sorted((item for item in root.rglob("*") if item.is_file()), key=str)
    ]
