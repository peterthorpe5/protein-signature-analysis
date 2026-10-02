"""Read-only DuckDB access and result validation for the application."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

from protein_signatures.errors import InputValidationError
from protein_signatures.exports import dataframe_to_tsv_bytes, dataframe_to_xlsx_bytes
from protein_signatures.publication import verify_completed_result
from protein_signatures.schemas import table_schemas

_RELATION_PATTERN = re.compile(
    r'\b(?:FROM|JOIN)\s+(?:"([A-Za-z_][A-Za-z0-9_]*)"|([A-Za-z_][A-Za-z0-9_]*))',
    flags=re.IGNORECASE,
)
_MAX_DOWNLOAD_BYTES = 512 * 1024 * 1024
_VERIFICATION_CACHE_VERSION = 1
LOGGER = logging.getLogger(__name__)
_DUCKDB_READ_ONLY_CONFIG = {
    "allow_community_extensions": "false",
    "allow_persistent_secrets": "false",
    "allow_unsigned_extensions": "false",
    "autoinstall_known_extensions": "false",
    "autoload_known_extensions": "false",
    "enable_external_access": "false",
    "lock_configuration": "true",
}


@dataclass(frozen=True)
class DownloadAsset:
    """One verified, ready-to-download result asset."""

    file_format: str
    relative_path: str
    row_count: int
    payload: bytes
    mime_type: str
    complete: bool = True


def result_inventory_identity(*, resource: Path) -> str:
    """Fingerprint result file metadata for interactive verification caching.

    The fingerprint is deliberately cheap to calculate. A new or modified file
    changes the key and triggers a fresh full checksum verification; the CLI's
    explicit ``verify`` command continues to verify every byte on every call.

    Args:
        resource: Completed result directory or physical DuckDB path.

    Returns:
        Deterministic digest of path names, types, sizes and change times.

    Raises:
        InputValidationError: If the result directory cannot be inspected.
    """

    candidate = Path(resource).expanduser().resolve()
    root = candidate.parent if candidate.is_file() else candidate
    digest = hashlib.sha256()
    walk_errors: list[OSError] = []
    try:
        for directory, names, files in os.walk(root, followlinks=False, onerror=walk_errors.append):
            names.sort()
            for name in sorted((*names, *files)):
                path = Path(directory) / name
                stat = path.lstat()
                digest.update(str(path.relative_to(root)).encode("utf-8"))
                identity = (
                    f"\0{stat.st_dev}:{stat.st_ino}:{stat.st_mode}:{stat.st_size}:"
                    f"{stat.st_mtime_ns}:{stat.st_ctime_ns}\n"
                )
                digest.update(identity.encode("ascii"))
    except (OSError, ValueError) as error:
        raise InputValidationError(f"Could not inspect result file inventory: {error}") from error
    if walk_errors:
        raise InputValidationError(f"Could not inspect result file inventory: {walk_errors[0]}")
    return digest.hexdigest()


def app_database_path(*, resource: Path) -> Path:
    """Check the local result paths before starting the viewer.

    Args:
        resource: Completed result directory or its physical DuckDB file.

    Returns:
        Canonical DuckDB path; file contents are not verified at this step.

    Raises:
        InputValidationError: If required paths are missing or unsupported.
    """

    candidate = Path(resource).expanduser().resolve()
    root = candidate.parent if candidate.is_file() else candidate
    if not root.is_dir():
        raise InputValidationError(f"Completed result directory does not exist: {root}")
    for filename in ("COMPLETED.json", "manifest.json"):
        if not (root / filename).is_file():
            raise InputValidationError(f"Result lacks completion marker or manifest: {root}")
    database = root / "protein_signatures.duckdb"
    if not database.is_file() or database.stat().st_size == 0:
        raise InputValidationError(f"Completed result lacks its DuckDB database: {database}")
    if candidate.is_file() and candidate != database:
        raise InputValidationError(f"Unsupported resource file; expected {database}")
    return database


def resolve_database(*, resource: Path) -> Path:
    """Verify every published output and resolve the result database.

    Args:
        resource: Completed result directory or its physical DuckDB file.

    Returns:
        Verified DuckDB path.

    Raises:
        InputValidationError: If any completion invariant fails.
    """

    try:
        database = app_database_path(resource=resource)
        verify_completed_result(result_dir=database.parent)
    except Exception as error:
        raise InputValidationError(
            f"App resource is not a valid completed result: {error}"
        ) from error
    return database


def resolve_database_for_app(
    *, resource: Path, inventory_identity: str, cache_dir: Path | None = None
) -> Path:
    """Reuse a full verification only while every result file remains unchanged.

    The cache lives outside the immutable result. A change to a file name,
    inode, type, size, modification time or change time invalidates its entry. The
    explicit CLI verify command always performs a full checksum verification.

    Args:
        resource: Completed result directory or canonical database file.
        inventory_identity: Fresh metadata fingerprint of the entire result.
        cache_dir: Optional verification-cache directory for isolated tests.

    Returns:
        Verified result database path.

    Raises:
        InputValidationError: If verification fails or the result changes
            while its checksums are being checked.
    """

    database = app_database_path(resource=resource)
    if not inventory_identity:
        raise InputValidationError("The result inventory fingerprint is missing.")
    cache_root = (
        Path(cache_dir)
        if cache_dir is not None
        else Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
        / "protein_signature_analysis"
    )
    cache_path = cache_root / (
        hashlib.sha256(str(database.parent).encode("utf-8")).hexdigest() + ".json"
    )
    expected = {
        "version": _VERIFICATION_CACHE_VERSION,
        "result_dir": str(database.parent),
        "inventory_identity": inventory_identity,
    }
    try:
        if not cache_root.is_symlink() and not cache_path.is_symlink():
            if cache_path.is_file() and cache_path.stat().st_size <= 4096:
                if json.loads(cache_path.read_text(encoding="utf-8")) == expected:
                    if result_inventory_identity(resource=database) == inventory_identity:
                        LOGGER.debug("Reused unchanged result verification for %s", database.parent)
                        return database
    except (OSError, UnicodeError, ValueError) as error:
        LOGGER.debug("Ignored unusable result verification cache: %s", error)

    verified = resolve_database(resource=resource)
    if result_inventory_identity(resource=verified) != inventory_identity:
        raise InputValidationError("Result files changed during checksum verification.")
    try:
        if cache_root.is_symlink():
            raise OSError(f"Verification cache directory is a symbolic link: {cache_root}")
        cache_root.mkdir(mode=0o700, parents=True, exist_ok=True)
        temporary_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=cache_root,
                prefix=".verification-",
                suffix=".json",
                delete=False,
            ) as handle:
                temporary_path = Path(handle.name)
                os.chmod(temporary_path, 0o600)
                json.dump(expected, handle, sort_keys=True)
            os.replace(temporary_path, cache_path)
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)
    except OSError as error:
        LOGGER.warning("Could not save reusable result verification: %s", error)
    return verified


def query_dataframe(*, database: Path, sql: str, parameters: tuple[Any, ...] = ()) -> pd.DataFrame:
    """Run one read-only query and return a Pandas data frame.

    Args:
        database: Verified physical DuckDB path.
        sql: Internal single-statement SELECT query over canonical result tables.
        parameters: Bound query parameters.

    Returns:
        Query result data frame.

    Raises:
        InputValidationError: If SQL is not read-only or execution fails.
    """

    _validate_canonical_select(sql=sql)
    try:
        with duckdb.connect(
            str(database),
            read_only=True,
            config=_DUCKDB_READ_ONLY_CONFIG,
        ) as connection:
            return connection.execute(sql, parameters).fetchdf()
    except duckdb.Error as error:
        raise InputValidationError(f"DuckDB query failed: {error}") from error


def _validate_canonical_select(*, sql: str) -> None:
    """Reject statements that can reach beyond canonical result relations.

    Args:
        sql: Candidate SQL statement used by the application.

    Raises:
        InputValidationError: If the statement is not one simple canonical SELECT.
    """

    statement = str(sql).strip()
    if not re.match(r"^SELECT\b", statement, flags=re.IGNORECASE):
        raise InputValidationError("The app backend accepts read-only SELECT queries only.")
    if any(token in statement for token in (";", "--", "/*", "*/", "\x00")):
        raise InputValidationError("The app backend accepts one comment-free SELECT statement.")
    relations = {
        next(value for value in match.groups() if value)
        for match in _RELATION_PATTERN.finditer(statement)
    }
    allowed = set(table_schemas()) | {"signature_evidence"}
    if not relations:
        raise InputValidationError("The app query must read a canonical result relation.")
    unknown = sorted(relations - allowed)
    if unknown:
        raise InputValidationError(
            "The app query references a non-canonical relation: " + ", ".join(unknown)
        )


def table_count(*, database: Path, table_name: str) -> int:
    """Return the row count for one canonical table.

    Args:
        database: Verified physical DuckDB path.
        table_name: Canonical table name.

    Returns:
        Non-negative row count.

    Raises:
        InputValidationError: If the name is not canonical.
    """

    if table_name not in table_schemas():
        raise InputValidationError(f"Unknown canonical result table: {table_name!r}")
    frame = query_dataframe(database=database, sql=f'SELECT count(*) AS n FROM "{table_name}"')
    return int(frame.iloc[0]["n"])


def distinct_values(*, database: Path, table_name: str, column_name: str) -> tuple[str, ...]:
    """Return ordered non-empty text values from a canonical table column.

    Args:
        database: Verified DuckDB path.
        table_name: Canonical table name.
        column_name: Canonical column name.

    Returns:
        Ordered values.

    Raises:
        InputValidationError: If table or column is not canonical.
    """

    schemas = table_schemas()
    if table_name not in schemas or column_name not in schemas[table_name].names:
        raise InputValidationError(
            f"Unknown canonical table/column: {table_name!r}/{column_name!r}"
        )
    frame = query_dataframe(
        database=database,
        sql=(
            f'SELECT DISTINCT "{column_name}" AS value FROM "{table_name}" '
            f'WHERE "{column_name}" IS NOT NULL AND "{column_name}" <> \'\' ORDER BY value'
        ),
    )
    return tuple(str(value) for value in frame["value"].tolist())


def load_metadata(*, database: Path) -> dict[str, Any]:
    """Load run metadata beside a verified database.

    Args:
        database: Verified DuckDB path.

    Returns:
        Metadata mapping.

    Raises:
        InputValidationError: If metadata is malformed.
    """

    path = database.parent / "run_metadata.json"
    try:
        with path.open(mode="r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError(f"Could not load run metadata {path}: {error}") from error
    if not isinstance(value, dict):
        raise InputValidationError(f"Run metadata must be a JSON object: {path}")
    return value


def canonical_table_names() -> tuple[str, ...]:
    """Return every canonical result-table name in publication order.

    Returns:
        Stable table names shared by publication and the application.
    """

    return tuple(table_schemas())


def canonical_table_preview(
    *, database: Path, table_name: str, limit: int = 500, offset: int = 0
) -> pd.DataFrame:
    """Load one bounded page from a canonical result table.

    Args:
        database: Verified physical DuckDB path.
        table_name: Canonical table name.
        limit: Number of rows to return, from 1 through 5,000.
        offset: Zero-based row offset.

    Returns:
        Canonically ordered result columns and the requested row page.

    Raises:
        InputValidationError: If the table or page bounds are invalid.
    """

    if table_name not in table_schemas():
        raise InputValidationError(f"Unknown canonical result table: {table_name!r}")
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 5_000:
        raise InputValidationError("Canonical table preview limit must be from 1 through 5,000.")
    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise InputValidationError("Canonical table preview offset must be a non-negative integer.")
    return query_dataframe(
        database=database,
        sql=f'SELECT * FROM "{table_name}" LIMIT ? OFFSET ?',
        parameters=(limit, offset),
    )


def load_canonical_table_assets(*, database: Path, table_name: str) -> tuple[DownloadAsset, ...]:
    """Load practical published downloads for one canonical table.

    Args:
        database: Verified physical DuckDB path.
        table_name: Canonical table name.

    Returns:
        Complete TSV and XLSX for manageable tables. Large tables return their
        compact XLSX index; their complete TSV.GZ remains available by path and
        through filtered exports without being loaded into app memory.

    Raises:
        InputValidationError: If inventory rows or assets violate the result contract.
    """

    if table_name not in table_schemas():
        raise InputValidationError(f"Unknown canonical result table: {table_name!r}")
    inventory = _read_report_inventory(database=database)
    rows = [
        row
        for row in inventory
        if row.get("content_id") == table_name and row.get("section") != "99_final_results"
    ]
    expected_count = table_count(database=database, table_name=table_name)
    assets: list[DownloadAsset] = []
    complete_rows = [row for row in rows if row.get("asset_kind") == "TABLE"]
    if complete_rows:
        selected = complete_rows
    else:
        selected = [row for row in rows if row.get("asset_kind") == "TABLE_INDEX"]
    for row in selected:
        file_format = str(row.get("file_format", ""))
        mime_type = (
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            if file_format == "XLSX"
            else "text/tab-separated-values"
        )
        try:
            row_count = int(row.get("row_count", ""))
        except (TypeError, ValueError) as error:
            raise InputValidationError(
                f"Canonical table {table_name!r} has an invalid inventory row count."
            ) from error
        complete = row.get("asset_kind") == "TABLE"
        if complete and row_count != expected_count:
            raise InputValidationError(
                f"Canonical table {table_name!r} inventory count differs from DuckDB."
            )
        relative_path = str(row.get("relative_path", ""))
        payload = _read_result_asset(database=database, relative_path=relative_path)
        assets.append(
            DownloadAsset(
                file_format=file_format,
                relative_path=relative_path,
                row_count=row_count,
                payload=payload,
                mime_type=mime_type,
                complete=complete,
            )
        )
    if not assets:
        raise InputValidationError(
            f"Canonical table {table_name!r} has no practical download or summary asset."
        )
    return tuple(assets)


def canonical_table_storage(*, database: Path, table_name: str) -> tuple[dict[str, Any], ...]:
    """Describe complete canonical files without reading them into memory.

    Args:
        database: Verified physical DuckDB path.
        table_name: Canonical table name.

    Returns:
        Parquet and TSV/TSV.GZ path, size and row-count records.

    Raises:
        InputValidationError: If the table or files are invalid.
    """

    if table_name not in table_schemas():
        raise InputValidationError(f"Unknown canonical result table: {table_name!r}")
    root = Path(database).expanduser().resolve().parent
    candidates = sorted((root / "tables").glob(f"{table_name}.*"))
    row_count = table_count(database=database, table_name=table_name)
    records: list[dict[str, Any]] = []
    for path in candidates:
        if path.name.endswith(".parquet"):
            file_format = "PARQUET"
            if path.is_dir():
                parts = tuple(sorted(path.glob("part-*.parquet")))
                if not parts:
                    raise InputValidationError(
                        f"Canonical Parquet dataset is empty for {table_name!r}."
                    )
                size_bytes = sum(part.stat().st_size for part in parts)
            else:
                size_bytes = path.stat().st_size
        elif path.name.endswith(".tsv.gz"):
            file_format = "TSV.GZ"
            size_bytes = path.stat().st_size
        elif path.name.endswith(".tsv"):
            file_format = "TSV"
            size_bytes = path.stat().st_size
        else:
            continue
        records.append(
            {
                "file_format": file_format,
                "relative_path": str(path.relative_to(root)),
                "size_bytes": size_bytes,
                "row_count": row_count,
            }
        )
    if {record["file_format"] for record in records} & {"TSV", "TSV.GZ"} == set() or not any(
        record["file_format"] == "PARQUET" for record in records
    ):
        raise InputValidationError(f"Canonical storage files are incomplete for {table_name!r}.")
    return tuple(records)


def filtered_feature_exports(
    *,
    database: Path,
    feature_types: tuple[str, ...] = (),
    protein_id: str = "",
    feature_id_contains: str = "",
    maximum_rows: int = 100_000,
) -> tuple[DownloadAsset, ...]:
    """Build bounded TSV and Excel downloads from filtered feature memberships.

    Args:
        database: Verified physical DuckDB path.
        feature_types: Optional exact feature-type filters.
        protein_id: Optional exact protein identifier.
        feature_id_contains: Optional literal feature-identifier substring.
        maximum_rows: Hard result cap from 1 through 250,000 rows.

    Returns:
        Matching TSV and formatted XLSX payloads.

    Raises:
        InputValidationError: If filters are unsafe, absent or too broad.
    """

    if (
        isinstance(maximum_rows, bool)
        or not isinstance(maximum_rows, int)
        or not 1 <= maximum_rows <= 250_000
    ):
        raise InputValidationError("Filtered feature export limit must be 1 through 250,000.")
    types = tuple(
        dict.fromkeys(str(value).strip() for value in feature_types if str(value).strip())
    )
    protein = str(protein_id).strip()
    feature_text = str(feature_id_contains).strip()
    if not types and not protein and not feature_text:
        raise InputValidationError("Select at least one feature filter before exporting.")
    clauses: list[str] = []
    parameters: list[Any] = []
    if types:
        placeholders = ",".join("?" for _value in types)
        clauses.append(f"feature_type IN ({placeholders})")
        parameters.extend(types)
    if protein:
        clauses.append("protein_id = ?")
        parameters.append(protein)
    if feature_text:
        clauses.append("position(? in feature_id) > 0")
        parameters.append(feature_text)
    where_sql = " AND ".join(clauses)
    count_frame = query_dataframe(
        database=database,
        sql=f"SELECT count(*) AS row_count FROM features WHERE {where_sql}",
        parameters=tuple(parameters),
    )
    matched_rows = int(count_frame.iloc[0]["row_count"])
    if matched_rows > maximum_rows:
        raise InputValidationError(
            f"Filtered feature selection contains {matched_rows:,} rows and exceeds the "
            f"{maximum_rows:,}-row export limit; narrow the filters."
        )
    parameters.append(maximum_rows)
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT * FROM features WHERE "
            + where_sql
            + ' ORDER BY feature_type, feature_id, protein_id, "start" LIMIT ?'
        ),
        parameters=tuple(parameters),
    )
    if len(frame) != matched_rows:
        raise InputValidationError(
            "Filtered feature export row count changed during its immutable DuckDB query."
        )
    tsv = dataframe_to_tsv_bytes(frame=frame)
    xlsx = dataframe_to_xlsx_bytes(frame=frame, title="Filtered Feature Memberships")
    return (
        DownloadAsset(
            file_format="TSV",
            relative_path="filtered_features.tsv",
            row_count=len(frame),
            payload=tsv,
            mime_type="text/tab-separated-values",
        ),
        DownloadAsset(
            file_format="XLSX",
            relative_path="filtered_features.xlsx",
            row_count=len(frame),
            payload=xlsx,
            mime_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        ),
    )


def load_report_inventory_assets(*, database: Path) -> tuple[DownloadAsset, ...]:
    """Load the report inventory itself as TSV and formatted XLSX.

    Args:
        database: Verified physical DuckDB path.

    Returns:
        TSV then XLSX inventory assets.
    """

    rows = _read_report_inventory(database=database)
    base = "analysis/00_run_information/report_inventory"
    return (
        DownloadAsset(
            file_format="TSV",
            relative_path=f"{base}.tsv",
            row_count=len(rows),
            payload=_read_result_asset(database=database, relative_path=f"{base}.tsv"),
            mime_type="text/tab-separated-values",
        ),
        DownloadAsset(
            file_format="XLSX",
            relative_path=f"{base}.xlsx",
            row_count=len(rows),
            payload=_read_result_asset(database=database, relative_path=f"{base}.xlsx"),
            mime_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        ),
    )


def _read_report_inventory(*, database: Path) -> tuple[dict[str, str], ...]:
    """Parse the verified human-report inventory beside one database.

    Args:
        database: Verified physical DuckDB path.

    Returns:
        Strict inventory rows.

    Raises:
        InputValidationError: If the inventory is absent or malformed.
    """

    relative_path = "analysis/00_run_information/report_inventory.tsv"
    payload = _read_result_asset(database=database, relative_path=relative_path)
    try:
        text = payload.decode("utf-8")
        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        required = {
            "section",
            "asset_kind",
            "content_id",
            "file_format",
            "relative_path",
            "row_count",
            "status",
            "description",
        }
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise InputValidationError("Report inventory headings violate the output contract.")
        rows = tuple(dict(row) for row in reader)
    except UnicodeError as error:
        raise InputValidationError("Report inventory is not valid UTF-8.") from error
    if not rows or any(row.get("status") != "COMPLETE" for row in rows):
        raise InputValidationError("Report inventory is empty or contains incomplete assets.")
    return rows


def _read_result_asset(*, database: Path, relative_path: str) -> bytes:
    """Read one bounded result-relative regular file without path escape.

    Args:
        database: Canonical database within the completed result.
        relative_path: Portable result-relative path.

    Returns:
        Complete asset bytes.

    Raises:
        InputValidationError: If the path or file is unsafe or too large.
    """

    root = Path(database).expanduser().resolve().parent
    relative = Path(relative_path)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise InputValidationError(f"Unsafe result asset path: {relative_path!r}")
    candidate = (root / relative).resolve()
    if root not in candidate.parents or not candidate.is_file() or candidate.is_symlink():
        raise InputValidationError(f"Result asset is missing or unsafe: {relative_path!r}")
    size = candidate.stat().st_size
    if size > _MAX_DOWNLOAD_BYTES:
        raise InputValidationError(
            f"Result asset exceeds the 512 MiB application download limit: {relative_path!r}"
        )
    try:
        return candidate.read_bytes()
    except OSError as error:
        raise InputValidationError(f"Could not read result asset: {relative_path!r}") from error
