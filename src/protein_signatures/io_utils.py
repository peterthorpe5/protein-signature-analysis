"""Safe TSV, JSON and text input/output helpers."""

from __future__ import annotations

import csv
import gzip
import json
import os
import tempfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, TextIO

from .errors import InputValidationError, PublicationError


@contextmanager
def open_text(*, path: Path) -> Iterator[TextIO]:
    """Open plain or gzip-compressed UTF-8 text for reading.

    Args:
        path: Existing text path, optionally ending in ``.gz``.

    Yields:
        Open text handle.

    Raises:
        InputValidationError: If the input is missing or empty.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty input file: {source}")
    opener = gzip.open if source.suffix == ".gz" else open
    try:
        with opener(source, mode="rt", encoding="utf-8", newline="") as handle:
            yield handle
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Could not read UTF-8 input {source}: {error}") from error


def iter_tsv(
    *, path: Path, required_fields: Sequence[str], allow_empty: bool = False
) -> Iterator[dict[str, str]]:
    """Yield strictly shaped TSV records.

    Args:
        path: Plain or gzip-compressed TSV.
        required_fields: Required unique column names.
        allow_empty: Whether a header-only table is valid.

    Yields:
        Ordered row dictionaries.

    Raises:
        InputValidationError: If the header or any row is malformed.
    """

    required = tuple(required_fields)
    if not required or len(required) != len(set(required)):
        raise InputValidationError("required_fields must contain unique column names.")
    with open_text(path=path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headings = reader.fieldnames
        if headings is None:
            raise InputValidationError(f"TSV has no header: {Path(path).resolve()}")
        if len(headings) != len(set(headings)) or any(not heading for heading in headings):
            raise InputValidationError(f"TSV has duplicate or empty headings: {headings}")
        missing = [field for field in required if field not in headings]
        if missing:
            raise InputValidationError(
                f"TSV is missing required fields {missing}: {Path(path).resolve()}"
            )
        count = 0
        for row_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise InputValidationError(
                    f"Malformed field count at row {row_number}: {Path(path).resolve()}"
                )
            count += 1
            yield {key: value for key, value in row.items()}
        if count == 0 and not allow_empty:
            raise InputValidationError(f"TSV contains no data rows: {Path(path).resolve()}")


def write_tsv_atomic(
    *,
    path: Path,
    fieldnames: Sequence[str],
    records: Iterable[Mapping[str, Any]],
) -> int:
    """Write a TSV through a same-directory temporary file and atomic rename.

    Args:
        path: Final table destination.
        fieldnames: Ordered output columns.
        records: Row mappings.

    Returns:
        Number of written data rows.

    Raises:
        PublicationError: If output fields or records are invalid.
    """

    fields = tuple(fieldnames)
    if not fields or len(fields) != len(set(fields)):
        raise PublicationError("fieldnames must be a non-empty unique sequence.")
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        text=True,
    )
    count = 0
    try:
        with os.fdopen(descriptor, mode="w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                extrasaction="raise",
                lineterminator="\n",
            )
            writer.writeheader()
            for record in records:
                missing = [field for field in fields if field not in record]
                if missing:
                    raise PublicationError(f"Output record is missing fields: {missing}")
                writer.writerow(record)
                count += 1
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    except (OSError, ValueError, PublicationError) as error:
        Path(temporary_name).unlink(missing_ok=True)
        if isinstance(error, PublicationError):
            raise
        raise PublicationError(f"Could not publish TSV {destination}: {error}") from error
    return count


def write_json_atomic(*, path: Path, value: Any) -> None:
    """Write formatted JSON atomically.

    Args:
        path: Final JSON destination.
        value: JSON-compatible data.

    Raises:
        PublicationError: If serialisation or publication fails.
    """

    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        text=True,
    )
    try:
        with os.fdopen(descriptor, mode="w", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    except (OSError, TypeError, ValueError) as error:
        Path(temporary_name).unlink(missing_ok=True)
        raise PublicationError(f"Could not publish JSON {destination}: {error}") from error


def write_text_atomic(*, path: Path, text: str) -> None:
    """Write UTF-8 text through a same-directory atomic rename.

    Args:
        path: Final text destination.
        text: Unicode text to publish.

    Raises:
        PublicationError: If the value is not text or publication fails.
    """

    if not isinstance(text, str):
        raise PublicationError("Atomic text output must be a string.")
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        text=True,
    )
    try:
        with os.fdopen(descriptor, mode="w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    except OSError as error:
        Path(temporary_name).unlink(missing_ok=True)
        raise PublicationError(f"Could not publish text {destination}: {error}") from error


def read_json(*, path: Path) -> Any:
    """Read one UTF-8 JSON document.

    Args:
        path: Existing JSON input.

    Returns:
        Decoded JSON value.

    Raises:
        InputValidationError: If the JSON is invalid.
    """

    try:
        with open_text(path=path) as handle:
            return json.load(handle)
    except json.JSONDecodeError as error:
        raise InputValidationError(f"Invalid JSON in {Path(path).resolve()}: {error}") from error
