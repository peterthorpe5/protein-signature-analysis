"""Portable DeepClust result import, separate from the mounted source release."""

from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from protein_signatures.checksums import sha256_file
from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic

from .catalogue import identifiers, normalise_sequence
from .jobs import MAX_UPLOAD_BYTES

LOGGER = logging.getLogger(__name__)


def result_reader(*, source: Path | bytes) -> pq.ParquetFile:
    """Open a bounded self-contained Parquet search export.

    Args:
        source: Uploaded bytes or the private uploaded snapshot.

    Returns:
        Reader with bounded metadata and no external column files.
    """
    size = source.stat().st_size if isinstance(source, Path) else len(source)
    if not 1 <= size <= MAX_UPLOAD_BYTES:
        raise InputValidationError("DeepClust export must be at most 20 MiB.")
    try:
        reader = pq.ParquetFile(
            source=source if isinstance(source, Path) else pa.BufferReader(source),
            thrift_string_size_limit=16_000_000,
            thrift_container_size_limit=1_000_000,
        )
    except (pa.ArrowException, OSError) as error:
        raise InputValidationError("Cannot open this DeepClust Parquet export.") from error
    metadata = reader.metadata
    names = reader.schema_arrow.names
    if (
        not 1 <= metadata.num_rows <= 200_000
        or not 1 <= metadata.num_columns <= 128
        or (len(set(names)) != len(names) or metadata.num_row_groups > 4096)
    ):
        reader.close()
        raise InputValidationError("Export 1–200,000 rows with at most 128 unique fields.")
    decoded = 0
    for index in range(metadata.num_row_groups):
        group = metadata.row_group(index)
        decoded += group.total_byte_size
        if decoded > 256 * 1024 * 1024 or any(
            group.column(column).file_path for column in range(group.num_columns)
        ):
            reader.close()
            raise InputValidationError("Use a smaller self-contained DeepClust export.")
    return reader


def import_result(*, source: Path, parameters: dict, output: Path) -> None:
    """Validate exact identifiers and extract one complete portable family.

    Args:
        source: Small search-result Parquet file, not a database partition.
        parameters: Identifier, sequence and optional exact cluster columns.
        output: Owned job destination.
    """
    id_column = parameters.get("id_column", "sequence_id")
    sequence_column = parameters.get("sequence_column", "sequence")
    cluster_column = parameters.get("cluster_column", "")
    cluster_id = parameters.get("cluster_id", "")
    columns = [id_column, sequence_column] + ([cluster_column] if cluster_column else [])
    if not all(isinstance(column, str) and column for column in columns) or (
        len(set(columns)) != len(columns) or bool(cluster_column) != bool(cluster_id)
    ):
        raise InputValidationError(
            "Choose distinct sequence/ID fields and an optional exact cluster."
        )
    records, size, duplicates = {}, 0, 0
    with result_reader(source=source) as reader:
        if not set(columns).issubset(reader.schema_arrow.names):
            raise InputValidationError("Selected columns are absent from this DeepClust export.")
        for batch in reader.iter_batches(batch_size=1024, columns=columns):
            for row in batch.to_pylist():
                if cluster_column and str(row[cluster_column]) != cluster_id:
                    continue
                identifier = identifiers(values=[row[id_column]])[0]
                sequence = normalise_sequence(value=row[sequence_column])
                if identifier in records:
                    if records[identifier] != sequence:
                        raise InputValidationError(f"Conflicting sequences for {identifier}.")
                    duplicates += 1
                    continue
                size += len(identifier.encode()) + len(sequence) + 3
                if size > MAX_UPLOAD_BYTES:
                    raise InputValidationError(
                        "Select a smaller family; FASTA is limited to 20 MiB."
                    )
                records[identifier] = sequence
    if not records:
        raise InputValidationError("No exact proteins matched the selected export rows.")
    (output / "family.fasta").write_text(
        data="".join(f">{key}\n{seq}\n" for key, seq in records.items()), encoding="utf-8"
    )
    write_json_atomic(
        path=output / "import_summary.json",
        value={
            "source_sha256": sha256_file(path=source),
            "parameters": parameters,
            "proteins": len(records),
            "duplicate_rows_collapsed": duplicates,
            "scope": "Uploaded result subset; not verified complete DeepClust cluster membership",
        },
    )
    LOGGER.info("Imported %d proteins from a portable DeepClust result", len(records))
