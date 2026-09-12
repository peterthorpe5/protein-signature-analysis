"""Streaming protein FASTA parsing and amino-acid feature generation."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Iterator
from pathlib import Path

from .errors import InputValidationError
from .io_utils import open_text
from .models import SequenceRecord
from .validation import validate_identifier, validate_text

LOGGER = logging.getLogger(__name__)
_VALID_AMINO_ACIDS = frozenset("ABCDEFGHIKLMNOPQRSTUVWXYZ*-?")


def read_protein_fasta(*, path: Path) -> tuple[SequenceRecord, ...]:
    """Read and validate a protein FASTA without altering identifiers.

    Args:
        path: Plain or gzip-compressed FASTA file.

    Returns:
        Ordered immutable sequence records.

    Raises:
        InputValidationError: If records are malformed or duplicated.
    """

    records = tuple(iter_protein_fasta(path=path))
    LOGGER.info("Loaded %d protein sequences from %s", len(records), Path(path).resolve())
    return records


def iter_protein_fasta(*, path: Path) -> Iterator[SequenceRecord]:
    """Yield validated protein records from a FASTA file.

    Args:
        path: Plain or gzip-compressed FASTA file.

    Yields:
        Validated sequence records.

    Raises:
        InputValidationError: If FASTA syntax or protein content is invalid.
    """

    seen: set[str] = set()
    current_header: str | None = None
    parts: list[str] = []
    with open_text(path=path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    record = _make_record(header=current_header, sequence="".join(parts))
                    _reject_duplicate(protein_id=record.protein_id, seen=seen)
                    yield record
                current_header = line[1:].strip()
                if not current_header:
                    raise InputValidationError(f"Empty FASTA header at line {line_number}.")
                parts = []
            else:
                if current_header is None:
                    raise InputValidationError(
                        f"Sequence text precedes the first FASTA header at line {line_number}."
                    )
                parts.append("".join(line.split()).upper())
    if current_header is None:
        raise InputValidationError(f"FASTA contains no records: {Path(path).resolve()}")
    record = _make_record(header=current_header, sequence="".join(parts))
    _reject_duplicate(protein_id=record.protein_id, seen=seen)
    yield record


def sequence_kmers(*, sequence: str, length: int) -> frozenset[str]:
    """Return unique overlapping amino-acid k-mers.

    Args:
        sequence: Validated protein sequence.
        length: Positive k-mer length.

    Returns:
        Unique k-mers in the sequence.

    Raises:
        InputValidationError: If the length is invalid.
    """

    if length < 1:
        raise InputValidationError(f"K-mer length must be positive: {length}")
    if length > len(sequence):
        return frozenset()
    return frozenset(
        sequence[index : index + length] for index in range(len(sequence) - length + 1)
    )


def _make_record(*, header: str, sequence: str) -> SequenceRecord:
    """Build one validated sequence record.

    Args:
        header: FASTA header without ``>``.
        sequence: Concatenated upper-case sequence.

    Returns:
        Validated sequence record.

    Raises:
        InputValidationError: If the sequence or header is invalid.
    """

    protein_id, separator, description = header.partition(" ")
    protein_id = validate_identifier(value=protein_id, field_name="FASTA protein identifier")
    description = validate_text(
        value=description if separator else "",
        field_name=f"FASTA description for {protein_id}",
        allow_empty=True,
    )
    if not sequence:
        raise InputValidationError(f"Protein {protein_id!r} has an empty sequence.")
    invalid = sorted(set(sequence) - _VALID_AMINO_ACIDS)
    if invalid:
        raise InputValidationError(
            f"Protein {protein_id!r} contains unsupported residue symbols: {invalid}"
        )
    if "*" in sequence[:-1]:
        raise InputValidationError(f"Protein {protein_id!r} contains an internal stop symbol.")
    clean_sequence = sequence[:-1] if sequence.endswith("*") else sequence
    digest = hashlib.sha256(clean_sequence.encode("ascii")).hexdigest()
    return SequenceRecord(
        protein_id=protein_id,
        description=description,
        sequence=clean_sequence,
        sequence_length=len(clean_sequence),
        sequence_sha256=digest,
    )


def _reject_duplicate(*, protein_id: str, seen: set[str]) -> None:
    """Reject and otherwise register a FASTA identifier.

    Args:
        protein_id: Candidate identifier.
        seen: Mutable set scoped to one FASTA read.

    Raises:
        InputValidationError: If the identifier was already observed.
    """

    if protein_id in seen:
        raise InputValidationError(f"Duplicate FASTA identifier: {protein_id!r}")
    seen.add(protein_id)
