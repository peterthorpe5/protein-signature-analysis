"""Prepare auditable starter authorities from a generic sequence catalogue TSV."""

from __future__ import annotations

import hashlib
import logging
import os
import re
import shutil
import tempfile
from pathlib import Path

from .checksums import sha256_file
from .errors import InputValidationError, PublicationError
from .fasta import read_protein_fasta
from .io_utils import iter_tsv, write_json_atomic, write_text_atomic, write_tsv_atomic
from .validation import validate_identifier, validate_text

LOGGER = logging.getLogger(__name__)
_UNIPROT_ACCESSION = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]"
    r"(?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
_MAX_CATALOGUE_TEXT_LENGTH = 32_767


def prepare_catalogue(
    *,
    catalogue_path: Path,
    output_dir: Path,
    id_column: str,
    sequence_column: str,
    name_column: str = "",
    proposed_category_column: str = "",
    starter_label_id: str = "e3:associated:unknown",
) -> Path:
    """Convert a TSV sequence catalogue into conservative starter authorities.

    No catalogue annotation is promoted to a reviewed positive label. Every
    generated assignment is ``UNMAPPED`` and therefore excluded from training
    until a curator edits the copied template deliberately.

    Args:
        catalogue_path: Plain or gzip-compressed source TSV.
        output_dir: New output bundle directory.
        id_column: Source column containing protein identifiers.
        sequence_column: Source column containing full protein sequences.
        name_column: Optional human-readable description column.
        proposed_category_column: Optional unreviewed category retained for audit.
        starter_label_id: Valid profile label used only as a curation placeholder.

    Returns:
        Atomically published starter bundle directory.

    Raises:
        InputValidationError: If source rows are missing, duplicated or malformed.
        PublicationError: If the destination already exists or publication fails.
    """

    source = Path(catalogue_path).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise PublicationError(f"Catalogue starter output already exists: {destination}")
    columns = tuple(
        column
        for column in (id_column, sequence_column, name_column, proposed_category_column)
        if column
    )
    if len(columns) != len(set(columns)):
        raise InputValidationError("Catalogue source column names must be distinct.")
    if not id_column or not sequence_column:
        raise InputValidationError("id_column and sequence_column are required.")
    label_id = validate_identifier(value=starter_label_id, field_name="starter_label_id")
    source_checksum = sha256_file(path=source)
    rows = _read_catalogue_rows(
        path=source,
        id_column=id_column,
        sequence_column=sequence_column,
        name_column=name_column,
        proposed_category_column=proposed_category_column,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        _write_starter_files(
            staging=staging,
            rows=rows,
            starter_label_id=label_id,
            source_name=source.name,
            source_checksum=source_checksum,
        )
        outputs = _output_inventory(root=staging)
        write_json_atomic(
            path=staging / "PREPARED.json",
            value={
                "schema_version": 1,
                "status": "COMPLETE",
                "source_path": str(source),
                "source_sha256": source_checksum,
                "protein_count": len(rows),
                "training_eligible_assignment_count": 0,
                "outputs": outputs,
            },
        )
        os.replace(staging, destination)
    except (OSError, UnicodeError, InputValidationError, PublicationError) as error:
        shutil.rmtree(staging, ignore_errors=True)
        if isinstance(error, (InputValidationError, PublicationError)):
            raise
        raise PublicationError(f"Could not publish catalogue starter bundle: {error}") from error
    LOGGER.info(
        "Prepared %d conservative catalogue records at %s; manual label review is required",
        len(rows),
        destination,
    )
    return destination


def _read_catalogue_rows(
    *,
    path: Path,
    id_column: str,
    sequence_column: str,
    name_column: str,
    proposed_category_column: str,
) -> tuple[dict[str, str], ...]:
    """Read and validate sequence-bearing catalogue rows.

    Args:
        path: Source TSV.
        id_column: Protein identifier field.
        sequence_column: Protein sequence field.
        name_column: Optional description field.
        proposed_category_column: Optional unreviewed category field.

    Returns:
        Ordered normalised catalogue rows.
    """

    required = tuple(
        column
        for column in (id_column, sequence_column, name_column, proposed_category_column)
        if column
    )
    result: list[dict[str, str]] = []
    seen: dict[str, str] = {}
    for source_row, raw in enumerate(iter_tsv(path=path, required_fields=required), start=2):
        protein_id = validate_identifier(value=raw[id_column], field_name=id_column)
        sequence = "".join(raw[sequence_column].split()).upper().rstrip("*")
        if not sequence:
            raise InputValidationError(
                f"Catalogue protein {protein_id!r} has no sequence at row {source_row}."
            )
        previous = seen.get(protein_id)
        if previous is not None:
            message = "conflicting sequences" if previous != sequence else "a duplicate identifier"
            raise InputValidationError(
                f"Catalogue has {message} for {protein_id!r} at row {source_row}."
            )
        seen[protein_id] = sequence
        description = (
            validate_text(
                value=raw[name_column],
                field_name=name_column,
                allow_empty=True,
                maximum_length=_MAX_CATALOGUE_TEXT_LENGTH,
            )
            if name_column
            else ""
        )
        proposed_category = (
            validate_text(
                value=raw[proposed_category_column],
                field_name=proposed_category_column,
                allow_empty=True,
                maximum_length=_MAX_CATALOGUE_TEXT_LENGTH,
            )
            if proposed_category_column
            else ""
        )
        result.append(
            {
                "protein_id": protein_id,
                "description": description,
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
                "source_row": str(source_row),
                "proposed_category": proposed_category,
            }
        )
    return tuple(sorted(result, key=lambda item: item["protein_id"]))


def _write_starter_files(
    *,
    staging: Path,
    rows: tuple[dict[str, str], ...],
    starter_label_id: str,
    source_name: str,
    source_checksum: str,
) -> None:
    """Write FASTA, curation, acquisition and provenance authorities.

    Args:
        staging: New staging directory.
        rows: Validated catalogue rows.
        starter_label_id: Placeholder profile label.
        source_name: Source filename for provenance.
        source_checksum: Source SHA-256 digest.
    """

    fasta = "".join(
        f">{row['protein_id']}{' ' + row['description'] if row['description'] else ''}\n"
        f"{_wrap_sequence(sequence=row['sequence'])}"
        for row in rows
    )
    fasta_path = staging / "proteins.faa"
    write_text_atomic(path=fasta_path, text=fasta)
    read_protein_fasta(path=fasta_path)
    write_tsv_atomic(
        path=staging / "label_assignments.REVIEW_REQUIRED.tsv",
        fieldnames=(
            "protein_id",
            "label_id",
            "curation_status",
            "evidence_status",
            "evidence_source",
            "evidence_reference",
            "component_role",
            "curation_reason",
        ),
        records=(
            {
                "protein_id": row["protein_id"],
                "label_id": starter_label_id,
                "curation_status": "UNMAPPED",
                "evidence_status": "UNREVIEWED_SOURCE",
                "evidence_source": source_name,
                "evidence_reference": source_checksum,
                "component_role": "UNKNOWN",
                "curation_reason": "Manual class and role review required before training.",
            }
            for row in rows
        ),
    )
    write_tsv_atomic(
        path=staging / "alphafold_accessions.REVIEW_REQUIRED.tsv",
        fieldnames=("protein_id", "uniprot_accession"),
        records=(
            {"protein_id": row["protein_id"], "uniprot_accession": row["protein_id"]}
            for row in rows
            if _is_uniprot_accession(value=row["protein_id"])
        ),
    )
    write_tsv_atomic(
        path=staging / "catalogue_audit.tsv",
        fieldnames=(
            "protein_id",
            "source_row",
            "sequence_sha256",
            "proposed_category",
            "promotion_status",
        ),
        records=(
            {
                "protein_id": row["protein_id"],
                "source_row": row["source_row"],
                "sequence_sha256": row["sequence_sha256"],
                "proposed_category": row["proposed_category"],
                "promotion_status": "NOT_PROMOTED_REVIEW_REQUIRED",
            }
            for row in rows
        ),
    )


def _wrap_sequence(*, sequence: str, width: int = 80) -> str:
    """Wrap a sequence deterministically for FASTA publication.

    Args:
        sequence: Non-empty sequence text.
        width: Positive line width.

    Returns:
        Wrapped sequence ending with a newline.

    Raises:
        InputValidationError: If the value or width is invalid.
    """

    if not sequence or width < 1:
        raise InputValidationError("FASTA wrapping requires a sequence and positive width.")
    return "".join(
        f"{sequence[index : index + width]}\n" for index in range(0, len(sequence), width)
    )


def _is_uniprot_accession(*, value: str) -> bool:
    """Return whether text has a canonical UniProt accession shape.

    Args:
        value: Candidate accession.

    Returns:
        ``True`` for canonical six- or ten-character accessions.
    """

    return _UNIPROT_ACCESSION.fullmatch(value.upper()) is not None


def _output_inventory(*, root: Path) -> list[dict[str, object]]:
    """Build a deterministic checksum inventory before the completion marker.

    Args:
        root: Staging bundle root.

    Returns:
        Output checksum records.
    """

    return [
        {
            "relative_path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in sorted((item for item in root.rglob("*") if item.is_file()), key=str)
    ]
