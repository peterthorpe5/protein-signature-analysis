"""Bounded seed selection from immutable protein-signature results."""

from __future__ import annotations

from pathlib import Path

from protein_signature_app.backend import query_dataframe
from protein_signatures.errors import InputValidationError

from .catalogue import identifiers
from .jobs import MAX_UPLOAD_BYTES


def export_seeds(
    *, database: Path, selection: str, values: list[str], maximum: int = 1000
) -> bytes:
    """Export exact proteins, local HOG members, labelled proteins or feature members.

    Args:
        database: Verified read-only reference database.
        selection: protein, group, label or feature.
        values: Exact canonical IDs to select.
        maximum: Maximum complete seed set; exceeding it is an error.

    Returns:
        Complete FASTA snapshot suitable for submitting a job.

    Raises:
        InputValidationError: If IDs are missing, ambiguous in scope or exceed bounds.
    """
    requested = identifiers(values=values)
    clauses = {
        "protein": ("proteins", "protein_id"),
        "group": ("orthofinder_memberships", "group_id"),
        "label": ("label_memberships", "label_id"),
        "feature": ("features", "feature_id"),
    }
    if selection not in clauses or isinstance(maximum, bool) or not 1 <= maximum <= 10_000:
        raise InputValidationError("Choose a supported seed selector and a limit from 1 to 10000.")
    table, field = clauses[selection]
    placeholders = ",".join("?" for _ in requested)
    found = query_dataframe(
        database=database,
        sql=f"SELECT DISTINCT {field} FROM {table} WHERE {field} IN ({placeholders})",
        parameters=requested,
    )
    missing = set(requested) - set(found[field].astype(str))
    if missing:
        raise InputValidationError(f"Seed selector IDs were not found: {sorted(missing)[:10]}")
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT protein_id,sequence FROM proteins WHERE protein_id IN ("
            f"SELECT protein_id FROM {table} WHERE {field} IN ({placeholders})) "
            "ORDER BY protein_id LIMIT ?"
        ),
        parameters=(*requested, maximum + 1),
    )
    if frame.empty or len(frame) > maximum:
        raise InputValidationError(
            "Seed selection is empty or exceeds the limit; choose fewer groups or exact proteins."
        )
    payload = "".join(
        f">{row.protein_id}\n{row.sequence}\n" for row in frame.itertuples(index=False)
    ).encode()
    if len(payload) > MAX_UPLOAD_BYTES:
        raise InputValidationError("Seed FASTA exceeds the 20 MiB submission limit.")
    return payload
