"""Search published HOGs and snapshot exact sequence and annotation evidence."""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from protein_signature_app.backend import query_dataframe
from protein_signatures.errors import InputValidationError

from .catalogue import identifiers, normalise_sequence, open_catalogue
from .config import WorkbenchConfig
from .context import context_bytes, reviewed_species, reviewed_taxonomy
from .jobs import MAX_UPLOAD_BYTES

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True, kw_only=True)
class HogScope:
    """Identify one published HOG, including its run and hierarchy node."""

    group_id: str
    run_id: str
    hierarchy_node: str


def search_hogs(*, database: Path, text: str, maximum: int = 50) -> list[HogScope]:
    """Search HOG identifiers literally, keeping same-named runs distinct."""
    if len(text) > 200 or not 1 <= maximum <= 200:
        raise InputValidationError("HOG search exceeds its display limits.")
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT DISTINCT group_id,run_id,hierarchy_node FROM orthofinder_memberships "
            "WHERE group_type='HOG' AND contains(lower(group_id),lower(?)) "
            "ORDER BY group_id,run_id,hierarchy_node LIMIT ?"
        ),
        parameters=(text.strip(), maximum),
    )
    return [HogScope(**record) for record in frame.fillna("").to_dict(orient="records")]


def hog_members(
    *,
    database: Path,
    scope: HogScope,
    maximum: int,
    selected_ids: list[str] | None = None,
) -> pd.DataFrame:
    """Read one explicit HOG scope and retain members without available sequences."""
    if not 1 <= maximum <= 100_000:
        raise InputValidationError("Invalid HOG member limit.")
    values = identifiers(values=selected_ids, maximum=maximum) if selected_ids else ()
    filter_ids = "AND m.protein_id IN (" + ",".join("?" for _ in values) + ") " if values else ""
    frame = query_dataframe(
        database=database,
        sql=(
            "SELECT DISTINCT m.protein_id,m.species_label,p.sequence,p.sequence_sha256 "
            "FROM orthofinder_memberships m LEFT JOIN proteins p USING(protein_id) "
            "WHERE m.group_type='HOG' AND m.group_id=? AND coalesce(m.run_id,'')=? "
            "AND coalesce(m.hierarchy_node,'')=? "
            f"{filter_ids}ORDER BY m.protein_id LIMIT ?"
        ),
        parameters=(scope.group_id, scope.run_id, scope.hierarchy_node, *values, maximum + 1),
    ).fillna("")
    if values and set(values) != set(frame.protein_id):
        raise InputValidationError("An explicitly selected protein is outside the released HOG.")
    if frame.empty or len(frame) > maximum:
        raise InputValidationError(
            "HOG selection is empty or exceeds the limit. Enter a smaller exact member selection."
        )
    if frame.protein_id.duplicated().any():
        raise InputValidationError("This HOG has conflicting rows for a member identifier.")
    return frame


def hog_snapshot(
    *,
    config: WorkbenchConfig,
    scope: HogScope,
    members: pd.DataFrame,
) -> tuple[bytes, bytes]:
    """Export released proteins and metadata bound to each exact sequence hash."""
    sequences, rows = [], []
    available = members.loc[members.sequence.ne("")]
    if available.empty:
        raise InputValidationError("No selected HOG member has a published sequence.")
    with open_catalogue(config=config, with_index=False) as connection:
        reviewed_taxonomy(connection=connection, path=config.taxonomy)
        reviewed_species(connection=connection, path=config.species_taxonomy)
        taxonomy = {
            row[0]: row[1:]
            for row in connection.execute(
                "SELECT * FROM reviewed_taxonomy WHERE sequence_id IN (SELECT unnest(?))",
                [available.protein_id.tolist()],
            ).fetchall()
        }
        species = {
            row[0]: row[1:]
            for row in connection.execute(
                "SELECT * FROM reviewed_species WHERE released_species IN (SELECT unnest(?))",
                [available.species_label.unique().tolist()],
            ).fetchall()
        }
    byte_count = 0
    for record in available.itertuples(index=False):
        identifier = identifiers(values=[record.protein_id])[0]
        sequence = normalise_sequence(value=record.sequence)
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if record.sequence_sha256 and record.sequence_sha256 != digest:
            raise InputValidationError(f"Published sequence checksum mismatch for {identifier}.")
        name_by_id, taxon_by_id = taxonomy.get(identifier, ("", None))
        name_by_label, taxon_by_label = species.get(record.species_label, ("", None))
        if taxon_by_id and taxon_by_label and taxon_by_id != taxon_by_label:
            reviewed_name, taxon = "", None
            LOGGER.warning("Reviewed taxon sources disagree for %s", identifier)
        else:
            reviewed_name, taxon = name_by_id or name_by_label, taxon_by_id or taxon_by_label
        row = f">{identifier}\n{sequence}\n"
        byte_count += len(row.encode())
        if byte_count > MAX_UPLOAD_BYTES:
            raise InputValidationError("HOG FASTA exceeds 20 MiB. Select fewer exact members.")
        sequences.append(row)
        rows.append(
            {
                "sequence_id": identifier,
                "sequence_sha256": digest,
                "source_group": scope.group_id,
                "group_type": "HOG",
                "run_id": scope.run_id,
                "hierarchy_node": scope.hierarchy_node,
                "released_species": record.species_label,
                "reviewed_species": reviewed_name,
                "reviewed_taxon_id": taxon,
            }
        )
    context = context_bytes(rows=rows)
    if len(context) > MAX_UPLOAD_BYTES:
        raise InputValidationError("HOG member metadata exceeds 20 MiB. Select fewer members.")
    LOGGER.info("Snapshotted %d published members of %s", len(rows), scope.group_id)
    return "".join(sequences).encode(), context


def recorded_domains(*, database: Path, identifier: str, maximum: int = 2_000) -> list[dict]:
    """Read recorded domain intervals for one selected canonical protein ID."""
    value = identifiers(values=[identifier])[0]
    if not 1 <= maximum <= 10_000:
        raise InputValidationError("Invalid published domain display limit.")
    frame = query_dataframe(
        database=database,
        sql=(
            'SELECT domain_authority,domain_id,domain_name,"start","end",'
            "evidence_source,evidence_reference FROM domain_hits WHERE protein_id=? "
            'ORDER BY "start","end",domain_id LIMIT ?'
        ),
        parameters=(value, maximum + 1),
    )
    if len(frame) > maximum:
        raise InputValidationError("This protein exceeds the domain feature display limit.")
    return frame.to_dict(orient="records")
