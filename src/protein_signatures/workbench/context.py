"""Sequence-bound group and taxonomy evidence for workbench handoffs."""

from __future__ import annotations

import csv
import io
import logging
from collections.abc import Iterable
from pathlib import Path

from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic

from .analysis import index_fasta
from .catalogue import column, open_catalogue, source_view
from .config import WorkbenchConfig

LOGGER = logging.getLogger(__name__)
CONTEXT_FILE = "family_context.tsv"
CONTEXT_FIELDS = (
    "sequence_id",
    "sequence_sha256",
    "source_group",
    "group_type",
    "run_id",
    "hierarchy_node",
    "released_species",
    "reviewed_species",
    "reviewed_taxon_id",
)
CONTEXT_FORMAT = ",".join(CONTEXT_FIELDS)


def context_bytes(*, rows: Iterable[dict]) -> bytes:
    """Encode sequence context as UTF-8 TSV with an explicit ordered header."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=CONTEXT_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


def reviewed_taxonomy(*, connection: object, path: str) -> None:
    """Expose only consistent positive reviewed IDs for exact sequence identifiers."""
    if not path:
        connection.execute(
            "CREATE TABLE reviewed_taxonomy(sequence_id VARCHAR,reviewed_species VARCHAR,"
            "reviewed_taxon_id BIGINT)"
        )
        return
    columns = set(source_view(connection=connection, path=path, name="taxonomy_source"))
    if not {"sequence_id", "species", "taxon_id"} <= columns:
        raise InputValidationError("Reviewed taxonomy needs sequence_id, species and taxon_id.")
    reviewed = "WHERE mapping_status='REVIEWED'" if "mapping_status" in columns else ""
    # Accept a legacy operator-supplied reviewed table without a status column.
    connection.execute(
        "CREATE VIEW reviewed_taxonomy AS SELECT CAST(sequence_id AS VARCHAR) AS sequence_id,"
        "CASE WHEN count(DISTINCT species)=1 THEN min(species) END AS reviewed_species,"
        "min(taxon) AS reviewed_taxon_id FROM (SELECT sequence_id,species,"
        "CASE WHEN regexp_full_match(CAST(taxon_id AS VARCHAR),'[1-9][0-9]*') "
        "THEN try_cast(taxon_id AS BIGINT) END AS taxon FROM taxonomy_source "
        f"{reviewed}) GROUP BY sequence_id "
        "HAVING count(DISTINCT taxon)=1 AND count(taxon)=count(*)"
    )


def open_context(*, connection: object, path: Path) -> None:
    """Load a single, unambiguous metadata record per exact sequence and identifier."""
    columns = source_view(connection=connection, path=str(path), name="family_context")
    if tuple(columns) != CONTEXT_FIELDS:
        raise InputValidationError("Family context TSV does not have the expected header.")
    invalid = connection.execute(
        "SELECT 1 FROM family_context WHERE sequence_id IS NULL "
        "OR NOT coalesce(regexp_full_match(sequence_sha256,'[0-9a-f]{64}'),false) "
        "OR (reviewed_taxon_id IS NOT NULL AND ("
        "NOT regexp_full_match(reviewed_taxon_id,'[1-9][0-9]*') "
        "OR try_cast(reviewed_taxon_id AS BIGINT) IS NULL)) "
        "UNION ALL SELECT 1 FROM family_context GROUP BY sequence_id HAVING count(*)>1 LIMIT 1"
    ).fetchone()
    if invalid:
        raise InputValidationError("Family context contains malformed or ambiguous evidence.")


def reviewed_species(*, connection: object, path: str) -> None:
    """Read exact names from the operator's E3-compatible reviewed species mapping."""
    if not path:
        connection.execute(
            "CREATE TABLE reviewed_species(released_species VARCHAR,reviewed_species VARCHAR,"
            "reviewed_taxon_id BIGINT)"
        )
        return
    columns = source_view(connection=connection, path=path, name="species_taxonomy_source")
    label = column(columns=columns, choices=("source_species_name", "workflow_species_label"))
    taxon = column(columns=columns, choices=("taxon_id", "ncbi_taxon_id"))
    name = column(columns=columns, choices=("canonical_species_name", "accepted_species_name"))
    reviewed = "AND upper(trim(mapping_status))='REVIEWED'" if "mapping_status" in columns else ""
    connection.execute(
        "CREATE VIEW reviewed_species AS SELECT released_species,"
        "CASE WHEN count(DISTINCT species)=1 THEN min(species) END AS reviewed_species,"
        "min(taxon_id) AS reviewed_taxon_id FROM ("
        f"SELECT {label} AS released_species,{name} AS species,"
        f"try_cast({taxon} AS BIGINT) AS taxon_id FROM species_taxonomy_source "
        f"WHERE regexp_full_match(CAST({taxon} AS VARCHAR),'[1-9][0-9]*') {reviewed}) "
        "GROUP BY released_species HAVING count(DISTINCT taxon_id)=1 AND count(taxon_id)=count(*)"
    )


def save_context(*, connection: object, query: str, output: Path, sequence_count: int) -> None:
    """Save bounded display labels and a streaming sequence context ledger."""
    connection.execute("CREATE TABLE matched_context AS " + query)
    connection.sql(f"SELECT {CONTEXT_FORMAT} FROM matched_context ORDER BY sequence_id").write_csv(
        file_name=str(output / CONTEXT_FILE), sep="\t", header=True
    )
    groups = connection.execute(
        "SELECT DISTINCT source_group FROM matched_context WHERE source_group <> '' "
        "ORDER BY source_group LIMIT 12"
    ).fetchall()
    matched = connection.execute("SELECT count(*) FROM matched_context").fetchone()[0]
    write_json_atomic(
        path=output / "family_context.json",
        value={
            "source_groups": [row[0] for row in groups],
            "metadata_matched": matched,
            "output_sequences": sequence_count,
            "match_policy": "Exact original protein identifier and normalised sequence SHA-256.",
        },
    )
    LOGGER.info("Retained family context for %d of %d sequences", matched, sequence_count)


def carry_context(
    *,
    config: WorkbenchConfig,
    context: Path,
    sequences: Path,
    output: Path,
) -> None:
    """Carry only the exact sequence and ID matches into the next job's ledger."""
    with open_catalogue(config=config, with_index=False) as connection:
        count = index_fasta(
            connection=connection, path=sequences, name="family", maximum=config.max_sequences
        )
        open_context(connection=connection, path=context)
        save_context(
            connection=connection,
            query="SELECT c.* FROM family_context c JOIN family USING(sequence_id,sequence_sha256)",
            output=output,
            sequence_count=count,
        )
