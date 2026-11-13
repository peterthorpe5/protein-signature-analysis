"""Inspect complete audited partitions and select exact unaligned members for handoffs."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from protein_signatures.errors import InputValidationError

from .analysis import index_fasta, write_fasta
from .catalogue import open_catalogue, source_view
from .config import WorkbenchConfig
from .context import open_context

SUBCLUSTER_ID = re.compile(r"sub_[0-9a-f]{24}")


def validate_partition_link(*, sources: dict) -> str:
    """Require the membership and original FASTA of the same completed subcluster job."""
    family, membership = sources.get("family", {}), sources.get("subcluster_membership", {})
    identifier = membership.get("subcluster_id")
    if (
        not isinstance(identifier, str)
        or not SUBCLUSTER_ID.fullmatch(identifier)
        or family.get("name") != "sequences.fasta"
        or membership.get("name") != "membership.parquet"
        or family.get("job_id") != membership.get("job_id")
    ):
        raise InputValidationError("Choose an audited subcluster from the same completed family.")
    return identifier


def filter_subcluster(
    *,
    config: WorkbenchConfig,
    fasta: Path,
    membership: Path,
    subcluster_id: str,
    destination: Path,
) -> None:
    """Export the full selected partition and check every ID and sequence length."""
    with open_catalogue(config=config, with_index=False) as connection:
        source_view(connection=connection, path=str(membership), name="partition_members")
        connection.execute(
            "CREATE TABLE selected AS SELECT sequence_id,sequence_length FROM partition_members "
            "WHERE subcluster_id=?",
            [subcluster_id],
        )
        count = connection.execute(
            "SELECT count(*),count(DISTINCT sequence_id) FROM selected"
        ).fetchone()
        if not count[0] or count[0] != count[1]:
            raise InputValidationError("The selected subcluster is missing or has duplicate IDs.")
        index_fasta(
            connection=connection,
            path=fasta,
            name="family",
            maximum=config.max_subcluster_sequences,
        )
        invalid = connection.execute(
            "SELECT 1 FROM selected LEFT JOIN family USING(sequence_id) "
            "WHERE sequence IS NULL OR sequence_length<>length(sequence) LIMIT 1"
        ).fetchone()
        if invalid:
            raise InputValidationError("Subcluster membership disagrees with the saved sequences.")
        write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM family JOIN selected USING(sequence_id) "
            "ORDER BY sequence_id",
            path=destination,
        )


def compare_partition(
    *,
    config: WorkbenchConfig,
    membership: Path,
    context: Path | None,
    subcluster_id: str,
) -> dict[str, pd.DataFrame]:
    """Compare the full selected subcluster with every saved input member."""
    with open_catalogue(config=config, with_index=False) as connection:
        source_view(connection=connection, path=str(membership), name="partition_members")
        connection.execute(
            "CREATE TABLE scopes AS SELECT *, 'All saved input' AS scope FROM partition_members "
            "UNION ALL SELECT *, 'Selected subcluster' AS scope FROM partition_members "
            "WHERE subcluster_id=?",
            [subcluster_id],
        )
        size = connection.execute("""
            SELECT scope,count(*) AS proteins,count(DISTINCT subcluster_id) AS audited_groups,
                count(DISTINCT native_representative_id) AS native_groups,
                min(sequence_length) AS shortest,median(sequence_length) AS median_length,
                avg(sequence_length) AS mean_length,max(sequence_length) AS longest,
                count(*) FILTER (WHERE alignment_status NOT IN
                    ('CENTROID_SELF','EXACT_THRESHOLDS_PASS')) AS rescued_singletons
            FROM scopes GROUP BY scope ORDER BY scope
        """).fetchdf()
        if len(size) != 2:
            raise InputValidationError("Select one existing audited subcluster.")
        evidence = connection.execute("""
            SELECT scope,alignment_status,count(*) AS proteins,
                count(pident) AS measured_alignments,min(pident) AS minimum_identity_percent,
                median(pident) AS median_identity_percent,
                median(representative_coverage) AS median_centroid_coverage_percent,
                median(member_coverage) AS median_member_coverage_percent,
                median(bitscore) AS median_bitscore,median(evalue) AS median_evalue
            FROM scopes GROUP BY scope,alignment_status ORDER BY scope,alignment_status
        """).fetchdf()
        if context is not None:
            open_context(connection=connection, path=context)
            taxonomy = connection.execute("""
                SELECT scope,coalesce(released_species,'') AS released_species,
                    coalesce(reviewed_species,'') AS reviewed_species,
                    coalesce(reviewed_taxon_id,'') AS reviewed_taxon_id,count(*) AS proteins
                FROM scopes LEFT JOIN family_context USING(sequence_id)
                GROUP BY scope,released_species,reviewed_species,reviewed_taxon_id
                ORDER BY scope,proteins DESC,released_species,reviewed_taxon_id
            """).fetchdf()
        else:
            taxonomy = pd.DataFrame()
    return {
        "Size and sequence lengths": size,
        "Exact audit evidence": evidence,
        "Taxonomy": taxonomy,
    }


def reviewed_taxon_sets(
    *,
    config: WorkbenchConfig,
    membership: Path,
    context: Path,
    groups: list[str],
) -> dict[str, set[str]]:
    """Collect exact reviewed taxon IDs per selected audited group, keeping unknowns out."""
    if not 1 <= len(groups) <= 12 or len(set(groups)) != len(groups):
        raise InputValidationError("Choose 1–12 distinct subclusters for taxon comparison.")
    result: dict[str, set[str]] = {key: set() for key in groups}
    with open_catalogue(config=config, with_index=False) as connection:
        source_view(connection=connection, path=str(membership), name="partition_members")
        open_context(connection=connection, path=context)
        rows = connection.execute(
            "SELECT DISTINCT subcluster_id,reviewed_taxon_id FROM partition_members "
            "JOIN family_context USING(sequence_id) WHERE reviewed_taxon_id IS NOT NULL "
            "AND subcluster_id IN (SELECT unnest(?))",
            [groups],
        ).fetchall()
    for group, taxon in rows:
        result[group].add(taxon)
    return result


def taxon_overlap(*, groups: dict[str, set[str]]) -> tuple[list[dict], list[dict]]:
    """Compute pairwise Jaccard and up to three-group exclusive reviewed-taxon regions."""
    pairwise, regions = [], []
    for left, a in groups.items():
        for right, b in groups.items():
            union, shared = a | b, a & b
            pairwise.append(
                {
                    "subcluster_a": left,
                    "subcluster_b": right,
                    "shared_reviewed_taxa": len(shared),
                    "union_reviewed_taxa": len(union),
                    "jaccard": len(shared) / len(union) if union else None,
                }
            )
    if 2 <= len(groups) <= 3:
        names = list(groups)
        for mask in range(1, 1 << len(names)):
            inside = [name for i, name in enumerate(names) if mask & (1 << i)]
            outside = [name for name in names if name not in inside]
            taxa = set.intersection(*(groups[name] for name in inside))
            for name in outside:
                taxa -= groups[name]
            regions.append({"groups": " + ".join(inside), "exclusive_reviewed_taxa": len(taxa)})
    return pairwise, regions
