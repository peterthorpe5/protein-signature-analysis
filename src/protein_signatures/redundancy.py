"""Exact and imported near-redundancy evidence for leakage-safe partitions."""

from __future__ import annotations

import logging
from pathlib import Path

from .errors import InputValidationError
from .io_utils import iter_tsv
from .models import RedundancyClusterMembership, SequenceRecord
from .validation import parse_optional_float, validate_identifier, validate_text

LOGGER = logging.getLogger(__name__)
REDUNDANCY_FIELDS = (
    "protein_id",
    "cluster_id",
    "method",
    "method_version",
    "identity_threshold",
    "coverage_threshold",
    "evidence_reference",
)


def derive_exact_sequence_clusters(
    *, sequences: tuple[SequenceRecord, ...]
) -> tuple[RedundancyClusterMembership, ...]:
    """Create a checksum-defined exact-sequence cluster for every protein.

    Args:
        sequences: Authoritative sequences.

    Returns:
        Ordered exact-sequence memberships, including singleton clusters.
    """

    return tuple(
        RedundancyClusterMembership(
            protein_id=record.protein_id,
            cluster_id=f"EXACT_{record.sequence_sha256}",
            cluster_type="EXACT_SEQUENCE",
            method="SHA256_EXACT_SEQUENCE",
            method_version="1",
            identity_threshold=1.0,
            coverage_threshold=1.0,
            evidence_reference=record.sequence_sha256,
        )
        for record in sorted(sequences, key=lambda item: item.protein_id)
    )


def read_redundancy_clusters(
    *, path: Path, protein_ids: frozenset[str]
) -> tuple[RedundancyClusterMembership, ...]:
    """Read a versioned near-redundancy membership authority.

    Args:
        path: Membership TSV.
        protein_ids: Authoritative campaign protein identifiers.

    Returns:
        Ordered near-redundancy memberships.

    Raises:
        InputValidationError: If rows are malformed, unknown or duplicated.
    """

    rows: list[RedundancyClusterMembership] = []
    observed: set[str] = set()
    for row in iter_tsv(path=path, required_fields=REDUNDANCY_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        if protein_id not in protein_ids:
            raise InputValidationError(
                f"Redundancy cluster references protein absent from FASTA: {protein_id!r}"
            )
        if protein_id in observed:
            raise InputValidationError(
                f"Protein occurs more than once in the near-redundancy authority: {protein_id!r}"
            )
        observed.add(protein_id)
        identity = parse_optional_float(
            value=row["identity_threshold"],
            field_name="identity_threshold",
            minimum=0.0,
            maximum=1.0,
        )
        coverage = parse_optional_float(
            value=row["coverage_threshold"],
            field_name="coverage_threshold",
            minimum=0.0,
            maximum=1.0,
        )
        if identity is None or coverage is None:
            raise InputValidationError(
                "Near-redundancy identity_threshold and coverage_threshold must be populated."
            )
        rows.append(
            RedundancyClusterMembership(
                protein_id=protein_id,
                cluster_id=validate_identifier(value=row["cluster_id"], field_name="cluster_id"),
                cluster_type="NEAR_REDUNDANCY",
                method=validate_identifier(value=row["method"], field_name="method"),
                method_version=validate_text(
                    value=row["method_version"], field_name="method_version"
                ),
                identity_threshold=float(identity),
                coverage_threshold=float(coverage),
                evidence_reference=validate_text(
                    value=row["evidence_reference"],
                    field_name="evidence_reference",
                ),
            )
        )
    result = tuple(sorted(rows, key=lambda item: (item.cluster_id, item.protein_id)))
    LOGGER.info(
        "Loaded %d proteins across %d supplied near-redundancy clusters",
        len(result),
        len({item.cluster_id for item in result}),
    )
    return result
