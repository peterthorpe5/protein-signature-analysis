"""Deterministic homology-aware discovery and validation partitions."""

from __future__ import annotations

import hashlib

from .errors import InputValidationError
from .models import (
    GroupMembership,
    PartitionAssignment,
    RedundancyClusterMembership,
    SequenceRecord,
)


def assign_partitions(
    *,
    sequences: tuple[SequenceRecord, ...],
    memberships: tuple[GroupMembership, ...],
    validation_fraction: float,
    random_seed: int,
    redundancy_memberships: tuple[RedundancyClusterMembership, ...] = (),
) -> tuple[PartitionAssignment, ...]:
    """Assign whole homology units to discovery or held-out validation.

    OrthoFinder groups, exact-sequence hashes and supplied near-redundancy
    clusters are joined into connected partition blocks. This prevents a copy
    linked by any authority from leaking across partitions.

    Args:
        sequences: Authoritative protein records.
        memberships: Selected completed-OrthoFinder memberships.
        validation_fraction: Fraction of partition units reserved for validation.
        random_seed: Stable integer seed incorporated into hash assignments.
        redundancy_memberships: Optional externally generated near-redundancy clusters.

    Returns:
        One deterministic assignment per protein.

    Raises:
        InputValidationError: If bounds or memberships are invalid.
    """

    if not 0.0 <= validation_fraction <= 0.9:
        raise InputValidationError("validation_fraction must be between 0.0 and 0.9.")
    if random_seed < 0:
        raise InputValidationError("random_seed must be non-negative.")
    by_protein = {item.protein_id: item for item in sequences}
    parent = {protein_id: protein_id for protein_id in by_protein}

    def find(protein_id: str) -> str:
        """Return and compress the current disjoint-set root."""

        while parent[protein_id] != protein_id:
            parent[protein_id] = parent[parent[protein_id]]
            protein_id = parent[protein_id]
        return protein_id

    def union(left: str, right: str) -> None:
        """Join two proteins using a deterministic root identifier."""

        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        smaller, larger = sorted((left_root, right_root))
        parent[larger] = smaller

    exact_groups: dict[str, list[str]] = {}
    for sequence in sequences:
        exact_groups.setdefault(sequence.sequence_sha256, []).append(sequence.protein_id)
    for members in exact_groups.values():
        for protein_id in members[1:]:
            union(members[0], protein_id)

    group_by_protein: dict[str, str] = {}
    members_by_group: dict[str, list[str]] = {}
    for membership in memberships:
        if membership.protein_id not in by_protein:
            raise InputValidationError(
                "OrthoFinder membership references protein absent from the FASTA: "
                f"{membership.protein_id!r}"
            )
        composite_group = "|".join(
            (
                membership.run_id,
                membership.group_type,
                membership.hierarchy_node,
                membership.group_id,
            )
        )
        previous = group_by_protein.get(membership.protein_id)
        if previous is not None and previous != composite_group:
            raise InputValidationError(
                f"Protein {membership.protein_id!r} occurs in multiple partition groups."
            )
        group_by_protein[membership.protein_id] = composite_group
        members_by_group.setdefault(composite_group, []).append(membership.protein_id)
    for members in members_by_group.values():
        for protein_id in members[1:]:
            union(members[0], protein_id)

    redundancy_by_protein: dict[str, str] = {}
    members_by_redundancy_cluster: dict[str, list[str]] = {}
    for membership in redundancy_memberships:
        if membership.protein_id not in by_protein:
            raise InputValidationError(
                "Redundancy membership references protein absent from the FASTA: "
                f"{membership.protein_id!r}"
            )
        previous = redundancy_by_protein.get(membership.protein_id)
        if previous is not None and previous != membership.cluster_id:
            raise InputValidationError(
                f"Protein {membership.protein_id!r} occurs in multiple redundancy clusters."
            )
        redundancy_by_protein[membership.protein_id] = membership.cluster_id
        members_by_redundancy_cluster.setdefault(membership.cluster_id, []).append(
            membership.protein_id
        )
    for members in members_by_redundancy_cluster.values():
        for protein_id in members[1:]:
            union(members[0], protein_id)

    components: dict[str, list[str]] = {}
    for protein_id in sorted(by_protein):
        components.setdefault(find(protein_id), []).append(protein_id)
    rows: list[PartitionAssignment] = []
    for members in components.values():
        group_tokens = sorted(
            {group_by_protein[item] for item in members if item in group_by_protein}
        )
        redundancy_tokens = sorted(
            {redundancy_by_protein[item] for item in members if item in redundancy_by_protein}
        )
        exact_tokens = sorted({by_protein[item].sequence_sha256 for item in members})
        digest_input = "\n".join(
            [
                *(f"ORTHOFINDER:{item}" for item in group_tokens),
                *(f"REDUNDANCY:{item}" for item in redundancy_tokens),
                *(f"SEQUENCE:{item}" for item in exact_tokens),
            ]
        )
        if not group_tokens and not redundancy_tokens:
            key = exact_tokens[0]
            unit = "EXACT_SEQUENCE"
        elif len(group_tokens) == 1 and not redundancy_tokens:
            key = group_tokens[0]
            unit = "ORTHOFINDER_GROUP"
        else:
            key = f"PB_{hashlib.sha256(digest_input.encode('utf-8')).hexdigest()}"
            unit = "REDUNDANCY_BLOCK"
        partition = _partition_for_key(
            key=f"{unit}:{key}",
            validation_fraction=validation_fraction,
            random_seed=random_seed,
        )
        for protein_id in members:
            rows.append(
                PartitionAssignment(
                    protein_id=protein_id,
                    partition=partition,
                    partition_unit=unit,
                    partition_key=key,
                )
            )
    return tuple(sorted(rows, key=lambda item: item.protein_id))


def _partition_for_key(*, key: str, validation_fraction: float, random_seed: int) -> str:
    """Map one stable unit key into a named partition.

    Args:
        key: Stable partition-unit key.
        validation_fraction: Reserved validation fraction.
        random_seed: Stable integer seed.

    Returns:
        ``DISCOVERY`` or ``VALIDATION``.
    """

    if validation_fraction == 0.0:
        return "DISCOVERY"
    digest = hashlib.sha256(f"{random_seed}:{key}".encode("utf-8")).digest()
    proportion = int.from_bytes(digest[:8], byteorder="big") / float(2**64)
    return "VALIDATION" if proportion < validation_fraction else "DISCOVERY"
