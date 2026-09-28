"""Build comparison-specific analysis cohorts from the control-match audit."""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Any

from .errors import InputValidationError
from .models import ComparisonDefinition, PartitionAssignment

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class MatchedComparisonCohort:
    """Target proteins with controls and controls matched to those target units."""

    target_protein_ids: frozenset[str]
    background_protein_ids: frozenset[str]
    target_unit_count: int
    unmatched_target_unit_count: int
    unmatched_target_protein_count: int
    control_unit_count: int

    def to_record(self) -> dict[str, int]:
        """Summarise the analysed subset for the campaign provenance record."""

        return {
            "target_protein_count": len(self.target_protein_ids),
            "control_protein_count": len(self.background_protein_ids),
            "matched_target_unit_count": self.target_unit_count,
            "excluded_unmatched_target_unit_count": self.unmatched_target_unit_count,
            "excluded_unmatched_target_protein_count": self.unmatched_target_protein_count,
            "matched_control_unit_count": self.control_unit_count,
        }


def derive_matched_comparison_cohorts(
    *,
    comparisons: tuple[ComparisonDefinition, ...],
    label_memberships: tuple[dict[str, str], ...],
    partitions: tuple[PartitionAssignment, ...],
    control_matching_audit: tuple[dict[str, Any], ...],
) -> dict[str, MatchedComparisonCohort]:
    """Restrict each audited comparison to supported target and control units.

    The matching policy can cover several target classes with one background.
    Controls for a particular comparison therefore come only from the units
    matched to target units carrying that comparison's target label. A target
    without any matched control remains in the label authority and audit, but
    cannot enter the association or prediction cohort for this comparison.

    Args:
        comparisons: Explicit analysis comparisons.
        label_memberships: Expanded positive label memberships.
        partitions: Analysis homology/redundancy blocks.
        control_matching_audit: Verified evidence-bundle matching rows.

    Returns:
        Audited comparison IDs mapped to their restricted cohorts. Comparisons
        with no audited background retain their standard reviewed-label cohort.

    Raises:
        InputValidationError: If the audit contradicts the campaign authorities
            or only some requested backgrounds have match provenance.
    """

    unit_by_protein = {item.protein_id: item.partition_key for item in partitions}
    labels_by_protein: dict[str, set[str]] = defaultdict(set)
    for row in label_memberships:
        labels_by_protein[row["protein_id"]].add(row["label_id"])
    matches_by_background: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    seen_control_units: dict[str, set[str]] = defaultdict(set)
    for row in control_matching_audit:
        background = str(row["background_label_id"])
        target_id = str(row["target_protein_id"])
        target_unit = str(row["target_unit_id"])
        target_labels = set(str(row["target_label_ids"]).split("|"))
        if not background or not target_unit or not target_labels - {""}:
            raise InputValidationError("Control audit has a missing target or background identity.")
        if unit_by_protein.get(target_id) != target_unit:
            raise InputValidationError(
                f"Control audit target block differs from campaign partitions: {target_id!r}."
            )
        if not labels_by_protein[target_id].intersection(target_labels):
            raise InputValidationError(
                f"Control audit target has no corresponding positive label: {target_id!r}."
            )
        matched_units = matches_by_background[background][target_unit]
        status = row["status"]
        if status == "UNMATCHED":
            if row["control_unit_id"] or row["control_protein_id"]:
                raise InputValidationError("Unmatched control audit row names a control.")
            continue
        if status != "MATCHED":
            raise InputValidationError(f"Unknown control audit status: {status!r}.")
        control_id = str(row["control_protein_id"])
        control_unit = str(row["control_unit_id"])
        if not control_unit or unit_by_protein.get(control_id) != control_unit:
            raise InputValidationError(
                f"Control audit block differs from campaign partitions: {control_id!r}."
            )
        if background not in labels_by_protein[control_id]:
            raise InputValidationError(
                f"Matched control lacks its background label: {control_id!r}."
            )
        if control_unit == target_unit or control_unit in seen_control_units[background]:
            raise InputValidationError(
                f"Control audit reused a block within {background!r}: {control_unit!r}."
            )
        seen_control_units[background].add(control_unit)
        matched_units.add(control_unit)

    audited_backgrounds = set(matches_by_background)
    result: dict[str, MatchedComparisonCohort] = {}
    for comparison in comparisons:
        requested_backgrounds = set(comparison.background_label_ids)
        audited = requested_backgrounds & audited_backgrounds
        if not audited:
            LOGGER.info(
                "Comparison %s has no audited matched background; using reviewed cohorts",
                comparison.comparison_id,
            )
            continue
        if audited != requested_backgrounds:
            raise InputValidationError(
                f"Comparison {comparison.comparison_id!r} mixes audited and unaudited "
                "background labels."
            )
        target_proteins = {
            protein_id
            for protein_id, labels in labels_by_protein.items()
            if labels.intersection(comparison.target_label_ids)
        }
        target_units = {unit_by_protein[protein_id] for protein_id in target_proteins}
        controls_by_target: dict[str, set[str]] = defaultdict(set)
        for background in sorted(audited):
            for target_unit, control_units in matches_by_background[background].items():
                controls_by_target[target_unit].update(control_units)
        matched_target_units = {unit_id for unit_id in target_units if controls_by_target[unit_id]}
        matched_control_units = (
            set().union(*(controls_by_target[unit_id] for unit_id in matched_target_units))
            if matched_target_units
            else set()
        )
        target_subset = frozenset(
            protein_id
            for protein_id in target_proteins
            if unit_by_protein[protein_id] in matched_target_units
        )
        background_subset = frozenset(
            protein_id
            for protein_id, labels in labels_by_protein.items()
            if labels.intersection(requested_backgrounds)
            and unit_by_protein[protein_id] in matched_control_units
        )
        cohort = MatchedComparisonCohort(
            target_protein_ids=target_subset,
            background_protein_ids=background_subset,
            target_unit_count=len(matched_target_units),
            unmatched_target_unit_count=len(target_units - matched_target_units),
            unmatched_target_protein_count=len(target_proteins - target_subset),
            control_unit_count=len(matched_control_units),
        )
        LOGGER.info(
            "Matched cohort comparison=%s target_units=%d excluded_unmatched_units=%d "
            "control_units=%d",
            comparison.comparison_id,
            cohort.target_unit_count,
            cohort.unmatched_target_unit_count,
            cohort.control_unit_count,
        )
        result[comparison.comparison_id] = cohort
    return result
