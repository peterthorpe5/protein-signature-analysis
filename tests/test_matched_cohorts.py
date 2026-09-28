"""Matched controls must be tied to the target blocks of each comparison."""

from __future__ import annotations

import pytest

from protein_signatures.associations import analyse_feature_associations
from protein_signatures.errors import InputValidationError
from protein_signatures.matched_cohorts import derive_matched_comparison_cohorts
from protein_signatures.models import (
    AnalysisSettings,
    ComparisonDefinition,
    FeatureRecord,
    PartitionAssignment,
)


def _inputs() -> dict[str, object]:
    """Supply two target classes sharing one pooled control label."""

    proteins = ("a1", "a2", "b1", "c1", "c2")
    return {
        "comparisons": (
            ComparisonDefinition("a_vs_c", "A versus control", ("a",), ("control",), ""),
            ComparisonDefinition("b_vs_c", "B versus control", ("b",), ("control",), ""),
        ),
        "label_memberships": tuple(
            {"protein_id": protein_id, "label_id": label}
            for protein_id, label in zip(
                proteins, ("a", "a", "b", "control", "control"), strict=True
            )
        ),
        "partitions": tuple(
            PartitionAssignment(protein_id, "DISCOVERY", "EXACT_SEQUENCE", protein_id)
            for protein_id in proteins
        ),
        "control_matching_audit": (
            _match(target="a1", control="c1", status="MATCHED"),
            _match(target="a2", control="", status="UNMATCHED"),
            _match(target="b1", control="c2", status="MATCHED"),
        ),
    }


def _match(*, target: str, control: str, status: str) -> dict[str, str]:
    """Build an audit row with stable synthetic partition identities."""

    return {
        "background_label_id": "control",
        "target_label_ids": "a|b",
        "target_protein_id": target,
        "target_unit_id": target,
        "control_protein_id": control,
        "control_unit_id": control,
        "status": status,
    }


def test_audited_comparisons_exclude_unmatched_and_unrelated_controls() -> None:
    """The label remains positive, but cannot enter its unmatched comparison."""

    inputs = _inputs()
    cohorts = derive_matched_comparison_cohorts(**inputs)
    assert cohorts["a_vs_c"].target_protein_ids == frozenset({"a1"})
    assert cohorts["a_vs_c"].background_protein_ids == frozenset({"c1"})
    assert cohorts["a_vs_c"].unmatched_target_unit_count == 1
    assert cohorts["b_vs_c"].target_protein_ids == frozenset({"b1"})
    assert cohorts["b_vs_c"].background_protein_ids == frozenset({"c2"})
    settings = AnalysisSettings((2,), 1, 1, 1, 20, 0.2, 0.05, 7, 0.5, 0.5)
    associations = analyse_feature_associations(
        comparisons=inputs["comparisons"],
        label_memberships=inputs["label_memberships"],
        partitions=inputs["partitions"],
        features=(
            FeatureRecord("a1", "TYPE", "signal", "Signal", None, None, "DERIVED", "test", ""),
            FeatureRecord("b1", "TYPE", "signal", "Signal", None, None, "DERIVED", "test", ""),
        ),
        settings=settings,
        matched_cohorts=cohorts,
    )
    discovery = {row.comparison_id: row for row in associations if row.partition == "DISCOVERY"}
    assert discovery["a_vs_c"].target_protein_count == 1
    assert discovery["a_vs_c"].background_protein_count == 1
    assert discovery["b_vs_c"].target_protein_count == 1
    assert discovery["b_vs_c"].background_protein_count == 1


def test_audit_rejects_reused_control_and_partition_disagreement() -> None:
    """A stale or internally inconsistent matching authority cannot enter inference."""

    inputs = _inputs()
    reused = (
        *inputs["control_matching_audit"],
        _match(target="a2", control="c1", status="MATCHED"),
    )
    with pytest.raises(InputValidationError, match="reused"):
        derive_matched_comparison_cohorts(**{**inputs, "control_matching_audit": reused})
    stale = dict(inputs["control_matching_audit"][0], target_unit_id="other")
    with pytest.raises(InputValidationError, match="partitions"):
        derive_matched_comparison_cohorts(
            **{**inputs, "control_matching_audit": (stale, *inputs["control_matching_audit"][1:])}
        )


def test_custom_comparison_cannot_mix_audited_and_other_backgrounds() -> None:
    """An explicit comparison cannot quietly mix matched and unmatched controls."""

    inputs = _inputs()
    mixed = ComparisonDefinition(
        "mixed", "Mixed controls", ("a",), ("control", "other_control"), ""
    )
    with pytest.raises(InputValidationError, match="mixes audited and unaudited"):
        derive_matched_comparison_cohorts(**{**inputs, "comparisons": (mixed,)})
    other = ComparisonDefinition("other", "Other", ("a",), ("other_control",), "")
    assert derive_matched_comparison_cohorts(**{**inputs, "comparisons": (other,)}) == {}
