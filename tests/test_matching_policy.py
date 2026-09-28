"""Focused coverage and stability tests for outcome-blind control matching."""

from __future__ import annotations

import importlib.util
import sys
import types
import unittest
from dataclasses import replace

# The minimal local test runtime does not include the optional DuckDB adapter.
# Matching does not use it; CI with project dependencies imports the real module.
if importlib.util.find_spec("duckdb") is None:
    sys.modules.setdefault("duckdb", types.ModuleType("duckdb"))

from protein_signatures.errors import InputValidationError  # noqa: E402
from protein_signatures.evidence_labels import (  # noqa: E402
    EvidenceSettings,
    ProteinEvidenceContext,
    _match_background_units,
    _maximum_coverage_assignment,
    _summarise_control_coverage,
)


def _context(*, protein_id: str, unit_id: str) -> ProteinEvidenceContext:
    """Make a clean, equally eligible representative for a matching test.

    Args:
        protein_id: Unique sequence identifier.
        unit_id: Independent homology block identifier.

    Returns:
        Minimal covariate-complete protein context.
    """

    return ProteinEvidenceContext(
        protein_id=protein_id,
        sequence_length=100,
        sequence_sha256="0" * 64,
        species=("plant",),
        domain_ids=(),
        domain_names={},
        domain_count=0,
        pfam_assessed=True,
        structure_eligible=True,
        mean_confidence=80.0,
        input_candidate=False,
        text_evidence=(),
        direct_labels=(),
        independence_unit=unit_id,
    )


def _settings(*, strategy: str) -> EvidenceSettings:
    """Make a validated synthetic evidence policy with three control slots.

    Args:
        strategy: Versioned matching allocation name.

    Returns:
        Policy with exact species/structure and simple numerical calipers.
    """

    return EvidenceSettings(
        annotation_score=2.0,
        specific_annotation_score=4.0,
        pfam_score=2.0,
        supporting_pfam_score=0.5,
        orthology_score=2.0,
        trusted_label_score=5.0,
        minimum_acceptance_score=3.5,
        minimum_evidence_groups=2,
        minimum_score_margin=1.0,
        require_target_structure_eligible=False,
        orthology_propagation_enabled=False,
        minimum_orthology_anchor_proteins=2,
        control_units_per_target_unit=3,
        minimum_control_units_per_background=3,
        require_species_match=True,
        require_structure_match=True,
        maximum_log2_length_difference=0.75,
        maximum_domain_count_difference=2,
        exclude_input_candidates_from_controls=True,
        matching_strategy=strategy,
    )


class MatchingPolicyTests(unittest.TestCase):
    """Ensure allocation improves target coverage without changing calipers."""

    def test_coverage_first_reaches_all_targets_before_extra_slots(self) -> None:
        """Three scarce controls should cover three targets, rather than one."""

        targets = {
            key: _context(protein_id=key, unit_id=f"target_{key}") for key in ("A", "B", "C")
        }
        controls = tuple(
            _context(protein_id=f"control_{index}", unit_id=f"control_{index}")
            for index in range(3)
        )
        inputs = {
            "background_label": "control:reference",
            "target_label_ids": ("protein:class",),
            "target_ids": frozenset(targets),
            "contexts": targets,
            "candidates": controls,
        }
        legacy, _ = _match_background_units(**inputs, settings=_settings(strategy="TARGET_GREEDY"))
        revised, selected = _match_background_units(
            **inputs, settings=_settings(strategy="COVERAGE_FIRST_BOUNDED")
        )
        legacy_coverage, _ = _summarise_control_coverage(rows=legacy, requested_control_count=3)
        coverage, summary = _summarise_control_coverage(rows=revised, requested_control_count=3)
        self.assertEqual(sorted(row["matched_control_count"] for row in legacy_coverage), [0, 0, 3])
        self.assertEqual(sorted(row["matched_control_count"] for row in coverage), [1, 1, 1])
        self.assertEqual(summary[0]["target_coverage_fraction"], 1.0)
        self.assertEqual(len(selected), 3)
        self.assertEqual(
            len({row["control_unit_id"] for row in revised if row["status"] == "MATCHED"}), 3
        )
        reversed_rows, _ = _match_background_units(
            **{
                **inputs,
                "candidates": tuple(reversed(controls)),
                "contexts": dict(reversed(tuple(targets.items()))),
            },
            settings=_settings(strategy="COVERAGE_FIRST_BOUNDED"),
        )
        self.assertEqual(revised, reversed_rows)

    def test_incompatible_controls_remain_unmatched_and_visible(self) -> None:
        """Coverage-first cannot fabricate compatible controls."""

        target = _context(protein_id="target", unit_id="target_unit")
        wrong_species = replace(
            _context(protein_id="control", unit_id="control_unit"), species=("animal",)
        )
        rows, selected = _match_background_units(
            background_label="control:reference",
            target_label_ids=("protein:class",),
            target_ids=frozenset({"target"}),
            contexts={"target": target},
            candidates=(wrong_species,),
            settings=_settings(strategy="COVERAGE_FIRST_BOUNDED"),
        )
        coverage, summary = _summarise_control_coverage(rows=rows, requested_control_count=3)
        self.assertEqual(selected, frozenset())
        self.assertEqual(coverage[0]["coverage_status"], "UNMATCHED")
        self.assertEqual(summary[0]["target_without_control_count"], 1)

    def test_coverage_summary_rejects_control_reuse(self) -> None:
        """The exported audit cannot silently count one control twice."""

        rows = (
            {
                "background_label_id": "control:a",
                "target_unit_id": "t1",
                "target_protein_id": "p1",
                "control_unit_id": "c1",
                "status": "MATCHED",
            },
            {
                "background_label_id": "control:a",
                "target_unit_id": "t2",
                "target_protein_id": "p2",
                "control_unit_id": "c1",
                "status": "MATCHED",
            },
        )
        with self.assertRaisesRegex(InputValidationError, "reused"):
            _summarise_control_coverage(rows=rows, requested_control_count=1)

    def test_augmenting_path_recovers_from_an_early_choice(self) -> None:
        """Candidate reallocation should cover a target that greedy selection misses."""

        controls = {name: _context(protein_id=name, unit_id=name) for name in ("A", "B", "C", "D")}
        choices = {
            "t1": (controls["A"], controls["B"]),
            "t2": (controls["A"], controls["C"]),
            "t3": (controls["A"], controls["C"]),
            "t4": (controls["C"], controls["D"]),
        }
        assigned = _maximum_coverage_assignment(preferences=choices)
        self.assertEqual(set(assigned), set(choices))
        self.assertEqual({value.independence_unit for value in assigned.values()}, set(controls))


if __name__ == "__main__":
    unittest.main()
