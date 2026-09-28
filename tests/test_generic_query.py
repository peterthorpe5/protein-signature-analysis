"""Focused tests for applying reference signatures to new protein collections."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.io_utils import iter_tsv, write_tsv_atomic
from protein_signatures.query import (
    HIT_FIELDS,
    STRUCTURAL_PROJECTION_FIELDS,
    SUMMARY_FIELDS,
    project_structural_signatures,
    query_signature_evidence,
)


class GenericQueryTests(unittest.TestCase):
    """Verify feature state semantics for proteins and overlapping groups."""

    def setUp(self) -> None:
        """Create a small, independent reference and query dataset."""

        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fasta = self.root / "query.faa"
        self.fasta.write_text(">p1\nMABCDE\n>p2\nMABCDF\n>p3\nMZZZZZ\n", encoding="utf-8")
        self.comparisons = self.root / "comparisons.tsv"
        write_tsv_atomic(
            path=self.comparisons,
            fieldnames=("comparison_id", "target_label_ids", "background_label_ids"),
            records=(
                {
                    "comparison_id": "enzyme_a",
                    "target_label_ids": "enzyme:a",
                    "background_label_ids": "control:a",
                },
                {
                    "comparison_id": "enzyme_b",
                    "target_label_ids": "enzyme:b",
                    "background_label_ids": "control:b",
                },
            ),
        )
        self.signatures = self.root / "signatures.tsv"
        write_tsv_atomic(
            path=self.signatures,
            fieldnames=(
                "comparison_id",
                "feature_type",
                "feature_id",
                "feature_name",
                "validation_q_value",
                "validation_study_q_value",
                "validation_prevalence_difference",
                "evidence_class",
                "status",
            ),
            records=(
                {
                    "comparison_id": "enzyme_a",
                    "feature_type": "STRUCTURE_CLUSTER",
                    "feature_id": "SC_1",
                    "feature_name": "Fold one",
                    "validation_q_value": 0.01,
                    "validation_study_q_value": 0.02,
                    "validation_prevalence_difference": 0.4,
                    "evidence_class": "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE",
                    "status": "COMPLETE",
                },
                {
                    "comparison_id": "enzyme_a",
                    "feature_type": "AMINO_ACID_KMER",
                    "feature_id": "k3:ABC",
                    "feature_name": "ABC",
                    "validation_q_value": 0.04,
                    "validation_study_q_value": 0.12,
                    "validation_prevalence_difference": 0.2,
                    "evidence_class": "DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON",
                    "status": "COMPLETE",
                },
                {
                    "comparison_id": "enzyme_a",
                    "feature_type": "DOMAIN",
                    "feature_id": "PF_test",
                    "feature_name": "Test domain",
                    "validation_q_value": 0.04,
                    "validation_study_q_value": 0.12,
                    "validation_prevalence_difference": 0.2,
                    "evidence_class": "DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON",
                    "status": "COMPLETE",
                },
                {
                    "comparison_id": "enzyme_a",
                    "feature_type": "STRUCTURE_CLUSTER",
                    "feature_id": "SC_reverse",
                    "feature_name": "Reverse",
                    "validation_q_value": 0.01,
                    "validation_study_q_value": 0.02,
                    "validation_prevalence_difference": -0.4,
                    "evidence_class": "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE",
                    "status": "COMPLETE",
                },
            ),
        )
        self.features = self.root / "query_features.tsv"
        write_tsv_atomic(
            path=self.features,
            fieldnames=("protein_id", "feature_type", "feature_id", "evidence_status"),
            records=(
                {
                    "protein_id": "p1",
                    "feature_type": "STRUCTURE_CLUSTER",
                    "feature_id": "SC_1",
                    "evidence_status": "ASSESSED_WITH_FEATURE",
                },
                {
                    "protein_id": "p2",
                    "feature_type": "STRUCTURE_CLUSTER",
                    "feature_id": "SC_1",
                    "evidence_status": "ASSESSED_NO_FEATURE",
                },
                {
                    "protein_id": "p3",
                    "feature_type": "STRUCTURE_CLUSTER",
                    "feature_id": "SC_1",
                    "evidence_status": "ASSESSED_NO_FEATURE",
                },
                {
                    "protein_id": "p1",
                    "feature_type": "DOMAIN",
                    "feature_id": "PF_test",
                    "evidence_status": "ASSESSED_NO_FEATURE",
                },
                {
                    "protein_id": "p2",
                    "feature_type": "DOMAIN",
                    "feature_id": "PF_test",
                    "evidence_status": "ASSESSED_NO_FEATURE",
                },
            ),
        )

    def test_protein_queries_preserve_unknown_and_frozen_scope(self) -> None:
        """Absent, missing and unvalidated evidence must remain distinct."""

        output = query_signature_evidence(
            signatures_path=self.signatures,
            comparisons_path=self.comparisons,
            query_fasta=self.fasta,
            query_features_path=self.features,
            output_dir=self.root / "protein_query",
        )
        rows = {
            (row["unit_id"], row["comparison_id"]): row
            for row in iter_tsv(
                path=output / "candidate_summary.tsv", required_fields=SUMMARY_FIELDS
            )
        }
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[("p1", "enzyme_a")]["evidence_state"], "FEATURE_EVIDENCE_PRESENT")
        self.assertEqual(rows[("p2", "enzyme_a")]["evidence_state"], "ASSESSED_NO_FEATURE")
        self.assertEqual(
            rows[("p3", "enzyme_b")]["evidence_state"], "NO_VALIDATED_REFERENCE_FEATURES"
        )
        hits = tuple(
            iter_tsv(path=output / "candidate_feature_hits.tsv", required_fields=HIT_FIELDS)
        )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["feature_id"], "SC_1")
        self.assertEqual(hits[0]["validation_fdr_scope"], "STUDY_WIDE")
        marker = json.loads((output / "QUERY_COMPLETED.json").read_text(encoding="utf-8"))
        self.assertEqual(marker["unit_count"], 3)
        self.assertEqual(
            marker["interpretation"], "UNREVIEWED_CANDIDATE_EVIDENCE_NOT_CLASS_PROBABILITY"
        )
        with self.assertRaises(PublicationError):
            query_signature_evidence(
                signatures_path=self.signatures,
                comparisons_path=self.comparisons,
                query_fasta=self.fasta,
                query_features_path=self.features,
                output_dir=output,
            )

    def test_groups_use_any_positive_and_require_all_negative_for_absence(self) -> None:
        """Overlapping groups must aggregate assessed status conservatively."""

        units = self.root / "units.tsv"
        write_tsv_atomic(
            path=units,
            fieldnames=("protein_id", "unit_id"),
            records=(
                {"protein_id": "p1", "unit_id": "group_1"},
                {"protein_id": "p2", "unit_id": "group_1"},
                {"protein_id": "p2", "unit_id": "group_2"},
                {"protein_id": "p3", "unit_id": "group_2"},
            ),
        )
        output = query_signature_evidence(
            signatures_path=self.signatures,
            comparisons_path=self.comparisons,
            query_fasta=self.fasta,
            query_features_path=self.features,
            query_units_path=units,
            output_dir=self.root / "group_query",
            include_within_comparison=True,
        )
        rows = {
            (row["unit_id"], row["comparison_id"]): row
            for row in iter_tsv(
                path=output / "candidate_summary.tsv", required_fields=SUMMARY_FIELDS
            )
        }
        self.assertEqual(rows[("group_1", "enzyme_a")]["present_feature_count"], "2")
        self.assertEqual(rows[("group_1", "enzyme_a")]["assessed_absent_feature_count"], "1")
        self.assertEqual(rows[("group_2", "enzyme_a")]["unknown_feature_count"], "1")
        self.assertEqual(
            rows[("group_2", "enzyme_a")]["evidence_state"], "FEATURE_EVIDENCE_PRESENT"
        )
        self.assertEqual(rows[("group_1", "enzyme_a")]["unit_type"], "GROUP")
        hits = tuple(
            iter_tsv(path=output / "candidate_feature_hits.tsv", required_fields=HIT_FIELDS)
        )
        structural_hit = next(row for row in hits if row["feature_id"] == "SC_1")
        self.assertEqual(structural_hit["positive_protein_ids"], "p1")

    def test_exact_kmers_are_assessed_from_fasta_without_external_tables(self) -> None:
        """A frozen exact sequence feature is present or absent for every protein."""

        output = query_signature_evidence(
            signatures_path=self.signatures,
            comparisons_path=self.comparisons,
            query_fasta=self.fasta,
            query_features_path=None,
            output_dir=self.root / "native_query",
            include_within_comparison=True,
        )
        rows = {
            (row["unit_id"], row["comparison_id"]): row
            for row in iter_tsv(
                path=output / "candidate_summary.tsv", required_fields=SUMMARY_FIELDS
            )
        }
        self.assertEqual(rows[("p1", "enzyme_a")]["present_feature_count"], "1")
        self.assertEqual(rows[("p2", "enzyme_a")]["present_feature_count"], "1")
        self.assertEqual(rows[("p3", "enzyme_a")]["assessed_absent_feature_count"], "1")
        self.assertEqual(rows[("p3", "enzyme_a")]["unknown_feature_count"], "2")

    def test_supplied_kmer_cannot_disagree_with_sequence(self) -> None:
        """A false negative for an exact query k-mer must fail closed."""

        with self.features.open(mode="a", encoding="utf-8") as handle:
            handle.write("p1\tAMINO_ACID_KMER\tk3:ABC\tASSESSED_NO_FEATURE\n")
        with self.assertRaisesRegex(InputValidationError, "disagrees with query FASTA"):
            query_signature_evidence(
                signatures_path=self.signatures,
                comparisons_path=self.comparisons,
                query_fasta=self.fasta,
                query_features_path=self.features,
                output_dir=self.root / "invalid_kmer_query",
                include_within_comparison=True,
            )

    def test_multiple_feature_tables_combine_without_losing_unknowns(self) -> None:
        """Separate feature producers should contribute to one frozen query."""

        extra = self.root / "additional_features.tsv"
        write_tsv_atomic(
            path=extra,
            fieldnames=("protein_id", "feature_type", "feature_id", "evidence_status"),
            records=(
                {
                    "protein_id": "p3",
                    "feature_type": "DOMAIN",
                    "feature_id": "PF_test",
                    "evidence_status": "ASSESSED_WITH_FEATURE",
                },
            ),
        )
        output = query_signature_evidence(
            signatures_path=self.signatures,
            comparisons_path=self.comparisons,
            query_fasta=self.fasta,
            query_features_path=self.features,
            additional_query_features_paths=(extra,),
            output_dir=self.root / "combined_query",
            include_within_comparison=True,
        )
        rows = {
            (row["unit_id"], row["comparison_id"]): row
            for row in iter_tsv(
                path=output / "candidate_summary.tsv", required_fields=SUMMARY_FIELDS
            )
        }
        self.assertEqual(rows[("p3", "enzyme_a")]["present_feature_count"], "1")
        self.assertEqual(rows[("p3", "enzyme_a")]["unknown_feature_count"], "0")
        marker = json.loads((output / "QUERY_COMPLETED.json").read_text(encoding="utf-8"))
        self.assertEqual(len(marker["inputs"]["query_feature_tables"]), 2)

    def test_one_frozen_comparison_can_be_queried_by_id(self) -> None:
        """A requested class should not produce all other reference comparisons."""

        output = query_signature_evidence(
            signatures_path=self.signatures,
            comparisons_path=self.comparisons,
            query_fasta=self.fasta,
            query_features_path=self.features,
            output_dir=self.root / "selected_query",
            comparison_ids=("enzyme_a",),
        )
        rows = tuple(
            iter_tsv(path=output / "candidate_summary.tsv", required_fields=SUMMARY_FIELDS)
        )
        self.assertEqual(len(rows), 3)
        self.assertEqual({row["comparison_id"] for row in rows}, {"enzyme_a"})
        marker = json.loads((output / "QUERY_COMPLETED.json").read_text(encoding="utf-8"))
        self.assertEqual(marker["selected_comparison_ids"], ["enzyme_a"])
        with self.assertRaisesRegex(InputValidationError, "Unknown query comparison"):
            query_signature_evidence(
                signatures_path=self.signatures,
                comparisons_path=self.comparisons,
                query_fasta=self.fasta,
                query_features_path=self.features,
                output_dir=self.root / "unknown_class",
                comparison_ids=("not_in_reference",),
            )

    def test_conflicting_feature_states_fail(self) -> None:
        """An ambiguous assessment must stop publication."""

        with self.features.open(mode="a", encoding="utf-8") as handle:
            handle.write("p1\tSTRUCTURE_CLUSTER\tSC_1\tASSESSED_NO_FEATURE\n")
        with self.assertRaisesRegex(InputValidationError, "Conflicting"):
            query_signature_evidence(
                signatures_path=self.signatures,
                comparisons_path=self.comparisons,
                query_fasta=self.fasta,
                query_features_path=self.features,
                output_dir=self.root / "invalid_query",
            )
        self.assertFalse((self.root / "invalid_query").exists())

    def test_group_membership_cannot_drop_fasta_proteins(self) -> None:
        """A partial orthogroup mapping must be rejected rather than silently filtered."""

        units = self.root / "partial_units.tsv"
        write_tsv_atomic(
            path=units,
            fieldnames=("protein_id", "unit_id"),
            records=({"protein_id": "p1", "unit_id": "group_1"},),
        )
        with self.assertRaisesRegex(InputValidationError, "omits FASTA proteins"):
            query_signature_evidence(
                signatures_path=self.signatures,
                comparisons_path=self.comparisons,
                query_fasta=self.fasta,
                query_features_path=self.features,
                query_units_path=units,
                output_dir=self.root / "partial_query",
            )


class StructuralProjectionTests(unittest.TestCase):
    """Check leakage-safe positive projection onto original discovery members."""

    def setUp(self) -> None:
        """Create a frozen reference and a small new query collection."""

        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fasta = self.root / "query.faa"
        self.fasta.write_text(">q1\nMABCDE\n>q2\nMABCDF\n", encoding="utf-8")
        self.clusters = self.root / "structure_clusters.tsv"
        fields = (
            "cluster_id",
            "protein_id",
            "reference_partition",
            "membership_method",
            "tm_score_threshold",
            "minimum_coverage",
            "coverage_scope",
            "comparison_tool",
            "comparison_tool_version",
        )
        write_tsv_atomic(
            path=self.clusters,
            fieldnames=fields,
            records=(
                {
                    "cluster_id": "SC_frozen",
                    "protein_id": protein_id,
                    "reference_partition": "DISCOVERY",
                    "membership_method": method,
                    "tm_score_threshold": 0.5,
                    "minimum_coverage": 0.5,
                    "coverage_scope": "STRUCTURE_MODEL_RESIDUES",
                    "comparison_tool": "Foldseek",
                    "comparison_tool_version": "test_1",
                }
                for protein_id, method in (
                    ("r1", "DISCOVERY_COMPONENT"),
                    ("r2", "DISCOVERY_COMPONENT"),
                    ("v1", "VALIDATION_PROJECTION"),
                )
            ),
        )
        self.comparisons = self.root / "query_comparisons.tsv"
        self.fields = (
            "protein_a_id",
            "protein_b_id",
            "comparison_tool",
            "comparison_tool_version",
            "tm_score",
            "coverage_a",
            "coverage_b",
            "comparison_status",
            "comparison_universe_id",
            "coverage_scope",
        )

    def test_projection_uses_only_frozen_members_and_never_infers_absence(self) -> None:
        """Passing edges from either orientation support one positive cluster call."""

        records = (
            ("q1", "r1", 0.8, 0.7, 0.6),
            ("r2", "q1", 0.9, 0.6, 0.8),
            ("q2", "r1", 0.9, 0.4, 0.9),
            ("q2", "v1", 0.9, 0.9, 0.9),
        )
        write_tsv_atomic(
            path=self.comparisons,
            fieldnames=self.fields,
            records=(
                {
                    "protein_a_id": a,
                    "protein_b_id": b,
                    "comparison_tool": "Foldseek",
                    "comparison_tool_version": "test_1",
                    "tm_score": score,
                    "coverage_a": coverage_a,
                    "coverage_b": coverage_b,
                    "comparison_status": "COMPLETE",
                    "comparison_universe_id": "new_search",
                    "coverage_scope": "STRUCTURE_MODEL_RESIDUES",
                }
                for a, b, score, coverage_a, coverage_b in records
            ),
        )
        output = project_structural_signatures(
            reference_clusters_path=self.clusters,
            query_comparisons_path=self.comparisons,
            query_fasta=self.fasta,
            output_dir=self.root / "projection",
        )
        rows = tuple(
            iter_tsv(
                path=output / "projected_structural_features.tsv",
                required_fields=STRUCTURAL_PROJECTION_FIELDS,
            )
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["protein_id"], "q1")
        self.assertEqual(rows[0]["supporting_reference_count"], "2")
        self.assertEqual(float(rows[0]["best_tm_score"]), 0.9)
        self.assertEqual(rows[0]["evidence_status"], "ASSESSED_WITH_FEATURE")
        marker = json.loads((output / "PROJECTION_COMPLETED.json").read_text(encoding="utf-8"))
        self.assertEqual(marker["interpretation"], "NO_HIT_IS_UNKNOWN_NOT_ASSESSED_NEGATIVE")

    def test_mismatched_tool_version_fails_before_publication(self) -> None:
        """Structural definitions cannot silently cross tool versions."""

        write_tsv_atomic(
            path=self.comparisons,
            fieldnames=self.fields,
            records=(
                {
                    "protein_a_id": "q1",
                    "protein_b_id": "r1",
                    "comparison_tool": "Foldseek",
                    "comparison_tool_version": "other",
                    "tm_score": 0.8,
                    "coverage_a": 0.8,
                    "coverage_b": 0.8,
                    "comparison_status": "COMPLETE",
                    "comparison_universe_id": "new_search",
                    "coverage_scope": "STRUCTURE_MODEL_RESIDUES",
                },
            ),
        )
        with self.assertRaisesRegex(InputValidationError, "incompatible"):
            project_structural_signatures(
                reference_clusters_path=self.clusters,
                query_comparisons_path=self.comparisons,
                query_fasta=self.fasta,
                output_dir=self.root / "invalid_projection",
            )


if __name__ == "__main__":
    unittest.main()
