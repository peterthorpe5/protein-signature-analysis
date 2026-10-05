"""Scientific dictionary and campaign-specific interpretation regression tests."""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
import zipfile
from io import BytesIO

import numpy as np
import pandas as pd
import pytest

from protein_signature_app.help_content import (
    GRAPH_HELP,
    LIMIT_SPECS,
    METRIC_HELP,
    PAGE_METHODS,
    PAGE_TERMS,
    STATUS_DEFINITIONS,
    campaign_limit_rows,
    feature_explanation,
    glossary_rows,
    graph_explanation,
    metric_reading,
    recorded_setting,
    status_explanation,
)
from protein_signature_app.viewer_help import GLOSSARY, PAGE_HELP
from protein_signatures.errors import InputValidationError
from protein_signatures.exports import dataframe_to_xlsx_bytes
from protein_signatures.feature_provenance import (
    FEATURE_ASSESSMENT_STATUSES,
    FEATURE_DERIVATION_SCOPES,
)
from protein_signatures.models import (
    AnalysisStatus,
    CurationStatus,
    DomainAssessmentStatus,
    FoldEvidenceStatus,
    StructureAnalysisEligibility,
    StructureComparisonStatus,
    StructureCoverageScope,
)
from protein_signatures.result_help import COLUMN_DEFINITIONS, column_definition
from protein_signatures.schemas import table_schemas


def test_every_canonical_field_and_page_term_has_an_explicit_glossary_entry() -> None:
    """New fields or pages must fail tests until scientifically defined."""
    rows = glossary_rows(base_rows=GLOSSARY)
    definitions = {(category, term): definition for category, term, definition in rows}
    terms = {term for _, term, _ in rows}
    fields = {name for schema in table_schemas().values() for name in schema.names}
    missing = fields - {term for category, term in definitions if category == "Table fields"}
    assert not missing, f"Missing canonical definitions: {sorted(missing)}"
    assert set(PAGE_METHODS) == set(PAGE_TERMS) == set(PAGE_HELP)
    assert len(PAGE_HELP) == 11
    assert all(set(page_terms).issubset(terms) for page_terms in PAGE_TERMS.values())
    assert all(column_definition(column_name=name) for name in fields)
    assert len(rows) == len(definitions)
    assert rows == tuple(sorted(rows))
    assert (
        "Exports",
        "Completed result",
        dict((term, text) for _, term, text in GLOSSARY)["Completed result"],
    ) in rows
    for section, key, _, _ in LIMIT_SPECS:
        assert ("Campaign settings", f"{section}.{key}") in definitions
    for prefix in ("discovery_", "validation_"):
        assert ("Table fields", f"{prefix}background_prevalence") in definitions


@pytest.mark.parametrize(
    "enum_type",
    [
        AnalysisStatus,
        CurationStatus,
        DomainAssessmentStatus,
        FoldEvidenceStatus,
        StructureAnalysisEligibility,
        StructureComparisonStatus,
        StructureCoverageScope,
    ],
)
def test_controlled_states_are_defined(enum_type: type) -> None:
    """Every canonical controlled state must be searchable in the glossary."""
    assert {state.value for state in enum_type}.issubset(STATUS_DEFINITIONS)


def test_feature_states_and_model_summary_semantics_are_explicit() -> None:
    """Coefficient ranks, penalty strengths and permutation scores must not be conflated."""
    assert FEATURE_ASSESSMENT_STATUSES.issubset(STATUS_DEFINITIONS)
    assert FEATURE_DERIVATION_SCOPES.issubset(STATUS_DEFINITIONS)
    assert "C = 1/lambda" in column_definition(column_name="regularisation_strength")
    assert "absolute fitted coefficient" in column_definition(column_name="importance_rank")
    assert "balanced-accuracy drop" in column_definition(
        column_name="validation_permutation_importance_mean"
    )
    assert "validation proteins" in column_definition(
        column_name="mean_absolute_validation_contribution"
    )
    assert "bounded preview" in column_definition(column_name="row_count")
    assert "ascending discovery q-value" in column_definition(column_name="within_comparison_rank")


def test_combined_viewer_summaries_define_validation_and_matching_scope() -> None:
    """Merged comparison summaries must preserve validation tiers and actual cohort scope."""
    fields = {
        "observed_proteins",
        "feature_label",
        "has_validation",
        "validated_within_count",
        "validated_study_count",
        "matched_target_units",
        "excluded_unmatched_target_units",
        "analysed_target_proteins",
        "excluded_unmatched_target_proteins",
        "matched_control_proteins",
        "minimum_coverage_bin",
        "tm_score_bin",
        "comparison_count",
    }
    assert fields.issubset(COLUMN_DEFINITIONS)
    assert "decision-candidate" in column_definition(column_name="enriched_count")
    assert "comparison or pooled audit" in column_definition(column_name="matched_control_units")
    assert "matched and unmatched" in column_definition(column_name="target_units")
    assert "provisional" in column_definition(column_name="observed_proteins")
    assert "held-out" in STATUS_DEFINITIONS["DISCOVERY_ENRICHED"]


def test_column_definition_preserves_scope_and_logs_unknown_fields(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Unknown producers must not acquire invented definitions."""
    assert column_definition(column_name=" protein_id ") == COLUMN_DEFINITIONS["protein_id"]
    assert column_definition(column_name="discovery_sequence_length").startswith(
        "Discovery partition:"
    )
    assert column_definition(column_name="validation_sequence_length").startswith(
        "Held-out validation:"
    )
    with caplog.at_level(level=logging.DEBUG, logger="protein_signatures.result_help"):
        assert "not registered" in column_definition(column_name="unknown_producer_field")
    assert "unknown_producer_field" in caplog.text


@pytest.mark.parametrize("name", [None, 1, "", "   "])
def test_column_help_rejects_invalid_names(name: object) -> None:
    """Malformed headings should be rejected rather than silently mislabelled."""
    with pytest.raises(InputValidationError, match="non-empty text"):
        column_definition(column_name=name)


def test_excel_definitions_match_tooltips_and_glossary() -> None:
    """Generated Excel dictionaries must share the authoritative field meanings."""
    fields = ("validation_roc_auc", "target_assessed_unit_count", "domain_architecture_jaccard")
    payload = dataframe_to_xlsx_bytes(frame=pd.DataFrame(columns=fields), title="Evidence")
    with zipfile.ZipFile(file=BytesIO(payload)) as archive:
        xml = ET.fromstring(text=archive.read(name="xl/sharedStrings.xml"))
    text = "".join(xml.itertext())
    for name in fields:
        assert column_definition(column_name=name) in text


@pytest.mark.parametrize("token", GRAPH_HELP)
def test_each_graph_has_specific_interpretation(token: str) -> None:
    """Figure identities must resolve to their own axes and evidence limits."""
    assert graph_explanation(graph_name=f"campaign_{token}") == GRAPH_HELP[token]


def test_graph_aliases_and_unknown_graphs_have_safe_help() -> None:
    """Model aliases and future plots must retain defensible interpretation."""
    assert graph_explanation(graph_name="three_dimensional_model") == GRAPH_HELP["trace"]
    assert graph_explanation(graph_name="top_enriched_features") == GRAPH_HELP["ranked_association"]
    assert "labelled axes" in graph_explanation(graph_name="future_producer_plot")


@pytest.mark.parametrize("name", [None, 1, "", "  "])
def test_graph_help_rejects_invalid_identity(name: object) -> None:
    """A missing plot identity must produce a useful validation error."""
    with pytest.raises(InputValidationError, match="non-empty text"):
        graph_explanation(graph_name=name)


@pytest.mark.parametrize(
    "feature_type, feature_id, expected",
    [
        ("AMINO_ACID_KMER", "k3:LPD", "3-residue amino-acid word LPD"),
        ("AMINO_ACID_KMER", "k10:AAAAAAAAAA", "10-residue"),
        ("AMINO_ACID_KMER", "k3:LP", "Unrecognised"),
        ("AMINO_ACID_KMER", "LPD", "Unrecognised"),
        ("STRUCTURE_CLUSTER", "SC_123", "no named fold or local residue interval"),
        ("PFAM_ARCHITECTURE", "PF1|PF2", "Ordered detected domain"),
        ("DOMAIN_ARCHITECTURE", "PF1|PF2", "not a contiguous"),
        ("PFAM_DOMAIN", "PF00001", "scan authority"),
        ("DOMAIN", "PF00001", "residue intervals"),
        ("FOLD", "fold1", "explicit authority"),
        ("STRUCTURAL_POCKET", "pocket1", "no enrichment-tested residue interval"),
        ("STRUCTURE_AVAILABILITY", "available", "Technical model"),
        ("STRUCTURE_QUALITY", "high", "not a biochemical"),
        ("STRUCTURE_AVAILABLE", "available", "Technical model"),
        ("EVIDENCE_AVAILABILITY", "ASSESSED", "not a biochemical"),
        ("FUTURE_PRODUCER", "id1", "declared producer"),
    ],
)
def test_feature_keys_have_readable_bounded_claims(
    feature_type: str,
    feature_id: str,
    expected: str,
) -> None:
    """Exact sequence words must not be described as demonstrated functional motifs."""
    assert expected in feature_explanation(feature_type=feature_type, feature_id=feature_id)


@pytest.mark.parametrize("feature_type, feature_id", [(None, "id"), ("FOLD", 1)])
def test_feature_explanations_reject_non_text_keys(
    feature_type: object,
    feature_id: object,
) -> None:
    """Invalid feature identifiers must not be guessed or string-coerced."""
    with pytest.raises(InputValidationError, match="text identifiers"):
        feature_explanation(feature_type=feature_type, feature_id=feature_id)


@pytest.mark.parametrize("value", [None, pd.NA, "0.9", True, float("nan"), float("inf")])
def test_unavailable_metrics_are_not_interpreted(value: object) -> None:
    """Missing, invalid and non-finite values cannot acquire a performance grade."""
    assert "Not calculated" in metric_reading(metric_name="validation_roc_auc", value=value)


@pytest.mark.parametrize(
    "name, value, expected",
    [
        ("cv_roc_auc", 0.5, "50.0%"),
        ("validation_roc_auc", np.float64(0.875), "87.5%"),
        ("validation_matthews_correlation", -0.604, "-0.604 correlation"),
        ("validation_balanced_accuracy", 0.8, "80.0%"),
        ("validation_average_precision", 0.9, "0.250 target prevalence"),
        ("validation_brier_score", 0.04, "reference 0.188"),
    ],
)
def test_metric_values_are_read_against_the_correct_baseline(
    name: str,
    value: float,
    expected: str,
) -> None:
    """Ranking, classification and probability-error values must remain distinct."""
    assert expected in metric_reading(metric_name=name, value=value, baseline=0.25)
    assert set(METRIC_HELP) == {
        "cv_roc_auc",
        "validation_roc_auc",
        "validation_average_precision",
        "validation_matthews_correlation",
        "validation_balanced_accuracy",
        "validation_brier_score",
    }


@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_out_of_range_metric_values_are_reported(value: float) -> None:
    """Corrupt stored values must not be presented as meaningful scores."""
    assert "valid range" in metric_reading(metric_name="validation_roc_auc", value=value)
    assert "valid range" in metric_reading(
        metric_name="validation_matthews_correlation", value=-1.01
    )


@pytest.mark.parametrize("baseline", [None, True, "0.2", float("nan"), float("inf"), -0.1, 1.1])
def test_metric_reference_is_not_invented(baseline: object) -> None:
    """Invalid class fractions must produce explicit unavailable baselines."""
    assert "unrecorded" in metric_reading(
        metric_name="validation_average_precision", value=0.9, baseline=baseline
    )
    assert "not available" in metric_reading(
        metric_name="validation_brier_score", value=0.04, baseline=baseline
    )


def test_unknown_metric_identity_is_rejected() -> None:
    """An unregistered statistic must not inherit another metric's meaning."""
    with pytest.raises(InputValidationError, match="Unknown model metric"):
        metric_reading(metric_name="future_metric", value=0.5)


@pytest.mark.parametrize(
    "key, value, expected",
    [
        ("minimum_target_proteins", 7, "7"),
        ("minimum_target_proteins", 0, "Not recorded"),
        ("minimum_target_proteins", 2.5, "Not recorded"),
        ("cross_validation_folds", 4, "4"),
        ("minimum_mean_plddt", 70.5, "70.5"),
        ("minimum_mean_plddt", 101, "Not recorded"),
        ("fdr_threshold", 0.01, "0.01"),
        ("fdr_threshold", -0.1, "Not recorded"),
        ("minimum_target_proteins", True, "Not recorded"),
        ("minimum_target_proteins", "10", "Not recorded"),
        ("minimum_target_proteins", float("inf"), "Not recorded"),
    ],
)
def test_recorded_limits_are_validated_without_defaulting(
    key: str,
    value: object,
    expected: str,
) -> None:
    """Saved non-default settings and malformed metadata must be distinguished."""
    metadata = {"campaign": {"analysis": {key: value}}}
    assert recorded_setting(metadata=metadata, section="analysis", key=key) == expected


@pytest.mark.parametrize("metadata", [{}, None, {"campaign": []}, {"campaign": {"analysis": []}}])
def test_missing_campaign_metadata_does_not_use_current_defaults(metadata: object) -> None:
    """Historical resources lacking settings must display Not recorded."""
    assert (
        recorded_setting(metadata=metadata, section="analysis", key="minimum_target_proteins")
        == "Not recorded"
    )
    assert all(row["limit"] == "Not recorded" for row in campaign_limit_rows(metadata=metadata))


def test_sample_status_explains_actual_limits_and_evaluation_units() -> None:
    """Association and model minimums must report different counted units."""
    metadata = {
        "campaign": {
            "analysis": {"minimum_target_proteins": 7, "minimum_background_proteins": 9},
            "explainable_ml": {"minimum_samples_per_class": 17, "minimum_groups_per_class": 6},
        }
    }
    association = status_explanation(status="INSUFFICIENT_SAMPLE_SIZE", metadata=metadata)
    assert "7 target and 9 background independent blocks" in association
    assert "Mixed blocks and unknown absences do not count" in association
    for status in (
        "INSUFFICIENT_SAMPLE_SIZE",
        "INSUFFICIENT_INDEPENDENT_GROUPS",
        "COMPLETE_VALIDATION_UNDERPOWERED",
    ):
        model = status_explanation(status=status, metadata=metadata, context="model")
        assert "17 protein samples AND 6 pure independent groups" in model
        assert "same two per-class limits" in model
    warned = status_explanation(
        status="COMPLETE_VALIDATION_UNDERPOWERED_CONVERGENCE_WARNING",
        metadata=metadata,
        context="model",
    )
    assert "iteration limit" in warned
    assert "17 protein samples" in warned
    assert "Producer-specific" in status_explanation(status="UNSEEN_UPSTREAM_CODE", metadata={})
    assert "Not recorded target" in status_explanation(
        status="INSUFFICIENT_SAMPLE_SIZE", metadata={}
    )
    limits = {row["setting"]: row for row in campaign_limit_rows(metadata=metadata)}
    assert limits["analysis.minimum_target_proteins"]["limit"] == "7"
    assert "blocks" in limits["analysis.minimum_target_proteins"]["unit"]


@pytest.mark.parametrize(
    "status, context", [(None, "model"), ("", "model"), ("COMPLETE", "unknown")]
)
def test_status_help_rejects_invalid_inputs(status: object, context: str) -> None:
    """Invalid codes and evaluation contexts must fail explicitly."""
    with pytest.raises(InputValidationError, match="valid context"):
        status_explanation(status=status, metadata={}, context=context)
