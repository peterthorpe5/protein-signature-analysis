"""Tests for explicit feature-assessment universes and modality composition."""

from __future__ import annotations

import logging
from dataclasses import replace

import pytest

from protein_signatures.assessment import (
    compose_feature_assessment_universes,
    normalise_feature_assessment_universes,
    normalise_feature_keys,
    validate_feature_inference_scopes,
)
from protein_signatures.errors import InputValidationError
from protein_signatures.feature_provenance import (
    ALL_DATA_EXPLORATORY,
    DISCOVERY_SUPERVISED,
)
from protein_signatures.models import (
    DomainAssessment,
    DomainAssessmentStatus,
    FeatureRecord,
)


def test_assessment_universe_normalisation_is_strict() -> None:
    """Observed positives must be assessed and all identifiers must be known."""

    positives = {("TYPE", "x"): {"p1"}}
    assert (
        normalise_feature_assessment_universes(
            universes=None,
            known_protein_ids={"p1", "p2"},
            feature_proteins=positives,
        )
        is None
    )
    normalised = normalise_feature_assessment_universes(
        universes={("TYPE", "x"): {"p1", "p2"}},
        known_protein_ids={"p1", "p2"},
        feature_proteins=positives,
    )
    assert normalised == {("TYPE", "x"): frozenset({"p1", "p2"})}
    with pytest.raises(InputValidationError, match="Missing assessment universe"):
        normalise_feature_assessment_universes(
            universes={},
            known_protein_ids={"p1"},
            feature_proteins=positives,
        )
    with pytest.raises(InputValidationError, match="not declared assessed"):
        normalise_feature_assessment_universes(
            universes={("TYPE", "x"): set()},
            known_protein_ids={"p1"},
            feature_proteins=positives,
        )
    with pytest.raises(InputValidationError, match="unknown proteins"):
        normalise_feature_assessment_universes(
            universes={("TYPE", "x"): {"p1", "missing"}},
            known_protein_ids={"p1"},
            feature_proteins=positives,
        )
    with pytest.raises(InputValidationError, match="outside the known universe"):
        normalise_feature_assessment_universes(
            universes=None,
            known_protein_ids={"p2"},
            feature_proteins=positives,
        )


def test_feature_key_and_identifier_collections_are_validated() -> None:
    """Malformed scalar collections and feature keys should fail clearly."""

    assert normalise_feature_keys(keys={("TYPE", "x")}, context="excluded features") == frozenset(
        {("TYPE", "x")}
    )
    with pytest.raises(InputValidationError, match="feature-key collection"):
        normalise_feature_keys(keys="TYPE:x", context="excluded features")  # type: ignore[arg-type]
    with pytest.raises(InputValidationError, match=r"must be \(feature_type"):
        normalise_feature_keys(
            keys={"TYPE:x"},  # type: ignore[arg-type]
            context="excluded features",
        )


def test_inference_scope_validation_fails_closed() -> None:
    """Explicit unsafe scopes must not enter inference through a caller omission."""

    base = _feature("p1", "TYPE", "x")
    all_data = replace(base, derivation_scope=ALL_DATA_EXPLORATORY)
    with pytest.raises(InputValidationError, match="requires every ALL_DATA_EXPLORATORY"):
        validate_feature_inference_scopes(
            features=(all_data,),
            exploratory_feature_keys=(),
            context="association",
        )
    assert validate_feature_inference_scopes(
        features=(all_data,),
        exploratory_feature_keys={("TYPE", "x")},
        context="association",
    ) == frozenset({("TYPE", "x")})
    supervised = replace(base, derivation_scope=DISCOVERY_SUPERVISED)
    with pytest.raises(InputValidationError, match="audit-only"):
        validate_feature_inference_scopes(
            features=(supervised,),
            exploratory_feature_keys=(),
            context="explainable ML",
        )
    with pytest.raises(InputValidationError, match="do not occur"):
        validate_feature_inference_scopes(
            features=(base,),
            exploratory_feature_keys={("TYPE", "missing")},
            context="association",
        )
    with pytest.raises(InputValidationError, match="collection of strings"):
        normalise_feature_assessment_universes(
            universes={("TYPE", "x"): "p1"},  # type: ignore[dict-item]
            known_protein_ids={"p1"},
            feature_proteins={("TYPE", "x"): {"p1"}},
        )


def test_composer_uses_sequence_domain_and_explicit_assessment_ledgers(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Each modality should receive only its defensible assessed protein universe."""

    features = (
        _feature("p1", "AMINO_ACID_KMER", "k3:AAA"),
        _feature("p1", "PFAM_DOMAIN", "Pfam:PF00646"),
        _feature("p1", "GENERIC", "curated"),
        _feature("p1", "UNSUPPORTED", "presence_only"),
    )
    assessments = (
        _domain_assessment("p1", DomainAssessmentStatus.ASSESSED_WITH_HIT),
        _domain_assessment("p2", DomainAssessmentStatus.ASSESSED_NO_HIT),
        _domain_assessment("p3", DomainAssessmentStatus.NOT_ASSESSED),
    )
    with caplog.at_level(logging.WARNING):
        universes = compose_feature_assessment_universes(
            protein_ids={"p1", "p2", "p3"},
            features=features,
            domain_assessments=assessments,
            explicit_assessment_universes=({("GENERIC", "curated"): {"p1", "p2"}},),
        )
    assert universes[("AMINO_ACID_KMER", "k3:AAA")] == frozenset({"p1", "p2", "p3"})
    assert universes[("PFAM_DOMAIN", "Pfam:PF00646")] == frozenset({"p1", "p2"})
    assert universes[("GENERIC", "curated")] == frozenset({"p1", "p2"})
    assert universes[("UNSUPPORTED", "presence_only")] == frozenset({"p1"})
    assert "1 feature assessment universes contain observed positives only" in caplog.text
    caplog.clear()
    compose_feature_assessment_universes(
        protein_ids={"p1"},
        features=(_feature("p1", "AMINO_ACID_KMER", "k3:ALL"),),
    )
    assert not caplog.records
    with pytest.raises(InputValidationError, match="Universal feature types"):
        compose_feature_assessment_universes(
            protein_ids={"p1"},
            features=features[:1],
            universally_assessed_feature_types="AMINO_ACID_KMER",  # type: ignore[arg-type]
        )


def test_composer_requires_explicit_complete_structure_universes() -> None:
    """Retained structural hits must not imply a no-hit assessment universe."""

    features = (
        _feature("p1", "STRUCTURE_AVAILABLE", "AlphaFoldDB"),
        _feature("p1", "FOLD", "CATH:1.10"),
        _feature(
            "p1",
            "STRUCTURE_CLUSTER",
            "SC_one",
            reference="U1|STRUCTURE_MODEL_RESIDUES|Foldseek|v1",
        ),
    )
    explicit = compose_feature_assessment_universes(
        protein_ids={"p1", "p2"},
        features=features,
        explicit_assessment_universes=(
            {
                ("STRUCTURE_AVAILABLE", "AlphaFoldDB"): {"p1", "p2"},
                ("FOLD", "CATH:1.10"): {"p1", "p2"},
                ("STRUCTURE_CLUSTER", "SC_one"): {"p1", "p2"},
            },
        ),
    )
    assert explicit[("STRUCTURE_AVAILABLE", "AlphaFoldDB")] == frozenset({"p1", "p2"})
    assert explicit[("FOLD", "CATH:1.10")] == frozenset({"p1", "p2"})
    assert explicit[("STRUCTURE_CLUSTER", "SC_one")] == frozenset({"p1", "p2"})
    sparse = compose_feature_assessment_universes(
        protein_ids={"p1", "p2"},
        features=features,
    )
    assert sparse[("STRUCTURE_CLUSTER", "SC_one")] == frozenset({"p1"})
    with pytest.raises(InputValidationError, match="contains unknown proteins"):
        compose_feature_assessment_universes(
            protein_ids={"p1", "p2"},
            features=features,
            explicit_assessment_universes=({("STRUCTURE_CLUSTER", "SC_one"): {"p1", "missing"}},),
        )


def _feature(
    protein_id: str,
    feature_type: str,
    feature_id: str,
    *,
    reference: str = "test",
) -> FeatureRecord:
    """Build one compact positive feature."""

    return FeatureRecord(
        protein_id=protein_id,
        feature_type=feature_type,
        feature_id=feature_id,
        feature_name=feature_id,
        start=None,
        end=None,
        evidence_status="ASSESSED_WITH_HIT",
        evidence_source="test",
        evidence_reference=reference,
    )


def _domain_assessment(protein_id: str, status: DomainAssessmentStatus) -> DomainAssessment:
    """Build one Pfam assessment row."""

    return DomainAssessment(
        protein_id=protein_id,
        domain_authority="Pfam",
        assessment_status=status,
        hit_count=1 if status == DomainAssessmentStatus.ASSESSED_WITH_HIT else 0,
        evidence_source="test",
        evidence_reference="test",
    )
