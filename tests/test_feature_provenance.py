"""Tests for leakage-resistant generic feature provenance."""

from __future__ import annotations

from dataclasses import replace

import pytest

from protein_signatures.checksums import sha256_text
from protein_signatures.errors import InputValidationError
from protein_signatures.feature_provenance import (
    ALL_DATA_EXPLORATORY,
    ASSESSED_NO_FEATURE,
    ASSESSED_WITH_FEATURE,
    DISCOVERY_DERIVED,
    DISCOVERY_SUPERVISED,
    FIXED_EXTERNAL,
    prepare_imported_feature_evidence,
    sequence_cohort_sha256,
)
from protein_signatures.models import FeatureRecord, SequenceRecord


def test_sequence_cohort_digest_is_exact_and_order_independent() -> None:
    """Cohort identity should bind sorted IDs to exact sequence digests."""

    sequences = (_sequence("p2", "b" * 64), _sequence("p1", "a" * 64))
    expected = sha256_text(text=f"p1\t{'a' * 64}\np2\t{'b' * 64}\n")

    assert (
        sequence_cohort_sha256(
            sequences=sequences,
            protein_ids=frozenset({"p2", "p1"}),
        )
        == expected
    )


def test_sequence_cohort_digest_rejects_invalid_authorities() -> None:
    """Cohort hashing must reject scalar, empty, unknown and duplicate IDs."""

    sequences = (_sequence("p1", "a" * 64),)
    for protein_ids, message in (
        ("p1", "must be a collection"),
        (frozenset(), "requires one or more"),
        (frozenset({"missing"}), "absent"),
    ):
        with pytest.raises(InputValidationError, match=message):
            sequence_cohort_sha256(
                sequences=sequences,
                protein_ids=protein_ids,  # type: ignore[arg-type]
            )
    with pytest.raises(InputValidationError, match="must be unique"):
        sequence_cohort_sha256(
            sequences=(sequences[0], sequences[0]),
            protein_ids=frozenset({"p1"}),
        )


def test_imported_feature_preparation_excludes_nonconfirmatory_definitions() -> None:
    """Only fixed or label-blind discovery-derived calls may enter inference."""

    sequences = (_sequence("d1", "a" * 64), _sequence("v1", "b" * 64))
    discovery_ids = frozenset({"d1"})
    discovery_digest = sequence_cohort_sha256(sequences=sequences, protein_ids=discovery_ids)
    full_digest = sequence_cohort_sha256(sequences=sequences, protein_ids=frozenset({"d1", "v1"}))
    records = (
        _feature("d1", "fixed", ASSESSED_WITH_FEATURE, FIXED_EXTERNAL, ""),
        _feature("v1", "fixed", ASSESSED_NO_FEATURE, FIXED_EXTERNAL, ""),
        _feature(
            "d1",
            "learned",
            ASSESSED_WITH_FEATURE,
            DISCOVERY_DERIVED,
            discovery_digest,
        ),
        _feature(
            "v1",
            "learned",
            ASSESSED_NO_FEATURE,
            DISCOVERY_DERIVED,
            discovery_digest,
        ),
        _feature(
            "v1",
            "exploratory",
            ASSESSED_WITH_FEATURE,
            ALL_DATA_EXPLORATORY,
            full_digest,
        ),
        _feature(
            "d1",
            "supervised",
            ASSESSED_WITH_FEATURE,
            DISCOVERY_SUPERVISED,
            discovery_digest,
        ),
        _feature("v1", "unknown", "NOT_ASSESSED", FIXED_EXTERNAL, ""),
    )

    prepared = prepare_imported_feature_evidence(
        records=records,
        sequences=sequences,
        discovery_protein_ids=discovery_ids,
    )

    assert {(item.protein_id, item.feature_id) for item in prepared.confirmatory_features} == {
        ("d1", "fixed"),
        ("d1", "learned"),
    }
    assert prepared.assessment_records == records
    assert prepared.assessment_universes == {
        ("MOTIF", "fixed"): frozenset({"d1", "v1"}),
        ("MOTIF", "learned"): frozenset({"d1", "v1"}),
        ("MOTIF", "unknown"): frozenset(),
    }
    assert prepared.exploratory_feature_keys == frozenset(
        {("MOTIF", "exploratory"), ("MOTIF", "supervised")}
    )


def test_imported_feature_preparation_rejects_false_provenance() -> None:
    """Campaign-derived and fixed declarations must match their cohort policy."""

    sequences = (_sequence("d1", "a" * 64), _sequence("v1", "b" * 64))
    discovery_ids = frozenset({"d1"})
    valid_fixed = _feature("d1", "fixed", ASSESSED_WITH_FEATURE, FIXED_EXTERNAL, "")
    cases = (
        (
            replace(valid_fixed, derivation_cohort_sha256="c" * 64),
            "must leave",
        ),
        (
            _feature(
                "d1",
                "learned",
                ASSESSED_WITH_FEATURE,
                DISCOVERY_DERIVED,
                "c" * 64,
            ),
            "exact discovery cohort",
        ),
        (
            _feature(
                "d1",
                "supervised",
                ASSESSED_WITH_FEATURE,
                DISCOVERY_SUPERVISED,
                "c" * 64,
            ),
            "exact discovery cohort",
        ),
        (
            _feature(
                "d1",
                "exploratory",
                ASSESSED_WITH_FEATURE,
                ALL_DATA_EXPLORATORY,
                "c" * 64,
            ),
            "full campaign cohort",
        ),
        (
            replace(valid_fixed, derivation_scope="UNKNOWN"),
            "Unknown derivation scope",
        ),
    )
    for record, message in cases:
        with pytest.raises(InputValidationError, match=message):
            prepare_imported_feature_evidence(
                records=(record,),
                sequences=sequences,
                discovery_protein_ids=discovery_ids,
            )

    conflicting = replace(valid_fixed, derivation_scope=DISCOVERY_DERIVED)
    with pytest.raises(InputValidationError, match="Conflicting derivation metadata"):
        prepare_imported_feature_evidence(
            records=(valid_fixed, conflicting),
            sequences=sequences,
            discovery_protein_ids=discovery_ids,
        )


def _sequence(protein_id: str, digest: str) -> SequenceRecord:
    """Build one compact authoritative sequence fixture."""

    return SequenceRecord(protein_id, "", "AAAA", 4, digest)


def _feature(
    protein_id: str,
    feature_id: str,
    evidence_status: str,
    derivation_scope: str,
    derivation_cohort_sha256: str,
) -> FeatureRecord:
    """Build one generic feature-assessment fixture."""

    return FeatureRecord(
        protein_id=protein_id,
        feature_type="MOTIF",
        feature_id=feature_id,
        feature_name=feature_id.title(),
        start=None,
        end=None,
        evidence_status=evidence_status,
        evidence_source="test",
        evidence_reference="definition-v1",
        derivation_scope=derivation_scope,
        feature_definition_sha256="f" * 64,
        derivation_cohort_sha256=derivation_cohort_sha256,
    )
