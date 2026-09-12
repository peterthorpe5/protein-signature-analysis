"""Defensive parser tests for every canonical scientific input table."""

from __future__ import annotations

from pathlib import Path

import pytest

from protein_signatures.checksums import sha256_text
from protein_signatures.errors import InputValidationError
from protein_signatures.models import (
    DomainAssessment,
    DomainAssessmentStatus,
    DomainHit,
    StructureAnalysisEligibility,
    StructureCoverageScope,
)
from protein_signatures.tables import (
    DOMAIN_FIELDS,
    FEATURE_FIELDS,
    LABEL_FIELDS,
    STRUCTURE_COMPARISON_FIELDS,
    STRUCTURE_FIELDS,
    _parse_coordinates,
    _reject_nonblank_domain_payload,
    _resolve_coordinate_path,
    _validate_known_protein,
    complete_domain_assessments,
    derive_domain_evidence,
    read_domains,
    read_features,
    read_label_assignments,
    read_structure_comparisons,
    read_structures,
)


def test_label_assignment_rejects_unknown_duplicate_and_status(
    tmp_path: Path,
) -> None:
    """Label assignments must reference known proteins, labels and curation states."""

    valid = {
        "protein_id": "p1",
        "label_id": "label",
        "curation_status": "REVIEWED_POSITIVE",
        "evidence_status": "REVIEWED",
        "evidence_source": "test",
        "evidence_reference": "",
        "component_role": "ROLE",
        "curation_reason": "reviewed",
    }
    path = _table(tmp_path / "labels.tsv", LABEL_FIELDS, (valid, valid))
    with pytest.raises(InputValidationError, match="Duplicate"):
        read_label_assignments(
            path=path,
            protein_ids=frozenset({"p1"}),
            label_ids=frozenset({"label"}),
        )
    for field, value, message in (
        ("protein_id", "missing", "absent"),
        ("label_id", "missing", "Unknown profile"),
        ("curation_status", "invented", "Unknown curation"),
    ):
        row = {**valid, field: value}
        path = _table(tmp_path / f"label_{field}.tsv", LABEL_FIELDS, (row,))
        with pytest.raises(InputValidationError, match=message):
            read_label_assignments(
                path=path,
                protein_ids=frozenset({"p1"}),
                label_ids=frozenset({"label"}),
            )


def test_feature_parser_validates_coordinates_and_uniqueness(tmp_path: Path) -> None:
    """Features should support blank coordinates but reject incomplete or duplicate spans."""

    row = {
        "protein_id": "p1",
        "feature_type": "MOTIF",
        "feature_id": "m1",
        "feature_name": "Motif one",
        "start": "",
        "end": "",
        "evidence_status": "ASSESSED_WITH_FEATURE",
        "evidence_source": "test",
        "evidence_reference": "",
        "derivation_scope": "FIXED_EXTERNAL",
        "feature_definition_sha256": "a" * 64,
        "derivation_cohort_sha256": "",
    }
    path = _table(tmp_path / "features.tsv", FEATURE_FIELDS, (row,))
    assert read_features(path=path, sequences=_sequences())[0].start is None
    duplicate = _table(tmp_path / "duplicate.tsv", FEATURE_FIELDS, (row, row))
    with pytest.raises(InputValidationError, match="Duplicate feature"):
        read_features(path=duplicate, sequences=_sequences())
    for start, end in (("1", ""), ("3", "2"), ("1", "9")):
        invalid = _table(
            tmp_path / f"feature_{start}_{end or 'blank'}.tsv",
            FEATURE_FIELDS,
            ({**row, "start": start, "end": end},),
        )
        with pytest.raises(InputValidationError):
            read_features(path=invalid, sequences=_sequences())
    with pytest.raises(InputValidationError, match="absent"):
        _validate_known_protein(protein_id="missing", protein_ids=frozenset({"p1"}))
    assert _parse_coordinates(start="", end="", protein_id="p1", length=4) == (None, None)


def test_feature_parser_enforces_assessment_and_derivation_contract(tmp_path: Path) -> None:
    """Generic features must carry consistent machine-checkable provenance."""

    base = {
        "protein_id": "p1",
        "feature_type": "MOTIF",
        "feature_id": "m1",
        "feature_name": "Motif one",
        "start": "",
        "end": "",
        "evidence_status": "ASSESSED_NO_FEATURE",
        "evidence_source": "test",
        "evidence_reference": "definition-v1",
        "derivation_scope": "FIXED_EXTERNAL",
        "feature_definition_sha256": "a" * 64,
        "derivation_cohort_sha256": "",
    }
    cases = (
        ({**base, "evidence_status": "PRESENT"}, "Unknown generic feature"),
        ({**base, "start": "1", "end": "2"}, "Only ASSESSED_WITH_FEATURE"),
        ({**base, "derivation_scope": "UNKNOWN"}, "Unknown derivation_scope"),
        ({**base, "feature_definition_sha256": "A" * 64}, "64 lower-case"),
        ({**base, "derivation_cohort_sha256": "b" * 64}, "must leave"),
        (
            {
                **base,
                "derivation_scope": "DISCOVERY_DERIVED",
                "derivation_cohort_sha256": "",
            },
            "Campaign-derived",
        ),
    )
    for index, (row, message) in enumerate(cases):
        path = _table(tmp_path / f"invalid_feature_{index}.tsv", FEATURE_FIELDS, (row,))
        with pytest.raises(InputValidationError, match=message):
            read_features(path=path, sequences=_sequences())

    conflicting_state = _table(
        tmp_path / "feature_state.tsv",
        FEATURE_FIELDS,
        (
            base,
            {**base, "evidence_status": "ASSESSED_WITH_FEATURE", "start": "1", "end": "2"},
        ),
    )
    with pytest.raises(InputValidationError, match="Conflicting assessment states"):
        read_features(path=conflicting_state, sequences=_sequences())
    conflicting_definition = _table(
        tmp_path / "feature_definition.tsv",
        FEATURE_FIELDS,
        (
            base,
            {
                **base,
                "protein_id": "p2",
                "feature_definition_sha256": "b" * 64,
            },
        ),
    )
    with pytest.raises(InputValidationError, match="Conflicting definition metadata"):
        read_features(path=conflicting_definition, sequences=_sequences())


def test_domain_parser_distinguishes_hits_no_hits_and_failures(tmp_path: Path) -> None:
    """Domain state, payload, provenance and hit identities must remain consistent."""

    hit = _domain_row()
    no_hit = {
        **_domain_row(),
        "protein_id": "p2",
        "assessment_status": "ASSESSED_NO_HIT",
        "domain_id": "",
        "domain_name": "",
        "start": "",
        "end": "",
        "score": "",
        "e_value": "",
    }
    path = _table(tmp_path / "domains.tsv", DOMAIN_FIELDS, (hit, no_hit))
    hits, assessments = read_domains(path=path, sequences=_sequences())
    assert len(hits) == 1
    assert {item.assessment_status for item in assessments} == {
        DomainAssessmentStatus.ASSESSED_WITH_HIT,
        DomainAssessmentStatus.ASSESSED_NO_HIT,
    }
    cases = (
        ({**hit, "assessment_status": "invented"}, "Unknown domain"),
        ({**hit, "domain_id": ""}, "lacks domain_id"),
        ({**no_hit, "domain_id": "PF1"}, "Non-hit"),
        ({**no_hit, "domain_name": "payload"}, "contains hit payload"),
        ({**hit, "start": ""}, "coordinates"),
    )
    for index, (row, message) in enumerate(cases):
        invalid = _table(tmp_path / f"domain_bad_{index}.tsv", DOMAIN_FIELDS, (row,))
        with pytest.raises(InputValidationError, match=message):
            read_domains(path=invalid, sequences=_sequences())
    duplicate = _table(tmp_path / "domain_duplicate.tsv", DOMAIN_FIELDS, (hit, hit))
    with pytest.raises(InputValidationError, match="Duplicate domain"):
        read_domains(path=duplicate, sequences=_sequences())
    conflicting_state = _table(
        tmp_path / "domain_state.tsv",
        DOMAIN_FIELDS,
        (hit, {**no_hit, "protein_id": "p1"}),
    )
    with pytest.raises(InputValidationError, match="Conflicting domain assessment"):
        read_domains(path=conflicting_state, sequences=_sequences())
    second_hit = {**hit, "domain_id": "PF2", "start": "2", "evidence_reference": "other"}
    conflicting_provenance = _table(
        tmp_path / "domain_provenance.tsv", DOMAIN_FIELDS, (hit, second_hit)
    )
    with pytest.raises(InputValidationError, match="Conflicting domain provenance"):
        read_domains(path=conflicting_provenance, sequences=_sequences())


def test_domain_derivation_and_completion_cover_other_authorities() -> None:
    """Domain evidence should include extracted sequence, architecture and missing states."""

    hit = DomainHit("p1", "InterPro", "IPR1", "Domain", 2, 3, 10.0, 1e-4, "test", "x")
    features, sequences = derive_domain_evidence(hits=(hit,), sequences=_sequences())
    assert {item.feature_type for item in features} == {"DOMAIN", "DOMAIN_ARCHITECTURE"}
    assert sequences[0]["domain_sequence"] == "CD"
    supplied = DomainAssessment(
        "p1",
        "InterPro",
        DomainAssessmentStatus.ASSESSED_WITH_HIT,
        1,
        "test",
        "x",
    )
    completed = complete_domain_assessments(
        protein_ids=frozenset({"p1", "p2"}), assessments=(supplied,)
    )
    assert len(completed) == 4
    assert any(
        item.domain_authority == "Pfam"
        and item.assessment_status == DomainAssessmentStatus.NOT_ASSESSED
        for item in completed
    )
    with pytest.raises(InputValidationError, match="hit payload"):
        _reject_nonblank_domain_payload(
            row={
                "domain_name": "x",
                "start": "",
                "end": "",
                "score": "",
                "e_value": "",
            },
            key=("p1", "Pfam"),
        )


def test_structure_inventory_checks_paths_digests_and_ids(tmp_path: Path) -> None:
    """Structure rows must bind coordinate bytes to unique identifiers."""

    coordinate = tmp_path / "model.pdb"
    coordinate.write_text("ATOM\n", encoding="utf-8")
    row = _structure_row(coordinate_path="model.pdb")
    path = _table(tmp_path / "structures.tsv", STRUCTURE_FIELDS, (row,))
    structures = read_structures(path=path, sequences=_sequences())
    assert structures[0].coordinate_path == coordinate
    assert structures[0].coordinate_sha256 == sha256_text(text="ATOM\n")
    assert structures[0].is_coordinate_analysis_eligible is True
    assert structures[0].comparison_universe_ids == ("campaign_1", "campaign_2")
    assert structures[0].analysis_eligibility_status == StructureAnalysisEligibility.ELIGIBLE
    assert structures[0].fold_authority_version == "4.3.0"
    assert structures[0].fold_evidence_reference == "CATH:4.3.0:fold1"
    for name, rows, message in (
        ("duplicate", (row, row), "Duplicate structure"),
        ("digest", ({**row, "coordinate_sha256": "0" * 64},), "checksum mismatch"),
        (
            "orphan_digest",
            ({**row, "coordinate_path": "", "coordinate_sha256": "0" * 64},),
            "without coordinate_path",
        ),
        ("missing", ({**row, "coordinate_path": "missing.pdb"},), "Missing or empty"),
    ):
        invalid = _table(tmp_path / f"structure_{name}.tsv", STRUCTURE_FIELDS, rows)
        with pytest.raises(InputValidationError, match=message):
            read_structures(path=invalid, sequences=_sequences())
    assert _resolve_coordinate_path(value="", table_path=path) is None

    duplicate_universe = {
        **row,
        "comparison_universe_ids": "campaign_1|campaign_1",
    }
    invalid = _table(tmp_path / "duplicate_universe.tsv", STRUCTURE_FIELDS, (duplicate_universe,))
    with pytest.raises(InputValidationError, match="duplicate identifier"):
        read_structures(path=invalid, sequences=_sequences())

    external_row = {
        **row,
        "coordinate_path": "",
        "availability_status": "INPUT_UNAVAILABLE",
        "analysis_eligibility_status": "NOT_APPLICABLE_EXTERNAL_EVIDENCE",
    }
    external = _table(tmp_path / "external_fold.tsv", STRUCTURE_FIELDS, (external_row,))
    assert read_structures(path=external, sequences=_sequences())[0].coordinate_path is None


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"availability_status": "UNKNOWN"}, "Unsupported availability_status"),
        ({"availability_status": "FAILED"}, "inconsistent"),
        (
            {"analysis_eligibility_status": "NOT_APPLICABLE_EXTERNAL_EVIDENCE"},
            "external-evidence eligibility",
        ),
        (
            {"analysis_eligibility_status": "INELIGIBLE_COORDINATE_UNAVAILABLE"},
            "declares unavailable coordinates",
        ),
        (
            {
                "analysis_eligibility_status": "INELIGIBLE_LOW_CONFIDENCE",
                "mean_confidence": "",
            },
            "low confidence",
        ),
        (
            {"analysis_eligibility_status": "INELIGIBLE_CONFIDENCE_UNAVAILABLE"},
            "confidence unavailable",
        ),
        ({"analysis_eligibility_status": "UNCONTROLLED"}, "must be one of"),
    ],
)
def test_structure_inventory_rejects_inconsistent_controlled_states(
    tmp_path: Path, updates: dict[str, str], message: str
) -> None:
    """Eligibility, physical availability and confidence must not contradict."""

    coordinate = tmp_path / "model.pdb"
    coordinate.write_text("ATOM\n", encoding="utf-8")
    row = {**_structure_row(coordinate_path="model.pdb"), **updates}
    path = _table(tmp_path / "invalid_structure.tsv", STRUCTURE_FIELDS, (row,))
    with pytest.raises(InputValidationError, match=message):
        read_structures(path=path, sequences=_sequences())

    no_coordinate = {
        **_structure_row(coordinate_path=""),
        "availability_status": "INPUT_UNAVAILABLE",
        "analysis_eligibility_status": "ELIGIBLE",
        "comparison_universe_ids": "campaign_2|campaign_1",
    }
    path = _table(tmp_path / "eligible_without_coordinate.tsv", STRUCTURE_FIELDS, (no_coordinate,))
    with pytest.raises(InputValidationError, match="cannot be ELIGIBLE"):
        read_structures(path=path, sequences=_sequences())


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"fold_evidence_status": "UNKNOWN"}, "must be one of"),
        ({"fold_id": ""}, "requires fold_id"),
        ({"fold_authority": ""}, "requires fold_id"),
        ({"fold_authority_version": ""}, "requires fold_id"),
        ({"fold_evidence_reference": ""}, "requires fold_id"),
        ({"fold_evidence_status": "ASSESSED_NO_HIT"}, "must not carry"),
        (
            {
                "fold_evidence_status": "FAILED",
                "fold_id": "",
                "fold_name": "",
                "fold_authority": "",
            },
            "requires fold_authority",
        ),
        (
            {
                "fold_evidence_status": "NOT_ASSESSED",
                "fold_id": "",
                "fold_name": "",
            },
            "must not carry fold authority",
        ),
    ],
)
def test_structure_inventory_rejects_inconsistent_fold_payloads(
    tmp_path: Path, updates: dict[str, str], message: str
) -> None:
    """Fold hit, no-hit and failure states must carry compatible fields."""

    coordinate = tmp_path / "model.pdb"
    coordinate.write_text("ATOM\n", encoding="utf-8")
    row = {**_structure_row(coordinate_path="model.pdb"), **updates}
    path = _table(tmp_path / "invalid_fold.tsv", STRUCTURE_FIELDS, (row,))
    with pytest.raises(InputValidationError, match=message):
        read_structures(path=path, sequences=_sequences())

    assessed_no_hit = {
        **_structure_row(coordinate_path="model.pdb"),
        "fold_evidence_status": "ASSESSED_NO_HIT",
        "fold_id": "",
        "fold_name": "",
    }
    path = _table(tmp_path / "fold_no_hit.tsv", STRUCTURE_FIELDS, (assessed_no_hit,))
    assert read_structures(path=path, sequences=_sequences())[0].fold_id == ""

    mixed_release = {
        **_structure_row(coordinate_path="model.pdb"),
        "protein_id": "p2",
        "structure_id": "s2",
        "fold_authority_version": "4.2.0",
    }
    path = _table(
        tmp_path / "mixed_fold_release.tsv",
        STRUCTURE_FIELDS,
        (_structure_row(coordinate_path="model.pdb"), mixed_release),
    )
    with pytest.raises(InputValidationError, match="multiple releases"):
        read_structures(path=path, sequences=_sequences())


def test_structure_comparisons_reject_self_duplicate_and_unknown(tmp_path: Path) -> None:
    """Pairwise structural records must describe unique known unordered pairs."""

    row = _comparison_row()
    path = _table(tmp_path / "comparisons.tsv", STRUCTURE_COMPARISON_FIELDS, (row,))
    parsed = read_structure_comparisons(path=path, protein_ids=frozenset({"p1", "p2"}))
    assert parsed[0].protein_a_id == "p1"
    assert parsed[0].coverage_a == pytest.approx(0.9)
    assert parsed[0].coverage_b == pytest.approx(0.7)
    assert parsed[0].coverage_scope == StructureCoverageScope.FULL_SEQUENCE
    assert parsed[0].comparison_universe_id == "campaign_1"
    for name, rows, message in (
        ("self", ({**row, "protein_a_id": "p1", "protein_b_id": "p1"},), "Self"),
        ("duplicate", (row, {**row, "protein_a_id": "p2", "protein_b_id": "p1"}), "Duplicate"),
        ("unknown", ({**row, "protein_b_id": "missing"},), "absent"),
        ("scope", ({**row, "coverage_scope": "UNKNOWN"},), "coverage_scope"),
        ("universe", ({**row, "comparison_universe_id": ""},), "comparison_universe_id"),
    ):
        invalid = _table(tmp_path / f"comparison_{name}.tsv", STRUCTURE_COMPARISON_FIELDS, rows)
        with pytest.raises(InputValidationError, match=message):
            read_structure_comparisons(
                path=invalid,
                protein_ids=frozenset({"p1", "p2"}),
            )


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"comparison_status": "UNKNOWN"}, "comparison_status must be one of"),
        ({"tm_score": ""}, "requires tm_score"),
        ({"coverage_a": ""}, "requires tm_score"),
        ({"comparison_status": "FAILED"}, "must not carry metrics"),
    ],
)
def test_structure_comparisons_enforce_status_payload_consistency(
    tmp_path: Path, updates: dict[str, str], message: str
) -> None:
    """Successful rows require clustering metrics and failures must carry none."""

    row = {**_comparison_row(), **updates}
    path = _table(tmp_path / "bad_comparison_payload.tsv", STRUCTURE_COMPARISON_FIELDS, (row,))
    with pytest.raises(InputValidationError, match=message):
        read_structure_comparisons(path=path, protein_ids=frozenset({"p1", "p2"}))

    failed = {
        **_comparison_row(),
        "comparison_status": "FAILED",
        "tm_score": "",
        "rmsd_angstrom": "",
        "aligned_residue_count": "",
        "coverage_a": "",
        "coverage_b": "",
    }
    path = _table(tmp_path / "failed_comparison.tsv", STRUCTURE_COMPARISON_FIELDS, (failed,))
    assert (
        read_structure_comparisons(path=path, protein_ids=frozenset({"p1", "p2"}))[
            0
        ].comparison_status.value
        == "FAILED"
    )


def _table(path: Path, fields: tuple[str, ...], rows: tuple[dict[str, str], ...]) -> Path:
    """Write a small TSV using the exact canonical field order."""

    lines = ["\t".join(fields)]
    lines.extend("\t".join(row.get(field, "") for field in fields) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _sequences() -> tuple[object, ...]:
    """Return two short authoritative protein records."""

    from protein_signatures.models import SequenceRecord

    return (
        SequenceRecord("p1", "", "ACDE", 4, "a" * 64),
        SequenceRecord("p2", "", "AAAA", 4, "b" * 64),
    )


def _domain_row() -> dict[str, str]:
    """Return one valid Pfam hit row."""

    return {
        "protein_id": "p1",
        "domain_authority": "Pfam",
        "assessment_status": "ASSESSED_WITH_HIT",
        "domain_id": "PF1",
        "domain_name": "Domain",
        "start": "1",
        "end": "2",
        "score": "10",
        "e_value": "0.001",
        "evidence_source": "test",
        "evidence_reference": "x",
    }


def _structure_row(*, coordinate_path: str) -> dict[str, str]:
    """Return one valid structure row."""

    return {
        "protein_id": "p1",
        "structure_id": "s1",
        "structure_source": "USER",
        "structure_version": "v1",
        "coordinate_path": coordinate_path,
        "coordinate_sha256": "",
        "availability_status": "AVAILABLE",
        "mean_confidence": "80",
        "fold_id": "fold1",
        "fold_name": "Fold one",
        "fold_authority": "CATH",
        "fold_authority_version": "4.3.0",
        "fold_evidence_reference": "CATH:4.3.0:fold1",
        "fold_evidence_status": "ASSESSED_WITH_HIT",
        "analysis_eligibility_status": "ELIGIBLE",
        "comparison_universe_ids": "campaign_2|campaign_1",
    }


def _comparison_row() -> dict[str, str]:
    """Return one valid pairwise structural comparison row."""

    return {
        "protein_a_id": "p2",
        "protein_b_id": "p1",
        "comparison_tool": "Foldseek",
        "comparison_tool_version": "1",
        "tm_score": "0.8",
        "rmsd_angstrom": "1.2",
        "aligned_residue_count": "3",
        "coverage_a": "0.7",
        "coverage_b": "0.9",
        "comparison_status": "COMPLETE",
        "source_record_id": "row1",
        "comparison_universe_id": "campaign_1",
        "coverage_scope": "FULL_SEQUENCE",
    }
