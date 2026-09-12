"""Strict parsers and derivations for campaign TSV input tables."""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from pathlib import Path

from .checksums import sha256_file
from .errors import InputValidationError
from .feature_provenance import (
    ASSESSED_WITH_FEATURE,
    FEATURE_ASSESSMENT_STATUSES,
    FEATURE_DERIVATION_SCOPES,
    FIXED_EXTERNAL,
)
from .io_utils import iter_tsv
from .models import (
    CurationStatus,
    DomainAssessment,
    DomainAssessmentStatus,
    DomainHit,
    FeatureRecord,
    FoldEvidenceStatus,
    LabelAssignment,
    PairwiseStructureComparison,
    SequenceRecord,
    StructureAnalysisEligibility,
    StructureComparisonStatus,
    StructureCoverageScope,
    StructureRecord,
)
from .validation import (
    parse_optional_float,
    parse_optional_integer,
    validate_identifier,
    validate_text,
)

LOGGER = logging.getLogger(__name__)

LABEL_FIELDS = (
    "protein_id",
    "label_id",
    "curation_status",
    "evidence_status",
    "evidence_source",
    "evidence_reference",
    "component_role",
    "curation_reason",
)
FEATURE_FIELDS = (
    "protein_id",
    "feature_type",
    "feature_id",
    "feature_name",
    "start",
    "end",
    "evidence_status",
    "evidence_source",
    "evidence_reference",
    "derivation_scope",
    "feature_definition_sha256",
    "derivation_cohort_sha256",
)
DOMAIN_FIELDS = (
    "protein_id",
    "domain_authority",
    "assessment_status",
    "domain_id",
    "domain_name",
    "start",
    "end",
    "score",
    "e_value",
    "evidence_source",
    "evidence_reference",
)
STRUCTURE_FIELDS = (
    "protein_id",
    "structure_id",
    "structure_source",
    "structure_version",
    "coordinate_path",
    "coordinate_sha256",
    "availability_status",
    "mean_confidence",
    "fold_id",
    "fold_name",
    "fold_authority",
    "fold_authority_version",
    "fold_evidence_reference",
    "fold_evidence_status",
    "analysis_eligibility_status",
    "comparison_universe_ids",
)
STRUCTURE_COMPARISON_FIELDS = (
    "protein_a_id",
    "protein_b_id",
    "comparison_tool",
    "comparison_tool_version",
    "tm_score",
    "rmsd_angstrom",
    "aligned_residue_count",
    "coverage_a",
    "coverage_b",
    "comparison_status",
    "source_record_id",
    "comparison_universe_id",
    "coverage_scope",
)

_AVAILABLE_COORDINATE_STATUSES = frozenset({"AVAILABLE", "COMPLETE"})
_STRUCTURE_AVAILABILITY_STATUSES = frozenset(
    {
        *_AVAILABLE_COORDINATE_STATUSES,
        "FAILED",
        "INPUT_UNAVAILABLE",
        "MODEL_NOT_AVAILABLE",
        "NOT_ASSESSED",
        "NOT_SELECTED",
    }
)
_SUCCESSFUL_STRUCTURE_COMPARISON_STATUSES = frozenset(
    {
        StructureComparisonStatus.COMPLETE,
        StructureComparisonStatus.PASS,
        StructureComparisonStatus.SUCCESS,
    }
)


def read_label_assignments(
    *, path: Path, protein_ids: frozenset[str], label_ids: frozenset[str]
) -> tuple[LabelAssignment, ...]:
    """Read curated protein-to-label assignments.

    Args:
        path: Assignment TSV.
        protein_ids: Authoritative FASTA identifiers.
        label_ids: Canonical profile identifiers.

    Returns:
        Deterministically ordered assignments.

    Raises:
        InputValidationError: If a row is invalid, unknown or duplicated.
    """

    assignments: list[LabelAssignment] = []
    observed: set[tuple[str, str]] = set()
    for row in iter_tsv(path=path, required_fields=LABEL_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        label_id = validate_identifier(value=row["label_id"], field_name="label_id")
        _validate_known_protein(protein_id=protein_id, protein_ids=protein_ids)
        if label_id not in label_ids:
            raise InputValidationError(f"Unknown profile label in assignments: {label_id!r}")
        key = (protein_id, label_id)
        if key in observed:
            raise InputValidationError(f"Duplicate protein/label assignment: {key!r}")
        observed.add(key)
        try:
            status = CurationStatus(row["curation_status"].strip().upper())
        except ValueError as error:
            raise InputValidationError(
                f"Unknown curation_status for {protein_id!r}: {row['curation_status']!r}"
            ) from error
        assignments.append(
            LabelAssignment(
                protein_id=protein_id,
                label_id=label_id,
                curation_status=status,
                evidence_status=validate_text(
                    value=row["evidence_status"], field_name="evidence_status"
                ),
                evidence_source=validate_text(
                    value=row["evidence_source"], field_name="evidence_source"
                ),
                evidence_reference=validate_text(
                    value=row["evidence_reference"],
                    field_name="evidence_reference",
                    allow_empty=True,
                ),
                component_role=validate_text(
                    value=row["component_role"], field_name="component_role"
                ),
                curation_reason=validate_text(
                    value=row["curation_reason"], field_name="curation_reason"
                ),
            )
        )
    result = tuple(sorted(assignments, key=lambda item: (item.protein_id, item.label_id)))
    LOGGER.info("Loaded %d curated protein-label assignments", len(result))
    return result


def read_features(
    *, path: Path, sequences: tuple[SequenceRecord, ...]
) -> tuple[FeatureRecord, ...]:
    """Read explicit assessment rows for externally generated features.

    Args:
        path: Feature TSV.
        sequences: Authoritative protein inventory.

    Returns:
        Deterministically ordered presence, absence and unknown assessment rows.

    Raises:
        InputValidationError: If assessment, coordinates or derivation metadata
            violate the generic feature contract.
    """
    lengths = {item.protein_id: item.sequence_length for item in sequences}
    features: list[FeatureRecord] = []
    observed: set[tuple[str, str, str, int | None, int | None]] = set()
    status_by_protein_feature: dict[tuple[str, str, str], str] = {}
    definition_by_feature: dict[tuple[str, str], tuple[str, str, str, str]] = {}
    for row in iter_tsv(path=path, required_fields=FEATURE_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        _validate_known_protein(protein_id=protein_id, protein_ids=frozenset(lengths))
        start, end = _parse_coordinates(
            start=row["start"], end=row["end"], protein_id=protein_id, length=lengths[protein_id]
        )
        feature_type = validate_identifier(value=row["feature_type"], field_name="feature_type")
        feature_id = validate_identifier(value=row["feature_id"], field_name="feature_id")
        feature_name = validate_text(value=row["feature_name"], field_name="feature_name")
        evidence_status = validate_identifier(
            value=row["evidence_status"], field_name="evidence_status"
        ).upper()
        if evidence_status not in FEATURE_ASSESSMENT_STATUSES:
            raise InputValidationError(
                f"Unknown generic feature evidence_status for {protein_id!r}: {evidence_status!r}."
            )
        if evidence_status != ASSESSED_WITH_FEATURE and (start is not None or end is not None):
            raise InputValidationError(
                "Only ASSESSED_WITH_FEATURE rows may contain feature coordinates: "
                f"{(protein_id, feature_type, feature_id)!r}."
            )
        derivation_scope = validate_identifier(
            value=row["derivation_scope"], field_name="derivation_scope"
        ).upper()
        if derivation_scope not in FEATURE_DERIVATION_SCOPES:
            raise InputValidationError(
                f"Unknown derivation_scope for generic feature "
                f"{(feature_type, feature_id)!r}: {derivation_scope!r}."
            )
        definition_digest = row["feature_definition_sha256"].strip()
        if len(definition_digest) != 64 or any(
            character not in "0123456789abcdef" for character in definition_digest
        ):
            raise InputValidationError(
                "feature_definition_sha256 must be exactly 64 lower-case hexadecimal "
                f"characters for {(feature_type, feature_id)!r}."
            )
        cohort_digest = row["derivation_cohort_sha256"].strip()
        if derivation_scope == FIXED_EXTERNAL:
            if cohort_digest:
                raise InputValidationError(
                    "FIXED_EXTERNAL generic features must leave "
                    f"derivation_cohort_sha256 blank: {(feature_type, feature_id)!r}."
                )
        elif len(cohort_digest) != 64 or any(
            character not in "0123456789abcdef" for character in cohort_digest
        ):
            raise InputValidationError(
                "Campaign-derived generic features require a 64-character lower-case "
                f"derivation_cohort_sha256: {(feature_type, feature_id)!r}."
            )
        definition_key = (feature_type, feature_id)
        definition = (
            feature_name,
            derivation_scope,
            definition_digest,
            cohort_digest,
        )
        previous_definition = definition_by_feature.setdefault(definition_key, definition)
        if previous_definition != definition:
            raise InputValidationError(
                f"Conflicting definition metadata for generic feature {definition_key!r}."
            )
        protein_feature_key = (protein_id, feature_type, feature_id)
        previous_status = status_by_protein_feature.setdefault(protein_feature_key, evidence_status)
        if previous_status != evidence_status:
            raise InputValidationError(
                f"Conflicting assessment states for generic feature {protein_feature_key!r}."
            )
        key = (protein_id, feature_type, feature_id, start, end)
        if key in observed:
            raise InputValidationError(f"Duplicate feature record: {key!r}")
        observed.add(key)
        features.append(
            FeatureRecord(
                protein_id=protein_id,
                feature_type=feature_type,
                feature_id=feature_id,
                feature_name=feature_name,
                start=start,
                end=end,
                evidence_status=evidence_status,
                evidence_source=validate_text(
                    value=row["evidence_source"], field_name="evidence_source"
                ),
                evidence_reference=validate_text(
                    value=row["evidence_reference"],
                    field_name="evidence_reference",
                    allow_empty=True,
                ),
                derivation_scope=derivation_scope,
                feature_definition_sha256=definition_digest,
                derivation_cohort_sha256=cohort_digest,
            )
        )
    result = tuple(
        sorted(
            features,
            key=lambda item: (
                item.protein_id,
                item.feature_type,
                item.feature_id,
                item.start or 0,
                item.end or 0,
            ),
        )
    )
    LOGGER.info("Loaded %d external protein-feature assessment rows", len(result))
    return result


def read_domains(
    *, path: Path, sequences: tuple[SequenceRecord, ...]
) -> tuple[tuple[DomainHit, ...], tuple[DomainAssessment, ...]]:
    """Read coordinate-resolved domain hits and explicit assessment states.

    Args:
        path: Domain TSV. Hit rows and no-hit/failure sentinel rows share one schema.
        sequences: Authoritative protein inventory.

    Returns:
        Domain hits and aggregated protein/authority assessments.

    Raises:
        InputValidationError: If state, coordinate or uniqueness contracts fail.
    """

    lengths = {item.protein_id: item.sequence_length for item in sequences}
    hits: list[DomainHit] = []
    states: dict[tuple[str, str], DomainAssessmentStatus] = {}
    provenance: dict[tuple[str, str], tuple[str, str]] = {}
    observed_hits: set[tuple[str, str, str, int, int]] = set()
    for row in iter_tsv(path=path, required_fields=DOMAIN_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        _validate_known_protein(protein_id=protein_id, protein_ids=frozenset(lengths))
        authority = validate_identifier(
            value=row["domain_authority"], field_name="domain_authority"
        )
        try:
            status = DomainAssessmentStatus(row["assessment_status"].strip().upper())
        except ValueError as error:
            raise InputValidationError(
                f"Unknown domain assessment_status for {protein_id!r}: {row['assessment_status']!r}"
            ) from error
        key = (protein_id, authority)
        previous = states.get(key)
        if previous is not None and previous != status:
            raise InputValidationError(f"Conflicting domain assessment states for {key!r}")
        states[key] = status
        source = validate_text(value=row["evidence_source"], field_name="evidence_source")
        reference = validate_text(
            value=row["evidence_reference"],
            field_name="evidence_reference",
            allow_empty=True,
        )
        if key in provenance and provenance[key] != (source, reference):
            raise InputValidationError(f"Conflicting domain provenance for {key!r}")
        provenance[key] = (source, reference)
        has_domain = bool(row["domain_id"].strip())
        if status == DomainAssessmentStatus.ASSESSED_WITH_HIT and not has_domain:
            raise InputValidationError(f"Domain hit state lacks domain_id for {key!r}")
        if status != DomainAssessmentStatus.ASSESSED_WITH_HIT and has_domain:
            raise InputValidationError(f"Non-hit domain state contains domain_id for {key!r}")
        if not has_domain:
            _reject_nonblank_domain_payload(row=row, key=key)
            continue
        start, end = _parse_coordinates(
            start=row["start"], end=row["end"], protein_id=protein_id, length=lengths[protein_id]
        )
        if start is None or end is None:
            raise InputValidationError(f"Domain hit lacks coordinates for {key!r}")
        domain_id = validate_identifier(value=row["domain_id"], field_name="domain_id")
        hit_key = (protein_id, authority, domain_id, start, end)
        if hit_key in observed_hits:
            raise InputValidationError(f"Duplicate domain hit: {hit_key!r}")
        observed_hits.add(hit_key)
        hits.append(
            DomainHit(
                protein_id=protein_id,
                domain_authority=authority,
                domain_id=domain_id,
                domain_name=validate_text(value=row["domain_name"], field_name="domain_name"),
                start=start,
                end=end,
                score=parse_optional_float(value=row["score"], field_name="score"),
                e_value=parse_optional_float(
                    value=row["e_value"], field_name="e_value", minimum=0.0
                ),
                evidence_source=source,
                evidence_reference=reference,
            )
        )
    hit_counts = Counter((item.protein_id, item.domain_authority) for item in hits)
    assessments = tuple(
        DomainAssessment(
            protein_id=key[0],
            domain_authority=key[1],
            assessment_status=status,
            hit_count=hit_counts[key],
            evidence_source=provenance[key][0],
            evidence_reference=provenance[key][1],
        )
        for key, status in sorted(states.items())
    )
    for assessment in assessments:
        expected_hits = assessment.assessment_status == DomainAssessmentStatus.ASSESSED_WITH_HIT
        if expected_hits != (assessment.hit_count > 0):
            raise InputValidationError(
                "Domain assessment/hit-count mismatch for "
                f"{assessment.protein_id!r}/{assessment.domain_authority!r}."
            )
    ordered_hits = tuple(
        sorted(
            hits,
            key=lambda item: (
                item.protein_id,
                item.domain_authority,
                item.start,
                item.end,
            ),
        )
    )
    LOGGER.info("Loaded %d domain hits and %d explicit assessments", len(hits), len(assessments))
    return ordered_hits, assessments


def complete_domain_assessments(
    *,
    protein_ids: frozenset[str],
    assessments: tuple[DomainAssessment, ...],
    default_authorities: tuple[str, ...] = ("Pfam",),
) -> tuple[DomainAssessment, ...]:
    """Add explicit ``NOT_ASSESSED`` rows for missing protein/authority pairs.

    Args:
        protein_ids: Authoritative protein identifiers.
        assessments: Explicitly supplied assessments.
        default_authorities: Authorities that must always appear in the coverage matrix.

    Returns:
        Complete, deterministic assessment matrix.
    """

    authorities = {
        validate_identifier(value=item, field_name="domain authority")
        for item in default_authorities
    }
    authorities.update(item.domain_authority for item in assessments)
    by_key = {(item.protein_id, item.domain_authority): item for item in assessments}
    result: list[DomainAssessment] = []
    for protein_id in sorted(protein_ids):
        for authority in sorted(authorities, key=str.casefold):
            result.append(
                by_key.get(
                    (protein_id, authority),
                    DomainAssessment(
                        protein_id=protein_id,
                        domain_authority=authority,
                        assessment_status=DomainAssessmentStatus.NOT_ASSESSED,
                        hit_count=0,
                        evidence_source="",
                        evidence_reference="",
                    ),
                )
            )
    return tuple(result)


def derive_domain_evidence(
    *, hits: tuple[DomainHit, ...], sequences: tuple[SequenceRecord, ...]
) -> tuple[tuple[FeatureRecord, ...], tuple[dict[str, object], ...]]:
    """Derive domain-presence, architecture and coordinate-sequence evidence.

    Args:
        hits: Validated coordinate-resolved domain hits.
        sequences: Authoritative protein inventory.

    Returns:
        Derived categorical features and extracted domain-sequence rows.
    """

    sequence_by_id = {item.protein_id: item.sequence for item in sequences}
    by_protein_authority: dict[tuple[str, str], list[DomainHit]] = defaultdict(list)
    features: list[FeatureRecord] = []
    domain_sequences: list[dict[str, object]] = []
    for hit in hits:
        by_protein_authority[(hit.protein_id, hit.domain_authority)].append(hit)
        feature_type = "PFAM_DOMAIN" if hit.domain_authority.casefold() == "pfam" else "DOMAIN"
        features.append(
            FeatureRecord(
                protein_id=hit.protein_id,
                feature_type=feature_type,
                feature_id=f"{hit.domain_authority}:{hit.domain_id}",
                feature_name=f"{hit.domain_authority} {hit.domain_name}",
                start=hit.start,
                end=hit.end,
                evidence_status="ASSESSED_WITH_HIT",
                evidence_source=hit.evidence_source,
                evidence_reference=hit.evidence_reference,
            )
        )
        domain_sequences.append(
            {
                "protein_id": hit.protein_id,
                "domain_authority": hit.domain_authority,
                "domain_id": hit.domain_id,
                "start": hit.start,
                "end": hit.end,
                "domain_sequence": sequence_by_id[hit.protein_id][hit.start - 1 : hit.end],
            }
        )
    for (protein_id, authority), ordered_hits in sorted(by_protein_authority.items()):
        ordered_hits.sort(key=lambda item: (item.start, item.end, item.domain_id))
        architecture = ">".join(item.domain_id for item in ordered_hits)
        feature_type = (
            "PFAM_ARCHITECTURE" if authority.casefold() == "pfam" else "DOMAIN_ARCHITECTURE"
        )
        features.append(
            FeatureRecord(
                protein_id=protein_id,
                feature_type=feature_type,
                feature_id=f"{authority}:{architecture}",
                feature_name=f"{authority} architecture {architecture}",
                start=None,
                end=None,
                evidence_status="DERIVED",
                evidence_source="protein-signature-analysis",
                evidence_reference="ordered_coordinate_architecture",
            )
        )
    ordered_features = tuple(
        sorted(
            features,
            key=lambda item: (
                item.protein_id,
                item.feature_type,
                item.feature_id,
                item.start or 0,
            ),
        )
    )
    ordered_sequences = tuple(
        sorted(
            domain_sequences,
            key=lambda item: (
                str(item["protein_id"]),
                str(item["domain_authority"]),
                int(item["start"]),
            ),
        )
    )
    return ordered_features, ordered_sequences


def read_structures(
    *, path: Path, sequences: tuple[SequenceRecord, ...]
) -> tuple[StructureRecord, ...]:
    """Read structure inventory rows and verify local coordinate checksums.

    Args:
        path: Structure TSV.
        sequences: Authoritative protein inventory.

    Returns:
        Validated structure records.
    """

    protein_ids = frozenset(item.protein_id for item in sequences)
    structures: list[StructureRecord] = []
    observed: set[str] = set()
    fold_versions_by_authority: dict[str, str] = {}
    for row in iter_tsv(path=path, required_fields=STRUCTURE_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        _validate_known_protein(protein_id=protein_id, protein_ids=protein_ids)
        structure_id = validate_identifier(value=row["structure_id"], field_name="structure_id")
        if structure_id in observed:
            raise InputValidationError(f"Duplicate structure_id: {structure_id!r}")
        observed.add(structure_id)
        coordinate_path = _resolve_coordinate_path(value=row["coordinate_path"], table_path=path)
        supplied_digest = row["coordinate_sha256"].strip().lower()
        if coordinate_path is not None:
            actual_digest = sha256_file(path=coordinate_path)
            if supplied_digest and supplied_digest != actual_digest:
                raise InputValidationError(
                    f"Coordinate checksum mismatch for structure {structure_id!r}."
                )
            supplied_digest = actual_digest
        elif supplied_digest:
            raise InputValidationError(
                f"coordinate_sha256 supplied without coordinate_path for {structure_id!r}."
            )
        availability_status = validate_identifier(
            value=row["availability_status"], field_name="availability_status"
        ).upper()
        mean_confidence = parse_optional_float(
            value=row["mean_confidence"],
            field_name="mean_confidence",
            minimum=0.0,
            maximum=100.0,
        )
        analysis_eligibility = _parse_structure_analysis_eligibility(
            value=row["analysis_eligibility_status"]
        )
        _validate_structure_state(
            structure_id=structure_id,
            coordinate_path=coordinate_path,
            availability_status=availability_status,
            mean_confidence=mean_confidence,
            analysis_eligibility=analysis_eligibility,
        )
        fold_id = validate_text(value=row["fold_id"], field_name="fold_id", allow_empty=True)
        fold_name = validate_text(value=row["fold_name"], field_name="fold_name", allow_empty=True)
        fold_authority = validate_text(
            value=row["fold_authority"],
            field_name="fold_authority",
            allow_empty=True,
        )
        fold_authority_version = validate_text(
            value=row["fold_authority_version"],
            field_name="fold_authority_version",
            allow_empty=True,
        )
        fold_evidence_reference = validate_text(
            value=row["fold_evidence_reference"],
            field_name="fold_evidence_reference",
            allow_empty=True,
        )
        fold_evidence_status = _parse_fold_evidence_status(value=row["fold_evidence_status"])
        comparison_universe_ids = _parse_comparison_universe_ids(
            value=row["comparison_universe_ids"]
        )
        _validate_fold_evidence(
            structure_id=structure_id,
            fold_id=fold_id,
            fold_name=fold_name,
            fold_authority=fold_authority,
            fold_authority_version=fold_authority_version,
            fold_evidence_reference=fold_evidence_reference,
            status=fold_evidence_status,
        )
        if fold_authority:
            previous_version = fold_versions_by_authority.setdefault(
                fold_authority, fold_authority_version
            )
            if previous_version != fold_authority_version:
                raise InputValidationError(
                    f"Fold authority {fold_authority!r} uses multiple releases in one "
                    f"campaign: {previous_version!r} and {fold_authority_version!r}."
                )
        structures.append(
            StructureRecord(
                protein_id=protein_id,
                structure_id=structure_id,
                structure_source=validate_text(
                    value=row["structure_source"], field_name="structure_source"
                ),
                structure_version=validate_text(
                    value=row["structure_version"],
                    field_name="structure_version",
                    allow_empty=True,
                ),
                coordinate_path=coordinate_path,
                coordinate_sha256=supplied_digest,
                availability_status=availability_status,
                mean_confidence=mean_confidence,
                fold_id=fold_id,
                fold_name=fold_name,
                fold_authority=fold_authority,
                fold_authority_version=fold_authority_version,
                fold_evidence_reference=fold_evidence_reference,
                fold_evidence_status=fold_evidence_status,
                analysis_eligibility_status=analysis_eligibility,
                comparison_universe_ids=comparison_universe_ids,
            )
        )
    result = tuple(sorted(structures, key=lambda item: item.structure_id))
    LOGGER.info("Loaded %d structure inventory records", len(result))
    return result


def read_structure_comparisons(
    *, path: Path, protein_ids: frozenset[str]
) -> tuple[PairwiseStructureComparison, ...]:
    """Read validated pairwise structure-alignment evidence.

    Args:
        path: Pairwise comparison TSV.
        protein_ids: Authoritative FASTA identifiers.

    Returns:
        Deterministically ordered comparison rows.
    """

    comparisons: list[PairwiseStructureComparison] = []
    observed: set[tuple[str, str, str, str]] = set()
    for row in iter_tsv(path=path, required_fields=STRUCTURE_COMPARISON_FIELDS):
        protein_a = validate_identifier(value=row["protein_a_id"], field_name="protein_a_id")
        protein_b = validate_identifier(value=row["protein_b_id"], field_name="protein_b_id")
        _validate_known_protein(protein_id=protein_a, protein_ids=protein_ids)
        _validate_known_protein(protein_id=protein_b, protein_ids=protein_ids)
        if protein_a == protein_b:
            raise InputValidationError(f"Self structure comparison is not allowed: {protein_a!r}")
        tool = validate_identifier(value=row["comparison_tool"], field_name="comparison_tool")
        record_id = validate_identifier(
            value=row["source_record_id"], field_name="source_record_id"
        )
        universe_id = validate_identifier(
            value=row["comparison_universe_id"],
            field_name="comparison_universe_id",
        )
        coverage_scope = _parse_structure_coverage_scope(value=row["coverage_scope"])
        pair = tuple(sorted((protein_a, protein_b)))
        key = (pair[0], pair[1], tool, record_id)
        if key in observed:
            raise InputValidationError(f"Duplicate structure comparison: {key!r}")
        observed.add(key)
        raw_coverage_a = parse_optional_float(
            value=row["coverage_a"],
            field_name="coverage_a",
            minimum=0.0,
            maximum=1.0,
        )
        raw_coverage_b = parse_optional_float(
            value=row["coverage_b"],
            field_name="coverage_b",
            minimum=0.0,
            maximum=1.0,
        )
        coverage_a, coverage_b = (
            (raw_coverage_a, raw_coverage_b)
            if protein_a == pair[0]
            else (raw_coverage_b, raw_coverage_a)
        )
        comparison_status = _parse_structure_comparison_status(value=row["comparison_status"])
        tm_score = parse_optional_float(
            value=row["tm_score"], field_name="tm_score", minimum=0.0, maximum=1.0
        )
        rmsd = parse_optional_float(
            value=row["rmsd_angstrom"], field_name="rmsd_angstrom", minimum=0.0
        )
        aligned_residue_count = parse_optional_integer(
            value=row["aligned_residue_count"],
            field_name="aligned_residue_count",
            minimum=1,
        )
        _validate_structure_comparison_payload(
            source_record_id=record_id,
            status=comparison_status,
            tm_score=tm_score,
            rmsd_angstrom=rmsd,
            aligned_residue_count=aligned_residue_count,
            coverage_a=coverage_a,
            coverage_b=coverage_b,
        )
        comparisons.append(
            PairwiseStructureComparison(
                protein_a_id=pair[0],
                protein_b_id=pair[1],
                comparison_tool=tool,
                comparison_tool_version=validate_text(
                    value=row["comparison_tool_version"],
                    field_name="comparison_tool_version",
                    allow_empty=True,
                ),
                tm_score=tm_score,
                rmsd_angstrom=rmsd,
                aligned_residue_count=aligned_residue_count,
                coverage_a=coverage_a,
                coverage_b=coverage_b,
                comparison_status=comparison_status,
                source_record_id=record_id,
                comparison_universe_id=universe_id,
                coverage_scope=coverage_scope,
            )
        )
    result = tuple(
        sorted(
            comparisons,
            key=lambda item: (
                item.protein_a_id,
                item.protein_b_id,
                item.comparison_tool,
                item.source_record_id,
            ),
        )
    )
    LOGGER.info("Loaded %d pairwise structure comparisons", len(result))
    return result


def _parse_structure_analysis_eligibility(*, value: str) -> StructureAnalysisEligibility:
    """Parse one controlled local-coordinate eligibility state.

    Args:
        value: Raw eligibility value.

    Returns:
        Controlled eligibility enum member.

    Raises:
        InputValidationError: If the state is unsupported.
    """

    normalised = validate_identifier(
        value=value,
        field_name="analysis_eligibility_status",
    ).upper()
    try:
        return StructureAnalysisEligibility(normalised)
    except ValueError as error:
        allowed = ", ".join(item.value for item in StructureAnalysisEligibility)
        raise InputValidationError(
            f"analysis_eligibility_status must be one of {allowed}; received {normalised!r}."
        ) from error


def _parse_fold_evidence_status(*, value: str) -> FoldEvidenceStatus:
    """Parse one controlled fold-assessment state.

    Args:
        value: Raw fold-assessment value.

    Returns:
        Controlled fold-evidence enum member.

    Raises:
        InputValidationError: If the state is unsupported.
    """

    normalised = validate_identifier(value=value, field_name="fold_evidence_status").upper()
    try:
        return FoldEvidenceStatus(normalised)
    except ValueError as error:
        allowed = ", ".join(item.value for item in FoldEvidenceStatus)
        raise InputValidationError(
            f"fold_evidence_status must be one of {allowed}; received {normalised!r}."
        ) from error


def _parse_structure_comparison_status(*, value: str) -> StructureComparisonStatus:
    """Parse one controlled pairwise-comparison state.

    Args:
        value: Raw comparison state.

    Returns:
        Controlled comparison-status enum member.

    Raises:
        InputValidationError: If the state is unsupported.
    """

    normalised = validate_identifier(value=value, field_name="comparison_status").upper()
    try:
        return StructureComparisonStatus(normalised)
    except ValueError as error:
        allowed = ", ".join(item.value for item in StructureComparisonStatus)
        raise InputValidationError(
            f"comparison_status must be one of {allowed}; received {normalised!r}."
        ) from error


def _parse_structure_coverage_scope(*, value: str) -> StructureCoverageScope:
    """Parse one controlled structural-coverage denominator.

    Args:
        value: Raw coverage-scope value.

    Returns:
        Controlled coverage-scope enum member.

    Raises:
        InputValidationError: If the scope is unsupported.
    """

    normalised = validate_identifier(value=value, field_name="coverage_scope").upper()
    try:
        return StructureCoverageScope(normalised)
    except ValueError as error:
        allowed = ", ".join(item.value for item in StructureCoverageScope)
        raise InputValidationError(
            f"coverage_scope must be one of {allowed}; received {normalised!r}."
        ) from error


def _parse_comparison_universe_ids(*, value: str) -> tuple[str, ...]:
    """Parse pipe-separated comparison-universe memberships.

    Args:
        value: Raw membership list; blank means no declared membership.

    Returns:
        Sorted unique comparison-universe identifiers.

    Raises:
        InputValidationError: If an identifier is empty, unsafe or repeated.
    """

    text = str(value or "").strip()
    if not text:
        return ()
    raw_identifiers = text.split("|")
    if any(not item.strip() for item in raw_identifiers):
        raise InputValidationError(
            "comparison_universe_ids contains an empty pipe-separated identifier."
        )
    identifiers = tuple(
        validate_identifier(value=item, field_name="comparison_universe_ids")
        for item in raw_identifiers
    )
    if len(set(identifiers)) != len(identifiers):
        raise InputValidationError("comparison_universe_ids contains a duplicate identifier.")
    return tuple(sorted(identifiers))


def _validate_structure_state(
    *,
    structure_id: str,
    coordinate_path: Path | None,
    availability_status: str,
    mean_confidence: float | None,
    analysis_eligibility: StructureAnalysisEligibility,
) -> None:
    """Require internally consistent structure availability and eligibility.

    Args:
        structure_id: Structure identifier for diagnostics.
        coordinate_path: Resolved coordinate file, when physically available.
        availability_status: Controlled physical-availability state.
        mean_confidence: Optional source-specific mean confidence.
        analysis_eligibility: Declared local-coordinate eligibility.

    Raises:
        InputValidationError: If availability, coordinates or eligibility conflict.
    """

    if availability_status not in _STRUCTURE_AVAILABILITY_STATUSES:
        allowed = ", ".join(sorted(_STRUCTURE_AVAILABILITY_STATUSES))
        raise InputValidationError(
            f"Unsupported availability_status for {structure_id!r}; expected one of {allowed}."
        )
    coordinate_available = coordinate_path is not None
    status_claims_available = availability_status in _AVAILABLE_COORDINATE_STATUSES
    if coordinate_available != status_claims_available:
        raise InputValidationError(
            f"Structure {structure_id!r} has inconsistent coordinate_path and "
            f"availability_status {availability_status!r}."
        )
    if analysis_eligibility == StructureAnalysisEligibility.ELIGIBLE and not coordinate_available:
        raise InputValidationError(
            f"Structure {structure_id!r} cannot be ELIGIBLE without a coordinate file."
        )
    if (
        analysis_eligibility == StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE
        and coordinate_available
    ):
        raise InputValidationError(
            f"Structure {structure_id!r} uses external-evidence eligibility despite "
            "having local coordinates."
        )
    if (
        analysis_eligibility == StructureAnalysisEligibility.INELIGIBLE_COORDINATE_UNAVAILABLE
        and coordinate_available
    ):
        raise InputValidationError(
            f"Structure {structure_id!r} declares unavailable coordinates but has a file."
        )
    if (
        analysis_eligibility == StructureAnalysisEligibility.INELIGIBLE_LOW_CONFIDENCE
        and mean_confidence is None
    ):
        raise InputValidationError(
            f"Structure {structure_id!r} declares low confidence without mean_confidence."
        )
    if (
        analysis_eligibility == StructureAnalysisEligibility.INELIGIBLE_CONFIDENCE_UNAVAILABLE
        and mean_confidence is not None
    ):
        raise InputValidationError(
            f"Structure {structure_id!r} declares confidence unavailable but supplies a value."
        )


def _validate_fold_evidence(
    *,
    structure_id: str,
    fold_id: str,
    fold_name: str,
    fold_authority: str,
    fold_authority_version: str,
    fold_evidence_reference: str,
    status: FoldEvidenceStatus,
) -> None:
    """Require fold payload to agree with its controlled assessment state.

    Args:
        structure_id: Structure identifier for diagnostics.
        fold_id: Optional authority-specific fold identifier.
        fold_name: Optional human-readable fold name.
        fold_authority: Fold classification authority.
        fold_authority_version: Exact authority release or database version.
        fold_evidence_reference: Reproducible annotation record or method reference.
        status: Controlled fold-assessment state.

    Raises:
        InputValidationError: If a hit lacks identifiers or a non-hit carries them.
    """

    if status == FoldEvidenceStatus.ASSESSED_WITH_HIT:
        if not all(
            (
                fold_id,
                fold_authority,
                fold_authority_version,
                fold_evidence_reference,
            )
        ):
            raise InputValidationError(
                f"Fold hit for structure {structure_id!r} requires fold_id, "
                "fold_authority, fold_authority_version and fold_evidence_reference."
            )
        return
    if fold_id or fold_name:
        raise InputValidationError(
            f"Non-hit fold status {status.value!r} for structure {structure_id!r} "
            "must not carry fold_id or fold_name."
        )
    assessed_without_hit = {
        FoldEvidenceStatus.ASSESSED_NO_HIT,
        FoldEvidenceStatus.FAILED,
    }
    if status in assessed_without_hit and not all(
        (fold_authority, fold_authority_version, fold_evidence_reference)
    ):
        raise InputValidationError(
            f"Fold status {status.value!r} for structure {structure_id!r} requires "
            "fold_authority, fold_authority_version and fold_evidence_reference."
        )
    if status in {FoldEvidenceStatus.NOT_ASSESSED, FoldEvidenceStatus.INPUT_UNAVAILABLE} and any(
        (fold_authority, fold_authority_version, fold_evidence_reference)
    ):
        raise InputValidationError(
            f"Fold status {status.value!r} for structure {structure_id!r} must not carry "
            "fold authority or evidence provenance."
        )


def _validate_structure_comparison_payload(
    *,
    source_record_id: str,
    status: StructureComparisonStatus,
    tm_score: float | None,
    rmsd_angstrom: float | None,
    aligned_residue_count: int | None,
    coverage_a: float | None,
    coverage_b: float | None,
) -> None:
    """Require pairwise metrics to agree with comparison completion state.

    Args:
        source_record_id: Source identifier for diagnostics.
        status: Controlled comparison completion state.
        tm_score: Optional conservative TM score.
        rmsd_angstrom: Optional RMSD.
        aligned_residue_count: Optional aligned length.
        coverage_a: Optional first-protein coverage.
        coverage_b: Optional second-protein coverage.

    Raises:
        InputValidationError: If completed evidence lacks clustering metrics or an
            unsuccessful row carries numerical evidence.
    """

    metrics = (tm_score, rmsd_angstrom, aligned_residue_count, coverage_a, coverage_b)
    if status in _SUCCESSFUL_STRUCTURE_COMPARISON_STATUSES:
        if tm_score is None or coverage_a is None or coverage_b is None:
            raise InputValidationError(
                f"Successful structure comparison {source_record_id!r} requires tm_score "
                "and bilateral coverage."
            )
        return
    if any(value is not None for value in metrics):
        raise InputValidationError(
            f"Unsuccessful structure comparison {source_record_id!r} must not carry metrics."
        )


def _parse_coordinates(
    *, start: str, end: str, protein_id: str, length: int
) -> tuple[int | None, int | None]:
    """Parse paired, one-based inclusive feature coordinates.

    Args:
        start: Raw start value.
        end: Raw end value.
        protein_id: Protein identifier for diagnostics.
        length: Protein sequence length.

    Returns:
        Parsed coordinate pair or two ``None`` values.
    """

    parsed_start = parse_optional_integer(value=start, field_name="start", minimum=1)
    parsed_end = parse_optional_integer(value=end, field_name="end", minimum=1)
    if (parsed_start is None) != (parsed_end is None):
        raise InputValidationError(f"Feature coordinates must be paired for {protein_id!r}.")
    if parsed_start is not None and parsed_end is not None:
        if parsed_start > parsed_end or parsed_end > length:
            raise InputValidationError(
                f"Invalid coordinates {parsed_start}-{parsed_end} for {protein_id!r} "
                f"of length {length}."
            )
    return parsed_start, parsed_end


def _validate_known_protein(*, protein_id: str, protein_ids: frozenset[str]) -> None:
    """Require a protein identifier to exist in the FASTA authority.

    Args:
        protein_id: Candidate identifier.
        protein_ids: Authoritative identifiers.

    Raises:
        InputValidationError: If the protein is unknown.
    """

    if protein_id not in protein_ids:
        raise InputValidationError(f"Input references protein absent from FASTA: {protein_id!r}")


def _reject_nonblank_domain_payload(*, row: dict[str, str], key: tuple[str, str]) -> None:
    """Require sentinel domain rows to leave hit-specific fields blank.

    Args:
        row: Raw TSV row.
        key: Protein and authority diagnostic key.

    Raises:
        InputValidationError: If a no-hit state carries hit data.
    """

    fields = ("domain_name", "start", "end", "score", "e_value")
    if any(row[field].strip() for field in fields):
        raise InputValidationError(f"Non-hit domain row contains hit payload for {key!r}")


def _resolve_coordinate_path(*, value: str, table_path: Path) -> Path | None:
    """Resolve an optional structure coordinate path relative to its TSV.

    Args:
        value: Raw path text.
        table_path: Structure inventory path.

    Returns:
        Existing coordinate file or ``None``.
    """

    if not value.strip():
        return None
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = Path(table_path).expanduser().resolve().parent / candidate
    candidate = candidate.resolve()
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty structure coordinate file: {candidate}")
    return candidate
