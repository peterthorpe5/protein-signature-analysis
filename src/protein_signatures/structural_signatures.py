"""Fold and pairwise-alignment feature derivation for structural signatures."""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from collections.abc import Mapping

from .checksums import sha256_json
from .errors import InputValidationError
from .feature_provenance import ALL_DATA_EXPLORATORY, DISCOVERY_DERIVED, FIXED_EXTERNAL
from .models import (
    AlphaFoldAcquisition,
    FeatureRecord,
    FoldEvidenceStatus,
    PairwiseStructureComparison,
    PartitionAssignment,
    StructureAnalysisEligibility,
    StructureComparisonStatus,
    StructureRecord,
)

LOGGER = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}")


def derive_structure_features(
    *,
    structures: tuple[StructureRecord, ...],
    comparisons: tuple[PairwiseStructureComparison, ...],
    tm_score_threshold: float,
    minimum_coverage: float,
    partitions: tuple[PartitionAssignment, ...] = (),
    derivation_cohort_sha256: str = "",
) -> tuple[tuple[FeatureRecord, ...], tuple[dict[str, object], ...]]:
    """Derive fold, model-availability and alignment-cluster features.

    Args:
        structures: Structure inventory including optional fold assignments.
        comparisons: Pairwise structural-alignment evidence.
        tm_score_threshold: Inclusive TM-score edge threshold.
        minimum_coverage: Inclusive coverage required for both proteins.
        partitions: Optional leakage-safe discovery/validation assignments. When supplied,
            structural clusters are defined only from discovery-to-discovery edges and
            validation proteins are projected onto those frozen components.
        derivation_cohort_sha256: Exact sequence-cohort digest for discovery, or
            for all proteins when no partitions are supplied. Required whenever
            passing comparisons could define a structural cluster.

    Returns:
        Derived protein features and transparent structural-cluster membership rows.

    Raises:
        InputValidationError: If thresholds or cluster derivation provenance are invalid.
    """

    if not 0.0 <= tm_score_threshold <= 1.0 or not 0.0 <= minimum_coverage <= 1.0:
        raise InputValidationError("Structural score and coverage thresholds must be 0.0 to 1.0.")
    features: list[FeatureRecord] = []
    fold_versions_by_authority: dict[str, str] = {}
    fold_definitions: dict[tuple[str, str], tuple[str, str]] = {}
    for structure in structures:
        if structure.is_coordinate_analysis_eligible:
            features.append(
                _feature(
                    protein_id=structure.protein_id,
                    feature_type="STRUCTURE_AVAILABLE",
                    feature_id=structure.structure_source,
                    feature_name=f"Structure available from {structure.structure_source}",
                    reference=structure.structure_id,
                    derivation_scope=FIXED_EXTERNAL,
                    feature_definition_sha256=sha256_json(
                        value={
                            "feature_type": "STRUCTURE_AVAILABLE",
                            "feature_id": structure.structure_source,
                            "rule": "EXPLICITLY_ELIGIBLE_AVAILABLE_COORDINATE",
                            "rule_version": 1,
                        }
                    ),
                )
            )
        fold_eligible_states = {
            StructureAnalysisEligibility.ELIGIBLE,
            StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE,
        }
        if (
            structure.analysis_eligibility_status in fold_eligible_states
            and structure.fold_evidence_status == FoldEvidenceStatus.ASSESSED_WITH_HIT
        ):
            if not all(
                (
                    structure.fold_id,
                    structure.fold_authority,
                    structure.fold_authority_version,
                    structure.fold_evidence_reference,
                )
            ):
                raise InputValidationError(
                    f"Fold hit for structure {structure.structure_id!r} lacks exact "
                    "authority release or evidence provenance."
                )
            previous_version = fold_versions_by_authority.setdefault(
                structure.fold_authority, structure.fold_authority_version
            )
            if previous_version != structure.fold_authority_version:
                raise InputValidationError(
                    f"Fold authority {structure.fold_authority!r} uses multiple releases "
                    "within one structural feature derivation."
                )
            fold_key = (structure.fold_authority, structure.fold_id)
            fold_definition = (
                structure.fold_authority_version,
                structure.fold_evidence_reference,
            )
            previous_definition = fold_definitions.setdefault(fold_key, fold_definition)
            if previous_definition != fold_definition:
                raise InputValidationError(
                    f"Fold feature {fold_key!r} has conflicting authority-release or "
                    "evidence-reference definitions."
                )
            features.append(
                _feature(
                    protein_id=structure.protein_id,
                    feature_type="FOLD",
                    feature_id=f"{structure.fold_authority}:{structure.fold_id}",
                    feature_name=(
                        f"{structure.fold_authority} {structure.fold_name or structure.fold_id}"
                    ),
                    reference=structure.fold_evidence_reference,
                    evidence_source=structure.fold_authority,
                    derivation_scope=FIXED_EXTERNAL,
                    feature_definition_sha256=sha256_json(
                        value={
                            "feature_type": "FOLD",
                            "feature_id": f"{structure.fold_authority}:{structure.fold_id}",
                            "fold_authority": structure.fold_authority,
                            "fold_authority_version": structure.fold_authority_version,
                            "fold_evidence_reference": structure.fold_evidence_reference,
                            "rule": "ASSESSED_WITH_HIT_BY_DECLARED_AUTHORITY",
                            "rule_version": 2,
                        }
                    ),
                )
            )
    passing = tuple(
        comparison
        for comparison in comparisons
        if _comparison_passes(
            comparison=comparison,
            tm_score_threshold=tm_score_threshold,
            minimum_coverage=minimum_coverage,
        )
    )
    if passing and _DIGEST.fullmatch(derivation_cohort_sha256) is None:
        scope_name = "discovery" if partitions else "all-data"
        raise InputValidationError(
            f"Structural cluster derivation requires the exact {scope_name} "
            "sequence-cohort SHA-256."
        )
    cluster_derivation_scope = DISCOVERY_DERIVED if partitions else ALL_DATA_EXPLORATORY
    partition_by_protein = {item.protein_id: item.partition for item in partitions}
    reference_partition = "DISCOVERY" if partitions else "ALL_DATA"
    comparisons_by_universe: dict[tuple[str, str, str, str], list[PairwiseStructureComparison]] = (
        defaultdict(list)
    )
    for comparison in passing:
        key = (
            comparison.comparison_universe_id,
            comparison.coverage_scope.value,
            comparison.comparison_tool,
            comparison.comparison_tool_version,
        )
        comparisons_by_universe[key].append(comparison)
    cluster_rows: list[dict[str, object]] = []
    for universe_key in sorted(comparisons_by_universe):
        universe_comparisons = tuple(comparisons_by_universe[universe_key])
        adjacency: dict[str, set[str]] = defaultdict(set)
        for comparison in universe_comparisons:
            if partitions and not (
                partition_by_protein.get(comparison.protein_a_id) == "DISCOVERY"
                and partition_by_protein.get(comparison.protein_b_id) == "DISCOVERY"
            ):
                continue
            adjacency[comparison.protein_a_id].add(comparison.protein_b_id)
            adjacency[comparison.protein_b_id].add(comparison.protein_a_id)
        components = _connected_components(adjacency=adjacency)
        universe_id, coverage_scope, comparison_tool, comparison_tool_version = universe_key
        evidence_reference = "|".join(universe_key)
        for component in components:
            if len(component) < 2:
                continue
            definition_digest = sha256_json(
                value={
                    "feature_type": "STRUCTURE_CLUSTER",
                    "comparison_universe_id": universe_id,
                    "coverage_scope": coverage_scope,
                    "comparison_tool": comparison_tool,
                    "comparison_tool_version": comparison_tool_version,
                    "reference_members": component,
                    "reference_partition": reference_partition,
                    "tm_score_threshold": tm_score_threshold,
                    "minimum_coverage": minimum_coverage,
                    "clustering_method": "CONNECTED_COMPONENTS_SINGLE_LINKAGE",
                }
            )
            cluster_id = f"SC_{definition_digest[:16]}"
            member_set = set(component)
            number_of_edges = sum(
                1
                for source in component
                for target in adjacency[source]
                if target in member_set and source < target
            )
            members: dict[str, tuple[str, int, float | None]] = {}
            for protein_id in component:
                incident = [
                    item.tm_score
                    for item in universe_comparisons
                    if item.tm_score is not None
                    and {item.protein_a_id, item.protein_b_id} <= member_set
                    and protein_id in {item.protein_a_id, item.protein_b_id}
                ]
                method = "DISCOVERY_COMPONENT" if partitions else "ALL_DATA_COMPONENT"
                members[protein_id] = (method, len(incident), max(incident) if incident else None)
            if partitions:
                for protein_id, partition in sorted(partition_by_protein.items()):
                    if partition != "VALIDATION":
                        continue
                    incident = [
                        item.tm_score
                        for item in universe_comparisons
                        if item.tm_score is not None
                        and (
                            item.protein_a_id == protein_id
                            and item.protein_b_id in member_set
                            or item.protein_b_id == protein_id
                            and item.protein_a_id in member_set
                        )
                    ]
                    if incident:
                        members[protein_id] = (
                            "VALIDATION_PROJECTION",
                            len(incident),
                            max(incident),
                        )
            for protein_id, (membership_method, support_count, best_score) in sorted(
                members.items()
            ):
                features.append(
                    _feature(
                        protein_id=protein_id,
                        feature_type="STRUCTURE_CLUSTER",
                        feature_id=cluster_id,
                        feature_name=f"{comparison_tool} alignment cluster {cluster_id}",
                        reference=evidence_reference,
                        derivation_scope=cluster_derivation_scope,
                        feature_definition_sha256=definition_digest,
                        derivation_cohort_sha256=derivation_cohort_sha256,
                    )
                )
                cluster_rows.append(
                    {
                        "cluster_id": cluster_id,
                        "protein_id": protein_id,
                        "member_count": len(members),
                        "reference_member_count": len(component),
                        "edge_count": number_of_edges,
                        "reference_partition": reference_partition,
                        "membership_method": membership_method,
                        "supporting_edge_count": support_count,
                        "best_tm_score": best_score,
                        "tm_score_threshold": tm_score_threshold,
                        "minimum_coverage": minimum_coverage,
                        "clustering_method": "CONNECTED_COMPONENTS_SINGLE_LINKAGE",
                        "comparison_universe_id": universe_id,
                        "coverage_scope": coverage_scope,
                        "comparison_tool": comparison_tool,
                        "comparison_tool_version": comparison_tool_version,
                    }
                )
    unique_features = {
        (item.protein_id, item.feature_type, item.feature_id, item.evidence_reference): item
        for item in features
    }
    ordered_features = tuple(
        unique_features[key]
        for key in sorted(unique_features, key=lambda item: (item[0], item[1], item[2], item[3]))
    )
    ordered_clusters = tuple(
        sorted(cluster_rows, key=lambda item: (str(item["cluster_id"]), str(item["protein_id"])))
    )
    return ordered_features, ordered_clusters


def _comparison_passes(
    *,
    comparison: PairwiseStructureComparison,
    tm_score_threshold: float,
    minimum_coverage: float,
) -> bool:
    """Return whether a pairwise structural comparison forms a cluster edge.

    Args:
        comparison: Validated comparison.
        tm_score_threshold: Inclusive TM-score threshold.
        minimum_coverage: Inclusive bilateral coverage threshold.

    Returns:
        ``True`` only for complete, numerically sufficient evidence.
    """

    return (
        comparison.comparison_status
        in {
            StructureComparisonStatus.COMPLETE,
            StructureComparisonStatus.SUCCESS,
            StructureComparisonStatus.PASS,
        }
        and comparison.tm_score is not None
        and comparison.coverage_a is not None
        and comparison.coverage_b is not None
        and comparison.tm_score >= tm_score_threshold
        and comparison.coverage_a >= minimum_coverage
        and comparison.coverage_b >= minimum_coverage
    )


def derive_structure_assessment_universes(
    *,
    structures: tuple[StructureRecord, ...],
    structure_clusters: tuple[dict[str, object], ...],
    comparison_universe_members: Mapping[str, frozenset[str]] | None = None,
    alphafold_acquisitions: tuple[AlphaFoldAcquisition, ...] = (),
) -> dict[tuple[str, str], frozenset[str]]:
    """Build only scientifically proved structural assessment universes.

    Absence is meaningful only among proteins that were successfully assessed for
    the same feature definition. Fold universes therefore use explicit hit/no-hit
    rows from the same authority. Structure-availability universes use terminal
    source-specific model eligibility assessments. Structural-cluster universes are
    emitted only when the caller supplies the complete membership of the declared
    comparison universe; endpoints observed in retained hits are not sufficient.

    Args:
        structures: Validated structure and fold-assessment records, including any
            explicit supplied comparison-universe memberships.
        structure_clusters: Derived cluster membership rows containing cluster and
            comparison-universe identifiers.
        comparison_universe_members: Optional complete protein membership keyed by
            comparison-universe identifier, such as Foldseek's eligible query set.
        alphafold_acquisitions: Complete AlphaFold request outcomes used to retain
            known model absence and sequence/quality exclusions in the denominator.

    Returns:
        Feature keys mapped to proteins with known presence or qualifying absence.
        Missing keys deliberately mean that the assessment universe is unproved.

    Raises:
        InputValidationError: If a supplied universe is empty or malformed.
    """

    declared_members: dict[str, set[str]] = defaultdict(set)
    for structure in structures:
        for universe_id in structure.comparison_universe_ids:
            declared_members[universe_id].add(structure.protein_id)
    known_universes: dict[str, frozenset[str]] = {
        universe_id: frozenset(members) for universe_id, members in declared_members.items()
    }
    for universe_id, members in (comparison_universe_members or {}).items():
        if not isinstance(universe_id, str) or not universe_id.strip():
            raise InputValidationError("Comparison-universe identifiers must be non-empty text.")
        if (
            not isinstance(members, frozenset)
            or not members
            or any(
                not isinstance(protein_id, str) or not protein_id.strip() for protein_id in members
            )
        ):
            raise InputValidationError(
                f"Comparison universe {universe_id!r} requires a non-empty frozenset "
                "of protein identifiers."
            )
        declared = known_universes.get(universe_id)
        if declared is not None and declared != members:
            raise InputValidationError(
                f"Comparison universe {universe_id!r} conflicts between structure "
                "membership and completed-search provenance."
            )
        known_universes[universe_id] = members

    result: dict[tuple[str, str], frozenset[str]] = {}
    terminal_availability = {
        StructureAnalysisEligibility.ELIGIBLE,
        StructureAnalysisEligibility.INELIGIBLE_CONFIDENCE_UNAVAILABLE,
        StructureAnalysisEligibility.INELIGIBLE_COORDINATE_UNAVAILABLE,
        StructureAnalysisEligibility.INELIGIBLE_LOW_CONFIDENCE,
        StructureAnalysisEligibility.INELIGIBLE_SEQUENCE_MISMATCH,
        StructureAnalysisEligibility.INELIGIBLE_SEQUENCE_UNVERIFIED,
        StructureAnalysisEligibility.INELIGIBLE_USER_EXCLUDED,
    }
    sources = sorted(
        {
            structure.structure_source
            for structure in structures
            if structure.analysis_eligibility_status in terminal_availability
        }
    )
    for source in sources:
        assessed = frozenset(
            structure.protein_id
            for structure in structures
            if structure.structure_source == source
            and structure.analysis_eligibility_status in terminal_availability
        )
        if assessed:
            result[("STRUCTURE_AVAILABLE", source)] = assessed
    terminal_alphafold_statuses = {
        "ACQUIRED",
        "ACQUIRED_CONFIDENCE_UNAVAILABLE",
        "ACQUIRED_LOW_CONFIDENCE",
        "ACQUIRED_SEQUENCE_UNVERIFIED",
        "MODEL_NOT_AVAILABLE",
        "SEQUENCE_MISMATCH",
    }
    alphafold_assessed = frozenset(
        outcome.protein_id
        for outcome in alphafold_acquisitions
        if outcome.acquisition_status in terminal_alphafold_statuses
    )
    if alphafold_assessed:
        key = ("STRUCTURE_AVAILABLE", "AlphaFoldDB")
        result[key] = result.get(key, frozenset()) | alphafold_assessed

    fold_eligible_states = {
        StructureAnalysisEligibility.ELIGIBLE,
        StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE,
    }
    assessed_fold_statuses = {
        FoldEvidenceStatus.ASSESSED_NO_HIT,
        FoldEvidenceStatus.ASSESSED_WITH_HIT,
    }
    positive_folds = {
        (
            structure.fold_authority,
            structure.fold_authority_version,
            structure.fold_id,
        )
        for structure in structures
        if structure.analysis_eligibility_status in fold_eligible_states
        and structure.fold_evidence_status == FoldEvidenceStatus.ASSESSED_WITH_HIT
    }
    authority_versions: dict[str, set[str]] = defaultdict(set)
    for structure in structures:
        if structure.fold_evidence_status in assessed_fold_statuses:
            authority_versions[structure.fold_authority].add(structure.fold_authority_version)
    conflicting_authorities = {
        authority: versions
        for authority, versions in authority_versions.items()
        if len(versions) > 1
    }
    if conflicting_authorities:
        raise InputValidationError(
            "Fold assessment universes cannot mix authority releases: "
            + "; ".join(
                f"{authority}={sorted(versions)!r}"
                for authority, versions in sorted(conflicting_authorities.items())
            )
        )
    for authority, authority_version, fold_id in sorted(positive_folds):
        assessed = frozenset(
            structure.protein_id
            for structure in structures
            if structure.fold_authority == authority
            and structure.fold_authority_version == authority_version
            and structure.analysis_eligibility_status in fold_eligible_states
            and structure.fold_evidence_status in assessed_fold_statuses
        )
        if assessed:
            result[("FOLD", f"{authority}:{fold_id}")] = assessed

    cluster_universes: dict[str, str] = {}
    for row in structure_clusters:
        cluster_id = str(row.get("cluster_id", "")).strip()
        universe_id = str(row.get("comparison_universe_id", "")).strip()
        if not cluster_id or not universe_id:
            raise InputValidationError(
                "Structure-cluster rows require cluster_id and comparison_universe_id."
            )
        previous = cluster_universes.setdefault(cluster_id, universe_id)
        if previous != universe_id:
            raise InputValidationError(
                f"Structure cluster {cluster_id!r} spans multiple comparison universes."
            )
        assessed = known_universes.get(universe_id)
        protein_id = str(row.get("protein_id", "")).strip()
        if assessed is not None and protein_id not in assessed:
            raise InputValidationError(
                f"Structure cluster {cluster_id!r} contains protein {protein_id!r} "
                f"outside declared comparison universe {universe_id!r}."
            )
    for cluster_id, universe_id in sorted(cluster_universes.items()):
        assessed = known_universes.get(universe_id)
        if assessed is not None:
            result[("STRUCTURE_CLUSTER", cluster_id)] = assessed

    LOGGER.info(
        "Derived %d proved structural feature-assessment universes; %d cluster "
        "definitions remain unknown without complete universe membership",
        len(result),
        sum(universe not in known_universes for universe in cluster_universes.values()),
    )
    return dict(sorted(result.items()))


def _connected_components(*, adjacency: dict[str, set[str]]) -> tuple[tuple[str, ...], ...]:
    """Find deterministic connected components in an undirected graph.

    Args:
        adjacency: Symmetric protein-neighbour mapping.

    Returns:
        Sorted components containing sorted protein identifiers.
    """

    remaining = set(adjacency)
    components: list[tuple[str, ...]] = []
    while remaining:
        start = min(remaining)
        pending = [start]
        observed: set[str] = set()
        while pending:
            current = pending.pop()
            if current in observed:
                continue
            observed.add(current)
            pending.extend(sorted(adjacency.get(current, set()) - observed, reverse=True))
        remaining -= observed
        components.append(tuple(sorted(observed)))
    return tuple(sorted(components))


def _feature(
    *,
    protein_id: str,
    feature_type: str,
    feature_id: str,
    feature_name: str,
    reference: str,
    derivation_scope: str,
    feature_definition_sha256: str,
    derivation_cohort_sha256: str = "",
    evidence_source: str = "protein-signature-analysis",
) -> FeatureRecord:
    """Build one consistently attributed structural feature.

    Args:
        protein_id: Campaign protein identifier.
        feature_type: Controlled feature category.
        feature_id: Stable feature identifier.
        feature_name: Human-readable name.
        reference: Evidence record or method reference.
        derivation_scope: Controlled feature-definition scope.
        feature_definition_sha256: Deterministic feature-definition digest.
        derivation_cohort_sha256: Exact data-derivation cohort digest, if applicable.
        evidence_source: Authority or package responsible for the evidence.

    Returns:
        Derived feature record.
    """

    return FeatureRecord(
        protein_id=protein_id,
        feature_type=feature_type,
        feature_id=feature_id,
        feature_name=feature_name,
        start=None,
        end=None,
        evidence_status="DERIVED",
        evidence_source=evidence_source,
        evidence_reference=reference,
        derivation_scope=derivation_scope,
        feature_definition_sha256=feature_definition_sha256,
        derivation_cohort_sha256=derivation_cohort_sha256,
    )
