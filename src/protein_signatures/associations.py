"""Comparison-aware feature enrichment and held-out signature synthesis."""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Collection, Mapping
from dataclasses import replace

from .assessment import (
    TECHNICAL_FEATURE_TYPES,
    FeatureAssessmentUniverses,
    FeatureKey,
    NormalisedAssessmentUniverses,
    normalise_feature_assessment_universes,
    normalise_feature_keys,
    validate_feature_inference_scopes,
)
from .errors import InputValidationError
from .models import (
    AnalysisSettings,
    AnalysisStatus,
    AssociationResult,
    ComparisonDefinition,
    FeatureRecord,
    PartitionAssignment,
    SignatureSummary,
)
from .statistics import benjamini_hochberg, fisher_exact_two_sided, wilson_score_interval

LOGGER = logging.getLogger(__name__)


def analyse_feature_associations(
    *,
    comparisons: tuple[ComparisonDefinition, ...],
    label_memberships: tuple[dict[str, str], ...],
    partitions: tuple[PartitionAssignment, ...],
    features: tuple[FeatureRecord, ...],
    settings: AnalysisSettings,
    feature_assessment_universes: FeatureAssessmentUniverses | None = None,
    exploratory_feature_keys: Collection[FeatureKey] = (),
) -> tuple[AssociationResult, ...]:
    """Test categorical features in each explicit comparison and partition.

    Args:
        comparisons: Target/background definitions.
        label_memberships: Expanded reviewed-positive protein memberships.
        partitions: Leakage-aware discovery/validation assignments.
        features: Protein-level feature records.
        settings: Sample, prevalence and FDR settings.
        feature_assessment_universes: Proteins successfully assessed for each
            feature. Proteins omitted from a supplied feature universe are
            unknown, not feature-negative. ``None`` is only appropriate when
            every feature was assessed for every partitioned protein.
        exploratory_feature_keys: Features derived from all data rather than a
            discovery-only reference. These remain eligible for exploratory
            discovery tests but are excluded from held-out validation tests.

    Returns:
        Exact association results with within-feature-type FDR correction.

    Raises:
        InputValidationError: If label sets overlap or feature names conflict.
    """

    labels_by_protein: dict[str, set[str]] = defaultdict(set)
    for row in label_memberships:
        labels_by_protein[row["protein_id"]].add(row["label_id"])
    partition_by_protein = {item.protein_id: item.partition for item in partitions}
    unit_by_protein = {item.protein_id: item.partition_key for item in partitions}
    feature_proteins: dict[tuple[str, str], set[str]] = defaultdict(set)
    feature_names: dict[tuple[str, str], str] = {}
    for feature in features:
        key = (feature.feature_type, feature.feature_id)
        previous = feature_names.get(key)
        if previous is not None and previous != feature.feature_name:
            raise InputValidationError(f"Conflicting names for feature {key!r}.")
        feature_names[key] = feature.feature_name
        feature_proteins[key].add(feature.protein_id)
    assessment_universes = normalise_feature_assessment_universes(
        universes=feature_assessment_universes,
        known_protein_ids=partition_by_protein,
        feature_proteins=feature_proteins,
    )
    exploratory_keys = validate_feature_inference_scopes(
        features=features,
        exploratory_feature_keys=exploratory_feature_keys,
        context="association",
    )
    if assessment_universes is None:
        LOGGER.warning(
            "No feature-assessment universes were supplied; association absence "
            "assumes every partitioned protein was assessed for every feature."
        )
    LOGGER.info(
        "Starting feature association analysis comparisons=%d feature_definitions=%d "
        "feature_rows=%d",
        len(comparisons),
        len(feature_proteins),
        len(features),
    )
    results: list[AssociationResult] = []
    for comparison_index, comparison in enumerate(comparisons, start=1):
        LOGGER.info(
            "Analysing feature comparison %d/%d comparison_id=%s",
            comparison_index,
            len(comparisons),
            comparison.comparison_id,
        )
        validation_candidate_keys: frozenset[FeatureKey] | None = None
        for partition in ("DISCOVERY", "VALIDATION"):
            eligible = {
                protein_id
                for protein_id, assigned_partition in partition_by_protein.items()
                if assigned_partition == partition
            }
            target = {
                protein_id
                for protein_id in eligible
                if labels_by_protein[protein_id] & set(comparison.target_label_ids)
            }
            background = {
                protein_id
                for protein_id in eligible
                if labels_by_protein[protein_id] & set(comparison.background_label_ids)
            }
            overlap = target & background
            if overlap:
                raise InputValidationError(
                    f"Comparison {comparison.comparison_id!r} has proteins in both classes: "
                    f"{sorted(overlap)[:10]}"
                )
            subset = _analyse_partition(
                comparison=comparison,
                partition=partition,
                target=target,
                background=background,
                feature_proteins=feature_proteins,
                feature_names=feature_names,
                settings=settings,
                unit_by_protein=unit_by_protein,
                feature_assessment_universes=assessment_universes,
                exploratory_feature_keys=exploratory_keys,
                candidate_feature_keys=(
                    validation_candidate_keys if partition == "VALIDATION" else None
                ),
            )
            adjusted_subset = _apply_fdr(results=subset)
            results.extend(adjusted_subset)
            if partition == "DISCOVERY":
                validation_candidate_keys = frozenset(
                    (item.feature_type, item.feature_id)
                    for item in adjusted_subset
                    if item.status == AnalysisStatus.COMPLETE
                    and item.q_value is not None
                    and item.q_value <= settings.fdr_threshold
                    and item.prevalence_difference not in {None, 0.0}
                )
    results = list(_apply_study_fdr(results=tuple(results)))
    LOGGER.info(
        "Completed feature association analysis comparisons=%d result_rows=%d",
        len(comparisons),
        len(results),
    )
    return tuple(
        sorted(
            results,
            key=lambda item: (
                item.comparison_id,
                item.partition,
                item.feature_type,
                item.q_value if item.q_value is not None else 2.0,
                item.feature_id,
            ),
        )
    )


def summarise_signatures(
    *,
    associations: tuple[AssociationResult, ...],
    fdr_threshold: float,
    exploratory_feature_keys: Collection[FeatureKey] = (),
) -> tuple[SignatureSummary, ...]:
    """Combine discovery associations with held-out validation evidence.

    Args:
        associations: Partition-specific association results.
        fdr_threshold: Inclusive significance threshold.
        exploratory_feature_keys: All-data-derived feature keys that cannot be
            promoted to decision candidates.

    Returns:
        Candidate signatures or explicit no-signature rows per comparison.
    """

    if not 0.0 <= fdr_threshold <= 1.0:
        raise InputValidationError("fdr_threshold must be between 0.0 and 1.0.")
    exploratory_keys = normalise_feature_keys(
        keys=exploratory_feature_keys,
        context="exploratory feature keys",
    )
    by_key = {
        (item.comparison_id, item.partition, item.feature_type, item.feature_id): item
        for item in associations
    }
    comparison_ids = sorted({item.comparison_id for item in associations})
    summaries: list[SignatureSummary] = []
    for comparison_id in comparison_ids:
        discoveries = [
            item
            for item in associations
            if item.comparison_id == comparison_id
            and item.partition == "DISCOVERY"
            and item.status == AnalysisStatus.COMPLETE
            and item.q_value is not None
            and item.q_value <= fdr_threshold
            and item.prevalence_difference not in {None, 0.0}
        ]
        if not discoveries:
            discovery_rows = [
                item
                for item in associations
                if item.comparison_id == comparison_id and item.partition == "DISCOVERY"
            ]
            summary_status = _no_discovery_status(discovery_rows=discovery_rows)
            completed_without_hit = summary_status == AnalysisStatus.NO_SIGNIFICANT_SIGNATURE
            summaries.append(
                SignatureSummary(
                    comparison_id=comparison_id,
                    feature_type="ANALYSIS",
                    feature_id=summary_status.value,
                    feature_name=(
                        "No discovery signature passed the configured FDR threshold"
                        if completed_without_hit
                        else "Discovery signature testing was not completed: "
                        f"{summary_status.value}"
                    ),
                    discovery_q_value=None,
                    discovery_study_q_value=None,
                    discovery_prevalence_difference=None,
                    validation_q_value=None,
                    validation_study_q_value=None,
                    validation_prevalence_difference=None,
                    evidence_class=(
                        "NO_DISCOVERY_SIGNATURE"
                        if completed_without_hit
                        else f"DISCOVERY_NOT_TESTED__{summary_status.value}"
                    ),
                    status=summary_status,
                )
            )
            continue
        for discovery in discoveries:
            validation = by_key.get(
                (
                    comparison_id,
                    "VALIDATION",
                    discovery.feature_type,
                    discovery.feature_id,
                )
            )
            evidence_class = _classify_discovery_scope(
                discovery=discovery,
                validation=validation,
                fdr_threshold=fdr_threshold,
                exploratory=(discovery.feature_type, discovery.feature_id) in exploratory_keys,
            )
            summaries.append(
                SignatureSummary(
                    comparison_id=comparison_id,
                    feature_type=discovery.feature_type,
                    feature_id=discovery.feature_id,
                    feature_name=discovery.feature_name,
                    discovery_q_value=discovery.q_value,
                    discovery_study_q_value=discovery.study_q_value,
                    discovery_prevalence_difference=discovery.prevalence_difference,
                    validation_q_value=validation.q_value if validation else None,
                    validation_study_q_value=(validation.study_q_value if validation else None),
                    validation_prevalence_difference=(
                        validation.prevalence_difference if validation else None
                    ),
                    evidence_class=evidence_class,
                    status=AnalysisStatus.COMPLETE,
                )
            )
    return tuple(
        sorted(
            summaries,
            key=lambda item: (
                item.comparison_id,
                item.discovery_q_value if item.discovery_q_value is not None else 2.0,
                item.feature_type,
                item.feature_id,
            ),
        )
    )


def _no_discovery_status(*, discovery_rows: Collection[AssociationResult]) -> AnalysisStatus:
    """Distinguish a negative completed screen from a screen not performed.

    Args:
        discovery_rows: Association rows from one discovery comparison.

    Returns:
        ``NO_SIGNIFICANT_SIGNATURE`` only when at least one eligible test was
        completed; otherwise the highest-priority explicit non-test status.
    """

    observed = {row.status for row in discovery_rows}
    if AnalysisStatus.COMPLETE in observed or AnalysisStatus.NO_SIGNIFICANT_SIGNATURE in observed:
        return AnalysisStatus.NO_SIGNIFICANT_SIGNATURE
    priority = (
        AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE,
        AnalysisStatus.NO_ELIGIBLE_FEATURES,
        AnalysisStatus.INPUT_UNAVAILABLE,
        AnalysisStatus.FAILED,
        AnalysisStatus.EXCLUDED,
        AnalysisStatus.AMBIGUOUS_CLASS,
        AnalysisStatus.MAPPING_FAILED,
        AnalysisStatus.NOT_SELECTED,
    )
    for status in priority:
        if status in observed:
            return status
    return AnalysisStatus.NO_ELIGIBLE_FEATURES


def _analyse_partition(
    *,
    comparison: ComparisonDefinition,
    partition: str,
    target: set[str],
    background: set[str],
    feature_proteins: dict[tuple[str, str], set[str]],
    feature_names: dict[tuple[str, str], str],
    settings: AnalysisSettings,
    unit_by_protein: dict[str, str] | None = None,
    feature_assessment_universes: NormalisedAssessmentUniverses | None = None,
    exploratory_feature_keys: Collection[FeatureKey] = (),
    candidate_feature_keys: Collection[FeatureKey] | None = None,
) -> tuple[AssociationResult, ...]:
    """Calculate unadjusted associations for one comparison partition.

    Args:
        comparison: Comparison definition.
        partition: Partition name.
        target: Target protein identifiers.
        background: Background protein identifiers.
        feature_proteins: Proteins carrying each feature.
        feature_names: Human-readable feature names.
        settings: Analysis thresholds.
        unit_by_protein: Proteins mapped to indivisible homology/redundancy blocks.
        feature_assessment_universes: Validated proteins successfully assessed
            for each feature. ``None`` means all proteins were assessed.
        exploratory_feature_keys: All-data-derived features excluded from
            validation inference.
        candidate_feature_keys: Optional frozen feature vocabulary for this
            partition. The public orchestrator uses locally significant,
            non-zero-effect discovery features for held-out validation.

    Returns:
        Unadjusted results or one explicit status row.
    """

    eligible = target | background
    effective_units = unit_by_protein or {protein_id: protein_id for protein_id in eligible}
    missing_units = eligible - effective_units.keys()
    if missing_units:
        raise InputValidationError(
            "Association samples are missing independence-block assignments: "
            f"{sorted(missing_units)[:10]}"
        )
    target_units = {effective_units[protein_id] for protein_id in target}
    background_units = {effective_units[protein_id] for protein_id in background}
    mixed_units = target_units & background_units
    target_units -= mixed_units
    background_units -= mixed_units
    target_members = _proteins_by_unit(
        protein_ids=target,
        unit_by_protein=effective_units,
        excluded_units=mixed_units,
    )
    background_members = _proteins_by_unit(
        protein_ids=background,
        unit_by_protein=effective_units,
        excluded_units=mixed_units,
    )
    if (
        len(target_units) < settings.minimum_target_proteins
        or len(background_units) < settings.minimum_background_proteins
    ):
        return (
            _status_result(
                comparison_id=comparison.comparison_id,
                partition=partition,
                target_count=len(target),
                background_count=len(background),
                target_unit_count=len(target_units),
                background_unit_count=len(background_units),
                excluded_mixed_unit_count=len(mixed_units),
                status=AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE,
            ),
        )
    exploratory_keys = frozenset(exploratory_feature_keys)
    if candidate_feature_keys is not None:
        candidates = set(candidate_feature_keys)
        unknown_candidates = candidates - feature_proteins.keys()
        if unknown_candidates:
            raise InputValidationError(
                "Candidate feature keys do not occur in the feature records: "
                f"{sorted(unknown_candidates)[:10]}"
            )
    elif partition == "VALIDATION":
        candidates = set(feature_proteins)
    else:
        candidates = {
            key
            for key, proteins in feature_proteins.items()
            if len(proteins & eligible) >= settings.minimum_feature_proteins
        }
    if not candidates:
        return (
            _status_result(
                comparison_id=comparison.comparison_id,
                partition=partition,
                target_count=len(target),
                background_count=len(background),
                target_unit_count=len(target_units),
                background_unit_count=len(background_units),
                excluded_mixed_unit_count=len(mixed_units),
                status=(
                    AnalysisStatus.NOT_SELECTED
                    if candidate_feature_keys is not None
                    else AnalysisStatus.NO_ELIGIBLE_FEATURES
                ),
            ),
        )
    rows: list[AssociationResult] = []
    for feature_type, feature_id in sorted(candidates):
        feature_key = (feature_type, feature_id)
        proteins = feature_proteins[feature_key]
        target_with = len(target & proteins)
        background_with = len(background & proteins)
        if feature_assessment_universes is None:
            assessed_proteins = frozenset(eligible)
            analysis_unit = "INDEPENDENCE_BLOCK"
            interval_method = "WILSON_95_PERCENT_INDEPENDENCE_BLOCKS"
        else:
            try:
                assessed_proteins = feature_assessment_universes[feature_key]
            except KeyError as error:
                raise InputValidationError(
                    f"Missing assessment universe for feature {feature_key!r}."
                ) from error
            unassessed_positives = (proteins & eligible) - assessed_proteins
            if unassessed_positives:
                raise InputValidationError(
                    f"Observed positives for feature {feature_key!r} are not assessed: "
                    f"{sorted(unassessed_positives)[:10]}"
                )
            analysis_unit = "FEATURE_ASSESSED_INDEPENDENCE_BLOCK"
            interval_method = "WILSON_95_PERCENT_FEATURE_ASSESSED_BLOCKS"
        target_denominator, target_units_with = _feature_block_counts(
            unit_members=target_members,
            positive_proteins=proteins,
            assessed_proteins=assessed_proteins,
        )
        background_denominator, background_units_with = _feature_block_counts(
            unit_members=background_members,
            positive_proteins=proteins,
            assessed_proteins=assessed_proteins,
        )
        if partition == "VALIDATION" and feature_key in exploratory_keys:
            rows.append(
                _feature_status_result(
                    comparison_id=comparison.comparison_id,
                    partition=partition,
                    analysis_unit=analysis_unit,
                    feature_type=feature_type,
                    feature_id=feature_id,
                    feature_name=feature_names[feature_key],
                    target_count=len(target),
                    background_count=len(background),
                    target_with_feature=target_with,
                    background_with_feature=background_with,
                    target_unit_count=len(target_units),
                    background_unit_count=len(background_units),
                    target_assessed_unit_count=len(target_denominator),
                    background_assessed_unit_count=len(background_denominator),
                    target_units_with_feature=len(target_units_with),
                    background_units_with_feature=len(background_units_with),
                    excluded_mixed_unit_count=len(mixed_units),
                    status=AnalysisStatus.EXCLUDED,
                    prevalence_ci_method="NOT_ESTIMATED_ALL_DATA_DERIVATION",
                )
            )
            continue
        if (
            len(target_denominator) < settings.minimum_target_proteins
            or len(background_denominator) < settings.minimum_background_proteins
        ):
            rows.append(
                _feature_status_result(
                    comparison_id=comparison.comparison_id,
                    partition=partition,
                    analysis_unit=analysis_unit,
                    feature_type=feature_type,
                    feature_id=feature_id,
                    feature_name=feature_names[feature_key],
                    target_count=len(target),
                    background_count=len(background),
                    target_with_feature=target_with,
                    background_with_feature=background_with,
                    target_unit_count=len(target_units),
                    background_unit_count=len(background_units),
                    target_assessed_unit_count=len(target_denominator),
                    background_assessed_unit_count=len(background_denominator),
                    target_units_with_feature=len(target_units_with),
                    background_units_with_feature=len(background_units_with),
                    excluded_mixed_unit_count=len(mixed_units),
                    status=AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE,
                    prevalence_ci_method=("NOT_ESTIMATED_FEATURE_ASSESSMENT_UNDERPOWERED"),
                )
            )
            continue
        target_prevalence = len(target_units_with) / len(target_denominator)
        background_prevalence = len(background_units_with) / len(background_denominator)
        target_interval = wilson_score_interval(
            successes=len(target_units_with), trials=len(target_denominator)
        )
        background_interval = wilson_score_interval(
            successes=len(background_units_with), trials=len(background_denominator)
        )
        difference = target_prevalence - background_prevalence
        odds_ratio, p_value = fisher_exact_two_sided(
            a=len(target_units_with),
            b=len(target_denominator) - len(target_units_with),
            c=len(background_units_with),
            d=len(background_denominator) - len(background_units_with),
        )
        rows.append(
            AssociationResult(
                comparison_id=comparison.comparison_id,
                partition=partition,
                analysis_unit=analysis_unit,
                feature_type=feature_type,
                feature_id=feature_id,
                feature_name=feature_names[feature_key],
                target_protein_count=len(target),
                background_protein_count=len(background),
                target_with_feature=target_with,
                background_with_feature=background_with,
                target_unit_count=len(target_units),
                background_unit_count=len(background_units),
                target_assessed_unit_count=len(target_denominator),
                background_assessed_unit_count=len(background_denominator),
                target_unknown_unit_count=len(target_units) - len(target_denominator),
                background_unknown_unit_count=(len(background_units) - len(background_denominator)),
                target_units_with_feature=len(target_units_with),
                background_units_with_feature=len(background_units_with),
                excluded_mixed_unit_count=len(mixed_units),
                target_prevalence=target_prevalence,
                background_prevalence=background_prevalence,
                target_prevalence_ci_lower=target_interval[0],
                target_prevalence_ci_upper=target_interval[1],
                background_prevalence_ci_lower=background_interval[0],
                background_prevalence_ci_upper=background_interval[1],
                prevalence_ci_method=interval_method,
                prevalence_difference=difference,
                odds_ratio=odds_ratio,
                p_value=p_value,
                q_value=None,
                study_q_value=None,
                direction=(
                    "TARGET_ENRICHED"
                    if difference > 0
                    else "BACKGROUND_ENRICHED"
                    if difference < 0
                    else "NO_DIFFERENCE"
                ),
                status=AnalysisStatus.COMPLETE,
            )
        )
    return tuple(rows)


def _proteins_by_unit(
    *,
    protein_ids: Collection[str],
    unit_by_protein: Mapping[str, str],
    excluded_units: Collection[str],
) -> dict[str, frozenset[str]]:
    """Group proteins into retained inference units.

    Args:
        protein_ids: Proteins in one comparison class.
        unit_by_protein: Protein-to-independence-block mapping.
        excluded_units: Mixed-class blocks that cannot enter inference.

    Returns:
        Retained unit identifiers mapped to their class-member proteins.
    """

    excluded = frozenset(excluded_units)
    grouped: dict[str, set[str]] = defaultdict(set)
    for protein_id in protein_ids:
        unit = unit_by_protein[protein_id]
        if unit not in excluded:
            grouped[unit].add(protein_id)
    return {unit: frozenset(members) for unit, members in grouped.items()}


def _feature_block_counts(
    *,
    unit_members: Mapping[str, Collection[str]],
    positive_proteins: Collection[str],
    assessed_proteins: Collection[str],
) -> tuple[frozenset[str], frozenset[str]]:
    """Return assessed denominator blocks and feature-present blocks.

    Presence is known if at least one block member has the feature. Absence is
    known only if every relevant member of that block was assessed and no member
    was positive. A partly assessed, no-hit block is therefore unknown.

    Args:
        unit_members: Retained inference blocks and their class members.
        positive_proteins: Proteins observed with the feature.
        assessed_proteins: Proteins successfully assessed for the feature.

    Returns:
        Pair of denominator block identifiers and feature-present block
        identifiers.
    """

    positive = frozenset(positive_proteins)
    assessed = frozenset(assessed_proteins)
    present_units: set[str] = set()
    assessed_absent_units: set[str] = set()
    for unit, raw_members in unit_members.items():
        members = frozenset(raw_members)
        if members & positive:
            present_units.add(unit)
        elif members <= assessed:
            assessed_absent_units.add(unit)
    denominator = present_units | assessed_absent_units
    return frozenset(denominator), frozenset(present_units)


def _feature_status_result(
    *,
    comparison_id: str,
    partition: str,
    analysis_unit: str,
    feature_type: str,
    feature_id: str,
    feature_name: str,
    target_count: int,
    background_count: int,
    target_with_feature: int,
    background_with_feature: int,
    target_unit_count: int,
    background_unit_count: int,
    target_assessed_unit_count: int,
    background_assessed_unit_count: int,
    target_units_with_feature: int,
    background_units_with_feature: int,
    excluded_mixed_unit_count: int,
    status: AnalysisStatus,
    prevalence_ci_method: str,
) -> AssociationResult:
    """Construct a feature-specific underpowered assessment row.

    Args:
        comparison_id: Comparison identifier.
        partition: Discovery or validation partition.
        analysis_unit: Denominator semantics for the feature.
        feature_type: Feature evidence family.
        feature_id: Feature identifier.
        feature_name: Human-readable feature name.
        target_count: Raw target protein count.
        background_count: Raw background protein count.
        target_with_feature: Raw target positive count.
        background_with_feature: Raw background positive count.
        target_unit_count: Total retained target block count.
        background_unit_count: Total retained background block count.
        target_assessed_unit_count: Assessed target block count.
        background_assessed_unit_count: Assessed background block count.
        target_units_with_feature: Target block positive count.
        background_units_with_feature: Background block positive count.
        excluded_mixed_unit_count: Mixed-class blocks excluded from inference.
        status: Controlled reason inference was not performed.
        prevalence_ci_method: Explicit non-estimation reason text.

    Returns:
        Association-shaped row with feature identity and assessed counts.
    """

    base = _status_result(
        comparison_id=comparison_id,
        partition=partition,
        target_count=target_count,
        background_count=background_count,
        target_unit_count=target_unit_count,
        background_unit_count=background_unit_count,
        excluded_mixed_unit_count=excluded_mixed_unit_count,
        status=status,
    )
    return replace(
        base,
        analysis_unit=analysis_unit,
        feature_type=feature_type,
        feature_id=feature_id,
        feature_name=feature_name,
        target_with_feature=target_with_feature,
        background_with_feature=background_with_feature,
        target_assessed_unit_count=target_assessed_unit_count,
        background_assessed_unit_count=background_assessed_unit_count,
        target_unknown_unit_count=target_unit_count - target_assessed_unit_count,
        background_unknown_unit_count=(background_unit_count - background_assessed_unit_count),
        target_units_with_feature=target_units_with_feature,
        background_units_with_feature=background_units_with_feature,
        prevalence_ci_method=prevalence_ci_method,
    )


def _apply_fdr(*, results: tuple[AssociationResult, ...]) -> tuple[AssociationResult, ...]:
    """Apply Benjamini-Hochberg correction within feature types.

    Args:
        results: One comparison/partition result set.

    Returns:
        Results with q-values populated where tests were complete.
    """

    positions: dict[str, list[int]] = defaultdict(list)
    for index, result in enumerate(results):
        if result.status == AnalysisStatus.COMPLETE and result.p_value is not None:
            positions[result.feature_type].append(index)
    adjusted = list(results)
    for indices in positions.values():
        q_values = benjamini_hochberg(
            p_values=tuple(float(results[index].p_value) for index in indices)
        )
        for index, q_value in zip(indices, q_values, strict=True):
            adjusted[index] = replace(results[index], q_value=q_value)
    return tuple(adjusted)


def _apply_study_fdr(*, results: tuple[AssociationResult, ...]) -> tuple[AssociationResult, ...]:
    """Correct across comparisons within each partition and evidence family.

    Args:
        results: Locally corrected association rows from the complete campaign.

    Returns:
        Rows with a second, study-wide q-value populated for complete tests.

    Notes:
        The feature type remains the statistical evidence family. This prevents a very
        large amino-acid k-mer search from consuming the error budget for a smaller,
        prespecified Pfam or structural-fold family while still correcting for inspection
        across all configured protein-type comparisons.
    """

    positions: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, result in enumerate(results):
        if result.status == AnalysisStatus.COMPLETE and result.p_value is not None:
            positions[(result.partition, result.feature_type)].append(index)
    adjusted = list(results)
    for indices in positions.values():
        q_values = benjamini_hochberg(
            p_values=tuple(float(results[index].p_value) for index in indices)
        )
        for index, q_value in zip(indices, q_values, strict=True):
            adjusted[index] = replace(results[index], study_q_value=q_value)
    return tuple(adjusted)


def _status_result(
    *,
    comparison_id: str,
    partition: str,
    target_count: int,
    background_count: int,
    target_unit_count: int,
    background_unit_count: int,
    excluded_mixed_unit_count: int,
    status: AnalysisStatus,
) -> AssociationResult:
    """Construct one explicit non-test status row.

    Args:
        comparison_id: Comparison identifier.
        partition: Partition name.
        target_count: Observed target proteins.
        background_count: Observed background proteins.
        target_unit_count: Independent target blocks retained for inference.
        background_unit_count: Independent background blocks retained for inference.
        excluded_mixed_unit_count: Blocks containing both target and background proteins.
        status: Controlled reason that tests were not produced.

    Returns:
        Association-shaped status record.
    """

    return AssociationResult(
        comparison_id=comparison_id,
        partition=partition,
        analysis_unit="INDEPENDENCE_BLOCK",
        feature_type="ANALYSIS",
        feature_id=status.value,
        feature_name=status.value.replace("_", " ").title(),
        target_protein_count=target_count,
        background_protein_count=background_count,
        target_with_feature=0,
        background_with_feature=0,
        target_unit_count=target_unit_count,
        background_unit_count=background_unit_count,
        target_assessed_unit_count=target_unit_count,
        background_assessed_unit_count=background_unit_count,
        target_unknown_unit_count=0,
        background_unknown_unit_count=0,
        target_units_with_feature=0,
        background_units_with_feature=0,
        excluded_mixed_unit_count=excluded_mixed_unit_count,
        target_prevalence=None,
        background_prevalence=None,
        target_prevalence_ci_lower=None,
        target_prevalence_ci_upper=None,
        background_prevalence_ci_lower=None,
        background_prevalence_ci_upper=None,
        prevalence_ci_method="NOT_ESTIMATED",
        prevalence_difference=None,
        odds_ratio=None,
        p_value=None,
        q_value=None,
        study_q_value=None,
        direction="NOT_TESTED",
        status=status,
    )


def _classify_discovery_scope(
    *,
    discovery: AssociationResult,
    validation: AssociationResult | None,
    fdr_threshold: float,
    exploratory: bool = False,
) -> str:
    """Classify discovery scope before appending held-out support.

    Args:
        discovery: Locally significant discovery association.
        validation: Matching held-out association, when available.
        fdr_threshold: Inclusive FDR decision threshold.
        exploratory: Whether the feature definition used all campaign data and
            therefore cannot support held-out confirmation.

    Returns:
        Evidence class separating study-wide decision candidates, local-only
        exploratory findings and technical availability quality-control rows.
    """

    validation_class = _classify_validation(
        discovery=discovery,
        validation=validation,
        fdr_threshold=fdr_threshold,
    )
    if discovery.feature_type in TECHNICAL_FEATURE_TYPES:
        return f"QC_TECHNICAL_NON_BIOLOGICAL__{validation_class}"
    if exploratory:
        return f"EXPLORATORY_ALL_DATA_DERIVATION__{validation_class}"
    if discovery.study_q_value is not None and discovery.study_q_value <= fdr_threshold:
        return f"DECISION_CANDIDATE__{validation_class}"
    return f"EXPLORATORY_LOCAL_ONLY__{validation_class}"


def _classify_validation(
    *,
    discovery: AssociationResult,
    validation: AssociationResult | None,
    fdr_threshold: float,
) -> str:
    """Classify held-out support for one discovery feature.

    Args:
        discovery: Significant discovery association.
        validation: Matching validation result, when available.
        fdr_threshold: Inclusive validation FDR threshold.

    Returns:
        Controlled evidence-class text.
    """

    if validation is None or validation.status != AnalysisStatus.COMPLETE:
        return "NO_VALIDATION_DATA"
    discovery_effect = discovery.prevalence_difference or 0.0
    validation_effect = validation.prevalence_difference or 0.0
    if discovery_effect * validation_effect < 0:
        return "DIRECTION_DISCORDANT"
    if validation.q_value is not None and validation.q_value <= fdr_threshold:
        if validation.study_q_value is not None and validation.study_q_value <= fdr_threshold:
            return "VALIDATED_STUDY_WIDE"
        return "VALIDATED_WITHIN_COMPARISON"
    return "DISCOVERY_ONLY"
