"""Unit tests for sequence, structural and statistical signature derivation."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest

from protein_signatures.associations import (
    _analyse_partition,
    _apply_fdr,
    _apply_study_fdr,
    _classify_discovery_scope,
    _classify_validation,
    _feature_block_counts,
    _no_discovery_status,
    _proteins_by_unit,
    _status_result,
    analyse_feature_associations,
    summarise_signatures,
)
from protein_signatures.errors import InputValidationError
from protein_signatures.models import (
    AlphaFoldAcquisition,
    AnalysisSettings,
    AnalysisStatus,
    AssociationResult,
    ComparisonDefinition,
    FeatureRecord,
    FoldEvidenceStatus,
    PairwiseStructureComparison,
    PartitionAssignment,
    SequenceRecord,
    StructureAnalysisEligibility,
    StructureComparisonStatus,
    StructureCoverageScope,
    StructureRecord,
)
from protein_signatures.sequence_signatures import build_kmer_features
from protein_signatures.statistics import (
    _hypergeometric_probability,
    _log_combination,
    benjamini_hochberg,
    fisher_exact_two_sided,
    wilson_score_interval,
)
from protein_signatures.structural_signatures import (
    _comparison_passes,
    _connected_components,
    derive_structure_assessment_universes,
    derive_structure_features,
)


def test_exact_statistics_cover_boundaries_and_known_result() -> None:
    """Exact tests should validate counts and return stable numerical results."""

    odds, p_value = fisher_exact_two_sided(a=8, b=2, c=1, d=9)
    assert odds == pytest.approx(36.0)
    assert p_value == pytest.approx(0.005477494641581322)
    assert fisher_exact_two_sided(a=0, b=0, c=0, d=0) == (None, 1.0)
    assert fisher_exact_two_sided(a=1, b=0, c=0, d=1)[0] == math.inf
    with pytest.raises(InputValidationError):
        fisher_exact_two_sided(a=True, b=0, c=0, d=1)
    assert benjamini_hochberg(p_values=()) == ()
    assert benjamini_hochberg(p_values=(0.01, 0.04, 0.03)) == pytest.approx((0.03, 0.04, 0.04))
    assert wilson_score_interval(successes=5, trials=10) == pytest.approx(
        (0.236593090512564, 0.7634069094874361)
    )
    for values in ((math.nan,), (-0.1,), (1.1,)):
        with pytest.raises(InputValidationError):
            benjamini_hochberg(p_values=values)
    for successes, trials, z_score in ((-1, 2, 1.96), (3, 2, 1.96), (1, 0, 1.96), (1, 2, 0.0)):
        with pytest.raises(InputValidationError):
            wilson_score_interval(
                successes=successes,
                trials=trials,
                z_score=z_score,
            )


def test_probability_helpers_reject_impossible_tables() -> None:
    """Probability helpers should return zero or negative infinity off support."""

    assert _hypergeometric_probability(x=-1, row_target=2, column_present=2, total=4) == 0.0
    assert _hypergeometric_probability(x=3, row_target=2, column_present=3, total=4) == 0.0
    assert _hypergeometric_probability(x=0, row_target=3, column_present=2, total=4) == 0.0
    assert _hypergeometric_probability(x=1, row_target=2, column_present=2, total=4) > 0.0
    assert _log_combination(n=2, k=-1) == -math.inf
    assert _log_combination(n=2, k=3) == -math.inf
    assert _log_combination(n=4, k=2) == pytest.approx(math.log(6))


def test_kmer_generation_and_safeguards() -> None:
    """K-mer extraction should filter prevalence and enforce cardinality limits."""

    sequences = (
        _sequence(protein_id="p1", sequence="AAAA"),
        _sequence(protein_id="p2", sequence="AAAT"),
    )
    rows = build_kmer_features(
        sequences=sequences,
        discovery_protein_ids=frozenset({"p1", "p2"}),
        lengths=(2,),
        minimum_proteins=2,
        maximum_features=10,
    )
    assert {(row.protein_id, row.feature_id) for row in rows} == {
        ("p1", "k2:AA"),
        ("p2", "k2:AA"),
    }
    for lengths, minimum, maximum in (((), 1, 1), ((0,), 1, 1), ((1,), 0, 1)):
        with pytest.raises(InputValidationError):
            build_kmer_features(
                sequences=sequences,
                discovery_protein_ids=frozenset({"p1", "p2"}),
                lengths=lengths,
                minimum_proteins=minimum,
                maximum_features=maximum,
            )
    with pytest.raises(InputValidationError, match="safeguard"):
        build_kmer_features(
            sequences=sequences,
            discovery_protein_ids=frozenset({"p1", "p2"}),
            lengths=(2,),
            minimum_proteins=1,
            maximum_features=1,
        )


def test_kmer_vocabulary_is_discovery_only_and_projects_to_validation() -> None:
    """Validation-only k-mers must not define features or affect safeguards."""

    sequences = (
        _sequence(protein_id="d1", sequence="AAAA"),
        _sequence(protein_id="d2", sequence="AAAT"),
        _sequence(protein_id="v1", sequence="AATT"),
        _sequence(protein_id="v2", sequence="CCGG"),
    )
    rows = build_kmer_features(
        sequences=sequences,
        discovery_protein_ids=frozenset({"d1", "d2"}),
        lengths=(2,),
        minimum_proteins=2,
        maximum_features=2,
    )

    assert {(row.protein_id, row.feature_id) for row in rows} == {
        ("d1", "k2:AA"),
        ("d2", "k2:AA"),
        ("v1", "k2:AA"),
    }
    assert all(row.protein_id != "v2" for row in rows)
    assert all(
        row.evidence_reference == "discovery_defined_exact_protein_level_kmer_presence"
        for row in rows
    )
    assert all(row.derivation_scope == "DISCOVERY_DERIVED" for row in rows)
    assert all(len(row.feature_definition_sha256) == 64 for row in rows)
    assert len({row.derivation_cohort_sha256 for row in rows}) == 1
    repeated = build_kmer_features(
        sequences=tuple(reversed(sequences)),
        discovery_protein_ids=frozenset({"d1", "d2"}),
        lengths=(2,),
        minimum_proteins=2,
        maximum_features=2,
    )
    assert repeated == rows


def test_kmer_generation_rejects_invalid_discovery_sets() -> None:
    """Discovery IDs must be non-empty, known and backed by unique records."""

    sequences = (_sequence(protein_id="p1", sequence="AAAA"),)
    for discovery_ids, message in (
        (frozenset(), "At least one discovery protein"),
        (frozenset({"missing"}), "absent from the supplied sequences"),
    ):
        with pytest.raises(InputValidationError, match=message):
            build_kmer_features(
                sequences=sequences,
                discovery_protein_ids=discovery_ids,
                lengths=(2,),
                minimum_proteins=1,
                maximum_features=10,
            )
    with pytest.raises(InputValidationError, match="must be unique"):
        build_kmer_features(
            sequences=(sequences[0], sequences[0]),
            discovery_protein_ids=frozenset({"p1"}),
            lengths=(2,),
            minimum_proteins=1,
            maximum_features=10,
        )


def test_partition_association_statuses_and_directions() -> None:
    """Partition analysis should expose sample, feature and direction states."""

    comparison = _comparison()
    strict = _settings(minimum_target=2, minimum_background=2, minimum_feature=2)
    insufficient = _analyse_partition(
        comparison=comparison,
        partition="DISCOVERY",
        target={"t1"},
        background={"b1", "b2"},
        feature_proteins={},
        feature_names={},
        settings=strict,
    )
    assert insufficient[0].status == AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE
    no_features = _analyse_partition(
        comparison=comparison,
        partition="DISCOVERY",
        target={"t1", "t2"},
        background={"b1", "b2"},
        feature_proteins={("TYPE", "rare"): {"t1"}},
        feature_names={("TYPE", "rare"): "Rare"},
        settings=strict,
    )
    assert no_features[0].status == AnalysisStatus.NO_ELIGIBLE_FEATURES
    feature_sets = {
        ("TYPE", "target"): {"t1", "t2"},
        ("TYPE", "background"): {"b1", "b2"},
        ("TYPE", "equal"): {"t1", "b1"},
    }
    complete = _analyse_partition(
        comparison=comparison,
        partition="DISCOVERY",
        target={"t1", "t2"},
        background={"b1", "b2"},
        feature_proteins=feature_sets,
        feature_names={key: key[1].title() for key in feature_sets},
        settings=_settings(minimum_target=1, minimum_background=1, minimum_feature=1),
    )
    assert {row.direction for row in complete} == {
        "TARGET_ENRICHED",
        "BACKGROUND_ENRICHED",
        "NO_DIFFERENCE",
    }
    adjusted = _apply_fdr(results=(*complete, insufficient[0]))
    assert all(row.q_value is not None for row in adjusted[:3])
    assert adjusted[-1].q_value is None
    study_adjusted = _apply_study_fdr(results=adjusted)
    assert all(row.study_q_value is not None for row in study_adjusted[:3])
    assert study_adjusted[-1].study_q_value is None
    status = _status_result(
        comparison_id="cmp",
        partition="VALIDATION",
        target_count=0,
        background_count=0,
        target_unit_count=0,
        background_unit_count=0,
        excluded_mixed_unit_count=0,
        status=AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE,
    )
    assert status.direction == "NOT_TESTED"


def test_association_inference_excludes_mixed_independence_blocks() -> None:
    """A homology block containing both classes must not count twice in Fisher tests."""

    comparison = _comparison()
    rows = _analyse_partition(
        comparison=comparison,
        partition="DISCOVERY",
        target={"t1", "t2"},
        background={"b1", "b2"},
        feature_proteins={("TYPE", "x"): {"t1", "t2"}},
        feature_names={("TYPE", "x"): "X"},
        settings=_settings(),
        unit_by_protein={"t1": "mixed", "b1": "mixed", "t2": "t2", "b2": "b2"},
    )
    assert rows[0].excluded_mixed_unit_count == 1
    assert rows[0].target_protein_count == 2
    assert rows[0].target_unit_count == 1
    assert rows[0].target_units_with_feature == 1


def test_feature_assessment_universes_exclude_unknown_blocks() -> None:
    """Only assessed absence, never unknown assay state, should enter denominators."""

    grouped = _proteins_by_unit(
        protein_ids={"t1", "t2", "t3"},
        unit_by_protein={"t1": "shared", "t2": "shared", "t3": "single"},
        excluded_units=(),
    )
    denominator, present = _feature_block_counts(
        unit_members=grouped,
        positive_proteins={"t1"},
        assessed_proteins={"t1", "t3"},
    )
    assert denominator == frozenset({"shared", "single"})
    assert present == frozenset({"shared"})

    key = ("TYPE", "x")
    arguments = {
        "comparison": _comparison(),
        "partition": "DISCOVERY",
        "target": {"t1", "t2", "t3"},
        "background": {"b1", "b2"},
        "feature_proteins": {key: {"t1"}},
        "feature_names": {key: "X"},
        "unit_by_protein": {
            "t1": "shared",
            "t2": "shared",
            "t3": "single",
            "b1": "b1",
            "b2": "b2",
        },
        "feature_assessment_universes": {key: frozenset({"t1", "t3", "b1"})},
    }
    underpowered = _analyse_partition(
        **arguments,
        settings=_settings(minimum_target=1, minimum_background=2),
    )[0]
    assert underpowered.status == AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE
    assert underpowered.feature_id == "x"
    assert underpowered.target_unit_count == 2
    assert underpowered.target_assessed_unit_count == 2
    assert underpowered.target_unknown_unit_count == 0
    assert underpowered.background_unit_count == 2
    assert underpowered.background_assessed_unit_count == 1
    assert underpowered.background_unknown_unit_count == 1
    complete = _analyse_partition(
        **arguments,
        settings=_settings(minimum_target=1, minimum_background=1),
    )[0]
    assert complete.status == AnalysisStatus.COMPLETE
    assert complete.analysis_unit == "FEATURE_ASSESSED_INDEPENDENCE_BLOCK"
    assert complete.target_prevalence == pytest.approx(0.5)
    assert complete.background_prevalence == pytest.approx(0.0)


def test_all_data_features_are_excluded_from_validation_inference() -> None:
    """All-data feature definitions should retain ledgers without pretending validation."""

    key = ("TYPE", "all_data")
    rows = _analyse_partition(
        comparison=_comparison(),
        partition="VALIDATION",
        target={"t1"},
        background={"b1"},
        feature_proteins={key: {"t1"}},
        feature_names={key: "All-data feature"},
        settings=_settings(),
        feature_assessment_universes={key: frozenset({"t1", "b1"})},
        exploratory_feature_keys={key},
    )
    assert rows[0].status == AnalysisStatus.EXCLUDED
    assert rows[0].feature_id == "all_data"
    assert rows[0].prevalence_ci_method == "NOT_ESTIMATED_ALL_DATA_DERIVATION"


def test_validation_tests_the_frozen_feature_vocabulary_without_prevalence_filtering() -> None:
    """An assessed discovery feature absent in validation is non-replication, not missing data."""

    key = ("TYPE", "discovery_only_positive")
    rows = _analyse_partition(
        comparison=_comparison(),
        partition="VALIDATION",
        target={"t1"},
        background={"b1"},
        feature_proteins={key: {"d1"}},
        feature_names={key: "Frozen discovery feature"},
        settings=_settings(minimum_feature=2),
        feature_assessment_universes={key: frozenset({"d1", "t1", "b1"})},
    )
    assert len(rows) == 1
    assert rows[0].status == AnalysisStatus.COMPLETE
    assert rows[0].target_units_with_feature == 0
    assert rows[0].background_units_with_feature == 0
    assert rows[0].prevalence_difference == pytest.approx(0.0)
    assert rows[0].p_value == pytest.approx(1.0)


def test_association_orchestrator_rejects_ambiguous_evidence() -> None:
    """The public association function should fail on overlap and name conflicts."""

    partitions = (
        PartitionAssignment("p1", "DISCOVERY", "EXACT_SEQUENCE", "a"),
        PartitionAssignment("p2", "DISCOVERY", "EXACT_SEQUENCE", "b"),
    )
    feature_a = _feature(protein_id="p1", feature_id="x", name="One")
    feature_b = _feature(protein_id="p2", feature_id="x", name="Two")
    with pytest.raises(InputValidationError, match="Conflicting names"):
        analyse_feature_associations(
            comparisons=(_comparison(),),
            label_memberships=(),
            partitions=partitions,
            features=(feature_a, feature_b),
            settings=_settings(),
        )
    with pytest.raises(InputValidationError, match="both classes"):
        analyse_feature_associations(
            comparisons=(_comparison(),),
            label_memberships=(
                {"protein_id": "p1", "label_id": "target"},
                {"protein_id": "p1", "label_id": "background"},
            ),
            partitions=partitions,
            features=(feature_a,),
            settings=_settings(),
        )
    results = analyse_feature_associations(
        comparisons=(_comparison(),),
        label_memberships=(
            {"protein_id": "p1", "label_id": "target"},
            {"protein_id": "p2", "label_id": "background"},
        ),
        partitions=partitions,
        features=(feature_a,),
        settings=_settings(),
    )
    assert {row.partition for row in results} == {"DISCOVERY", "VALIDATION"}


def test_orchestrator_validates_only_frozen_discovery_candidates() -> None:
    """Held-out tests should be selected only by discovery and retain zero-positive rows."""

    discovery_targets = tuple(f"dt{index}" for index in range(6))
    discovery_backgrounds = tuple(f"db{index}" for index in range(6))
    validation_targets = ("vt1", "vt2")
    validation_backgrounds = ("vb1", "vb2")
    protein_ids = (
        *discovery_targets,
        *discovery_backgrounds,
        *validation_targets,
        *validation_backgrounds,
    )
    partitions = tuple(
        PartitionAssignment(
            protein_id,
            "VALIDATION" if protein_id.startswith("v") else "DISCOVERY",
            "GROUP",
            protein_id,
        )
        for protein_id in protein_ids
    )
    memberships = tuple(
        {
            "protein_id": protein_id,
            "label_id": "target"
            if protein_id in {*discovery_targets, *validation_targets}
            else "background",
        }
        for protein_id in protein_ids
    )
    candidate = ("TYPE", "candidate")
    nonsignificant = ("TYPE", "nonsignificant")
    features = tuple(
        [
            FeatureRecord(
                protein_id,
                candidate[0],
                candidate[1],
                "Candidate",
                None,
                None,
                "DERIVED",
                "test",
                "test",
            )
            for protein_id in discovery_targets
        ]
        + [
            FeatureRecord(
                protein_id,
                nonsignificant[0],
                nonsignificant[1],
                "Nonsignificant",
                None,
                None,
                "DERIVED",
                "test",
                "test",
            )
            for protein_id in (*discovery_targets[:3], *discovery_backgrounds[:3])
        ]
    )
    universes = {candidate: set(protein_ids), nonsignificant: set(protein_ids)}
    rows = analyse_feature_associations(
        comparisons=(_comparison(),),
        label_memberships=memberships,
        partitions=partitions,
        features=features,
        settings=_settings(minimum_target=2, minimum_background=2),
        feature_assessment_universes=universes,
    )
    validation_rows = [row for row in rows if row.partition == "VALIDATION"]
    assert {(row.feature_type, row.feature_id) for row in validation_rows} == {candidate}
    assert validation_rows[0].status == AnalysisStatus.COMPLETE
    assert validation_rows[0].target_units_with_feature == 0
    assert validation_rows[0].background_units_with_feature == 0
    assert validation_rows[0].p_value == pytest.approx(1.0)


def test_signature_synthesis_covers_each_validation_class() -> None:
    """Discovery evidence should be classified against held-out direction and FDR."""

    discovery = _association(partition="DISCOVERY", q_value=0.01, prevalence_difference=0.8)
    assert (
        _classify_validation(
            discovery=discovery,
            validation=None,
            fdr_threshold=0.05,
        )
        == "NO_VALIDATION_DATA"
    )
    incomplete = _association(
        partition="VALIDATION",
        q_value=None,
        prevalence_difference=None,
        status=AnalysisStatus.INSUFFICIENT_SAMPLE_SIZE,
    )
    assert (
        _classify_validation(
            discovery=discovery,
            validation=incomplete,
            fdr_threshold=0.05,
        )
        == "NO_VALIDATION_DATA"
    )
    discordant = _association(partition="VALIDATION", q_value=0.01, prevalence_difference=-0.5)
    assert (
        _classify_validation(
            discovery=discovery,
            validation=discordant,
            fdr_threshold=0.05,
        )
        == "DIRECTION_DISCORDANT"
    )
    validated = _association(partition="VALIDATION", q_value=0.05, prevalence_difference=0.5)
    assert (
        _classify_validation(
            discovery=discovery,
            validation=validated,
            fdr_threshold=0.05,
        )
        == "VALIDATED_STUDY_WIDE"
    )
    within_comparison = _association(
        partition="VALIDATION",
        q_value=0.05,
        study_q_value=0.2,
        prevalence_difference=0.5,
    )
    assert (
        _classify_validation(
            discovery=discovery,
            validation=within_comparison,
            fdr_threshold=0.05,
        )
        == "VALIDATED_WITHIN_COMPARISON"
    )
    weak = _association(partition="VALIDATION", q_value=0.2, prevalence_difference=0.5)
    assert (
        _classify_validation(
            discovery=discovery,
            validation=weak,
            fdr_threshold=0.05,
        )
        == "DISCOVERY_ONLY"
    )
    summaries = summarise_signatures(
        associations=(discovery, validated),
        fdr_threshold=0.05,
    )
    assert summaries[0].evidence_class == "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE"
    local_discovery = _association(
        partition="DISCOVERY",
        q_value=0.01,
        study_q_value=0.2,
    )
    assert (
        _classify_discovery_scope(
            discovery=local_discovery,
            validation=None,
            fdr_threshold=0.05,
        )
        == "EXPLORATORY_LOCAL_ONLY__NO_VALIDATION_DATA"
    )
    exploratory = summarise_signatures(
        associations=(discovery, validated),
        fdr_threshold=0.05,
        exploratory_feature_keys={("TYPE", "feature")},
    )
    assert exploratory[0].evidence_class == (
        "EXPLORATORY_ALL_DATA_DERIVATION__VALIDATED_STUDY_WIDE"
    )
    technical_discovery = _association(
        partition="DISCOVERY",
        q_value=0.01,
        feature_type="STRUCTURE_AVAILABLE",
    )
    assert (
        _classify_discovery_scope(
            discovery=technical_discovery,
            validation=None,
            fdr_threshold=0.05,
        )
        == "QC_TECHNICAL_NON_BIOLOGICAL__NO_VALIDATION_DATA"
    )
    no_discovery = summarise_signatures(
        associations=(_association(partition="DISCOVERY", q_value=0.5),),
        fdr_threshold=0.05,
    )
    assert no_discovery[0].status == AnalysisStatus.NO_SIGNIFICANT_SIGNATURE
    not_tested_row = _status_result(
        comparison_id="cmp",
        partition="DISCOVERY",
        target_count=5,
        background_count=5,
        target_unit_count=5,
        background_unit_count=5,
        excluded_mixed_unit_count=0,
        status=AnalysisStatus.NO_ELIGIBLE_FEATURES,
    )
    assert _no_discovery_status(discovery_rows=(not_tested_row,)) == (
        AnalysisStatus.NO_ELIGIBLE_FEATURES
    )
    not_tested = summarise_signatures(
        associations=(not_tested_row,),
        fdr_threshold=0.05,
    )
    assert not_tested[0].status == AnalysisStatus.NO_ELIGIBLE_FEATURES
    assert not_tested[0].evidence_class == "DISCOVERY_NOT_TESTED__NO_ELIGIBLE_FEATURES"
    with pytest.raises(InputValidationError):
        summarise_signatures(associations=(), fdr_threshold=1.1)


def test_structural_features_cover_fold_filtering_and_graph_paths() -> None:
    """Structural signatures should require supported models, folds and pairwise edges."""

    structures = (
        _structure("p1", "s1", availability="AVAILABLE", fold_status="ASSESSED_WITH_HIT"),
        _structure("p2", "s2", availability="FAILED", fold_status="FAILED"),
        _structure("p3", "s3", availability="COMPLETE", fold_status="ASSESSED_WITH_HIT"),
    )
    passing = _structure_comparison("p1", "p2", "one", score=0.7)
    bridge = _structure_comparison("p2", "p3", "two", score=0.8)
    failing = _structure_comparison("p1", "p3", "three", score=0.2)
    features, clusters = derive_structure_features(
        structures=structures,
        comparisons=(passing, bridge, failing),
        tm_score_threshold=0.5,
        minimum_coverage=0.5,
        derivation_cohort_sha256="a" * 64,
    )
    assert len([row for row in features if row.feature_type == "STRUCTURE_CLUSTER"]) == 3
    assert len(clusters) == 3
    assert {row.feature_type for row in features} >= {"STRUCTURE_AVAILABLE", "FOLD"}
    fold_features = [row for row in features if row.feature_type == "FOLD"]
    assert {row.evidence_source for row in fold_features} == {"CATH"}
    assert {row.evidence_reference for row in fold_features} == {"CATH:4.3.0:fold_1"}
    cluster_features = [row for row in features if row.feature_type == "STRUCTURE_CLUSTER"]
    assert {row.derivation_scope for row in cluster_features} == {"ALL_DATA_EXPLORATORY"}
    assert all(row.derivation_cohort_sha256 == "a" * 64 for row in cluster_features)
    assert all(len(row.feature_definition_sha256) == 64 for row in features)
    updated_release, _ = derive_structure_features(
        structures=(replace(structures[0], fold_authority_version="4.4.0"),),
        comparisons=(),
        tm_score_threshold=0.5,
        minimum_coverage=0.5,
    )
    original_fold_digest = next(
        row.feature_definition_sha256 for row in fold_features if row.protein_id == "p1"
    )
    updated_fold_digest = next(
        row.feature_definition_sha256 for row in updated_release if row.feature_type == "FOLD"
    )
    assert updated_fold_digest != original_fold_digest
    with pytest.raises(InputValidationError, match="conflicting authority-release"):
        derive_structure_features(
            structures=(
                structures[0],
                replace(
                    structures[0],
                    protein_id="p4",
                    structure_id="s4",
                    fold_evidence_reference="different-method-reference",
                ),
            ),
            comparisons=(),
            tm_score_threshold=0.5,
            minimum_coverage=0.5,
        )
    with pytest.raises(InputValidationError, match="cohort SHA-256"):
        derive_structure_features(
            structures=(),
            comparisons=(passing,),
            tm_score_threshold=0.5,
            minimum_coverage=0.5,
        )
    assert (
        _comparison_passes(
            comparison=passing,
            tm_score_threshold=0.7,
            minimum_coverage=0.8,
        )
        is True
    )
    assert (
        _comparison_passes(
            comparison=_structure_comparison("p1", "p2", "failed", status="FAILED"),
            tm_score_threshold=0.5,
            minimum_coverage=0.5,
        )
        is False
    )
    cyclic = {"p1": {"p2"}, "p2": {"p1"}}
    assert _connected_components(adjacency=cyclic) == (("p1", "p2"),)
    with pytest.raises(InputValidationError):
        derive_structure_features(
            structures=(),
            comparisons=(),
            tm_score_threshold=-0.1,
            minimum_coverage=0.5,
        )


def test_structural_clusters_are_frozen_before_validation_projection() -> None:
    """Validation-to-validation edges must not define a discovery structural cluster."""

    comparisons = (
        _structure_comparison("d1", "d2", "discovery"),
        _structure_comparison("v1", "d1", "projection"),
        _structure_comparison("v1", "v2", "validation_only"),
    )
    partitions = (
        PartitionAssignment("d1", "DISCOVERY", "GROUP", "d1"),
        PartitionAssignment("d2", "DISCOVERY", "GROUP", "d2"),
        PartitionAssignment("v1", "VALIDATION", "GROUP", "v1"),
        PartitionAssignment("v2", "VALIDATION", "GROUP", "v2"),
    )
    _, rows = derive_structure_features(
        structures=(),
        comparisons=comparisons,
        tm_score_threshold=0.5,
        minimum_coverage=0.5,
        partitions=partitions,
        derivation_cohort_sha256="b" * 64,
    )
    assert {row["protein_id"] for row in rows} == {"d1", "d2", "v1"}
    projected = next(row for row in rows if row["protein_id"] == "v1")
    assert projected["membership_method"] == "VALIDATION_PROJECTION"
    assert projected["reference_member_count"] == 2
    assert all(row["reference_partition"] == "DISCOVERY" for row in rows)


def test_structural_clusters_do_not_bridge_comparison_universes_or_tools() -> None:
    """Disconnected evidence campaigns and tools must define separate components."""

    comparisons = (
        _structure_comparison("p1", "p2", "u1_edge", universe="universe_1"),
        _structure_comparison("p2", "p3", "u2_edge", universe="universe_2"),
        _structure_comparison(
            "p3",
            "p4",
            "other_tool_edge",
            universe="universe_2",
            tool="TM-align",
        ),
    )
    features, rows = derive_structure_features(
        structures=(),
        comparisons=comparisons,
        tm_score_threshold=0.5,
        minimum_coverage=0.5,
        derivation_cohort_sha256="c" * 64,
    )
    cluster_members = {
        cluster_id: {str(row["protein_id"]) for row in rows if row["cluster_id"] == cluster_id}
        for cluster_id in {str(row["cluster_id"]) for row in rows}
    }
    assert set(map(frozenset, cluster_members.values())) == {
        frozenset({"p1", "p2"}),
        frozenset({"p2", "p3"}),
        frozenset({"p3", "p4"}),
    }
    assert len({item.feature_id for item in features}) == 3
    assert {str(row["coverage_scope"]) for row in rows} == {"FULL_SEQUENCE"}


def test_structural_assessment_universes_require_proved_coverage() -> None:
    """Model, fold and cluster absences should use only explicit assessment ledgers."""

    eligible = _structure("p1", "s1", availability="AVAILABLE", fold_status="ASSESSED_WITH_HIT")
    low_confidence = replace(
        eligible,
        protein_id="p2",
        structure_id="s2",
        coordinate_path=Path("/models/s2.pdb"),
        mean_confidence=20.0,
        fold_id="",
        fold_name="",
        fold_evidence_status=FoldEvidenceStatus.ASSESSED_NO_HIT,
        analysis_eligibility_status=StructureAnalysisEligibility.INELIGIBLE_LOW_CONFIDENCE,
    )
    external_fold = replace(
        eligible,
        protein_id="p3",
        structure_id="s3",
        coordinate_path=None,
        coordinate_sha256="",
        availability_status="INPUT_UNAVAILABLE",
        analysis_eligibility_status=(StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE),
    )
    external_no_hit = replace(
        external_fold,
        protein_id="p4",
        structure_id="s4",
        fold_id="",
        fold_name="",
        fold_evidence_status=FoldEvidenceStatus.ASSESSED_NO_HIT,
    )
    features, clusters = derive_structure_features(
        structures=(eligible, low_confidence, external_fold, external_no_hit),
        comparisons=(_structure_comparison("p1", "p3", "edge", universe="U1"),),
        tm_score_threshold=0.5,
        minimum_coverage=0.5,
        derivation_cohort_sha256="d" * 64,
    )
    universes = derive_structure_assessment_universes(
        structures=(eligible, low_confidence, external_fold, external_no_hit),
        structure_clusters=clusters,
        comparison_universe_members={"U1": frozenset({"p1", "p2", "p3", "p4"})},
        alphafold_acquisitions=(
            AlphaFoldAcquisition(
                "p5",
                "P5",
                "MODEL_NOT_AVAILABLE",
                "",
                "",
                "https://alphafold.ebi.ac.uk/api/prediction/P5",
                None,
                "",
                None,
                None,
                "not available",
            ),
            AlphaFoldAcquisition(
                "p6",
                "P6",
                "FAILED",
                "",
                "",
                "https://alphafold.ebi.ac.uk/api/prediction/P6",
                None,
                "",
                None,
                None,
                "network failure",
            ),
        ),
    )
    assert universes[("STRUCTURE_AVAILABLE", "SOURCE")] == frozenset({"p1", "p2"})
    assert universes[("STRUCTURE_AVAILABLE", "AlphaFoldDB")] == frozenset({"p5"})
    assert universes[("FOLD", "CATH:fold_1")] == frozenset({"p1", "p3", "p4"})
    cluster_id = next(
        item.feature_id for item in features if item.feature_type == "STRUCTURE_CLUSTER"
    )
    assert universes[("STRUCTURE_CLUSTER", cluster_id)] == frozenset({"p1", "p2", "p3", "p4"})
    declared_structures = tuple(
        replace(item, comparison_universe_ids=("U1",))
        for item in (eligible, low_confidence, external_fold, external_no_hit)
    )
    declared = derive_structure_assessment_universes(
        structures=declared_structures,
        structure_clusters=clusters,
    )
    assert declared[("STRUCTURE_CLUSTER", cluster_id)] == frozenset({"p1", "p2", "p3", "p4"})
    unproved = derive_structure_assessment_universes(
        structures=(eligible,),
        structure_clusters=clusters,
    )
    assert ("STRUCTURE_CLUSTER", cluster_id) not in unproved
    with pytest.raises(InputValidationError, match="non-empty frozenset"):
        derive_structure_assessment_universes(
            structures=(),
            structure_clusters=(),
            comparison_universe_members={"U1": frozenset()},
        )
    with pytest.raises(InputValidationError, match="conflicts"):
        derive_structure_assessment_universes(
            structures=declared_structures,
            structure_clusters=clusters,
            comparison_universe_members={"U1": frozenset({"p1", "p3"})},
        )


def _sequence(*, protein_id: str, sequence: str) -> SequenceRecord:
    """Create a compact sequence record for unit tests."""

    return SequenceRecord(protein_id, "", sequence, len(sequence), protein_id * 32)


def _feature(*, protein_id: str, feature_id: str, name: str) -> FeatureRecord:
    """Create a categorical feature for unit tests."""

    return FeatureRecord(protein_id, "TYPE", feature_id, name, None, None, "DERIVED", "test", "x")


def _comparison() -> ComparisonDefinition:
    """Create one target-versus-background definition."""

    return ComparisonDefinition("cmp", "Comparison", ("target",), ("background",), "")


def _settings(
    *, minimum_target: int = 1, minimum_background: int = 1, minimum_feature: int = 1
) -> AnalysisSettings:
    """Create permissive analysis settings for unit tests."""

    return AnalysisSettings(
        (2,),
        minimum_target,
        minimum_background,
        minimum_feature,
        100,
        0.2,
        0.05,
        1,
        0.5,
        0.5,
    )


def _association(
    *,
    partition: str,
    q_value: float | None,
    study_q_value: float | None = None,
    prevalence_difference: float | None = 0.8,
    status: AnalysisStatus = AnalysisStatus.COMPLETE,
    feature_type: str = "TYPE",
) -> AssociationResult:
    """Create one association result for synthesis tests."""

    return AssociationResult(
        comparison_id="cmp",
        partition=partition,
        analysis_unit="INDEPENDENCE_BLOCK",
        feature_type=feature_type,
        feature_id="feature",
        feature_name="Feature",
        target_protein_count=10,
        background_protein_count=10,
        target_with_feature=8,
        background_with_feature=1,
        target_unit_count=10,
        background_unit_count=10,
        target_assessed_unit_count=10,
        background_assessed_unit_count=10,
        target_unknown_unit_count=0,
        background_unknown_unit_count=0,
        target_units_with_feature=8,
        background_units_with_feature=1,
        excluded_mixed_unit_count=0,
        target_prevalence=0.8,
        background_prevalence=0.1,
        target_prevalence_ci_lower=0.49,
        target_prevalence_ci_upper=0.94,
        background_prevalence_ci_lower=0.02,
        background_prevalence_ci_upper=0.40,
        prevalence_ci_method="WILSON_95_PERCENT_INDEPENDENCE_BLOCKS",
        prevalence_difference=prevalence_difference,
        odds_ratio=36.0,
        p_value=0.01,
        q_value=q_value,
        study_q_value=q_value if study_q_value is None else study_q_value,
        direction="TARGET_ENRICHED",
        status=status,
    )


def _structure(
    protein_id: str,
    structure_id: str,
    *,
    availability: str,
    fold_status: str,
) -> StructureRecord:
    """Create one structure record for derivation tests."""

    available = availability in {"AVAILABLE", "COMPLETE"}
    return StructureRecord(
        protein_id,
        structure_id,
        "SOURCE",
        "v1",
        Path(f"/models/{structure_id}.pdb") if available else None,
        "a" * 64 if available else "",
        availability,
        80.0,
        "fold_1",
        "Fold one",
        "CATH",
        "4.3.0",
        "CATH:4.3.0:fold_1",
        FoldEvidenceStatus(fold_status),
        (
            StructureAnalysisEligibility.ELIGIBLE
            if available
            else StructureAnalysisEligibility.INELIGIBLE_VALIDATION_FAILED
        ),
        (),
    )


def _structure_comparison(
    protein_a: str,
    protein_b: str,
    record_id: str,
    *,
    score: float | None = 0.8,
    status: str = "COMPLETE",
    universe: str = "campaign_1",
    tool: str = "Foldseek",
) -> PairwiseStructureComparison:
    """Create one structural comparison for graph tests."""

    return PairwiseStructureComparison(
        protein_a,
        protein_b,
        tool,
        "v1",
        score,
        1.0,
        80,
        0.8,
        0.8,
        StructureComparisonStatus(status),
        record_id,
        universe,
        StructureCoverageScope.FULL_SEQUENCE,
    )
