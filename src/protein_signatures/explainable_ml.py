"""Mandatory group-aware explainable classification over protein evidence features.

The module deliberately keeps prediction separate from hypothesis-tested
association. Feature selection is label-blind, homology/redundancy groups never
cross validation folds, and SHAP explanations use an interventional discovery
background in the fitted linear model's additive log-odds space.
"""

from __future__ import annotations

import hashlib
import logging
import math
import os
import re
import tempfile
import warnings
from collections import Counter, defaultdict
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .assessment import (
    TECHNICAL_FEATURE_TYPES,
    FeatureAssessmentUniverses,
    FeatureKey,
    NormalisedAssessmentUniverses,
    normalise_feature_assessment_universes,
    validate_feature_inference_scopes,
)
from .errors import ConfigurationError, InputValidationError, PublicationError
from .models import (
    ComparisonDefinition,
    ExplainableMLSettings,
    FeatureRecord,
    PartitionAssignment,
)

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class ExplainableMLResult:
    """Canonical tables and status produced by explainable classification."""

    models: tuple[dict[str, Any], ...]
    feature_importance: tuple[dict[str, Any], ...]
    predictions: tuple[dict[str, Any], ...]
    explanations: tuple[dict[str, Any], ...]
    plots: tuple[dict[str, Any], ...]
    plot_assets: tuple[tuple[str, Path], ...]
    status: str
    implementation_version: str
    shap_version: str
    matplotlib_version: str


def run_explainable_models(
    *,
    comparisons: tuple[ComparisonDefinition, ...],
    label_memberships: tuple[dict[str, str], ...],
    partitions: tuple[PartitionAssignment, ...],
    features: tuple[FeatureRecord, ...],
    settings: ExplainableMLSettings,
    random_seed: int,
    feature_assessment_universes: FeatureAssessmentUniverses | None = None,
    exploratory_feature_keys: Collection[FeatureKey] = (),
) -> ExplainableMLResult:
    """Fit independently validated elastic-net models for eligible comparisons.

    Args:
        comparisons: Explicit target-versus-background definitions.
        label_memberships: Expanded reviewed-positive label memberships.
        partitions: Leakage-safe protein partition and group assignments.
        features: Protein-level biological and technical features.
        settings: Validated explainable-model settings.
        random_seed: Stable campaign random seed.
        feature_assessment_universes: Proteins successfully assessed for each
            candidate feature. Proteins outside a supplied universe are unknown,
            not zero. ``None`` is only appropriate when every feature was
            assessed for every model sample.
        exploratory_feature_keys: All-data-derived features excluded from model
            selection to protect held-out evaluation and SHAP interpretation.

    Returns:
        Canonical model, importance, prediction and local-explanation rows.

    Raises:
        ConfigurationError: If required modelling dependencies are unavailable.
        InputValidationError: If memberships or feature names are contradictory.
    """

    if not settings.enabled:
        raise InputValidationError("Explainable ML is a mandatory campaign stage.")
    dependencies = _load_dependencies()
    exploratory_keys = validate_feature_inference_scopes(
        features=features,
        exploratory_feature_keys=exploratory_feature_keys,
        context="explainable ML",
    )
    labels_by_protein: dict[str, set[str]] = defaultdict(set)
    for row in label_memberships:
        labels_by_protein[str(row["protein_id"])].add(str(row["label_id"]))
    partition_by_protein = {item.protein_id: item for item in partitions}
    features_by_protein: dict[str, set[tuple[str, str]]] = defaultdict(set)
    feature_proteins: dict[tuple[str, str], set[str]] = defaultdict(set)
    feature_names: dict[tuple[str, str], str] = {}
    excluded_technical = 0
    for feature in features:
        key = (feature.feature_type, feature.feature_id)
        previous = feature_names.get(key)
        if previous is not None and previous != feature.feature_name:
            raise InputValidationError(f"Conflicting names for ML feature {key!r}.")
        feature_names[key] = feature.feature_name
        if key in exploratory_keys:
            continue
        if settings.exclude_technical_features and feature.feature_type in (
            TECHNICAL_FEATURE_TYPES
        ):
            excluded_technical += 1
            continue
        features_by_protein[feature.protein_id].add(key)
        feature_proteins[key].add(feature.protein_id)
    assessment_universes = normalise_feature_assessment_universes(
        universes=feature_assessment_universes,
        known_protein_ids=partition_by_protein,
        feature_proteins=feature_proteins,
    )
    if assessment_universes is None:
        LOGGER.warning(
            "No feature-assessment universes were supplied; explainable models "
            "assume every partitioned protein was assessed for every feature."
        )
    model_rows: list[dict[str, Any]] = []
    importance_rows: list[dict[str, Any]] = []
    prediction_rows: list[dict[str, Any]] = []
    explanation_rows: list[dict[str, Any]] = []
    plot_rows: list[dict[str, Any]] = []
    plot_assets: dict[str, Path] = {}
    for comparison_index, comparison in enumerate(comparisons):
        result = _fit_one_comparison(
            comparison=comparison,
            labels_by_protein=labels_by_protein,
            partition_by_protein=partition_by_protein,
            features_by_protein=features_by_protein,
            feature_names=feature_names,
            settings=settings,
            random_seed=random_seed + comparison_index,
            excluded_technical_feature_rows=excluded_technical,
            dependencies=dependencies,
            feature_assessment_universes=assessment_universes,
        )
        model_rows.append(result[0])
        importance_rows.extend(result[1])
        prediction_rows.extend(result[2])
        explanation_rows.extend(result[3])
        if result[4] is not None:
            rendered_rows, rendered_assets = _render_shap_graphics(
                comparison_id=comparison.comparison_id,
                plot_context=result[4],
                cache_dir=settings.plot_cache_dir,
                maximum_waterfalls=settings.maximum_waterfall_plots,
                maximum_display=settings.top_local_explanations,
                dependencies=dependencies,
            )
            plot_rows.extend(rendered_rows)
            for relative_path, source_path in rendered_assets:
                if relative_path in plot_assets:
                    raise PublicationError(f"Duplicate SHAP plot asset path: {relative_path!r}")
                plot_assets[relative_path] = source_path
    complete = sum(row["status"].startswith("COMPLETE") for row in model_rows)
    status = "COMPLETE" if complete else "NO_ELIGIBLE_MODELS"
    return ExplainableMLResult(
        models=tuple(sorted(model_rows, key=lambda row: str(row["comparison_id"]))),
        feature_importance=tuple(
            sorted(
                importance_rows,
                key=lambda row: (str(row["comparison_id"]), int(row["importance_rank"])),
            )
        ),
        predictions=tuple(
            sorted(
                prediction_rows,
                key=lambda row: (
                    str(row["comparison_id"]),
                    str(row["partition"]),
                    str(row["protein_id"]),
                ),
            )
        ),
        explanations=tuple(
            sorted(
                explanation_rows,
                key=lambda row: (
                    str(row["comparison_id"]),
                    str(row["protein_id"]),
                    int(row["absolute_rank"]),
                ),
            )
        ),
        plots=tuple(
            sorted(
                plot_rows,
                key=lambda row: (
                    str(row["comparison_id"]),
                    str(row["plot_type"]),
                    str(row["protein_id"]),
                    str(row["file_format"]),
                ),
            )
        ),
        plot_assets=tuple(sorted(plot_assets.items())),
        status=status,
        implementation_version=str(dependencies["sklearn_version"]),
        shap_version=str(dependencies["shap_version"]),
        matplotlib_version=str(dependencies["matplotlib_version"]),
    )


def _fit_one_comparison(
    *,
    comparison: ComparisonDefinition,
    labels_by_protein: dict[str, set[str]],
    partition_by_protein: dict[str, PartitionAssignment],
    features_by_protein: dict[str, set[tuple[str, str]]],
    feature_names: dict[tuple[str, str], str],
    settings: ExplainableMLSettings,
    random_seed: int,
    excluded_technical_feature_rows: int,
    dependencies: dict[str, Any],
    feature_assessment_universes: NormalisedAssessmentUniverses | None,
) -> tuple[
    dict[str, Any],
    tuple[dict[str, Any], ...],
    tuple[dict[str, Any], ...],
    tuple[dict[str, Any], ...],
    dict[str, Any] | None,
]:
    """Fit and explain one binary comparison without crossing group boundaries."""

    target_labels = set(comparison.target_label_ids)
    background_labels = set(comparison.background_label_ids)
    samples: dict[str, int] = {}
    for protein_id, assignment in sorted(partition_by_protein.items()):
        labels = labels_by_protein.get(protein_id, set())
        is_target = bool(labels & target_labels)
        is_background = bool(labels & background_labels)
        if is_target and is_background:
            raise InputValidationError(
                f"Comparison {comparison.comparison_id!r} has an ML sample in both classes: "
                f"{protein_id!r}."
            )
        if is_target or is_background:
            samples[protein_id] = 1 if is_target else 0
    discovery = tuple(
        protein_id
        for protein_id in sorted(samples)
        if partition_by_protein[protein_id].partition == "DISCOVERY"
    )
    validation = tuple(
        protein_id
        for protein_id in sorted(samples)
        if partition_by_protein[protein_id].partition == "VALIDATION"
    )
    discovery_group_classes: dict[str, set[int]] = defaultdict(set)
    for protein_id in discovery:
        discovery_group_classes[partition_by_protein[protein_id].partition_key].add(
            samples[protein_id]
        )
    mixed_discovery_groups = frozenset(
        group for group, classes in discovery_group_classes.items() if len(classes) > 1
    )
    model_discovery = tuple(
        protein_id
        for protein_id in discovery
        if partition_by_protein[protein_id].partition_key not in mixed_discovery_groups
    )
    class_counts = Counter(samples[protein_id] for protein_id in model_discovery)
    group_counts = {
        class_value: len(
            {
                partition_by_protein[protein_id].partition_key
                for protein_id in model_discovery
                if samples[protein_id] == class_value
            }
        )
        for class_value in (0, 1)
    }
    base = _model_row(
        comparison_id=comparison.comparison_id,
        discovery=model_discovery,
        validation=validation,
        samples=samples,
        partition_by_protein=partition_by_protein,
        settings=settings,
        excluded_technical_feature_rows=excluded_technical_feature_rows,
    )
    if any(class_counts[value] < settings.minimum_samples_per_class for value in (0, 1)):
        return (
            {
                **base,
                "status": "INSUFFICIENT_SAMPLE_SIZE",
                "status_message": (
                    "Discovery partition lacks the configured samples per class after "
                    f"excluding {len(mixed_discovery_groups)} mixed-class blocks."
                ),
            },
            (),
            (),
            (),
            None,
        )
    if any(group_counts[value] < settings.minimum_groups_per_class for value in (0, 1)):
        return (
            {
                **base,
                "status": "INSUFFICIENT_INDEPENDENT_GROUPS",
                "status_message": (
                    "Discovery partition lacks pure independent groups per class after "
                    f"excluding {len(mixed_discovery_groups)} mixed-class blocks."
                ),
            },
            (),
            (),
            (),
            None,
        )
    all_ids = (*discovery, *validation)
    all_id_set = frozenset(all_ids)
    eligible_feature_keys = (
        None
        if feature_assessment_universes is None
        else frozenset(
            key
            for key, assessed_proteins in feature_assessment_universes.items()
            if all_id_set <= assessed_proteins
        )
    )
    feature_keys = _select_feature_keys(
        protein_ids=model_discovery,
        features_by_protein=features_by_protein,
        maximum_features=settings.maximum_features,
        eligible_feature_keys=eligible_feature_keys,
    )
    if not feature_keys:
        return (
            {
                **base,
                "status": "NO_ELIGIBLE_FEATURES",
                "status_message": (
                    "No non-constant discovery features with complete assessment "
                    "across every prediction sample were available."
                    if feature_assessment_universes is not None
                    else "No non-constant discovery features were available."
                ),
            },
            (),
            (),
            (),
            None,
        )
    np = dependencies["numpy"]
    x_discovery = _feature_matrix(
        protein_ids=model_discovery,
        feature_keys=feature_keys,
        features_by_protein=features_by_protein,
        numpy_module=np,
    )
    y_discovery = np.asarray([samples[item] for item in model_discovery], dtype=int)
    discovery_groups = np.asarray(
        [partition_by_protein[item].partition_key for item in model_discovery], dtype=object
    )
    discovery_block_matrix, discovery_block_labels, unexpected_mixed_blocks = (
        _aggregate_pure_blocks(
            matrix=x_discovery,
            labels=y_discovery,
            groups=discovery_groups,
            numpy_module=np,
        )
    )
    if unexpected_mixed_blocks:
        raise InputValidationError(
            "Mixed discovery blocks remained after the pre-fit exclusion safeguard."
        )
    selected_strength, cv_score, completed_folds = _select_regularisation(
        matrix=x_discovery,
        labels=y_discovery,
        groups=discovery_groups,
        settings=settings,
        random_seed=random_seed,
        dependencies=dependencies,
    )
    if selected_strength is None:
        return (
            {
                **base,
                "status": "CROSS_VALIDATION_FAILED",
                "status_message": "No group-aware fold contained both classes.",
                "feature_count": len(feature_keys),
            },
            (),
            (),
            (),
            None,
        )
    model, converged = _fit_model(
        matrix=x_discovery,
        labels=y_discovery,
        settings=settings,
        strength=selected_strength,
        random_seed=random_seed,
        dependencies=dependencies,
        sample_weights=_group_balanced_sample_weights(
            labels=y_discovery,
            groups=discovery_groups,
            numpy_module=np,
        ),
    )
    x_all = _feature_matrix(
        protein_ids=all_ids,
        feature_keys=feature_keys,
        features_by_protein=features_by_protein,
        numpy_module=np,
    )
    probabilities = model.predict_proba(x_all)[:, 1]
    predictions = tuple(
        {
            "comparison_id": comparison.comparison_id,
            "protein_id": protein_id,
            "partition": partition_by_protein[protein_id].partition,
            "partition_key": partition_by_protein[protein_id].partition_key,
            "true_class": "TARGET" if samples[protein_id] else "BACKGROUND",
            "predicted_probability": float(probability),
            "predicted_class": "TARGET" if probability >= 0.5 else "BACKGROUND",
            "correct": bool((probability >= 0.5) == bool(samples[protein_id])),
        }
        for protein_id, probability in zip(all_ids, probabilities, strict=True)
    )
    x_validation = x_all[len(discovery) :]
    y_validation = np.asarray([samples[item] for item in validation], dtype=int)
    validation_groups = np.asarray(
        [partition_by_protein[item].partition_key for item in validation], dtype=object
    )
    validation_probabilities = probabilities[len(discovery) :]
    validation_block_matrix, validation_block_labels, mixed_validation_blocks = (
        _aggregate_pure_blocks(
            matrix=x_validation,
            labels=y_validation,
            groups=validation_groups,
            numpy_module=np,
        )
    )
    validation_group_classes: dict[str, set[int]] = defaultdict(set)
    for label, group in zip(y_validation.tolist(), validation_groups.tolist(), strict=True):
        validation_group_classes[str(group)].add(int(label))
    pure_validation_groups = frozenset(
        group for group, classes in validation_group_classes.items() if len(classes) == 1
    )
    validation_sample_counts = Counter(
        int(label)
        for label, group in zip(y_validation.tolist(), validation_groups.tolist(), strict=True)
        if str(group) in pure_validation_groups
    )
    validation_block_counts = Counter(validation_block_labels.tolist())
    validation_is_adequate = bool(validation) and all(
        validation_sample_counts[class_value] >= settings.minimum_samples_per_class
        and validation_block_counts[class_value] >= settings.minimum_groups_per_class
        for class_value in (0, 1)
    )
    if validation_is_adequate:
        validation_block_probabilities = model.predict_proba(validation_block_matrix)[:, 1]
        metrics = _classification_metrics(
            labels=validation_block_labels,
            probabilities=validation_block_probabilities,
            dependencies=dependencies,
        )
    else:
        validation_block_probabilities = np.asarray([], dtype=float)
        metrics = _classification_metrics(
            labels=np.asarray([], dtype=int),
            probabilities=validation_block_probabilities,
            dependencies=dependencies,
        )
    coefficients = model.coef_[0]
    permutation = (
        _held_out_permutation_importance(
            matrix=validation_block_matrix,
            labels=validation_block_labels,
            coefficients=coefficients,
            intercept=float(model.intercept_[0]),
            repeats=settings.permutation_repeats,
            maximum_features=settings.maximum_permutation_features,
            random_seed=random_seed,
            dependencies=dependencies,
        )
        if validation_is_adequate
        else {}
    )
    means = np.asarray(discovery_block_matrix.mean(axis=0), dtype=float)
    explanation_ids = validation or discovery
    explanation_partition = "VALIDATION" if validation else "DISCOVERY"
    explanation_matrix = x_validation if validation else x_all[: len(discovery)]
    explanation_labels = (
        y_validation if validation else np.asarray([samples[item] for item in discovery], dtype=int)
    )
    explanation_probabilities = (
        validation_probabilities if validation else probabilities[: len(discovery)]
    )
    shap_explanation = _build_shap_explanation(
        model=model,
        background=discovery_block_matrix,
        matrix=explanation_matrix,
        feature_keys=feature_keys,
        feature_names=feature_names,
        dependencies=dependencies,
    )
    validation_contributions = (
        (x_validation - means) * coefficients if len(validation) else np.empty((0, len(means)))
    )
    target_prevalence = discovery_block_matrix[discovery_block_labels == 1].mean(axis=0)
    background_prevalence = discovery_block_matrix[discovery_block_labels == 0].mean(axis=0)
    order = sorted(
        range(len(feature_keys)),
        key=lambda index: (-abs(float(coefficients[index])), feature_keys[index]),
    )
    importance = tuple(
        {
            "comparison_id": comparison.comparison_id,
            "feature_type": feature_keys[index][0],
            "feature_id": feature_keys[index][1],
            "feature_name": feature_names[feature_keys[index]],
            "coefficient_log_odds": float(coefficients[index]),
            "odds_multiplier": float(math.exp(max(-700.0, min(700.0, coefficients[index])))),
            "discovery_target_prevalence": float(target_prevalence[index]),
            "discovery_background_prevalence": float(background_prevalence[index]),
            "discovery_prevalence_difference": float(
                target_prevalence[index] - background_prevalence[index]
            ),
            "validation_permutation_importance_mean": permutation.get(index, (None, None))[0],
            "validation_permutation_importance_stddev": permutation.get(index, (None, None))[1],
            "mean_absolute_validation_contribution": (
                float(np.abs(validation_contributions[:, index]).mean())
                if len(validation)
                else None
            ),
            "importance_rank": rank,
        }
        for rank, index in enumerate(order, start=1)
    )
    explanations = _local_explanations(
        comparison_id=comparison.comparison_id,
        partition=explanation_partition,
        protein_ids=explanation_ids,
        labels=explanation_labels,
        probabilities=explanation_probabilities,
        matrix=explanation_matrix,
        means=means,
        coefficients=coefficients,
        shap_values=shap_explanation.values,
        feature_keys=feature_keys,
        feature_names=feature_names,
        maximum_per_protein=settings.top_local_explanations,
        numpy_module=np,
    )
    training_note = (
        "Discovery fitting used inverse-block-size and class-balanced sample weights over "
        f"{len(discovery_block_labels)} pure blocks; {len(mixed_discovery_groups)} "
        "mixed-class blocks were excluded. "
    )
    if not validation:
        status = "COMPLETE_NO_VALIDATION"
        status_message = (
            training_note + "No held-out samples were available; predictions and SHAP explanations "
            "use discovery samples and no held-out metrics were calculated."
        )
    elif not validation_is_adequate:
        status = "COMPLETE_VALIDATION_UNDERPOWERED"
        status_message = (
            training_note
            + "Held-out metrics and permutation importance were withheld because validation "
            "did not meet both per-class thresholds after mixed blocks were excluded: "
            f"target proteins={validation_sample_counts[1]}, background proteins="
            f"{validation_sample_counts[0]}, target pure blocks={validation_block_counts[1]}, "
            f"background pure blocks={validation_block_counts[0]}, mixed blocks="
            f"{mixed_validation_blocks}; predictions and SHAP explanations remain descriptive."
        )
    else:
        status = "COMPLETE"
        status_message = (
            training_note
            + "Held-out metrics and permutation importance use pure independence blocks "
            f"(target={validation_block_counts[1]}, background={validation_block_counts[0]}, "
            f"mixed blocks excluded={mixed_validation_blocks})."
        )
    if not converged:
        status += "_CONVERGENCE_WARNING"
    model_row = {
        **base,
        "status": status,
        "status_message": status_message,
        "feature_count": len(feature_keys),
        "cross_validation_folds": completed_folds,
        "regularisation_strength": selected_strength,
        "l1_ratio": settings.l1_ratio,
        "cv_roc_auc": cv_score,
        "shap_background_partition": "DISCOVERY",
        "shap_explained_partition": explanation_partition,
        "shap_explained_sample_count": len(explanation_ids),
        "shap_expected_value_log_odds": float(shap_explanation.base_values[0]),
        **metrics,
    }
    plot_context = {
        "explanation": shap_explanation,
        "protein_ids": explanation_ids,
        "partition": explanation_partition,
        "probabilities": explanation_probabilities,
        "feature_count": len(feature_keys),
    }
    return model_row, importance, predictions, explanations, plot_context


def _model_row(
    *,
    comparison_id: str,
    discovery: tuple[str, ...],
    validation: tuple[str, ...],
    samples: dict[str, int],
    partition_by_protein: dict[str, PartitionAssignment],
    settings: ExplainableMLSettings,
    excluded_technical_feature_rows: int,
) -> dict[str, Any]:
    """Create a schema-complete model status row before fitting."""

    return {
        "comparison_id": comparison_id,
        "model_type": "ELASTIC_NET_LOGISTIC_REGRESSION",
        "status": "NOT_RUN",
        "status_message": "",
        "analysis_unit": "INVERSE_BLOCK_WEIGHTED_PROTEIN_FIT_BLOCK_LEVEL_VALIDATION",
        "discovery_target_count": sum(samples[item] == 1 for item in discovery),
        "discovery_background_count": sum(samples[item] == 0 for item in discovery),
        "discovery_group_count": len(
            {partition_by_protein[item].partition_key for item in discovery}
        ),
        "validation_target_count": sum(samples[item] == 1 for item in validation),
        "validation_background_count": sum(samples[item] == 0 for item in validation),
        "validation_group_count": len(
            {partition_by_protein[item].partition_key for item in validation}
        ),
        "feature_count": 0,
        "excluded_technical_feature_rows": excluded_technical_feature_rows,
        "cross_validation_folds": 0,
        "regularisation_strength": None,
        "l1_ratio": settings.l1_ratio,
        "cv_roc_auc": None,
        "validation_roc_auc": None,
        "validation_average_precision": None,
        "validation_balanced_accuracy": None,
        "validation_matthews_correlation": None,
        "validation_brier_score": None,
        "decision_threshold": 0.5,
        "explanation_method": "SHAP_LINEAR_INTERVENTIONAL_LOG_ODDS",
        "shap_background_partition": "DISCOVERY",
        "shap_explained_partition": "",
        "shap_explained_sample_count": 0,
        "shap_expected_value_log_odds": None,
    }


def _select_feature_keys(
    *,
    protein_ids: tuple[str, ...],
    features_by_protein: dict[str, set[tuple[str, str]]],
    maximum_features: int,
    eligible_feature_keys: Collection[tuple[str, str]] | None = None,
) -> tuple[tuple[str, str], ...]:
    """Select non-constant features without inspecting class labels.

    Features are ranked by discovery document frequency within type, then taken
    round-robin across types so abundant k-mers cannot crowd out Pfam, fold or
    pocket evidence.

    Args:
        protein_ids: Discovery proteins used for label-blind prevalence.
        features_by_protein: Observed feature-presence mapping.
        maximum_features: Maximum number of selected features.
        eligible_feature_keys: Optional keys known to have complete assessment
            across every model sample.

    Returns:
        Deterministically ordered non-constant feature keys.
    """

    eligible = None if eligible_feature_keys is None else frozenset(eligible_feature_keys)
    prevalence: Counter[tuple[str, str]] = Counter()
    for protein_id in protein_ids:
        observed = features_by_protein.get(protein_id, set())
        prevalence.update(observed if eligible is None else observed & eligible)
    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for key, count in prevalence.items():
        if 0 < count < len(protein_ids):
            grouped[key[0]].append(key)
    for feature_type, keys in grouped.items():
        grouped[feature_type] = sorted(
            keys,
            key=lambda key: (-prevalence[key], key[1]),
        )
    selected: list[tuple[str, str]] = []
    offset = 0
    feature_types = sorted(grouped)
    while len(selected) < maximum_features:
        appended = False
        for feature_type in feature_types:
            keys = grouped[feature_type]
            if offset < len(keys):
                selected.append(keys[offset])
                appended = True
                if len(selected) == maximum_features:
                    break
        if not appended:
            break
        offset += 1
    return tuple(selected)


def _feature_matrix(
    *,
    protein_ids: tuple[str, ...],
    feature_keys: tuple[tuple[str, str], ...],
    features_by_protein: dict[str, set[tuple[str, str]]],
    numpy_module: Any,
) -> Any:
    """Build a deterministic dense binary matrix for a bounded feature set."""

    key_index = {key: index for index, key in enumerate(feature_keys)}
    matrix = numpy_module.zeros((len(protein_ids), len(feature_keys)), dtype=float)
    for row_index, protein_id in enumerate(protein_ids):
        for key in features_by_protein.get(protein_id, set()):
            column = key_index.get(key)
            if column is not None:
                matrix[row_index, column] = 1.0
    return matrix


def _aggregate_pure_blocks(
    *,
    matrix: Any,
    labels: Any,
    groups: Any,
    numpy_module: Any,
) -> tuple[Any, Any, int]:
    """Aggregate protein rows into pure independence-block feature means.

    Args:
        matrix: Protein-by-feature matrix.
        labels: Binary protein labels.
        groups: Independence-block identifiers aligned with rows.
        numpy_module: Imported NumPy module.

    Returns:
        Block-by-feature matrix, block labels and number of mixed-class blocks
        excluded from the result.

    Raises:
        InputValidationError: If shapes, labels or group identifiers are invalid.
    """

    values = numpy_module.asarray(matrix, dtype=float)
    label_values = numpy_module.asarray(labels, dtype=int)
    group_values = numpy_module.asarray(groups, dtype=object)
    if values.ndim != 2:
        raise InputValidationError("Block aggregation requires a two-dimensional matrix.")
    if label_values.ndim != 1 or group_values.ndim != 1:
        raise InputValidationError("Block labels and groups must be one-dimensional.")
    if not (len(values) == len(label_values) == len(group_values)):
        raise InputValidationError("Block matrix, labels and groups must have matching rows.")
    if not set(label_values.tolist()) <= {0, 1}:
        raise InputValidationError("Block aggregation labels must be binary zero/one values.")
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, raw_group in enumerate(group_values.tolist()):
        if not isinstance(raw_group, str) or not raw_group.strip():
            raise InputValidationError("Independence-block identifiers must be non-empty strings.")
        grouped[raw_group].append(index)
    rows: list[Any] = []
    block_labels: list[int] = []
    mixed_count = 0
    for group in sorted(grouped):
        indices = grouped[group]
        classes = set(label_values[indices].tolist())
        if len(classes) != 1:
            mixed_count += 1
            continue
        rows.append(values[indices].mean(axis=0))
        block_labels.append(int(next(iter(classes))))
    block_matrix = (
        numpy_module.vstack(rows) if rows else numpy_module.empty((0, values.shape[1]), dtype=float)
    )
    return block_matrix, numpy_module.asarray(block_labels, dtype=int), mixed_count


def _group_balanced_sample_weights(*, labels: Any, groups: Any, numpy_module: Any) -> Any:
    """Give each pure block equal influence within balanced outcome classes.

    Each protein first receives inverse-block-size weight. Weights are then
    rescaled within each class so target and background have equal total weight,
    while proteins in duplicated or unusually large blocks cannot dominate fit.

    Args:
        labels: Binary protein labels.
        groups: Pure independence-block identifiers aligned with labels.
        numpy_module: Imported NumPy module.

    Returns:
        Positive sample weights summing to the number of protein rows.

    Raises:
        InputValidationError: If inputs are empty, misaligned, non-binary or a
            block contains members from both classes.
    """

    label_values = numpy_module.asarray(labels, dtype=int)
    group_values = numpy_module.asarray(groups, dtype=object)
    if label_values.ndim != 1 or group_values.ndim != 1 or len(label_values) != len(group_values):
        raise InputValidationError("Sample-weight labels and groups must be aligned vectors.")
    if len(label_values) == 0 or set(label_values.tolist()) != {0, 1}:
        raise InputValidationError("Sample weights require both binary outcome classes.")
    group_sizes = Counter(group_values.tolist())
    if any(not isinstance(group, str) or not group.strip() for group in group_sizes):
        raise InputValidationError("Sample-weight group identifiers must be non-empty strings.")
    group_classes: dict[str, set[int]] = defaultdict(set)
    for label, group in zip(label_values.tolist(), group_values.tolist(), strict=True):
        group_classes[group].add(label)
    mixed = sorted(group for group, classes in group_classes.items() if len(classes) > 1)
    if mixed:
        raise InputValidationError(
            f"Sample weights require pure outcome blocks; mixed blocks: {mixed[:10]}"
        )
    weights = numpy_module.asarray(
        [1.0 / group_sizes[group] for group in group_values.tolist()],
        dtype=float,
    )
    class_target = len(label_values) / 2.0
    for class_value in (0, 1):
        mask = label_values == class_value
        class_total = float(weights[mask].sum())
        if class_total <= 0.0:
            raise InputValidationError("Sample weights require positive class totals.")
        weights[mask] *= class_target / class_total
    return weights


def _select_regularisation(
    *,
    matrix: Any,
    labels: Any,
    groups: Any,
    settings: ExplainableMLSettings,
    random_seed: int,
    dependencies: dict[str, Any],
) -> tuple[float | None, float | None, int]:
    """Choose regularisation through stratified, non-overlapping group folds."""

    group_classes: dict[str, set[int]] = defaultdict(set)
    for label, group in zip(labels.tolist(), groups.tolist(), strict=True):
        group_classes[str(group)].add(int(label))
    if any(len(classes) > 1 for classes in group_classes.values()):
        raise InputValidationError("Cross-validation requires pure outcome blocks.")
    class_groups = [
        sum(classes == {class_value} for classes in group_classes.values())
        for class_value in (0, 1)
    ]
    folds = min(settings.cross_validation_folds, *class_groups)
    if folds < 2:
        return None, None, 0
    splitter = dependencies["StratifiedGroupKFold"](
        n_splits=folds,
        shuffle=True,
        random_state=random_seed,
    )
    scores: dict[float, list[float]] = defaultdict(list)
    for strength in settings.regularisation_strengths:
        for train_index, test_index in splitter.split(matrix, labels, groups):
            if (
                len(set(labels[train_index].tolist())) < 2
                or len(set(labels[test_index].tolist())) < 2
            ):
                continue
            model, _ = _fit_model(
                matrix=matrix[train_index],
                labels=labels[train_index],
                settings=settings,
                strength=strength,
                random_seed=random_seed,
                dependencies=dependencies,
                sample_weights=_group_balanced_sample_weights(
                    labels=labels[train_index],
                    groups=groups[train_index],
                    numpy_module=dependencies["numpy"],
                ),
            )
            block_matrix, block_labels, _ = _aggregate_pure_blocks(
                matrix=matrix[test_index],
                labels=labels[test_index],
                groups=groups[test_index],
                numpy_module=dependencies["numpy"],
            )
            if len(set(block_labels.tolist())) < 2:
                continue
            probability = model.predict_proba(block_matrix)[:, 1]
            scores[strength].append(float(dependencies["roc_auc_score"](block_labels, probability)))
    eligible = {strength: values for strength, values in scores.items() if len(values) == folds}
    if not eligible:
        return None, None, 0
    selected = max(
        eligible,
        key=lambda strength: (
            sum(eligible[strength]) / len(eligible[strength]),
            strength,
        ),
    )
    mean_score = sum(eligible[selected]) / len(eligible[selected])
    return selected, float(mean_score), folds


def _fit_model(
    *,
    matrix: Any,
    labels: Any,
    settings: ExplainableMLSettings,
    strength: float,
    random_seed: int,
    dependencies: dict[str, Any],
    sample_weights: Any | None = None,
) -> tuple[Any, bool]:
    """Fit one deterministic block-aware elastic-net logistic model.

    Args:
        matrix: Protein-by-feature training matrix.
        labels: Binary training labels.
        settings: Validated explainable-model settings.
        strength: Positive regularisation strength.
        random_seed: Stable fitting seed.
        dependencies: Lazily imported modelling dependencies.
        sample_weights: Optional precomputed block/class-balanced row weights.
            When omitted, scikit-learn's class-balanced fallback is retained for
            backwards-compatible helper use.

    Returns:
        Fitted model and whether the optimiser converged before its iteration cap.

    Raises:
        InputValidationError: If supplied sample weights are invalid.
    """

    if sample_weights is not None:
        np = dependencies["numpy"]
        checked_weights = np.asarray(sample_weights, dtype=float)
        if (
            checked_weights.ndim != 1
            or len(checked_weights) != len(labels)
            or not np.all(np.isfinite(checked_weights))
            or np.any(checked_weights <= 0.0)
        ):
            raise InputValidationError(
                "Model sample weights must be finite, positive and aligned to training rows."
            )
    else:
        checked_weights = None

    model = dependencies["LogisticRegression"](
        C=1.0 / strength,
        solver="saga",
        l1_ratio=settings.l1_ratio,
        class_weight=None if checked_weights is not None else "balanced",
        random_state=random_seed,
        max_iter=settings.maximum_iterations,
        tol=settings.convergence_tolerance,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", dependencies["ConvergenceWarning"])
        model.fit(matrix, labels, sample_weight=checked_weights)
    converged = bool(max(model.n_iter_) < settings.maximum_iterations)
    return model, converged


def _classification_metrics(
    *, labels: Any, probabilities: Any, dependencies: dict[str, Any]
) -> dict[str, float | None]:
    """Calculate held-out classification metrics when both classes are present."""

    if len(labels) == 0 or len(set(labels.tolist())) < 2:
        return {
            "validation_roc_auc": None,
            "validation_average_precision": None,
            "validation_balanced_accuracy": None,
            "validation_matthews_correlation": None,
            "validation_brier_score": None,
        }
    predicted = probabilities >= 0.5
    return {
        "validation_roc_auc": float(dependencies["roc_auc_score"](labels, probabilities)),
        "validation_average_precision": float(
            dependencies["average_precision_score"](labels, probabilities)
        ),
        "validation_balanced_accuracy": float(
            dependencies["balanced_accuracy_score"](labels, predicted)
        ),
        "validation_matthews_correlation": float(
            dependencies["matthews_corrcoef"](labels, predicted)
        ),
        "validation_brier_score": float(dependencies["brier_score_loss"](labels, probabilities)),
    }


def _held_out_permutation_importance(
    *,
    matrix: Any,
    labels: Any,
    coefficients: Any,
    intercept: float,
    repeats: int,
    maximum_features: int,
    random_seed: int,
    dependencies: dict[str, Any],
) -> dict[int, tuple[float, float]]:
    """Measure held-out balanced-accuracy loss under column permutation."""

    if len(labels) == 0 or len(set(labels.tolist())) < 2:
        return {}
    np = dependencies["numpy"]
    logits = intercept + matrix @ coefficients
    baseline = float(
        dependencies["balanced_accuracy_score"](labels, _sigmoid(logits, numpy_module=np) >= 0.5)
    )
    selected = sorted(
        range(len(coefficients)),
        key=lambda index: (-abs(float(coefficients[index])), index),
    )[:maximum_features]
    generator = np.random.default_rng(random_seed)
    output: dict[int, tuple[float, float]] = {}
    for index in selected:
        values = matrix[:, index]
        changes: list[float] = []
        for _ in range(repeats):
            permuted = generator.permutation(values)
            changed_logits = logits + coefficients[index] * (permuted - values)
            changed = float(
                dependencies["balanced_accuracy_score"](
                    labels,
                    _sigmoid(changed_logits, numpy_module=np) >= 0.5,
                )
            )
            changes.append(baseline - changed)
        output[index] = (float(np.mean(changes)), float(np.std(changes, ddof=0)))
    return output


def _local_explanations(
    *,
    comparison_id: str,
    partition: str,
    protein_ids: tuple[str, ...],
    labels: Any,
    probabilities: Any,
    matrix: Any,
    means: Any,
    coefficients: Any,
    shap_values: Any,
    feature_keys: tuple[tuple[str, str], ...],
    feature_names: dict[tuple[str, str], str],
    maximum_per_protein: int,
    numpy_module: Any,
) -> tuple[dict[str, Any], ...]:
    """Return top genuine linear SHAP log-odds contributions for explained proteins."""

    if not protein_ids:
        return ()
    contributions = numpy_module.asarray(shap_values, dtype=float)
    if contributions.shape != matrix.shape:
        raise InputValidationError(
            "SHAP contribution matrix does not match the explained feature matrix."
        )
    rows: list[dict[str, Any]] = []
    for row_index, protein_id in enumerate(protein_ids):
        order = numpy_module.argsort(-numpy_module.abs(contributions[row_index]), kind="stable")
        for rank, feature_index in enumerate(order[:maximum_per_protein], start=1):
            index = int(feature_index)
            key = feature_keys[index]
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "protein_id": protein_id,
                    "partition": partition,
                    "true_class": "TARGET" if int(labels[row_index]) else "BACKGROUND",
                    "predicted_probability": float(probabilities[row_index]),
                    "feature_type": key[0],
                    "feature_id": key[1],
                    "feature_name": feature_names[key],
                    "feature_value": int(matrix[row_index, index]),
                    "discovery_background_mean": float(means[index]),
                    "coefficient_log_odds": float(coefficients[index]),
                    "log_odds_contribution": float(contributions[row_index, index]),
                    "absolute_rank": rank,
                    "explanation_method": "SHAP_LINEAR_INTERVENTIONAL_LOG_ODDS",
                }
            )
    return tuple(rows)


def _build_shap_explanation(
    *,
    model: Any,
    background: Any,
    matrix: Any,
    feature_keys: tuple[tuple[str, str], ...],
    feature_names: dict[tuple[str, str], str],
    dependencies: dict[str, Any],
) -> Any:
    """Build a named interventional SHAP explanation in model log-odds space.

    Args:
        model: Fitted scikit-learn linear classifier.
        background: Discovery feature matrix defining the intervention baseline.
        matrix: Discovery or held-out matrix to explain.
        feature_keys: Ordered feature identifiers matching matrix columns.
        feature_names: Human-readable feature names.
        dependencies: Lazily imported numerical and plotting dependencies.

    Returns:
        SHAP ``Explanation`` with stable, unique feature labels.

    Raises:
        InputValidationError: If matrix dimensions do not match the feature contract.
    """

    np = dependencies["numpy"]
    if (
        np.asarray(background).ndim != 2
        or np.asarray(matrix).ndim != 2
        or background.shape[1] != len(feature_keys)
        or matrix.shape[1] != len(feature_keys)
        or len(background) == 0
        or len(matrix) == 0
    ):
        raise InputValidationError(
            "SHAP background and explained matrices must be non-empty, two-dimensional, "
            "and match the selected feature count."
        )
    shap_module = dependencies["shap"]
    masker = shap_module.maskers.Independent(
        background,
        max_samples=len(background),
    )
    explainer = shap_module.LinearExplainer(model, masker=masker)
    raw = explainer(matrix)
    labels = [f"{key[0]} | {feature_names[key]} [{key[1]}]" for key in feature_keys]
    return shap_module.Explanation(
        values=np.asarray(raw.values, dtype=float),
        base_values=np.asarray(raw.base_values, dtype=float),
        data=np.asarray(matrix, dtype=float),
        feature_names=labels,
    )


def _render_shap_graphics(
    *,
    comparison_id: str,
    plot_context: dict[str, Any],
    cache_dir: Path,
    maximum_waterfalls: int,
    maximum_display: int,
    dependencies: dict[str, Any],
) -> tuple[tuple[dict[str, Any], ...], tuple[tuple[str, Path], ...]]:
    """Render global and local SHAP graphics to deterministic cache paths.

    Args:
        comparison_id: Stable comparison identifier.
        plot_context: Explanation, sample identifiers and display metadata.
        cache_dir: Local generated-plot cache.
        maximum_waterfalls: Maximum local waterfall plots to publish.
        maximum_display: Maximum named features displayed by each plot.
        dependencies: Lazily imported plotting dependencies.

    Returns:
        Plot-inventory rows and portable-result-path/source-path pairs.

    Raises:
        InputValidationError: If plotting bounds or context are invalid.
        PublicationError: If a graphic cannot be written atomically.
    """

    if maximum_waterfalls < 0 or maximum_display < 1:
        raise InputValidationError("SHAP plot limits are outside their valid bounds.")
    explanation = plot_context.get("explanation")
    protein_ids = tuple(str(value) for value in plot_context.get("protein_ids", ()))
    partition = str(plot_context.get("partition", ""))
    probabilities = dependencies["numpy"].asarray(
        plot_context.get("probabilities", ()), dtype=float
    )
    feature_count = int(plot_context.get("feature_count", 0))
    if (
        explanation is None
        or not protein_ids
        or len(probabilities) != len(protein_ids)
        or feature_count < 1
    ):
        raise InputValidationError("SHAP plot context is incomplete or inconsistent.")
    shap_module = dependencies["shap"]
    pyplot = dependencies["pyplot"]
    token = _safe_asset_token(value=comparison_id)
    relative_root = Path("analysis") / "06_explainable_models" / "figures" / "shap" / token
    source_root = Path(cache_dir).expanduser().resolve() / token
    rows: list[dict[str, Any]] = []
    assets: list[tuple[str, Path]] = []

    for plot_type, renderer, stem in (
        (
            "SHAP_BEESWARM",
            lambda: shap_module.plots.beeswarm(
                explanation,
                max_display=maximum_display,
                show=False,
            ),
            "beeswarm",
        ),
        (
            "SHAP_GLOBAL_BAR",
            lambda: shap_module.plots.bar(
                explanation,
                max_display=maximum_display,
                show=False,
            ),
            "global_bar",
        ),
    ):
        random_state = dependencies["numpy"].random.get_state()
        plot_seed = int(
            hashlib.sha256(f"{comparison_id}:{plot_type}".encode("utf-8")).hexdigest()[:8],
            16,
        )
        try:
            dependencies["numpy"].random.seed(plot_seed)
            axis = renderer()
        finally:
            dependencies["numpy"].random.set_state(random_state)
        figure = axis.figure
        _register_figure_formats(
            figure=figure,
            comparison_id=comparison_id,
            plot_type=plot_type,
            partition=partition,
            protein_id="",
            stem=stem,
            relative_root=relative_root,
            source_root=source_root,
            sample_count=len(protein_ids),
            feature_count=feature_count,
            rows=rows,
            assets=assets,
        )
        pyplot.close(figure)

    waterfall_order = sorted(
        range(len(protein_ids)),
        key=lambda index: (-abs(float(probabilities[index]) - 0.5), protein_ids[index]),
    )[:maximum_waterfalls]
    for index in waterfall_order:
        protein_id = protein_ids[index]
        axis = shap_module.plots.waterfall(
            explanation[index],
            max_display=maximum_display,
            show=False,
        )
        figure = axis.figure
        _register_figure_formats(
            figure=figure,
            comparison_id=comparison_id,
            plot_type="SHAP_WATERFALL",
            partition=partition,
            protein_id=protein_id,
            stem=f"waterfall_{_safe_asset_token(value=protein_id)}",
            relative_root=relative_root / "waterfalls",
            source_root=source_root / "waterfalls",
            sample_count=1,
            feature_count=feature_count,
            rows=rows,
            assets=assets,
        )
        pyplot.close(figure)
    return tuple(rows), tuple(assets)


def _register_figure_formats(
    *,
    figure: Any,
    comparison_id: str,
    plot_type: str,
    partition: str,
    protein_id: str,
    stem: str,
    relative_root: Path,
    source_root: Path,
    sample_count: int,
    feature_count: int,
    rows: list[dict[str, Any]],
    assets: list[tuple[str, Path]],
) -> None:
    """Save PNG, SVG and PDF forms of one SHAP figure and register every asset.

    Args:
        figure: Matplotlib figure returned by SHAP.
        comparison_id: Stable model comparison.
        plot_type: Controlled SHAP plot type.
        partition: Explained data partition.
        protein_id: Local protein identifier or blank for global plots.
        stem: Safe filename stem.
        relative_root: Portable result-relative parent.
        source_root: Local cache parent.
        sample_count: Number of samples represented.
        feature_count: Number of model features.
        rows: Mutable plot-inventory accumulator.
        assets: Mutable result-path/source-path accumulator.
    """

    for file_format in ("png", "svg", "pdf"):
        relative_path = relative_root / f"{stem}.{file_format}"
        source_path = source_root / f"{stem}.{file_format}"
        _save_figure_atomic(
            figure=figure,
            path=source_path,
            file_format=file_format,
        )
        rows.append(
            {
                "comparison_id": comparison_id,
                "plot_type": plot_type,
                "partition": partition,
                "protein_id": protein_id,
                "file_format": file_format.upper(),
                "asset_path": relative_path.as_posix(),
                "sample_count": sample_count,
                "feature_count": feature_count,
                "explanation_method": "SHAP_LINEAR_INTERVENTIONAL_LOG_ODDS",
            }
        )
        assets.append((relative_path.as_posix(), source_path))


def _save_figure_atomic(*, figure: Any, path: Path, file_format: str) -> None:
    """Save one Matplotlib figure through an atomic same-directory rename.

    Args:
        figure: Matplotlib figure.
        path: Final cache path.
        file_format: ``png``, ``svg`` or ``pdf``.

    Raises:
        PublicationError: If the format is unsupported or writing fails.
    """

    if file_format not in {"png", "svg", "pdf"}:
        raise PublicationError(f"Unsupported SHAP figure format: {file_format!r}")
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=f".{file_format}",
    )
    os.close(descriptor)
    if file_format == "png":
        metadata = {"Software": "protein-signature-analysis"}
    elif file_format == "svg":
        metadata = {"Creator": "protein-signature-analysis", "Date": None}
    else:
        metadata = {
            "Creator": "protein-signature-analysis",
            "CreationDate": None,
            "ModDate": None,
        }
    try:
        figure.savefig(
            temporary_name,
            format=file_format,
            dpi=180,
            bbox_inches="tight",
            metadata=metadata,
        )
        os.replace(temporary_name, destination)
    except (OSError, TypeError, ValueError) as error:
        Path(temporary_name).unlink(missing_ok=True)
        raise PublicationError(
            f"Could not write SHAP {file_format.upper()} graphic {destination}: {error}"
        ) from error


def _safe_asset_token(*, value: str) -> str:
    """Return a bounded, collision-resistant token for a result asset path.

    Args:
        value: Stable scientific identifier.

    Returns:
        Readable filesystem-safe token with a checksum suffix.
    """

    text = str(value).strip()
    if not text:
        raise InputValidationError("A SHAP asset identifier must not be empty.")
    readable = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._-") or "item"
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()[:10]
    return f"{readable[:80]}_{digest}"


def _sigmoid(values: Any, *, numpy_module: Any) -> Any:
    """Calculate a numerically stable vector sigmoid."""

    clipped = numpy_module.clip(values, -709.0, 709.0)
    return 1.0 / (1.0 + numpy_module.exp(-clipped))


def _load_dependencies() -> dict[str, Any]:
    """Load required numerical, SHAP and non-interactive plotting dependencies."""

    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
        import numpy as np
        import shap
        import sklearn
        from sklearn.exceptions import ConvergenceWarning
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import (
            average_precision_score,
            balanced_accuracy_score,
            brier_score_loss,
            matthews_corrcoef,
            roc_auc_score,
        )
        from sklearn.model_selection import StratifiedGroupKFold
    except ImportError as error:
        raise ConfigurationError(
            "Explainable ML and SHAP graphics require the package's required "
            "numerical and plotting dependencies; reinstall protein-signature-analysis."
        ) from error
    matplotlib.rcParams["svg.hashsalt"] = "protein-signature-analysis"
    return {
        "matplotlib_version": matplotlib.__version__,
        "numpy": np,
        "pyplot": plt,
        "shap": shap,
        "shap_version": shap.__version__,
        "sklearn_version": sklearn.__version__,
        "ConvergenceWarning": ConvergenceWarning,
        "LogisticRegression": LogisticRegression,
        "StratifiedGroupKFold": StratifiedGroupKFold,
        "average_precision_score": average_precision_score,
        "balanced_accuracy_score": balanced_accuracy_score,
        "brier_score_loss": brier_score_loss,
        "matthews_corrcoef": matthews_corrcoef,
        "roc_auc_score": roc_auc_score,
    }
