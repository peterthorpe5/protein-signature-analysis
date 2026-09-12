"""Unit tests for the mandatory group-aware explainable prediction layer."""

from __future__ import annotations

import builtins
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from protein_signatures.errors import (
    ConfigurationError,
    InputValidationError,
    PublicationError,
)
from protein_signatures.explainable_ml import (
    _aggregate_pure_blocks,
    _build_shap_explanation,
    _classification_metrics,
    _feature_matrix,
    _fit_model,
    _group_balanced_sample_weights,
    _held_out_permutation_importance,
    _load_dependencies,
    _local_explanations,
    _model_row,
    _render_shap_graphics,
    _safe_asset_token,
    _save_figure_atomic,
    _select_feature_keys,
    _select_regularisation,
    _sigmoid,
    run_explainable_models,
)
from protein_signatures.models import (
    ComparisonDefinition,
    ExplainableMLSettings,
    FeatureRecord,
    PartitionAssignment,
)


def test_mandatory_and_complete_explainable_models(tmp_path: Path) -> None:
    """The mandatory layer should produce genuine held-out SHAP evidence and graphics."""

    inputs = _dataset()
    with pytest.raises(InputValidationError, match="mandatory"):
        run_explainable_models(
            **inputs,
            settings=replace(_settings(), enabled=False),
            random_seed=7,
        )
    result = run_explainable_models(
        **inputs,
        settings=_settings(plot_cache_dir=tmp_path / "plots"),
        random_seed=7,
    )
    assert result.status == "COMPLETE"
    assert result.models[0]["status"] == "COMPLETE"
    assert result.models[0]["validation_roc_auc"] == pytest.approx(1.0)
    assert result.models[0]["analysis_unit"] == (
        "INVERSE_BLOCK_WEIGHTED_PROTEIN_FIT_BLOCK_LEVEL_VALIDATION"
    )
    assert len(result.predictions) == 16
    assert len(result.explanations) == 8
    assert result.feature_importance[0]["importance_rank"] == 1
    assert result.implementation_version
    assert result.shap_version
    assert result.matplotlib_version
    assert {row["plot_type"] for row in result.plots} == {
        "SHAP_BEESWARM",
        "SHAP_GLOBAL_BAR",
        "SHAP_WATERFALL",
    }
    assert all(path.is_file() for _, path in result.plot_assets)


def test_model_statuses_are_explicit_for_weak_datasets(tmp_path: Path) -> None:
    """Small, dependent and constant-feature datasets should not yield fitted models."""

    inputs = _dataset()
    too_small = run_explainable_models(
        **inputs,
        settings=replace(
            _settings(plot_cache_dir=tmp_path / "plots"),
            minimum_samples_per_class=20,
        ),
        random_seed=1,
    )
    assert too_small.models[0]["status"] == "INSUFFICIENT_SAMPLE_SIZE"
    same_groups = tuple(
        replace(item, partition_key=f"class_{'t' if item.protein_id.startswith('t') else 'b'}")
        for item in inputs["partitions"]
    )
    dependent = run_explainable_models(
        **{**inputs, "partitions": same_groups},
        settings=_settings(plot_cache_dir=tmp_path / "plots"),
        random_seed=1,
    )
    assert dependent.models[0]["status"] == "INSUFFICIENT_INDEPENDENT_GROUPS"
    constant = tuple(
        _feature(protein_id=item.protein_id, feature_id="constant") for item in inputs["partitions"]
    )
    no_features = run_explainable_models(
        **{**inputs, "features": constant},
        settings=_settings(plot_cache_dir=tmp_path / "plots"),
        random_seed=1,
    )
    assert no_features.models[0]["status"] == "NO_ELIGIBLE_FEATURES"


def test_ml_excludes_incompletely_assessed_and_all_data_features(tmp_path: Path) -> None:
    """Unknown values and all-data feature definitions must never become model zeros."""

    inputs = _dataset()
    protein_ids = {item.protein_id for item in inputs["partitions"]}
    incomplete = {
        ("TYPE", "signal"): protein_ids - {"b2"},
        ("TYPE", "noise"): protein_ids - {"b2"},
    }
    unavailable = run_explainable_models(
        **inputs,
        settings=_settings(plot_cache_dir=tmp_path / "incomplete"),
        random_seed=4,
        feature_assessment_universes=incomplete,
    )
    assert unavailable.models[0]["status"] == "NO_ELIGIBLE_FEATURES"
    assert "complete assessment" in unavailable.models[0]["status_message"]

    complete = {key: protein_ids for key in incomplete}
    fitted = run_explainable_models(
        **inputs,
        settings=_settings(plot_cache_dir=tmp_path / "complete"),
        random_seed=4,
        feature_assessment_universes=complete,
        exploratory_feature_keys={("TYPE", "noise")},
    )
    assert fitted.models[0]["status"] == "COMPLETE"
    assert {row["feature_id"] for row in fitted.feature_importance} == {"signal"}
    with pytest.raises(InputValidationError, match="do not occur"):
        run_explainable_models(
            **inputs,
            settings=_settings(plot_cache_dir=tmp_path / "invalid"),
            random_seed=4,
            exploratory_feature_keys={("TYPE", "missing")},
        )


def test_underpowered_validation_retains_predictions_and_shap(tmp_path: Path) -> None:
    """Small held-out sets should be described but must not emit inferential metrics."""

    inputs = _dataset()
    partitions = tuple(
        replace(item, partition="DISCOVERY") if item.protein_id in {"tv2", "bv2"} else item
        for item in inputs["partitions"]
    )
    result = run_explainable_models(
        **{**inputs, "partitions": partitions},
        settings=_settings(plot_cache_dir=tmp_path / "underpowered"),
        random_seed=5,
    )
    model = result.models[0]
    assert model["status"] == "COMPLETE_VALIDATION_UNDERPOWERED"
    assert model["validation_roc_auc"] is None
    assert model["validation_average_precision"] is None
    assert "metrics and permutation importance were withheld" in model["status_message"]
    assert len([row for row in result.predictions if row["partition"] == "VALIDATION"]) == 2
    assert {row["partition"] for row in result.explanations} == {"VALIDATION"}
    assert all(
        row["validation_permutation_importance_mean"] is None for row in result.feature_importance
    )

    mixed_partitions = tuple(
        replace(item, partition_key="mixed") if item.protein_id in {"tv2", "bv2"} else item
        for item in inputs["partitions"]
    )
    mixed = run_explainable_models(
        **{**inputs, "partitions": mixed_partitions},
        settings=replace(
            _settings(plot_cache_dir=tmp_path / "mixed"),
            minimum_groups_per_class=1,
        ),
        random_seed=5,
    )
    assert mixed.models[0]["status"] == "COMPLETE_VALIDATION_UNDERPOWERED"
    assert "target proteins=1, background proteins=1" in mixed.models[0]["status_message"]
    assert "mixed blocks=1" in mixed.models[0]["status_message"]


def test_model_rejects_overlap_and_conflicting_feature_names(tmp_path: Path) -> None:
    """Contradictory labels and feature identities should fail before fitting."""

    inputs = _dataset()
    overlapping = (*inputs["label_memberships"], {"protein_id": "t1", "label_id": "bg"})
    with pytest.raises(InputValidationError, match="both classes"):
        run_explainable_models(
            **{**inputs, "label_memberships": overlapping},
            settings=_settings(plot_cache_dir=tmp_path / "plots"),
            random_seed=1,
        )
    conflict = (
        _feature(protein_id="t1", feature_id="same", name="First"),
        _feature(protein_id="t2", feature_id="same", name="Second"),
    )
    with pytest.raises(InputValidationError, match="Conflicting names"):
        run_explainable_models(
            **{**inputs, "features": conflict},
            settings=_settings(plot_cache_dir=tmp_path / "plots"),
            random_seed=1,
        )


def test_feature_selection_matrix_and_model_helpers() -> None:
    """Feature helpers should remain deterministic, bounded and label blind."""

    feature_map = {
        "p1": {("KMER", "a"), ("PFAM", "x")},
        "p2": {("KMER", "a"), ("PFAM", "y")},
        "p3": {("KMER", "b"), ("PFAM", "x")},
    }
    selected = _select_feature_keys(
        protein_ids=("p1", "p2", "p3"),
        features_by_protein=feature_map,
        maximum_features=3,
    )
    assert len(selected) == 3
    assert {key[0] for key in selected[:2]} == {"KMER", "PFAM"}
    matrix = _feature_matrix(
        protein_ids=("p1", "p3"),
        feature_keys=selected,
        features_by_protein=feature_map,
        numpy_module=np,
    )
    assert matrix.shape == (2, 3)
    assert set(np.unique(matrix)) <= {0.0, 1.0}
    assert (
        _select_feature_keys(
            protein_ids=("p1",),
            features_by_protein={"p1": {("TYPE", "constant")}},
            maximum_features=2,
        )
        == ()
    )


def test_block_aggregation_and_inverse_size_weights() -> None:
    """Blocks should be independent metric rows and equally weighted in fitting."""

    weights = _group_balanced_sample_weights(
        labels=np.asarray([0, 0, 0, 1, 1]),
        groups=np.asarray(["b1", "b1", "b2", "t1", "t2"], dtype=object),
        numpy_module=np,
    )
    assert weights.sum() == pytest.approx(5.0)
    assert weights[:2].sum() == pytest.approx(weights[2])
    assert weights[3] == pytest.approx(weights[4])
    matrix, labels, mixed = _aggregate_pure_blocks(
        matrix=np.asarray([[0.0], [1.0], [1.0], [0.0]]),
        labels=np.asarray([0, 0, 1, 0]),
        groups=np.asarray(["pure", "pure", "mixed", "mixed"], dtype=object),
        numpy_module=np,
    )
    assert matrix.shape == (1, 1)
    assert matrix[0, 0] == pytest.approx(0.5)
    assert labels.tolist() == [0]
    assert mixed == 1
    with pytest.raises(InputValidationError, match="pure outcome blocks"):
        _group_balanced_sample_weights(
            labels=np.asarray([0, 1]),
            groups=np.asarray(["mixed", "mixed"], dtype=object),
            numpy_module=np,
        )
    with pytest.raises(InputValidationError, match="matching rows"):
        _aggregate_pure_blocks(
            matrix=np.ones((2, 1)),
            labels=np.asarray([0]),
            groups=np.asarray(["a", "b"], dtype=object),
            numpy_module=np,
        )


def test_cross_validation_fitting_metrics_and_convergence() -> None:
    """Numerical helpers should tune, score and report convergence deterministically."""

    dependencies = _load_dependencies()
    matrix = np.asarray([[0], [0], [1], [1], [0], [0], [1], [1]], dtype=float)
    labels = np.asarray([0, 0, 1, 1, 0, 0, 1, 1], dtype=int)
    groups = np.asarray(["b1", "b2", "t1", "t2", "b3", "b4", "t3", "t4"])
    strength, score, folds = _select_regularisation(
        matrix=matrix,
        labels=labels,
        groups=groups,
        settings=_settings(),
        random_seed=3,
        dependencies=dependencies,
    )
    assert strength in _settings().regularisation_strengths
    assert score == pytest.approx(1.0)
    assert folds == 2
    no_strength = _select_regularisation(
        matrix=matrix[:4],
        labels=labels[:4],
        groups=np.asarray(["b", "b", "t", "t"]),
        settings=_settings(),
        random_seed=3,
        dependencies=dependencies,
    )
    assert no_strength == (None, None, 0)
    model, converged = _fit_model(
        matrix=matrix,
        labels=labels,
        settings=_settings(),
        strength=0.1,
        random_seed=3,
        dependencies=dependencies,
    )
    assert converged is True
    assert model.predict_proba([[1]])[0, 1] > model.predict_proba([[0]])[0, 1]
    weighted_model, _ = _fit_model(
        matrix=matrix,
        labels=labels,
        settings=_settings(),
        strength=0.1,
        random_seed=3,
        dependencies=dependencies,
        sample_weights=_group_balanced_sample_weights(
            labels=labels,
            groups=groups,
            numpy_module=np,
        ),
    )
    assert weighted_model.predict_proba([[1]])[0, 1] > 0.5
    with pytest.raises(InputValidationError, match="sample weights"):
        _fit_model(
            matrix=matrix,
            labels=labels,
            settings=_settings(),
            strength=0.1,
            random_seed=3,
            dependencies=dependencies,
            sample_weights=np.zeros(len(labels)),
        )
    _, short_convergence = _fit_model(
        matrix=matrix,
        labels=labels,
        settings=replace(_settings(), maximum_iterations=1, convergence_tolerance=1e-12),
        strength=0.1,
        random_seed=3,
        dependencies=dependencies,
    )
    assert short_convergence is False
    metrics = _classification_metrics(
        labels=np.asarray([0, 1]),
        probabilities=np.asarray([0.1, 0.9]),
        dependencies=dependencies,
    )
    assert all(value is not None for value in metrics.values())
    empty = _classification_metrics(
        labels=np.asarray([], dtype=int),
        probabilities=np.asarray([], dtype=float),
        dependencies=dependencies,
    )
    assert set(empty.values()) == {None}


def test_permutation_and_linear_shap_helpers() -> None:
    """Held-out importance and additive local log-odds explanations should be auditable."""

    dependencies = _load_dependencies()
    matrix = np.asarray([[0.0, 1.0], [1.0, 0.0], [0.0, 0.0], [1.0, 1.0]])
    labels = np.asarray([0, 1, 0, 1])
    coefficients = np.asarray([4.0, 0.0])
    importance = _held_out_permutation_importance(
        matrix=matrix,
        labels=labels,
        coefficients=coefficients,
        intercept=-2.0,
        repeats=4,
        maximum_features=1,
        random_seed=2,
        dependencies=dependencies,
    )
    assert set(importance) == {0}
    assert importance[0][0] >= 0.0
    assert (
        _held_out_permutation_importance(
            matrix=matrix[:1],
            labels=labels[:1],
            coefficients=coefficients,
            intercept=-2.0,
            repeats=2,
            maximum_features=1,
            random_seed=1,
            dependencies=dependencies,
        )
        == {}
    )
    explanations = _local_explanations(
        comparison_id="cmp",
        partition="VALIDATION",
        protein_ids=("p1",),
        labels=np.asarray([1]),
        probabilities=np.asarray([0.9]),
        matrix=matrix[1:2],
        means=matrix.mean(axis=0),
        coefficients=coefficients,
        shap_values=(matrix[1:2] - matrix.mean(axis=0)) * coefficients,
        feature_keys=(("TYPE", "a"), ("TYPE", "b")),
        feature_names={("TYPE", "a"): "A", ("TYPE", "b"): "B"},
        maximum_per_protein=1,
        numpy_module=np,
    )
    assert explanations[0]["feature_id"] == "a"
    assert explanations[0]["explanation_method"].startswith("SHAP_LINEAR")
    assert explanations[0]["partition"] == "VALIDATION"
    assert (
        _local_explanations(
            comparison_id="cmp",
            partition="VALIDATION",
            protein_ids=(),
            labels=np.asarray([]),
            probabilities=np.asarray([]),
            matrix=np.empty((0, 1)),
            means=np.asarray([0.0]),
            coefficients=np.asarray([1.0]),
            shap_values=np.empty((0, 1)),
            feature_keys=(("TYPE", "a"),),
            feature_names={("TYPE", "a"): "A"},
            maximum_per_protein=1,
            numpy_module=np,
        )
        == ()
    )
    with pytest.raises(InputValidationError, match="does not match"):
        _local_explanations(
            comparison_id="cmp",
            partition="VALIDATION",
            protein_ids=("p1",),
            labels=np.asarray([1]),
            probabilities=np.asarray([0.9]),
            matrix=np.ones((1, 1)),
            means=np.asarray([0.0]),
            coefficients=np.asarray([1.0]),
            shap_values=np.ones((1, 2)),
            feature_keys=(("TYPE", "a"),),
            feature_names={("TYPE", "a"): "A"},
            maximum_per_protein=1,
            numpy_module=np,
        )
    sigmoid = _sigmoid(np.asarray([-1000.0, 0.0, 1000.0]), numpy_module=np)
    assert sigmoid[0] < 1e-300 and sigmoid[1] == 0.5 and sigmoid[2] == 1.0


def test_model_row_and_missing_dependency_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Status templates should count groups and required imports should fail clearly."""

    partitions = {
        "p1": PartitionAssignment("p1", "DISCOVERY", "GROUP", "g1"),
        "p2": PartitionAssignment("p2", "VALIDATION", "GROUP", "g2"),
    }
    row = _model_row(
        comparison_id="cmp",
        discovery=("p1",),
        validation=("p2",),
        samples={"p1": 1, "p2": 0},
        partition_by_protein=partitions,
        settings=_settings(),
        excluded_technical_feature_rows=2,
    )
    assert row["discovery_group_count"] == 1
    assert row["excluded_technical_feature_rows"] == 2
    original_import = builtins.__import__

    def failing_import(name: str, *args: object, **kwargs: object) -> object:
        """Reject only the required numerical dependency import."""

        if name == "numpy":
            raise ImportError("simulated")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", failing_import)
    with pytest.raises(ConfigurationError, match="required"):
        _load_dependencies()


def test_shap_explanation_and_graphics_are_named_and_atomic(tmp_path: Path) -> None:
    """SHAP helpers should publish named global and local PNG/SVG graphics."""

    dependencies = _load_dependencies()
    matrix = np.asarray([[0.0], [0.0], [1.0], [1.0]])
    labels = np.asarray([0, 0, 1, 1])
    model, _ = _fit_model(
        matrix=matrix,
        labels=labels,
        settings=_settings(),
        strength=0.1,
        random_seed=3,
        dependencies=dependencies,
    )
    explanation = _build_shap_explanation(
        model=model,
        background=matrix,
        matrix=matrix,
        feature_keys=(("PFAM_DOMAIN", "PF00646"),),
        feature_names={("PFAM_DOMAIN", "PF00646"): "F-box"},
        dependencies=dependencies,
    )
    assert explanation.values.shape == matrix.shape
    assert "F-box" in explanation.feature_names[0]
    rows, assets = _render_shap_graphics(
        comparison_id="fbox:comparison",
        plot_context={
            "explanation": explanation,
            "protein_ids": ("p1", "p2", "p3", "p4"),
            "partition": "VALIDATION",
            "probabilities": model.predict_proba(matrix)[:, 1],
            "feature_count": 1,
        },
        cache_dir=tmp_path / "plots",
        maximum_waterfalls=1,
        maximum_display=1,
        dependencies=dependencies,
    )
    assert len(rows) == 9
    assert len(assets) == 9
    assert {row["file_format"] for row in rows} == {"PDF", "PNG", "SVG"}
    assert all(path.is_file() and path.stat().st_size > 0 for _, path in assets)
    repeated_rows, repeated_assets = _render_shap_graphics(
        comparison_id="fbox:comparison",
        plot_context={
            "explanation": explanation,
            "protein_ids": ("p1", "p2", "p3", "p4"),
            "partition": "VALIDATION",
            "probabilities": model.predict_proba(matrix)[:, 1],
            "feature_count": 1,
        },
        cache_dir=tmp_path / "plots-repeated",
        maximum_waterfalls=1,
        maximum_display=1,
        dependencies=dependencies,
    )
    assert repeated_rows == rows
    assert [path.read_bytes() for _, path in repeated_assets] == [
        path.read_bytes() for _, path in assets
    ]
    assert _safe_asset_token(value="a/b") == _safe_asset_token(value="a/b")
    assert _safe_asset_token(value="///").startswith("item_")
    with pytest.raises(InputValidationError):
        _safe_asset_token(value="")
    with pytest.raises(InputValidationError, match="non-empty"):
        _build_shap_explanation(
            model=model,
            background=np.empty((0, 1)),
            matrix=matrix,
            feature_keys=(("TYPE", "x"),),
            feature_names={("TYPE", "x"): "X"},
            dependencies=dependencies,
        )
    with pytest.raises(InputValidationError, match="limits"):
        _render_shap_graphics(
            comparison_id="cmp",
            plot_context={},
            cache_dir=tmp_path,
            maximum_waterfalls=-1,
            maximum_display=1,
            dependencies=dependencies,
        )
    with pytest.raises(InputValidationError, match="incomplete"):
        _render_shap_graphics(
            comparison_id="cmp",
            plot_context={"protein_ids": ("p1",), "feature_count": 1},
            cache_dir=tmp_path,
            maximum_waterfalls=1,
            maximum_display=1,
            dependencies=dependencies,
        )
    figure = dependencies["pyplot"].figure()
    with pytest.raises(PublicationError):
        _save_figure_atomic(
            figure=figure,
            path=tmp_path / "bad.eps",
            file_format="eps",
        )
    dependencies["pyplot"].close(figure)

    class BrokenFigure:
        """Figure double that rejects output writes."""

        def savefig(self, *_args: object, **_kwargs: object) -> None:
            """Simulate a filesystem-level graphics failure."""

            raise OSError("simulated")

    with pytest.raises(PublicationError, match="Could not write"):
        _save_figure_atomic(
            figure=BrokenFigure(),
            path=tmp_path / "broken.png",
            file_format="png",
        )


def _dataset() -> dict[str, object]:
    """Return a separable discovery/validation dataset with unique group blocks."""

    discovery_targets = tuple(f"t{index}" for index in range(1, 7))
    discovery_backgrounds = tuple(f"b{index}" for index in range(1, 7))
    validation_targets = ("tv1", "tv2")
    validation_backgrounds = ("bv1", "bv2")
    proteins = (
        *discovery_targets,
        *discovery_backgrounds,
        *validation_targets,
        *validation_backgrounds,
    )
    target_ids = {*discovery_targets, *validation_targets}
    partitions = tuple(
        PartitionAssignment(
            protein_id,
            "VALIDATION" if protein_id.endswith("v1") or protein_id.endswith("v2") else "DISCOVERY",
            "REDUNDANCY_BLOCK",
            f"group:{protein_id}",
        )
        for protein_id in proteins
    )
    memberships = tuple(
        {"protein_id": protein_id, "label_id": "target" if protein_id in target_ids else "bg"}
        for protein_id in proteins
    )
    features = []
    for index, protein_id in enumerate(proteins):
        if protein_id in target_ids:
            features.append(_feature(protein_id=protein_id, feature_id="signal"))
        if index % 2 == 0:
            features.append(_feature(protein_id=protein_id, feature_id="noise"))
        features.append(
            _feature(
                protein_id=protein_id,
                feature_id="available",
                feature_type="STRUCTURE_AVAILABLE",
            )
        )
    return {
        "comparisons": (ComparisonDefinition("cmp", "Comparison", ("target",), ("bg",), ""),),
        "label_memberships": memberships,
        "partitions": partitions,
        "features": tuple(features),
    }


def _settings(*, plot_cache_dir: Path | None = None) -> ExplainableMLSettings:
    """Return fast but fully enabled modelling settings for tests."""

    return ExplainableMLSettings(
        enabled=True,
        minimum_samples_per_class=2,
        minimum_groups_per_class=2,
        maximum_features=10,
        regularisation_strengths=(0.1, 1.0),
        l1_ratio=0.2,
        cross_validation_folds=2,
        maximum_iterations=2000,
        convergence_tolerance=1e-4,
        permutation_repeats=3,
        maximum_permutation_features=2,
        top_local_explanations=2,
        maximum_waterfall_plots=2,
        plot_cache_dir=plot_cache_dir or Path("unused-test-shap-cache"),
        exclude_technical_features=True,
    )


def _feature(
    *, protein_id: str, feature_id: str, name: str = "Feature", feature_type: str = "TYPE"
) -> FeatureRecord:
    """Create one binary feature record."""

    return FeatureRecord(
        protein_id,
        feature_type,
        feature_id,
        name,
        None,
        None,
        "DERIVED",
        "test",
        "test",
    )
