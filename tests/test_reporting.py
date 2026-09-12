"""Unit tests for numbered file-first analysis reports."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.reporting import (
    _association_figure,
    _build_static_figures,
    _draw_no_data,
    _final_shap_relative_path,
    _inventory_row,
    _model_importance_figure,
    _prediction_figure,
    _prevalence_figure,
    _register_figure_set,
    _register_table_exports,
    _safe_report_token,
    _save_figure_atomic,
    _table_frame,
    _validate_fdr_threshold,
    _write_binary_atomic,
    build_human_reports,
)
from protein_signatures.schemas import table_schemas


def test_numbered_reports_publish_tables_figures_and_final_copies(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One populated campaign should produce the complete file-first report contract."""

    tables = {name: () for name in table_schemas()}
    tables.update(
        {
            "comparisons": (
                {
                    "comparison_id": "family/a vs control",
                    "display_name": "Family A",
                },
            ),
            "label_assignments": ({"protein_id": "p1", "label_id": "family:a"},),
            "partitions": (
                {
                    "protein_id": "p1",
                    "partition": "DISCOVERY",
                    "partition_unit": "CONNECTED_COMPONENT",
                    "partition_key": "block1",
                },
            ),
            "orthofinder_group_context": ({"member_count": 4},),
            "features": (
                {
                    "protein_id": "p1",
                    "feature_type": "PFAM_DOMAIN",
                    "feature_id": "PF00646",
                },
            ),
            "domain_assessments": (
                {
                    "protein_id": "p1",
                    "domain_authority": "Pfam",
                    "assessment_status": "ASSESSED_WITH_HIT",
                },
            ),
            "structures": (
                {
                    "protein_id": "p1",
                    "structure_id": "s1",
                    "structure_source": "AlphaFoldDB",
                    "availability_status": "COMPLETE",
                },
            ),
            "associations": (
                {
                    "comparison_id": "family/a vs control",
                    "partition": "DISCOVERY",
                    "feature_type": "PFAM_DOMAIN",
                    "feature_id": "PF00646",
                    "feature_name": "F-box",
                    "target_prevalence": 0.8,
                    "background_prevalence": 0.2,
                    "prevalence_difference": 0.6,
                    "q_value": 0.01,
                    "study_q_value": 0.02,
                    "status": "COMPLETE",
                },
            ),
            "signatures": (
                {
                    "comparison_id": "family/a vs control",
                    "feature_type": "PFAM_DOMAIN",
                    "feature_id": "PF00646",
                    "evidence_class": "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE",
                },
            ),
            "ml_models": (
                {
                    "comparison_id": "family/a vs control",
                    "status": "COMPLETE_VALIDATED",
                },
            ),
            "ml_feature_importance": (
                {
                    "comparison_id": "family/a vs control",
                    "feature_type": "PFAM_DOMAIN",
                    "feature_id": "PF00646",
                    "feature_name": "F-box",
                    "coefficient_log_odds": 1.2,
                },
            ),
            "ml_predictions": (
                {
                    "comparison_id": "family/a vs control",
                    "protein_id": "p1",
                    "partition": "VALIDATION",
                    "true_class": "TARGET",
                    "predicted_probability": 0.9,
                },
            ),
        }
    )
    native_root = Path("analysis") / "06_explainable_models" / "figures" / "shap" / "family_a_token"
    shap_assets = []
    for suffix, payload in (("png", b"png"), ("svg", b"svg"), ("pdf", b"pdf")):
        source = tmp_path / f"shap.{suffix}"
        source.write_bytes(payload)
        shap_assets.append(((native_root / "waterfalls" / f"p1.{suffix}").as_posix(), source))

    observed_association_arguments: list[tuple[str, float]] = []
    original_association_figure = _association_figure

    def record_association_arguments(
        *,
        frame: pd.DataFrame,
        comparison_display_name: str,
        fdr_threshold: float,
    ) -> object:
        observed_association_arguments.append((comparison_display_name, fdr_threshold))
        return original_association_figure(
            frame=frame,
            comparison_display_name=comparison_display_name,
            fdr_threshold=fdr_threshold,
        )

    monkeypatch.setattr(
        "protein_signatures.reporting._association_figure",
        record_association_arguments,
    )

    first = build_human_reports(
        tables=tables,
        cache_dir=tmp_path / "reports",
        fdr_threshold=0.2,
        existing_plot_assets=tuple(shap_assets),
    )
    assets = dict(first.assets)
    assert first.figure_count >= 13
    assert first.workbook_count >= len(table_schemas())
    assert "analysis/01_proteins_and_curation/tables/proteins.tsv" in assets
    assert "analysis/99_final_results/tables/signatures.xlsx" in assets
    assert "analysis/05_association_statistics/figures/00_signature_evidence_classes.pdf" in assets
    final_shap = "analysis/99_final_results/figures/shap/family_a_token/waterfalls/p1.pdf"
    assert assets[final_shap].read_bytes() == b"pdf"
    assert (
        assets["analysis/00_run_information/report_inventory.xlsx"].read_bytes().startswith(b"PK")
    )
    static_pdf = assets[
        "analysis/05_association_statistics/figures/00_signature_evidence_classes.pdf"
    ].read_bytes()
    assert static_pdf.startswith(b"%PDF-") and b"%%EOF" in static_pdf[-1024:]
    assert any(row["section"] == "99_final_results" for row in first.inventory)
    comparison_token = _safe_report_token(value="family/a vs control")
    association_path = (
        f"analysis/05_association_statistics/figures/{comparison_token}_association_landscape.pdf"
    )
    assert association_path in assets
    assert any(
        row["content_id"] == "family/a vs control:association_landscape" for row in first.inventory
    )
    assert observed_association_arguments == [("Family A", 0.2)]


def test_report_validation_and_low_level_writers_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Malformed routing, identifiers, formats and writes should be rejected."""

    with pytest.raises(InputValidationError, match="routing mismatch"):
        build_human_reports(tables={}, cache_dir=tmp_path, fdr_threshold=0.05)
    with pytest.raises(InputValidationError, match="routing mismatch"):
        build_human_reports(
            tables={"unexpected": ()},
            cache_dir=tmp_path,
            fdr_threshold=0.05,
        )
    assert _validate_fdr_threshold(fdr_threshold=0.05) == pytest.approx(0.05)
    for invalid_threshold in (True, None, "invalid", 0.0, 1.01, float("inf")):
        with pytest.raises(InputValidationError, match="FDR threshold"):
            _validate_fdr_threshold(fdr_threshold=invalid_threshold)  # type: ignore[arg-type]
    assert _table_frame(table_name="proteins", records=()).columns.tolist()[0] == "protein_id"
    assert _safe_report_token(value="family/a") == _safe_report_token(value="family/a")
    with pytest.raises(InputValidationError, match="must not be blank"):
        _safe_report_token(value=" ")
    assert _safe_report_token(value="///").startswith("item_")
    valid = "analysis/06_explainable_models/figures/shap/token/waterfalls/p1.pdf"
    assert _final_shap_relative_path(relative_path=valid).endswith(
        "figures/shap/token/waterfalls/p1.pdf"
    )
    with pytest.raises(InputValidationError, match="outside"):
        _final_shap_relative_path(relative_path="assets/ml_shap/p1.pdf")
    with pytest.raises(InputValidationError, match="invalid"):
        _final_shap_relative_path(relative_path="analysis/06_explainable_models/figures/shap")
    with pytest.raises(PublicationError, match="empty"):
        _write_binary_atomic(path=tmp_path / "empty.xlsx", payload=b"")
    _write_binary_atomic(path=tmp_path / "value.bin", payload=b"value")
    assert (tmp_path / "value.bin").read_bytes() == b"value"
    monkeypatch.setattr(
        "protein_signatures.reporting.os.replace",
        lambda *_args: (_ for _ in ()).throw(OSError("rename failure")),
    )
    with pytest.raises(PublicationError, match="Could not write report"):
        _write_binary_atomic(path=tmp_path / "failed.bin", payload=b"value")

    figure, axis = plt.subplots()
    _draw_no_data(axis=axis, message="Nothing available")
    with pytest.raises(PublicationError, match="Unsupported"):
        _save_figure_atomic(figure=figure, path=tmp_path / "x.bad", file_format="bad")

    def fail_savefig(*_args: object, **_kwargs: object) -> None:
        raise OSError("disk failure")

    monkeypatch.setattr(figure, "savefig", fail_savefig)
    with pytest.raises(PublicationError, match="Could not write"):
        _save_figure_atomic(figure=figure, path=tmp_path / "x.pdf", file_format="pdf")
    plt.close(figure)


def test_empty_static_reports_and_duplicate_destinations_fail_closed(tmp_path: Path) -> None:
    """Unavailable evidence should be plotted and repeated report paths rejected."""

    frames = {name: _table_frame(table_name=name, records=()) for name in table_schemas()}
    assets: dict[str, Path] = {}
    inventory: list[dict[str, object]] = []
    assert (
        _build_static_figures(
            frames=frames,
            fdr_threshold=0.05,
            root=tmp_path,
            assets=assets,
            inventory=inventory,
        )
        == 8
    )
    table_assets = {"analysis/01_proteins_and_curation/tables/proteins.tsv": tmp_path / "old"}
    with pytest.raises(PublicationError, match="table report path is duplicated"):
        _register_table_exports(
            frame=frames["proteins"],
            table_name="proteins",
            section="01_proteins_and_curation",
            root=tmp_path,
            assets=table_assets,
            inventory=[],
        )

    figure, _axis = plt.subplots()
    final_assets = {"analysis/99_final_results/figures/repeated.png": tmp_path / "old"}
    with pytest.raises(PublicationError, match="Final figure path is duplicated"):
        _register_figure_set(
            figure=figure,
            stem="repeated",
            section="05_association_statistics",
            description="Description",
            content_id="content",
            root=tmp_path,
            assets=final_assets,
            inventory=[],
            copy_to_final=True,
        )
    plt.close(figure)


def test_individual_report_figures_cover_empty_and_populated_states(tmp_path: Path) -> None:
    """Scientific plotting helpers should render both available and unavailable evidence."""

    association = pd.DataFrame(
        {
            "feature_type": ["PFAM_DOMAIN", "STRUCTURAL_CLUSTER"],
            "feature_id": ["PF00646", "SC1"],
            "feature_name": ["F-box", "Fold one"],
            "prevalence_difference": [0.5, None],
            "target_prevalence": [0.8, None],
            "background_prevalence": [0.3, None],
            "q_value": [0.01, None],
            "study_q_value": [0.02, None],
        }
    )
    importance = pd.DataFrame(
        {
            "feature_type": ["PFAM_DOMAIN", "KMER"],
            "feature_id": ["PF00646", "ABC"],
            "feature_name": ["F-box", "ABC"],
            "coefficient_log_odds": [1.5, -0.5],
        }
    )
    predictions = pd.DataFrame(
        {
            "partition": ["DISCOVERY", "VALIDATION", "OTHER"],
            "true_class": ["TARGET", "BACKGROUND", "UNKNOWN"],
            "predicted_probability": [0.9, 0.2, 0.5],
        }
    )
    association_figure = _association_figure(
        frame=association,
        comparison_display_name="Human-readable comparison",
        fdr_threshold=0.2,
    )
    prevalence_figure = _prevalence_figure(
        frame=association,
        comparison_display_name="Human-readable comparison",
    )
    model_figure = _model_importance_figure(
        frame=importance,
        comparison_display_name="Human-readable comparison",
    )
    prediction_figure = _prediction_figure(
        frame=predictions,
        comparison_display_name="Human-readable comparison",
    )
    figures = (
        association_figure,
        _association_figure(
            frame=association.iloc[1:],
            comparison_display_name="Empty association",
            fdr_threshold=0.2,
        ),
        prevalence_figure,
        _prevalence_figure(
            frame=association.iloc[1:],
            comparison_display_name="Empty prevalence",
        ),
        model_figure,
        prediction_figure,
        _prediction_figure(
            frame=predictions.iloc[2:],
            comparison_display_name="Empty prediction",
        ),
    )
    assert association_figure.axes[0].get_title().endswith("Human-readable comparison")
    threshold_lines = [
        line for line in association_figure.axes[0].lines if line.get_linestyle() == "--"
    ]
    assert len(threshold_lines) == 1
    assert threshold_lines[0].get_ydata() == pytest.approx([-math.log10(0.2)] * 2)
    assert prevalence_figure.axes[0].get_title().endswith("Human-readable comparison")
    assert model_figure.axes[0].get_title().endswith("Human-readable comparison")
    assert prediction_figure.axes[0].get_title().endswith("Human-readable comparison")
    prevalence_figure.canvas.draw()
    legend = prevalence_figure.axes[0].get_legend()
    assert legend is not None
    assert legend.get_window_extent().x0 >= prevalence_figure.axes[0].get_window_extent().x1
    assets: dict[str, Path] = {}
    inventory: list[dict[str, object]] = []
    _register_figure_set(
        figure=figures[0],
        stem="figure",
        section="05_association_statistics",
        description="Description",
        content_id="content",
        root=tmp_path,
        assets=assets,
        inventory=inventory,
        copy_to_final=False,
    )
    assert len(assets) == 3
    with pytest.raises(PublicationError, match="duplicated"):
        _register_figure_set(
            figure=figures[0],
            stem="figure",
            section="05_association_statistics",
            description="Description",
            content_id="content",
            root=tmp_path,
            assets=assets,
            inventory=inventory,
            copy_to_final=False,
        )
    repeated_figure = _association_figure(
        frame=association,
        comparison_display_name="Human-readable comparison",
        fdr_threshold=0.2,
    )
    repeated_assets: dict[str, Path] = {}
    _register_figure_set(
        figure=repeated_figure,
        stem="figure",
        section="05_association_statistics",
        description="Description",
        content_id="content",
        root=tmp_path / "repeated",
        assets=repeated_assets,
        inventory=[],
        copy_to_final=False,
    )
    for file_format in ("png", "svg", "pdf"):
        relative_path = f"analysis/05_association_statistics/figures/figure.{file_format}"
        assert assets[relative_path].read_bytes() == repeated_assets[relative_path].read_bytes()
    plt.close(repeated_figure)
    for figure in figures:
        plt.close(figure)
    assert (
        _inventory_row(
            section="01",
            asset_kind="TABLE",
            content_id="x",
            file_format="TSV",
            relative_path="x.tsv",
            row_count=1,
            description="x",
        )["status"]
        == "COMPLETE"
    )
