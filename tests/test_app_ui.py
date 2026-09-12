"""Headless interaction tests for every application page."""

from __future__ import annotations

import sys
from pathlib import Path
from types import TracebackType
from typing import Any

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import protein_signature_app.app as application_module
from protein_signatures.errors import PublicationError

APP_TEST_TIMEOUT_SECONDS = 60


class _Block:
    """Minimal context block used by the direct rendering tests."""

    def __enter__(self) -> _Block:
        """Enter the no-op block."""

        return self

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        """Leave the no-op block without suppressing exceptions."""

        return False

    def metric(self, *_args: object, **_kwargs: object) -> None:
        """Accept a metric rendering call."""


class _Sidebar:
    """Minimal sidebar that returns a configured page."""

    def __init__(self, *, page: str) -> None:
        """Store the navigation choice."""

        self.page = page

    def title(self, *_args: object, **_kwargs: object) -> None:
        """Accept a title call."""

    def caption(self, *_args: object, **_kwargs: object) -> None:
        """Accept a caption call."""

    def radio(self, *_args: object, **_kwargs: object) -> str:
        """Return the configured page."""

        return self.page


class _FakeStreamlit:
    """Small Streamlit surface for direct branch coverage."""

    def __init__(self, *, page: str = "Overview", multiselect_empty: bool = False) -> None:
        """Initialise a fake navigation and widget state."""

        self.sidebar = _Sidebar(page=page)
        self.multiselect_empty = multiselect_empty
        self.messages: list[str] = []

    def __getattr__(self, name: str) -> Any:
        """Return a no-op renderer for ordinary output methods."""

        if name in {
            "caption",
            "dataframe",
            "error",
            "info",
            "image",
            "json",
            "markdown",
            "plotly_chart",
            "set_page_config",
            "subheader",
            "success",
            "title",
            "warning",
            "download_button",
        }:
            return lambda *args, **_kwargs: self.messages.append(str(args[0]) if args else name)
        raise AttributeError(name)

    def columns(self, specification: int | tuple[int, ...]) -> tuple[_Block, ...]:
        """Return the requested number of context blocks."""

        count = specification if isinstance(specification, int) else len(specification)
        return tuple(_Block() for _ in range(count))

    def tabs(self, labels: tuple[str, ...]) -> tuple[_Block, ...]:
        """Return one context block per tab label."""

        return tuple(_Block() for _ in labels)

    def expander(self, *_args: object, **_kwargs: object) -> _Block:
        """Return one context block for an expander."""

        return _Block()

    def selectbox(self, _label: str, options: tuple[str, ...]) -> str:
        """Select the first deterministic option."""

        return options[0]

    def multiselect(
        self, _label: str, options: tuple[str, ...], *, default: tuple[str, ...]
    ) -> tuple[str, ...]:
        """Return all defaults or an explicit empty selection."""

        return () if self.multiselect_empty else default

    def stop(self) -> None:
        """Model Streamlit's terminating stop call."""

        raise RuntimeError("STREAMLIT_STOP")


def test_every_application_page_renders_and_stays_synchronised(
    completed_result: Path,
) -> None:
    """Every navigation choice should render against the same verified resource."""

    application = Path("src/protein_signature_app/app.py").resolve()
    previous = list(sys.argv)
    try:
        sys.argv = [str(application), "--resource", str(completed_result)]
        test_app = AppTest.from_file(
            application,
            default_timeout=APP_TEST_TIMEOUT_SECONDS,
        ).run()
        assert not test_app.exception
        assert [item.value for item in test_app.title if item.value == "Protein signature analysis"]
        assert [(item.label, item.value) for item in test_app.metric][:2] == [
            ("Proteins", "24"),
            ("Candidate signatures", "195"),
        ]
        for page, title in (
            ("Signature explorer", "Signature explorer"),
            ("Explainable prediction", "Explainable prediction"),
            ("Protein & Pfam", "Protein & Pfam explorer"),
            ("Classes & roles", "Classes & component roles"),
            ("Structures & folds", "Structures & folds"),
            ("Orthology & partitions", "Orthology & partitions"),
            ("Canonical data & downloads", "Canonical data & downloads"),
            ("Data quality & provenance", "Data quality & provenance"),
        ):
            test_app.sidebar.radio[0].set_value(page).run()
            assert not test_app.exception
            assert title in [item.value for item in test_app.title]
        test_app.sidebar.radio[0].set_value("Signature explorer").run()
        test_app.multiselect[0].set_value([]).run()
        assert any("Select at least one" in item.value for item in test_app.info)
    finally:
        sys.argv = previous


def test_application_rejects_an_unverified_resource(tmp_path: Path) -> None:
    """The user interface should show a contextual error for an invalid resource."""

    application = Path("src/protein_signature_app/app.py").resolve()
    previous = list(sys.argv)
    try:
        sys.argv = [str(application), "--resource", str(tmp_path / "missing")]
        test_app = AppTest.from_file(
            application,
            default_timeout=APP_TEST_TIMEOUT_SECONDS,
        ).run()
        assert not test_app.exception
        assert any("Could not open the result" in item.value for item in test_app.error)
    finally:
        sys.argv = previous


@pytest.mark.parametrize(
    "page",
    [
        "Overview",
        "Signature explorer",
        "Explainable prediction",
        "Protein & Pfam",
        "Classes & roles",
        "Structures & folds",
        "Orthology & partitions",
        "Canonical data & downloads",
        "Data quality & provenance",
    ],
)
def test_main_routes_directly_to_every_page(
    completed_result: Path,
    monkeypatch: pytest.MonkeyPatch,
    page: str,
) -> None:
    """The Python entry point should route every supported page directly."""

    fake = _FakeStreamlit(page=page)
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(
        sys,
        "argv",
        ["protein-signature-app", "--resource", str(completed_result)],
    )
    application_module.main()
    assert fake.messages


def test_direct_renderers_cover_controlled_empty_states(
    completed_result: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Empty evidence states should remain informative instead of failing."""

    database = completed_result / "protein_signatures.duckdb"
    fake = _FakeStreamlit(multiselect_empty=True)
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_signatures(database=database)
    assert any("Select at least one" in message for message in fake.messages)
    monkeypatch.setattr(application_module, "distinct_values", lambda **_kwargs: ())
    application_module._render_signatures(database=database)
    application_module._render_proteins(database=database)
    empty_features = pd.DataFrame(columns=("feature_type", "feature_count", "protein_count"))
    monkeypatch.setattr(application_module, "query_dataframe", lambda **_kwargs: empty_features)
    application_module._render_overview(database=database, metadata={"evidence_availability": []})
    assert any("No positive feature" in message for message in fake.messages)


def test_explainable_page_renders_complete_and_non_fitted_models(
    completed_result: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The XAI page should expose model metrics, importance, predictions and explanations."""

    database = completed_result / "protein_signatures.duckdb"
    fake = _FakeStreamlit(page="Explainable prediction")
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(application_module, "distinct_values", lambda **_kwargs: ("cmp",))

    def model_frames(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return the canonical frame requested by each XAI query."""

        if "FROM ml_models" in sql:
            return pd.DataFrame(
                [
                    {
                        "status": "COMPLETE",
                        "cv_roc_auc": 0.8,
                        "validation_roc_auc": 0.75,
                        "validation_average_precision": 0.7,
                        "validation_matthews_correlation": 0.5,
                        "shap_explained_partition": "VALIDATION",
                    }
                ]
            )
        if "FROM ml_feature_importance" in sql:
            return pd.DataFrame(
                [
                    {
                        "feature_type": "FOLD",
                        "feature_id": "fold1",
                        "feature_name": "Fold one",
                        "coefficient_log_odds": 2.0,
                        "validation_permutation_importance_mean": 0.2,
                    }
                ]
            )
        if "FROM ml_predictions" in sql:
            return pd.DataFrame(
                [
                    {
                        "protein_id": "p1",
                        "true_class": "TARGET",
                        "predicted_probability": 0.9,
                    },
                    {
                        "protein_id": "p2",
                        "true_class": "BACKGROUND",
                        "predicted_probability": 0.1,
                    },
                ]
            )
        if "FROM ml_plot_inventory" in sql:
            return pd.DataFrame(
                columns=(
                    "comparison_id",
                    "plot_type",
                    "protein_id",
                    "file_format",
                    "asset_path",
                )
            )
        return pd.DataFrame(
            [
                {
                    "protein_id": "p1",
                    "feature_id": "fold1",
                    "log_odds_contribution": 1.2,
                }
            ]
        )

    monkeypatch.setattr(application_module, "query_dataframe", model_frames)
    application_module._render_explainable_ml(database=database)
    assert any("Validation predictions" in message for message in fake.messages)
    assert application_module._metric_text(0.12345) == "0.123"
    assert application_module._metric_text(None) == "NA"
    assert application_module._metric_text(float("nan")) == "NA"

    def incomplete_model(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return an explicit non-fitted model state."""

        assert "FROM ml_models" in sql
        return pd.DataFrame([{"status": "INSUFFICIENT_SAMPLE_SIZE"}])

    monkeypatch.setattr(application_module, "query_dataframe", incomplete_model)
    application_module._render_explainable_ml(database=database)
    assert any("non-fitted" in message for message in fake.messages)


def test_shap_asset_rendering_is_result_scoped(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """SHAP graphics should render and download without accepting escaped paths."""

    result = tmp_path / "result"
    asset_dir = result / "analysis" / "06_explainable_models" / "figures" / "shap" / "cmp"
    asset_dir.mkdir(parents=True)
    database = result / "protein_signatures.duckdb"
    database.touch()
    png = asset_dir / "beeswarm.png"
    svg = asset_dir / "beeswarm.svg"
    pdf = asset_dir / "beeswarm.pdf"
    png.write_bytes(b"\x89PNG\r\n\x1a\nimage")
    svg.write_bytes(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
    pdf.write_bytes(b"%PDF-1.7\n%%EOF\n")
    inventory = pd.DataFrame(
        [
            {
                "comparison_id": "cmp",
                "plot_type": "SHAP_BEESWARM",
                "protein_id": "p1",
                "file_format": "PNG",
                "asset_path": ("analysis/06_explainable_models/figures/shap/cmp/beeswarm.png"),
            },
            {
                "comparison_id": "cmp",
                "plot_type": "SHAP_BEESWARM",
                "protein_id": "p1",
                "file_format": "SVG",
                "asset_path": ("analysis/06_explainable_models/figures/shap/cmp/beeswarm.svg"),
            },
            {
                "comparison_id": "cmp",
                "plot_type": "SHAP_BEESWARM",
                "protein_id": "p1",
                "file_format": "PDF",
                "asset_path": ("analysis/06_explainable_models/figures/shap/cmp/beeswarm.pdf"),
            },
        ]
    )
    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_shap_assets(
        database=database,
        inventory=inventory,
        heading="SHAP figures",
    )
    assert (
        application_module._result_asset_path(
            database=database,
            relative_path="analysis/06_explainable_models/figures/shap/cmp/beeswarm.png",
        )
        == png
    )
    assert any("PNG" in message for message in fake.messages)
    assert fake.messages.count("download_button") == 3
    with pytest.raises(application_module.InputValidationError, match="result-relative"):
        application_module._result_asset_path(
            database=database,
            relative_path=str(png),
        )
    with pytest.raises(application_module.InputValidationError, match="outside"):
        application_module._result_asset_path(
            database=database,
            relative_path="../escape.png",
        )
    application_module._render_shap_assets(
        database=database,
        inventory=pd.DataFrame(columns=inventory.columns),
        heading="Empty figures",
    )
    missing = inventory.copy()
    missing["asset_path"] = "assets/missing.png"
    application_module._render_shap_assets(
        database=database,
        inventory=missing,
        heading="Missing figures",
    )
    assert any("outside" in message for message in fake.messages)


def test_shap_assets_require_safe_content_and_a_pdf_companion(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Static previews must be bounded, inert and accompanied by a PDF export."""

    result = tmp_path / "result"
    asset_dir = result / "assets"
    asset_dir.mkdir(parents=True)
    database = result / "protein_signatures.duckdb"
    database.touch()
    (asset_dir / "plot.png").write_bytes(b"not a PNG")
    (asset_dir / "plot.svg").write_text("<svg onload='alert(1)'></svg>", encoding="utf-8")
    (asset_dir / "valid.png").write_bytes(b"\x89PNG\r\n\x1a\nimage")
    inventory = pd.DataFrame(
        [
            {
                "comparison_id": "cmp",
                "plot_type": "SHAP_GLOBAL_BAR",
                "protein_id": "",
                "file_format": "PNG",
                "asset_path": "assets/plot.png",
            },
            {
                "comparison_id": "cmp",
                "plot_type": "SHAP_GLOBAL_BAR",
                "protein_id": "",
                "file_format": "SVG",
                "asset_path": "assets/plot.svg",
            },
            {
                "comparison_id": "cmp",
                "plot_type": "SHAP_BEESWARM",
                "protein_id": "",
                "file_format": "PNG",
                "asset_path": "assets/valid.png",
            },
        ]
    )
    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_shap_assets(
        database=database,
        inventory=inventory,
        heading="SHAP figures",
    )
    assert any("invalid signature" in message for message in fake.messages)
    assert any("active content" in message for message in fake.messages)
    assert any("PDF companion" in message for message in fake.messages)
    assert not any(message.startswith("b'\\x89PNG") for message in fake.messages)

    monkeypatch.setattr(application_module, "_MAX_APP_ASSET_BYTES", 4)
    with pytest.raises(application_module.InputValidationError, match="allowed range"):
        application_module._read_result_asset(
            database=database,
            relative_path="assets/valid.png",
            file_format="PNG",
        )
    with pytest.raises(application_module.InputValidationError, match="Unsupported"):
        application_module._read_result_asset(
            database=database,
            relative_path="assets/valid.png",
            file_format="HTML",
        )


def test_table_and_plot_renderers_always_offer_declared_downloads(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every rendered table and interactive plot should expose its portable downloads."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_downloadable_table(
        frame=pd.DataFrame({"protein_id": ["p1"], "score": [0.8]}),
        download_name="Protein table",
        height=250,
    )
    assert fake.messages.count("download_button") == 2
    assert any("protein_id" in message and "score" in message for message in fake.messages)

    class Figure:
        """Minimal interactive-figure test double."""

    monkeypatch.setattr(
        application_module,
        "plotly_figure_to_pdf_bytes",
        lambda **_kwargs: b"%PDF-1.7\n%%EOF\n",
    )
    application_module._render_plotly_figure(
        figure=Figure(),
        download_name="Association plot",
    )
    assert fake.messages.count("download_button") == 3
    rendered_plot_count = sum("Figure object" in message for message in fake.messages)
    assert rendered_plot_count == 1

    def fail_pdf(**_kwargs: object) -> bytes:
        """Model an unavailable Plotly PDF runtime."""

        raise PublicationError("PDF runtime unavailable")

    monkeypatch.setattr(application_module, "plotly_figure_to_pdf_bytes", fail_pdf)
    application_module._render_plotly_figure(
        figure=Figure(),
        download_name="Failed plot",
    )
    assert any("PDF runtime unavailable" in message for message in fake.messages)
    assert sum("Figure object" in message for message in fake.messages) == rendered_plot_count


def test_page_renderers_cover_sparse_and_imported_evidence_branches(
    completed_result: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Sparse results and an imported structural summary should remain browsable."""

    database = completed_result / "protein_signatures.duckdb"
    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(application_module, "_render_plotly_figure", lambda **_kwargs: None)

    monkeypatch.setattr(application_module, "distinct_values", lambda **_kwargs: ())
    application_module._render_explainable_ml(database=database)
    assert any("No configured comparison" in message for message in fake.messages)

    monkeypatch.setattr(
        application_module,
        "distinct_values",
        lambda **kwargs: ("cmp",) if kwargs["column_name"] == "comparison_id" else ("TYPE",),
    )

    def sparse_signature_query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return a signature with no plottable effect and an empty ledger."""

        if "FROM signatures" in sql:
            return pd.DataFrame(
                {
                    "comparison_id": ["cmp"],
                    "feature_type": ["TYPE"],
                    "feature_id": ["f1"],
                    "feature_name": ["Feature"],
                    "evidence_class": ["DECISION_CANDIDATE__DISCOVERY_ONLY"],
                    "discovery_prevalence_difference": [None],
                    "discovery_q_value": [None],
                }
            )
        return pd.DataFrame({"partition": [], "feature_id": []})

    monkeypatch.setattr(application_module, "query_dataframe", sparse_signature_query)
    application_module._render_signatures(database=database)

    def sparse_model_query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return a fitted model with no importance or explained samples."""

        if "FROM ml_models" in sql:
            return pd.DataFrame(
                [
                    {
                        "status": "COMPLETE_NO_VALIDATION",
                        "cv_roc_auc": 0.6,
                        "validation_roc_auc": None,
                        "validation_average_precision": None,
                        "validation_matthews_correlation": None,
                        "shap_explained_partition": "VALIDATION",
                    }
                ]
            )
        if "FROM ml_feature_importance" in sql:
            return pd.DataFrame(
                columns=(
                    "feature_type",
                    "feature_id",
                    "feature_name",
                    "coefficient_log_odds",
                    "validation_permutation_importance_mean",
                )
            )
        if "FROM ml_predictions" in sql:
            return pd.DataFrame(columns=("protein_id", "predicted_probability", "true_class"))
        return pd.DataFrame(
            columns=(
                "comparison_id",
                "plot_type",
                "protein_id",
                "file_format",
                "asset_path",
            )
        )

    monkeypatch.setattr(application_module, "query_dataframe", sparse_model_query)
    application_module._render_explainable_ml(database=database)

    def empty_class_query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return an unpopulated controlled vocabulary and empty role table."""

        if "FROM profile_labels" in sql:
            return pd.DataFrame(
                [
                    {
                        "label_id": "protein",
                        "display_name": "Protein",
                        "parent_label_id": "",
                        "level": "root",
                        "system_class": "",
                        "mechanistic_class": "",
                        "component_role": "UNKNOWN",
                        "family": "",
                        "active_site_expected": "UNKNOWN",
                        "active_site_residue": "",
                        "reviewed_proteins": 0,
                    }
                ]
            )
        return pd.DataFrame(columns=("component_role", "proteins"))

    monkeypatch.setattr(application_module, "query_dataframe", empty_class_query)
    application_module._render_classes(database=database)
    assert any("No reviewed-positive" in message for message in fake.messages)

    def structural_query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return non-plottable pairwise rows and one imported summary."""

        if "GROUP BY structure_source" in sql:
            return pd.DataFrame({"structure_source": ["TEST"], "availability_status": ["COMPLETE"]})
        if "WHERE fold_id" in sql:
            return pd.DataFrame(columns=("fold_authority", "fold_id", "fold_name", "proteins"))
        if "FROM structure_clusters" in sql:
            return pd.DataFrame(columns=("cluster_id", "total_members"))
        if "FROM structure_comparisons" in sql:
            return pd.DataFrame(
                {
                    "protein_a_id": ["p1"],
                    "protein_b_id": ["p2"],
                    "comparison_tool": ["test"],
                    "tm_score": [None],
                    "coverage_a": [None],
                    "coverage_b": [None],
                    "rmsd_angstrom": [None],
                }
            )
        return pd.DataFrame([{"cluster_id": "c1", "group_support_fraction": 0.8}])

    monkeypatch.setattr(application_module, "query_dataframe", structural_query)
    application_module._render_structures(database=database)
    assert any("Imported within-group" in message for message in fake.messages)

    def orthology_query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return count frames, a partition row and non-empty group context."""

        if "count(DISTINCT" in sql:
            return pd.DataFrame({"n": [0]})
        if "FROM partitions" in sql:
            return pd.DataFrame(
                {
                    "partition": ["DISCOVERY"],
                    "partition_unit": ["CONNECTED_COMPONENT"],
                    "proteins": [1],
                    "blocks": [1],
                }
            )
        return pd.DataFrame([{"group_id": "HOG1", "member_count": 1}])

    monkeypatch.setattr(application_module, "query_dataframe", orthology_query)
    application_module._render_orthology(database=database)

    monkeypatch.setattr(
        application_module,
        "dataframe_to_xlsx_bytes",
        lambda **_kwargs: (_ for _ in ()).throw(PublicationError("Excel failed")),
    )
    application_module._render_downloadable_table(
        frame=pd.DataFrame({"protein_id": ["p1"]}),
        download_name="failed_table",
    )
    assert any("Excel failed" in message for message in fake.messages)
    assert not any("protein_id" in message and "p1" in message for message in fake.messages[-2:])


def test_main_error_boundary_stops_after_reporting(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The direct application boundary should report and stop on invalid input."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(
        sys,
        "argv",
        ["protein-signature-app", "--resource", str(tmp_path / "missing")],
    )
    with pytest.raises(RuntimeError, match="STREAMLIT_STOP"):
        application_module.main()
    assert any("Could not open" in message for message in fake.messages)
