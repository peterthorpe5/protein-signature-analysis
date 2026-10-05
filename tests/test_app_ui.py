"""Headless interaction tests for every application page."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from types import TracebackType
from typing import Any

import duckdb
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import protein_signature_app.app as application_module
from protein_signature_app.viewer_help import PAGE_HELP, PAGE_METHODS
from protein_signatures.errors import PublicationError

APP_TEST_TIMEOUT_SECONDS = max(
    60,
    int(os.environ.get("PROTEIN_SIGNATURE_APP_TEST_TIMEOUT_SECONDS", "180")),
)


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

    def expander(self, *_args: object, **_kwargs: object) -> _Block:
        """Accept contextual sidebar help."""

        return _Block()


class _FakeStreamlit:
    """Small Streamlit surface for direct branch coverage."""

    def __init__(self, *, page: str = "Overview", multiselect_empty: bool = False) -> None:
        """Initialise a fake navigation and widget state."""

        self.sidebar = _Sidebar(page=page)
        self.multiselect_empty = multiselect_empty
        self.button_result = False
        self.messages: list[str] = []
        self.session_state: dict[str, object] = {}

    def __getattr__(self, name: str) -> Any:
        """Return a no-op renderer for ordinary output methods."""

        if name in {
            "caption",
            "dataframe",
            "error",
            "info",
            "image",
            "json",
            "link_button",
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

    def selectbox(self, _label: str, options: tuple[str, ...], **_kwargs: object) -> str:
        """Select the first deterministic option."""

        return options[0]

    def radio(self, _label: str, options: tuple[str, ...], **_kwargs: object) -> str:
        """Select the first evidence section."""

        return options[0]

    def checkbox(self, *_args: object, **_kwargs: object) -> bool:
        """Leave optional, expensive evidence panels closed by default."""

        return False

    def text_input(self, *_args: object, **_kwargs: object) -> str:
        """Leave optional exact-ID overrides unset."""

        return ""

    def multiselect(
        self,
        _label: str,
        options: tuple[str, ...],
        *,
        default: tuple[str, ...],
        **_kwargs: object,
    ) -> tuple[str, ...]:
        """Return all defaults or an explicit empty selection."""

        return () if self.multiselect_empty else default

    def button(self, *_args: object, **_kwargs: object) -> bool:
        """Record a button and return its configured activation state."""

        self.messages.append("button")
        return self.button_result

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
            ("Signature records", "195"),
        ]
        for page, title in (
            ("Signature explorer", "Signature explorer"),
            ("Explainable prediction", "Explainable prediction"),
            ("Protein & Pfam", "Protein & Pfam explorer"),
            ("Classes & roles", "Classes & component roles"),
            ("Structures & folds", "Structures & folds"),
            ("Model & alignment explorer", "Model & alignment explorer"),
            ("Orthology & partitions", "Orthology & partitions"),
            ("Canonical data & downloads", "Canonical data & downloads"),
            ("Data quality & provenance", "Data quality & provenance"),
            ("Glossary & help", "Glossary & help"),
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
    assert any("No completed signature" in message for message in fake.messages)


def test_overview_uses_published_counts_and_defers_full_feature_scan(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Opening the landing page must not scan a very large feature table."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(
        application_module,
        "table_count",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("Unexpected count query")),
    )
    queries: list[str] = []

    def query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return only the small signature summary unless coverage is requested."""

        queries.append(sql)
        if "FROM comparisons c LEFT JOIN signatures s" in sql:
            return pd.DataFrame(
                [
                    {
                        "comparison_id": "fbox",
                        "display_name": "F-box versus controls",
                        "target_label_ids": "fbox",
                        "enriched_count": 97,
                        "validated_within_count": 12,
                        "validated_study_count": 7,
                        "complete_count": 97,
                        "insufficient_count": 0,
                        "no_signature_count": 0,
                    },
                    {
                        "comparison_id": "unresolved",
                        "display_name": "Unresolved versus controls",
                        "target_label_ids": "unresolved",
                        "enriched_count": 0,
                        "validated_within_count": 0,
                        "validated_study_count": 0,
                        "complete_count": 0,
                        "insufficient_count": 1,
                        "no_signature_count": 0,
                    },
                ]
            )
        if "FROM signatures" in sql:
            return pd.DataFrame({"feature_type": ["AMINO_ACID_KMER"], "signature_count": [97]})
        return pd.DataFrame(
            {"feature_type": ["AMINO_ACID_KMER"], "feature_count": [120], "protein_count": [24]}
        )

    monkeypatch.setattr(application_module, "query_dataframe", query)
    figures: list[object] = []
    monkeypatch.setattr(
        application_module,
        "_render_plotly_figure",
        lambda *, figure, **_kwargs: figures.append(figure),
    )
    monkeypatch.setattr(application_module, "_render_downloadable_table", lambda **_kwargs: None)
    metadata: dict[str, object] = {
        "counts": {"proteins": 24, "signatures": 97, "structures": 12, "domain_hits": 28},
        "evidence_availability": {},
    }
    application_module._render_overview(
        database=Path("example.duckdb"), metadata=metadata, inventory_identity="snapshot"
    )
    assert len(queries) == 2 and "FROM signatures" in queries[1]
    assert len(figures) == 2
    assert any("insufficient" in message.lower() for message in fake.messages)
    assert (
        application_module._overview_count(
            database=Path("unused"), metadata=metadata, table_name="proteins"
        )
        == 24
    )

    fake.button_result = True
    application_module._render_overview(
        database=Path("example.duckdb"), metadata=metadata, inventory_identity="snapshot"
    )
    assert len(queries) == 5 and "FROM features" in queries[-1]
    assert len(figures) == 5


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

        if "FROM ml_models m" in sql:
            return pd.DataFrame(
                [
                    {
                        "comparison_id": "cmp",
                        "display_name": "Modelled class",
                        "status": "COMPLETE",
                        "validation_target_count": 2,
                    }
                ]
            )
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

        if "FROM ml_models m" in sql:
            return pd.DataFrame(
                [
                    {
                        "comparison_id": "cmp",
                        "display_name": "Underpowered class",
                        "status": "INSUFFICIENT_SAMPLE_SIZE",
                        "validation_target_count": 0,
                    }
                ]
            )
        assert "FROM ml_models" in sql
        return pd.DataFrame([{"status": "INSUFFICIENT_SAMPLE_SIZE"}])

    monkeypatch.setattr(application_module, "query_dataframe", incomplete_model)
    application_module._render_explainable_ml(database=database)
    assert any("No fitted model" in message for message in fake.messages)


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
    assert any(
        "Each point represents one explained protein" in message for message in fake.messages
    )
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
    """TSV is immediate; formatted workbooks and plot files are prepared on demand."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_downloadable_table(
        frame=pd.DataFrame({"protein_id": ["p1"], "score": [0.8]}),
        download_name="Protein table",
        height=250,
    )
    assert fake.messages.count("download_button") == 1
    assert any("protein_id" in message and "score" in message for message in fake.messages)
    fake.button_result = True
    application_module._render_downloadable_table(
        frame=pd.DataFrame({"protein_id": ["p1"], "score": [0.8]}),
        download_name="Protein table",
    )
    assert fake.messages.count("download_button") == 3
    fake.button_result = False

    class Figure:
        """Minimal interactive-figure test double."""

    def unexpected_pdf(**_kwargs: object) -> bytes:
        """Fail if PDF rendering happens before the user requests it."""

        raise AssertionError("PDF rendering must be deferred")

    monkeypatch.setattr(
        application_module,
        "plotly_figure_to_pdf_bytes",
        unexpected_pdf,
    )
    application_module._render_plotly_figure(
        figure=Figure(),
        download_name="Deferred association plot",
    )
    assert fake.messages.count("download_button") == 3
    assert fake.messages.count("button") == 5

    fake.button_result = True
    monkeypatch.setattr(
        application_module,
        "plotly_figure_to_pdf_bytes",
        lambda **_kwargs: b"%PDF-1.7\n%%EOF\n",
    )
    application_module._render_plotly_figure(
        figure=Figure(),
        download_name="Association plot",
    )
    assert fake.messages.count("download_button") == 4
    rendered_plot_count = sum("Figure object" in message for message in fake.messages)
    assert rendered_plot_count == 2

    def fail_pdf(**_kwargs: object) -> bytes:
        """Model an unavailable Plotly PDF runtime."""

        raise PublicationError("PDF runtime unavailable")

    monkeypatch.setattr(application_module, "plotly_figure_to_pdf_bytes", fail_pdf)
    application_module._render_plotly_figure(
        figure=Figure(),
        download_name="Failed plot",
    )
    assert any("PDF runtime unavailable" in message for message in fake.messages)
    assert sum("Figure object" in message for message in fake.messages) == rendered_plot_count + 1


def test_signature_table_formats_screen_without_changing_numeric_export(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Readable scientific notation must not alter raw TSV statistics."""

    fake = _FakeStreamlit()
    visible: list[pd.DataFrame] = []
    downloads: list[bytes] = []
    fake.dataframe = lambda frame, **_kwargs: visible.append(frame)
    fake.download_button = lambda **kwargs: downloads.append(kwargs["data"])
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_downloadable_table(
        frame=pd.DataFrame(
            {
                "feature_id": ["k3:LPD"],
                "discovery_q_value": [0.0000000123],
                "discovery_prevalence_difference": [0.25],
            }
        ),
        download_name="scientific_evidence",
    )
    assert visible[0].iloc[0]["discovery_q_value"] == "1.23e-08"
    assert visible[0].iloc[0]["discovery_prevalence_difference"] == "+0.250"
    assert b"+0.250" not in downloads[0]
    assert b"0.25" in downloads[0]


def test_packaged_mmcif_model_reaches_the_interactive_trace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A verified local mmCIF should render without requiring an AFDB download."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)

    def records(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return the selected local structure and no remote acquisition."""

        if "FROM structures" in sql:
            return pd.DataFrame(
                [
                    {
                        "structure_id": "s1",
                        "structure_source": "packaged",
                        "coordinate_path": "assets/structures/digest.cif",
                        "coordinate_sha256": "digest",
                        "mean_confidence": 90.0,
                        "analysis_eligibility_status": "ELIGIBLE",
                    }
                ]
            )
        return pd.DataFrame(columns=("uniprot_accession",))

    monkeypatch.setattr(application_module, "query_dataframe", records)
    monkeypatch.setattr(application_module, "read_published_model", lambda **_kwargs: b"CIF")
    parsed = []

    class Trace:
        """Record exact-sequence validation on the locally parsed model."""

        def matches_sequence(self, *, sequence: str) -> bool:
            """Accept the fixture's full-length protein sequence."""

            return sequence == "ACG"

    def parse(*, payload: bytes) -> Trace:
        """Record use of the mmCIF parser for packaged coordinates."""

        parsed.append(payload)
        return Trace()

    monkeypatch.setattr(application_module, "parse_mmcif_trace", parse)
    rendered = []
    monkeypatch.setattr(
        application_module, "_render_trace_figure", lambda **kwargs: rendered.append(kwargs)
    )
    application_module._render_model_for_protein(
        database=Path("unused"),
        protein_id="local-p1",
        sequence="ACG",
        scores=[0.1, 0.2, 0.3],
        descriptions=["one", "two", "three"],
        comparison_id="cmp",
    )
    assert parsed == [b"CIF"]
    assert rendered[0]["scores"] == [0.1, 0.2, 0.3]
    assert any("Download selected coordinate model" in message for message in fake.messages)


def test_ranked_associations_shows_fbox_u_box_and_ring_together(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Per-comparison ranking must include every selected enriched class."""

    database = tmp_path / "signatures.duckdb"
    with duckdb.connect(str(database)) as connection:
        connection.execute(
            "CREATE TABLE signatures (comparison_id VARCHAR, feature_type VARCHAR, "
            "feature_id VARCHAR, feature_name VARCHAR, discovery_q_value DOUBLE, "
            "discovery_prevalence_difference DOUBLE, evidence_class VARCHAR, "
            "status VARCHAR)"
        )
        rows = [
            (
                "FBOX",
                "AMINO_ACID_KMER",
                f"f{i}",
                f"F-box {i}",
                0.001,
                0.3,
                "DECISION_CANDIDATE__DISCOVERY_ONLY",
                "COMPLETE",
            )
            for i in range(25)
        ]
        rows += [
            (
                "UBOX",
                "AMINO_ACID_KMER",
                "u1",
                "U-box",
                0.002,
                0.2,
                "DECISION_CANDIDATE__DISCOVERY_ONLY",
                "COMPLETE",
            ),
            (
                "RING",
                "AMINO_ACID_KMER",
                "r1",
                "RING",
                0.003,
                0.1,
                "DECISION_CANDIDATE__DISCOVERY_ONLY",
                "COMPLETE",
            ),
            (
                "RING",
                "AMINO_ACID_KMER",
                "r2",
                "Not significant",
                0.5,
                0.2,
                "DECISION_CANDIDATE__DISCOVERY_ONLY",
                "COMPLETE",
            ),
            (
                "UBOX",
                "AMINO_ACID_KMER",
                "u2",
                "Depleted",
                0.001,
                -0.2,
                "DECISION_CANDIDATE__DISCOVERY_ONLY",
                "COMPLETE",
            ),
        ]
        connection.executemany("INSERT INTO signatures VALUES (?, ?, ?, ?, ?, ?, ?, ?)", rows)

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    tables: list[pd.DataFrame] = []
    figures: list[object] = []
    monkeypatch.setattr(
        application_module,
        "_render_downloadable_table",
        lambda *, frame, **_kwargs: tables.append(frame),
    )
    monkeypatch.setattr(
        application_module,
        "_render_plotly_figure",
        lambda *, figure, **_kwargs: figures.append(figure),
    )
    application_module._render_ranked_associations(
        database=database,
        enriched_comparisons=("FBOX", "UBOX", "RING"),
        descriptions={"FBOX": "F-box", "UBOX": "U-box", "RING": "RING"},
    )
    assert tables[0].groupby("comparison_id").size().to_dict() == {
        "FBOX": 20,
        "UBOX": 1,
        "RING": 1,
    }
    assert tables[0]["within_comparison_rank"].max() == 20
    assert len(figures) == 1

    fake.multiselect_empty = True
    application_module._render_ranked_associations(
        database=database,
        enriched_comparisons=("FBOX", "UBOX", "RING"),
        descriptions={"FBOX": "F-box", "UBOX": "U-box", "RING": "RING"},
    )
    assert len(tables) == 1
    assert any("Select one or more comparisons" in message for message in fake.messages)


def test_comparison_choices_use_display_values_and_disambiguate_names() -> None:
    """Visible widget options must map back to unique published IDs."""

    assert application_module._comparison_choices(
        descriptions={"fbox": "F-box · 97 enriched", "ring": "RING · 30 enriched"}
    ) == {"F-box · 97 enriched": "fbox", "RING · 30 enriched": "ring"}
    assert application_module._comparison_choices(
        descriptions={"first": "Shared label", "second": "Shared label"}
    ) == {
        "Shared label [first]": "first",
        "Shared label [second]": "second",
    }


def test_comparison_outcomes_keep_unanalysed_separate_from_non_significant(
    tmp_path: Path,
) -> None:
    """Discovery, validation and inadequate samples must remain distinct."""

    database = tmp_path / "comparison_outcomes.duckdb"
    with duckdb.connect(str(database)) as connection:
        connection.execute(
            "CREATE TABLE comparisons (comparison_id VARCHAR, display_name VARCHAR, "
            "target_label_ids VARCHAR)"
        )
        connection.execute(
            "CREATE TABLE signatures (comparison_id VARCHAR, status VARCHAR, "
            "discovery_prevalence_difference DOUBLE, discovery_q_value DOUBLE, "
            "evidence_class VARCHAR)"
        )
        connection.executemany(
            "INSERT INTO comparisons VALUES (?, ?, ?)",
            [
                ("signal", "Signal class", "signal"),
                ("weak", "Weak class", "weak"),
                ("small", "Underpowered class", "small"),
            ],
        )
        connection.executemany(
            "INSERT INTO signatures VALUES (?, ?, ?, ?, ?)",
            [
                ("signal", "COMPLETE", 0.3, 0.001, "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE"),
                ("signal", "COMPLETE", 0.1, 0.03, "DECISION_CANDIDATE__DISCOVERY_ONLY"),
                ("signal", "COMPLETE", 0.5, 0.0001, "QC_TECHNICAL_NON_BIOLOGICAL__DISCOVERY_ONLY"),
                ("weak", "COMPLETE", -0.2, 0.04, "DECISION_CANDIDATE__DISCOVERY_ONLY"),
                ("small", "INSUFFICIENT_SAMPLE_SIZE", None, None, ""),
            ],
        )
    outcomes = application_module._comparison_outcomes(database=database).set_index("comparison_id")
    assert outcomes.loc["signal", "enriched_count"] == 2
    assert outcomes.loc["signal", "validated_study_count"] == 1
    assert outcomes.loc["weak", "analysis_status"] == "COMPLETE_NO_POSITIVE_ENRICHMENT"
    assert outcomes.loc["small", "analysis_status"] == "INSUFFICIENT_SAMPLE_SIZE"


def test_matched_control_coverage_counts_targets_not_only_controls(tmp_path: Path) -> None:
    """Multiple controls for one target must not conceal an unmatched target."""

    database = tmp_path / "matched_controls.duckdb"
    with duckdb.connect(str(database)) as connection:
        connection.execute(
            "CREATE TABLE control_matching_audit (background_label_id VARCHAR, "
            "target_unit_id VARCHAR, control_unit_id VARCHAR, status VARCHAR)"
        )
        connection.executemany(
            "INSERT INTO control_matching_audit VALUES (?, ?, ?, ?)",
            [
                ("background", "target_1", "control_a", "MATCHED"),
                ("background", "target_1", "control_b", "MATCHED"),
                ("background", "target_2", "", "UNMATCHED"),
            ],
        )
    coverage = application_module._matched_control_coverage(database=database).iloc[0]
    assert coverage["target_units"] == 2
    assert coverage["covered_target_units"] == 1
    assert coverage["matched_control_units"] == 2
    assert coverage["target_coverage_fraction"] == 0.5


def test_comparison_coverage_uses_published_analysis_cohorts(tmp_path: Path) -> None:
    """Shared background pools must retain distinct comparison denominators."""

    database = tmp_path / "cohorts.duckdb"
    with duckdb.connect(str(database)) as connection:
        connection.execute("CREATE TABLE comparisons (comparison_id VARCHAR, display_name VARCHAR)")
        connection.executemany(
            "INSERT INTO comparisons VALUES (?, ?)",
            [("fbox", "F-box class"), ("other", "Other receptor class")],
        )
    counts = {
        "target_protein_count": 4,
        "control_protein_count": 8,
        "matched_target_unit_count": 2,
        "excluded_unmatched_target_unit_count": 3,
        "excluded_unmatched_target_protein_count": 6,
        "control_unit_count": 4,
    }
    other = {**counts, "matched_target_unit_count": 1, "excluded_unmatched_target_unit_count": 0}
    metadata = {
        "automated_label_evidence": {"matched_comparison_cohorts": {"fbox": counts, "other": other}}
    }
    coverage = application_module._matched_comparison_coverage(
        database=database, metadata=metadata
    ).set_index("comparison_id")
    assert coverage.loc["fbox", "target_units"] == 5
    assert coverage.loc["fbox", "target_coverage_fraction"] == 0.4
    assert coverage.loc["other", "target_units"] == 1
    assert coverage.loc["other", "target_coverage_fraction"] == 1.0
    assert coverage.loc["fbox", "display_name"] == "F-box class"
    assert application_module._matched_comparison_coverage(database=database, metadata={}).empty


def test_feature_meaning_and_page_methods_are_visible() -> None:
    """Opaque feature hashes and sequence words need different explanations."""

    assert "Leu–Pro–Asp" in application_module._feature_explanation(
        feature_type="AMINO_ACID_KMER", feature_id="k3:LPD", feature_name="k3 LPD"
    )
    assert "not a named fold" in application_module._feature_explanation(
        feature_type="STRUCTURE_CLUSTER", feature_id="SC_123", feature_name="Foldseek group"
    )
    assert set(PAGE_METHODS) == set(PAGE_HELP)
    assert all("**Method.**" in note and "**Limit.**" in note for note in PAGE_METHODS.values())


def test_signature_evidence_respects_held_out_and_study_wide_tiers() -> None:
    """A strong discovery cannot silently become a validated candidate."""

    frame = pd.DataFrame(
        [
            {
                "feature_type": "AMINO_ACID_KMER",
                "feature_id": "k3:AAB",
                "status": "COMPLETE",
                "discovery_q_value": 0.001,
                "discovery_prevalence_difference": 0.4,
                "validation_q_value": 0.01,
                "validation_study_q_value": 0.2,
                "validation_prevalence_difference": 0.1,
                "evidence_class": "DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON",
            },
            {
                "feature_type": "STRUCTURE_CLUSTER",
                "feature_id": "sc:strong",
                "status": "COMPLETE",
                "discovery_q_value": 0.02,
                "discovery_prevalence_difference": 0.3,
                "validation_q_value": 0.01,
                "validation_study_q_value": 0.04,
                "validation_prevalence_difference": 0.2,
                "evidence_class": "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE",
            },
            {
                "feature_type": "STRUCTURE_AVAILABLE",
                "feature_id": "technical",
                "status": "COMPLETE",
                "discovery_q_value": 0.0001,
                "discovery_prevalence_difference": 0.5,
                "validation_q_value": 0.001,
                "validation_study_q_value": 0.001,
                "validation_prevalence_difference": 0.5,
                "evidence_class": "QC_TECHNICAL_NON_BIOLOGICAL__VALIDATED_STUDY_WIDE",
            },
        ]
    )
    assert len(application_module._signature_evidence(frame=frame, tier="Discovery")) == 2
    assert (
        len(application_module._signature_evidence(frame=frame, tier="Validated within comparison"))
        == 2
    )
    study = application_module._signature_evidence(frame=frame, tier="Validated study-wide")
    assert study["feature_id"].tolist() == ["sc:strong"]
    with pytest.raises(application_module.InputValidationError, match="Unknown signature"):
        application_module._signature_evidence(frame=frame, tier="unsupported")


def test_alignment_residue_strip_shows_letters_without_injecting_identifiers() -> None:
    """The alignment should be readable while treating supplied IDs as plain text."""

    window = [
        {
            "reference_residue": "A",
            "reference_position": 1,
            "reference_enrichment": 1.0,
            "comparison_residue": "A",
            "comparison_position": 4,
            "comparison_enrichment": 0.0,
            "identity": True,
        },
        {
            "reference_residue": "-",
            "reference_position": None,
            "reference_enrichment": 0.0,
            "comparison_residue": "G",
            "comparison_position": 5,
            "comparison_enrichment": 0.2,
            "identity": False,
        },
    ]
    markup = application_module._alignment_residue_strip(
        rows=window, reference_id="<script>", comparison_id="partner"
    )
    assert "&lt;script&gt;" in markup
    assert "<script>" not in markup
    assert "Identity" in markup and "gap" in markup
    assert "hsl(0, 68%, 78%)" in markup


def test_structure_summary_distinguishes_absent_named_folds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An empty fold vocabulary should not hide published model coverage."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(application_module, "_render_downloadable_table", lambda **_kwargs: None)
    figures = []
    monkeypatch.setattr(
        application_module, "_render_plotly_figure", lambda **kwargs: figures.append(kwargs)
    )

    def records(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return model availability and no explicit named fold assignments."""

        if "count(*) AS models" in sql:
            return pd.DataFrame(
                [
                    {
                        "structure_source": "packaged",
                        "availability_status": "AVAILABLE",
                        "analysis_eligibility_status": "ELIGIBLE",
                        "fold_evidence_status": "NOT_ASSESSED",
                        "models": 12,
                        "mean_confidence": 80.0,
                    }
                ]
            )
        return pd.DataFrame()

    monkeypatch.setattr(application_module, "query_dataframe", records)
    application_module._render_structures(database=Path("unused"))
    assert any("No named fold assignments" in message for message in fake.messages)
    assert len(figures) == 1
    assert figures[0]["figure"].layout.title.text == "Published model coverage and eligibility"


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

        if "FROM comparisons c LEFT JOIN signatures s" in sql:
            return pd.DataFrame(
                [
                    {
                        "comparison_id": "cmp",
                        "display_name": "Sparse comparison",
                        "target_label_ids": "target",
                        "enriched_count": 0,
                        "validated_within_count": 0,
                        "validated_study_count": 0,
                        "complete_count": 1,
                        "insufficient_count": 0,
                        "no_signature_count": 0,
                    }
                ]
            )
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

        if "FROM ml_models m" in sql:
            return pd.DataFrame(
                [
                    {
                        "comparison_id": "cmp",
                        "display_name": "No validation class",
                        "status": "COMPLETE_NO_VALIDATION",
                        "validation_target_count": 0,
                    }
                ]
            )
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
                        "observed_proteins": 0,
                    }
                ]
            )
        return pd.DataFrame(columns=("component_role", "proteins"))

    monkeypatch.setattr(application_module, "query_dataframe", empty_class_query)
    application_module._render_classes(database=database)
    assert any("No positive profile" in message for message in fake.messages)

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
                columns=(
                    "comparison_tool",
                    "minimum_coverage_bin",
                    "tm_score_bin",
                    "comparison_count",
                )
            )
        return pd.DataFrame([{"cluster_id": "c1", "group_support_fraction": 0.8}])

    monkeypatch.setattr(application_module, "query_dataframe", structural_query)
    application_module._render_structures(database=database)
    assert any("Imported exploratory" in message for message in fake.messages)

    def orthology_query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return count frames, a partition row and non-empty group context."""

        if "SELECT count(DISTINCT" in sql:
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
    fake.button_result = True
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
