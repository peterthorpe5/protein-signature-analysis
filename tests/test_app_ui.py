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
import streamlit
from streamlit.testing.v1 import AppTest

import protein_signature_app.app as application_module
from protein_signature_app.viewer_help import PAGE_HELP, PAGE_METHODS
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.result_help import column_definition

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
        self.column_config = streamlit.column_config

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
            return lambda *args, **kwargs: self.messages.append(
                str(args[0]) if args else str(kwargs.get("body", name))
            )
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

    def selectbox(self, label: str, options: tuple[str, ...], **_kwargs: object) -> str:
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
            ("Signature outcome rows", "195"),
        ]
        common_help = {
            "❓ What can I answer on this page?",
            "❓ Methods, evidence and limitations",
            "❓ Terms used on this page",
            "❓ Sample-size limits and thresholds for this result",
        }
        assert common_help.issubset({item.label for item in test_app.expander})
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
            assert common_help.issubset({item.label for item in test_app.expander})
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


def test_contextual_help_rejects_unknown_pages_and_explains_status_limits(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Help must use saved model and association limits instead of hard-coded defaults."""
    fake = _FakeStreamlit()
    fake.session_state["help-result-metadata"] = {
        "campaign": {
            "analysis": {"minimum_target_proteins": 7, "minimum_background_proteins": 9},
            "explainable_ml": {"minimum_samples_per_class": 17, "minimum_groups_per_class": 6},
        }
    }
    monkeypatch.setattr(application_module, "st", fake)
    application_module._render_page_help(page="Signature explorer")
    assert any("Association target minimum: 7" in message for message in fake.messages)
    with pytest.raises(InputValidationError, match="No contextual help"):
        application_module._render_page_help(page="Unknown page")
    application_module._render_table_help(
        frame=pd.DataFrame({"status": ["INSUFFICIENT_SAMPLE_SIZE", None]})
    )
    assert any("7 target and 9 background" in message for message in fake.messages)
    application_module._render_table_help(
        frame=pd.DataFrame(
            {"status": ["INSUFFICIENT_SAMPLE_SIZE"], "model_type": ["ELASTIC_NET_LOGISTIC"]}
        )
    )
    assert any("17 protein samples AND 6 pure" in message for message in fake.messages)
    application_module._render_graph_help(graph_name="SHAP_GLOBAL_BAR")
    assert any("not a single discovered motif" in message for message in fake.messages)


def test_table_headers_preserve_small_values_and_define_all_column_types() -> None:
    """Question-mark headings must preserve useful q-values and probability precision."""
    frame = pd.DataFrame(
        {
            "q_value": [1e-12],
            "discovery_p_value": [1e-8],
            "e_value": [1e-20],
            "predicted_probability": [0.9999],
            "protein_id": ["p1"],
            "species_match": [True],
            "target_blocks": [7],
        }
    )
    columns = application_module._table_column_configuration(frame=frame)
    for name in frame.columns:
        assert columns[name]["help"] == f"{name}: {column_definition(column_name=name)}"
        assert columns[name]["label"] == name.replace("_", " ")
    for name in ("q_value", "discovery_p_value", "e_value"):
        assert columns[name]["type_config"]["format"] == "%.3g"
    assert columns["predicted_probability"]["type_config"]["format"] == "%.6f"
    assert columns["species_match"].get("type_config") is None
    assert columns["target_blocks"]["type_config"]["format"] is None


def test_feature_explanation_table_is_detached_and_preserves_producer_text() -> None:
    """Feature descriptions may be added without changing the published rows."""
    source = pd.DataFrame(
        {
            "feature_type": ["AMINO_ACID_KMER", "AMINO_ACID_KMER", "STRUCTURE_CLUSTER", None],
            "feature_id": ["k3:LPD", "k3:LPD", "SC_123", ""],
        }
    )
    augmented = application_module._with_feature_explanations(frame=source)
    assert "feature_explanation" not in source.columns
    assert augmented.loc[0, "feature_explanation"] == augmented.loc[1, "feature_explanation"]
    assert "3-residue" in augmented.loc[0, "feature_explanation"]
    assert "no named fold" in augmented.loc[2, "feature_explanation"]
    assert "No feature identifier" in augmented.loc[3, "feature_explanation"]
    authoritative = source.assign(feature_explanation="Producer definition")
    assert application_module._with_feature_explanations(frame=authoritative).equals(authoritative)
    plain = pd.DataFrame({"protein_id": ["p1"]})
    detached = application_module._with_feature_explanations(frame=plain)
    detached.loc[0, "protein_id"] = "changed"
    assert plain.loc[0, "protein_id"] == "p1"


def test_metric_guide_counts_pure_blocks_and_handles_no_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AP reference counts must exclude mixed blocks and duplicate protein predictions."""
    database = tmp_path / "metrics.duckdb"
    with duckdb.connect(database=str(database)) as connection:
        connection.execute(
            "CREATE TABLE ml_predictions (comparison_id TEXT, partition TEXT, "
            "partition_key TEXT, true_class TEXT)"
        )
        connection.executemany(
            query="INSERT INTO ml_predictions VALUES (?, ?, ?, ?)",
            parameters=[
                ("cmp", "VALIDATION", "t1", "TARGET"),
                ("cmp", "VALIDATION", "t1", "TARGET"),
                ("cmp", "VALIDATION", "b1", "BACKGROUND"),
                ("cmp", "VALIDATION", "b2", "BACKGROUND"),
                ("cmp", "VALIDATION", "b3", "BACKGROUND"),
                ("cmp", "VALIDATION", "mixed", "TARGET"),
                ("cmp", "VALIDATION", "mixed", "BACKGROUND"),
                ("cmp", "DISCOVERY", "d1", "TARGET"),
            ],
        )
    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    record = {"validation_roc_auc": 0.875, "validation_average_precision": 0.9}
    application_module._render_metric_guide(database=database, comparison_id="cmp", record=record)
    assert any("**1 target**, **3 background**" in message for message in fake.messages)
    assert any("0.250 target prevalence" in message for message in fake.messages)
    assert any("87.5%" in message for message in fake.messages)
    fake.messages.clear()
    application_module._render_metric_guide(database=database, comparison_id="absent", record={})
    assert any("Not calculated" in message for message in fake.messages)
    assert not any("Held-out pure blocks:" in message for message in fake.messages)


@pytest.mark.parametrize(
    "selection, table_name",
    [
        ("Decisions", "label_evidence_audit"),
        ("Matched controls", "control_matching_audit"),
        ("Excluded label features", "label_definition_features"),
        ("Abstentions", "unresolved_assignments"),
    ],
)
def test_audit_preview_queries_only_the_selected_bounded_dataset(
    selection: str,
    table_name: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Hidden audits must not be loaded or offered as falsely complete downloads."""
    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(fake, "selectbox", lambda **_kwargs: selection)
    queries = []
    rendered = []

    def query(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Record the exact query while returning a minimal preview."""
        queries.append(sql)
        return pd.DataFrame({"protein_id": ["p1"]})

    monkeypatch.setattr(application_module, "query_dataframe", query)
    monkeypatch.setattr(
        application_module, "_render_downloadable_table", lambda **kwargs: rendered.append(kwargs)
    )
    application_module._render_label_evidence_audit(database=Path("unused"))
    assert queries == rendered == []
    fake.button_result = True
    application_module._render_label_evidence_audit(database=Path("unused"))
    assert len(queries) == 1
    assert queries[0].startswith(f"SELECT * FROM {table_name} ORDER BY ")
    assert queries[0].endswith("LIMIT 5000")
    assert rendered[0]["download_name"] == f"{table_name}_preview"
    assert any("at most 5,000 rows" in message for message in fake.messages)


def test_matching_coverage_counts_distinct_targets_and_handles_empty_pools(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Coverage must measure matched targets, including targets with multiple controls."""
    database = tmp_path / "matching.duckdb"
    with duckdb.connect(database=str(database)) as connection:
        connection.execute(
            "CREATE TABLE control_matching_audit (background_label_id TEXT, "
            "target_unit_id TEXT, control_unit_id TEXT, status TEXT)"
        )
        connection.executemany(
            query="INSERT INTO control_matching_audit VALUES (?, ?, ?, ?)",
            parameters=[
                ("pool", "t1", "c1", "MATCHED"),
                ("pool", "t1", "c2", "MATCHED"),
                ("pool", "t2", "", "UNMATCHED"),
                ("pool", "t3", "", "UNMATCHED"),
            ],
        )
    fake = _FakeStreamlit()
    monkeypatch.setattr(application_module, "st", fake)
    tables = []
    figures = []
    monkeypatch.setattr(
        application_module,
        "_render_downloadable_table",
        lambda **kwargs: tables.append(kwargs["frame"]),
    )
    monkeypatch.setattr(
        application_module,
        "_render_plotly_figure",
        lambda **kwargs: figures.append(kwargs["figure"]),
    )
    application_module._render_control_matching_coverage(database=database)
    record = tables[0].iloc[0]
    assert record["target_units"] == 3
    assert record["covered_target_units"] == 1
    assert record["matched_control_units"] == 2
    assert record["target_coverage_fraction"] == pytest.approx(expected=1 / 3)
    assert len(figures) == 1
    assert any("not test denominators" in message for message in fake.messages)
    with duckdb.connect(database=str(database)) as connection:
        connection.execute("DELETE FROM control_matching_audit")
    application_module._render_control_matching_coverage(database=database)
    assert tables[-1].empty
    assert len(figures) == 1


def test_help_query_failures_preserve_the_backend_error_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unreadable dataset must not be substituted with invented zero coverage."""
    fake = _FakeStreamlit()
    fake.button_result = True
    monkeypatch.setattr(application_module, "st", fake)
    with pytest.raises(InputValidationError, match="DuckDB query failed"):
        application_module._render_label_evidence_audit(database=Path("missing.duckdb"))
    with pytest.raises(InputValidationError, match="DuckDB query failed"):
        application_module._render_control_matching_coverage(database=Path("missing.duckdb"))


@pytest.mark.parametrize(
    "near_count, availability, expected",
    [(0, "INPUT_UNAVAILABLE", "Not assessed"), (0, "COMPLETE", "0"), (2, "COMPLETE", "2")],
)
def test_orthology_distinguishes_unassessed_redundancy_and_plots_blocks(
    near_count: int,
    availability: str,
    expected: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Zero imported rows must not imply an absence assessment when the stage was unavailable."""
    fake = _FakeStreamlit()
    fake.session_state["help-result-metadata"] = {
        "evidence_availability": {"near_redundancy": availability}
    }
    metrics = []
    figures = []

    class MetricBlock(_Block):
        """Capture metric values while preserving rendering context semantics."""

        def metric(self, *args: object, **_kwargs: object) -> None:
            """Capture the displayed label and value."""
            metrics.append(args)

    def records(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Return distinct count authorities and a valid partition allocation."""
        if "FROM redundancy_clusters" in sql:
            return pd.DataFrame({"n": [near_count]})
        if "SELECT count(DISTINCT" in sql:
            return pd.DataFrame({"n": [2]})
        if "GROUP BY ALL" in sql:
            return pd.DataFrame(
                {
                    "partition": ["DISCOVERY", "VALIDATION"],
                    "partition_unit": ["ORTHOFINDER_GROUP", "ORTHOFINDER_GROUP"],
                    "proteins": [20, 5],
                    "blocks": [2, 1],
                }
            )
        return pd.DataFrame()

    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(fake, "columns", lambda _count: tuple(MetricBlock() for _ in range(4)))
    monkeypatch.setattr(application_module, "query_dataframe", records)
    monkeypatch.setattr(application_module, "table_count", lambda **_kwargs: 25)
    monkeypatch.setattr(application_module, "_render_downloadable_table", lambda **_kwargs: None)
    monkeypatch.setattr(
        application_module, "_render_plotly_figure", lambda **kwargs: figures.append(kwargs)
    )
    application_module._render_orthology(database=Path("unused"))
    assert ("Near-redundancy clusters", expected) in metrics
    assert figures[0]["download_name"] == "discovery_validation_block_allocation"
    assert figures[0]["figure"].layout.yaxis.title.text == "Independent partition blocks"
    assert "proteins" in figures[0]["figure"].data[0].hovertemplate
    assert any("No optional OrthoFinder group-context" in message for message in fake.messages)


def test_signature_plot_coordinates_preserve_recorded_values_and_flag_underflow() -> None:
    """Display caps must not change published q-values or admit invalid coordinates."""
    frame = pd.DataFrame(
        data={
            "feature_id": [f"f{index}" for index in range(10)],
            "discovery_q_value": [1.0, 1e-8, 0.0, 1e-320, -0.1, 1.1, None, float("inf"), 0.1, 0.1],
            "discovery_prevalence_difference": [0.2] * 8 + [float("inf"), 1.1],
        }
    )
    original = frame.copy(deep=True)
    plotted = application_module._signature_plot_data(frame=frame)
    assert plotted["feature_id"].tolist() == ["f0", "f1", "f2", "f3"]
    assert plotted["plot_negative_log10_q_value"].tolist() == [0.0, 8.0, 300.0, 300.0]
    assert plotted["plot_q_value_clipped"].tolist() == [False, False, True, True]
    assert plotted["discovery_q_value"].tolist() == [1.0, 1e-8, 0.0, 1e-320]
    pd.testing.assert_frame_equal(left=frame, right=original)


def test_signature_plot_accepts_nullable_fields_and_empty_valid_sets() -> None:
    """Historical unavailable values should give an empty plot rather than an exception."""
    frame = pd.DataFrame(
        data={
            "discovery_q_value": pd.Series(data=[pd.NA, -1.0], dtype="Float64"),
            "discovery_prevalence_difference": pd.Series(data=[pd.NA, 0.1], dtype="Float64"),
        }
    )
    assert application_module._signature_plot_data(frame=frame).empty


@pytest.mark.parametrize("frame", [None, [], pd.DataFrame(data={"unrelated": [1]})])
def test_signature_plot_rejects_missing_coordinate_authorities(frame: object) -> None:
    """A plotting request must not invent q-values or prevalence differences."""
    with pytest.raises(InputValidationError, match="recorded q-values and effects"):
        application_module._signature_plot_data(frame=frame)


@pytest.mark.parametrize("axis", ["−log10 q-value", "Recorded q-value"])
def test_signature_chart_axis_choice_keeps_zero_and_missing_outcomes_explicit(
    axis: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both axes must retain original hover values and disclose omitted status rows."""
    fake = _FakeStreamlit()
    figures = []
    frame = pd.DataFrame(
        data={
            "feature_id": ["k3:AAA", "k3:BBB", "unfitted"],
            "feature_name": ["AAA", "BBB", "Unavailable"],
            "feature_type": ["AMINO_ACID_KMER"] * 3,
            "evidence_class": ["DECISION_CANDIDATE__DISCOVERY_ONLY"] * 3,
            "discovery_q_value": [0.0, 1e-8, None],
            "discovery_prevalence_difference": [0.2] * 3,
        }
    )

    def signature_records(*, sql: str, **_kwargs: object) -> pd.DataFrame:
        """Supply two completed outcomes and one unavailable statistical record."""
        if "AS enriched_count" in sql:
            return pd.DataFrame(
                data={
                    "comparison_id": ["cmp"],
                    "display_name": ["Example target"],
                    "enriched_count": [2],
                    "complete_count": [2],
                    "insufficient_count": [0],
                    "no_signature_count": [0],
                }
            )
        return frame if "SELECT * FROM signatures" in sql else pd.DataFrame()

    def axis_selection(*, label: str, options: tuple[str, ...], **_kwargs: object) -> str:
        """Choose the requested display without changing the comparison."""
        return axis if label == "q-value axis" else options[0]

    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(
        application_module,
        "distinct_values",
        lambda **kwargs: (
            ("cmp",) if kwargs["column_name"] == "comparison_id" else ("AMINO_ACID_KMER",)
        ),
    )
    monkeypatch.setattr(application_module, "query_dataframe", signature_records)
    monkeypatch.setattr(application_module, "_render_downloadable_table", lambda **_kwargs: None)
    monkeypatch.setattr(
        fake,
        "selectbox",
        lambda *args, **kwargs: axis_selection(
            label=kwargs.get("label", args[0] if args else ""),
            options=kwargs.get("options", args[1] if len(args) > 1 else ()),
        ),
    )
    monkeypatch.setattr(
        application_module, "_render_plotly_figure", lambda **kwargs: figures.append(kwargs)
    )
    application_module._render_signatures(database=Path("unused"))
    assert len(figures) == 2
    assert figures[0]["download_name"] == "cmp_Discovery_shortlist_chart"
    assert figures[1]["download_name"] == "signature_scatter"
    figure = figures[1]["figure"]
    logarithmic = axis == "−log10 q-value"
    assert list(figure.data[0].y) == ([300.0, 8.0] if logarithmic else [0.0, 1e-8])
    assert len(figure.layout.shapes) == int(logarithmic)
    if logarithmic:
        assert "q = 0.05" in str(figure.layout.annotations)
    assert any("1 rows without finite" in message for message in fake.messages)
    assert any("display cap of 300" in message for message in fake.messages) == logarithmic


def test_individual_shap_query_uses_full_validation_ranking_and_a_bounded_preview(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The chart must rank individual SHAP values, excluding invalid or other-model rows."""
    database = tmp_path / "shap.duckdb"
    with duckdb.connect(database=str(database)) as connection:
        connection.execute(
            query=(
                "CREATE TABLE ml_feature_importance (comparison_id VARCHAR, feature_type VARCHAR, "
                "feature_id VARCHAR, feature_name VARCHAR, "
                "mean_absolute_validation_contribution DOUBLE, coefficient_log_odds DOUBLE)"
            )
        )
        connection.executemany(
            query="INSERT INTO ml_feature_importance VALUES (?, ?, ?, ?, ?, ?)",
            parameters=[
                ("cmp", "FOLD", f"f{index}", f"Fold {index}", index / 40, 40 - index)
                for index in range(41)
            ]
            + [
                ("cmp", "FOLD", "null", "Unavailable", None, 99),
                ("cmp", "FOLD", "infinite", "Invalid", float("inf"), 99),
                ("cmp", "FOLD", "negative", "Invalid", -1.0, 99),
                ("other", "FOLD", "other", "Other comparison", 999, 99),
            ],
        )
    fake = _FakeStreamlit()
    figures = []
    monkeypatch.setattr(application_module, "st", fake)
    monkeypatch.setattr(
        application_module, "_render_plotly_figure", lambda **kwargs: figures.append(kwargs)
    )
    application_module._render_individual_shap_summary(database=database, comparison_id="cmp")
    assert len(figures) == 1
    figure = figures[0]["figure"]
    assert len(figure.data[0].x) == 30
    assert list(figure.data[0].x) == [index / 40 for index in range(40, 10, -1)]
    assert figure.data[0].y[0] == "Fold 40 · f40"
    assert "log-odds" in figure.layout.xaxis.title.text
    assert any("remainder is excluded" in message for message in fake.messages)
    application_module._render_individual_shap_summary(database=database, comparison_id="missing")
    assert len(figures) == 1
    assert any("No individual validation SHAP summary" in message for message in fake.messages)
