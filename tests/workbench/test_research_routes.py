"""Exercise the new workbench routes and viewer handoffs through Streamlit."""

import io
import json
import zipfile
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from streamlit.testing.v1 import AppTest

from protein_signature_app import research_tools
from protein_signature_app.research_tools import portable_model_request, trace_figure
from protein_signature_app.research_viewers import (
    EXAMPLE_ALIGNMENT,
    EXAMPLE_TREE,
    alignment_payload,
)
from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.jobs import list_jobs, run_next, submit_job


def app_for(*, path: Path) -> AppTest:
    """Open the actual workbench with a temporary profile and no production data."""
    return AppTest.from_string(
        "from pathlib import Path\n"
        "from protein_signature_app.workbench import render_workbench\n"
        f"render_workbench(config_path=Path({str(path)!r}))\n",
        default_timeout=15,
    ).run()


def test_synthetic_msa_is_visible_and_compact(config_path):
    """The JBrowse route opens on a clearly labelled 6-by-56 synthetic example."""
    data = alignment_payload(fasta=EXAMPLE_ALIGNMENT, newick=EXAMPLE_TREE)
    assert data["height"] < 400 and data["alignment_length"] == 56
    app = app_for(path=config_path)
    app.radio(key="wb-stage").set_value("5. Research tools").run()
    assert not app.exception
    assert any("Synthetic example" in item.value for item in app.info)


def test_seed_to_esmfold_handoff_downloads_exact_input_without_api(config_path):
    """A seed handoff remains usable on a server without a configured API."""
    app = app_for(path=config_path)
    app.session_state["wb-seeds"] = b">one\nAAAA\n>two\nACDE\n"
    app.radio(key="wb-stage").set_value("1. Choose seeds").run()
    next(
        button for button in app.button if button.label == "Open ESMFold for these seeds"
    ).click().run()
    assert app.selectbox(key="wb-research-tool").value == "ESMFold & models"
    assert app.selectbox(key="wb-family-source").value == "Current seed snapshot"
    app.text_area(key="wb-selected-esmfold").set_value("two").run()
    app.button(key="wb-model-prepare").click().run()
    assert not app.exception
    _, fasta, archive = app.session_state["wb-model-request"]
    assert fasta == b">two\nACDE\n"
    with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
        assert bundle.read("model_inputs.fasta") == fasta
        assert json.loads(bundle.read("parameters.json")) == {"remote_allowed": False}
        assert "--operation predict_models" in bundle.read("README.txt").decode()


def test_retrieved_family_portable_request_uses_verified_artifact(workbench_config):
    """Reuse exact proteins from the mounted release without another upload."""
    job = submit_job(
        config=workbench_config, operation="retrieve", parameters={"cluster_ids": ["UNIPROT_P1"]}
    )
    assert run_next(config=workbench_config) == job
    fasta, _ = portable_model_request(
        config=workbench_config,
        selected=["B"],
        uploads={},
        sources={"family": {"job_id": job, "name": "family.fasta"}},
    )
    assert fasta == b">B\nACDEFGHIK\n"


def test_completed_mafft_job_opens_in_jbrowse(config_path, workbench_config, fake_tools):
    """Open the verified alignment and its matching tree from the job directly."""
    job = submit_job(
        config=workbench_config,
        operation="align",
        uploads={"family": b">a\nAAA\n>b\nACA\n"},
    )
    assert run_next(config=workbench_config) == job
    app = app_for(path=config_path)
    app.radio(key="wb-stage").set_value("Jobs & results").run()
    next(button for button in app.button if button.label == "Open in JBrowse MSA").click().run()
    assert app.radio(key="wb-jb-source").value == "Completed MAFFT job"
    assert app.checkbox(key="wb-jb-open").value is True
    assert not app.exception and not app.error


@pytest.mark.parametrize(
    "tool",
    [
        "Predict annotations",
        "Nightingale features",
        "LoL-align structures",
        "Import DeepClust results",
    ],
)
def test_research_routes_render_without_optional_tools(config_path, tool):
    """Optional providers and absent completed jobs must not break navigation."""
    app = app_for(path=config_path)
    app.radio(key="wb-stage").set_value("5. Research tools").run()
    app.selectbox(key="wb-research-tool").set_value(tool).run()
    assert not app.exception and not app.error


def test_trace_uses_transform_and_explicit_confidence_label():
    """Use the reported transform for overlay and only colour confirmed predictions."""
    row = {
        "x": 1,
        "y": 2,
        "z": 3,
        "bfactor": 78,
        "residue": "A",
        "model_residue": "1",
        "chain": "A",
    }
    model = trace_figure(models={"a": [row]}, prediction_confidence=True)
    assert list(model.data[0].marker.color) == [78]
    paired = trace_figure(
        models={"query": [row], "target": [row]},
        transform={"rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "translation": [1, 2, 3]},
    )
    assert list(paired.data[1].x) == [2] and list(paired.data[1].z) == [6]


def upload(*, data: bytes) -> SimpleNamespace:
    """Represent a Streamlit upload using complete bytes and their size."""
    return SimpleNamespace(size=len(data), getvalue=lambda: data)


@pytest.fixture(name="surface")
def fixture_surface(monkeypatch):
    """Replace widgets while exercising the real viewer and import handlers."""
    widgets = MagicMock()
    widgets.selectbox.side_effect = lambda _label, options, **kwargs: options[
        kwargs.get("index", 0)
    ]
    widgets.radio.side_effect = lambda _label, options, **_kwargs: options[-1]
    widgets.checkbox.return_value = True
    widgets.button.return_value = False
    widgets.text_input.return_value = ""
    monkeypatch.setattr(research_tools, "st", widgets)
    monkeypatch.setattr(research_tools, "components", MagicMock())
    return widgets


def test_upload_alignment_and_nightingale_tsv(workbench_config, surface):
    """Consume complete uploaded viewer inputs and provide offline downloads."""
    surface.file_uploader.side_effect = [
        upload(data=b">p1\nAC-D\n>p2\nACED\n"),
        upload(data=b"(p1:0.1,p2:0.2);"),
    ]
    research_tools.render_alignment_viewer(config=workbench_config)
    assert surface.download_button.call_args.kwargs["file_name"] == "jbrowse_protein_viewer.html"
    fasta = b">p1\nACDE\n"
    tsv = (
        "sequence_id\tsequence_sha256\tstart\tend\taccession\tdescription\n"
        f"p1\t{sha256(b'ACDE').hexdigest()}\t2\t4\tPF1\tPortable feature\n"
    ).encode()
    surface.file_uploader.side_effect = [upload(data=fasta), upload(data=tsv)]
    research_tools.render_feature_viewer(config=workbench_config)
    download = surface.download_button.call_args.kwargs
    assert download["file_name"] == "nightingale_protein_viewer.html"
    assert "Portable feature" in download["data"]
    surface.file_uploader.side_effect = [
        upload(data=fasta.replace(b"ACDE", b"ACDF")),
        upload(data=tsv),
    ]
    with pytest.raises(InputValidationError, match="different sequence"):
        research_tools.render_feature_viewer(config=workbench_config)


def test_imported_model_requires_the_original_sequence(workbench_config, surface):
    """Show a portable PDB only when its C-alpha sequence matches the provided FASTA."""
    model = (
        b"ATOM      1  CA  ALA A   1       2.000   1.000   2.000  1.00 85.00           C\n"
        b"ATOM      2  CA  ALA A   2       4.000   1.000   2.000  1.00 86.00           C\nEND\n"
    )
    surface.file_uploader.side_effect = [upload(data=b">p1\nAA\n"), upload(data=model)]
    research_tools.render_imported_model(config=workbench_config)
    assert surface.plotly_chart.call_count == 1
    surface.file_uploader.side_effect = [upload(data=b">p1\nAG\n"), upload(data=model)]
    with pytest.raises(InputValidationError, match="does not match"):
        research_tools.render_imported_model(config=workbench_config)


def test_deepclust_upload_queues_the_selected_columns(workbench_config, surface):
    """Preview the portable file and queue the same selected cluster and fields."""
    target = pa.BufferOutputStream()
    pq.write_table(pa.table({"id": ["a", "b"], "aa": ["AAA", "ACA"]}), target)
    surface.file_uploader.return_value = upload(data=target.getvalue().to_pybytes())
    surface.button.return_value = True
    research_tools.render_deepclust_import(config=workbench_config)
    queued = list_jobs(config=workbench_config)[0]
    assert queued["operation"] == "import_deepclust" and queued["status"] == "QUEUED"
    assert run_next(config=workbench_config) == queued["job_id"]
    assert list_jobs(config=workbench_config)[0]["status"] == "COMPLETE"


def test_structural_report_filter_shows_recorded_method_and_status(surface):
    """Inspect a report by actual method/status and keep the raw upload intact."""
    tsv = (
        b"previous_alignment_tool\tprevious_cluster_id\tlol_status\n"
        b"US-align\tC1\tALIGNED\nTM-align\tC1\tNOT_ASSESSED\n"
        b"TM-align\tC2\tALIGNED\n"
    )
    surface.file_uploader.return_value = upload(data=tsv)
    surface.selectbox.side_effect = ["TM-align", "ALIGNED"]
    research_tools.render_comparison_report()
    table = surface.dataframe.call_args.args[0]
    assert table["previous_cluster_id"].tolist() == ["C2"]
    assert surface.file_uploader.return_value.getvalue() == tsv
