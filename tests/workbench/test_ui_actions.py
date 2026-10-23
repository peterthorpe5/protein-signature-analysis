"""Direct UI action contracts alongside the real headless navigation tests."""

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from protein_signature_app import workbench as ui
from protein_signatures.workbench.jobs import read_job, run_next, submit_job


@pytest.fixture(name="ui_surface")
def fixture_ui_surface(monkeypatch):
    """Provide configurable widgets while retaining the production rendering functions."""
    surface = MagicMock()
    surface.session_state = {}
    surface.choices = {}
    surface.texts = {}
    surface.clicks = set()
    surface.ticks = set()
    surface.button.side_effect = lambda *, label, **_kwargs: label in surface.clicks
    surface.checkbox.side_effect = lambda *, label, **kwargs: (
        label in surface.ticks or kwargs.get("value", False)
    )
    surface.selectbox.side_effect = lambda *, label, options, **_kwargs: surface.choices.get(
        label, list(options)[0]
    )
    surface.radio.side_effect = lambda *, label, options, **_kwargs: surface.choices.get(
        label, list(options)[0]
    )
    surface.text_area.side_effect = lambda *, label, **kwargs: surface.texts.get(
        label, kwargs.get("value", "")
    )
    surface.text_input.side_effect = lambda *, label, **kwargs: surface.texts.get(
        label, kwargs.get("value", "")
    )
    surface.number_input.side_effect = lambda *, label, **kwargs: kwargs.get("value", 1)
    surface.file_uploader.return_value = None
    monkeypatch.setattr(ui, "st", surface)
    return surface


def test_seed_reference_export_and_missing_reference_help(
    workbench_config, ui_surface, monkeypatch
):
    """The user-selected canonical group and bounds reach the actual seed export API."""
    ui_surface.choices["Seed source"] = "Select from the completed result"
    ui.render_seeds(config=workbench_config, database=None)
    assert ui_surface.info.called
    ui_surface.clicks.add("Prepare seed FASTA from this selection")
    ui_surface.texts["Exact selector IDs, one per line"] = "HOG1"
    observed = []
    monkeypatch.setattr(
        ui, "export_seeds", lambda **kwargs: observed.append(kwargs) or b">a\nACD\n"
    )
    ui.render_seeds(config=workbench_config, database=Path("reference.duckdb"))
    assert observed[0]["values"] == ["HOG1"]
    assert ui_surface.session_state["wb-seeds"] == b">a\nACD\n"
    ui_surface.choices["Seed source"] = "Upload or paste protein FASTA"
    ui_surface.clicks = {"Use these seed proteins"}
    with pytest.raises(ValueError):
        ui.render_seeds(config=workbench_config, database=None)
    ui_surface.file_uploader.return_value = SimpleNamespace(getvalue=lambda: b">upload\nACD\n")
    ui.render_seeds(config=workbench_config, database=None)
    assert b"upload" in ui_surface.session_state["wb-seeds"]


def test_representative_search_and_build_buttons_queue_named_operations(
    config_path, workbench_config, ui_surface, monkeypatch
):
    """GUI discovery never substitutes identifier mapping for sequence-based search."""
    captured = []
    monkeypatch.setattr(ui, "queue_job", lambda **kwargs: captured.append(kwargs))
    ui_surface.choices["Discovery route"] = "Search representative sequences"
    ui_surface.clicks.add("Queue representative search")
    with pytest.raises(ValueError, match="seed FASTA"):
        ui.render_discovery(config=workbench_config)
    ui_surface.session_state["wb-seeds"] = b">a\nACD\n"
    ui.render_discovery(config=replace(workbench_config, centroid_database="existing.dmnd"))
    assert captured[-1]["operation"] == "search"
    assert captured[-1]["parameters"]["maximum_hits"] == 50
    ui_surface.clicks = {"Queue representative database build"}
    ui.render_workbench(config_path=config_path)
    assert captured[-1]["operation"] == "build_search"


@pytest.mark.parametrize(
    "task,operation",
    [
        ("Screen members against seeds", "screen"),
        ("Alignment and tree", "align"),
        ("Local Pfam domains", "annotate"),
        ("Frozen signature evidence", "score"),
        ("Prepare expanded campaign inputs", "prepare"),
    ],
)
def test_analysis_submissions_preserve_input_and_scientific_controls(
    workbench_config, ui_surface, monkeypatch, task, operation
):
    """Each displayed task queues its corresponding operation and complete FASTA snapshot."""
    captured = []
    monkeypatch.setattr(ui, "queue_job", lambda **kwargs: captured.append(kwargs))
    ui_surface.choices["Task"] = task
    ui_surface.clicks.add(f"Queue {task.lower()}")
    with pytest.raises(ValueError, match="family FASTA"):
        ui.render_analysis(config=workbench_config)
    ui_surface.file_uploader.return_value = SimpleNamespace(getvalue=lambda: b">a\nACD\n>b\nACE\n")
    if operation in {"screen", "prepare"}:
        with pytest.raises(ValueError, match="original seed"):
            ui.render_analysis(config=workbench_config)
    ui_surface.session_state["wb-seeds"] = b">seed\nACD\n"
    ui.render_analysis(config=workbench_config)
    assert captured[-1]["operation"] == operation
    assert "family" in captured[-1]["uploads"]


def test_completed_family_sources_and_verified_job_download(workbench_config, ui_surface):
    """Selected job FASTA is reused on the server and downloads verify the stored manifest."""
    identifier = submit_job(
        config=workbench_config, operation="retrieve", parameters={"cluster_ids": ["UNIPROT_P1"]}
    )
    run_next(config=workbench_config)
    record = read_job(config=workbench_config, job_id=identifier)
    label = f"retrieve · family.fasta · {record['created_at']} · {identifier}"
    ui_surface.choices["Family FASTA from a completed job"] = label
    uploads, sources = ui.family_input(config=workbench_config)
    assert not uploads and sources["family"]["job_id"] == identifier
    ui_surface.choices["Result file"] = "clusters.tsv"
    ui_surface.clicks.add("Prepare verified download")
    ui_surface.ticks.add("Show table preview (first 1000 rows)")
    ui.render_jobs(config=workbench_config)
    assert ui_surface.dataframe.called and ui_surface.download_button.called
    assert ui_surface.session_state["wb-download"][:2] == (identifier, "clusters.tsv")


def test_alignment_job_view_shows_alignment_conservation_and_phylogram(
    workbench_config, ui_surface, fake_tools
):
    """The result viewer renders the exact completed alignment and its safe tip mapping."""
    identifier = submit_job(
        config=workbench_config,
        operation="align",
        uploads={"family": b">a\nACDEFGHIK\n>b\nACDEFGHIL\n>c\nACDEFGHIM\n"},
        parameters={"tree": True},
    )
    run_next(config=workbench_config)
    assert read_job(config=workbench_config, job_id=identifier)["status"] == "COMPLETE"
    ui_surface.choices["Result file"] = "conservation.tsv"
    ui_surface.ticks.update({"Show interactive phylogram", "Show table preview (first 1000 rows)"})
    ui.render_jobs(config=workbench_config)
    assert ui_surface.plotly_chart.call_count == 2
    assert ui_surface.code.called


def test_queued_failed_and_large_artifact_states(workbench_config, ui_surface, monkeypatch):
    """Cancellation, failure diagnostics and direct large-file transfer remain explicit."""
    identifier = submit_job(
        config=workbench_config, operation="search", uploads={"seeds": b">a\nACD\n"}
    )
    ui_surface.clicks.add("Cancel this waiting job")
    ui.render_jobs(config=workbench_config)
    assert read_job(config=workbench_config, job_id=identifier)["status"] == "CANCELLED"
    failed = submit_job(
        config=workbench_config, operation="search", uploads={"seeds": b">a\nACD\n"}
    )
    run_next(config=workbench_config)
    ui.render_jobs(config=workbench_config)
    assert ui_surface.error.called and ui_surface.code.called
    record = read_job(config=workbench_config, job_id=failed)
    record.update(status="COMPLETE", artifacts={"family.fasta": {"bytes": 200 * 1024 * 1024}})
    monkeypatch.setattr(ui, "list_jobs", lambda **_kwargs: [record])
    ui_surface.choices["Result file"] = "family.fasta"
    ui.render_jobs(config=workbench_config)
    assert "rsync" in ui_surface.code.call_args.kwargs["body"]
    record["artifacts"] = {}
    ui.render_jobs(config=workbench_config)


def test_bad_profile_is_a_contextual_error(config_path, ui_surface):
    """A bad operator file is reported in the page without crashing Streamlit."""
    config_path.write_text(data="invalid TOML", encoding="utf-8")
    ui.render_workbench(config_path=config_path)
    assert ui_surface.error.called
