"""Headless interaction coverage for the guided workbench and optional launch mode."""

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from streamlit.testing.v1 import AppTest

from protein_signature_app import launcher
from protein_signature_app.workbench import parse_ids
from protein_signatures.workbench.config import load_config
from protein_signatures.workbench.jobs import list_jobs, run_next


def workbench_test(*, config_path: Path | None) -> AppTest:
    """Build the real workbench rendering script for a temporary operator profile.

    Args:
        config_path: Temporary trusted profile, or None for disabled-help mode.

    Returns:
        Headless Streamlit test application.
    """
    argument = f"Path({str(config_path)!r})" if config_path else "None"
    return AppTest.from_string(
        script=(
            "from pathlib import Path\n"
            "from protein_signature_app.workbench import render_workbench\n"
            f"render_workbench(config_path={argument}, database=None)\n"
        ),
        default_timeout=30,
    ).run()


@pytest.mark.parametrize(
    "stage",
    [
        "Start here",
        "1. Choose seeds",
        "2. Find families",
        "3. Retrieve sequences",
        "4. Analyse family",
        "Jobs & results",
    ],
)
def test_all_guided_stages_render_without_a_reference_result(config_path, stage):
    """Discovery remains usable when the server has no copied completed campaign."""
    app = workbench_test(config_path=config_path)
    app.radio(key="wb-stage").set_value(stage).run()
    assert not app.exception


@pytest.mark.parametrize(
    "task",
    [
        "Screen members against seeds",
        "Alignment and tree",
        "Local Pfam domains",
        "Frozen signature evidence",
        "Prepare expanded campaign inputs",
    ],
)
def test_every_analysis_task_has_usable_help_and_controls(config_path, task):
    """Optional dependencies and required inputs must be explained before submission."""
    app = workbench_test(config_path=config_path)
    app.radio(key="wb-stage").set_value("4. Analyse family").run()
    app.selectbox(key="wb-analysis-task").set_value(task).run()
    assert not app.exception
    assert any(button.label.startswith("Queue ") for button in app.button)


def test_seed_mapping_job_and_completed_retrieval_controls(config_path):
    """A user can start with pasted proteins, queue mapping and inspect completed families."""
    app = workbench_test(config_path=config_path)
    app.radio(key="wb-stage").set_value("1. Choose seeds").run()
    app.text_area(key="wb-seed-paste").set_value(">B\nACDEFGHIK\n").run()
    next(button for button in app.button if button.label == "Use these seed proteins").click().run()
    assert not app.exception and app.session_state["wb-seeds"] == b">B\nACDEFGHIK\n"
    app.radio(key="wb-stage").set_value("2. Find families").run()
    next(
        button for button in app.button if button.label == "Queue identifier mapping"
    ).click().run()
    config = load_config(path=config_path)
    assert list_jobs(config=config)[0]["status"] == "QUEUED"
    run_next(config=config)
    app.radio(key="wb-stage").set_value("3. Retrieve sequences").run()
    app.multiselect(key="wb-clusters").set_value(["UNIPROT_P1"]).run()
    next(
        button for button in app.button if button.label == "Queue complete cluster retrieval"
    ).click().run()
    assert not app.exception
    run_next(config=config)
    assert list_jobs(config=config)[0]["status"] == "COMPLETE"
    app.radio(key="wb-stage").set_value("Jobs & results").run()
    assert not app.exception and len(app.dataframe) >= 1


def test_disabled_profile_help_and_id_parsing():
    """The result-only application still explains how to enable the workbench."""
    app = workbench_test(config_path=None)
    assert not app.exception and app.info
    assert parse_ids(text="\n NR_A.1 \nUNIPROT_P1\n") == ["NR_A.1", "UNIPROT_P1"]


def test_launcher_supports_workbench_only_mode(config_path, monkeypatch):
    """The launcher passes the trusted profile without requiring a result directory."""
    command = launcher.build_streamlit_command(
        database=None, port=8505, address="127.0.0.1", workbench_config=config_path
    )
    assert command[-2:] == ("--workbench-config", str(config_path))
    assert "--resource" not in command
    monkeypatch.setattr(
        launcher.subprocess, "run", lambda *_args, **_kwargs: SimpleNamespace(returncode=0)
    )
    assert launcher.main(argv=["--workbench-config", str(config_path)]) == 0
    assert launcher.main(argv=[]) == 2


def test_main_application_routes_to_workbench_without_result(config_path, monkeypatch):
    """The actual entry file opens the workbench immediately in standalone mode."""
    import protein_signature_app.app as application

    monkeypatch.setattr(sys, "argv", ["app.py", "--workbench-config", str(config_path)])
    app = AppTest.from_file(script_path=str(Path(application.__file__)), default_timeout=30).run()
    assert not app.exception
    assert app.title[0].value == "Protein research workbench"
