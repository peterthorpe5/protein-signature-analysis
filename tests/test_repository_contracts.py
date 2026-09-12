"""Tests for shipped configuration, documentation and release metadata."""

from __future__ import annotations

import json
import os
import re
import stat
import subprocess
import textwrap
import tomllib
from pathlib import Path

import jsonschema
import yaml

from protein_signature_app import __version__ as app_version
from protein_signatures import __version__ as package_version
from protein_signatures.profiles import default_profile_comparisons, load_profile


def test_example_campaign_matches_the_published_json_schema() -> None:
    """The generic example should remain valid against the declared public schema."""

    root = Path(__file__).parents[1]
    schema = json.loads((root / "configs/schema.json").read_text(encoding="utf-8"))
    campaign = yaml.safe_load((root / "configs/campaign.example.yaml").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.validate(instance=campaign, schema=schema)


def test_example_profile_and_e3_defaults_are_executable() -> None:
    """The generic profile and the advertised E3 comparison count should stay valid."""

    root = Path(__file__).parents[1]
    profile_schema = json.loads((root / "configs/profile.schema.json").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(profile_schema)
    generic_document = yaml.safe_load(
        (root / "configs/profile.example.yaml").read_text(encoding="utf-8")
    )
    e3_document = yaml.safe_load(
        (root / "src/protein_signatures/data/profiles/e3.yaml").read_text(encoding="utf-8")
    )
    jsonschema.validate(instance=generic_document, schema=profile_schema)
    jsonschema.validate(instance=e3_document, schema=profile_schema)
    generic = load_profile(source=root / "configs/profile.example.yaml")
    assert generic.profile_id == "kinase_example"
    generic_comparisons = default_profile_comparisons(profile=generic)
    assert len(generic_comparisons) == 2
    assert generic.require_structural_evidence is False
    assert {item.background_label_ids for item in generic_comparisons} == {
        ("control:matched_non_kinase",),
        ("control:matched_non_tyrosine_kinase",),
    }
    e3 = load_profile(source="e3")
    comparisons = default_profile_comparisons(profile=e3)
    assert len(e3.labels) == 118
    assert len(comparisons) == 73
    assert any(
        comparison.target_label_ids == ("e3:ubiquitin:crl:crl1_scf:f_box",)
        for comparison in comparisons
    )


def test_documentation_local_links_and_release_metadata_exist() -> None:
    """Every relative Markdown link and release-facing root file should resolve."""

    root = Path(__file__).parents[1]
    markdown_paths = (root / "README.md", *sorted((root / "docs").glob("*.md")))
    pattern = re.compile(r"\[[^]]+\]\((?!https?://|#)([^)]+)\)")
    for markdown_path in markdown_paths:
        text = markdown_path.read_text(encoding="utf-8")
        for match in pattern.finditer(text):
            target = match.group(1).split("#", maxsplit=1)[0]
            assert (markdown_path.parent / target).resolve().exists(), (
                markdown_path,
                target,
            )
    citation = yaml.safe_load((root / "CITATION.cff").read_text(encoding="utf-8"))
    assert citation["version"] == "0.1.0"
    assert (root / "CONTRIBUTING.md").is_file()
    assert (root / "SECURITY.md").is_file()


def test_release_versions_and_launcher_modes_are_consistent() -> None:
    """Release metadata should agree and advertised launchers should be executable."""

    root = Path(__file__).parents[1]
    with (root / "pyproject.toml").open(mode="rb") as handle:
        project = tomllib.load(handle)["project"]
    project_version = project["version"]
    citation = yaml.safe_load((root / "CITATION.cff").read_text(encoding="utf-8"))
    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    alpha_fold = (root / "src/protein_signatures/alphafold.py").read_text(encoding="utf-8")

    assert project_version == package_version == app_version == citation["version"]
    assert re.search(rf"^## {re.escape(package_version)} - ", changelog, flags=re.MULTILINE)
    assert f"protein-signature-analysis/{package_version}" in alpha_fold
    assert project["urls"]["Repository"] == citation["repository-code"]
    assert project["scripts"] == {
        "protein-signatures": "protein_signatures.cli:main",
        "protein-signature-app": "protein_signature_app.launcher:main",
    }

    launchers = (
        "run_tests.sh",
        "run_protein_signature_analysis.sh",
        "run_protein_signature_app.sh",
        "start_from_inputs.sh",
    )
    for launcher in launchers:
        assert (root / launcher).stat().st_mode & stat.S_IXUSR


def test_source_distribution_manifest_excludes_generated_and_sensitive_files() -> None:
    """The source manifest must prune caches and common local secret formats."""

    root = Path(__file__).parents[1]
    manifest = (root / "MANIFEST.in").read_text(encoding="utf-8")
    assert "prune examples/minimal_e3/.protein_signature_cache" in manifest
    assert "global-exclude .env .env.* *.pem *.key *.p12 *.pfx" in manifest


def test_launchers_reject_missing_option_values() -> None:
    """Every value-taking launcher option should fail cleanly when truncated."""

    root = Path(__file__).parents[1]
    cases = (
        ("run_protein_signature_analysis.sh", "--config"),
        ("run_protein_signature_app.sh", "--resource"),
        ("run_protein_signature_app.sh", "--address"),
        ("start_from_inputs.sh", "--work-dir"),
        ("start_from_inputs.sh", "--orthofinder-results"),
    )
    for launcher_name, option in cases:
        result = subprocess.run(
            ("bash", str(root / launcher_name), option),
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 2
        assert f"Option {option} requires a non-empty value." in result.stderr


def test_start_launcher_guards_new_and_existing_config_modes(tmp_path: Path) -> None:
    """New campaigns need inputs and existing campaigns cannot ignore new inputs."""

    root = Path(__file__).parents[1]
    launcher = root / "start_from_inputs.sh"
    missing_new_inputs = subprocess.run(
        ("bash", str(launcher), "--work-dir", str(tmp_path / "new_campaign")),
        capture_output=True,
        text=True,
        check=False,
    )
    assert missing_new_inputs.returncode == 2
    assert "A new campaign requires" in missing_new_inputs.stderr

    work_dir = tmp_path / "campaign"
    work_dir.mkdir()
    (work_dir / "campaign.yaml").write_text("schema_version: 1\n", encoding="utf-8")

    without_resume = subprocess.run(
        ("bash", str(launcher), "--work-dir", str(work_dir)),
        capture_output=True,
        text=True,
        check=False,
    )
    assert without_resume.returncode == 2
    assert "pass --resume" in without_resume.stderr

    with_new_input = subprocess.run(
        (
            "bash",
            str(launcher),
            "--work-dir",
            str(work_dir),
            "--resume",
            "--profile",
            "e3",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert with_new_input.returncode == 2
    assert "sole configuration authority" in with_new_input.stderr
    assert "--profile" in with_new_input.stderr


def test_mutating_launchers_reject_broad_or_option_like_destinations(
    tmp_path: Path,
) -> None:
    """Launchers should reject dangerous roots before invoking Conda or analysis."""

    root = Path(__file__).parents[1]
    config = tmp_path / "campaign.yaml"
    config.write_text("schema_version: 1\n", encoding="utf-8")
    cases = (
        ("start_from_inputs.sh", ("--work-dir", "/"), "unsafe --work-dir"),
        ("start_from_inputs.sh", ("--work-dir", "-campaign"), "must not begin"),
        (
            "run_protein_signature_analysis.sh",
            ("--config", str(config), "--output-dir", "/"),
            "unsafe --output-dir",
        ),
        (
            "run_protein_signature_analysis.sh",
            ("--config", str(config), "--output-dir", "-result"),
            "must not begin",
        ),
    )
    for launcher_name, arguments, message in cases:
        result = subprocess.run(
            ("bash", str(root / launcher_name), *arguments),
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 2
        assert message in result.stderr


def test_start_launcher_initialise_only_updates_existing_environment(tmp_path: Path) -> None:
    """Initialise-only should synchronise Conda, validate YAML and skip analysis."""

    root = Path(__file__).parents[1]
    launcher = root / "start_from_inputs.sh"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_conda = fake_bin / "conda"
    fake_conda.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env bash
            set -Eeuo pipefail
            printf '%s\\n' "$*" >> "${FAKE_CONDA_LOG}"
            if [[ " $* " == *" protein-signatures initialise "* ]]; then
                previous=""
                for argument in "$@"; do
                    if [[ "${previous}" == "--config" ]]; then
                        mkdir -p "$(dirname "${argument}")"
                        printf 'schema_version: 1\\n' > "${argument}"
                        break
                    fi
                    previous="${argument}"
                done
            fi
            exit 0
            """
        ),
        encoding="utf-8",
    )
    fake_conda.chmod(fake_conda.stat().st_mode | stat.S_IXUSR)
    log_path = tmp_path / "conda.log"
    sequences = tmp_path / "proteins.faa"
    assignments = tmp_path / "labels.tsv"
    sequences.write_text(">protein01\nAAAA\n", encoding="utf-8")
    assignments.write_text("header\nrow\n", encoding="utf-8")
    work_dir = tmp_path / "campaign"
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["FAKE_CONDA_LOG"] = str(log_path)
    environment["BASH_COMPAT"] = "3.2"

    result = subprocess.run(
        (
            "bash",
            str(launcher),
            "--work-dir",
            str(work_dir),
            "--campaign-id",
            "campaign_001",
            "--sequences-fasta",
            str(sequences),
            "--label-assignments",
            str(assignments),
            "--initialise-only",
        ),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )

    assert result.returncode == 0, result.stderr
    assert "analysis was not started" in result.stdout
    log = log_path.read_text(encoding="utf-8")
    assert "env update --file" in log
    assert "--prune" in log
    assert "protein-signatures initialise" in log
    assert "protein-signatures validate" in log
    assert "protein-signatures run-all" not in log
