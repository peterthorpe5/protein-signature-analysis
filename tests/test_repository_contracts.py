"""Tests for shipped configuration, documentation and release metadata."""

from __future__ import annotations

import json
import os
import re
import shutil
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
        "run_completed_e3_workflow.sh",
        "run_evidence_signature_workflow.sh",
        "submit_protein_signature_workflow_slurm.sh",
        "slurm/protein_signature_workflow_controller.sbatch",
        "slurm/run_completed_e3_workflow.sbatch",
    )
    for launcher in launchers:
        assert (root / launcher).stat().st_mode & stat.S_IXUSR


def test_source_distribution_manifest_excludes_generated_and_sensitive_files() -> None:
    """The source manifest must prune caches and common local secret formats."""

    root = Path(__file__).parents[1]
    manifest = (root / "MANIFEST.in").read_text(encoding="utf-8")
    assert "prune examples/minimal_e3/.protein_signature_cache" in manifest
    assert "global-exclude .env .env.* *.pem *.key *.p12 *.pfx" in manifest
    assert "recursive-include workflow Snakefile E3Snakefile" in manifest
    assert "EvidenceSnakefile" in manifest
    assert "include run_evidence_signature_workflow.sh" in manifest
    assert "recursive-include profiles *.yaml" in manifest
    assert "recursive-include slurm *.sbatch" in manifest


def test_snakemake_workflow_profiles_and_environment_are_consistent() -> None:
    """Packaged workflow rules, executor profile and runtime dependencies should agree."""

    root = Path(__file__).parents[1]
    snakefile = (root / "workflow/Snakefile").read_text(encoding="utf-8")
    e3_snakefile = (root / "workflow/E3Snakefile").read_text(encoding="utf-8")
    evidence_snakefile = (root / "workflow/EvidenceSnakefile").read_text(encoding="utf-8")
    assert "rule validate_campaign:" in snakefile
    assert "rule run_campaign:" in snakefile
    assert "rule verify_campaign:" in snakefile
    assert "protein-signatures workflow-validate" in snakefile
    assert "protein-signatures run-all" in snakefile
    assert "protein-signatures workflow-verify" in snakefile
    assert "mem_mb=ANALYSIS_MEMORY_MB" in snakefile
    assert "ANALYSIS_MARKER" in snakefile
    assert "RESULT_COMPLETION" not in snakefile
    assert "rule ensure_e3_preparation:" in e3_snakefile
    assert "rule stage_e3_label_review:" in e3_snakefile
    assert "rule create_automated_test_labels:" in e3_snakefile
    assert "rule approve_automated_test_labels:" in e3_snakefile
    assert "rule verify_e3_label_review:" in e3_snakefile
    assert "rule initialise_e3_campaign:" in e3_snakefile
    assert "rule run_e3_campaign:" in e3_snakefile
    assert "rule verify_e3_campaign:" in e3_snakefile
    assert "protein-signatures workflow-prepare-e3" in e3_snakefile
    assert "protein-signatures workflow-stage-e3-review" in e3_snakefile
    assert "protein-signatures create-automated-test-labels" in e3_snakefile
    assert "--allow-empty-output-dir" in e3_snakefile
    assert "--automated-test-marker" in e3_snakefile
    assert "protein-signatures workflow-verify-e3-review" in e3_snakefile
    assert "protein-signatures workflow-initialise-e3" in e3_snakefile
    assert "PREPARED_DIR in REVIEWED_LABELS.parents" in e3_snakefile
    assert "--resume" in e3_snakefile
    assert "rule create_evidence_labels:" in evidence_snakefile
    assert "rule initialise_campaign:" in evidence_snakefile
    assert "rule run_campaign:" in evidence_snakefile
    assert "rule verify_campaign:" in evidence_snakefile
    assert "--allow-empty-output-dir" in evidence_snakefile
    assert "--label-definition-features" in evidence_snakefile
    assert "Explicit provisional-evidence acceptance is required" in evidence_snakefile

    local_profile = yaml.safe_load(
        (root / "profiles/local/config.v8+.yaml").read_text(encoding="utf-8")
    )
    slurm_profile = yaml.safe_load(
        (root / "profiles/slurm/config.v8+.yaml").read_text(encoding="utf-8")
    )
    environment = yaml.safe_load((root / "environment.yml").read_text(encoding="utf-8"))
    assert local_profile["rerun-incomplete"] is True
    assert slurm_profile["executor"] == "slurm"
    assert slurm_profile["default-resources"][:2] == [
        "slurm_account=barton",
        "slurm_partition=barton",
    ]
    dependencies = environment["dependencies"]
    assert "gawk" in dependencies
    assert "snakemake>=9,<10" in dependencies
    assert "snakemake-executor-plugin-slurm>=2.7.1,<3" in dependencies
    assert environment["channels"][-1] == "nodefaults"


def test_completed_e3_automated_smoke_route_is_explicit_and_submit_only(
    tmp_path: Path,
) -> None:
    """One Slurm command should request all classes without creating local state."""

    root = Path(__file__).parents[1]
    launcher = root / "run_completed_e3_workflow.sh"
    predecessor = tmp_path / "completed_e3_run"
    predecessor.mkdir()
    work = tmp_path / "all_classes_smoke_test"
    result = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "all",
            "--run-root",
            str(predecessor),
            "--work-dir",
            str(work),
            "--campaign-id",
            "all_classes_smoke_test",
            "--automated-test-labels",
            "--submit-slurm",
            "--slurm-account",
            "barton",
            "--slurm-partition",
            "barton",
            "--slurm-dry-run",
            "--threads",
            "24",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    command = result.stdout.replace("\\ ", " ")
    assert "nothing was submitted" in command
    assert "--automated-test-labels" in command
    assert "--test-target-label ALL" in command
    assert "--test-samples-per-class 20" in command
    assert "--account=barton" in command
    assert "--partition=barton" in command
    assert not work.exists()

    unsafe = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "all",
            "--run-root",
            str(predecessor),
            "--work-dir",
            str(tmp_path / "scientific_campaign"),
            "--campaign-id",
            "scientific_campaign",
            "--automated-test-labels",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert unsafe.returncode == 2
    assert "require 'smoke' or 'test'" in unsafe.stderr


def test_launchers_reject_missing_option_values() -> None:
    """Every value-taking launcher option should fail cleanly when truncated."""

    root = Path(__file__).parents[1]
    cases = (
        ("run_protein_signature_analysis.sh", "--config"),
        ("run_protein_signature_app.sh", "--resource"),
        ("run_protein_signature_app.sh", "--address"),
        ("start_from_inputs.sh", "--work-dir"),
        ("start_from_inputs.sh", "--orthofinder-results"),
        ("start_from_inputs.sh", "--foldseek-maximum-hits"),
        ("run_completed_e3_workflow.sh", "--run-root"),
        ("run_completed_e3_workflow.sh", "--minimum-mean-plddt"),
        ("run_completed_e3_workflow.sh", "--slurm-memory"),
        ("run_completed_e3_workflow.sh", "--slurm-scratch-base"),
        ("run_completed_e3_workflow.sh", "--slurm-min-scratch-free-gib"),
        ("run_protein_signature_analysis.sh", "--memory-mb"),
        ("submit_protein_signature_workflow_slurm.sh", "--controller-memory"),
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


def test_completed_e3_launcher_enforces_phases_and_review_gate(tmp_path: Path) -> None:
    """The predecessor launcher should reject unsafe or unreviewed transitions."""

    root = Path(__file__).parents[1]
    launcher = root / "run_completed_e3_workflow.sh"
    missing_phase = subprocess.run(
        ("bash", str(launcher), "--work-dir", str(tmp_path / "work")),
        capture_output=True,
        text=True,
        check=False,
    )
    assert missing_phase.returncode == 2
    assert "--phase and --work-dir are required" in missing_phase.stderr
    invalid_phase = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "guess",
            "--work-dir",
            str(tmp_path / "work"),
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid_phase.returncode == 2
    assert "prepare, approve, initialise, all, run or verify" in invalid_phase.stderr

    invalid_confidence = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--work-dir",
            str(tmp_path / "work"),
            "--minimum-mean-plddt",
            "101",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid_confidence.returncode == 2
    assert "number from 0 to 100" in invalid_confidence.stderr

    run_root = tmp_path / "run"
    run_root.mkdir()
    work = tmp_path / "work"
    prepared = work / "prepared_inputs"
    prepared.mkdir(parents=True)
    (prepared / "PREPARED.json").write_text("{}\n", encoding="utf-8")
    labels = prepared / "label_assignments.REVIEW_REQUIRED.tsv"
    labels.write_text("header\nrow\n", encoding="utf-8")
    rejected = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "initialise",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work),
            "--campaign-id",
            "campaign",
            "--label-assignments",
            str(labels),
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert rejected.returncode == 2
    assert "Checksum-bound label approval is missing" in rejected.stderr


def test_completed_e3_launcher_submits_bounded_slurm_worker(tmp_path: Path) -> None:
    """Prepare submission should preserve arguments without recursive submission."""

    root = Path(__file__).parents[1]
    launcher = root / "run_completed_e3_workflow.sh"
    run_root = tmp_path / "completed_run"
    run_root.mkdir()
    work_dir = tmp_path / "campaign"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    capture = tmp_path / "sbatch_arguments.txt"
    environment_capture = tmp_path / "sbatch_environment.txt"
    scratch_environment_capture = tmp_path / "sbatch_scratch_environment.txt"
    sbatch = fake_bin / "sbatch"
    sbatch.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$@" > "${FAKE_SBATCH_LOG}"\n'
        'printf \'%s\\n\' "${SLURM_CPUS_PER_TASK-unset}" > "${FAKE_SBATCH_ENV}"\n'
        "printf '%s\\n' \"${PROTEIN_SIGNATURE_SCRATCH_BASE-unset}\" "
        '"${PROTEIN_SIGNATURE_WORK_DIR-unset}" '
        '"${PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB-unset}" '
        '> "${FAKE_SBATCH_SCRATCH_ENV}"\n'
        "printf '98765;cluster\\n'\n",
        encoding="utf-8",
    )
    sbatch.chmod(sbatch.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["FAKE_SBATCH_LOG"] = str(capture)
    environment["FAKE_SBATCH_ENV"] = str(environment_capture)
    environment["FAKE_SBATCH_SCRATCH_ENV"] = str(scratch_environment_capture)
    environment["SLURM_CPUS_PER_TASK"] = "99"
    environment["BASH_COMPAT"] = "3.2"

    result = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--submit-slurm",
            "--slurm-account",
            "barton",
            "--slurm-partition",
            "barton",
            "--slurm-memory",
            "64G",
            "--slurm-time",
            "04:00:00",
            "--slurm-scratch-base",
            str(tmp_path / "preferred_scratch"),
            "--slurm-min-scratch-free-gib",
            "7",
            "--threads",
            "4",
        ),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )

    assert result.returncode == 0, result.stderr
    assert "Submitted Slurm job 98765" in result.stdout
    arguments = capture.read_text(encoding="utf-8").splitlines()
    assert "--mem=64G" in arguments
    assert "--time=04:00:00" in arguments
    assert "--cpus-per-task=4" in arguments
    assert "--account=barton" in arguments
    assert "--partition=barton" in arguments
    assert str(root / "slurm/run_completed_e3_workflow.sbatch") in arguments
    assert "--submit-slurm" not in arguments
    assert environment_capture.read_text(encoding="utf-8").strip() == "unset"
    assert scratch_environment_capture.read_text(encoding="utf-8").splitlines() == [
        str(tmp_path / "preferred_scratch"),
        str(work_dir),
        "7",
    ]
    assert (work_dir / "slurm_logs").is_dir()


def test_completed_e3_prepare_is_owned_by_snakemake_and_stages_review(
    tmp_path: Path,
) -> None:
    """The adapter launcher should express preparation and review as one DAG target."""

    root = Path(__file__).parents[1]
    launcher = root / "run_completed_e3_workflow.sh"
    run_root = tmp_path / "completed_run"
    run_root.mkdir()
    work_dir = tmp_path / "campaign"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    conda_log = tmp_path / "conda.log"
    conda = fake_bin / "conda"
    conda.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$*" >> "${FAKE_CONDA_LOG}"\nexit 0\n',
        encoding="utf-8",
    )
    conda.chmod(conda.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["FAKE_CONDA_LOG"] = str(conda_log)
    environment["BASH_COMPAT"] = "3.2"

    result = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--campaign-id",
            "fixture_campaign",
            "--threads",
            "6",
        ),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert result.returncode == 0, result.stderr
    commands = conda_log.read_text(encoding="utf-8")
    assert f"--snakefile {root / 'workflow/E3Snakefile'}" in commands
    assert "review_ready --config" in commands
    assert f"e3_run_root={run_root}" in commands
    assert f"e3_work_dir={work_dir}" in commands
    assert f"e3_reviewed_labels={work_dir / 'reviewed_label_assignments.tsv'}" in commands
    assert "e3_threads=6" in commands
    assert "Curate that file" in result.stdout


def test_completed_e3_launcher_validates_slurm_options_and_dry_run(tmp_path: Path) -> None:
    """Malformed scheduler requests should fail before submission and dry-run should be inert."""

    launcher = Path(__file__).parents[1] / "run_completed_e3_workflow.sh"
    run_root = tmp_path / "completed_run"
    run_root.mkdir()
    work_dir = tmp_path / "campaign"
    without_submit = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--work-dir",
            str(work_dir),
            "--slurm-memory",
            "64G",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert without_submit.returncode == 2
    assert "require --submit-slurm" in without_submit.stderr
    invalid_memory = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--submit-slurm",
            "--slurm-memory",
            "lots",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid_memory.returncode == 2
    assert "positive Slurm size" in invalid_memory.stderr
    invalid_time = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--submit-slurm",
            "--slurm-time",
            "04:99:00",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid_time.returncode == 2
    assert "HH:MM:SS" in invalid_time.stderr
    invalid_scratch = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--submit-slurm",
            "--slurm-scratch-base",
            "relative/scratch",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid_scratch.returncode == 2
    assert "--slurm-scratch-base must be absolute" in invalid_scratch.stderr
    invalid_scratch_space = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--submit-slurm",
            "--slurm-min-scratch-free-gib",
            "0",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid_scratch_space.returncode == 2
    assert "must be a positive integer" in invalid_scratch_space.stderr
    dry_run = subprocess.run(
        (
            "bash",
            str(launcher),
            "--phase",
            "prepare",
            "--run-root",
            str(run_root),
            "--work-dir",
            str(work_dir),
            "--submit-slurm",
            "--slurm-dry-run",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert dry_run.returncode == 0, dry_run.stderr
    assert "nothing was submitted" in dry_run.stdout
    assert "--mem=128G" in dry_run.stdout
    assert not work_dir.exists()


def test_completed_e3_slurm_worker_checks_allocation_and_executes(tmp_path: Path) -> None:
    """The direct worker should validate CPUs, scratch and exact launcher arguments."""

    root = Path(__file__).parents[1]
    worker = root / "slurm/run_completed_e3_workflow.sbatch"
    runner_log = tmp_path / "runner.log"
    temporary_environment_log = tmp_path / "temporary_environment.log"
    runner = tmp_path / "runner.sh"
    runner.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$@" > "${FAKE_RUNNER_LOG}"\n'
        'printf \'%s\\n\' "${TMPDIR}" "${TMP}" "${TEMP}" '
        '> "${FAKE_TEMPORARY_ENVIRONMENT_LOG}"\n',
        encoding="utf-8",
    )
    runner.chmod(runner.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment.update(
        {
            "SLURM_JOB_ID": "1122",
            "SLURM_CPUS_PER_TASK": "3",
            "PROTEIN_SIGNATURE_REQUESTED_CPUS": "4",
            "PROTEIN_SIGNATURE_SCRATCH_BASE": str(tmp_path / "scratch"),
            "PROTEIN_SIGNATURE_WORK_DIR": str(tmp_path / "work"),
            "PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB": "1",
            "FAKE_RUNNER_LOG": str(runner_log),
            "FAKE_TEMPORARY_ENVIRONMENT_LOG": str(temporary_environment_log),
        }
    )
    mismatch = subprocess.run(
        ("bash", str(worker), str(runner), "--phase", "prepare"),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert mismatch.returncode == 2
    assert "fewer than requested" in mismatch.stderr
    environment["SLURM_CPUS_PER_TASK"] = "5"
    extra_allocation = subprocess.run(
        ("bash", str(worker), str(runner), "--phase", "prepare"),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert extra_allocation.returncode == 0, extra_allocation.stderr
    environment["SLURM_CPUS_PER_TASK"] = "4"
    success = subprocess.run(
        ("bash", str(worker), str(runner), "--phase", "prepare"),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert success.returncode == 0, success.stderr
    assert "Job ID: 1122" in success.stdout
    assert "Scratch source: PROTEIN_SIGNATURE_SCRATCH_BASE" in success.stdout
    assert runner_log.read_text(encoding="utf-8").splitlines() == ["--phase", "prepare"]
    temporary_paths = temporary_environment_log.read_text(encoding="utf-8").splitlines()
    assert len(set(temporary_paths)) == 1
    assert temporary_paths[0].endswith("/scratch/protein_signature_1122/generic_tmp")
    assert not (tmp_path / "scratch/protein_signature_1122").exists()
    provenance = tmp_path / "work/slurm_logs/protein_signature_scratch_1122.tsv"
    assert "source\tPROTEIN_SIGNATURE_SCRATCH_BASE" in provenance.read_text(encoding="utf-8")


def test_completed_e3_slurm_worker_falls_back_and_retains_failed_scratch(
    tmp_path: Path,
) -> None:
    """Invalid scheduler scratch should fall back, and failed-job evidence should remain."""

    root = Path(__file__).parents[1]
    worker = root / "slurm/run_completed_e3_workflow.sbatch"
    unavailable_slurm_tmp = tmp_path / "not_a_directory"
    unavailable_slurm_tmp.write_text("occupied\n", encoding="utf-8")
    fallback_tmp = tmp_path / "fallback_tmp"
    temporary_environment_log = tmp_path / "failed_temporary_environment.log"
    runner = tmp_path / "failing_runner.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        'printf \'%s\\n\' "${TMPDIR}" > "${FAKE_TEMPORARY_ENVIRONMENT_LOG}"\n'
        "exit 7\n",
        encoding="utf-8",
    )
    runner.chmod(runner.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment.pop("PROTEIN_SIGNATURE_SCRATCH_BASE", None)
    environment.update(
        {
            "SLURM_JOB_ID": "2233",
            "SLURM_CPUS_PER_TASK": "2",
            "PROTEIN_SIGNATURE_REQUESTED_CPUS": "2",
            "PROTEIN_SIGNATURE_WORK_DIR": str(tmp_path / "work"),
            "PROTEIN_SIGNATURE_MIN_SCRATCH_FREE_GIB": "1",
            "SLURM_TMPDIR": str(unavailable_slurm_tmp),
            "TMPDIR": str(fallback_tmp),
            "FAKE_TEMPORARY_ENVIRONMENT_LOG": str(temporary_environment_log),
        }
    )

    result = subprocess.run(
        ("bash", str(worker), str(runner)),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )

    expected_scratch = fallback_tmp / "protein_signature_2233"
    assert result.returncode == 7
    assert "Skipping unavailable SLURM_TMPDIR" in result.stderr
    assert "Job scratch retained after exit 7" in result.stderr
    assert "Scratch source: TMPDIR" in result.stdout
    assert temporary_environment_log.read_text(encoding="utf-8").strip() == str(
        expected_scratch / "generic_tmp"
    )
    assert expected_scratch.is_dir()


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


def test_generic_analysis_launcher_routes_through_snakemake(tmp_path: Path) -> None:
    """The normal input-driven runner should use the packaged workflow and local profile."""

    root = Path(__file__).parents[1]
    launcher = root / "run_protein_signature_analysis.sh"
    config = tmp_path / "campaign.yaml"
    config.write_text("schema_version: 1\n", encoding="utf-8")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    log = tmp_path / "conda.log"
    conda = fake_bin / "conda"
    conda.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$*" >> "${FAKE_CONDA_LOG}"\nexit 0\n',
        encoding="utf-8",
    )
    conda.chmod(conda.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["FAKE_CONDA_LOG"] = str(log)
    environment["BASH_COMPAT"] = "3.2"
    output = tmp_path / "campaign" / "result"
    state = tmp_path / "campaign" / "workflow_state"

    result = subprocess.run(
        (
            "bash",
            str(launcher),
            "--config",
            str(config),
            "--output-dir",
            str(output),
            "--workflow-state-dir",
            str(state),
            "--threads",
            "7",
            "--memory-mb",
            "70000",
            "--runtime-minutes",
            "360",
        ),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )

    assert result.returncode == 0, result.stderr
    commands = log.read_text(encoding="utf-8")
    assert "env update --file" in commands
    assert "snakemake --version" in commands
    assert f"--snakefile {root / 'workflow/Snakefile'}" in commands
    assert f"--profile {root / 'profiles/local'}" in commands
    assert f"workflow_state_dir={state}" in commands
    assert f"workflow_result_dir={output}" in commands
    assert "workflow_threads=7" in commands
    assert "workflow_memory_mb=70000" in commands
    assert "workflow_runtime_minutes=360" in commands
    assert "--forcerun validate_campaign run_campaign verify_campaign" in commands


def test_generic_controller_dry_run_and_worker_contract(tmp_path: Path) -> None:
    """Controller submission should be inspectable and its worker should enforce Slurm."""

    root = Path(__file__).parents[1]
    submitter = root / "submit_protein_signature_workflow_slurm.sh"
    config = tmp_path / "campaign.yaml"
    config.write_text("schema_version: 1\n", encoding="utf-8")
    work_dir = tmp_path / "campaign"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    conda = fake_bin / "conda"
    conda.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    conda.chmod(conda.stat().st_mode | stat.S_IXUSR)
    flock = fake_bin / "flock"
    flock.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    flock.chmod(flock.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["BASH_COMPAT"] = "3.2"

    result = subprocess.run(
        (
            "bash",
            str(submitter),
            "--config",
            str(config),
            "--work-dir",
            str(work_dir),
            "--threads",
            "24",
            "--memory-mb",
            "128000",
            "--dry-run",
        ),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert result.returncode == 0, result.stderr
    assert "nothing was submitted" in result.stdout
    assert "--mem=4G" in result.stdout
    assert "--profile slurm" in result.stdout.replace("\\ ", " ")
    assert "--memory-mb 128000" in result.stdout.replace("\\ ", " ")
    assert not work_dir.exists()

    worker = root / "slurm/protein_signature_workflow_controller.sbatch"
    true_executable = shutil.which("true")
    assert true_executable is not None
    outside_slurm = subprocess.run(
        ("bash", str(worker), "--state-dir", str(work_dir), "--", true_executable),
        capture_output=True,
        text=True,
        check=False,
    )
    assert outside_slurm.returncode == 2
    assert "must be launched by sbatch" in outside_slurm.stderr
    # The worker runs on Linux Slurm, where flock is a required dependency.
    # Supply the fixture implementation so this contract remains portable to macOS.
    worker_environment = environment.copy()
    worker_environment["SLURM_JOB_ID"] = "24680"
    inside_slurm = subprocess.run(
        ("bash", str(worker), "--state-dir", str(work_dir), "--", true_executable),
        capture_output=True,
        text=True,
        check=False,
        env=worker_environment,
    )
    assert inside_slurm.returncode == 0, inside_slurm.stderr
    assert "Controller job ID: 24680" in inside_slurm.stdout


def test_generic_controller_submitter_syncs_environment_and_submits(tmp_path: Path) -> None:
    """Real controller mode should sync dependencies and preserve one exact worker argv."""

    root = Path(__file__).parents[1]
    submitter = root / "submit_protein_signature_workflow_slurm.sh"
    config = tmp_path / "campaign.yaml"
    config.write_text("schema_version: 1\n", encoding="utf-8")
    work_dir = tmp_path / "campaign"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    conda_log = tmp_path / "conda.log"
    sbatch_log = tmp_path / "sbatch.log"
    conda = fake_bin / "conda"
    conda.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$*" >> "${FAKE_CONDA_LOG}"\nexit 0\n',
        encoding="utf-8",
    )
    conda.chmod(conda.stat().st_mode | stat.S_IXUSR)
    sbatch = fake_bin / "sbatch"
    sbatch.write_text(
        "#!/usr/bin/env bash\nprintf '%s\\n' \"$@\" > \"${FAKE_SBATCH_LOG}\"\nprintf '13579\\n'\n",
        encoding="utf-8",
    )
    sbatch.chmod(sbatch.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["FAKE_CONDA_LOG"] = str(conda_log)
    environment["FAKE_SBATCH_LOG"] = str(sbatch_log)
    environment["BASH_COMPAT"] = "3.2"

    result = subprocess.run(
        (
            "bash",
            str(submitter),
            "--config",
            str(config),
            "--work-dir",
            str(work_dir),
            "--threads",
            "12",
            "--memory-mb",
            "96000",
        ),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert result.returncode == 0, result.stderr
    assert "controller job 13579" in result.stdout
    assert "env update --file" in conda_log.read_text(encoding="utf-8")
    arguments = sbatch_log.read_text(encoding="utf-8").splitlines()
    assert "--job-name=protein_signature_controller" in arguments
    assert "--account=barton" in arguments
    assert "--partition=barton" in arguments
    assert str(root / "slurm/protein_signature_workflow_controller.sbatch") in arguments
    assert "--profile" in arguments
    assert "slurm" in arguments
    assert "--memory-mb" in arguments
    assert "96000" in arguments
    assert (work_dir / "workflow_state/slurm").is_dir()
