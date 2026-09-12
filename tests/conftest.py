"""Shared deterministic fixtures for the protein-signature test suite."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from protein_signatures.pipeline import run_campaign


@pytest.fixture(name="example_dir")
def fixture_example_dir() -> Path:
    """Return the packaged offline example directory.

    Returns:
        Absolute example directory.
    """

    return (Path(__file__).parents[1] / "examples" / "minimal_e3").resolve()


@pytest.fixture(name="completed_result")
def fixture_completed_result(tmp_path: Path, example_dir: Path) -> Path:
    """Run and return one complete offline campaign.

    Args:
        tmp_path: Pytest temporary directory.
        example_dir: Packaged example directory.

    Returns:
        Completed result directory.
    """

    configuration = yaml.safe_load((example_dir / "campaign.yaml").read_text(encoding="utf-8"))
    for field in (
        "sequences_fasta",
        "label_assignments",
        "domains",
        "structures",
        "structure_comparisons",
    ):
        value = configuration["inputs"].get(field)
        if value:
            configuration["inputs"][field] = str(example_dir / value)
    configuration["explainable_ml"]["plot_cache_dir"] = str(tmp_path / "shap-cache")
    copied_config = tmp_path / "campaign.yaml"
    copied_config.write_text(yaml.safe_dump(configuration, sort_keys=False), encoding="utf-8")
    return run_campaign(
        config_path=copied_config,
        output_dir=tmp_path / "result",
        threads=1,
    )
