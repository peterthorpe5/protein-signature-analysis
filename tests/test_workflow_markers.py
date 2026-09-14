"""Tests for deterministic Snakemake boundary markers."""

from __future__ import annotations

from pathlib import Path

import pytest

import protein_signatures.workflow_markers as marker_module
from protein_signatures.checksums import sha256_file
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.workflow_markers import (
    publish_validation_marker,
    publish_verification_marker,
    read_workflow_marker,
)


def test_validation_marker_records_config_and_validation_summary(
    tmp_path: Path,
    example_dir: Path,
) -> None:
    """Campaign validation should publish a checksum-bound deterministic marker."""

    config = example_dir / "campaign.yaml"
    marker = publish_validation_marker(
        config_path=config,
        marker_path=tmp_path / "state" / "VALIDATED.json",
    )
    document = read_workflow_marker(
        marker_path=marker,
        action="CAMPAIGN_VALIDATION",
    )
    assert document["config_path"] == str(config.resolve())
    assert document["config_sha256"] == sha256_file(path=config)
    assert document["validation_summary"]["status"] == "VALID"
    with pytest.raises(PublicationError, match="must not replace"):
        publish_validation_marker(config_path=config, marker_path=config)
    with pytest.raises(InputValidationError, match="action differs"):
        read_workflow_marker(marker_path=marker, action="WRONG_ACTION")


def test_verification_marker_checks_result_and_original_inputs(
    tmp_path: Path,
    completed_result: Path,
) -> None:
    """Post-run verification should bind both immutable outputs and live inputs."""

    marker = publish_verification_marker(
        result_dir=completed_result,
        marker_path=tmp_path / "workflow_state" / "VERIFIED.json",
    )
    document = read_workflow_marker(
        marker_path=marker,
        action="RESULT_AND_INPUT_VERIFICATION",
    )
    assert document["result_dir"] == str(completed_result.resolve())
    assert document["completion_sha256"] == sha256_file(path=completed_result / "COMPLETED.json")
    assert document["manifest_sha256"] == sha256_file(path=completed_result / "manifest.json")
    with pytest.raises(PublicationError, match="outside the immutable result"):
        publish_verification_marker(
            result_dir=completed_result,
            marker_path=completed_result / "workflow_marker.json",
        )


def test_workflow_marker_reader_rejects_invalid_shape_and_status(tmp_path: Path) -> None:
    """Marker consumers should fail closed on malformed or incomplete documents."""

    marker = tmp_path / "marker.json"
    marker.write_text("[]\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="contain an object"):
        read_workflow_marker(marker_path=marker, action="CAMPAIGN_VALIDATION")
    marker.write_text('{"schema_version": 1, "status": "RUNNING"}\n', encoding="utf-8")
    with pytest.raises(InputValidationError, match="not valid and complete"):
        read_workflow_marker(marker_path=marker, action="CAMPAIGN_VALIDATION")


def test_validation_marker_rejects_configuration_changed_in_flight(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The validation boundary must not publish against two different config digests."""

    config = tmp_path / "campaign.yaml"
    config.write_text("schema_version: 1\n", encoding="utf-8")
    digests = iter(("a" * 64, "b" * 64))
    monkeypatch.setattr(marker_module, "sha256_file", lambda **_kwargs: next(digests))
    monkeypatch.setattr(
        marker_module,
        "validate_campaign",
        lambda **_kwargs: {"status": "VALID"},
    )
    destination = tmp_path / "VALIDATED.json"
    with pytest.raises(InputValidationError, match="changed during"):
        publish_validation_marker(config_path=config, marker_path=destination)
    assert not destination.exists()
