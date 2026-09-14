"""Checksummed boundary markers for external workflow orchestration."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from . import __version__
from .checksums import sha256_file
from .errors import InputValidationError, PublicationError
from .io_utils import read_json, write_json_atomic
from .pipeline import validate_campaign
from .publication import verify_completed_result, verify_input_authorities


def publish_validation_marker(*, config_path: Path, marker_path: Path) -> Path:
    """Validate a campaign and publish a deterministic orchestration marker.

    Args:
        config_path: Campaign YAML authority.
        marker_path: Marker path outside the campaign configuration.

    Returns:
        Published absolute marker path.

    Raises:
        InputValidationError: If configuration validation fails or changes in flight.
        PublicationError: If the marker destination is unsafe or cannot be written.
    """

    config = Path(config_path).expanduser().resolve()
    marker = Path(marker_path).expanduser().resolve()
    if marker == config:
        raise PublicationError("Workflow validation marker must not replace campaign YAML.")
    initial_digest = sha256_file(path=config)
    summary = validate_campaign(config_path=config)
    final_digest = sha256_file(path=config)
    if initial_digest != final_digest:
        raise InputValidationError("Campaign YAML changed during workflow validation.")
    write_json_atomic(
        path=marker,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "CAMPAIGN_VALIDATION",
            "package_version": __version__,
            "config_path": str(config),
            "config_sha256": final_digest,
            "validation_summary": summary,
        },
    )
    return marker


def publish_verification_marker(*, result_dir: Path, marker_path: Path) -> Path:
    """Verify an immutable result and its inputs, then publish a marker.

    Args:
        result_dir: Completed protein-signature result directory.
        marker_path: Marker path outside the immutable result directory.

    Returns:
        Published absolute marker path.

    Raises:
        PublicationError: If verification fails or the marker destination is unsafe.
    """

    result = Path(result_dir).expanduser().resolve()
    marker = Path(marker_path).expanduser().resolve()
    if marker == result or result in marker.parents:
        raise PublicationError("Workflow verification marker must be outside the immutable result.")
    verify_completed_result(result_dir=result)
    verify_input_authorities(result_dir=result)
    completion_path = result / "COMPLETED.json"
    manifest_path = result / "manifest.json"
    write_json_atomic(
        path=marker,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "RESULT_AND_INPUT_VERIFICATION",
            "package_version": __version__,
            "result_dir": str(result),
            "completion_sha256": sha256_file(path=completion_path),
            "manifest_sha256": sha256_file(path=manifest_path),
        },
    )
    return marker


def read_workflow_marker(*, marker_path: Path, action: str) -> dict[str, Any]:
    """Read one workflow marker and require its expected action and status.

    Args:
        marker_path: Existing JSON marker.
        action: Expected controlled action value.

    Returns:
        Validated marker mapping.

    Raises:
        InputValidationError: If marker shape, status or action is invalid.
    """

    marker = Path(marker_path).expanduser().resolve()
    document = read_json(path=marker)
    if not isinstance(document, dict):
        raise InputValidationError(f"Workflow marker must contain an object: {marker}")
    if document.get("schema_version") != 1 or document.get("status") != "VALID":
        raise InputValidationError(f"Workflow marker is not valid and complete: {marker}")
    if document.get("action") != action:
        raise InputValidationError(f"Workflow marker action differs from {action!r}: {marker}")
    return document
