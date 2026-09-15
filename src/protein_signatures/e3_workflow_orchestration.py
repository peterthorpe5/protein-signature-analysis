"""Reproducible orchestration boundaries for completed-E3 workflow inputs."""

from __future__ import annotations

import logging
import math
import os
import re
import shutil
import tempfile
from collections import Counter
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .automated_test_labels import (
    AUTOMATED_TEST_APPROVER,
    AUTOMATED_TEST_EVIDENCE_STATUS,
    AUTOMATED_TEST_TOKEN,
)
from .checksums import sha256_file
from .config import load_config
from .e3_workflow_bridge import prepare_e3_workflow_inputs
from .errors import InputValidationError, PublicationError
from .evidence_labels import (
    EVIDENCE_APPROVER,
    EVIDENCE_CONTROL_STATUS,
    EVIDENCE_INTERPRETATION_SCOPE,
    EVIDENCE_POSITIVE_STATUS,
    EVIDENCE_TARGET_STATUS,
    verify_evidence_label_bundle,
)
from .fasta import iter_protein_fasta
from .io_utils import iter_tsv, read_json, write_json_atomic
from .pipeline import validate_campaign
from .profiles import (
    default_profile_comparisons,
    load_profile,
    validate_assignment_profile_compatibility,
)
from .starter import initialise_campaign
from .tables import LABEL_FIELDS, read_label_assignments
from .validation import validate_text
from .workflow_markers import read_workflow_marker

LOGGER = logging.getLogger(__name__)
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_PREPARED_OUTPUTS = frozenset(
    {
        "alphafold_accessions.MISSING_MODELS_REVIEW_REQUIRED.tsv",
        "domains.tsv",
        "e3_label_curation_review.tsv",
        "label_assignments.REVIEW_REQUIRED.tsv",
        "proteins.faa",
        "source_inventory.tsv",
        "structures.tsv",
    }
)


def verify_prepared_e3_bundle(
    *,
    prepared_dir: Path,
    run_root: Path,
    minimum_mean_plddt: float,
    verify_source_authorities: bool = True,
) -> dict[str, Any]:
    """Verify a completed-E3 review bundle and its source identity.

    Args:
        prepared_dir: Published preparation directory.
        run_root: Expected completed predecessor root.
        minimum_mean_plddt: Expected structural eligibility threshold.
        verify_source_authorities: Rehash predecessor authorities when true.

    Returns:
        Validated preparation marker.

    Raises:
        InputValidationError: If any identity, inventory or checksum differs.
    """

    prepared = Path(prepared_dir).expanduser().resolve()
    source_root = Path(run_root).expanduser().resolve()
    if not prepared.is_dir():
        raise InputValidationError(f"Prepared E3 bundle is not a directory: {prepared}")
    marker_path = prepared / "PREPARED.json"
    document = read_json(path=marker_path)
    if not isinstance(document, dict):
        raise InputValidationError(f"Prepared marker must contain an object: {marker_path}")
    if document.get("schema_version") != 1 or document.get("status") != "COMPLETE":
        raise InputValidationError(f"Prepared marker is not complete: {marker_path}")
    recorded_root = Path(str(document.get("source_run_root", ""))).expanduser().resolve()
    if recorded_root != source_root:
        raise InputValidationError(
            f"Prepared source root differs from the requested run: {recorded_root} != {source_root}"
        )
    for field in ("structural_alignment_resource", "orthofinder_results"):
        authority = (
            Path(validate_text(value=document.get(field), field_name=f"prepared {field}"))
            .expanduser()
            .resolve()
        )
        if not authority.is_dir() or source_root not in authority.parents:
            raise InputValidationError(
                f"Prepared {field} must be an existing directory below the source run: {authority}"
            )
    try:
        threshold = float(document.get("minimum_mean_plddt"))
        expected_threshold = float(minimum_mean_plddt)
    except (TypeError, ValueError) as error:
        raise InputValidationError("Prepared mean-pLDDT threshold is not numeric.") from error
    if (
        not math.isfinite(threshold)
        or not math.isfinite(expected_threshold)
        or not math.isclose(threshold, expected_threshold, abs_tol=1e-9)
    ):
        raise InputValidationError(
            f"Prepared mean-pLDDT threshold differs: {threshold!r} != {expected_threshold!r}"
        )
    for field in (
        "protein_count",
        "pfam_record_count",
        "structure_count",
        "foldseek_eligible_structure_count",
    ):
        _require_nonnegative_integer(value=document.get(field), field_name=field)
    outputs = document.get("outputs")
    if not isinstance(outputs, list):
        raise InputValidationError("Prepared output inventory must be a list.")
    observed: set[str] = set()
    for index, value in enumerate(outputs):
        if not isinstance(value, Mapping):
            raise InputValidationError(f"Malformed prepared output record at index {index}.")
        relative_text = str(value.get("relative_path", ""))
        relative = Path(relative_text)
        if (
            not relative_text
            or relative.is_absolute()
            or ".." in relative.parts
            or len(relative.parts) != 1
        ):
            raise InputValidationError(f"Unsafe prepared output path: {relative_text!r}")
        if relative_text in observed:
            raise InputValidationError(f"Duplicate prepared output record: {relative_text!r}")
        observed.add(relative_text)
        expected_size = _require_nonnegative_integer(
            value=value.get("size_bytes"),
            field_name=f"outputs[{index}].size_bytes",
        )
        expected_digest = _require_sha256(
            value=value.get("sha256"),
            field_name=f"outputs[{index}].sha256",
        )
        candidate = prepared / relative
        if not candidate.is_file() or candidate.stat().st_size != expected_size:
            raise InputValidationError(f"Prepared output size differs: {candidate}")
        if sha256_file(path=candidate) != expected_digest:
            raise InputValidationError(f"Prepared output checksum differs: {candidate}")
    if observed != _PREPARED_OUTPUTS:
        raise InputValidationError(
            "Prepared output inventory differs from the supported contract: "
            f"missing={sorted(_PREPARED_OUTPUTS - observed)} "
            f"unexpected={sorted(observed - _PREPARED_OUTPUTS)}"
        )
    if verify_source_authorities:
        _verify_source_inventory(
            inventory_path=prepared / "source_inventory.tsv",
            run_root=source_root,
        )
    LOGGER.info(
        "Verified completed-E3 prepared bundle proteins=%s structures=%s eligible=%s at %s",
        document["protein_count"],
        document["structure_count"],
        document["foldseek_eligible_structure_count"],
        prepared,
    )
    return document


def ensure_e3_preparation_marker(
    *,
    run_root: Path,
    prepared_dir: Path,
    minimum_mean_plddt: float,
    marker_path: Path,
) -> Path:
    """Create or verify preparation and publish a Snakemake-owned marker.

    Args:
        run_root: Completed E3 predecessor root.
        prepared_dir: Stable preparation directory.
        minimum_mean_plddt: Structural eligibility threshold.
        marker_path: Marker outside the immutable prepared bundle.

    Returns:
        Published marker path.
    """

    prepared = Path(prepared_dir).expanduser().resolve()
    marker = Path(marker_path).expanduser().resolve()
    if prepared == marker or prepared in marker.parents:
        raise PublicationError("Preparation workflow marker must be outside prepared inputs.")
    action = "REUSED_VERIFIED" if prepared.exists() else "CREATED"
    if not prepared.exists():
        prepare_e3_workflow_inputs(
            run_root=run_root,
            output_dir=prepared,
            minimum_mean_plddt=minimum_mean_plddt,
        )
    document = verify_prepared_e3_bundle(
        prepared_dir=prepared,
        run_root=run_root,
        minimum_mean_plddt=minimum_mean_plddt,
        verify_source_authorities=True,
    )
    write_json_atomic(
        path=marker,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "E3_PREPARATION_VERIFICATION",
            "package_version": __version__,
            "preparation_action": action,
            "prepared_dir": str(prepared),
            "prepared_marker_sha256": sha256_file(path=prepared / "PREPARED.json"),
            "source_run_root": str(Path(run_root).expanduser().resolve()),
            "minimum_mean_plddt": float(minimum_mean_plddt),
            "protein_count": document["protein_count"],
            "structure_count": document["structure_count"],
            "foldseek_eligible_structure_count": document["foldseek_eligible_structure_count"],
        },
    )
    LOGGER.info("Published E3 preparation boundary action=%s at %s", action, marker)
    return marker


def stage_e3_label_review(
    *, preparation_marker: Path, reviewed_labels: Path, marker_path: Path
) -> Path:
    """Safely stage or adopt a label-review file without overwriting it.

    Args:
        preparation_marker: Verified preparation boundary marker.
        reviewed_labels: Editable review authority outside prepared inputs.
        marker_path: Snakemake-owned review-ready marker.

    Returns:
        Published review-ready marker.
    """

    preparation = read_workflow_marker(
        marker_path=preparation_marker,
        action="E3_PREPARATION_VERIFICATION",
    )
    prepared = Path(str(preparation.get("prepared_dir", ""))).resolve()
    if sha256_file(path=prepared / "PREPARED.json") != preparation.get("prepared_marker_sha256"):
        raise InputValidationError("Prepared marker changed after workflow verification.")
    template = prepared / "label_assignments.REVIEW_REQUIRED.tsv"
    reviewed = Path(reviewed_labels).expanduser().resolve()
    marker = Path(marker_path).expanduser().resolve()
    if reviewed == template or prepared in reviewed.parents:
        raise PublicationError("Reviewed labels must be outside the immutable prepared bundle.")
    if marker == reviewed:
        raise PublicationError("Review-ready marker must not replace reviewed labels.")
    if reviewed.exists() and not reviewed.is_file():
        raise PublicationError(f"Reviewed-label destination is not a regular file: {reviewed}")
    staging_action = "ADOPTED_EXISTING"
    if not reviewed.exists():
        _copy_file_atomic(source=template, destination=reviewed)
        staging_action = "CREATED_FROM_TEMPLATE"
    row_count = sum(1 for _ in iter_tsv(path=reviewed, required_fields=LABEL_FIELDS))
    template_digest = sha256_file(path=template)
    reviewed_digest = sha256_file(path=reviewed)
    write_json_atomic(
        path=marker,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "E3_LABEL_REVIEW_STAGING",
            "package_version": __version__,
            "review_status": "AWAITING_CURATOR_APPROVAL",
            "staging_action": staging_action,
            "prepared_dir": str(prepared),
            "preparation_marker": str(Path(preparation_marker).expanduser().resolve()),
            "preparation_marker_sha256": sha256_file(
                path=Path(preparation_marker).expanduser().resolve()
            ),
            "template_path": str(template),
            "template_sha256": template_digest,
            "reviewed_labels": str(reviewed),
            "reviewed_sha256_at_staging": reviewed_digest,
            "reviewed_row_count_at_staging": row_count,
            "reviewed_differs_from_template": reviewed_digest != template_digest,
        },
    )
    LOGGER.info(
        "Published E3 label-review staging action=%s rows=%d differs_from_template=%s at %s",
        staging_action,
        row_count,
        reviewed_digest != template_digest,
        marker,
    )
    return marker


def approve_e3_label_review(
    *,
    preparation_marker: Path,
    review_marker: Path,
    reviewed_labels: Path,
    approval_marker: Path,
    curator: str,
    note: str = "",
    profile: str | Path = "e3",
    automated_test_mode: bool = False,
    automated_test_marker: Path | None = None,
    automated_evidence_mode: bool = False,
    evidence_label_marker: Path | None = None,
) -> Path:
    """Validate reviewed labels and publish an immutable curator approval.

    Args:
        preparation_marker: Verified preparation marker.
        review_marker: Review-staging marker.
        reviewed_labels: Completed label-assignment authority.
        approval_marker: New checksum-bound approval marker.
        curator: Named person accepting the reviewed authority.
        note: Optional bounded approval note.
        profile: Built-in or custom label profile.
        automated_test_mode: Accept only explicitly synthetic software-test labels.
        automated_test_marker: Required generation audit in automated test mode.
        automated_evidence_mode: Accept a checksummed provisional evidence bundle.
        evidence_label_marker: Required evidence-bundle marker in evidence mode.

    Returns:
        Published approval marker.

    Raises:
        InputValidationError: If labels remain unsuitable for analysis.
        PublicationError: If approval would overwrite an existing marker.
    """

    approval = Path(approval_marker).expanduser().resolve()
    if approval.exists():
        raise PublicationError(f"Label-review approval already exists: {approval}")
    curator_name = validate_text(
        value=curator,
        field_name="curator",
        maximum_length=256,
    )
    approval_note = validate_text(
        value=note,
        field_name="approval note",
        allow_empty=True,
        maximum_length=2_000,
    )
    context = _review_context(
        preparation_marker=preparation_marker,
        review_marker=review_marker,
        reviewed_labels=reviewed_labels,
    )
    if context["reviewed_sha256"] == context["template_sha256"]:
        raise InputValidationError(
            "Reviewed labels are byte-identical to the all-UNMAPPED template."
        )
    summary = _validate_reviewed_assignments(
        prepared_dir=context["prepared_dir"],
        reviewed_labels=context["reviewed_labels"],
        profile=profile,
    )
    synthetic_count = summary["synthetic_test_positive_count"]
    evidence_count = summary["evidence_supported_positive_count"]
    automated_authority: dict[str, Any] = {}
    if automated_test_mode and automated_evidence_mode:
        raise InputValidationError(
            "Automated smoke-test and evidence-led approval modes are mutually exclusive."
        )
    if automated_test_mode:
        if evidence_label_marker is not None:
            raise InputValidationError(
                "Evidence-label markers cannot enter automated smoke-test approval."
            )
        if curator_name != AUTOMATED_TEST_APPROVER:
            raise InputValidationError(
                "Automated test approval requires curator='AUTOMATED_TEST_MODE'."
            )
        if synthetic_count != summary["reviewed_positive_count"]:
            raise InputValidationError(
                "Automated test approval requires every positive assignment to carry "
                f"evidence_status={AUTOMATED_TEST_EVIDENCE_STATUS!r}."
            )
        test_marker_path, test_marker_sha256 = _validate_automated_test_marker(
            marker_path=automated_test_marker,
            prepared_dir=context["prepared_dir"],
            reviewed_labels=context["reviewed_labels"],
            reviewed_sha256=context["reviewed_sha256"],
            profile=profile,
        )
        automated_authority = {
            "automated_test_marker": str(test_marker_path),
            "automated_test_marker_sha256": test_marker_sha256,
        }
        approval_mode = "AUTOMATED_SMOKE_TEST"
        review_status = "AUTOMATED_TEST_ONLY"
        interpretation_allowed = False
        interpretation_scope = "SOFTWARE_EXECUTION_TEST_ONLY"
    elif automated_evidence_mode:
        if automated_test_marker is not None:
            raise InputValidationError("Automated test markers cannot enter evidence-led approval.")
        if curator_name != EVIDENCE_APPROVER:
            raise InputValidationError(
                f"Automated evidence approval requires curator={EVIDENCE_APPROVER!r}."
            )
        if synthetic_count:
            raise InputValidationError(
                "Synthetic test labels cannot enter an evidence-led approval."
            )
        if evidence_count != summary["reviewed_positive_count"]:
            raise InputValidationError(
                "Evidence-led approval requires every analysis-positive assignment to use "
                f"curation_status={EVIDENCE_POSITIVE_STATUS!r}."
            )
        if (
            summary["evidence_supported_target_count"] + summary["evidence_supported_control_count"]
            != evidence_count
        ):
            raise InputValidationError(
                "Evidence-led positives must carry the controlled automated target or "
                "matched-control evidence status."
            )
        marker_path, marker_digest, evidence_document = _validate_evidence_label_marker(
            marker_path=evidence_label_marker,
            prepared_dir=context["prepared_dir"],
            reviewed_labels=context["reviewed_labels"],
            reviewed_sha256=context["reviewed_sha256"],
            profile=profile,
        )
        automated_authority = {
            "evidence_label_marker": str(marker_path),
            "evidence_label_marker_sha256": marker_digest,
            "evidence_ruleset_id": evidence_document["ruleset_id"],
            "evidence_ruleset_version": evidence_document["ruleset_version"],
            "evidence_ruleset_sha256": evidence_document["ruleset_sha256"],
        }
        approval_mode = "AUTOMATED_EVIDENCE_PROPOSAL"
        review_status = "PROVISIONAL_EVIDENCE_SUPPORTED"
        interpretation_allowed = False
        interpretation_scope = EVIDENCE_INTERPRETATION_SCOPE
    else:
        if automated_test_marker is not None or evidence_label_marker is not None:
            raise InputValidationError(
                "Automated generation markers are valid only in their explicit approval mode."
            )
        if synthetic_count:
            raise InputValidationError(
                "Synthetic test labels cannot receive human-review approval; use the "
                "explicit automated test workflow in an isolated smoke-test campaign."
            )
        if evidence_count:
            raise InputValidationError(
                "Automated evidence-supported labels require the explicit provisional "
                "evidence approval mode, or human review must replace their curation status "
                "with REVIEWED_POSITIVE."
            )
        approval_mode = "HUMAN_REVIEW"
        review_status = "APPROVED"
        interpretation_allowed = True
        interpretation_scope = "HUMAN_REVIEWED_ANALYSIS"
    write_json_atomic(
        path=approval,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "E3_LABEL_REVIEW_APPROVAL",
            "package_version": __version__,
            "review_status": review_status,
            "approval_mode": approval_mode,
            "scientific_interpretation_allowed": interpretation_allowed,
            "interpretation_scope": interpretation_scope,
            "approved_at_utc": datetime.now(timezone.utc).isoformat(),
            "curator": curator_name,
            "approval_note": approval_note,
            "prepared_dir": str(context["prepared_dir"]),
            "preparation_marker": str(context["preparation_marker"]),
            "preparation_marker_sha256": context["preparation_marker_sha256"],
            "review_marker": str(context["review_marker"]),
            "review_marker_sha256": context["review_marker_sha256"],
            "template_sha256": context["template_sha256"],
            "reviewed_labels": str(context["reviewed_labels"]),
            "reviewed_labels_size_bytes": context["reviewed_labels"].stat().st_size,
            "reviewed_labels_sha256": context["reviewed_sha256"],
            **automated_authority,
            **summary,
        },
    )
    LOGGER.info(
        "Published checksum-bound E3 label approval curator=%s positives=%d at %s",
        curator_name,
        summary["reviewed_positive_count"],
        approval,
    )
    return approval


def verify_e3_label_review(
    *,
    preparation_marker: Path,
    review_marker: Path,
    reviewed_labels: Path,
    approval_marker: Path,
    marker_path: Path,
    profile: str | Path = "e3",
) -> Path:
    """Verify current reviewed labels against their immutable approval.

    Args:
        preparation_marker: Verified preparation marker.
        review_marker: Review-staging marker.
        reviewed_labels: Current label authority.
        approval_marker: Existing curator approval.
        marker_path: Snakemake-owned verification marker.
        profile: Built-in or custom profile.

    Returns:
        Published verification marker.
    """

    approval = read_workflow_marker(
        marker_path=approval_marker,
        action="E3_LABEL_REVIEW_APPROVAL",
    )
    context = _review_context(
        preparation_marker=preparation_marker,
        review_marker=review_marker,
        reviewed_labels=reviewed_labels,
    )
    expected = {
        "prepared_dir": str(context["prepared_dir"]),
        "preparation_marker_sha256": context["preparation_marker_sha256"],
        "review_marker_sha256": context["review_marker_sha256"],
        "template_sha256": context["template_sha256"],
        "reviewed_labels": str(context["reviewed_labels"]),
        "reviewed_labels_size_bytes": context["reviewed_labels"].stat().st_size,
        "reviewed_labels_sha256": context["reviewed_sha256"],
    }
    for field, value in expected.items():
        if approval.get(field) != value:
            raise InputValidationError(
                f"Label-review approval field {field!r} no longer matches its authority."
            )
    summary = _validate_reviewed_assignments(
        prepared_dir=context["prepared_dir"],
        reviewed_labels=context["reviewed_labels"],
        profile=profile,
    )
    for field, value in summary.items():
        if approval.get(field) != value:
            raise InputValidationError(
                f"Label-review approval summary {field!r} no longer matches reviewed labels."
            )
    approval_mode = approval.get("approval_mode", "HUMAN_REVIEW")
    if approval_mode not in {
        "HUMAN_REVIEW",
        "AUTOMATED_SMOKE_TEST",
        "AUTOMATED_EVIDENCE_PROPOSAL",
    }:
        raise InputValidationError(f"Unknown E3 label approval mode: {approval_mode!r}")
    interpretation_allowed = approval_mode == "HUMAN_REVIEW"
    if (
        approval.get("scientific_interpretation_allowed", interpretation_allowed)
        != interpretation_allowed
    ):
        raise InputValidationError("Label-review approval interpretation policy is inconsistent.")
    expected_scope = {
        "HUMAN_REVIEW": "HUMAN_REVIEWED_ANALYSIS",
        "AUTOMATED_SMOKE_TEST": "SOFTWARE_EXECUTION_TEST_ONLY",
        "AUTOMATED_EVIDENCE_PROPOSAL": EVIDENCE_INTERPRETATION_SCOPE,
    }[approval_mode]
    if approval.get("interpretation_scope", expected_scope) != expected_scope:
        raise InputValidationError("Label-review approval interpretation scope is inconsistent.")
    if approval_mode == "AUTOMATED_SMOKE_TEST":
        test_marker_path, test_marker_sha256 = _validate_automated_test_marker(
            marker_path=Path(str(approval.get("automated_test_marker", ""))),
            prepared_dir=context["prepared_dir"],
            reviewed_labels=context["reviewed_labels"],
            reviewed_sha256=context["reviewed_sha256"],
            profile=profile,
        )
        if approval.get("automated_test_marker_sha256") != test_marker_sha256:
            raise InputValidationError("Automated test-label marker changed after approval.")
        evidence_marker_path = None
        evidence_marker_sha256 = None
    elif approval_mode == "AUTOMATED_EVIDENCE_PROPOSAL":
        evidence_marker_path, evidence_marker_sha256, evidence_document = (
            _validate_evidence_label_marker(
                marker_path=Path(str(approval.get("evidence_label_marker", ""))),
                prepared_dir=context["prepared_dir"],
                reviewed_labels=context["reviewed_labels"],
                reviewed_sha256=context["reviewed_sha256"],
                profile=profile,
            )
        )
        if approval.get("evidence_label_marker_sha256") != evidence_marker_sha256:
            raise InputValidationError("Evidence-label marker changed after approval.")
        expected_rules = {
            "evidence_ruleset_id": evidence_document["ruleset_id"],
            "evidence_ruleset_version": evidence_document["ruleset_version"],
            "evidence_ruleset_sha256": evidence_document["ruleset_sha256"],
        }
        for field, expected_value in expected_rules.items():
            if approval.get(field) != expected_value:
                raise InputValidationError(
                    f"Evidence approval field {field!r} differs from its bundle."
                )
        test_marker_path = None
        test_marker_sha256 = None
    else:
        test_marker_path = None
        test_marker_sha256 = None
        evidence_marker_path = None
        evidence_marker_sha256 = None
    marker = Path(marker_path).expanduser().resolve()
    write_json_atomic(
        path=marker,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "E3_LABEL_REVIEW_VERIFICATION",
            "package_version": __version__,
            "approval_marker": str(Path(approval_marker).expanduser().resolve()),
            "approval_marker_sha256": sha256_file(
                path=Path(approval_marker).expanduser().resolve()
            ),
            "prepared_dir": str(context["prepared_dir"]),
            "preparation_marker": str(context["preparation_marker"]),
            "preparation_marker_sha256": context["preparation_marker_sha256"],
            "review_marker": str(context["review_marker"]),
            "review_marker_sha256": context["review_marker_sha256"],
            "reviewed_labels": str(context["reviewed_labels"]),
            "reviewed_labels_sha256": context["reviewed_sha256"],
            "approval_mode": approval_mode,
            "scientific_interpretation_allowed": interpretation_allowed,
            "automated_test_marker": (
                str(test_marker_path) if test_marker_path is not None else None
            ),
            "automated_test_marker_sha256": test_marker_sha256,
            "evidence_label_marker": (
                str(evidence_marker_path) if evidence_marker_path is not None else None
            ),
            "evidence_label_marker_sha256": evidence_marker_sha256,
            "interpretation_scope": approval.get(
                "interpretation_scope",
                "HUMAN_REVIEWED_ANALYSIS" if interpretation_allowed else "RESTRICTED",
            ),
            **summary,
        },
    )
    LOGGER.info(
        "Verified approved E3 labels assignments=%d checksum=%s at %s",
        summary["reviewed_assignment_count"],
        context["reviewed_sha256"],
        marker,
    )
    return marker


def ensure_e3_campaign_marker(
    *,
    preparation_marker: Path,
    review_verification_marker: Path,
    campaign_config: Path,
    campaign_id: str,
    profile: str | Path,
    marker_path: Path,
) -> Path:
    """Create or adopt the E3 campaign and publish its validation marker.

    Args:
        preparation_marker: Verified preparation boundary.
        review_verification_marker: Verified curator-approval boundary.
        campaign_config: Stable campaign YAML path.
        campaign_id: Expected campaign identifier.
        profile: Built-in or custom profile.
        marker_path: Snakemake-owned validation marker.

    Returns:
        Published campaign-validation marker.
    """

    preparation = read_workflow_marker(
        marker_path=preparation_marker,
        action="E3_PREPARATION_VERIFICATION",
    )
    review = read_workflow_marker(
        marker_path=review_verification_marker,
        action="E3_LABEL_REVIEW_VERIFICATION",
    )
    prepared = Path(str(preparation.get("prepared_dir", ""))).resolve()
    labels = Path(str(review.get("reviewed_labels", ""))).resolve()
    preparation_path = Path(preparation_marker).expanduser().resolve()
    if review.get("prepared_dir") != str(prepared):
        raise InputValidationError("Review verification refers to a different prepared bundle.")
    if review.get("preparation_marker") != str(preparation_path):
        raise InputValidationError("Review verification refers to a different preparation marker.")
    if review.get("preparation_marker_sha256") != sha256_file(path=preparation_path):
        raise InputValidationError("Preparation marker changed after label-review verification.")
    approval_path = Path(str(review.get("approval_marker", ""))).resolve()
    if review.get("approval_marker_sha256") != sha256_file(path=approval_path):
        raise InputValidationError("Label-review approval changed after verification.")
    prepared_marker = prepared / "PREPARED.json"
    if sha256_file(path=prepared_marker) != preparation.get("prepared_marker_sha256"):
        raise InputValidationError("Prepared marker changed after workflow verification.")
    if sha256_file(path=labels) != review.get("reviewed_labels_sha256"):
        raise InputValidationError("Reviewed labels changed after approval verification.")
    prepared_document = read_json(path=prepared_marker)
    if not isinstance(prepared_document, dict):
        raise InputValidationError("Prepared marker must contain an object.")
    eligible_count = _require_nonnegative_integer(
        value=prepared_document.get("foldseek_eligible_structure_count"),
        field_name="foldseek_eligible_structure_count",
    )
    if eligible_count < 2:
        raise InputValidationError("At least two Foldseek-eligible structures are required.")
    config_path = Path(campaign_config).expanduser().resolve()
    approval_mode = str(review.get("approval_mode", "HUMAN_REVIEW"))
    if approval_mode == "AUTOMATED_SMOKE_TEST":
        if re.search(r"(?:smoke|test)", campaign_id, flags=re.IGNORECASE) is None:
            raise InputValidationError(
                "An automated smoke-test approval requires 'smoke' or 'test' in campaign_id."
            )
        if re.search(r"(?:smoke|test)", config_path.parent.name, flags=re.IGNORECASE) is None:
            raise InputValidationError(
                "An automated smoke-test approval requires an isolated work directory whose "
                "name contains 'smoke' or 'test'."
            )
    analysis_domains = prepared / "domains.tsv"
    evidence_inputs: dict[str, Path | None] = {
        "label_evidence_marker": None,
        "label_evidence_audit": None,
        "control_matching_audit": None,
        "label_definition_features": None,
        "class_labelling_summary": None,
        "unresolved_assignments": None,
    }
    if approval_mode == "AUTOMATED_EVIDENCE_PROPOSAL":
        evidence_marker, evidence_digest, evidence_document = _validate_evidence_label_marker(
            marker_path=Path(str(review.get("evidence_label_marker", ""))),
            prepared_dir=prepared,
            reviewed_labels=labels,
            reviewed_sha256=str(review.get("reviewed_labels_sha256", "")),
            profile=profile,
        )
        if review.get("evidence_label_marker_sha256") != evidence_digest:
            raise InputValidationError("Evidence-label bundle changed after review verification.")
        recorded_domains = evidence_document.get("analysis_domains")
        if not recorded_domains:
            raise InputValidationError(
                "Completed-E3 evidence labelling requires a label-definition-safe "
                "domain projection."
            )
        analysis_domains = Path(str(recorded_domains)).resolve()
        evidence_inputs = {
            "label_evidence_marker": evidence_marker,
            "label_evidence_audit": evidence_marker.parent / "label_evidence_audit.tsv",
            "control_matching_audit": evidence_marker.parent / "control_matching_audit.tsv",
            "label_definition_features": (evidence_marker.parent / "label_definition_features.tsv"),
            "class_labelling_summary": (evidence_marker.parent / "class_labelling_summary.tsv"),
            "unresolved_assignments": (evidence_marker.parent / "unresolved_assignments.tsv"),
        }
    expected_profile = _normalise_profile_source(profile=profile)
    creation_action = "ADOPTED_EXISTING" if config_path.exists() else "CREATED"
    if not config_path.exists():
        initialise_campaign(
            config_path=config_path,
            campaign_id=campaign_id,
            profile=expected_profile,
            sequences_fasta=prepared / "proteins.faa",
            label_assignments=labels,
            domains=analysis_domains,
            **evidence_inputs,
            structures=prepared / "structures.tsv",
            structural_alignment_resource=Path(
                str(prepared_document.get("structural_alignment_resource", ""))
            ),
            orthofinder_results=Path(str(prepared_document.get("orthofinder_results", ""))),
            orthofinder_group_type="HOG",
            orthofinder_hierarchy_node="N0",
            enable_foldseek=True,
            foldseek_maximum_hits=eligible_count,
        )
    _validate_existing_e3_campaign(
        config_path=config_path,
        campaign_id=campaign_id,
        profile=expected_profile,
        prepared_dir=prepared,
        reviewed_labels=labels,
        domains=analysis_domains,
        evidence_inputs=evidence_inputs,
        structural_resource=Path(str(prepared_document.get("structural_alignment_resource", ""))),
        orthofinder_results=Path(str(prepared_document.get("orthofinder_results", ""))),
        foldseek_maximum_hits=eligible_count,
    )
    initial_digest = sha256_file(path=config_path)
    summary = validate_campaign(config_path=config_path)
    final_digest = sha256_file(path=config_path)
    if initial_digest != final_digest:
        raise InputValidationError("Campaign YAML changed during E3 workflow validation.")
    marker = Path(marker_path).expanduser().resolve()
    write_json_atomic(
        path=marker,
        value={
            "schema_version": 1,
            "status": "VALID",
            "action": "E3_CAMPAIGN_VALIDATION",
            "package_version": __version__,
            "campaign_action": creation_action,
            "campaign_config": str(config_path),
            "campaign_config_sha256": final_digest,
            "preparation_marker_sha256": sha256_file(path=preparation_path),
            "review_verification_marker_sha256": sha256_file(
                path=Path(review_verification_marker).expanduser().resolve()
            ),
            "approval_mode": approval_mode,
            "scientific_interpretation_allowed": review.get(
                "scientific_interpretation_allowed", True
            ),
            "interpretation_scope": review.get("interpretation_scope", "HUMAN_REVIEWED_ANALYSIS"),
            "validation_summary": summary,
        },
    )
    LOGGER.info(
        "Published E3 campaign boundary action=%s config=%s marker=%s",
        creation_action,
        config_path,
        marker,
    )
    return marker


def _verify_source_inventory(*, inventory_path: Path, run_root: Path) -> None:
    """Verify every predecessor authority recorded by preparation."""

    observed: set[str] = set()
    count = 0
    for row in iter_tsv(
        path=inventory_path,
        required_fields=("authority", "path", "size_bytes", "sha256"),
    ):
        authority = validate_text(value=row["authority"], field_name="authority")
        if authority in observed:
            raise InputValidationError(f"Duplicate prepared source authority: {authority!r}")
        observed.add(authority)
        source = Path(row["path"]).expanduser().resolve()
        if run_root not in source.parents:
            raise InputValidationError(f"Prepared source authority is outside run root: {source}")
        expected_size = _require_nonnegative_integer(
            value=row["size_bytes"],
            field_name=f"{authority}.size_bytes",
        )
        expected_digest = _require_sha256(
            value=row["sha256"],
            field_name=f"{authority}.sha256",
        )
        if not source.is_file() or source.stat().st_size != expected_size:
            raise InputValidationError(f"Prepared source authority size differs: {source}")
        if sha256_file(path=source) != expected_digest:
            raise InputValidationError(f"Prepared source authority checksum differs: {source}")
        count += 1
    if count < 6:
        raise InputValidationError("Prepared source inventory is unexpectedly incomplete.")


def _review_context(
    *, preparation_marker: Path, review_marker: Path, reviewed_labels: Path
) -> dict[str, Any]:
    """Validate linked preparation/review markers and return current identities."""

    preparation_path = Path(preparation_marker).expanduser().resolve()
    review_path = Path(review_marker).expanduser().resolve()
    preparation = read_workflow_marker(
        marker_path=preparation_path,
        action="E3_PREPARATION_VERIFICATION",
    )
    review = read_workflow_marker(
        marker_path=review_path,
        action="E3_LABEL_REVIEW_STAGING",
    )
    prepared = Path(str(preparation.get("prepared_dir", ""))).resolve()
    reviewed = Path(reviewed_labels).expanduser().resolve()
    if review.get("prepared_dir") != str(prepared):
        raise InputValidationError("Review marker refers to a different prepared bundle.")
    if review.get("preparation_marker") != str(preparation_path):
        raise InputValidationError("Review marker refers to a different preparation marker.")
    preparation_digest = sha256_file(path=preparation_path)
    if review.get("preparation_marker_sha256") != preparation_digest:
        raise InputValidationError("Preparation marker changed after review staging.")
    if review.get("reviewed_labels") != str(reviewed):
        raise InputValidationError("Review marker refers to a different reviewed-label file.")
    template = prepared / "label_assignments.REVIEW_REQUIRED.tsv"
    template_digest = sha256_file(path=template)
    if review.get("template_sha256") != template_digest:
        raise InputValidationError("Review template changed after staging.")
    if not reviewed.is_file() or reviewed.stat().st_size == 0:
        raise InputValidationError(f"Reviewed labels are missing or empty: {reviewed}")
    return {
        "prepared_dir": prepared,
        "preparation_marker": preparation_path,
        "preparation_marker_sha256": preparation_digest,
        "review_marker": review_path,
        "review_marker_sha256": sha256_file(path=review_path),
        "template_sha256": template_digest,
        "reviewed_labels": reviewed,
        "reviewed_sha256": sha256_file(path=reviewed),
    }


def _validate_automated_test_marker(
    *,
    marker_path: Path | None,
    prepared_dir: Path,
    reviewed_labels: Path,
    reviewed_sha256: str,
    profile: str | Path,
) -> tuple[Path, str]:
    """Validate and bind one synthetic-label generation authority.

    Args:
        marker_path: Explicit automated label-generation marker.
        prepared_dir: Verified completed-E3 input bundle.
        reviewed_labels: Synthetic label table being approved.
        reviewed_sha256: Current checksum of the synthetic label table.
        profile: Built-in or custom classification profile.

    Returns:
        Resolved marker path and its current SHA-256 digest.

    Raises:
        InputValidationError: If the audit is absent, inconsistent or mutable.
    """

    if marker_path is None:
        raise InputValidationError("Automated test approval requires --automated-test-marker.")
    marker = Path(marker_path).expanduser().resolve()
    if AUTOMATED_TEST_TOKEN not in marker.name.upper():
        raise InputValidationError(
            f"Automated test marker filename must contain {AUTOMATED_TEST_TOKEN!r}."
        )
    document = read_json(path=marker)
    if not isinstance(document, dict):
        raise InputValidationError(f"Automated test-label marker must contain an object: {marker}")
    if (
        document.get("schema_version") != 1
        or document.get("status") != "AUTOMATED_TEST_ONLY"
        or document.get("action") != "AUTOMATED_TEST_LABEL_GENERATION"
        or document.get("scientific_interpretation_allowed") is not False
    ):
        raise InputValidationError(
            f"Automated test-label marker is not a complete test-only authority: {marker}"
        )
    expected_paths = {
        "sequence_fasta": prepared_dir / "proteins.faa",
        "template_labels": prepared_dir / "label_assignments.REVIEW_REQUIRED.tsv",
        "output_labels": reviewed_labels,
    }
    for field, expected_path in expected_paths.items():
        recorded = Path(str(document.get(field, ""))).expanduser().resolve()
        if recorded != expected_path.resolve():
            raise InputValidationError(
                f"Automated test-label marker {field!r} refers to another authority."
            )
    expected_digests = {
        "sequence_fasta_sha256": sha256_file(path=expected_paths["sequence_fasta"]),
        "template_labels_sha256": sha256_file(path=expected_paths["template_labels"]),
        "output_labels_sha256": reviewed_sha256,
    }
    for field, expected_digest in expected_digests.items():
        if document.get(field) != expected_digest:
            raise InputValidationError(
                f"Automated test-label marker {field!r} no longer matches its authority."
            )
    loaded_profile = load_profile(source=profile)
    if (
        document.get("profile_id") != loaded_profile.profile_id
        or document.get("profile_version") != loaded_profile.profile_version
    ):
        raise InputValidationError(
            "Automated test-label marker refers to a different classification profile."
        )
    return marker, sha256_file(path=marker)


def _validate_evidence_label_marker(
    *,
    marker_path: Path | None,
    prepared_dir: Path,
    reviewed_labels: Path,
    reviewed_sha256: str,
    profile: str | Path,
) -> tuple[Path, str, Mapping[str, Any]]:
    """Validate and bind a provisional evidence-label bundle.

    Args:
        marker_path: Explicit ``EVIDENCE_LABELS.json`` authority.
        prepared_dir: Verified completed-E3 input bundle.
        reviewed_labels: Evidence-supported assignment table being approved.
        reviewed_sha256: Current checksum of the assignment table.
        profile: Built-in or custom classification profile.

    Returns:
        Marker path, marker digest and validated marker document.

    Raises:
        InputValidationError: If paths, inputs, profile or checksums disagree.
    """

    if marker_path is None:
        raise InputValidationError("Automated evidence approval requires --evidence-label-marker.")
    marker = Path(marker_path).expanduser().resolve()
    if marker.name != "EVIDENCE_LABELS.json":
        raise InputValidationError(
            f"Evidence-label authority must be named EVIDENCE_LABELS.json: {marker}"
        )
    document = verify_evidence_label_bundle(bundle_dir=marker.parent)
    recorded_labels = Path(str(document.get("label_assignments", ""))).resolve()
    if recorded_labels != reviewed_labels.resolve():
        raise InputValidationError("Evidence-label marker refers to a different assignment table.")
    if document.get("label_assignments_sha256") != reviewed_sha256:
        raise InputValidationError(
            "Evidence-label assignment checksum differs from the reviewed authority."
        )
    loaded_profile = load_profile(source=profile)
    if (
        document.get("profile_id") != loaded_profile.profile_id
        or document.get("profile_version") != loaded_profile.profile_version
    ):
        raise InputValidationError(
            "Evidence-label marker refers to a different classification profile."
        )
    expected_inputs = {
        prepared_dir / "proteins.faa",
        prepared_dir / "domains.tsv",
        prepared_dir / "structures.tsv",
        prepared_dir / "e3_label_curation_review.tsv",
        prepared_dir / "label_assignments.REVIEW_REQUIRED.tsv",
    }
    recorded_inputs = {
        Path(str(row.get("path", ""))).resolve(): row
        for row in document.get("inputs", ())
        if isinstance(row, Mapping)
    }
    missing = expected_inputs - set(recorded_inputs)
    if missing:
        raise InputValidationError(
            "Evidence-label marker is not bound to every required prepared authority: "
            f"{sorted(str(path) for path in missing)}"
        )
    for path in expected_inputs:
        row = recorded_inputs[path]
        if row.get("size_bytes") != path.stat().st_size:
            raise InputValidationError(f"Prepared evidence input size differs: {path}")
        if row.get("sha256") != sha256_file(path=path):
            raise InputValidationError(f"Prepared evidence input checksum differs: {path}")
    return marker, sha256_file(path=marker), document


def _validate_reviewed_assignments(
    *, prepared_dir: Path, reviewed_labels: Path, profile: str | Path
) -> dict[str, Any]:
    """Validate label syntax, complete coverage and minimum class evidence."""

    protein_ids = frozenset(
        record.protein_id for record in iter_protein_fasta(path=prepared_dir / "proteins.faa")
    )
    loaded_profile = load_profile(source=profile)
    assignments = read_label_assignments(
        path=reviewed_labels,
        protein_ids=protein_ids,
        label_ids=loaded_profile.label_ids(),
    )
    validate_assignment_profile_compatibility(
        assignments=assignments,
        profile=loaded_profile,
    )
    covered = frozenset(assignment.protein_id for assignment in assignments)
    if covered != protein_ids:
        raise InputValidationError(
            "Reviewed labels must retain an explicit row for every prepared protein: "
            f"missing={len(protein_ids - covered)}"
        )
    positives = tuple(assignment for assignment in assignments if assignment.is_eligible_positive)
    comparison_by_target = {
        comparison.target_label_ids[0]: comparison.background_label_ids[0]
        for comparison in default_profile_comparisons(profile=loaded_profile)
    }
    control_label_ids = frozenset(
        label.label_id
        for label in loaded_profile.labels
        if label.label_id.startswith("control:") or label.system_class.upper() == "CONTROL"
    )
    positive_label_counts = Counter(assignment.label_id for assignment in positives)
    positive_evidence_status_counts = Counter(
        assignment.evidence_status for assignment in positives
    )
    control_positive_count = sum(
        count for label_id, count in positive_label_counts.items() if label_id in control_label_ids
    )
    target_positive_count = sum(
        count
        for label_id, count in positive_label_counts.items()
        if label_id in comparison_by_target
    )
    if target_positive_count < 1 or control_positive_count < 1:
        raise InputValidationError(
            "Approved labels require at least one reviewed-positive target and one "
            "reviewed-positive control assignment."
        )
    missing_matched_backgrounds = {
        target_label: background_label
        for target_label, background_label in comparison_by_target.items()
        if positive_label_counts[target_label] > 0 and positive_label_counts[background_label] == 0
    }
    if missing_matched_backgrounds:
        raise InputValidationError(
            "Each reviewed-positive analysis target requires its profile-resolved matched "
            f"background label: {missing_matched_backgrounds}"
        )
    status_counts = Counter(assignment.curation_status.value for assignment in assignments)
    evidence_supported_count = status_counts[EVIDENCE_POSITIVE_STATUS]
    return {
        "prepared_protein_count": len(protein_ids),
        "reviewed_assignment_count": len(assignments),
        "covered_protein_count": len(covered),
        "reviewed_positive_count": len(positives),
        "target_reviewed_positive_count": target_positive_count,
        "control_reviewed_positive_count": control_positive_count,
        "curation_status_counts": dict(sorted(status_counts.items())),
        "reviewed_positive_label_counts": dict(sorted(positive_label_counts.items())),
        "reviewed_positive_evidence_status_counts": dict(
            sorted(positive_evidence_status_counts.items())
        ),
        "synthetic_test_positive_count": positive_evidence_status_counts[
            AUTOMATED_TEST_EVIDENCE_STATUS
        ],
        "evidence_supported_positive_count": evidence_supported_count,
        "evidence_supported_target_count": positive_evidence_status_counts[EVIDENCE_TARGET_STATUS],
        "evidence_supported_control_count": positive_evidence_status_counts[
            EVIDENCE_CONTROL_STATUS
        ],
    }


def _validate_existing_e3_campaign(
    *,
    config_path: Path,
    campaign_id: str,
    profile: str,
    prepared_dir: Path,
    reviewed_labels: Path,
    domains: Path,
    evidence_inputs: Mapping[str, Path | None],
    structural_resource: Path,
    orthofinder_results: Path,
    foldseek_maximum_hits: int,
) -> None:
    """Require an existing campaign to retain all orchestration authorities."""

    config = load_config(path=config_path)
    if not config.foldseek.enabled:
        raise InputValidationError("Existing E3 campaign must keep Foldseek enabled.")
    if config.alphafold.enabled:
        raise InputValidationError(
            "Existing E3 campaign must reuse prepared AlphaFold coordinates rather than "
            "downloading models."
        )
    if config.inputs.orthofinder_resource is not None:
        raise InputValidationError(
            "Existing E3 campaign must use the checksum-verified raw OrthoFinder results."
        )
    expected = {
        "campaign_id": (config.campaign_id, campaign_id),
        "profile": (config.profile_name, profile),
        "sequences_fasta": (config.inputs.sequences_fasta, prepared_dir / "proteins.faa"),
        "label_assignments": (config.inputs.label_assignments, reviewed_labels),
        "domains": (config.inputs.domains, domains),
        "structures": (config.inputs.structures, prepared_dir / "structures.tsv"),
        "structural_alignment_resource": (
            config.inputs.structural_alignment_resource,
            structural_resource.resolve(),
        ),
        "orthofinder_results": (
            config.inputs.orthofinder_results,
            orthofinder_results.resolve(),
        ),
        "orthofinder_group_type": (config.inputs.orthofinder_group_type, "HOG"),
        "orthofinder_hierarchy_node": (config.inputs.orthofinder_hierarchy_node, "N0"),
        "foldseek.maximum_hits": (config.foldseek.maximum_hits, foldseek_maximum_hits),
        **{
            field: (getattr(config.inputs, field), expected)
            for field, expected in evidence_inputs.items()
        },
    }
    for field, (observed, wanted) in expected.items():
        if observed != wanted:
            raise InputValidationError(
                f"Existing E3 campaign changed orchestration field {field!r}: "
                f"{observed!r} != {wanted!r}"
            )


def _normalise_profile_source(*, profile: str | Path) -> str:
    """Normalise a built-in name or an existing custom-profile path."""

    source = str(profile)
    if "/" not in source and "\\" not in source and not source.endswith((".yaml", ".yml")):
        return validate_text(value=source, field_name="profile")
    candidate = Path(source).expanduser().resolve()
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"Custom profile is missing or empty: {candidate}")
    return str(candidate)


def _copy_file_atomic(*, source: Path, destination: Path) -> None:
    """Copy one file through a same-directory temporary path without overwrite."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        shutil.copyfile(source, temporary)
        os.chmod(temporary, 0o600)
        with temporary.open(mode="rb") as handle:
            os.fsync(handle.fileno())
        try:
            os.link(temporary, destination)
        except FileExistsError as error:
            raise PublicationError(
                f"Reviewed-label destination appeared concurrently: {destination}"
            ) from error
        temporary.unlink()
    except (OSError, PublicationError) as error:
        temporary.unlink(missing_ok=True)
        if isinstance(error, PublicationError):
            raise
        raise PublicationError(f"Could not stage reviewed labels {destination}: {error}") from error


def _require_nonnegative_integer(*, value: Any, field_name: str) -> int:
    """Return a non-negative integer or fail with context."""

    if isinstance(value, bool):
        raise InputValidationError(f"{field_name} must be a non-negative integer.")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise InputValidationError(f"{field_name} must be a non-negative integer.") from error
    if parsed < 0 or str(parsed) != str(value):
        raise InputValidationError(f"{field_name} must be a non-negative integer.")
    return parsed


def _require_sha256(*, value: Any, field_name: str) -> str:
    """Return a canonical SHA-256 digest or fail with context."""

    digest = str(value)
    if not _SHA256_PATTERN.fullmatch(digest):
        raise InputValidationError(f"{field_name} must be a lower-case SHA-256 digest.")
    return digest
