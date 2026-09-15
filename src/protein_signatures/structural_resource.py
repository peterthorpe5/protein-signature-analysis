"""Adapter for completed US-align/TM-align structural evidence resources."""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import duckdb

from .checksums import sha256_file, sha256_json
from .errors import InputValidationError
from .feature_provenance import ALL_DATA_EXPLORATORY
from .io_utils import read_json
from .models import (
    FeatureRecord,
    PairwiseStructureComparison,
    SequenceRecord,
    StructuralAlignmentImport,
    StructureComparisonStatus,
    StructureCoverageScope,
)

LOGGER = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}")
_DATASET_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*")
_STANDALONE_MANIFEST = "standalone_outputs"
_AGGREGATE_MANIFEST = "end_to_end_datasets"
_REQUIRED_TABLES = (
    "structural_alignments.parquet",
    "pocket_comparisons.parquet",
    "structural_alignment_summary.parquet",
)


def import_structural_alignment_resource(
    *, resource_dir: Path, sequences: tuple[SequenceRecord, ...]
) -> StructuralAlignmentImport:
    """Import a checksum-valid structural-alignment result without recalculation.

    The adapter understands the published ``e3_structural_alignment`` table
    contract but returns generic structural comparisons and features.  Imported
    US-align/TM-align evidence and a new all-versus-all Foldseek screen therefore
    remain distinct and may coexist in one campaign.

    Args:
        resource_dir: Structural result root or end-to-end workflow root.
        sequences: Authoritative campaign sequences used to calculate coverage.

    Returns:
        Pairwise comparisons, pocket features, summaries and provenance inputs.

    Raises:
        InputValidationError: If the resource is incomplete or incompatible.
    """

    root = resolve_structural_resource_root(path=resource_dir)
    manifest_path = root / "provenance" / "run_manifest.json"
    manifest = read_json(path=manifest_path)
    if not isinstance(manifest, dict) or manifest.get("status") != "complete":
        raise InputValidationError(
            f"Structural resource manifest is not marked complete: {manifest_path}"
        )
    verify_structural_resource_outputs(resource_dir=root, manifest=manifest)
    _require_manifested_structural_tables(manifest=manifest)
    run_digest, package_version, manifest_paths = _structural_manifest_identity(
        resource_dir=root,
        manifest=manifest,
    )
    tables = root / "tables"
    required_paths = tuple(tables / name for name in _REQUIRED_TABLES)
    missing = [path.name for path in required_paths if not path.is_file()]
    if missing:
        raise InputValidationError(
            "Structural resource lacks required Parquet tables: " + "; ".join(missing)
        )
    lengths = {sequence.protein_id: sequence.sequence_length for sequence in sequences}
    summaries = _read_group_summaries(path=required_paths[2])
    expected_universe_sizes = {
        str(row["cluster_id"]): int(row["aligned_accession_count"]) for row in summaries
    }
    expected_reference_accessions = {
        str(row["cluster_id"]): str(row["reference_accession"]) for row in summaries
    }
    (
        comparisons,
        comparison_universe_members,
        reference_membership_row_count,
    ) = _read_global_alignments(
        path=required_paths[0],
        lengths=lengths,
        run_digest=run_digest,
        expected_universe_sizes=expected_universe_sizes,
        expected_reference_accessions=expected_reference_accessions,
    )
    features = _read_pocket_features(
        path=required_paths[1],
        protein_ids=frozenset(lengths),
        run_digest=run_digest,
    )
    LOGGER.info(
        "Imported structural resource comparisons=%d reference_membership_rows=%d "
        "features=%d groups=%d",
        len(comparisons),
        reference_membership_row_count,
        len(features),
        len(summaries),
    )
    return StructuralAlignmentImport(
        comparisons=comparisons,
        features=features,
        group_summaries=summaries,
        input_paths=(*manifest_paths, *required_paths),
        package_version=package_version,
        run_digest=run_digest,
        comparison_universe_members=comparison_universe_members,
        reference_membership_row_count=reference_membership_row_count,
    )


def resolve_structural_resource_root(*, path: Path) -> Path:
    """Resolve one explicit structural result beneath a supplied directory.

    Args:
        path: Structural result root or parent end-to-end run root.

    Returns:
        Directory containing ``provenance/run_manifest.json`` and ``tables``.

    Raises:
        InputValidationError: If no unique supported root can be resolved.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_dir():
        raise InputValidationError(f"Structural resource directory does not exist: {source}")
    candidates = [source]
    candidates.extend(
        (
            source / "09b_structural_alignment" / "structural_alignment",
            source / "structural_alignment",
        )
    )
    matches = tuple(
        candidate
        for candidate in candidates
        if (candidate / "provenance" / "run_manifest.json").is_file()
    )
    unique = tuple(dict.fromkeys(matches))
    if len(unique) != 1:
        raise InputValidationError(
            "Expected exactly one structural-alignment resource at the supplied path; "
            f"found {len(unique)}."
        )
    return unique[0]


def verify_structural_resource_outputs(
    *, resource_dir: Path, manifest: Mapping[str, Any] | None = None
) -> None:
    """Verify the scientific files declared by a structural resource manifest.

    Standalone component manifests checksum their files directly in ``outputs``.
    Aggregated Stage 09b manifests declare Parquet files in ``datasets``; those
    declarations are cross-checked against the checksum-bearing outer stage
    manifest before any file is accepted.

    Args:
        resource_dir: Resolved structural resource root.
        manifest: Optional already-decoded manifest.

    Raises:
        InputValidationError: If any declaration, size or checksum is invalid.
    """

    root = Path(resource_dir).expanduser().resolve()
    value = (
        read_json(path=root / "provenance" / "run_manifest.json") if manifest is None else manifest
    )
    if not isinstance(value, Mapping):
        raise InputValidationError("Structural run manifest must contain an object.")
    if value.get("status") != "complete":
        raise InputValidationError("Structural run manifest is not marked complete.")
    variant = _manifest_variant(manifest=value)
    if variant == _STANDALONE_MANIFEST:
        records = _output_records(
            manifest=value,
            missing_message="Structural run manifest has no output inventory.",
            record_name="Structural output record",
        )
        for relative, (index, row) in records.items():
            _verify_output_record(
                root=root,
                relative=relative,
                index=index,
                row=row,
                record_name="Structural output record",
            )
        LOGGER.info(
            "Verified standalone structural resource outputs=%d root=%s",
            len(records),
            root,
        )
        return
    _verify_aggregate_outputs(resource_dir=root, manifest=value)


def _require_manifested_structural_tables(*, manifest: Mapping[str, Any]) -> None:
    """Require every consumed structural table in the checksum inventory.

    Args:
        manifest: Already validated structural resource manifest.

    Raises:
        InputValidationError: If a required Parquet table is not declared exactly once
            beneath the resource ``tables`` directory.
    """

    variant = _manifest_variant(manifest=manifest)
    if variant == _STANDALONE_MANIFEST:
        outputs = manifest.get("outputs")
        if not isinstance(outputs, list):  # pragma: no cover - variant guard
            raise InputValidationError("Structural run manifest has no output inventory.")
        declared = {
            _safe_path(value=row.get("path"), index=index)
            for index, row in enumerate(outputs)
            if isinstance(row, Mapping)
        }
        required = {Path("tables") / name for name in _REQUIRED_TABLES}
        missing = sorted(str(path) for path in required - declared)
    else:
        datasets = manifest.get("datasets")
        if not isinstance(datasets, Mapping):  # pragma: no cover - variant guard
            raise InputValidationError("Structural aggregate has no dataset inventory.")
        required_names = {Path(name).stem for name in _REQUIRED_TABLES}
        missing = sorted(f"tables/{name}.parquet" for name in required_names - set(datasets))
    if missing:
        raise InputValidationError(
            "Required structural tables are not checksum-inventoried in the manifest: "
            + "; ".join(missing)
        )


def _manifest_variant(*, manifest: Mapping[str, Any]) -> str:
    """Identify one supported, unambiguous structural manifest contract.

    Args:
        manifest: Decoded structural run manifest.

    Returns:
        Internal identifier for the standalone or aggregate contract.

    Raises:
        InputValidationError: If inventories are absent, malformed or ambiguous.
    """

    has_outputs = "outputs" in manifest
    has_datasets = "datasets" in manifest
    if has_outputs and has_datasets:
        raise InputValidationError(
            "Structural run manifest mixes outputs and datasets inventories."
        )
    if has_outputs:
        outputs = manifest.get("outputs")
        if not isinstance(outputs, list) or not outputs:
            raise InputValidationError("Structural run manifest has no output inventory.")
        return _STANDALONE_MANIFEST
    if has_datasets:
        datasets = manifest.get("datasets")
        if not isinstance(datasets, Mapping) or not datasets:
            raise InputValidationError("Structural aggregate has no dataset inventory.")
        return _AGGREGATE_MANIFEST
    raise InputValidationError(
        "Structural run manifest has neither an outputs nor a datasets inventory."
    )


def _output_records(
    *,
    manifest: Mapping[str, Any],
    missing_message: str,
    record_name: str,
) -> dict[Path, tuple[int, Mapping[str, Any]]]:
    """Index a safe, duplicate-free manifest output inventory.

    Args:
        manifest: Manifest containing an ``outputs`` list.
        missing_message: Diagnostic used for a missing or empty list.
        record_name: Human-readable record name used in diagnostics.

    Returns:
        Output records keyed by safe relative path.

    Raises:
        InputValidationError: If the inventory shape or a path is invalid.
    """

    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise InputValidationError(missing_message)
    records: dict[Path, tuple[int, Mapping[str, Any]]] = {}
    for index, row in enumerate(outputs):
        if not isinstance(row, Mapping):
            raise InputValidationError(f"{record_name} {index} is not an object.")
        relative = _safe_path(
            value=row.get("path"),
            index=index,
            record_name=record_name,
        )
        if relative in records:
            raise InputValidationError(f"Structural output path is repeated: {relative}")
        records[relative] = (index, row)
    return records


def _verify_output_record(
    *,
    root: Path,
    relative: Path,
    index: int,
    row: Mapping[str, Any],
    record_name: str,
) -> str:
    """Verify one safe output record and return its declared SHA-256.

    Args:
        root: Directory against which the relative path is resolved.
        relative: Previously validated relative file path.
        index: Zero-based manifest record index.
        row: Output record containing size and checksum fields.
        record_name: Human-readable record name used in diagnostics.

    Returns:
        Verified lower-case SHA-256 digest.

    Raises:
        InputValidationError: If the file, size or checksum is invalid.
    """

    candidate = (root / relative).resolve()
    if root not in candidate.parents or not candidate.is_file():
        raise InputValidationError(f"Structural output is missing or unsafe: {relative}")
    size = _non_negative_integer(value=row.get("size_bytes"), field="size_bytes")
    if candidate.stat().st_size != size:
        raise InputValidationError(f"Structural output size mismatch: {relative}")
    expected = str(row.get("sha256", ""))
    if _DIGEST.fullmatch(expected) is None:
        raise InputValidationError(f"{record_name} {index} has an invalid SHA-256.")
    if sha256_file(path=candidate) != expected:
        raise InputValidationError(f"Structural output checksum mismatch: {relative}")
    return expected


def _aggregate_dataset_records(
    *, resource_dir: Path, manifest: Mapping[str, Any]
) -> dict[str, tuple[Path, str]]:
    """Resolve and validate aggregate dataset declarations.

    Absolute paths in the upstream manifest are treated as provenance text, not
    as read targets, so a checksum-valid copied run remains portable. The final
    ``tables/<dataset>.parquet`` suffix must still match exactly.

    Args:
        resource_dir: Resolved aggregate structural resource root.
        manifest: Decoded aggregate run manifest.

    Returns:
        Local dataset paths and inner-manifest checksums keyed by dataset name.

    Raises:
        InputValidationError: If a dataset name, path or checksum is invalid.
    """

    datasets = manifest.get("datasets")
    if not isinstance(datasets, Mapping) or not datasets:
        raise InputValidationError("Structural aggregate has no dataset inventory.")
    root = Path(resource_dir).expanduser().resolve()
    records: dict[str, tuple[Path, str]] = {}
    for name, row in sorted(datasets.items(), key=lambda item: str(item[0])):
        if not isinstance(name, str) or _DATASET_NAME.fullmatch(name) is None:
            raise InputValidationError(f"Structural aggregate dataset name is invalid: {name!r}")
        if not isinstance(row, Mapping):
            raise InputValidationError(f"Structural aggregate dataset {name!r} is not an object.")
        raw_path = row.get("path")
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise InputValidationError(f"Structural aggregate dataset {name!r} lacks a path.")
        recorded_path = Path(raw_path)
        expected_name = f"{name}.parquet"
        if len(recorded_path.parts) < 2 or recorded_path.parts[-2:] != (
            "tables",
            expected_name,
        ):
            raise InputValidationError(
                f"Structural aggregate dataset {name!r} has an incompatible path: {raw_path!r}"
            )
        expected_digest = str(row.get("sha256", ""))
        if _DIGEST.fullmatch(expected_digest) is None:
            raise InputValidationError(
                f"Structural aggregate dataset {name!r} has an invalid SHA-256."
            )
        records[name] = (root / "tables" / expected_name, expected_digest)
    return records


def _read_aggregate_stage_manifest(*, resource_dir: Path) -> tuple[Path, Mapping[str, Any]]:
    """Read the completed checksum authority for an aggregated Stage 09b result.

    Args:
        resource_dir: Resolved ``structural_alignment`` directory.

    Returns:
        Outer stage-manifest path and decoded object.

    Raises:
        InputValidationError: If the manifest is absent, malformed or incomplete.
    """

    path = Path(resource_dir).expanduser().resolve().parent / "stage_manifest.json"
    manifest = read_json(path=path)
    if not isinstance(manifest, Mapping) or manifest.get("status") != "complete":
        raise InputValidationError(
            f"Aggregated structural Stage 09b manifest is not marked complete: {path}"
        )
    return path, manifest


def _verify_aggregate_outputs(*, resource_dir: Path, manifest: Mapping[str, Any]) -> None:
    """Cross-check aggregate datasets against the outer Stage 09b inventory.

    Args:
        resource_dir: Resolved aggregate structural resource root.
        manifest: Decoded aggregate run manifest.

    Raises:
        InputValidationError: If either manifest or any scientific file differs.
    """

    root = Path(resource_dir).expanduser().resolve()
    stage_manifest_path, stage_manifest = _read_aggregate_stage_manifest(resource_dir=root)
    records = _output_records(
        manifest=stage_manifest,
        missing_message=f"Stage 09b manifest has no output inventory: {stage_manifest_path}",
        record_name="Stage 09b output record",
    )
    datasets = _aggregate_dataset_records(resource_dir=root, manifest=manifest)
    required_files: list[tuple[str, Path, str | None]] = [
        ("aggregate run manifest", root / "provenance" / "run_manifest.json", None)
    ]
    required_files.extend(
        (f"dataset {name!r}", path, digest) for name, (path, digest) in datasets.items()
    )
    stage_root = root.parent
    for label, path, inner_digest in required_files:
        try:
            relative = path.resolve().relative_to(stage_root)
        except ValueError as error:  # pragma: no cover - constructed local paths
            raise InputValidationError(
                f"Aggregated structural {label} is outside Stage 09b: {path}"
            ) from error
        record = records.get(relative)
        if record is None:
            raise InputValidationError(
                f"Aggregated structural {label} is not checksum-inventoried in "
                f"{stage_manifest_path}: {relative}"
            )
        index, row = record
        outer_digest = _verify_output_record(
            root=stage_root,
            relative=relative,
            index=index,
            row=row,
            record_name="Stage 09b output record",
        )
        if inner_digest is not None and inner_digest != outer_digest:
            raise InputValidationError(
                f"Structural aggregate checksum disagrees with Stage 09b for {label}."
            )
    LOGGER.info(
        "Verified aggregated Stage 09b structural datasets=%d authority=%s",
        len(datasets),
        stage_manifest_path,
    )


def _structural_manifest_identity(
    *, resource_dir: Path, manifest: Mapping[str, Any]
) -> tuple[str, str, tuple[Path, ...]]:
    """Return a stable run identity, package version and manifest authorities.

    Args:
        resource_dir: Resolved structural resource root.
        manifest: Verified structural run manifest.

    Returns:
        Run digest, producing package version and ordered manifest paths.

    Raises:
        InputValidationError: If the applicable upstream identity is invalid.
    """

    root = Path(resource_dir).expanduser().resolve()
    manifest_path = root / "provenance" / "run_manifest.json"
    variant = _manifest_variant(manifest=manifest)
    if variant == _STANDALONE_MANIFEST:
        run_digest = str(manifest.get("run_digest", ""))
        if _DIGEST.fullmatch(run_digest) is None:
            raise InputValidationError("Structural resource run_digest is not a SHA-256 value.")
        package_version = str(manifest.get("package_version") or "unknown").strip()
        return run_digest, package_version or "unknown", (manifest_path,)

    configuration_digest = str(manifest.get("configuration_digest", ""))
    if _DIGEST.fullmatch(configuration_digest) is None:
        raise InputValidationError(
            "Structural aggregate configuration_digest is not a SHA-256 value."
        )
    datasets = _aggregate_dataset_records(resource_dir=root, manifest=manifest)
    run_digest = sha256_json(
        value={
            "schema": "e3workflow_stage_09b_aggregate_v1",
            "configuration_digest": configuration_digest,
            "datasets": {name: digest for name, (_, digest) in datasets.items()},
        }
    )
    stage_manifest_path, stage_manifest = _read_aggregate_stage_manifest(resource_dir=root)
    package_version = str(stage_manifest.get("package_version") or "unknown").strip()
    return (
        run_digest,
        package_version or "unknown",
        (manifest_path, stage_manifest_path),
    )


def _read_global_alignments(
    *,
    path: Path,
    lengths: Mapping[str, int],
    run_digest: str,
    expected_universe_sizes: Mapping[str, int],
    expected_reference_accessions: Mapping[str, str],
) -> tuple[
    tuple[PairwiseStructureComparison, ...],
    dict[str, frozenset[str]],
    int,
]:
    """Convert published structural alignments to generic comparisons.

    Args:
        path: Published global-alignment Parquet table.
        lengths: Campaign protein lengths used as coverage denominators.
        run_digest: Verified upstream run identity.
        expected_universe_sizes: Aligned protein count keyed by upstream cluster.
        expected_reference_accessions: Declared reference protein keyed by upstream
            cluster.

    Returns:
        Validated comparisons, proved cluster-local assessment universes and the
        number of explicit reference-membership rows omitted from comparisons.
    """

    required = (
        "cluster_id",
        "reference_accession",
        "mobile_accession",
        "alignment_tool",
        "status",
        "tool_version",
        "aligned_length",
        "rmsd_angstrom",
        "minimum_tm_score",
    )
    rows = _read_parquet(path=path, required=required)
    comparisons: dict[tuple[str, str, str, str], PairwiseStructureComparison] = {}
    universe_members: dict[str, set[str]] = {}
    universe_by_cluster: dict[str, str] = {}
    reference_membership_rows: set[tuple[str, str, str]] = set()
    matching_row_count = 0
    for row_number, row in enumerate(rows, start=1):
        reference = str(row["reference_accession"] or "").strip()
        mobile = str(row["mobile_accession"] or "").strip()
        if reference not in lengths or mobile not in lengths:
            continue
        matching_row_count += 1
        tool = str(row["alignment_tool"] or "").strip()
        version = str(row["tool_version"] or "").strip()
        cluster_id = str(row["cluster_id"] or "").strip()
        if not tool or not version or not cluster_id:
            raise InputValidationError(
                f"Structural alignment row {row_number} lacks tool or cluster provenance."
            )
        raw_status = str(row["status"] or "").strip().upper()
        universe_id = "ES3A_" + sha256_json(value=(run_digest, cluster_id))[:24]
        previous_universe = universe_by_cluster.setdefault(cluster_id, universe_id)
        if previous_universe != universe_id:  # pragma: no cover - deterministic hash guard
            raise InputValidationError(
                f"Structural cluster {cluster_id!r} resolved to inconsistent universes."
            )
        expected_reference = expected_reference_accessions.get(cluster_id)
        if expected_reference is None:
            raise InputValidationError(
                f"Structural alignment cluster {cluster_id!r} lacks a group summary."
            )
        if reference != expected_reference:
            raise InputValidationError(
                f"Structural alignment row {row_number} names reference {reference!r}, "
                f"but cluster {cluster_id!r} declares {expected_reference!r}."
            )
        if reference == mobile:
            if raw_status != "REFERENCE":
                raise InputValidationError(
                    f"Structural alignment row {row_number} compares a protein with "
                    "itself without the explicit REFERENCE status."
                )
            aligned = _optional_integer(
                value=row["aligned_length"],
                field="aligned_length",
                row_number=row_number,
            )
            rmsd = _optional_number(
                value=row["rmsd_angstrom"], field="rmsd_angstrom", row_number=row_number
            )
            score = _optional_number(
                value=row["minimum_tm_score"],
                field="minimum_tm_score",
                row_number=row_number,
            )
            if aligned is not None or rmsd != 0.0 or score != 1.0:
                raise InputValidationError(
                    f"Structural reference row {row_number} does not contain the "
                    "expected aligned_length=NULL, RMSD=0 and minimum TM-score=1 sentinels."
                )
            reference_key = (cluster_id, tool, reference)
            if reference_key in reference_membership_rows:
                raise InputValidationError(
                    f"Duplicate imported structural reference row: {reference_key!r}"
                )
            reference_membership_rows.add(reference_key)
            universe_members.setdefault(universe_id, set()).add(reference)
            continue
        status = _parse_comparison_status(value=raw_status, row_number=row_number)
        aligned = _optional_integer(
            value=row["aligned_length"], field="aligned_length", row_number=row_number
        )
        rmsd = _optional_number(
            value=row["rmsd_angstrom"], field="rmsd_angstrom", row_number=row_number
        )
        score = _optional_number(
            value=row["minimum_tm_score"], field="minimum_tm_score", row_number=row_number
        )
        if score is not None and not 0.0 <= score <= 1.0:
            raise InputValidationError(
                f"minimum_tm_score outside zero to one at row {row_number}: {score}"
            )
        if rmsd is not None and rmsd < 0.0:
            raise InputValidationError(f"Negative RMSD at structural row {row_number}.")
        coverage_reference = min(1.0, aligned / lengths[reference]) if aligned is not None else None
        coverage_mobile = min(1.0, aligned / lengths[mobile]) if aligned is not None else None
        successful = status in {
            StructureComparisonStatus.COMPLETE,
            StructureComparisonStatus.PASS,
            StructureComparisonStatus.SUCCESS,
        }
        if successful and (score is None or coverage_reference is None or coverage_mobile is None):
            raise InputValidationError(
                f"Successful structural alignment at row {row_number} lacks TM score "
                "or bilateral coverage."
            )
        if not successful and any(value is not None for value in (aligned, rmsd, score)):
            raise InputValidationError(
                f"Unsuccessful structural alignment at row {row_number} carries metrics."
            )
        source_id = "|".join((cluster_id, tool, reference, mobile))
        universe_members.setdefault(universe_id, set()).update((reference, mobile))
        key = (cluster_id, tool, reference, mobile)
        if key in comparisons:
            raise InputValidationError(f"Duplicate imported structural comparison: {key!r}")
        comparisons[key] = PairwiseStructureComparison(
            protein_a_id=reference,
            protein_b_id=mobile,
            comparison_tool=tool,
            comparison_tool_version=version,
            tm_score=score,
            rmsd_angstrom=rmsd,
            aligned_residue_count=aligned,
            coverage_a=coverage_reference,
            coverage_b=coverage_mobile,
            comparison_status=status,
            source_record_id=source_id,
            comparison_universe_id=universe_id,
            coverage_scope=StructureCoverageScope.FULL_SEQUENCE,
        )
    if rows and matching_row_count == 0:
        raise InputValidationError(
            "No structural-alignment accessions matched the campaign FASTA identifiers."
        )
    for cluster_id, universe_id in sorted(universe_by_cluster.items()):
        expected_size = expected_universe_sizes.get(cluster_id)
        if expected_size is None:
            raise InputValidationError(
                f"Structural alignment cluster {cluster_id!r} lacks a group summary."
            )
        observed_size = len(universe_members[universe_id])
        if observed_size != expected_size:
            raise InputValidationError(
                f"Structural comparison universe {cluster_id!r} contains {observed_size} "
                f"campaign proteins but its summary declares {expected_size}; complete "
                "absence assessment is not possible."
            )
    frozen_members = {
        universe_id: frozenset(members) for universe_id, members in sorted(universe_members.items())
    }
    if reference_membership_rows:
        LOGGER.info(
            "Accepted %d explicit structural REFERENCE rows as assessment-universe "
            "membership and omitted them from pairwise comparisons",
            len(reference_membership_rows),
        )
    return (
        tuple(comparisons[key] for key in sorted(comparisons)),
        frozen_members,
        len(reference_membership_rows),
    )


def _read_pocket_features(
    *, path: Path, protein_ids: frozenset[str], run_digest: str
) -> tuple[FeatureRecord, ...]:
    """Derive cluster-namespaced exploratory pocket features.

    Args:
        path: Published pocket-comparison Parquet table.
        protein_ids: Authoritative campaign protein identifiers.
        run_digest: Verified predecessor run identity.

    Returns:
        Positive pocket features whose identifiers cannot collide across upstream
        structural clusters.
    """

    required = (
        "cluster_id",
        "reference_accession",
        "mobile_accession",
        "alignment_tool",
        "status",
        "same_pocket_position_supported",
        "pocket_structure_conserved",
    )
    features: dict[tuple[str, str, str], FeatureRecord] = {}
    for row in _read_parquet(path=path, required=required):
        reference = str(row["reference_accession"] or "").strip()
        mobile = str(row["mobile_accession"] or "").strip()
        if reference not in protein_ids or mobile not in protein_ids:
            continue
        status = str(row["status"] or "").strip().upper()
        if status not in {"ASSESSED", "COMPLETE", "PASS"}:
            continue
        cluster_id = str(row["cluster_id"] or "").strip()
        alignment_tool = str(row["alignment_tool"] or "").strip()
        if not cluster_id or not alignment_tool:
            raise InputValidationError(
                "Assessed pocket comparison lacks cluster or alignment-tool provenance."
            )
        feature_namespace = "ES3A_" + sha256_json(value=(run_digest, cluster_id))[:24]
        evidence_reference = "|".join(
            (
                cluster_id,
                alignment_tool,
                reference,
                mobile,
            )
        )
        for field, local_feature_id, feature_name in (
            (
                "same_pocket_position_supported",
                "SAME_3D_POCKET_POSITION",
                "Same three-dimensional pocket position supported",
            ),
            (
                "pocket_structure_conserved",
                "CONSERVED_3D_POCKET",
                "Conserved three-dimensional pocket supported",
            ),
        ):
            if not _boolean(value=row[field], field=field):
                continue
            feature_id = f"{feature_namespace}:{local_feature_id}"
            definition_digest = sha256_json(
                value={
                    "upstream_run_digest": run_digest,
                    "upstream_cluster_id": cluster_id,
                    "feature_type": "STRUCTURAL_POCKET",
                    "feature_id": feature_id,
                    "support_field": field,
                    "support_rule": "ANY_ASSESSED_SUPPORTED_PAIR",
                }
            )
            for protein_id in (reference, mobile):
                key = (protein_id, feature_id, evidence_reference)
                features[key] = FeatureRecord(
                    protein_id=protein_id,
                    feature_type="STRUCTURAL_POCKET",
                    feature_id=feature_id,
                    feature_name=f"{feature_name} in upstream cluster {cluster_id}",
                    start=None,
                    end=None,
                    evidence_status="IMPORTED_COMPLETE",
                    evidence_source="e3_structural_alignment",
                    evidence_reference=evidence_reference,
                    derivation_scope=ALL_DATA_EXPLORATORY,
                    feature_definition_sha256=definition_digest,
                    derivation_cohort_sha256="",
                )
    return tuple(features[key] for key in sorted(features))


def _read_group_summaries(*, path: Path) -> tuple[dict[str, Any], ...]:
    """Read a stable subset of structural group-summary fields."""

    fields = (
        "cluster_id",
        "primary_group_type",
        "primary_group_id",
        "reference_accession",
        "alignment_tools",
        "alignment_tool_count",
        "selected_accession_count",
        "model_available_accession_count",
        "aligned_accession_count",
        "supported_accession_count",
        "position_supported_accession_count",
        "group_support_fraction",
        "group_position_support_fraction",
        "mean_minimum_tm_score",
        "mean_pocket_overlap_fraction",
        "median_centroid_distance_angstrom",
        "position_alignment_status",
        "alignment_status",
        "interpretation",
    )
    rows = _read_parquet(path=path, required=fields)
    normalised: list[dict[str, Any]] = []
    identifiers: set[tuple[str, str, str]] = set()
    integer_fields = fields[5:11]
    number_fields = fields[11:16]
    for row_number, row in enumerate(rows, start=1):
        record = {field: row[field] for field in fields}
        key = (
            str(record["cluster_id"] or ""),
            str(record["primary_group_type"] or ""),
            str(record["primary_group_id"] or ""),
        )
        if not all(key) or key in identifiers:
            raise InputValidationError(
                f"Invalid or duplicate structural group identity at row {row_number}: {key!r}"
            )
        identifiers.add(key)
        for field in integer_fields:
            record[field] = _non_negative_integer(value=record[field], field=field)
        for field in number_fields:
            record[field] = _optional_number(
                value=record[field], field=field, row_number=row_number
            )
        for field in fields:
            if field not in {*integer_fields, *number_fields}:
                record[field] = str(record[field] or "")
        normalised.append(record)
    return tuple(
        sorted(
            normalised,
            key=lambda row: (
                str(row["primary_group_type"]),
                str(row["primary_group_id"]),
                str(row["cluster_id"]),
            ),
        )
    )


def _parse_comparison_status(*, value: Any, row_number: int) -> StructureComparisonStatus:
    """Parse an upstream pairwise-comparison completion state.

    Args:
        value: Raw upstream status.
        row_number: One-based row number for diagnostics.

    Returns:
        Controlled comparison status.

    Raises:
        InputValidationError: If the upstream status is unsupported.
    """

    normalised = str(value or "").strip().upper()
    try:
        return StructureComparisonStatus(normalised)
    except ValueError as error:
        allowed = ", ".join(item.value for item in StructureComparisonStatus)
        raise InputValidationError(
            f"Unsupported structural status at row {row_number}; expected one of "
            f"{allowed}, received {normalised!r}."
        ) from error


def _read_parquet(*, path: Path, required: tuple[str, ...]) -> list[dict[str, Any]]:
    """Read a Parquet table after validating its required unique columns."""

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty structural Parquet table: {source}")
    try:
        with duckdb.connect(":memory:") as connection:
            description = connection.execute(
                "DESCRIBE SELECT * FROM read_parquet(?)", [str(source)]
            ).fetchall()
            columns = tuple(str(row[0]) for row in description)
            missing = sorted(set(required) - set(columns))
            if missing:
                raise InputValidationError(
                    f"Structural table {source} lacks required fields: {missing}"
                )
            selected = ", ".join(f'"{field}"' for field in required)
            cursor = connection.execute(f"SELECT {selected} FROM read_parquet(?)", [str(source)])
            return [dict(zip(required, row, strict=True)) for row in cursor.fetchall()]
    except duckdb.Error as error:
        raise InputValidationError(
            f"Could not read structural Parquet table {source}: {error}"
        ) from error


def _optional_integer(*, value: Any, field: str, row_number: int) -> int | None:
    """Parse an optional positive integer from structural evidence."""

    if value is None or str(value).strip() == "":
        return None
    parsed = _non_negative_integer(value=value, field=field)
    if parsed == 0:
        raise InputValidationError(f"{field} must be positive at row {row_number}.")
    return parsed


def _non_negative_integer(*, value: Any, field: str) -> int:
    """Parse a non-negative non-Boolean integer with field context."""

    if isinstance(value, bool):
        raise InputValidationError(f"{field} must be an integer, not a Boolean.")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise InputValidationError(f"{field} must be an integer: {value!r}") from error
    if parsed < 0:
        raise InputValidationError(f"{field} must be non-negative: {parsed}")
    return parsed


def _optional_number(*, value: Any, field: str, row_number: int) -> float | None:
    """Parse an optional finite floating-point structural value."""

    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            f"{field} must be numeric at row {row_number}: {value!r}"
        ) from error
    if not math.isfinite(parsed):
        raise InputValidationError(f"{field} must be finite at row {row_number}.")
    return parsed


def _boolean(*, value: Any, field: str) -> bool:
    """Parse one explicit Boolean without treating missing evidence as false."""

    if isinstance(value, bool):
        return value
    text = str(value).strip().lower() if value is not None else ""
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise InputValidationError(f"{field} must contain an explicit Boolean: {value!r}")


def _safe_path(*, value: Any, index: int, record_name: str = "Structural output record") -> Path:
    """Return one safe non-empty manifest-relative path.

    Args:
        value: Raw manifest path.
        index: Zero-based record index.
        record_name: Human-readable record name used in diagnostics.

    Returns:
        Validated relative path.

    Raises:
        InputValidationError: If the path is empty, absolute or traverses parents.
    """

    if not isinstance(value, str) or not value:
        raise InputValidationError(f"{record_name} {index} lacks a path.")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path == Path("."):
        raise InputValidationError(f"{record_name} {index} has an unsafe path: {value!r}")
    return path
