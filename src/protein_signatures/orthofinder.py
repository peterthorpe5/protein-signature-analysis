"""Read-only adapters for completed OrthoFinder 2.5.5 and 3 results."""

from __future__ import annotations

import csv
import logging
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .checksums import sha256_file
from .errors import InputValidationError
from .io_utils import iter_tsv, open_text, read_json
from .models import GroupMembership, OrthoFinderLayout
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)
_VERSION_PATTERN = re.compile(r"OrthoFinder\s+(?:version\s+)?([0-9]+(?:\.[0-9A-Za-z]+)+)", re.I)
_COMPLETION_MARKER = "OrthoFinder run completed"
_DIGEST_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_RAW_SOURCE_MODE = "RAW_COMPLETED_RESULTS"
_WORKFLOW_STAGE_SOURCE_MODE = "CHECKSUMMED_WORKFLOW_STAGE"
_WORKFLOW_AUTHORITY_FIELDS = (
    "mode",
    "archive_path",
    "archive_size_bytes",
    "archive_sha256",
    "published_results",
    "orthofinder_version",
    "decision_basis",
)
_WORKFLOW_VALIDATION_FIELDS = ("relative_path", "size_bytes", "sha256", "status")
_WORKFLOW_REQUIRED_RESULTS = (
    "WorkingDirectory/SpeciesIDs.txt",
    "WorkingDirectory/SequenceIDs.txt",
    "Orthogroups/Orthogroups.tsv",
    "Phylogenetic_Hierarchical_Orthogroups/N0.tsv",
    "Species_Tree/SpeciesTree_rooted_node_labels.txt",
)


def discover_orthofinder_layout(*, results_dir: Path) -> OrthoFinderLayout:
    """Discover authorities in completed OrthoFinder output without running it.

    Args:
        results_dir: Completed OrthoFinder result directory.

    Returns:
        Version-aware, read-only layout description.

    Raises:
        InputValidationError: If the output is incomplete or unsupported.
    """

    root = Path(results_dir).expanduser().resolve()
    if not root.is_dir():
        raise InputValidationError(f"OrthoFinder results directory does not exist: {root}")
    orthogroup_candidates = sorted(root.rglob("Orthogroups.tsv"))
    hog_paths = tuple(
        sorted(
            (
                path
                for path in root.rglob("*.tsv")
                if "Hierarchical_Orthogroups" in str(path.parent)
                and re.fullmatch(r"N[^/]*\.tsv", path.name)
            ),
            key=lambda item: item.name,
        )
    )
    orthogroups_path = orthogroup_candidates[0] if orthogroup_candidates else None
    if orthogroups_path is None and not hog_paths:
        raise InputValidationError(
            f"No Orthogroups.tsv or hierarchical orthogroup tables found beneath {root}."
        )
    species_ids = _first_or_none(paths=sorted(root.rglob("SpeciesIDs.txt")))
    sequence_ids = _first_or_none(paths=sorted(root.rglob("SequenceIDs.txt")))
    log_candidates = sorted(root.rglob("Log.txt"), key=lambda item: (len(item.parts), str(item)))
    workflow_authorities = _workflow_authority_paths(results_dir=root)
    if log_candidates:
        log_path: Path | None = log_candidates[0]
        version = _read_version(log_path=log_path)
        _require_completion_marker(log_path=log_path)
        source_mode = _RAW_SOURCE_MODE
        completion_authorities: tuple[Path, ...] = ()
    elif any(path.exists() for path in workflow_authorities):
        log_path = None
        version = _validate_workflow_stage_completion(
            results_dir=root,
            discovered_paths=tuple(
                path
                for path in (
                    orthogroups_path,
                    *hog_paths,
                    species_ids,
                    sequence_ids,
                )
                if path is not None
            ),
        )
        source_mode = _WORKFLOW_STAGE_SOURCE_MODE
        completion_authorities = workflow_authorities
    else:
        raise InputValidationError(f"Completed OrthoFinder output lacks Log.txt: {root}")
    major = _validate_supported_version(version=version)
    primary = "HIERARCHICAL_ORTHOGROUP" if hog_paths else "ORTHOGROUP"
    layout = OrthoFinderLayout(
        results_dir=root,
        version=version,
        major_version=major,
        adapter_name="orthofinder_v3" if major == 3 else "orthofinder_v2_5_5",
        primary_group_authority=primary,
        source_mode=source_mode,
        log_path=log_path,
        completion_authority_paths=completion_authorities,
        orthogroups_path=orthogroups_path,
        hog_paths=hog_paths,
        species_ids_path=species_ids,
        sequence_ids_path=sequence_ids,
    )
    LOGGER.info(
        "Detected completed OrthoFinder %s output in %s mode with %d HOG tables",
        version,
        source_mode,
        len(hog_paths),
    )
    return layout


def _workflow_authority_paths(*, results_dir: Path) -> tuple[Path, ...]:
    """Return the ordered completion authorities beside workflow Stage 04 results.

    Args:
        results_dir: Published ``04_orthofinder/Results`` directory.

    Returns:
        Stage manifest, archive authority and validation table paths.
    """

    stage_root = results_dir.parent
    return (
        stage_root / "stage_manifest.json",
        stage_root / "orthofinder_authority.tsv",
        stage_root / "orthofinder_reuse_validation.tsv",
    )


def _validate_workflow_stage_completion(
    *, results_dir: Path, discovered_paths: Sequence[Path]
) -> str:
    """Validate the checksum-bound completion contract for a reused Stage 04.

    This is the only accepted alternative to OrthoFinder's own completed
    ``Log.txt``.  All three workflow authorities must exist, the stage must be
    complete, and every required or consumed result must match both the outer
    manifest and the reuse-validation table where applicable.

    Args:
        results_dir: Published ``Results`` directory.
        discovered_paths: Group and identifier files selected by the adapter.

    Returns:
        Validated OrthoFinder version declared by the workflow authority.

    Raises:
        InputValidationError: If any completion authority is incomplete or inconsistent.
    """

    stage_root = results_dir.parent.resolve()
    manifest_path, authority_path, validation_path = _workflow_authority_paths(
        results_dir=results_dir
    )
    for path in (manifest_path, authority_path, validation_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise InputValidationError(
                f"Log-less OrthoFinder workflow stage lacks a required completion authority: {path}"
            )
    manifest = read_json(path=manifest_path)
    if not isinstance(manifest, Mapping):
        raise InputValidationError(f"OrthoFinder stage manifest must be an object: {manifest_path}")
    if manifest.get("status") != "complete":
        raise InputValidationError(
            f"OrthoFinder stage manifest is not marked complete: {manifest_path}"
        )
    configuration_digest = manifest.get("configuration_digest")
    if not isinstance(configuration_digest, str) or not _DIGEST_PATTERN.fullmatch(
        configuration_digest
    ):
        raise InputValidationError(
            f"OrthoFinder stage manifest has an invalid configuration_digest: {manifest_path}"
        )
    output_records = _stage_output_records(stage_root=stage_root, manifest=manifest)
    required_paths = dict.fromkeys(
        (
            authority_path,
            validation_path,
            *discovered_paths,
            *(results_dir / relative for relative in _WORKFLOW_REQUIRED_RESULTS),
        )
    )
    for required_path in required_paths:
        _verify_stage_output(
            stage_root=stage_root,
            path=required_path,
            output_records=output_records,
        )
    version = _read_workflow_authority(path=authority_path)
    _validate_workflow_reuse_rows(
        path=validation_path,
        results_dir=results_dir,
        stage_output_records=output_records,
    )
    return version


def _stage_output_records(
    *, stage_root: Path, manifest: Mapping[str, Any]
) -> dict[str, tuple[int, str]]:
    """Index a workflow stage output inventory after strict shape validation.

    Args:
        stage_root: Directory containing the stage manifest and outputs.
        manifest: Decoded workflow stage manifest.

    Returns:
        Mapping from safe relative path to declared size and checksum.

    Raises:
        InputValidationError: If the inventory is empty, malformed or ambiguous.
    """

    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise InputValidationError("OrthoFinder stage manifest has no output inventory.")
    records: dict[str, tuple[int, str]] = {}
    for index, item in enumerate(outputs):
        if not isinstance(item, Mapping):
            raise InputValidationError(
                f"OrthoFinder stage output record {index} must be an object."
            )
        relative = _safe_workflow_relative_path(
            value=item.get("path"),
            root=stage_root,
            context=f"stage output record {index}",
        )
        if relative in records:
            raise InputValidationError(
                f"OrthoFinder stage output inventory repeats path: {relative}"
            )
        size = _manifest_nonnegative_integer(
            value=item.get("size_bytes"),
            context=f"size_bytes in stage output record {index}",
        )
        digest = item.get("sha256")
        if not isinstance(digest, str) or not _DIGEST_PATTERN.fullmatch(digest):
            raise InputValidationError(
                f"Invalid sha256 in OrthoFinder stage output record {index}."
            )
        records[relative] = (size, digest)
    return records


def _verify_stage_output(
    *,
    stage_root: Path,
    path: Path,
    output_records: Mapping[str, tuple[int, str]],
) -> None:
    """Require one file to match its workflow stage output declaration.

    Args:
        stage_root: Directory against which output paths are relative.
        path: Required output file.
        output_records: Validated stage output declarations.

    Raises:
        InputValidationError: If the file is absent, undeclared or altered.
    """

    candidate = path.resolve()
    try:
        relative = candidate.relative_to(stage_root.resolve()).as_posix()
    except ValueError as error:
        raise InputValidationError(
            f"OrthoFinder stage authority escapes its stage root: {path}"
        ) from error
    declared = output_records.get(relative)
    if declared is None:
        raise InputValidationError(
            f"OrthoFinder stage manifest does not declare required output: {relative}"
        )
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"OrthoFinder stage output is missing or empty: {candidate}")
    expected_size, expected_digest = declared
    if candidate.stat().st_size != expected_size:
        raise InputValidationError(
            f"OrthoFinder stage output size differs from its manifest: {relative}"
        )
    if sha256_file(path=candidate) != expected_digest:
        raise InputValidationError(
            f"OrthoFinder stage output checksum differs from its manifest: {relative}"
        )


def _read_workflow_authority(*, path: Path) -> str:
    """Validate the reviewed-archive authority and return its version.

    Args:
        path: ``orthofinder_authority.tsv`` path already bound to the stage manifest.

    Returns:
        Declared supported OrthoFinder version.

    Raises:
        InputValidationError: If the authority is malformed or unsupported.
    """

    rows = tuple(iter_tsv(path=path, required_fields=_WORKFLOW_AUTHORITY_FIELDS))
    if len(rows) != 1:
        raise InputValidationError(
            f"OrthoFinder workflow authority must contain exactly one row: {path}"
        )
    row = rows[0]
    if row["mode"] != "reused_reviewed_archive":
        raise InputValidationError(
            f"Unsupported OrthoFinder workflow authority mode: {row['mode']!r}"
        )
    if row["published_results"] != "Results":
        raise InputValidationError(
            "OrthoFinder workflow authority must publish the adjacent Results directory."
        )
    if not row["archive_path"].strip() or not row["decision_basis"].strip():
        raise InputValidationError(
            "OrthoFinder workflow authority lacks archive_path or decision_basis."
        )
    archive_size = _tsv_nonnegative_integer(
        value=row["archive_size_bytes"], context="archive_size_bytes"
    )
    if archive_size == 0:
        raise InputValidationError(
            "OrthoFinder workflow authority archive_size_bytes must be positive."
        )
    if not _DIGEST_PATTERN.fullmatch(row["archive_sha256"]):
        raise InputValidationError("OrthoFinder workflow authority has an invalid archive_sha256.")
    version = row["orthofinder_version"]
    _validate_supported_version(version=version)
    return version


def _validate_workflow_reuse_rows(
    *,
    path: Path,
    results_dir: Path,
    stage_output_records: Mapping[str, tuple[int, str]],
) -> None:
    """Validate the exact five-file reviewed-archive extraction record.

    Args:
        path: ``orthofinder_reuse_validation.tsv`` path.
        results_dir: Published workflow ``Results`` directory.
        stage_output_records: Validated outer stage declarations.

    Raises:
        InputValidationError: If rows are missing, duplicated or inconsistent.
    """

    rows = tuple(iter_tsv(path=path, required_fields=_WORKFLOW_VALIDATION_FIELDS))
    observed: set[str] = set()
    for index, row in enumerate(rows, start=2):
        relative = _safe_workflow_relative_path(
            value=row["relative_path"],
            root=results_dir,
            context=f"reuse validation row {index}",
        )
        if relative in observed:
            raise InputValidationError(f"OrthoFinder reuse validation repeats path: {relative}")
        observed.add(relative)
        if row["status"] != "VALID":
            raise InputValidationError(f"OrthoFinder reuse validation is not VALID for {relative}.")
        declared_size = _tsv_nonnegative_integer(
            value=row["size_bytes"], context=f"size_bytes for {relative}"
        )
        declared_digest = row["sha256"]
        if not _DIGEST_PATTERN.fullmatch(declared_digest):
            raise InputValidationError(
                f"OrthoFinder reuse validation has an invalid sha256 for {relative}."
            )
        stage_record = stage_output_records.get(f"Results/{relative}")
        if stage_record != (declared_size, declared_digest):
            raise InputValidationError(
                f"OrthoFinder reuse validation disagrees with the stage manifest for {relative}."
            )
    expected = set(_WORKFLOW_REQUIRED_RESULTS)
    if observed != expected:
        missing = sorted(expected - observed)
        unexpected = sorted(observed - expected)
        raise InputValidationError(
            "OrthoFinder reuse validation does not contain the exact required file set; "
            f"missing={missing}, unexpected={unexpected}."
        )


def _safe_workflow_relative_path(*, value: Any, root: Path, context: str) -> str:
    """Validate and normalise a manifest or TSV path beneath one root.

    Args:
        value: Candidate POSIX-style relative path.
        root: Directory the value must remain beneath after resolution.
        context: Human-readable source location for errors.

    Returns:
        Normalised POSIX relative path.

    Raises:
        InputValidationError: If the path is missing, non-canonical or unsafe.
    """

    if not isinstance(value, str) or not value or value != value.strip() or "\\" in value:
        raise InputValidationError(f"Invalid relative path in OrthoFinder {context}: {value!r}")
    relative = Path(value)
    if relative.is_absolute() or value in {".", ".."} or ".." in relative.parts:
        raise InputValidationError(f"Unsafe relative path in OrthoFinder {context}: {value!r}")
    normalised = relative.as_posix()
    if normalised != value:
        raise InputValidationError(
            f"Non-canonical relative path in OrthoFinder {context}: {value!r}"
        )
    base = root.resolve()
    candidate = (base / relative).resolve()
    if not candidate.is_relative_to(base):
        raise InputValidationError(f"Unsafe relative path in OrthoFinder {context}: {value!r}")
    return normalised


def _manifest_nonnegative_integer(*, value: Any, context: str) -> int:
    """Validate a JSON non-negative integer without accepting booleans.

    Args:
        value: Candidate JSON value.
        context: Human-readable field location.

    Returns:
        Validated integer.

    Raises:
        InputValidationError: If the value is not a non-negative integer.
    """

    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InputValidationError(
            f"Invalid non-negative integer in OrthoFinder {context}: {value!r}"
        )
    return value


def _tsv_nonnegative_integer(*, value: str, context: str) -> int:
    """Parse one canonical non-negative integer from a workflow TSV.

    Args:
        value: Candidate decimal text.
        context: Human-readable field location.

    Returns:
        Parsed integer.

    Raises:
        InputValidationError: If the text is not canonical non-negative decimal.
    """

    if not re.fullmatch(r"0|[1-9][0-9]*", value):
        raise InputValidationError(
            f"Invalid non-negative integer in OrthoFinder {context}: {value!r}"
        )
    return int(value)


def read_group_memberships(
    *,
    layout: OrthoFinderLayout,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    protein_ids: frozenset[str],
) -> tuple[GroupMembership, ...]:
    """Read memberships from the selected HOG node or legacy orthogroups.

    Args:
        layout: Discovered completed-output layout.
        run_id: Caller-defined immutable identity for this OrthoFinder run.
        group_type: ``HOG`` or ``LEGACY_ORTHOGROUP``.
        hierarchy_node: HOG node such as ``N0``; ignored when only legacy groups exist.
        protein_ids: Campaign protein identifiers to retain.

    Returns:
        Unique memberships for campaign proteins.

    Raises:
        InputValidationError: If the table is malformed or mappings are ambiguous.
    """

    authority = group_type.strip().upper()
    if authority not in {"HOG", "LEGACY_ORTHOGROUP"}:
        raise InputValidationError(f"Unsupported OrthoFinder group type: {group_type!r}.")
    node = (
        validate_identifier(value=hierarchy_node, field_name="hierarchy node")
        if authority == "HOG"
        else ""
    )
    if authority == "LEGACY_ORTHOGROUP" and hierarchy_node:
        raise InputValidationError("LEGACY_ORTHOGROUP membership requires an empty hierarchy node.")
    source_run_id = validate_identifier(value=run_id, field_name="OrthoFinder run_id")
    internal_to_original = _read_sequence_id_map(path=layout.sequence_ids_path)
    hog_path = next((path for path in layout.hog_paths if path.stem == node), None)
    if authority == "HOG" and hog_path is None:
        available = [path.stem for path in layout.hog_paths]
        raise InputValidationError(
            f"Requested OrthoFinder hierarchy node {node!r} is unavailable; found {available}."
        )
    if authority == "HOG" and hog_path is not None:
        memberships = _read_hog_table(
            path=hog_path,
            run_id=source_run_id,
            node=node,
            protein_ids=protein_ids,
            identifier_map=internal_to_original,
        )
    elif authority == "LEGACY_ORTHOGROUP" and layout.orthogroups_path is not None:
        memberships = _read_orthogroup_table(
            path=layout.orthogroups_path,
            run_id=source_run_id,
            protein_ids=protein_ids,
            identifier_map=internal_to_original,
        )
    else:
        raise InputValidationError("Discovered OrthoFinder layout has no readable group authority.")
    by_protein: dict[str, str] = {}
    for membership in memberships:
        previous = by_protein.get(membership.protein_id)
        if previous is not None and previous != membership.group_id:
            raise InputValidationError(
                f"Protein {membership.protein_id!r} maps to multiple selected groups: "
                f"{previous!r} and {membership.group_id!r}."
            )
        by_protein[membership.protein_id] = membership.group_id
    if protein_ids and not memberships:
        raise InputValidationError(
            "No campaign FASTA identifiers mapped to the selected OrthoFinder membership table."
        )
    LOGGER.info(
        "Mapped %d of %d campaign proteins to completed OrthoFinder groups",
        len(by_protein),
        len(protein_ids),
    )
    return memberships


def _read_version(*, log_path: Path) -> str:
    """Extract an OrthoFinder version from its log.

    Args:
        log_path: Completed-run log.

    Returns:
        Version text.
    """

    with open_text(path=log_path) as handle:
        for line in handle:
            match = _VERSION_PATTERN.search(line)
            if match:
                return match.group(1)
    raise InputValidationError(f"Could not identify OrthoFinder version in {log_path}")


def _require_completion_marker(*, log_path: Path) -> None:
    """Require the completion record written by OrthoFinder itself.

    OrthoFinder's official entry point writes ``OrthoFinder run completed`` to
    ``Log.txt`` only after result production and citation output finish. Merely
    finding result-shaped tables is therefore not accepted as evidence of a
    completed upstream run.

    Args:
        log_path: Candidate OrthoFinder run log.

    Raises:
        InputValidationError: If the authoritative completion record is absent.
    """

    expected = _COMPLETION_MARKER.casefold()
    with open_text(path=log_path) as handle:
        if any(line.strip().casefold() == expected for line in handle):
            return
    raise InputValidationError(
        f"OrthoFinder Log.txt lacks the official completion marker "
        f"{_COMPLETION_MARKER!r}: {log_path}"
    )


def _validate_supported_version(*, version: str) -> int:
    """Validate the deliberately narrow raw-results compatibility contract.

    Args:
        version: Version extracted from an OrthoFinder log.

    Returns:
        Supported OrthoFinder major version.

    Raises:
        InputValidationError: If the version is neither 2.5.5 nor a 3.x release.
    """

    if version == "2.5.5":
        return 2
    if version.startswith("3."):
        return 3
    raise InputValidationError(
        f"Unsupported OrthoFinder version {version!r}; supported raw results are "
        "exactly OrthoFinder 2.5.5 or OrthoFinder 3.x."
    )


def _read_hog_table(
    *,
    path: Path,
    run_id: str,
    node: str,
    protein_ids: frozenset[str],
    identifier_map: dict[str, str],
) -> tuple[GroupMembership, ...]:
    """Parse a hierarchical orthogroup table.

    Args:
        path: HOG TSV.
        run_id: Immutable caller-supplied run identity.
        node: HOG hierarchy node.
        protein_ids: Campaign identifiers to retain.
        identifier_map: Optional OrthoFinder internal-to-original identifier map.

    Returns:
        Ordered campaign memberships.
    """

    with open_text(path=path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        hog_field = _find_field(fields=fields, candidates=("HOG", "HOG ID", "HOG_ID"))
        og_field = _find_field(
            fields=fields, candidates=("OG", "Orthogroup", "OG ID"), required=False
        )
        clade_field = _find_field(
            fields=fields,
            candidates=("Gene Tree Parent Clade", "Gene_Tree_Parent_Clade"),
            required=False,
        )
        metadata_fields = {item for item in (hog_field, og_field, clade_field) if item}
        species_fields = tuple(field for field in fields if field not in metadata_fields)
        if not species_fields:
            raise InputValidationError(f"HOG table contains no species columns: {path}")
        rows: list[GroupMembership] = []
        for row_number, row in enumerate(reader, start=2):
            group_id = validate_identifier(
                value=row.get(hog_field, ""), field_name=f"HOG at row {row_number}"
            )
            legacy = row.get(og_field, "").strip() if og_field else ""
            clade = row.get(clade_field, "").strip() if clade_field else ""
            for species in species_fields:
                for identifier in _split_members(value=row.get(species, "")):
                    original = identifier_map.get(identifier, identifier)
                    if original in protein_ids:
                        rows.append(
                            GroupMembership(
                                run_id=run_id,
                                group_type="HOG",
                                hierarchy_node=node,
                                group_id=group_id,
                                legacy_orthogroup_id=legacy,
                                gene_tree_parent_clade=clade,
                                species_label=species,
                                protein_id=original,
                            )
                        )
    return tuple(
        sorted(rows, key=lambda item: (item.group_id, item.species_label, item.protein_id))
    )


def _read_orthogroup_table(
    *,
    path: Path,
    run_id: str,
    protein_ids: frozenset[str],
    identifier_map: dict[str, str],
) -> tuple[GroupMembership, ...]:
    """Parse a legacy Orthogroups.tsv membership table.

    Args:
        path: Legacy orthogroup TSV.
        run_id: Immutable caller-supplied run identity.
        protein_ids: Campaign identifiers to retain.
        identifier_map: Optional internal-to-original identifier map.

    Returns:
        Ordered campaign memberships.
    """

    with open_text(path=path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        group_field = _find_field(fields=fields, candidates=("Orthogroup", "OG"))
        species_fields = tuple(field for field in fields if field != group_field)
        if not species_fields:
            raise InputValidationError(f"Orthogroups table contains no species columns: {path}")
        rows: list[GroupMembership] = []
        for row_number, row in enumerate(reader, start=2):
            group_id = validate_identifier(
                value=row.get(group_field, ""), field_name=f"orthogroup at row {row_number}"
            )
            for species in species_fields:
                for identifier in _split_members(value=row.get(species, "")):
                    original = identifier_map.get(identifier, identifier)
                    if original in protein_ids:
                        rows.append(
                            GroupMembership(
                                run_id=run_id,
                                group_type="LEGACY_ORTHOGROUP",
                                hierarchy_node="",
                                group_id=group_id,
                                legacy_orthogroup_id=group_id,
                                gene_tree_parent_clade="",
                                species_label=species,
                                protein_id=original,
                            )
                        )
    return tuple(
        sorted(rows, key=lambda item: (item.group_id, item.species_label, item.protein_id))
    )


def _read_sequence_id_map(*, path: Path | None) -> dict[str, str]:
    """Read optional OrthoFinder internal-to-original sequence identifiers.

    Args:
        path: SequenceIDs.txt path or ``None``.

    Returns:
        Mapping keyed by internal identifiers.
    """

    if path is None:
        return {}
    result: dict[str, str] = {}
    with open_text(path=path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            internal, separator, original = line.partition(": ")
            if not separator or not internal or not original:
                raise InputValidationError(
                    f"Malformed SequenceIDs.txt record at line {line_number}: {line!r}"
                )
            if internal in result:
                raise InputValidationError(f"Duplicate internal sequence identifier: {internal!r}")
            result[internal] = original.split(maxsplit=1)[0]
    return result


def _find_field(
    *, fields: tuple[str, ...], candidates: tuple[str, ...], required: bool = True
) -> str | None:
    """Find the first accepted field spelling.

    Args:
        fields: Actual table headings.
        candidates: Accepted alternatives in preference order.
        required: Whether absence is an error.

    Returns:
        Matching heading or ``None``.
    """

    by_normalised = {field.strip().casefold(): field for field in fields}
    for candidate in candidates:
        match = by_normalised.get(candidate.strip().casefold())
        if match is not None:
            return match
    if required:
        raise InputValidationError(
            f"Table headings {list(fields)} lack any required field from {list(candidates)}."
        )
    return None


def _split_members(*, value: str | None) -> tuple[str, ...]:
    """Split one OrthoFinder membership cell defensively.

    Args:
        value: Comma-delimited cell text.

    Returns:
        Non-empty member identifiers.
    """

    if value is None or not value.strip():
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _first_or_none(*, paths: list[Path]) -> Path | None:
    """Return the first ordered path or ``None``.

    Args:
        paths: Ordered paths.

    Returns:
        First path when present.
    """

    return paths[0] if paths else None
