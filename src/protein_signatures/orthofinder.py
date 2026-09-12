"""Read-only adapters for completed OrthoFinder 2.5.5 and 3 results."""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

from .errors import InputValidationError
from .io_utils import open_text
from .models import GroupMembership, OrthoFinderLayout
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)
_VERSION_PATTERN = re.compile(r"OrthoFinder\s+(?:version\s+)?([0-9]+(?:\.[0-9A-Za-z]+)+)", re.I)
_COMPLETION_MARKER = "OrthoFinder run completed"


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
    log_candidates = sorted(root.rglob("Log.txt"), key=lambda item: (len(item.parts), str(item)))
    if not log_candidates:
        raise InputValidationError(f"Completed OrthoFinder output lacks Log.txt: {root}")
    log_path = log_candidates[0]
    version = _read_version(log_path=log_path)
    major = _validate_supported_version(version=version)
    _require_completion_marker(log_path=log_path)
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
    primary = "HIERARCHICAL_ORTHOGROUP" if hog_paths else "ORTHOGROUP"
    layout = OrthoFinderLayout(
        results_dir=root,
        version=version,
        major_version=major,
        adapter_name="orthofinder_v3" if major == 3 else "orthofinder_v2_5_5",
        primary_group_authority=primary,
        log_path=log_path,
        orthogroups_path=orthogroups_path,
        hog_paths=hog_paths,
        species_ids_path=species_ids,
        sequence_ids_path=sequence_ids,
    )
    LOGGER.info(
        "Detected completed OrthoFinder %s output with %d HOG tables",
        version,
        len(hog_paths),
    )
    return layout


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
