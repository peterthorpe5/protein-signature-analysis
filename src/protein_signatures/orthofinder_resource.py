"""Read-only adapter for completed ``orthofinder-results`` resources.

The companion project owns OrthoFinder parsing, taxonomy, trees, evolutionary
distances and focus-seed matching.  This module verifies that published
authority and imports only the relationships needed for leakage-safe signature
analysis.  It deliberately has no import-time dependency on that project.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import duckdb

from .checksums import sha256_file
from .errors import InputValidationError
from .io_utils import read_json
from .models import GroupMembership, OrthoFinderResource
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)
SUPPORTED_RESOURCE_SCHEMAS = frozenset({3, 4})
REQUIRED_RELATIONS = frozenset(
    {
        "distance_statistics",
        "group_statistics",
        "hog_memberships",
        "legacy_orthogroup_memberships",
        "resource_metadata",
        "sequences",
    }
)
FOCUS_RELATIONS = frozenset(
    {
        "e3_cluster_results",
        "e3_seed_catalogue_audit",
        "e3_seed_matches",
    }
)
_UNIPROT_IDENTIFIER = re.compile(r"^(?:sp|tr)\|([^|\s]+)\|([^|\s]+)$")
_DUCKDB_READ_ONLY_CONFIG = {
    "allow_community_extensions": "false",
    "allow_persistent_secrets": "false",
    "allow_unsigned_extensions": "false",
    "autoinstall_known_extensions": "false",
    "autoload_known_extensions": "false",
    "enable_external_access": "false",
    "lock_configuration": "true",
}


def discover_orthofinder_resource(
    *, resource_dir: Path, verify_checksums: bool = True
) -> OrthoFinderResource:
    """Validate and describe a completed ``orthofinder-results`` resource.

    Args:
        resource_dir: Published resource root containing ``run_manifest.json``.
        verify_checksums: Verify every output declared by the upstream manifest.

    Returns:
        Immutable compatible resource identity.

    Raises:
        InputValidationError: If completion, checksums or schemas are invalid.
    """

    root = Path(resource_dir).expanduser().resolve()
    if not root.is_dir():
        raise InputValidationError(f"OrthoFinder resource directory does not exist: {root}")
    manifest_path = root / "run_manifest.json"
    manifest = read_json(path=manifest_path)
    if not isinstance(manifest, dict) or manifest.get("status") != "complete":
        raise InputValidationError(
            f"OrthoFinder resource manifest is not marked complete: {manifest_path}"
        )
    schema_version = _manifest_integer(manifest=manifest, field="schema_version")
    if schema_version not in SUPPORTED_RESOURCE_SCHEMAS:
        supported = ", ".join(str(value) for value in sorted(SUPPORTED_RESOURCE_SCHEMAS))
        raise InputValidationError(
            f"Unsupported orthofinder-results schema {schema_version}; supported: {supported}."
        )
    _verify_required_database_output(resource_dir=root, manifest=manifest)
    if verify_checksums:
        verify_resource_outputs(resource_dir=root, manifest=manifest)
    database_path = root / "duckdb" / "orthofinder_results.duckdb"
    database_schema, relations, database_run_ids = _inspect_database(path=database_path)
    missing = sorted(REQUIRED_RELATIONS - relations)
    if missing:
        raise InputValidationError(
            "OrthoFinder resource lacks required DuckDB relations: " + "; ".join(missing)
        )
    if database_schema != schema_version:
        raise InputValidationError(
            "OrthoFinder manifest and DuckDB schema versions disagree: "
            f"{schema_version} versus {database_schema}."
        )
    if len(database_run_ids) != 1:
        raise InputValidationError(
            "OrthoFinder resource must contain exactly one run_id; observed: "
            + ("; ".join(database_run_ids) if database_run_ids else "none")
        )
    run_id = validate_identifier(
        value=manifest.get("run_id"), field_name="orthofinder resource run_id"
    )
    if database_run_ids[0] != run_id:
        raise InputValidationError(
            "OrthoFinder manifest and DuckDB run identifiers disagree: "
            f"{run_id!r} versus {database_run_ids[0]!r}."
        )
    resource = OrthoFinderResource(
        resource_dir=root,
        manifest_path=manifest_path,
        database_path=database_path,
        run_id=run_id,
        schema_version=schema_version,
        package_version=str(manifest.get("package_version", "unknown")),
        orthofinder_version=str(manifest.get("orthofinder_version", "unknown")),
        adapter_name=str(manifest.get("adapter_name", "unknown")),
        primary_group_authority=str(manifest.get("primary_group_authority", "unknown")),
        relations=frozenset(relations),
        focus_analysis_available=FOCUS_RELATIONS <= relations,
    )
    LOGGER.info(
        "Validated orthofinder-results resource run=%s schema=%d relations=%d",
        resource.run_id,
        resource.schema_version,
        len(resource.relations),
    )
    return resource


def verify_resource_outputs(
    *, resource_dir: Path, manifest: Mapping[str, Any] | None = None
) -> None:
    """Verify every upstream output against its manifest size and checksum.

    Args:
        resource_dir: Completed resource root.
        manifest: Optional already-decoded run manifest.

    Raises:
        InputValidationError: If the manifest or any declared output is invalid.
    """

    root = Path(resource_dir).expanduser().resolve()
    record = read_json(path=root / "run_manifest.json") if manifest is None else manifest
    if not isinstance(record, Mapping):
        raise InputValidationError("OrthoFinder run_manifest.json must contain an object.")
    outputs = record.get("outputs")
    if not isinstance(outputs, list):
        raise InputValidationError("OrthoFinder run manifest lacks an outputs list.")
    seen: set[Path] = set()
    for index, output in enumerate(outputs):
        if not isinstance(output, Mapping):
            raise InputValidationError(f"OrthoFinder output record {index} must contain an object.")
        relative = _safe_relative_path(value=output.get("path"), index=index)
        if relative in seen:
            raise InputValidationError(f"OrthoFinder output manifest repeats path: {relative}")
        seen.add(relative)
        path = (root / relative).resolve()
        if root not in path.parents or not path.is_file():
            raise InputValidationError(
                f"OrthoFinder resource is missing declared output: {relative}"
            )
        expected_size = _output_integer(output=output, field="size_bytes", index=index)
        if path.stat().st_size != expected_size:
            raise InputValidationError(f"OrthoFinder output size differs from manifest: {relative}")
        expected_sha = str(output.get("sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{64}", expected_sha):
            raise InputValidationError(f"OrthoFinder output record {index} has an invalid SHA-256.")
        if sha256_file(path=path) != expected_sha:
            raise InputValidationError(
                f"OrthoFinder output checksum differs from manifest: {relative}"
            )


def _verify_required_database_output(*, resource_dir: Path, manifest: Mapping[str, Any]) -> None:
    """Bind the exact preferred DuckDB to one valid manifest record.

    This verification is mandatory even when the caller skips the potentially
    expensive checksum pass over every upstream output.

    Args:
        resource_dir: Completed resource root.
        manifest: Decoded upstream run manifest.

    Raises:
        InputValidationError: If the database record is absent, repeated, unsafe,
            malformed or differs from the physical database.
    """

    outputs = manifest.get("outputs")
    if not isinstance(outputs, list):
        raise InputValidationError("OrthoFinder run manifest lacks an outputs list.")
    required = Path("duckdb/orthofinder_results.duckdb")
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, output in enumerate(outputs):
        if not isinstance(output, Mapping):
            raise InputValidationError(f"OrthoFinder output record {index} must contain an object.")
        relative = _safe_relative_path(value=output.get("path"), index=index)
        if relative == required:
            matches.append((index, output))
    if len(matches) != 1:
        raise InputValidationError(
            "OrthoFinder manifest must declare exactly one "
            "duckdb/orthofinder_results.duckdb output."
        )
    index, output = matches[0]
    candidate = resource_dir / required
    if candidate.is_symlink():
        raise InputValidationError("OrthoFinder resource DuckDB must not be a symbolic link.")
    try:
        source = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputValidationError("OrthoFinder resource DuckDB is missing or unsafe.") from error
    if resource_dir not in source.parents or not source.is_file():
        raise InputValidationError("OrthoFinder resource DuckDB is missing or unsafe.")
    expected_size = _output_integer(output=output, field="size_bytes", index=index)
    if source.stat().st_size != expected_size:
        raise InputValidationError("OrthoFinder resource DuckDB size differs from manifest.")
    expected_sha = str(output.get("sha256", ""))
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha):
        raise InputValidationError("OrthoFinder resource DuckDB manifest SHA-256 is invalid.")
    if sha256_file(path=source) != expected_sha:
        raise InputValidationError("OrthoFinder resource DuckDB checksum differs from manifest.")


def read_resource_memberships(
    *,
    resource: OrthoFinderResource,
    group_type: str,
    hierarchy_node: str,
    protein_ids: frozenset[str],
) -> tuple[GroupMembership, ...]:
    """Import exact campaign memberships from one published group collection.

    Controlled aliases from the upstream sequence authority are accepted:
    canonical member ID, OrthoFinder internal ID, and exact accession/entry
    components of a syntactically valid UniProt identifier.  Substrings are
    never matched.

    Args:
        resource: Validated upstream resource.
        group_type: ``HOG`` or ``LEGACY_ORTHOGROUP``.
        hierarchy_node: Exact HOG node, or blank for legacy orthogroups.
        protein_ids: Authoritative campaign FASTA identifiers.

    Returns:
        Unique memberships keyed by the upstream composite identity.

    Raises:
        InputValidationError: If selection or identifier mapping is ambiguous.
    """

    relation, node = _membership_selection(group_type=group_type, hierarchy_node=hierarchy_node)
    if not protein_ids:
        return ()
    aliases = _campaign_alias_map(resource=resource, protein_ids=protein_ids)
    source_members = sorted(aliases)
    if not source_members:
        raise InputValidationError(
            "No campaign FASTA identifiers matched the orthofinder-results sequence authority."
        )
    sql = (
        "SELECT run_id, group_type, hierarchy_node, group_id, "
        "legacy_orthogroup_id, gene_tree_parent_clade, species_label, member_id "
        f"FROM {relation} WHERE run_id = ? AND group_type = ? AND hierarchy_node = ? "
        "AND member_id IN (SELECT unnest(?::VARCHAR[])) "
        "ORDER BY group_id, species_label, member_id"
    )
    rows = _query_rows(
        database=resource.database_path,
        sql=sql,
        parameters=(resource.run_id, group_type, node, source_members),
    )
    memberships: dict[tuple[str, str, str, str], GroupMembership] = {}
    group_by_protein: dict[str, tuple[str, str, str, str]] = {}
    for row in rows:
        source_member = str(row[7])
        campaign_id = aliases[source_member]
        key = (str(row[0]), str(row[1]), str(row[2]), str(row[3]))
        previous = group_by_protein.get(campaign_id)
        if previous is not None and previous != key:
            raise InputValidationError(
                f"Campaign protein {campaign_id!r} maps to multiple selected groups: "
                f"{previous!r} and {key!r}."
            )
        group_by_protein[campaign_id] = key
        membership_key = (*key, campaign_id)
        memberships[membership_key] = GroupMembership(
            run_id=key[0],
            group_type=key[1],
            hierarchy_node=key[2],
            group_id=key[3],
            legacy_orthogroup_id=str(row[4] or ""),
            gene_tree_parent_clade=str(row[5] or ""),
            species_label=str(row[6]),
            protein_id=campaign_id,
        )
    if not memberships:
        raise InputValidationError(
            "Matched sequence identifiers had no membership at the requested group level."
        )
    LOGGER.info(
        "Imported %d campaign memberships from upstream run %s (%s/%s)",
        len(memberships),
        resource.run_id,
        group_type,
        node or "ROOT",
    )
    return tuple(memberships[key] for key in sorted(memberships))


def read_resource_group_context(
    *, resource: OrthoFinderResource, memberships: Iterable[GroupMembership]
) -> tuple[dict[str, Any], ...]:
    """Read upstream group and evolutionary summaries for selected memberships.

    Args:
        resource: Validated upstream resource.
        memberships: Imported campaign memberships.

    Returns:
        One portable context row per selected composite group identity.
    """

    keys = sorted(
        {
            (
                membership.run_id,
                membership.group_type,
                membership.hierarchy_node,
                membership.group_id,
            )
            for membership in memberships
        }
    )
    if not keys:
        return ()
    values_sql = ", ".join("(?, ?, ?, ?)" for _ in keys)
    parameters = tuple(value for key in keys for value in key)
    distance_join = (
        "LEFT JOIN distance_statistics AS d USING (run_id, group_type, hierarchy_node, group_id)"
        if "distance_statistics" in resource.relations
        else ""
    )
    distance_columns = (
        "d.distance_method, d.computation_status, d.total_member_count, "
        "d.sampled_member_count, d.distance_pair_count, d.unresolved_pair_count, "
        "d.minimum_distance, d.q25_distance, d.median_distance, d.mean_distance, "
        "d.q75_distance, d.maximum_distance, d.population_stddev_distance, "
        "d.failure_reason"
        if distance_join
        else "'' AS distance_method, '' AS computation_status, NULL, NULL, NULL, "
        "NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, '' AS failure_reason"
    )
    sql = (
        "WITH selected(run_id, group_type, hierarchy_node, group_id) AS (VALUES "
        f"{values_sql}) SELECT g.run_id, g.group_type, g.hierarchy_node, g.group_id, "
        "g.legacy_orthogroup_id, g.gene_tree_parent_clade, g.member_count, "
        "g.species_count, g.single_copy_species_count, g.max_copies_per_species, "
        f"g.mean_copies_per_species, g.is_singleton, {distance_columns} "
        "FROM group_statistics AS g JOIN selected AS s USING "
        "(run_id, group_type, hierarchy_node, group_id) "
        f"{distance_join} ORDER BY g.run_id, g.group_type, g.hierarchy_node, g.group_id"
    )
    rows = _query_dicts(
        database=resource.database_path,
        sql=sql,
        parameters=parameters,
    )
    if len(rows) != len(keys):
        raise InputValidationError(
            "One or more imported OrthoFinder memberships lack a unique group-statistics row."
        )
    return tuple(rows)


def resource_input_paths(*, resource: OrthoFinderResource) -> tuple[Path, ...]:
    """Return the minimal upstream authorities bound into this run's provenance.

    Args:
        resource: Validated upstream resource.

    Returns:
        Manifest, database and relevant analytical tables when present.
    """

    candidates = [
        resource.manifest_path,
        resource.database_path,
        resource.resource_dir / "tables" / "hog_memberships.tsv.gz",
        resource.resource_dir / "tables" / "legacy_orthogroup_memberships.tsv.gz",
        resource.resource_dir / "tables" / "group_statistics.tsv.gz",
        resource.resource_dir / "tables" / "distance_statistics.tsv.gz",
        resource.resource_dir / "tables" / "e3_cluster_results.tsv.gz",
        resource.resource_dir / "tables" / "e3_seed_matches.tsv.gz",
    ]
    return tuple(path for path in candidates if path.is_file())


def _campaign_alias_map(
    *, resource: OrthoFinderResource, protein_ids: frozenset[str]
) -> dict[str, str]:
    """Map upstream canonical member IDs to unique campaign identifiers."""

    identifiers = sorted(protein_ids)
    sql = (
        "SELECT DISTINCT member_id, internal_id FROM sequences WHERE run_id = ? AND "
        "(member_id IN (SELECT unnest(?::VARCHAR[])) OR "
        "internal_id IN (SELECT unnest(?::VARCHAR[])) OR "
        "regexp_extract(member_id, '^(?:sp|tr)\\|([^|[:space:]]+)\\|', 1) "
        "IN (SELECT unnest(?::VARCHAR[])) OR "
        "regexp_extract(member_id, '^(?:sp|tr)\\|[^|[:space:]]+\\|([^|[:space:]]+)$', 1) "
        "IN (SELECT unnest(?::VARCHAR[]))) ORDER BY member_id, internal_id"
    )
    rows = _query_rows(
        database=resource.database_path,
        sql=sql,
        parameters=(resource.run_id, identifiers, identifiers, identifiers, identifiers),
    )
    candidates: dict[str, set[str]] = {}
    for member_value, internal_value in rows:
        member_id = str(member_value)
        aliases = {member_id, str(internal_value or "")}
        match = _UNIPROT_IDENTIFIER.fullmatch(member_id)
        if match is not None:
            aliases.update(match.groups())
        matched = aliases & protein_ids
        if matched:
            candidates.setdefault(member_id, set()).update(matched)
    for identifier in protein_ids:
        candidates.setdefault(identifier, set()).add(identifier)
    ambiguous = {member: sorted(values) for member, values in candidates.items() if len(values) > 1}
    if ambiguous:
        first_member = sorted(ambiguous)[0]
        raise InputValidationError(
            "Campaign FASTA identifiers ambiguously match one upstream member: "
            f"{first_member!r} -> {ambiguous[first_member]!r}."
        )
    return {member: next(iter(values)) for member, values in candidates.items() if values}


def _membership_selection(*, group_type: str, hierarchy_node: str) -> tuple[str, str]:
    """Validate a group authority and return its static relation and node."""

    authority = group_type.strip().upper()
    if authority == "HOG":
        node = validate_identifier(value=hierarchy_node, field_name="hierarchy_node")
        return "hog_memberships", node
    if authority == "LEGACY_ORTHOGROUP":
        if hierarchy_node:
            raise InputValidationError(
                "LEGACY_ORTHOGROUP membership requires an empty hierarchy_node."
            )
        return "legacy_orthogroup_memberships", ""
    raise InputValidationError(
        f"Unsupported OrthoFinder group_type {group_type!r}; expected HOG or LEGACY_ORTHOGROUP."
    )


def _inspect_database(*, path: Path) -> tuple[int, frozenset[str], tuple[str, ...]]:
    """Return schema, relations and run IDs from an upstream DuckDB."""

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"OrthoFinder resource DuckDB is missing or empty: {source}")
    try:
        with duckdb.connect(
            str(source),
            read_only=True,
            config=_DUCKDB_READ_ONLY_CONFIG,
        ) as connection:
            relations = frozenset(
                str(row[0]) for row in connection.execute("SHOW TABLES").fetchall()
            )
            if "resource_metadata" not in relations or "group_statistics" not in relations:
                return -1, relations, ()
            schema_rows = connection.execute(
                "SELECT DISTINCT schema_version FROM resource_metadata"
            ).fetchall()
            run_rows = connection.execute(
                "SELECT DISTINCT run_id FROM group_statistics ORDER BY run_id"
            ).fetchall()
    except (duckdb.Error, OSError) as error:
        raise InputValidationError(
            f"Could not inspect OrthoFinder resource DuckDB {source}: {error}"
        ) from error
    if len(schema_rows) != 1:
        raise InputValidationError(
            "OrthoFinder resource metadata must contain exactly one schema version."
        )
    try:
        schema_version = int(schema_rows[0][0])
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            f"OrthoFinder resource schema version is invalid: {schema_rows[0][0]!r}."
        ) from error
    return schema_version, relations, tuple(str(row[0]) for row in run_rows)


def _query_rows(*, database: Path, sql: str, parameters: tuple[Any, ...]) -> list[tuple[Any, ...]]:
    """Execute one internal read-only query and return tuples."""

    try:
        with duckdb.connect(
            str(database),
            read_only=True,
            config=_DUCKDB_READ_ONLY_CONFIG,
        ) as connection:
            return connection.execute(sql, parameters).fetchall()
    except duckdb.Error as error:
        raise InputValidationError(
            f"Could not read OrthoFinder resource {database}: {error}"
        ) from error


def _query_dicts(*, database: Path, sql: str, parameters: tuple[Any, ...]) -> list[dict[str, Any]]:
    """Execute one internal read-only query and return named records."""

    try:
        with duckdb.connect(
            str(database),
            read_only=True,
            config=_DUCKDB_READ_ONLY_CONFIG,
        ) as connection:
            cursor = connection.execute(sql, parameters)
            headings = tuple(item[0] for item in cursor.description)
            return [dict(zip(headings, row, strict=True)) for row in cursor.fetchall()]
    except duckdb.Error as error:
        raise InputValidationError(
            f"Could not read OrthoFinder resource {database}: {error}"
        ) from error


def _manifest_integer(*, manifest: Mapping[str, Any], field: str) -> int:
    """Return one non-Boolean integer manifest field."""

    value = manifest.get(field)
    if isinstance(value, bool):
        raise InputValidationError(f"OrthoFinder manifest {field} must be an integer.")
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            f"OrthoFinder manifest {field} must be an integer: {value!r}."
        ) from error


def _output_integer(*, output: Mapping[str, Any], field: str, index: int) -> int:
    """Return one non-negative integer from an output manifest record."""

    value = output.get(field)
    if isinstance(value, bool):
        raise InputValidationError(
            f"OrthoFinder output record {index} field {field} must be an integer."
        )
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            f"OrthoFinder output record {index} field {field} is invalid: {value!r}."
        ) from error
    if parsed < 0:
        raise InputValidationError(
            f"OrthoFinder output record {index} field {field} must be non-negative."
        )
    return parsed


def _safe_relative_path(*, value: Any, index: int) -> Path:
    """Validate one portable output path from an upstream manifest."""

    if not isinstance(value, str) or not value:
        raise InputValidationError(f"OrthoFinder output record {index} lacks a path string.")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path == Path("."):
        raise InputValidationError(
            f"OrthoFinder output record {index} has an unsafe path: {value!r}."
        )
    return path
