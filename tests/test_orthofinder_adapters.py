"""Tests for raw and published-resource OrthoFinder adapters."""

from __future__ import annotations

import json
from pathlib import Path

import duckdb
import pytest

from protein_signatures.checksums import sha256_file
from protein_signatures.errors import InputValidationError
from protein_signatures.orthofinder import (
    _find_field,
    _first_or_none,
    _read_sequence_id_map,
    _read_version,
    _require_completion_marker,
    _split_members,
    _validate_supported_version,
    discover_orthofinder_layout,
    read_group_memberships,
)
from protein_signatures.orthofinder_resource import (
    _inspect_database,
    _manifest_integer,
    _membership_selection,
    _output_integer,
    _query_dicts,
    _query_rows,
    _safe_relative_path,
    discover_orthofinder_resource,
    read_resource_group_context,
    read_resource_memberships,
    resource_input_paths,
    verify_resource_outputs,
)


def test_raw_v2_hog_and_legacy_membership(tmp_path: Path) -> None:
    """Completed OrthoFinder 2.5.5 data should support both explicit authorities."""

    root = _raw_results(root=tmp_path / "results", version="2.5.5")
    layout = discover_orthofinder_layout(results_dir=root)
    assert layout.version == "2.5.5"
    assert layout.adapter_name == "orthofinder_v2_5_5"
    assert layout.to_record()["source_mode"] == "RAW_COMPLETED_RESULTS"
    hogs = read_group_memberships(
        layout=layout,
        run_id="source_run",
        group_type="HOG",
        hierarchy_node="N0",
        protein_ids=frozenset({"protA", "protB"}),
    )
    assert [(row.protein_id, row.group_id, row.run_id) for row in hogs] == [
        ("protA", "N0.HOG0001", "source_run"),
        ("protB", "N0.HOG0001", "source_run"),
    ]
    legacy = read_group_memberships(
        layout=layout,
        run_id="source_run",
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        protein_ids=frozenset({"protA", "protB"}),
    )
    assert {row.group_type for row in legacy} == {"LEGACY_ORTHOGROUP"}
    assert {row.group_id for row in legacy} == {"OG0001"}


def test_raw_v3_detection_and_identifier_helpers(tmp_path: Path) -> None:
    """OrthoFinder 3 should be recognised and parsing helpers should remain exact."""

    root = _raw_results(root=tmp_path / "results", version="3.1.0")
    layout = discover_orthofinder_layout(results_dir=root)
    assert layout.major_version == 3
    assert layout.primary_group_authority == "HIERARCHICAL_ORTHOGROUP"
    assert _read_version(log_path=root / "Log.txt") == "3.1.0"
    assert _read_sequence_id_map(path=None) == {}
    assert _split_members(value=" a, b ,, ") == ("a", "b")
    assert _split_members(value=None) == ()
    assert _find_field(fields=("x", "HOG"), candidates=("hog",)) == "HOG"
    assert _find_field(fields=("x",), candidates=("missing",), required=False) is None
    assert _first_or_none(paths=[Path("a"), Path("b")]) == Path("a")
    assert _first_or_none(paths=[]) is None
    _require_completion_marker(log_path=root / "Log.txt")
    assert _validate_supported_version(version="2.5.5") == 2
    assert _validate_supported_version(version="3.0.0-beta1") == 3


@pytest.mark.parametrize(
    ("mutation", "expected"),
    [
        ("missing_dir", "does not exist"),
        ("missing_log", "lacks Log.txt"),
        ("bad_version", "identify OrthoFinder version"),
        ("unsupported", "Unsupported OrthoFinder version"),
        ("unsupported_v2", "Unsupported OrthoFinder version"),
        ("unfinished", "official completion marker"),
        ("no_groups", "No Orthogroups.tsv"),
    ],
)
def test_raw_layout_rejects_incomplete_results(
    tmp_path: Path, mutation: str, expected: str
) -> None:
    """Layout discovery should reject unsupported or incomplete result trees."""

    root = tmp_path / mutation
    if mutation != "missing_dir":
        root = _raw_results(root=root, version="2.5.5")
    if mutation == "missing_log":
        (root / "Log.txt").unlink()
    elif mutation == "bad_version":
        (root / "Log.txt").write_text("unknown\n", encoding="utf-8")
    elif mutation == "unsupported":
        (root / "Log.txt").write_text(
            "OrthoFinder version 4.0.0\nOrthoFinder run completed\n", encoding="utf-8"
        )
    elif mutation == "unsupported_v2":
        (root / "Log.txt").write_text(
            "OrthoFinder version 2.6.0\nOrthoFinder run completed\n", encoding="utf-8"
        )
    elif mutation == "unfinished":
        (root / "Log.txt").write_text(
            "OrthoFinder version 2.5.5\nRun was not completed\n", encoding="utf-8"
        )
    elif mutation == "no_groups":
        (root / "Orthogroups" / "Orthogroups.tsv").unlink()
        (root / "Phylogenetic_Hierarchical_Orthogroups" / "N0.tsv").unlink()
    with pytest.raises(InputValidationError, match=expected):
        discover_orthofinder_layout(results_dir=root)


def test_raw_version_contract_rejects_other_v2_releases() -> None:
    """The tested v2 adapter must not silently claim support for other v2 releases."""

    for version in ("2.5.4", "2.6.0", "4.0.0"):
        with pytest.raises(InputValidationError, match="exactly OrthoFinder 2.5.5"):
            _validate_supported_version(version=version)


def test_raw_membership_rejects_bad_selection_and_mapping(tmp_path: Path) -> None:
    """Raw membership parsing should fail on unknown levels, types and ambiguity."""

    root = _raw_results(root=tmp_path / "results", version="2.5.5")
    layout = discover_orthofinder_layout(results_dir=root)
    with pytest.raises(InputValidationError, match="Unsupported OrthoFinder group type"):
        read_group_memberships(
            layout=layout,
            run_id="run",
            group_type="OTHER",
            hierarchy_node="N0",
            protein_ids=frozenset({"protA"}),
        )
    with pytest.raises(InputValidationError, match="requires an empty"):
        read_group_memberships(
            layout=layout,
            run_id="run",
            group_type="LEGACY_ORTHOGROUP",
            hierarchy_node="N0",
            protein_ids=frozenset({"protA"}),
        )
    with pytest.raises(InputValidationError, match="unavailable"):
        read_group_memberships(
            layout=layout,
            run_id="run",
            group_type="HOG",
            hierarchy_node="N9",
            protein_ids=frozenset({"protA"}),
        )
    with pytest.raises(InputValidationError, match="No campaign FASTA"):
        read_group_memberships(
            layout=layout,
            run_id="run",
            group_type="HOG",
            hierarchy_node="N0",
            protein_ids=frozenset({"absent"}),
        )
    sequence_ids = root / "WorkingDirectory" / "SequenceIDs.txt"
    sequence_ids.write_text("broken\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Malformed SequenceIDs"):
        _read_sequence_id_map(path=sequence_ids)
    sequence_ids.write_text("0_0: protA\n0_0: protB\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Duplicate internal"):
        _read_sequence_id_map(path=sequence_ids)
    with pytest.raises(InputValidationError, match="lack any required"):
        _find_field(fields=("a",), candidates=("b",))


def test_published_resource_imports_aliases_and_context(tmp_path: Path) -> None:
    """Schema-4 resources should preserve composite identities and explicit distances."""

    root = _published_resource(root=tmp_path / "resource", include_focus=True)
    resource = discover_orthofinder_resource(resource_dir=root)
    assert resource.schema_version == 4
    assert resource.focus_analysis_available is True
    assert resource.to_record()["source_mode"] == "ORTHOFINDER_RESULTS_RESOURCE"
    memberships = read_resource_memberships(
        resource=resource,
        group_type="HOG",
        hierarchy_node="N0",
        protein_ids=frozenset({"Q1", "protB"}),
    )
    assert {(row.protein_id, row.group_id) for row in memberships} == {
        ("Q1", "N0.HOG0001"),
        ("protB", "N0.HOG0001"),
    }
    context = read_resource_group_context(resource=resource, memberships=memberships)
    assert len(context) == 1
    assert context[0]["mean_distance"] == pytest.approx(0.25)
    assert context[0]["computation_status"] == "EXACT"
    inputs = resource_input_paths(resource=resource)
    assert resource.manifest_path in inputs
    assert resource.database_path in inputs
    verify_resource_outputs(resource_dir=root)


def test_published_resource_supports_internal_and_legacy_ids(tmp_path: Path) -> None:
    """Internal IDs and the separate legacy authority should map without substring use."""

    root = _published_resource(root=tmp_path / "resource", include_focus=False)
    resource = discover_orthofinder_resource(resource_dir=root, verify_checksums=False)
    internal = read_resource_memberships(
        resource=resource,
        group_type="HOG",
        hierarchy_node="N0",
        protein_ids=frozenset({"0_0"}),
    )
    assert internal[0].protein_id == "0_0"
    legacy = read_resource_memberships(
        resource=resource,
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        protein_ids=frozenset({"protB"}),
    )
    assert legacy[0].group_id == "OG0001"
    assert read_resource_group_context(resource=resource, memberships=()) == ()


def test_published_resource_rejects_unsafe_or_inconsistent_state(tmp_path: Path) -> None:
    """Manifest traversal, schema disagreement and ambiguous aliases must fail closed."""

    root = _published_resource(root=tmp_path / "resource", include_focus=False)
    manifest_path = root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"][0]["path"] = "../escape"
    with pytest.raises(InputValidationError, match="unsafe path"):
        verify_resource_outputs(resource_dir=root, manifest=manifest)
    manifest["outputs"][0]["path"] = "duckdb/orthofinder_results.duckdb"
    manifest["outputs"][0]["size_bytes"] += 1
    with pytest.raises(InputValidationError, match="size differs"):
        verify_resource_outputs(resource_dir=root, manifest=manifest)
    manifest["outputs"][0]["size_bytes"] -= 1
    manifest["outputs"][0]["sha256"] = "0" * 64
    with pytest.raises(InputValidationError, match="checksum differs"):
        verify_resource_outputs(resource_dir=root, manifest=manifest)
    manifest["schema_version"] = 2
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputValidationError, match="Unsupported"):
        discover_orthofinder_resource(resource_dir=root, verify_checksums=False)


def test_resource_helpers_validate_values_and_query_errors(tmp_path: Path) -> None:
    """Small resource helpers should reject invalid values and database operations."""

    assert _membership_selection(group_type="HOG", hierarchy_node="N0") == (
        "hog_memberships",
        "N0",
    )
    assert _membership_selection(group_type="legacy_orthogroup", hierarchy_node="") == (
        "legacy_orthogroup_memberships",
        "",
    )
    with pytest.raises(InputValidationError):
        _membership_selection(group_type="OTHER", hierarchy_node="")
    with pytest.raises(InputValidationError):
        _membership_selection(group_type="LEGACY_ORTHOGROUP", hierarchy_node="N0")
    assert _manifest_integer(manifest={"x": "4"}, field="x") == 4
    assert _output_integer(output={"x": "0"}, field="x", index=0) == 0
    assert _safe_relative_path(value="tables/a.tsv", index=0) == Path("tables/a.tsv")
    for value in (True, None, "bad"):
        with pytest.raises(InputValidationError):
            _manifest_integer(manifest={"x": value}, field="x")
    for value in (True, None, "bad", -1):
        with pytest.raises(InputValidationError):
            _output_integer(output={"x": value}, field="x", index=0)
    for value in (None, "", "/absolute", "../escape", "."):
        with pytest.raises(InputValidationError):
            _safe_relative_path(value=value, index=0)
    with pytest.raises(InputValidationError, match="missing or empty"):
        _inspect_database(path=tmp_path / "missing.duckdb")
    invalid = tmp_path / "invalid.duckdb"
    invalid.write_bytes(b"not duckdb")
    with pytest.raises(InputValidationError, match="Could not inspect"):
        _inspect_database(path=invalid)
    with pytest.raises(InputValidationError, match="Could not read"):
        _query_rows(database=invalid, sql="SELECT 1", parameters=())
    with pytest.raises(InputValidationError, match="Could not read"):
        _query_dicts(database=invalid, sql="SELECT 1", parameters=())


def test_resource_queries_cannot_execute_external_file_views(tmp_path: Path) -> None:
    """An upstream DuckDB view must not read files through trusted adapter SQL."""

    secret = tmp_path / "outside.txt"
    secret.write_text("must-not-be-read\n", encoding="utf-8")
    database = tmp_path / "malicious.duckdb"
    escaped = str(secret).replace("'", "''")
    with duckdb.connect(str(database)) as connection:
        connection.execute(
            f"CREATE VIEW sequences AS SELECT content AS member_id FROM read_text('{escaped}')"
        )
    with pytest.raises(InputValidationError, match="Could not read"):
        _query_rows(database=database, sql="SELECT * FROM sequences", parameters=())


@pytest.mark.parametrize("mutation", ["absent", "duplicate", "size", "checksum"])
def test_preferred_resource_database_is_always_manifest_bound(
    tmp_path: Path,
    mutation: str,
) -> None:
    """Skipping the full output pass must still authenticate the consumed DuckDB."""

    root = _published_resource(root=tmp_path / mutation, include_focus=False)
    manifest_path = root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if mutation == "absent":
        manifest["outputs"] = []
    elif mutation == "duplicate":
        manifest["outputs"].append(dict(manifest["outputs"][0]))
    elif mutation == "size":
        manifest["outputs"][0]["size_bytes"] += 1
    else:
        manifest["outputs"][0]["sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputValidationError, match="manifest|declare"):
        discover_orthofinder_resource(resource_dir=root, verify_checksums=False)


def _raw_results(*, root: Path, version: str) -> Path:
    """Create a tiny completed raw OrthoFinder result tree for tests."""

    (root / "WorkingDirectory").mkdir(parents=True)
    (root / "Orthogroups").mkdir()
    (root / "Phylogenetic_Hierarchical_Orthogroups").mkdir()
    (root / "Log.txt").write_text(
        f"OrthoFinder version {version}\nOrthoFinder run completed\n", encoding="utf-8"
    )
    (root / "WorkingDirectory" / "SequenceIDs.txt").write_text(
        "0_0: protA description\n1_0: protB\n", encoding="utf-8"
    )
    (root / "WorkingDirectory" / "SpeciesIDs.txt").write_text(
        "0: Species_A.faa\n1: Species_B.faa\n", encoding="utf-8"
    )
    (root / "Phylogenetic_Hierarchical_Orthogroups" / "N0.tsv").write_text(
        "HOG\tOG\tGene Tree Parent Clade\tSpecies_A\tSpecies_B\nN0.HOG0001\tOG0001\tN0\t0_0\t1_0\n",
        encoding="utf-8",
    )
    (root / "Orthogroups" / "Orthogroups.tsv").write_text(
        "Orthogroup\tSpecies_A\tSpecies_B\nOG0001\t0_0\t1_0\n",
        encoding="utf-8",
    )
    return root


def _published_resource(*, root: Path, include_focus: bool) -> Path:
    """Create a minimal schema-4 companion resource with a physical DuckDB."""

    database = root / "duckdb" / "orthofinder_results.duckdb"
    database.parent.mkdir(parents=True)
    with duckdb.connect(str(database)) as connection:
        connection.execute("CREATE TABLE resource_metadata AS SELECT 4::INTEGER AS schema_version")
        connection.execute(
            "CREATE TABLE group_statistics AS SELECT * FROM (VALUES "
            "('run_1','HOG','N0','N0.HOG0001','OG0001','N0',2,2,2,1,1.0,false),"
            "('run_1','LEGACY_ORTHOGROUP','','OG0001','OG0001','',2,2,2,1,1.0,false)) "
            "t(run_id,group_type,hierarchy_node,group_id,legacy_orthogroup_id,"
            "gene_tree_parent_clade,member_count,species_count,single_copy_species_count,"
            "max_copies_per_species,mean_copies_per_species,is_singleton)"
        )
        connection.execute(
            "CREATE TABLE hog_memberships AS SELECT * FROM (VALUES "
            "('run_1','HOG','N0','N0.HOG0001','OG0001','N0','Species_A',"
            "'sp|Q1|ENTRY1'),"
            "('run_1','HOG','N0','N0.HOG0001','OG0001','N0','Species_B','protB')) "
            "t(run_id,group_type,hierarchy_node,group_id,legacy_orthogroup_id,"
            "gene_tree_parent_clade,species_label,member_id)"
        )
        connection.execute(
            "CREATE TABLE legacy_orthogroup_memberships AS SELECT * FROM (VALUES "
            "('run_1','LEGACY_ORTHOGROUP','','OG0001','OG0001','','Species_A',"
            "'sp|Q1|ENTRY1'),"
            "('run_1','LEGACY_ORTHOGROUP','','OG0001','OG0001','','Species_B','protB')) "
            "t(run_id,group_type,hierarchy_node,group_id,legacy_orthogroup_id,"
            "gene_tree_parent_clade,species_label,member_id)"
        )
        connection.execute(
            "CREATE TABLE sequences AS SELECT * FROM (VALUES "
            "('run_1','0_0','Species_A','sp|Q1|ENTRY1'),"
            "('run_1','1_0','Species_B','protB')) "
            "t(run_id,internal_id,species_label,member_id)"
        )
        connection.execute(
            "CREATE TABLE distance_statistics AS SELECT * FROM (VALUES "
            "('run_1','HOG','N0','N0.HOG0001','RESOLVED_GENE_TREE','EXACT',2,2,1,0,"
            "0.25,0.25,0.25,0.25,0.25,0.25,0.0,''),"
            "('run_1','LEGACY_ORTHOGROUP','','OG0001','DISABLED','UNAVAILABLE',2,0,0,0,"
            "NULL,NULL,NULL,NULL,NULL,NULL,NULL,'disabled')) "
            "t(run_id,group_type,hierarchy_node,group_id,distance_method,"
            "computation_status,total_member_count,sampled_member_count,distance_pair_count,"
            "unresolved_pair_count,minimum_distance,q25_distance,median_distance,mean_distance,"
            "q75_distance,maximum_distance,population_stddev_distance,failure_reason)"
        )
        if include_focus:
            connection.execute("CREATE TABLE e3_cluster_results AS SELECT 1 AS x")
            connection.execute("CREATE TABLE e3_seed_catalogue_audit AS SELECT 1 AS x")
            connection.execute("CREATE TABLE e3_seed_matches AS SELECT 1 AS x")
        connection.execute("CHECKPOINT")
    output = {
        "path": "duckdb/orthofinder_results.duckdb",
        "size_bytes": database.stat().st_size,
        "sha256": sha256_file(path=database),
    }
    manifest = {
        "status": "complete",
        "schema_version": 4,
        "run_id": "run_1",
        "package_version": "0.9.1",
        "orthofinder_version": "2.5.5",
        "adapter_name": "orthofinder_2",
        "primary_group_authority": "HOG",
        "outputs": [output],
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root
