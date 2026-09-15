"""Tests for importing completed US-align/TM-align structural resources."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from protein_signatures.checksums import sha256_file, sha256_json, sha256_text
from protein_signatures.errors import InputValidationError
from protein_signatures.models import SequenceRecord
from protein_signatures.structural_resource import (
    _boolean,
    _non_negative_integer,
    _optional_integer,
    _optional_number,
    _parse_comparison_status,
    _read_global_alignments,
    _read_parquet,
    _safe_path,
    import_structural_alignment_resource,
    resolve_structural_resource_root,
    verify_structural_resource_outputs,
)


def _alignment_row(**changes: object) -> dict[str, object]:
    """Return one valid explicit structural-reference sentinel row."""

    row: dict[str, object] = {
        "cluster_id": "cluster_1",
        "reference_accession": "p1",
        "mobile_accession": "p1",
        "alignment_tool": "TM-align",
        "status": "REFERENCE",
        "tool_version": "20240303",
        "aligned_length": None,
        "rmsd_angstrom": 0.0,
        "minimum_tm_score": 1.0,
    }
    row.update(changes)
    return row


def test_import_structural_resource_converts_global_and_pocket_evidence(
    tmp_path: Path,
) -> None:
    """Completed upstream tables should yield coverage-aware generic evidence."""

    root = _structural_resource(root=tmp_path / "structural")
    imported = import_structural_alignment_resource(resource_dir=root, sequences=_sequences())
    assert imported.package_version == "0.6.0"
    assert imported.run_digest == "a" * 64
    assert imported.reference_membership_row_count == 2
    assert len(imported.comparisons) == 2
    usalign = imported.comparisons[0]
    assert usalign.tm_score == pytest.approx(0.8)
    assert usalign.coverage_a == pytest.approx(0.8)
    assert usalign.coverage_b == pytest.approx(0.5)
    assert usalign.comparison_universe_id.startswith("ES3A_")
    assert usalign.coverage_scope.value == "FULL_SEQUENCE"
    assert imported.comparison_universe_members == {
        usalign.comparison_universe_id: frozenset({"p1", "p2"})
    }
    assert {item.feature_id.rsplit(":", 1)[1] for item in imported.features} == {
        "CONSERVED_3D_POCKET",
        "SAME_3D_POCKET_POSITION",
    }
    assert len({item.feature_id.split(":", 1)[0] for item in imported.features}) == 2
    assert {item.derivation_scope for item in imported.features} == {"ALL_DATA_EXPLORATORY"}
    assert all(len(item.feature_definition_sha256) == 64 for item in imported.features)
    assert all(item.derivation_cohort_sha256 == "" for item in imported.features)
    assert {item.protein_id for item in imported.features} == {"p1", "p2"}
    assert imported.group_summaries[0]["alignment_status"] == ("CONSERVED_3D_POCKET_SUPPORTED")
    assert imported.group_summaries[0]["aligned_accession_count"] == 2
    verify_structural_resource_outputs(resource_dir=root)


def test_structural_resource_resolves_end_to_end_parent(tmp_path: Path) -> None:
    """An end-to-end root should resolve only its conventional structural child."""

    parent = tmp_path / "run"
    root = _structural_resource(root=parent / "09b_structural_alignment" / "structural_alignment")
    assert resolve_structural_resource_root(path=parent) == root
    assert resolve_structural_resource_root(path=root) == root
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(InputValidationError, match="found 0"):
        resolve_structural_resource_root(path=empty)
    with pytest.raises(InputValidationError, match="does not exist"):
        resolve_structural_resource_root(path=tmp_path / "missing")


def test_import_supports_checksum_bound_end_to_end_aggregate(tmp_path: Path) -> None:
    """A real Stage 09b datasets manifest should import via its outer inventory."""

    root = _structural_resource(root=tmp_path / "09b_structural_alignment" / "structural_alignment")
    aggregate = _convert_to_aggregate(root=root)
    imported = import_structural_alignment_resource(
        resource_dir=tmp_path,
        sequences=_sequences(),
    )
    expected_datasets = {
        name.removesuffix(".parquet"): sha256_file(path=root / "tables" / name)
        for name in _REQUIRED_NAMES
    }
    assert imported.run_digest == sha256_json(
        value={
            "schema": "e3workflow_stage_09b_aggregate_v1",
            "configuration_digest": "b" * 64,
            "datasets": expected_datasets,
        }
    )
    assert imported.package_version == "0.16.0"
    assert imported.input_paths[:2] == (
        aggregate,
        root.parent / "stage_manifest.json",
    )
    assert len(imported.comparisons) == 2
    verify_structural_resource_outputs(resource_dir=root)


def test_aggregate_manifest_validation_rejects_unbound_evidence(tmp_path: Path) -> None:
    """Stage 09b datasets must agree with the outer checksum authority."""

    root = _structural_resource(
        root=tmp_path / "first" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    stage_manifest = root.parent / "stage_manifest.json"
    stage_manifest.unlink()
    with pytest.raises(InputValidationError, match="Missing or empty input file"):
        verify_structural_resource_outputs(resource_dir=root)

    root = _structural_resource(
        root=tmp_path / "second" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    stage_manifest = root.parent / "stage_manifest.json"
    outer = json.loads(stage_manifest.read_text(encoding="utf-8"))
    outer["status"] = "running"
    stage_manifest.write_text(json.dumps(outer), encoding="utf-8")
    with pytest.raises(InputValidationError, match="not marked complete"):
        verify_structural_resource_outputs(resource_dir=root)

    root = _structural_resource(
        root=tmp_path / "third" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    stage_manifest = root.parent / "stage_manifest.json"
    outer = json.loads(stage_manifest.read_text(encoding="utf-8"))
    outer["outputs"] = [
        row
        for row in outer["outputs"]
        if row["path"] != "structural_alignment/provenance/run_manifest.json"
    ]
    stage_manifest.write_text(json.dumps(outer), encoding="utf-8")
    with pytest.raises(InputValidationError, match="not checksum-inventoried"):
        verify_structural_resource_outputs(resource_dir=root)

    root = _structural_resource(
        root=tmp_path / "fourth" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    manifest = json.loads(aggregate.read_text(encoding="utf-8"))
    manifest["datasets"]["structural_alignments"]["sha256"] = "0" * 64
    aggregate.write_text(json.dumps(manifest), encoding="utf-8")
    _write_aggregate_stage_manifest(root=root)
    with pytest.raises(InputValidationError, match="checksum disagrees"):
        verify_structural_resource_outputs(resource_dir=root)


def test_aggregate_manifest_rejects_bad_identity_and_required_dataset(
    tmp_path: Path,
) -> None:
    """Aggregate paths, configuration identity and required tables must be explicit."""

    root = _structural_resource(
        root=tmp_path / "first" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    manifest = json.loads(aggregate.read_text(encoding="utf-8"))
    manifest["datasets"]["structural_alignments"]["path"] = "/unsafe/wrong.parquet"
    aggregate.write_text(json.dumps(manifest), encoding="utf-8")
    _write_aggregate_stage_manifest(root=root)
    with pytest.raises(InputValidationError, match="incompatible path"):
        verify_structural_resource_outputs(resource_dir=root)

    root = _structural_resource(
        root=tmp_path / "second" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    manifest = json.loads(aggregate.read_text(encoding="utf-8"))
    del manifest["datasets"]["structural_alignments"]
    aggregate.write_text(json.dumps(manifest), encoding="utf-8")
    _write_aggregate_stage_manifest(root=root)
    with pytest.raises(InputValidationError, match="not checksum-inventoried"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())

    root = _structural_resource(
        root=tmp_path / "third" / "09b_structural_alignment" / "structural_alignment"
    )
    aggregate = _convert_to_aggregate(root=root)
    manifest = json.loads(aggregate.read_text(encoding="utf-8"))
    manifest["configuration_digest"] = "bad"
    aggregate.write_text(json.dumps(manifest), encoding="utf-8")
    _write_aggregate_stage_manifest(root=root)
    with pytest.raises(InputValidationError, match="configuration_digest"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())


def test_structural_manifest_validation_rejects_tampering(tmp_path: Path) -> None:
    """Unsafe, missing, duplicated and changed outputs should be rejected."""

    root = _structural_resource(root=tmp_path / "structural")
    manifest = json.loads((root / "provenance" / "run_manifest.json").read_text(encoding="utf-8"))
    unsafe = dict(manifest)
    unsafe["outputs"] = [dict(manifest["outputs"][0], path="../escape")]
    with pytest.raises(InputValidationError, match="unsafe path"):
        verify_structural_resource_outputs(resource_dir=root, manifest=unsafe)
    duplicate = dict(manifest)
    duplicate["outputs"] = [manifest["outputs"][0], manifest["outputs"][0]]
    with pytest.raises(InputValidationError, match="repeated"):
        verify_structural_resource_outputs(resource_dir=root, manifest=duplicate)
    changed = dict(manifest)
    changed["outputs"] = [dict(manifest["outputs"][0], size_bytes=0)]
    with pytest.raises(InputValidationError, match="size mismatch"):
        verify_structural_resource_outputs(resource_dir=root, manifest=changed)
    changed["outputs"] = [dict(manifest["outputs"][0], sha256="0" * 64)]
    with pytest.raises(InputValidationError, match="checksum mismatch"):
        verify_structural_resource_outputs(resource_dir=root, manifest=changed)
    with pytest.raises(InputValidationError, match="no output inventory"):
        verify_structural_resource_outputs(
            resource_dir=root,
            manifest={"status": "complete", "outputs": []},
        )
    with pytest.raises(InputValidationError, match="mixes outputs and datasets"):
        verify_structural_resource_outputs(
            resource_dir=root,
            manifest={
                "status": "complete",
                "outputs": manifest["outputs"],
                "datasets": {"structural_alignments": {}},
            },
        )


def test_structural_import_requires_every_consumed_table_in_manifest(tmp_path: Path) -> None:
    """A present but unmanifested Parquet file must not become evidence."""

    root = _structural_resource(root=tmp_path / "structural")
    manifest_path = root / "provenance" / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"] = [
        row for row in manifest["outputs"] if row["path"] != "tables/structural_alignments.parquet"
    ]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(InputValidationError, match="not checksum-inventoried"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())


def test_structural_import_rejects_incompatible_content(tmp_path: Path) -> None:
    """Incomplete manifests, bad run digests and unmatched accessions should fail."""

    root = _structural_resource(root=tmp_path / "structural")
    manifest_path = root / "provenance" / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "failed"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputValidationError, match="not marked complete"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())
    manifest["status"] = "complete"
    manifest["run_digest"] = "bad"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputValidationError, match="run_digest"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())
    manifest["run_digest"] = "a" * 64
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputValidationError, match="No structural-alignment"):
        import_structural_alignment_resource(
            resource_dir=root,
            sequences=(SequenceRecord("other", "", "AAAA", 4, sha256_text(text="AAAA")),),
        )


def test_structural_import_rejects_uncontrolled_self_alignment(tmp_path: Path) -> None:
    """Only the precursor's explicit reference sentinel may be a diagonal row."""

    root = _structural_resource(
        root=tmp_path / "structural",
        reference_status="COMPLETE",
    )
    with pytest.raises(InputValidationError, match="without the explicit REFERENCE status"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())


def test_structural_import_rejects_reference_summary_disagreement(tmp_path: Path) -> None:
    """Reference membership must agree with the checksum-bound group summary."""

    root = _structural_resource(
        root=tmp_path / "structural",
        summary_reference_accession="p2",
    )
    with pytest.raises(InputValidationError, match="declares 'p2'"):
        import_structural_alignment_resource(resource_dir=root, sequences=_sequences())


def test_global_alignment_reader_supports_reference_only_and_pair_only_tables(
    tmp_path: Path,
) -> None:
    """Reference sentinels and older pair-only resources should remain executable."""

    path = tmp_path / "alignments.parquet"
    reference = _alignment_row()
    pq.write_table(pa.Table.from_pylist([reference]), path)
    comparisons, universes, reference_count = _read_global_alignments(
        path=path,
        lengths={"p1": 100},
        run_digest="a" * 64,
        expected_universe_sizes={"cluster_1": 1},
        expected_reference_accessions={"cluster_1": "p1"},
    )
    assert comparisons == ()
    assert tuple(universes.values()) == (frozenset({"p1"}),)
    assert reference_count == 1

    pair = _alignment_row(
        mobile_accession="p2",
        status="COMPLETE",
        aligned_length=80,
        rmsd_angstrom=1.5,
        minimum_tm_score=0.8,
    )
    pq.write_table(pa.Table.from_pylist([pair]), path)
    comparisons, universes, reference_count = _read_global_alignments(
        path=path,
        lengths={"p1": 100, "p2": 160},
        run_digest="a" * 64,
        expected_universe_sizes={"cluster_1": 2},
        expected_reference_accessions={"cluster_1": "p1"},
    )
    assert len(comparisons) == 1
    assert tuple(universes.values()) == (frozenset({"p1", "p2"}),)
    assert reference_count == 0


@pytest.mark.parametrize(
    ("rows", "expected_sizes", "expected_references", "message"),
    (
        (
            (_alignment_row(alignment_tool=""),),
            {"cluster_1": 1},
            {"cluster_1": "p1"},
            "lacks tool or cluster provenance",
        ),
        (
            (_alignment_row(),),
            {"cluster_1": 1},
            {},
            "lacks a group summary",
        ),
        (
            (_alignment_row(aligned_length=100),),
            {"cluster_1": 1},
            {"cluster_1": "p1"},
            "does not contain the expected",
        ),
        (
            (_alignment_row(), _alignment_row()),
            {"cluster_1": 1},
            {"cluster_1": "p1"},
            "Duplicate imported structural reference row",
        ),
        (
            (_alignment_row(),),
            {"cluster_1": 2},
            {"cluster_1": "p1"},
            "contains 1 campaign proteins",
        ),
    ),
)
def test_global_alignment_reference_rows_fail_closed(
    tmp_path: Path,
    rows: tuple[dict[str, object], ...],
    expected_sizes: dict[str, int],
    expected_references: dict[str, str],
    message: str,
) -> None:
    """Malformed reference sentinels must not weaken imported evidence."""

    path = tmp_path / "alignments.parquet"
    pq.write_table(pa.Table.from_pylist(list(rows)), path)
    with pytest.raises(InputValidationError, match=message):
        _read_global_alignments(
            path=path,
            lengths={"p1": 100},
            run_digest="a" * 64,
            expected_universe_sizes=expected_sizes,
            expected_reference_accessions=expected_references,
        )


def test_structural_value_helpers_cover_valid_and_invalid_states(tmp_path: Path) -> None:
    """Scalar and Parquet helpers should preserve blanks and reject malformed values."""

    assert _boolean(value=True, field="flag") is True
    assert _boolean(value="yes", field="flag") is True
    assert _boolean(value="0", field="flag") is False
    with pytest.raises(InputValidationError):
        _boolean(value=None, field="flag")
    assert _non_negative_integer(value="2", field="n") == 2
    assert _optional_integer(value="", field="n", row_number=1) is None
    assert _optional_integer(value=2, field="n", row_number=1) == 2
    assert _optional_number(value=None, field="x", row_number=1) is None
    assert _optional_number(value="1.5", field="x", row_number=1) == 1.5
    for value in (True, "bad", -1):
        with pytest.raises(InputValidationError):
            _non_negative_integer(value=value, field="n")
    with pytest.raises(InputValidationError, match="positive"):
        _optional_integer(value=0, field="n", row_number=1)
    for value in ("bad", "nan"):
        with pytest.raises(InputValidationError):
            _optional_number(value=value, field="x", row_number=1)
    assert _parse_comparison_status(value="complete", row_number=1).value == "COMPLETE"
    with pytest.raises(InputValidationError, match="Unsupported structural status"):
        _parse_comparison_status(value="invented", row_number=1)
    assert _safe_path(value="tables/a.parquet", index=0) == Path("tables/a.parquet")
    for value in (None, "", "/absolute", "../escape", "."):
        with pytest.raises(InputValidationError):
            _safe_path(value=value, index=0)
    missing = tmp_path / "missing.parquet"
    with pytest.raises(InputValidationError, match="Missing or empty"):
        _read_parquet(path=missing, required=("x",))
    wrong = tmp_path / "wrong.parquet"
    pq.write_table(pa.Table.from_pylist([{"y": 1}]), wrong)
    with pytest.raises(InputValidationError, match="lacks required"):
        _read_parquet(path=wrong, required=("x",))


def _sequences() -> tuple[SequenceRecord, ...]:
    """Return small structures-compatible protein sequences."""

    return (
        SequenceRecord("p1", "", "A" * 100, 100, sha256_text(text="A" * 100)),
        SequenceRecord("p2", "", "C" * 160, 160, sha256_text(text="C" * 160)),
        SequenceRecord("p3", "", "D" * 120, 120, sha256_text(text="D" * 120)),
    )


def _structural_resource(
    *,
    root: Path,
    reference_status: str = "REFERENCE",
    summary_reference_accession: str = "p1",
) -> Path:
    """Create a checksum-manifested structural-alignment result fixture."""

    tables = root / "tables"
    tables.mkdir(parents=True)
    alignments = [
        {
            "cluster_id": "cluster_1",
            "reference_accession": "p1",
            "mobile_accession": "p1",
            "alignment_tool": "TM-align",
            "status": reference_status,
            "tool_version": "20240303",
            "aligned_length": None,
            "rmsd_angstrom": 0.0,
            "minimum_tm_score": 1.0,
        },
        {
            "cluster_id": "cluster_1",
            "reference_accession": "p1",
            "mobile_accession": "p1",
            "alignment_tool": "US-align",
            "status": "REFERENCE",
            "tool_version": "20241201",
            "aligned_length": None,
            "rmsd_angstrom": 0.0,
            "minimum_tm_score": 1.0,
        },
        {
            "cluster_id": "cluster_1",
            "reference_accession": "p1",
            "mobile_accession": "p2",
            "alignment_tool": "TM-align",
            "status": "COMPLETE",
            "tool_version": "20240303",
            "aligned_length": 80,
            "rmsd_angstrom": 1.5,
            "minimum_tm_score": 0.8,
        },
        {
            "cluster_id": "cluster_1",
            "reference_accession": "p1",
            "mobile_accession": "p2",
            "alignment_tool": "US-align",
            "status": "COMPLETE",
            "tool_version": "20241201",
            "aligned_length": 80,
            "rmsd_angstrom": 1.6,
            "minimum_tm_score": 0.79,
        },
    ]
    pockets = [
        {
            "cluster_id": "cluster_1",
            "reference_accession": "p1",
            "mobile_accession": "p2",
            "alignment_tool": "US-align",
            "status": "ASSESSED",
            "same_pocket_position_supported": True,
            "pocket_structure_conserved": True,
        },
        {
            "cluster_id": "cluster_2",
            "reference_accession": "p1",
            "mobile_accession": "p2",
            "alignment_tool": "US-align",
            "status": "ASSESSED",
            "same_pocket_position_supported": True,
            "pocket_structure_conserved": True,
        },
    ]
    summaries = [
        {
            "cluster_id": "cluster_1",
            "primary_group_type": "HOG",
            "primary_group_id": "N0.HOG1",
            "reference_accession": summary_reference_accession,
            "alignment_tools": "TM-align;US-align",
            "alignment_tool_count": 2,
            "selected_accession_count": 3,
            "model_available_accession_count": 2,
            "aligned_accession_count": 2,
            "supported_accession_count": 2,
            "position_supported_accession_count": 2,
            "group_support_fraction": 1.0,
            "group_position_support_fraction": 1.0,
            "mean_minimum_tm_score": 0.795,
            "mean_pocket_overlap_fraction": 0.8,
            "median_centroid_distance_angstrom": 2.0,
            "position_alignment_status": "SAME_3D_POCKET_POSITION_SUPPORTED",
            "alignment_status": "CONSERVED_3D_POCKET_SUPPORTED",
            "interpretation": "computational evidence",
        }
    ]
    pq.write_table(pa.Table.from_pylist(alignments), tables / _REQUIRED_NAMES[0])
    pq.write_table(pa.Table.from_pylist(pockets), tables / _REQUIRED_NAMES[1])
    pq.write_table(pa.Table.from_pylist(summaries), tables / _REQUIRED_NAMES[2])
    output_rows = [
        {
            "path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in sorted(tables.glob("*.parquet"))
    ]
    manifest = {
        "status": "complete",
        "package_version": "0.6.0",
        "run_digest": "a" * 64,
        "outputs": output_rows,
    }
    provenance = root / "provenance"
    provenance.mkdir()
    (provenance / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root


def _convert_to_aggregate(*, root: Path) -> Path:
    """Replace a standalone fixture manifest with the Stage 09b aggregate schema."""

    manifest_path = root / "provenance" / "run_manifest.json"
    datasets = {
        name.removesuffix(".parquet"): {
            "path": str((root / "tables" / name).resolve()),
            "sha256": sha256_file(path=root / "tables" / name),
        }
        for name in _REQUIRED_NAMES
    }
    manifest_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "configuration_digest": "b" * 64,
                "finished_at_utc": "2026-09-14T00:00:00Z",
                "task_count": 1,
                "summary_group_count": 1,
                "structural_evidence_counts": {},
                "datasets": datasets,
                "shards": [],
            }
        ),
        encoding="utf-8",
    )
    _write_aggregate_stage_manifest(root=root)
    return manifest_path


def _write_aggregate_stage_manifest(*, root: Path) -> None:
    """Checksum-bind aggregate fixture datasets through the outer stage manifest."""

    stage_root = root.parent
    paths = (
        root / "provenance" / "run_manifest.json",
        *(root / "tables" / name for name in _REQUIRED_NAMES),
    )
    outputs = [
        {
            "path": str(path.relative_to(stage_root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in paths
    ]
    (stage_root / "stage_manifest.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "package_version": "0.16.0",
                "outputs": outputs,
            }
        ),
        encoding="utf-8",
    )


_REQUIRED_NAMES = (
    "structural_alignments.parquet",
    "pocket_comparisons.parquet",
    "structural_alignment_summary.parquet",
)
