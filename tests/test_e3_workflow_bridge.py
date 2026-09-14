"""Tests for the completed E3 end-to-end workflow review bridge."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import protein_signatures.e3_workflow_bridge as bridge_module
from protein_signatures.checksums import sha256_file
from protein_signatures.e3_workflow_bridge import (
    _prepare_domains,
    _prepare_sequences,
    _prepare_structures,
    _read_parquet_records,
    _require_complete_manifest,
    _resolve_asset_path,
    _resolve_first_file,
    _verify_manifested_files,
    prepare_e3_workflow_inputs,
    resolve_e3_workflow_paths,
)
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.fasta import read_protein_fasta
from protein_signatures.io_utils import iter_tsv
from protein_signatures.tables import read_domains, read_structures


def _write_parquet(*, path: Path, rows: list[dict[str, object]], schema: pa.Schema) -> None:
    """Write one deterministic Parquet fixture with an explicit schema."""

    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)


def _write_stage_manifest(*, stage_root: Path, paths: tuple[Path, ...]) -> None:
    """Write a checksum inventory for selected stage fixture outputs."""

    records = [
        {
            "path": str(path.relative_to(stage_root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in paths
    ]
    (stage_root / "stage_manifest.json").write_text(
        json.dumps({"status": "complete", "outputs": records}),
        encoding="utf-8",
    )


def _completed_workflow(tmp_path: Path) -> Path:
    """Build the supported predecessor layout with two sequence-matched models."""

    root = tmp_path / "e3_run"
    final = root / "11_app_ready"
    final.mkdir(parents=True)
    (final / "stage_manifest.json").write_text(json.dumps({"status": "complete"}), encoding="utf-8")

    sequences = root / "05_orthology" / "orthology" / "tables"
    sequence_path = sequences / "candidate_group_member_sequences.parquet"
    sequence_schema = pa.schema(
        [
            ("cluster_id", pa.string()),
            ("group_id", pa.string()),
            ("orthogroup_id", pa.string()),
            ("species", pa.string()),
            ("parsed_accession", pa.string()),
            ("is_input_candidate", pa.bool_()),
            ("protein_sequence", pa.string()),
            ("sequence_length", pa.int64()),
            ("sequence_sha256", pa.string()),
        ]
    )
    sequence_rows: list[dict[str, object]] = []
    for cluster, group, protein_id, sequence, candidate in (
        ("cluster_1", "N0.HOG0001", "P00001", "MACDEFGH", True),
        ("cluster_2", "N0.HOG0002", "Q00002", "MAAAACCC", False),
        ("cluster_3", "N0.HOG0003", "P00001", "MACDEFGH", False),
    ):
        sequence_rows.append(
            {
                "cluster_id": cluster,
                "group_id": group,
                "orthogroup_id": "OG0001",
                "species": "Arabidopsis_thaliana",
                "parsed_accession": protein_id,
                "is_input_candidate": candidate,
                "protein_sequence": sequence,
                "sequence_length": len(sequence),
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
        )
    _write_parquet(path=sequence_path, rows=sequence_rows, schema=sequence_schema)
    _write_stage_manifest(stage_root=root / "05_orthology", paths=(sequence_path,))

    domain_root = root / "06_domains"
    hit_path = domain_root / "tables" / "domain_hits.parquet"
    hit_schema = pa.schema(
        [
            ("member_accession", pa.string()),
            ("source_database", pa.string()),
            ("entry_accession", pa.string()),
            ("entry_name", pa.string()),
            ("location_start", pa.int64()),
            ("location_end", pa.int64()),
            ("score", pa.float64()),
            ("e3_family", pa.string()),
            ("evidence_role", pa.string()),
        ]
    )
    _write_parquet(
        path=hit_path,
        rows=[
            {
                "member_accession": "P00001",
                "source_database": "Pfam",
                "entry_accession": "PF00646",
                "entry_name": "F-box domain",
                "location_start": 2,
                "location_end": 6,
                "score": 21.5,
                "e3_family": "F-box",
                "evidence_role": "SUBSTRATE_RECEPTOR",
            }
        ],
        schema=hit_schema,
    )
    summary_path = domain_root / "tables" / "domain_summary.parquet"
    summary_schema = pa.schema(
        [
            ("member_accession", pa.string()),
            ("annotation_availability_status", pa.string()),
            ("annotation_status_detail", pa.string()),
            ("pfam_hit_count", pa.int64()),
            ("interpro_version", pa.string()),
        ]
    )
    _write_parquet(
        path=summary_path,
        rows=[
            {
                "member_accession": "P00001",
                "annotation_availability_status": "AVAILABLE",
                "annotation_status_detail": "retrieved",
                "pfam_hit_count": 1,
                "interpro_version": "105.0",
            },
            {
                "member_accession": "Q00002",
                "annotation_availability_status": "AVAILABLE",
                "annotation_status_detail": "retrieved no Pfam hit",
                "pfam_hit_count": 0,
                "interpro_version": "105.0",
            },
        ],
        schema=summary_schema,
    )
    _write_stage_manifest(
        stage_root=domain_root,
        paths=(hit_path, summary_path),
    )

    asset_root = root / "09_ligandability"
    asset_rows: list[dict[str, object]] = []
    for protein_id in ("P00001", "Q00002"):
        coordinate = asset_root / "models" / protein_id / f"{protein_id}.cif"
        coordinate.parent.mkdir(parents=True, exist_ok=True)
        coordinate.write_text(f"data_{protein_id}\n", encoding="utf-8")
        asset_rows.append(
            {
                "accession": protein_id,
                "action": "REUSED",
                "bytes": coordinate.stat().st_size,
                "path": str(coordinate.relative_to(asset_root)),
                "sha256": sha256_file(path=coordinate),
                "url": f"https://example.invalid/{protein_id}.cif",
            }
        )
    asset_path = asset_root / "tables" / "reused_asset_manifest.parquet"
    _write_parquet(
        path=asset_path,
        rows=asset_rows,
        schema=pa.schema(
            [
                ("accession", pa.string()),
                ("action", pa.string()),
                ("bytes", pa.int64()),
                ("path", pa.string()),
                ("sha256", pa.string()),
                ("url", pa.string()),
            ]
        ),
    )
    quality_path = asset_root / "tables" / "reused_model_quality.parquet"
    _write_parquet(
        path=quality_path,
        rows=[
            {"accession": "P00001", "mean_plddt": 90.0},
            {"accession": "Q00002", "mean_plddt": 88.0},
        ],
        schema=pa.schema(
            [
                ("accession", pa.string()),
                ("mean_plddt", pa.float64()),
            ]
        ),
    )
    _write_stage_manifest(stage_root=asset_root, paths=(asset_path, quality_path))

    structural = root / "09b_structural_alignment" / "structural_alignment"
    structural_tables = []
    datasets = {}
    for name in (
        "structural_alignments.parquet",
        "pocket_comparisons.parquet",
        "structural_alignment_summary.parquet",
    ):
        path = structural / "tables" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"fixture {name}\n", encoding="utf-8")
        structural_tables.append(path)
        datasets[name.removesuffix(".parquet")] = {
            "path": str(path.resolve()),
            "sha256": sha256_file(path=path),
        }
    provenance = structural / "provenance"
    provenance.mkdir(parents=True)
    structural_manifest = provenance / "run_manifest.json"
    structural_manifest.write_text(
        json.dumps(
            {
                "status": "complete",
                "configuration_digest": "a" * 64,
                "task_count": 2,
                "summary_group_count": 2,
                "structural_evidence_counts": {},
                "datasets": datasets,
                "shards": [],
            }
        ),
        encoding="utf-8",
    )
    structural_stage = root / "09b_structural_alignment"
    _write_stage_manifest(
        stage_root=structural_stage,
        paths=(structural_manifest, *structural_tables),
    )
    orthofinder = root / "04_orthofinder" / "Results"
    orthofinder.mkdir(parents=True)
    (orthofinder / "Log.txt").write_text(
        "OrthoFinder version 3.0.1\nOrthoFinder run completed\n",
        encoding="utf-8",
    )
    hog = orthofinder / "Phylogenetic_Hierarchical_Orthogroups" / "N0.tsv"
    hog.parent.mkdir()
    hog.write_text(
        "HOG\tOG\tGene Tree Parent Clade\tArabidopsis_thaliana\n"
        "N0.HOG0001\tOG0001\tn0\tP00001, Q00002\n",
        encoding="utf-8",
    )
    return root


def test_completed_workflow_preparation_is_conservative_and_executable(
    tmp_path: Path,
) -> None:
    """The bridge should publish valid inputs without promoting family hints."""

    root = _completed_workflow(tmp_path)
    resolved = resolve_e3_workflow_paths(run_root=root)
    assert resolved.orthofinder_results == root / "04_orthofinder" / "Results"
    assert resolved.structural_stage_manifest == (
        root / "09b_structural_alignment" / "stage_manifest.json"
    )
    destination = prepare_e3_workflow_inputs(
        run_root=root,
        output_dir=tmp_path / "prepared",
    )
    marker = json.loads((destination / "PREPARED.json").read_text(encoding="utf-8"))
    assert marker["protein_count"] == 2
    assert marker["structure_count"] == 2
    assert marker["foldseek_eligible_structure_count"] == 2
    assert marker["minimum_mean_plddt"] == 50.0
    assert marker["training_eligible_assignment_count"] == 0
    assert marker["next_action"] == "CURATE_LABEL_ASSIGNMENTS_AND_CONTROLS"

    sequences = read_protein_fasta(path=destination / "proteins.faa")
    assert [record.protein_id for record in sequences] == ["P00001", "Q00002"]
    domain_hits, domain_assessments = read_domains(
        path=destination / "domains.tsv",
        sequences=sequences,
    )
    assert domain_hits[0].domain_id == "PF00646"
    assert [item.assessment_status.value for item in domain_assessments] == [
        "ASSESSED_WITH_HIT",
        "ASSESSED_NO_HIT",
    ]
    structures = read_structures(path=destination / "structures.tsv", sequences=sequences)
    assert len(structures) == 2
    assert all(item.is_coordinate_analysis_eligible for item in structures)
    assert {item.protein_id: item.mean_confidence for item in structures} == {
        "P00001": 90.0,
        "Q00002": 88.0,
    }
    assert {item.structure_version for item in structures} == {""}
    assert all("publication_action=REUSED" in item.structure_source for item in structures)
    labels = tuple(
        iter_tsv(
            path=destination / "label_assignments.REVIEW_REQUIRED.tsv",
            required_fields=("protein_id", "curation_status"),
        )
    )
    assert {row["curation_status"] for row in labels} == {"UNMAPPED"}
    review = tuple(
        iter_tsv(
            path=destination / "e3_label_curation_review.tsv",
            required_fields=("protein_id", "upstream_e3_families"),
        )
    )
    assert review[0]["upstream_e3_families"] == "F-box"
    assert "cluster_1|cluster_3" in review[0]["cluster_ids"]
    assert review[0]["input_candidate_states"] == "FALSE|TRUE"
    source_inventory = tuple(
        iter_tsv(
            path=destination / "source_inventory.tsv",
            required_fields=("authority", "sha256"),
        )
    )
    assert len(source_inventory) == 8
    assert {row["authority"] for row in source_inventory} >= {
        "structural_run_manifest",
        "structural_stage_manifest",
    }
    with pytest.raises(PublicationError, match="already exists"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=destination)


def test_workflow_resolution_rejects_incomplete_and_changed_sources(tmp_path: Path) -> None:
    """Final status and stage checksums should be hard compatibility boundaries."""

    root = _completed_workflow(tmp_path)
    final = root / "11_app_ready" / "stage_manifest.json"
    final.write_text(json.dumps({"status": "running"}), encoding="utf-8")
    with pytest.raises(InputValidationError, match="not marked complete"):
        resolve_e3_workflow_paths(run_root=root)

    root = _completed_workflow(tmp_path / "second")
    summary = root / "06_domains" / "tables" / "domain_summary.parquet"
    summary.write_bytes(summary.read_bytes() + b"changed")
    with pytest.raises(InputValidationError, match="size mismatch"):
        resolve_e3_workflow_paths(run_root=root)

    root = _completed_workflow(tmp_path / "third")
    stage_manifest = root / "09_ligandability" / "stage_manifest.json"
    stage_manifest.write_text(json.dumps({"status": "complete", "outputs": []}), encoding="utf-8")
    with pytest.raises(InputValidationError, match="no output inventory"):
        resolve_e3_workflow_paths(run_root=root)

    with pytest.raises(InputValidationError, match="does not exist"):
        resolve_e3_workflow_paths(run_root=tmp_path / "absent")

    root = _completed_workflow(tmp_path / "fourth")
    (root / "04_orthofinder" / "Results").rename(root / "04_orthofinder" / "moved")
    with pytest.raises(InputValidationError, match="OrthoFinder Results"):
        resolve_e3_workflow_paths(run_root=root)


def test_manifest_helpers_reject_unsafe_duplicate_and_invalid_records(tmp_path: Path) -> None:
    """Every consumed stage authority should have one safe valid checksum record."""

    manifest = tmp_path / "document.json"
    manifest.write_text("[]\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="contain an object"):
        _require_complete_manifest(path=manifest)

    stage = tmp_path / "05_stage"
    stage.mkdir()
    source = stage / "table.tsv"
    source.write_text("id\n1\n", encoding="utf-8")
    valid = {
        "path": "table.tsv",
        "size_bytes": source.stat().st_size,
        "sha256": sha256_file(path=source),
    }
    cases = (
        (["bad"], "Malformed output record"),
        ([{**valid, "path": "../table.tsv"}], "Unsafe output path"),
        ([valid, valid], "found 2"),
        ([{**valid, "size_bytes": "bad"}], "Invalid size_bytes"),
        ([{**valid, "sha256": "bad"}], "Invalid checksum record"),
        ([{**valid, "sha256": "0" * 64}], "checksum mismatch"),
    )
    for outputs, message in cases:
        (stage / "stage_manifest.json").write_text(
            json.dumps({"status": "complete", "outputs": outputs}),
            encoding="utf-8",
        )
        with pytest.raises(InputValidationError, match=message):
            _verify_manifested_files(stage_root=stage, paths=(source,))

    (stage / "stage_manifest.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "outputs": [{**valid, "path": "05_stage/table.tsv"}],
            }
        ),
        encoding="utf-8",
    )
    _verify_manifested_files(stage_root=stage, paths=(source,))
    with pytest.raises(InputValidationError, match="Could not find table"):
        _resolve_first_file(root=stage, candidates=("missing.tsv",), label="table")


def test_parquet_and_sequence_helpers_reject_malformed_authorities(tmp_path: Path) -> None:
    """Parquet shape and sequence identity failures should be diagnosed exactly."""

    table = tmp_path / "table.parquet"
    _write_parquet(
        path=table,
        rows=[{"id": "x"}],
        schema=pa.schema([("id", pa.string())]),
    )
    assert _read_parquet_records(path=table, required=("id",))[0]["id"] == "x"
    with pytest.raises(InputValidationError, match="non-empty and unique"):
        _read_parquet_records(path=table, required=("id", "id"))
    with pytest.raises(InputValidationError, match="missing required columns"):
        _read_parquet_records(path=table, required=("other",))
    with pytest.raises(InputValidationError, match="Missing or empty"):
        _read_parquet_records(path=tmp_path / "missing.parquet", required=("id",))
    empty = tmp_path / "empty.parquet"
    _write_parquet(
        path=empty,
        rows=[],
        schema=pa.schema([("id", pa.string())]),
    )
    assert _read_parquet_records(path=empty, required=("id",), allow_empty=True) == ()
    with pytest.raises(InputValidationError, match="contains no records"):
        _read_parquet_records(path=empty, required=("id",))
    corrupt = tmp_path / "corrupt.parquet"
    corrupt.write_text("not parquet\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Could not read Parquet"):
        _read_parquet_records(path=corrupt, required=("id",))

    sequence = "MACD"
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    base = {
        "parsed_accession": "P00001",
        "protein_sequence": sequence,
        "sequence_length": len(sequence),
        "sequence_sha256": digest,
    }
    prepared = _prepare_sequences(rows=({**base}, {**base}, {**base, "parsed_accession": ""}))
    assert prepared.skipped_unmapped_rows == 1
    with pytest.raises(InputValidationError, match="inconsistent length"):
        _prepare_sequences(rows=({**base, "sequence_length": 99},))
    with pytest.raises(InputValidationError, match="invalid sequence_length"):
        _prepare_sequences(rows=({**base, "sequence_length": "bad"},))
    with pytest.raises(InputValidationError, match="whitespace-bearing"):
        _prepare_sequences(rows=({**base, "protein_sequence": "MA CD"},))
    with pytest.raises(InputValidationError, match="Conflicting sequences"):
        _prepare_sequences(
            rows=(
                base,
                {
                    **base,
                    "protein_sequence": "MACE",
                    "sequence_sha256": hashlib.sha256(b"MACE").hexdigest(),
                },
            )
        )
    with pytest.raises(InputValidationError, match="non-ASCII"):
        _prepare_sequences(
            rows=(
                {
                    **base,
                    "protein_sequence": "MAÇD",
                    "sequence_sha256": hashlib.sha256("MAÇD".encode()).hexdigest(),
                },
            )
        )


def test_domain_helper_preserves_no_hit_unknown_and_rejects_bad_payloads() -> None:
    """Pfam conversion should preserve assessment state and coordinate integrity."""

    sequences = {"P1": "MACDE", "P2": "MAAAA", "P3": "MCCCC"}
    hit = {
        "member_accession": "P1",
        "source_database": "Pfam",
        "entry_accession": "PF1",
        "entry_name": "Domain one",
        "location_start": 2,
        "location_end": 4,
        "score": 2.5,
        "e3_family": "RING",
        "evidence_role": "CATALYTIC_E3",
    }
    summaries = (
        {
            "member_accession": "P1",
            "annotation_availability_status": "AVAILABLE",
            "annotation_status_detail": "ok",
            "pfam_hit_count": 1,
            "interpro_version": "105",
        },
        {
            "member_accession": "P2",
            "annotation_availability_status": "AVAILABLE",
            "annotation_status_detail": "no hit",
            "pfam_hit_count": 0,
            "interpro_version": "105",
        },
    )
    prepared = _prepare_domains(
        hit_rows=(hit,),
        summary_rows=summaries,
        sequences=sequences,
        evidence_reference="sha256-bound",
    )
    assert [row["assessment_status"] for row in prepared.domain_records] == [
        "ASSESSED_WITH_HIT",
        "ASSESSED_NO_HIT",
        "NOT_ASSESSED",
    ]
    assert prepared.audit_by_protein["P1"]["upstream_e3_families"] == "RING"
    with pytest.raises(InputValidationError, match="exceed sequence bounds"):
        _prepare_domains(
            hit_rows=({**hit, "location_end": 99},),
            summary_rows=summaries,
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    with pytest.raises(InputValidationError, match="Non-finite Pfam score"):
        _prepare_domains(
            hit_rows=({**hit, "score": float("nan")},),
            summary_rows=summaries,
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    with pytest.raises(InputValidationError, match="declares hits"):
        _prepare_domains(
            hit_rows=(),
            summary_rows=(summaries[0],),
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    with pytest.raises(InputValidationError, match="Invalid pfam_hit_count"):
        _prepare_domains(
            hit_rows=(),
            summary_rows=({**summaries[0], "pfam_hit_count": "bad"},),
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    with pytest.raises(InputValidationError, match="Negative pfam_hit_count"):
        _prepare_domains(
            hit_rows=(),
            summary_rows=({**summaries[0], "pfam_hit_count": -1},),
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    with pytest.raises(InputValidationError, match="Invalid Pfam coordinates"):
        _prepare_domains(
            hit_rows=({**hit, "location_start": "bad"},),
            summary_rows=summaries,
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    with pytest.raises(InputValidationError, match="Invalid Pfam score"):
        _prepare_domains(
            hit_rows=({**hit, "score": "bad"},),
            summary_rows=summaries,
            sequences=sequences,
            evidence_reference="sha256-bound",
        )
    non_pfam = _prepare_domains(
        hit_rows=({**hit, "source_database": "InterPro"},),
        summary_rows=({**summaries[0], "pfam_hit_count": 0}, summaries[1]),
        sequences=sequences,
        evidence_reference="sha256-bound",
    )
    assert non_pfam.audit_by_protein["P1"]["upstream_e3_families"] == "RING"
    assert non_pfam.domain_records[0]["assessment_status"] == "ASSESSED_NO_HIT"


def test_structure_helpers_verify_paths_checksums_sizes_and_duplicates(tmp_path: Path) -> None:
    """Only unique checksum-matching coordinate assets should be emitted."""

    run_root = tmp_path / "run"
    manifest = run_root / "09_ligandability" / "tables" / "assets.parquet"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("fixture\n", encoding="utf-8")
    coordinate = run_root / "09_ligandability" / "models" / "P1.cif"
    coordinate.parent.mkdir(parents=True)
    coordinate.write_text("data_P1\n", encoding="utf-8")
    digest = sha256_file(path=coordinate)
    row = {
        "accession": "P1",
        "action": "REUSED",
        "bytes": coordinate.stat().st_size,
        "path": "models/P1.cif",
        "sha256": digest,
    }
    assert (
        _resolve_asset_path(
            value=str(coordinate),
            run_root=run_root,
            asset_manifest_path=manifest,
        )
        == coordinate
    )
    prepared = _prepare_structures(
        rows=(row, {**row}, {**row, "accession": "P2"}, {**row, "path": "note.json"}),
        quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
        sequence_ids=frozenset({"P1"}),
        run_root=run_root,
        asset_manifest_path=manifest,
        minimum_mean_plddt=50.0,
    )
    assert len(prepared.structure_records) == 1
    assert prepared.skipped_unmatched_rows == 1
    assert prepared.skipped_non_coordinate_rows == 1
    with pytest.raises(InputValidationError, match="checksum mismatch"):
        _prepare_structures(
            rows=({**row, "sha256": "0" * 64},),
            quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt=50.0,
        )
    with pytest.raises(InputValidationError, match="Invalid coordinate SHA-256"):
        _prepare_structures(
            rows=({**row, "sha256": "bad"},),
            quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt=50.0,
        )
    with pytest.raises(InputValidationError, match="Invalid coordinate byte count"):
        _prepare_structures(
            rows=({**row, "bytes": "bad"},),
            quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt=50.0,
        )
    with pytest.raises(InputValidationError, match="byte-count mismatch"):
        _prepare_structures(
            rows=({**row, "bytes": 999},),
            quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt=50.0,
        )
    second = coordinate.with_name("P1_second.cif")
    second.write_text("different\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Multiple different"):
        _prepare_structures(
            rows=(
                row,
                {
                    **row,
                    "path": str(second),
                    "bytes": second.stat().st_size,
                    "sha256": sha256_file(path=second),
                },
            ),
            quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt=50.0,
        )
    with pytest.raises(InputValidationError, match="found 0"):
        _resolve_asset_path(
            value="missing.cif",
            run_root=run_root,
            asset_manifest_path=manifest,
        )

    preferred = coordinate.with_name("A.cif")
    preferred.write_bytes(coordinate.read_bytes())
    selected = _prepare_structures(
        rows=(
            row,
            {
                **row,
                "path": str(preferred),
                "bytes": preferred.stat().st_size,
            },
        ),
        quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
        sequence_ids=frozenset({"P1"}),
        run_root=run_root,
        asset_manifest_path=manifest,
        minimum_mean_plddt=50.0,
    )
    assert selected.structure_records[0]["coordinate_path"] == str(preferred)

    low_confidence = _prepare_structures(
        rows=(row,),
        quality_rows=({"accession": "P1", "mean_plddt": 49.9},),
        sequence_ids=frozenset({"P1"}),
        run_root=run_root,
        asset_manifest_path=manifest,
        minimum_mean_plddt=50.0,
    )
    assert low_confidence.eligible_protein_ids == frozenset()
    assert (
        low_confidence.structure_records[0]["analysis_eligibility_status"]
        == "INELIGIBLE_LOW_CONFIDENCE"
    )
    unavailable = _prepare_structures(
        rows=(row,),
        quality_rows=(),
        sequence_ids=frozenset({"P1"}),
        run_root=run_root,
        asset_manifest_path=manifest,
        minimum_mean_plddt=50.0,
    )
    assert unavailable.structure_records[0]["mean_confidence"] == ""
    assert (
        unavailable.structure_records[0]["analysis_eligibility_status"]
        == "INELIGIBLE_CONFIDENCE_UNAVAILABLE"
    )
    duplicate_quality = _prepare_structures(
        rows=(row,),
        quality_rows=(
            {"accession": "P1", "mean_plddt": 90.0},
            {"accession": "P1", "mean_plddt": 90.0},
            {"accession": "P2", "mean_plddt": 10.0},
        ),
        sequence_ids=frozenset({"P1"}),
        run_root=run_root,
        asset_manifest_path=manifest,
        minimum_mean_plddt=50.0,
    )
    assert duplicate_quality.eligible_protein_ids == frozenset({"P1"})
    with pytest.raises(InputValidationError, match="Conflicting model-quality"):
        _prepare_structures(
            rows=(row,),
            quality_rows=(
                {"accession": "P1", "mean_plddt": 90.0},
                {"accession": "P1", "mean_plddt": 80.0},
            ),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt=50.0,
        )
    for bad_confidence, message in (
        ("bad", "Invalid mean_plddt"),
        (True, "Invalid mean_plddt"),
        (float("nan"), "finite and between"),
        (101.0, "finite and between"),
    ):
        with pytest.raises(InputValidationError, match=message):
            _prepare_structures(
                rows=(row,),
                quality_rows=({"accession": "P1", "mean_plddt": bad_confidence},),
                sequence_ids=frozenset({"P1"}),
                run_root=run_root,
                asset_manifest_path=manifest,
                minimum_mean_plddt=50.0,
            )
    with pytest.raises(InputValidationError, match="minimum_mean_plddt"):
        _prepare_structures(
            rows=(row,),
            quality_rows=({"accession": "P1", "mean_plddt": 90.0},),
            sequence_ids=frozenset({"P1"}),
            run_root=run_root,
            asset_manifest_path=manifest,
            minimum_mean_plddt="bad",  # type: ignore[arg-type]
        )


def test_public_preparer_rejects_empty_or_structureless_cohorts(tmp_path: Path) -> None:
    """Preparation should require mapped sequences and a searchable model cohort."""

    root = _completed_workflow(tmp_path / "empty")
    sequence_path = (
        root / "05_orthology" / "orthology" / "tables" / "candidate_group_member_sequences.parquet"
    )
    table = pq.read_table(sequence_path)
    rows = table.to_pylist()
    for row in rows:
        row["parsed_accession"] = ""
    _write_parquet(path=sequence_path, rows=rows, schema=table.schema)
    _write_stage_manifest(stage_root=root / "05_orthology", paths=(sequence_path,))
    with pytest.raises(InputValidationError, match="No exact accession-bearing"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "no_sequences")

    root = _completed_workflow(tmp_path / "one_model")
    asset_path = root / "09_ligandability" / "tables" / "reused_asset_manifest.parquet"
    quality_path = root / "09_ligandability" / "tables" / "reused_model_quality.parquet"
    table = pq.read_table(asset_path)
    _write_parquet(path=asset_path, rows=table.to_pylist()[:1], schema=table.schema)
    _write_stage_manifest(
        stage_root=root / "09_ligandability",
        paths=(asset_path, quality_path),
    )
    with pytest.raises(InputValidationError, match="At least two"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "one_structure")

    root = _completed_workflow(tmp_path / "low_confidence")
    quality_path = root / "09_ligandability" / "tables" / "reused_model_quality.parquet"
    quality_table = pq.read_table(quality_path)
    low_rows = [{**row, "mean_plddt": 49.0} for row in quality_table.to_pylist()]
    _write_parquet(path=quality_path, rows=low_rows, schema=quality_table.schema)
    asset_path = root / "09_ligandability" / "tables" / "reused_asset_manifest.parquet"
    _write_stage_manifest(
        stage_root=root / "09_ligandability",
        paths=(asset_path, quality_path),
    )
    with pytest.raises(InputValidationError, match="confidence-eligible"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "low_structures")

    for invalid_threshold in (True, "bad", float("inf"), -1.0, 101.0):
        with pytest.raises(InputValidationError, match="minimum_mean_plddt"):
            prepare_e3_workflow_inputs(
                run_root=root,
                output_dir=tmp_path / f"invalid_{str(invalid_threshold).replace('/', '_')}",
                minimum_mean_plddt=invalid_threshold,  # type: ignore[arg-type]
            )


def test_publication_failures_are_cleaned_and_contextualised(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Failed bundle publication should remove staging directories and retain context."""

    root = _completed_workflow(tmp_path / "publication")
    destination = tmp_path / "failed_publication"
    monkeypatch.setattr(bridge_module, "read_protein_fasta", lambda **_kwargs: ())
    with pytest.raises(PublicationError, match="record count changed"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=destination)
    assert not destination.exists()
    assert not tuple(tmp_path.glob(".failed_publication.staging.*"))

    monkeypatch.undo()
    root = _completed_workflow(tmp_path / "os_error")
    monkeypatch.setattr(
        bridge_module,
        "_output_inventory",
        lambda **_kwargs: (_ for _ in ()).throw(OSError("fixture failure")),
    )
    with pytest.raises(PublicationError, match="fixture failure"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "os_failed")
