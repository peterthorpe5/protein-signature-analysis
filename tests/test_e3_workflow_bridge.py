"""Tests for the completed E3 end-to-end workflow review bridge."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

import protein_signatures.e3_workflow_bridge as bridge_module
import protein_signatures.e3_workflow_orchestration as orchestration_module
from protein_signatures.checksums import sha256_file
from protein_signatures.e3_workflow_bridge import (
    _iter_parquet_records,
    _prepare_domains,
    _prepare_sequences,
    _prepare_structures,
    _read_parquet_records,
    _require_complete_manifest,
    _resolve_asset_path,
    _resolve_first_file,
    _serialise_review_values,
    _verify_manifested_files,
    _write_prepared_fasta,
    prepare_e3_workflow_inputs,
    resolve_e3_workflow_paths,
)
from protein_signatures.e3_workflow_orchestration import (
    _copy_file_atomic,
    _normalise_profile_source,
    _require_nonnegative_integer,
    _require_sha256,
    _review_context,
    _validate_existing_e3_campaign,
    _verify_source_inventory,
    approve_e3_label_review,
    ensure_e3_campaign_marker,
    ensure_e3_preparation_marker,
    stage_e3_label_review,
    verify_e3_label_review,
    verify_prepared_e3_bundle,
)
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.fasta import read_protein_fasta
from protein_signatures.io_utils import iter_tsv, write_tsv_atomic
from protein_signatures.tables import LABEL_FIELDS, read_domains, read_structures


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
        json.dumps(
            {
                "status": "complete",
                "configuration_digest": "c" * 64,
                "package_version": "0.16.0",
                "outputs": records,
            }
        ),
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
        (
            "Arabidopsis_thaliana@@sp|B3H578|PHD1_ARATH",
            "N0.HOG0001",
            "P00001",
            "MACDEFGH",
            True,
        ),
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
    working = orthofinder / "WorkingDirectory"
    working.mkdir()
    species_ids = working / "SpeciesIDs.txt"
    species_ids.write_text("0: Arabidopsis_thaliana.faa\n", encoding="utf-8")
    sequence_ids = working / "SequenceIDs.txt"
    sequence_ids.write_text("0_0: P00001\n0_1: Q00002\n", encoding="utf-8")
    orthogroups = orthofinder / "Orthogroups" / "Orthogroups.tsv"
    orthogroups.parent.mkdir()
    orthogroups.write_text(
        "Orthogroup\tArabidopsis_thaliana\nOG0001\tP00001, Q00002\n",
        encoding="utf-8",
    )
    hog = orthofinder / "Phylogenetic_Hierarchical_Orthogroups" / "N0.tsv"
    hog.parent.mkdir()
    hog.write_text(
        "HOG\tOG\tGene Tree Parent Clade\tArabidopsis_thaliana\n"
        "N0.HOG0001\tOG0001\tn0\tP00001, Q00002\n",
        encoding="utf-8",
    )
    species_tree = orthofinder / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt"
    species_tree.parent.mkdir()
    species_tree.write_text("(Arabidopsis_thaliana)N0;\n", encoding="utf-8")
    orthofinder_stage = root / "04_orthofinder"
    authority = orthofinder_stage / "orthofinder_authority.tsv"
    authority.write_text(
        "mode\tarchive_path\tarchive_size_bytes\tarchive_sha256\t"
        "published_results\torthofinder_version\tdecision_basis\n"
        "reused_reviewed_archive\t/archive/Results_Feb26.tar.gz\t123\t"
        + "d" * 64
        + "\tResults\t2.5.5\tproject-reviewed Results_Feb26 phylogeny\n",
        encoding="utf-8",
    )
    validation = orthofinder_stage / "orthofinder_reuse_validation.tsv"
    required = (species_ids, sequence_ids, orthogroups, hog, species_tree)
    validation.write_text(
        "relative_path\tsize_bytes\tsha256\tstatus\n"
        + "".join(
            f"{path.relative_to(orthofinder)}\t{path.stat().st_size}\t"
            f"{sha256_file(path=path)}\tVALID\n"
            for path in required
        ),
        encoding="utf-8",
    )
    _write_stage_manifest(
        stage_root=orthofinder_stage,
        paths=(*required, authority, validation),
    )
    return root


def _curate_fixture_labels(*, path: Path, include_control: bool = True) -> None:
    """Replace fixture placeholders with one target and optional control assignment."""

    records = list(iter_tsv(path=path, required_fields=LABEL_FIELDS))
    records[0].update(
        {
            "label_id": "e3:ubiquitin:crl:crl1_scf:f_box",
            "curation_status": "REVIEWED_POSITIVE",
            "evidence_status": "REVIEWED",
            "component_role": "SUBSTRATE_RECEPTOR",
            "curation_reason": "Fixture target reviewed against source evidence.",
        }
    )
    if include_control:
        records[1].update(
            {
                "label_id": "control:matched_substrate_receptor_reference",
                "curation_status": "REVIEWED_POSITIVE",
                "evidence_status": "REVIEWED",
                "component_role": "CONTROL",
                "curation_reason": "Fixture control reviewed against source evidence.",
            }
        )
    write_tsv_atomic(path=path, fieldnames=LABEL_FIELDS, records=records)


def _approved_workflow(tmp_path: Path) -> dict[str, Path]:
    """Create a fully linked preparation, review and approval fixture."""

    root = _completed_workflow(tmp_path)
    work = tmp_path / "campaign"
    state = work / "workflow_state" / "e3"
    paths = {
        "root": root,
        "work": work,
        "prepared": work / "prepared_inputs",
        "preparation": state / "01_preparation" / "PREPARED_VERIFIED.json",
        "review": state / "02_label_review" / "REVIEW_READY.json",
        "approval": state / "02_label_review" / "REVIEW_APPROVED.json",
        "verification": state / "02_label_review" / "REVIEW_VERIFIED.json",
        "reviewed": work / "reviewed_label_assignments.tsv",
        "config": work / "campaign.yaml",
        "campaign_marker": state / "03_initialisation" / "CAMPAIGN_VALIDATED.json",
    }
    ensure_e3_preparation_marker(
        run_root=paths["root"],
        prepared_dir=paths["prepared"],
        minimum_mean_plddt=50,
        marker_path=paths["preparation"],
    )
    stage_e3_label_review(
        preparation_marker=paths["preparation"],
        reviewed_labels=paths["reviewed"],
        marker_path=paths["review"],
    )
    _curate_fixture_labels(path=paths["reviewed"])
    approve_e3_label_review(
        preparation_marker=paths["preparation"],
        review_marker=paths["review"],
        reviewed_labels=paths["reviewed"],
        approval_marker=paths["approval"],
        curator="Test Curator",
    )
    verify_e3_label_review(
        preparation_marker=paths["preparation"],
        review_marker=paths["review"],
        reviewed_labels=paths["reviewed"],
        approval_marker=paths["approval"],
        marker_path=paths["verification"],
    )
    return paths


def test_completed_workflow_preparation_is_conservative_and_executable(
    tmp_path: Path,
) -> None:
    """The bridge should publish valid inputs without promoting family hints."""

    root = _completed_workflow(tmp_path)
    resolved = resolve_e3_workflow_paths(run_root=root)
    assert resolved.orthofinder_results == root / "04_orthofinder" / "Results"
    assert resolved.orthofinder_version == "2.5.5"
    assert resolved.orthofinder_source_mode == "CHECKSUMMED_WORKFLOW_STAGE"
    assert resolved.orthofinder_log_path is None
    assert len(resolved.orthofinder_completion_authority_paths) == 3
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
    assert marker["orthofinder_version"] == "2.5.5"
    assert marker["orthofinder_source_mode"] == "CHECKSUMMED_WORKFLOW_STAGE"
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
    assert json.loads(review[0]["cluster_ids"]) == [
        "Arabidopsis_thaliana@@sp|B3H578|PHD1_ARATH",
        "cluster_3",
    ]
    assert review[0]["input_candidate_states"] == "FALSE|TRUE"
    source_inventory = tuple(
        iter_tsv(
            path=destination / "source_inventory.tsv",
            required_fields=("authority", "sha256"),
        )
    )
    assert len(source_inventory) == 11
    assert {row["authority"] for row in source_inventory} >= {
        "orthofinder_authority",
        "orthofinder_reuse_validation",
        "orthofinder_stage_manifest",
        "structural_run_manifest",
        "structural_stage_manifest",
    }
    with pytest.raises(PublicationError, match="already exists"):
        prepare_e3_workflow_inputs(run_root=root, output_dir=destination)


def test_e3_orchestration_happy_path_is_review_gated_and_resumable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every automated boundary should bind the authorities it consumes."""

    root = _completed_workflow(tmp_path)
    work = tmp_path / "campaign"
    prepared = work / "prepared_inputs"
    state = work / "workflow_state" / "e3"
    preparation_marker = state / "01_preparation" / "PREPARED_VERIFIED.json"
    review_marker = state / "02_label_review" / "REVIEW_READY.json"
    approval_marker = state / "02_label_review" / "REVIEW_APPROVED.json"
    verification_marker = state / "02_label_review" / "REVIEW_VERIFIED.json"
    campaign_marker = state / "03_initialisation" / "CAMPAIGN_VALIDATED.json"
    reviewed = work / "reviewed_label_assignments.tsv"

    assert (
        ensure_e3_preparation_marker(
            run_root=root,
            prepared_dir=prepared,
            minimum_mean_plddt=50.0,
            marker_path=preparation_marker,
        )
        == preparation_marker
    )
    preparation = json.loads(preparation_marker.read_text(encoding="utf-8"))
    assert preparation["preparation_action"] == "CREATED"
    ensure_e3_preparation_marker(
        run_root=root,
        prepared_dir=prepared,
        minimum_mean_plddt=50.0,
        marker_path=preparation_marker,
    )
    assert (
        json.loads(preparation_marker.read_text(encoding="utf-8"))["preparation_action"]
        == "REUSED_VERIFIED"
    )

    assert (
        stage_e3_label_review(
            preparation_marker=preparation_marker,
            reviewed_labels=reviewed,
            marker_path=review_marker,
        )
        == review_marker
    )
    assert (
        reviewed.read_bytes() == (prepared / "label_assignments.REVIEW_REQUIRED.tsv").read_bytes()
    )
    assert json.loads(review_marker.read_text(encoding="utf-8"))["staging_action"] == (
        "CREATED_FROM_TEMPLATE"
    )
    _curate_fixture_labels(path=reviewed)

    assert (
        approve_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            curator="Test Curator",
            note="Fixture approval",
        )
        == approval_marker
    )
    approval = json.loads(approval_marker.read_text(encoding="utf-8"))
    assert approval["review_status"] == "APPROVED"
    assert approval["target_reviewed_positive_count"] == 1
    assert approval["control_reviewed_positive_count"] == 1
    assert approval["reviewed_labels_sha256"] == sha256_file(path=reviewed)
    assert (
        verify_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            marker_path=verification_marker,
        )
        == verification_marker
    )

    monkeypatch.setattr(
        orchestration_module,
        "validate_campaign",
        lambda **_kwargs: {"status": "VALID", "protein_count": 2},
    )
    campaign_config = work / "campaign.yaml"
    assert (
        ensure_e3_campaign_marker(
            preparation_marker=preparation_marker,
            review_verification_marker=verification_marker,
            campaign_config=campaign_config,
            campaign_id="fixture_e3_campaign",
            profile="e3",
            marker_path=campaign_marker,
        )
        == campaign_marker
    )
    campaign = json.loads(campaign_marker.read_text(encoding="utf-8"))
    assert campaign["campaign_action"] == "CREATED"
    assert campaign["validation_summary"]["protein_count"] == 2
    assert "maximum_hits: 2" in campaign_config.read_text(encoding="utf-8")

    ensure_e3_campaign_marker(
        preparation_marker=preparation_marker,
        review_verification_marker=verification_marker,
        campaign_config=campaign_config,
        campaign_id="fixture_e3_campaign",
        profile="e3",
        marker_path=campaign_marker,
    )
    assert json.loads(campaign_marker.read_text(encoding="utf-8"))["campaign_action"] == (
        "ADOPTED_EXISTING"
    )


def test_e3_review_staging_adopts_but_never_overwrites_existing_work(
    tmp_path: Path,
) -> None:
    """A manually copied or partly curated review file must remain untouched."""

    root = _completed_workflow(tmp_path)
    prepared = tmp_path / "work" / "prepared_inputs"
    preparation_marker = tmp_path / "state" / "PREPARED_VERIFIED.json"
    ensure_e3_preparation_marker(
        run_root=root,
        prepared_dir=prepared,
        minimum_mean_plddt=50.0,
        marker_path=preparation_marker,
    )
    reviewed = tmp_path / "work" / "reviewed.tsv"
    reviewed.parent.mkdir(parents=True, exist_ok=True)
    reviewed.write_text("curator work must survive\n", encoding="utf-8")
    before = reviewed.read_bytes()
    marker = tmp_path / "state" / "REVIEW_READY.json"
    with pytest.raises(InputValidationError, match="required fields"):
        stage_e3_label_review(
            preparation_marker=preparation_marker,
            reviewed_labels=reviewed,
            marker_path=marker,
        )
    assert reviewed.read_bytes() == before

    reviewed.write_bytes((prepared / "label_assignments.REVIEW_REQUIRED.tsv").read_bytes())
    _curate_fixture_labels(path=reviewed)
    before = reviewed.read_bytes()
    stage_e3_label_review(
        preparation_marker=preparation_marker,
        reviewed_labels=reviewed,
        marker_path=marker,
    )
    assert reviewed.read_bytes() == before
    assert json.loads(marker.read_text(encoding="utf-8"))["staging_action"] == ("ADOPTED_EXISTING")


def test_e3_review_approval_rejects_placeholder_incomplete_and_changed_labels(
    tmp_path: Path,
) -> None:
    """Approval must require curated classes and remain invalid after any edit."""

    root = _completed_workflow(tmp_path)
    prepared = tmp_path / "work" / "prepared_inputs"
    preparation_marker = tmp_path / "state" / "PREPARED_VERIFIED.json"
    review_marker = tmp_path / "state" / "REVIEW_READY.json"
    approval_marker = tmp_path / "state" / "REVIEW_APPROVED.json"
    reviewed = tmp_path / "work" / "reviewed.tsv"
    ensure_e3_preparation_marker(
        run_root=root,
        prepared_dir=prepared,
        minimum_mean_plddt=50.0,
        marker_path=preparation_marker,
    )
    stage_e3_label_review(
        preparation_marker=preparation_marker,
        reviewed_labels=reviewed,
        marker_path=review_marker,
    )
    with pytest.raises(InputValidationError, match="byte-identical"):
        approve_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            curator="Test Curator",
        )

    _curate_fixture_labels(path=reviewed, include_control=False)
    with pytest.raises(InputValidationError, match="reviewed-positive control"):
        approve_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            curator="Test Curator",
        )

    records = list(iter_tsv(path=reviewed, required_fields=LABEL_FIELDS))
    records[1].update(
        {
            "label_id": "control:prespecified_non_e3_reference",
            "curation_status": "REVIEWED_POSITIVE",
            "evidence_status": "REVIEWED",
            "component_role": "CONTROL",
            "curation_reason": "Reviewed but not target-matched fixture control.",
        }
    )
    write_tsv_atomic(path=reviewed, fieldnames=LABEL_FIELDS, records=records)
    with pytest.raises(InputValidationError, match="profile-resolved matched background"):
        approve_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            curator="Test Curator",
        )

    _curate_fixture_labels(path=reviewed)
    approve_e3_label_review(
        preparation_marker=preparation_marker,
        review_marker=review_marker,
        reviewed_labels=reviewed,
        approval_marker=approval_marker,
        curator="Test Curator",
    )
    with pytest.raises(PublicationError, match="already exists"):
        approve_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            curator="Test Curator",
        )
    reviewed.write_bytes(reviewed.read_bytes() + b"\n")
    with pytest.raises(InputValidationError, match="no longer matches"):
        verify_e3_label_review(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
            approval_marker=approval_marker,
            marker_path=tmp_path / "state" / "REVIEW_VERIFIED.json",
        )


def test_e3_prepared_bundle_and_scalar_helpers_reject_tampering(
    tmp_path: Path,
) -> None:
    """Prepared outputs and marker scalars should fail closed on mutation."""

    root = _completed_workflow(tmp_path)
    prepared = prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "prepared")
    assert (
        verify_prepared_e3_bundle(
            prepared_dir=prepared,
            run_root=root,
            minimum_mean_plddt=50.0,
        )["protein_count"]
        == 2
    )
    (prepared / "domains.tsv").write_bytes((prepared / "domains.tsv").read_bytes() + b"changed")
    with pytest.raises(InputValidationError, match="size differs"):
        verify_prepared_e3_bundle(
            prepared_dir=prepared,
            run_root=root,
            minimum_mean_plddt=50.0,
        )

    assert _require_nonnegative_integer(value="0", field_name="count") == 0
    assert _require_sha256(value="a" * 64, field_name="digest") == "a" * 64
    for value in (True, -1, "1.5", None):
        with pytest.raises(InputValidationError, match="non-negative integer"):
            _require_nonnegative_integer(value=value, field_name="count")
    with pytest.raises(InputValidationError, match="lower-case SHA-256"):
        _require_sha256(value="A" * 64, field_name="digest")
    assert _normalise_profile_source(profile="e3") == "e3"
    custom_profile = tmp_path / "custom.yaml"
    custom_profile.write_text("profile_id: custom\n", encoding="utf-8")
    assert _normalise_profile_source(profile=custom_profile) == str(custom_profile.resolve())
    with pytest.raises(InputValidationError, match="missing or empty"):
        _normalise_profile_source(profile=tmp_path / "missing.yaml")


@pytest.mark.parametrize(
    ("case", "message"),
    (
        ("not_object", "contain an object"),
        ("incomplete", "not complete"),
        ("wrong_root", "source root differs"),
        ("nonnumeric_threshold", "not numeric"),
        ("different_threshold", "threshold differs"),
        ("bad_outputs", "must be a list"),
        ("malformed_output", "Malformed prepared output"),
        ("unsafe_output", "Unsafe prepared output"),
        ("duplicate_output", "Duplicate prepared output"),
        ("wrong_checksum", "checksum differs"),
        ("wrong_inventory", "inventory differs"),
    ),
)
def test_prepared_bundle_rejects_malformed_manifest_contracts(
    tmp_path: Path,
    case: str,
    message: str,
) -> None:
    """Each preparation-manifest boundary should fail with exact context."""

    root = _completed_workflow(tmp_path)
    prepared = prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "prepared")
    marker_path = prepared / "PREPARED.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if case == "not_object":
        marker = []
    elif case == "incomplete":
        marker["status"] = "RUNNING"
    elif case == "wrong_root":
        marker["source_run_root"] = str(tmp_path / "different")
    elif case == "nonnumeric_threshold":
        marker["minimum_mean_plddt"] = "bad"
    elif case == "different_threshold":
        marker["minimum_mean_plddt"] = 51
    elif case == "bad_outputs":
        marker["outputs"] = {}
    elif case == "malformed_output":
        marker["outputs"][0] = "bad"
    elif case == "unsafe_output":
        marker["outputs"][0]["relative_path"] = "../outside.tsv"
    elif case == "duplicate_output":
        marker["outputs"].append(dict(marker["outputs"][0]))
    elif case == "wrong_checksum":
        marker["outputs"][0]["sha256"] = "0" * 64
    else:
        marker["outputs"].pop()
    marker_path.write_text(json.dumps(marker), encoding="utf-8")
    with pytest.raises(InputValidationError, match=message):
        verify_prepared_e3_bundle(
            prepared_dir=prepared,
            run_root=root,
            minimum_mean_plddt=50,
            verify_source_authorities=False,
        )


@pytest.mark.parametrize(
    ("case", "message"),
    (
        ("duplicate", "Duplicate prepared source"),
        ("outside", "outside run root"),
        ("size", "size differs"),
        ("checksum", "checksum differs"),
        ("incomplete", "unexpectedly incomplete"),
    ),
)
def test_prepared_source_inventory_fails_closed(
    tmp_path: Path,
    case: str,
    message: str,
) -> None:
    """Predecessor authorities must remain unique, contained and checksum exact."""

    root = _completed_workflow(tmp_path)
    prepared = prepare_e3_workflow_inputs(run_root=root, output_dir=tmp_path / "prepared")
    inventory = prepared / "source_inventory.tsv"
    fields = ("authority", "path", "size_bytes", "sha256")
    records = list(iter_tsv(path=inventory, required_fields=fields))
    if case == "duplicate":
        records.append(dict(records[0]))
    elif case == "outside":
        outside = tmp_path / "outside.tsv"
        outside.write_text("outside\n", encoding="utf-8")
        records[0].update(
            {
                "path": str(outside),
                "size_bytes": str(outside.stat().st_size),
                "sha256": sha256_file(path=outside),
            }
        )
    elif case == "size":
        records[0]["size_bytes"] = str(int(records[0]["size_bytes"]) + 1)
    elif case == "checksum":
        records[0]["sha256"] = "0" * 64
    else:
        records = records[:5]
    write_tsv_atomic(path=inventory, fieldnames=fields, records=records)
    with pytest.raises(InputValidationError, match=message):
        _verify_source_inventory(inventory_path=inventory, run_root=root.resolve())


def test_review_staging_rejects_unsafe_destinations_and_changed_preparation(
    tmp_path: Path,
) -> None:
    """Review staging should protect both immutable inputs and curator destinations."""

    root = _completed_workflow(tmp_path)
    prepared = tmp_path / "work" / "prepared_inputs"
    with pytest.raises(PublicationError, match="outside prepared inputs"):
        ensure_e3_preparation_marker(
            run_root=root,
            prepared_dir=prepared,
            minimum_mean_plddt=50,
            marker_path=prepared / "marker.json",
        )
    preparation_marker = tmp_path / "state" / "prepared.json"
    ensure_e3_preparation_marker(
        run_root=root,
        prepared_dir=prepared,
        minimum_mean_plddt=50,
        marker_path=preparation_marker,
    )
    with pytest.raises(PublicationError, match="outside the immutable"):
        stage_e3_label_review(
            preparation_marker=preparation_marker,
            reviewed_labels=prepared / "reviewed.tsv",
            marker_path=tmp_path / "state" / "review.json",
        )
    same = tmp_path / "work" / "same.json"
    with pytest.raises(PublicationError, match="must not replace"):
        stage_e3_label_review(
            preparation_marker=preparation_marker,
            reviewed_labels=same,
            marker_path=same,
        )
    destination_directory = tmp_path / "work" / "directory"
    destination_directory.mkdir()
    with pytest.raises(PublicationError, match="not a regular file"):
        stage_e3_label_review(
            preparation_marker=preparation_marker,
            reviewed_labels=destination_directory,
            marker_path=tmp_path / "state" / "review.json",
        )
    (prepared / "PREPARED.json").write_bytes((prepared / "PREPARED.json").read_bytes() + b"\n")
    with pytest.raises(InputValidationError, match="changed after workflow"):
        stage_e3_label_review(
            preparation_marker=preparation_marker,
            reviewed_labels=tmp_path / "work" / "reviewed.tsv",
            marker_path=tmp_path / "state" / "review.json",
        )


@pytest.mark.parametrize(
    ("case", "message"),
    (
        ("prepared_dir", "different prepared bundle"),
        ("preparation_path", "different preparation marker"),
        ("preparation_digest", "changed after review staging"),
        ("reviewed_path", "different reviewed-label file"),
        ("template_digest", "template changed after staging"),
        ("missing_review", "missing or empty"),
    ),
)
def test_review_context_rejects_broken_provenance_links(
    tmp_path: Path,
    case: str,
    message: str,
) -> None:
    """Review markers must retain the complete preparation and file-identity chain."""

    root = _completed_workflow(tmp_path)
    prepared = tmp_path / "work" / "prepared_inputs"
    preparation_marker = tmp_path / "state" / "prepared.json"
    review_marker = tmp_path / "state" / "review.json"
    reviewed = tmp_path / "work" / "reviewed.tsv"
    ensure_e3_preparation_marker(
        run_root=root,
        prepared_dir=prepared,
        minimum_mean_plddt=50,
        marker_path=preparation_marker,
    )
    stage_e3_label_review(
        preparation_marker=preparation_marker,
        reviewed_labels=reviewed,
        marker_path=review_marker,
    )
    review = json.loads(review_marker.read_text(encoding="utf-8"))
    if case == "prepared_dir":
        review["prepared_dir"] = str(tmp_path / "other")
    elif case == "preparation_path":
        review["preparation_marker"] = str(tmp_path / "other.json")
    elif case == "preparation_digest":
        preparation_marker.write_bytes(preparation_marker.read_bytes() + b"\n")
    elif case == "reviewed_path":
        review["reviewed_labels"] = str(tmp_path / "other.tsv")
    elif case == "template_digest":
        review["template_sha256"] = "0" * 64
    else:
        reviewed.unlink()
    review_marker.write_text(json.dumps(review), encoding="utf-8")
    with pytest.raises(InputValidationError, match=message):
        _review_context(
            preparation_marker=preparation_marker,
            review_marker=review_marker,
            reviewed_labels=reviewed,
        )


def test_atomic_review_copy_rejects_races_and_io_failures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The no-clobber copy should translate race and I/O failures safely."""

    source = tmp_path / "source.tsv"
    source.write_text("header\n", encoding="utf-8")
    destination = tmp_path / "reviewed.tsv"
    monkeypatch.setattr(
        orchestration_module.os,
        "link",
        lambda *_args: (_ for _ in ()).throw(FileExistsError("race")),
    )
    with pytest.raises(PublicationError, match="appeared concurrently"):
        _copy_file_atomic(source=source, destination=destination)
    assert not destination.exists()

    monkeypatch.undo()
    monkeypatch.setattr(
        orchestration_module.shutil,
        "copyfile",
        lambda *_args: (_ for _ in ()).throw(OSError("copy failed")),
    )
    with pytest.raises(PublicationError, match="copy failed"):
        _copy_file_atomic(source=source, destination=destination)


def test_review_approval_rejects_incomplete_coverage_and_changed_summary(
    tmp_path: Path,
) -> None:
    """Coverage and recorded summary counts are mandatory approval invariants."""

    paths = _approved_workflow(tmp_path)
    approval = json.loads(paths["approval"].read_text(encoding="utf-8"))
    approval["reviewed_assignment_count"] = 999
    paths["approval"].write_text(json.dumps(approval), encoding="utf-8")
    with pytest.raises(InputValidationError, match="summary"):
        verify_e3_label_review(
            preparation_marker=paths["preparation"],
            review_marker=paths["review"],
            reviewed_labels=paths["reviewed"],
            approval_marker=paths["approval"],
            marker_path=tmp_path / "summary_verification.json",
        )

    second = _approved_workflow(tmp_path / "second")
    records = list(iter_tsv(path=second["reviewed"], required_fields=LABEL_FIELDS))
    write_tsv_atomic(
        path=second["reviewed"],
        fieldnames=LABEL_FIELDS,
        records=records[:1],
    )
    second["approval"].unlink()
    with pytest.raises(InputValidationError, match="explicit row for every"):
        approve_e3_label_review(
            preparation_marker=second["preparation"],
            review_marker=second["review"],
            reviewed_labels=second["reviewed"],
            approval_marker=second["approval"],
            curator="Test Curator",
        )


@pytest.mark.parametrize(
    ("case", "message"),
    (
        ("prepared", "Prepared marker changed"),
        ("labels", "Reviewed labels changed"),
        ("approval", "approval changed"),
        ("malformed_prepared", "must contain an object"),
        ("too_few_structures", "At least two Foldseek-eligible"),
    ),
)
def test_campaign_boundary_rejects_changed_approved_authorities(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    case: str,
    message: str,
) -> None:
    """Campaign creation should recheck every approved upstream identity."""

    paths = _approved_workflow(tmp_path)
    if case == "prepared":
        marker = paths["prepared"] / "PREPARED.json"
        marker.write_bytes(marker.read_bytes() + b"\n")
    elif case == "labels":
        paths["reviewed"].write_bytes(paths["reviewed"].read_bytes() + b"\n")
    elif case == "approval":
        paths["approval"].write_bytes(paths["approval"].read_bytes() + b"\n")
    elif case == "malformed_prepared":
        monkeypatch.setattr(orchestration_module, "read_json", lambda **_kwargs: [])
    else:
        prepared = json.loads((paths["prepared"] / "PREPARED.json").read_text(encoding="utf-8"))
        prepared["foldseek_eligible_structure_count"] = 1
        monkeypatch.setattr(orchestration_module, "read_json", lambda **_kwargs: prepared)
    with pytest.raises(InputValidationError, match=message):
        ensure_e3_campaign_marker(
            preparation_marker=paths["preparation"],
            review_verification_marker=paths["verification"],
            campaign_config=paths["config"],
            campaign_id="fixture_campaign",
            profile="e3",
            marker_path=paths["campaign_marker"],
        )


@pytest.mark.parametrize(
    ("case", "message"),
    (
        ("campaign_id", "changed orchestration field"),
        ("foldseek", "keep Foldseek enabled"),
        ("alphafold", "reuse prepared AlphaFold"),
        ("orthofinder_resource", "raw OrthoFinder results"),
    ),
)
def test_existing_campaign_cannot_change_bridge_authorities(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    case: str,
    message: str,
) -> None:
    """Adopted campaign YAML may tune analyses but not replace input authorities."""

    paths = _approved_workflow(tmp_path)
    monkeypatch.setattr(
        orchestration_module,
        "validate_campaign",
        lambda **_kwargs: {"status": "VALID"},
    )
    ensure_e3_campaign_marker(
        preparation_marker=paths["preparation"],
        review_verification_marker=paths["verification"],
        campaign_config=paths["config"],
        campaign_id="fixture_campaign",
        profile="e3",
        marker_path=paths["campaign_marker"],
    )
    document = yaml.safe_load(paths["config"].read_text(encoding="utf-8"))
    if case == "foldseek":
        document["foldseek"]["enabled"] = False
    elif case == "alphafold":
        document["alphafold"]["enabled"] = True
    elif case == "orthofinder_resource":
        document["inputs"]["orthofinder"]["resource_dir"] = str(paths["root"])
        document["inputs"]["orthofinder"]["results_dir"] = None
    paths["config"].write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    prepared_document = json.loads(
        (paths["prepared"] / "PREPARED.json").read_text(encoding="utf-8")
    )
    with pytest.raises(InputValidationError, match=message):
        _validate_existing_e3_campaign(
            config_path=paths["config"],
            campaign_id=("different_campaign" if case == "campaign_id" else "fixture_campaign"),
            profile="e3",
            prepared_dir=paths["prepared"],
            reviewed_labels=paths["reviewed"],
            structural_resource=Path(prepared_document["structural_alignment_resource"]),
            orthofinder_results=Path(prepared_document["orthofinder_results"]),
            foldseek_maximum_hits=2,
        )


def test_campaign_validation_rejects_in_flight_configuration_change(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The campaign checksum should remain stable throughout expensive validation."""

    paths = _approved_workflow(tmp_path)

    def mutate_campaign(*, config_path: Path) -> dict[str, str]:
        """Change the validated file to simulate a concurrent editor."""

        config_path.write_bytes(config_path.read_bytes() + b"\n# concurrent edit\n")
        return {"status": "VALID"}

    monkeypatch.setattr(orchestration_module, "validate_campaign", mutate_campaign)
    with pytest.raises(InputValidationError, match="changed during"):
        ensure_e3_campaign_marker(
            preparation_marker=paths["preparation"],
            review_verification_marker=paths["verification"],
            campaign_config=paths["config"],
            campaign_id="fixture_campaign",
            profile="e3",
            marker_path=paths["campaign_marker"],
        )


def test_prepared_verification_rejects_missing_directory(tmp_path: Path) -> None:
    """A missing prepared directory should fail before marker access."""

    with pytest.raises(InputValidationError, match="not a directory"):
        verify_prepared_e3_bundle(
            prepared_dir=tmp_path / "missing",
            run_root=tmp_path / "run",
            minimum_mean_plddt=50,
        )


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
    projected = tuple(
        _iter_parquet_records(
            path=table,
            required=("id",),
            optional=("absent_optional",),
            batch_size=1,
        )
    )
    assert projected == ({"id": "x"},)
    with pytest.raises(InputValidationError, match="non-empty and unique"):
        _read_parquet_records(path=table, required=("id", "id"))
    with pytest.raises(InputValidationError, match="Optional Parquet fields"):
        _read_parquet_records(path=table, required=("id",), optional=("id",))
    with pytest.raises(InputValidationError, match="batch_size"):
        _read_parquet_records(path=table, required=("id",), batch_size=0)
    with pytest.raises(InputValidationError, match="batch_size"):
        _read_parquet_records(path=table, required=("id",), batch_size=True)
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
    composite = "Arabidopsis_thaliana@@sp|B3H578|PHD1_ARATH"
    with_context = _prepare_sequences(rows=({**base, "cluster_id": composite},))
    assert json.loads(with_context.audit_records[0]["cluster_ids"]) == [composite]
    assert _serialise_review_values(values=("second", "first", "second")) == ('["first","second"]')
    with pytest.raises(InputValidationError, match="TSV-breaking"):
        _prepare_sequences(rows=({**base, "cluster_id": "bad\tcluster"},))
    with pytest.raises(InputValidationError, match="must match"):
        _prepare_sequences(rows=({**base, "group_id": "bad|group"},))
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


def test_prepared_fasta_writer_streams_sorted_validated_records(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Prepared FASTA publication should be atomic, sorted and independently readable."""

    destination = tmp_path / "nested" / "proteins.faa"
    count = _write_prepared_fasta(
        path=destination,
        sequences={"Q00002": "MAAA", "P00001": "MACD"},
    )
    assert count == 2
    records = read_protein_fasta(path=destination)
    assert [(record.protein_id, record.sequence) for record in records] == [
        ("P00001", "MACD"),
        ("Q00002", "MAAA"),
    ]
    with pytest.raises(PublicationError, match="at least one sequence"):
        _write_prepared_fasta(path=tmp_path / "empty.faa", sequences={})
    invalid = tmp_path / "invalid.faa"
    with pytest.raises(InputValidationError, match="requires a sequence"):
        _write_prepared_fasta(path=invalid, sequences={"P00001": ""})
    assert not invalid.exists()

    failed = tmp_path / "failed.faa"
    monkeypatch.setattr(
        bridge_module.os,
        "replace",
        lambda *_args: (_ for _ in ()).throw(OSError("replace failure")),
    )
    with pytest.raises(PublicationError, match="replace failure"):
        _write_prepared_fasta(path=failed, sequences={"P00001": "MACD"})
    assert not failed.exists()


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
    monkeypatch.setattr(bridge_module, "iter_protein_fasta", lambda **_kwargs: iter(()))
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
