"""Focused tests for explicitly synthetic automated smoke-test labels."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from protein_signatures.automated_test_labels import (
    AUTOMATED_TEST_EVIDENCE_STATUS,
    _direct_test_target_labels,
    _select_test_cohorts,
    create_automated_test_labels,
)
from protein_signatures.cli import build_parser, main
from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import iter_tsv
from protein_signatures.models import PartitionAssignment
from protein_signatures.profiles import default_profile_comparisons, load_profile
from protein_signatures.tables import LABEL_FIELDS


def _write_unique_fasta(*, path: Path, count: int) -> tuple[str, ...]:
    """Write deterministic proteins with distinct amino-acid sequences."""

    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    protein_ids: list[str] = []
    lines: list[str] = []
    for index in range(count):
        protein_id = f"protein_{index:04d}"
        value = index
        encoded: list[str] = []
        for _ in range(5):
            encoded.append(alphabet[value % len(alphabet)])
            value //= len(alphabet)
        sequence = f"M{''.join(encoded)}ACDEFGHIKLMNPQRSTVWY"
        protein_ids.append(protein_id)
        lines.extend((f">{protein_id}", sequence))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return tuple(protein_ids)


def test_all_generic_profile_defaults_create_disjoint_reproducible_cohorts(
    tmp_path: Path,
) -> None:
    """ALL should exercise every custom-profile comparison without unit leakage."""

    root = Path(__file__).parents[1]
    profile_path = root / "configs/profile.example.yaml"
    fasta = tmp_path / "proteins.faa"
    protein_ids = _write_unique_fasta(path=fasta, count=500)
    labels = tmp_path / "AUTOMATED_TEST_ONLY.labels.tsv"
    marker = tmp_path / "AUTOMATED_TEST_ONLY.labels.json"

    published = create_automated_test_labels(
        sequences_fasta=fasta,
        output_labels=labels,
        marker_path=marker,
        profile=profile_path,
        target_label="ALL",
        samples_per_class=20,
    )

    assert published == marker.resolve()
    rows = list(iter_tsv(path=labels, required_fields=LABEL_FIELDS))
    assert len(rows) == len(protein_ids)
    positives = [row for row in rows if row["curation_status"] == "REVIEWED_POSITIVE"]
    counts = Counter(row["label_id"] for row in positives)
    assert counts == {
        "protein:kinase:serine_threonine": 20,
        "protein:kinase:tyrosine": 20,
        "control:matched_non_kinase": 20,
        "control:matched_non_tyrosine_kinase": 20,
    }
    assert {row["evidence_status"] for row in positives} == {AUTOMATED_TEST_EVIDENCE_STATUS}

    document = json.loads(marker.read_text(encoding="utf-8"))
    assert document["status"] == "AUTOMATED_TEST_ONLY"
    assert document["scientific_interpretation_allowed"] is False
    assert document["comparison_count"] == 2
    assert document["cohort_count"] == 4
    assert all(
        comparison["target_protein_count"] == 20 and comparison["background_protein_count"] == 20
        for comparison in document["comparisons"]
    )
    unit_owners: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for label_id, selected in document["selected_cohorts"].items():
        assert Counter(row["partition"] for row in selected) == {
            "DISCOVERY": 16,
            "VALIDATION": 4,
        }
        assert (
            len(
                {
                    (row["partition_unit"], row["partition_key"])
                    for row in selected
                    if row["partition"] == "DISCOVERY"
                }
            )
            >= 3
        )
        for row in selected:
            unit_owners[
                (
                    row["partition"],
                    row["partition_unit"],
                    row["partition_key"],
                )
            ].add(label_id)
    assert all(len(owners) == 1 for owners in unit_owners.values())

    second_labels = tmp_path / "second.AUTOMATED_TEST_ONLY.labels.tsv"
    second_marker = tmp_path / "second.AUTOMATED_TEST_ONLY.labels.json"
    create_automated_test_labels(
        sequences_fasta=fasta,
        output_labels=second_labels,
        marker_path=second_marker,
        profile=profile_path,
        samples_per_class=20,
    )
    assert second_labels.read_bytes() == labels.read_bytes()


def test_structure_filter_and_alias_work_for_a_non_e3_profile(tmp_path: Path) -> None:
    """A generic target alias should draw both cohorts only from eligible models."""

    root = Path(__file__).parents[1]
    fasta = tmp_path / "proteins.faa"
    protein_ids = _write_unique_fasta(path=fasta, count=400)
    eligible = frozenset(protein_ids[:200])
    structures = tmp_path / "structures.tsv"
    lines = ["protein_id\tavailability_status\tcoordinate_path\tanalysis_eligibility_status"]
    for protein_id in protein_ids:
        is_eligible = protein_id in eligible
        lines.append(
            "\t".join(
                (
                    protein_id,
                    "AVAILABLE" if is_eligible else "MODEL_NOT_AVAILABLE",
                    f"/coordinates/{protein_id}.cif" if is_eligible else "",
                    "ELIGIBLE" if is_eligible else "INELIGIBLE_COORDINATE_UNAVAILABLE",
                )
            )
        )
    structures.write_text("\n".join(lines) + "\n", encoding="utf-8")
    labels = tmp_path / "AUTOMATED_TEST_ONLY.structure_labels.tsv"
    marker = tmp_path / "AUTOMATED_TEST_ONLY.structure_labels.json"

    create_automated_test_labels(
        sequences_fasta=fasta,
        output_labels=labels,
        marker_path=marker,
        profile=root / "configs/profile.example.yaml",
        target_label="STK",
        structures=structures,
        samples_per_class=20,
    )

    document = json.loads(marker.read_text(encoding="utf-8"))
    assert document["comparison_count"] == 1
    assert document["structure_eligible_only"] is True
    selected_ids = {
        row["protein_id"] for rows in document["selected_cohorts"].values() for row in rows
    }
    assert selected_ids <= eligible
    assert len(selected_ids) == 40


def test_partition_allocator_uses_whole_groups_and_fails_if_groups_are_short() -> None:
    """Multi-protein HOG-like units must remain exclusive to one synthetic label."""

    assignments: list[PartitionAssignment] = []
    eligible: set[str] = set()
    for partition, group_count, group_size in (
        ("DISCOVERY", 6, 6),
        ("VALIDATION", 2, 4),
    ):
        for group_index in range(group_count):
            group = f"{partition}_HOG_{group_index}"
            for protein_index in range(group_size):
                protein_id = f"{group}_protein_{protein_index}"
                eligible.add(protein_id)
                assignments.append(
                    PartitionAssignment(
                        protein_id=protein_id,
                        partition=partition,
                        partition_unit="ORTHOFINDER_GROUP",
                        partition_key=group,
                    )
                )
    cohorts = _select_test_cohorts(
        partitions=tuple(assignments),
        eligible_ids=frozenset(eligible),
        cohort_roles={"target": "TARGET", "control": "BACKGROUND"},
        samples_per_class=20,
        validation_fraction=0.2,
        random_seed=7,
    )
    assert {label_id: len(rows) for label_id, rows in cohorts.items()} == {
        "target": 20,
        "control": 20,
    }
    owners: dict[tuple[str, str], set[str]] = defaultdict(set)
    for label_id, rows in cohorts.items():
        for row in rows:
            owners[(row["partition"], row["partition_key"])].add(label_id)
    assert all(len(labels) == 1 for labels in owners.values())

    with pytest.raises(InputValidationError, match="Insufficient disjoint"):
        _select_test_cohorts(
            partitions=tuple(assignments[:-6]),
            eligible_ids=frozenset(eligible),
            cohort_roles={"target": "TARGET", "control": "BACKGROUND"},
            samples_per_class=20,
            validation_fraction=0.2,
            random_seed=7,
        )


def test_e3_all_contract_and_cli_defaults_cover_every_subclass() -> None:
    """The smoke default should cover all 73 E3 comparisons, including F-box."""

    profile = load_profile(source="e3")
    comparisons = default_profile_comparisons(profile=profile)
    direct_targets = _direct_test_target_labels(
        profile=profile,
        comparisons=comparisons,
    )
    backgrounds = {comparison.background_label_ids[0] for comparison in comparisons}
    assert len(comparisons) == 73
    assert len(direct_targets) == 65
    assert len(backgrounds) == 14
    assert "e3:ubiquitin:crl:crl1_scf:f_box" in {
        comparison.target_label_ids[0] for comparison in comparisons
    }

    arguments = build_parser().parse_args(
        (
            "create-automated-test-labels",
            "--sequences-fasta",
            "proteins.faa",
            "--output-labels",
            "AUTOMATED_TEST_ONLY.labels.tsv",
            "--marker",
            "AUTOMATED_TEST_ONLY.labels.json",
        )
    )
    assert arguments.target_label == "ALL"
    assert arguments.samples_per_class == 20


def test_e3_all_generation_populates_every_inherited_comparison(tmp_path: Path) -> None:
    """Terminal synthetic cohorts should populate all parent and F-box comparisons."""

    fasta = tmp_path / "proteins.faa"
    _write_unique_fasta(path=fasta, count=3_000)
    labels = tmp_path / "AUTOMATED_TEST_ONLY.e3_labels.tsv"
    marker = tmp_path / "AUTOMATED_TEST_ONLY.e3_labels.json"
    create_automated_test_labels(
        sequences_fasta=fasta,
        output_labels=labels,
        marker_path=marker,
        profile="e3",
        target_label="ALL",
        samples_per_class=20,
    )

    document = json.loads(marker.read_text(encoding="utf-8"))
    assert document["comparison_count"] == 73
    assert document["cohort_count"] == 79
    assert len(document["direct_target_label_ids"]) == 65
    assert len(document["background_label_ids"]) == 14
    assert all(
        comparison["target_protein_count"] >= 20 and comparison["background_protein_count"] >= 20
        for comparison in document["comparisons"]
    )
    f_box = next(
        comparison
        for comparison in document["comparisons"]
        if comparison["target_label_id"] == "e3:ubiquitin:crl:crl1_scf:f_box"
    )
    assert f_box["target_protein_count"] > 20


def test_cli_routes_generic_automated_label_arguments(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The public command should forward generic authorities and report test status."""

    observed: dict[str, object] = {}
    marker = tmp_path / "AUTOMATED_TEST_ONLY.labels.json"

    def fake_create(**kwargs: object) -> Path:
        """Record the parsed CLI contract without reading scientific inputs."""

        observed.update(kwargs)
        return marker.resolve()

    monkeypatch.setattr(
        "protein_signatures.cli.create_automated_test_labels",
        fake_create,
    )
    exit_code = main(
        (
            "create-automated-test-labels",
            "--sequences-fasta",
            str(tmp_path / "proteins.faa"),
            "--output-labels",
            str(tmp_path / "AUTOMATED_TEST_ONLY.labels.tsv"),
            "--marker",
            str(marker),
            "--profile",
            str(tmp_path / "profile.yaml"),
            "--orthofinder-resource",
            str(tmp_path / "orthofinder_resource"),
            "--redundancy-clusters",
            str(tmp_path / "redundancy.tsv"),
        )
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload == {
        "marker": str(marker.resolve()),
        "status": "AUTOMATED_TEST_ONLY",
    }
    assert observed["target_label"] == "ALL"
    assert observed["samples_per_class"] == 20
    assert observed["orthofinder_resource"] == tmp_path / "orthofinder_resource"
    assert observed["redundancy_clusters"] == tmp_path / "redundancy.tsv"


def test_synthetic_outputs_require_conspicuous_names_and_sufficient_data(
    tmp_path: Path,
) -> None:
    """Synthetic generation should fail closed on unsafe naming and tiny inputs."""

    root = Path(__file__).parents[1]
    fasta = tmp_path / "proteins.faa"
    _write_unique_fasta(path=fasta, count=30)
    with pytest.raises(InputValidationError, match="filename must contain"):
        create_automated_test_labels(
            sequences_fasta=fasta,
            output_labels=tmp_path / "looks_scientific.tsv",
            marker_path=tmp_path / "AUTOMATED_TEST_ONLY.marker.json",
            profile=root / "configs/profile.example.yaml",
        )
    with pytest.raises(InputValidationError, match="Insufficient disjoint"):
        create_automated_test_labels(
            sequences_fasta=fasta,
            output_labels=tmp_path / "AUTOMATED_TEST_ONLY.labels.tsv",
            marker_path=tmp_path / "AUTOMATED_TEST_ONLY.marker.json",
            profile=root / "configs/profile.example.yaml",
        )
