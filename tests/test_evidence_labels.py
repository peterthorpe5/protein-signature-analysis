"""Focused tests for evidence-led production label proposals."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest
import yaml

import protein_signatures.cli as cli_module
from protein_signatures.cli import build_parser
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.evidence_labels import (
    CLASS_SUMMARY_FIELDS,
    CONTROL_MATCH_FIELDS,
    EVIDENCE_AUDIT_FIELDS,
    EVIDENCE_BUNDLE_STATUS,
    EVIDENCE_CONTROL_STATUS,
    EVIDENCE_POSITIVE_STATUS,
    EVIDENCE_TARGET_STATUS,
    LABEL_DEFINITION_FIELDS,
    UNRESOLVED_FIELDS,
    create_evidence_label_bundle,
    load_evidence_rules,
    read_evidence_audit_table,
    verify_evidence_label_bundle,
)
from protein_signatures.io_utils import iter_tsv, write_tsv_atomic
from protein_signatures.pipeline import validate_campaign
from protein_signatures.profiles import default_profile_comparisons, load_profile
from protein_signatures.starter import initialise_campaign
from protein_signatures.tables import DOMAIN_FIELDS, LABEL_FIELDS


def _write_profile(*, path: Path) -> None:
    """Write a small non-E3 profile with two target classes."""

    document = {
        "profile_id": "enzyme_test",
        "profile_version": "1.0.0",
        "display_name": "Evidence labelling test enzymes",
        "require_structural_evidence": False,
        "default_comparison": {
            "target_root_label_id": "protein",
            "background_label_id": "control:matched_reference",
            "excluded_label_ids": [],
            "excluded_subtree_label_ids": [],
            "description": "Prespecified test comparison.",
        },
        "labels": [
            {
                "label_id": "protein",
                "display_name": "Protein",
                "parent_label_id": "",
                "level": "root",
                "reviewed_positive_allowed": False,
                "description": "Profile root.",
            },
            {
                "label_id": "protein:kinase",
                "display_name": "Kinase",
                "parent_label_id": "protein",
                "level": "family",
                "component_role": "CATALYTIC_PROTEIN",
                "default_analysis": True,
                "assignment_exclusivity_group": "primary_enzyme_type",
                "description": "Protein kinase target.",
            },
            {
                "label_id": "protein:phosphatase",
                "display_name": "Phosphatase",
                "parent_label_id": "protein",
                "level": "family",
                "component_role": "CATALYTIC_PROTEIN",
                "default_analysis": True,
                "assignment_exclusivity_group": "primary_enzyme_type",
                "description": "Protein phosphatase target.",
            },
            {
                "label_id": "control:matched_reference",
                "display_name": "Matched reference",
                "parent_label_id": "protein",
                "level": "control",
                "component_role": "CONTROL",
                "description": "Outcome-blind matched background.",
            },
        ],
    }
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")


def _write_rules(*, path: Path) -> None:
    """Write strict generic rules for the test profile."""

    settings = {
        "annotation_score": 2.0,
        "specific_annotation_score": 4.0,
        "pfam_score": 2.0,
        "supporting_pfam_score": 0.5,
        "orthology_score": 2.0,
        "trusted_label_score": 5.0,
        "minimum_acceptance_score": 3.5,
        "minimum_evidence_groups": 2,
        "minimum_score_margin": 1.0,
        "require_target_structure_eligible": False,
        "orthology_propagation_enabled": True,
        "minimum_orthology_anchor_proteins": 2,
        "control_units_per_target_unit": 1,
        "minimum_control_units_per_background": 2,
        "require_species_match": True,
        "require_structure_match": True,
        "maximum_log2_length_difference": 1.0,
        "maximum_domain_count_difference": 2,
        "exclude_input_candidates_from_controls": True,
    }
    rules = []
    for name, accession in (("kinase", "PF00069"), ("phosphatase", "PF00102")):
        rules.append(
            {
                "rule_id": f"{name}_rule",
                "label_id": f"protein:{name}",
                "priority": 100,
                "annotation_patterns": [rf"\b{name}\b"],
                "specific_annotation_patterns": [rf"\b{name} enzyme\b"],
                "exclusion_patterns": [],
                "required_domain_patterns": [rf"^{accession}\b"],
                "supporting_domain_patterns": [],
                "allow_specific_annotation_only": False,
                "propagate_by_orthology": True,
            }
        )
    document = {
        "schema_version": 1,
        "ruleset_id": "enzyme_test_rules",
        "ruleset_version": "1.0.0",
        "profile_id": "enzyme_test",
        "settings": settings,
        "label_rules": rules,
    }
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")


def _write_inputs(*, root: Path) -> tuple[Path, Path, Path]:
    """Write sequences, generic metadata and explicit Pfam assessment rows."""

    protein_ids = (
        "kinase_a",
        "kinase_b",
        "phosphatase_a",
        "phosphatase_b",
        "ambiguous",
        "weak_annotation",
        "control_a",
        "control_b",
        "control_c",
        "control_d",
    )
    fasta = root / "proteins.faa"
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    fasta.write_text(
        "".join(
            f">{protein_id}\nM{alphabet[index:]}{alphabet[:index]}ACDEFGH\n"
            for index, protein_id in enumerate(protein_ids)
        ),
        encoding="utf-8",
    )
    metadata = root / "protein_metadata.tsv"
    write_tsv_atomic(
        path=metadata,
        fieldnames=("protein_id", "species", "input_candidate"),
        records=(
            {
                "protein_id": protein_id,
                "species": "Test_species",
                "input_candidate": protein_id.startswith(("kinase", "phosphatase"))
                or protein_id in {"ambiguous", "weak_annotation"},
            }
            for protein_id in protein_ids
        ),
    )
    domains = root / "domains.tsv"
    accessions = {
        "kinase_a": ("PF00069",),
        "kinase_b": ("PF00069",),
        "phosphatase_a": ("PF00102",),
        "phosphatase_b": ("PF00102",),
        "ambiguous": ("PF00069", "PF00102"),
    }
    rows = []
    for protein_id in protein_ids:
        protein_accessions = accessions.get(protein_id, ())
        if not protein_accessions:
            rows.append(
                {
                    "protein_id": protein_id,
                    "domain_authority": "Pfam",
                    "assessment_status": "ASSESSED_NO_HIT",
                    "domain_id": "",
                    "domain_name": "",
                    "start": "",
                    "end": "",
                    "score": "",
                    "e_value": "",
                    "evidence_source": "Test Pfam release",
                    "evidence_reference": "test-pfam-1",
                }
            )
            continue
        for domain_index, accession in enumerate(protein_accessions):
            rows.append(
                {
                    "protein_id": protein_id,
                    "domain_authority": "Pfam",
                    "assessment_status": "ASSESSED_WITH_HIT",
                    "domain_id": accession,
                    "domain_name": (
                        "Protein kinase domain" if accession == "PF00069" else "Phosphatase"
                    ),
                    "start": 2 + domain_index * 5,
                    "end": 6 + domain_index * 5,
                    "score": 40.0,
                    "e_value": "1e-20",
                    "evidence_source": "Test Pfam release",
                    "evidence_reference": "test-pfam-1",
                }
            )
    write_tsv_atomic(path=domains, fieldnames=DOMAIN_FIELDS, records=rows)
    return fasta, metadata, domains


def _write_annotations(*, path: Path) -> None:
    """Write direct-protein annotations, including one conflict and one weak call."""

    annotations = {
        "kinase_a": "putative kinase",
        "kinase_b": "putative kinase",
        "phosphatase_a": "putative phosphatase",
        "phosphatase_b": "putative phosphatase",
        "ambiguous": "kinase phosphatase",
        "weak_annotation": "putative kinase",
    }
    write_tsv_atomic(
        path=path,
        fieldnames=(
            "protein_id",
            "label_id",
            "annotation_text",
            "evidence_status",
            "evidence_source",
            "evidence_reference",
            "annotation_scope",
        ),
        records=(
            {
                "protein_id": protein_id,
                "label_id": "",
                "annotation_text": text,
                "evidence_status": "CURATED",
                "evidence_source": "Independent annotation authority",
                "evidence_reference": f"annotation:{protein_id}",
                "annotation_scope": "EXACT_PROTEIN",
            }
            for protein_id, text in annotations.items()
        ),
    )


def test_evidence_bundle_is_generic_auditable_and_circularity_safe(
    tmp_path: Path,
) -> None:
    """Generic evidence should accept corroborated targets and match clean controls."""

    profile = tmp_path / "profile.yaml"
    rules = tmp_path / "rules.yaml"
    annotations = tmp_path / "annotations.tsv"
    _write_profile(path=profile)
    _write_rules(path=rules)
    fasta, metadata, domains = _write_inputs(root=tmp_path)
    _write_annotations(path=annotations)
    bundle = tmp_path / "evidence_bundle"

    marker = create_evidence_label_bundle(
        sequences_fasta=fasta,
        output_dir=bundle,
        profile=profile,
        evidence_rules=rules,
        protein_metadata=metadata,
        domains=domains,
        external_annotations=annotations,
    )

    assert marker == (bundle / "EVIDENCE_LABELS.json").resolve()
    document = verify_evidence_label_bundle(bundle_dir=bundle)
    assert document["status"] == EVIDENCE_BUNDLE_STATUS
    assert document["human_review_completed"] is False
    assert document["evidence_supported_target_count"] == 4
    assert document["matched_control_protein_count"] == 4

    labels = tuple(iter_tsv(path=bundle / "label_assignments.tsv", required_fields=LABEL_FIELDS))
    targets = [row for row in labels if row["evidence_status"] == EVIDENCE_TARGET_STATUS]
    controls = [row for row in labels if row["evidence_status"] == EVIDENCE_CONTROL_STATUS]
    assert len(targets) == 4
    assert len(controls) == 4
    assert {row["curation_status"] for row in targets + controls} == {EVIDENCE_POSITIVE_STATUS}
    unresolved = {
        row["protein_id"]: row
        for row in iter_tsv(
            path=bundle / "unresolved_assignments.tsv",
            required_fields=UNRESOLVED_FIELDS,
        )
    }
    assert unresolved["ambiguous"]["curation_status"] == "AMBIGUOUS"
    assert unresolved["weak_annotation"]["curation_status"] == "PROPOSED"

    filtered_domains = tuple(
        iter_tsv(
            path=bundle / "domains.for_signature_analysis.tsv",
            required_fields=DOMAIN_FIELDS,
        )
    )
    assert not {row["domain_id"] for row in filtered_domains} & {"PF00069", "PF00102"}
    definitions = read_evidence_audit_table(
        path=bundle / "label_definition_features.tsv",
        fields=LABEL_DEFINITION_FIELDS,
    )
    assert {row["feature_id"] for row in definitions} == {"PF00069", "PF00102"}
    assert {row["exclusion_scope"] for row in definitions} == {"GLOBAL_CONFIRMATORY_DOMAIN_AND_ML"}

    matches = read_evidence_audit_table(
        path=bundle / "control_matching_audit.tsv",
        fields=CONTROL_MATCH_FIELDS,
    )
    assert len(matches) == 4
    assert {row["status"] for row in matches} == {"MATCHED"}
    assert all(row["target_unit_id"] != row["control_unit_id"] for row in matches)
    assert read_evidence_audit_table(
        path=bundle / "label_evidence_audit.tsv",
        fields=EVIDENCE_AUDIT_FIELDS,
    )
    assert read_evidence_audit_table(
        path=bundle / "class_labelling_summary.tsv",
        fields=CLASS_SUMMARY_FIELDS,
    )

    config = initialise_campaign(
        config_path=tmp_path / "campaign.yaml",
        campaign_id="generic_evidence_test",
        profile=str(profile),
        sequences_fasta=fasta,
        label_assignments=bundle / "label_assignments.tsv",
        label_evidence_marker=marker,
        label_evidence_audit=bundle / "label_evidence_audit.tsv",
        control_matching_audit=bundle / "control_matching_audit.tsv",
        label_definition_features=bundle / "label_definition_features.tsv",
        class_labelling_summary=bundle / "class_labelling_summary.tsv",
        unresolved_assignments=bundle / "unresolved_assignments.tsv",
        domains=bundle / "domains.for_signature_analysis.tsv",
    )
    validation = validate_campaign(config_path=config)
    assert validation["automated_label_evidence"]["status"] == ("PROVISIONAL_EVIDENCE_SUPPORTED")
    metadata.write_text(metadata.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Evidence input size differs"):
        verify_evidence_label_bundle(bundle_dir=bundle)


def test_bundle_verification_detects_a_modified_audit(tmp_path: Path) -> None:
    """Any post-publication audit edit should invalidate the bundle checksum."""

    profile = tmp_path / "profile.yaml"
    rules = tmp_path / "rules.yaml"
    annotations = tmp_path / "annotations.tsv"
    _write_profile(path=profile)
    _write_rules(path=rules)
    fasta, metadata, domains = _write_inputs(root=tmp_path)
    _write_annotations(path=annotations)
    bundle = tmp_path / "evidence_bundle"
    create_evidence_label_bundle(
        sequences_fasta=fasta,
        output_dir=bundle,
        profile=profile,
        evidence_rules=rules,
        protein_metadata=metadata,
        domains=domains,
        external_annotations=annotations,
    )
    audit = bundle / "control_matching_audit.tsv"
    audit.write_text(audit.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(InputValidationError, match="size differs|checksum differs"):
        verify_evidence_label_bundle(bundle_dir=bundle)


def test_evidence_bundle_allows_only_explicit_empty_workflow_directory(
    tmp_path: Path,
) -> None:
    """Workflow compatibility may remove an empty directory but never existing data."""

    profile = tmp_path / "profile.yaml"
    rules = tmp_path / "rules.yaml"
    annotations = tmp_path / "annotations.tsv"
    _write_profile(path=profile)
    _write_rules(path=rules)
    fasta, metadata, domains = _write_inputs(root=tmp_path)
    _write_annotations(path=annotations)

    default_destination = tmp_path / "default_existing"
    default_destination.mkdir()
    with pytest.raises(PublicationError, match="already exists"):
        create_evidence_label_bundle(
            sequences_fasta=fasta,
            output_dir=default_destination,
            profile=profile,
            evidence_rules=rules,
        )

    workflow_destination = tmp_path / "workflow_created"
    workflow_destination.mkdir()
    marker = create_evidence_label_bundle(
        sequences_fasta=fasta,
        output_dir=workflow_destination,
        profile=profile,
        evidence_rules=rules,
        protein_metadata=metadata,
        domains=domains,
        external_annotations=annotations,
        allow_empty_output_dir=True,
    )
    assert marker == workflow_destination / "EVIDENCE_LABELS.json"

    protected_destination = tmp_path / "protected"
    protected_destination.mkdir()
    protected_file = protected_destination / "retain.txt"
    protected_file.write_text("retain\n", encoding="utf-8")
    with pytest.raises(PublicationError, match="not empty"):
        create_evidence_label_bundle(
            sequences_fasta=fasta,
            output_dir=protected_destination,
            profile=profile,
            evidence_rules=rules,
            allow_empty_output_dir=True,
        )
    assert protected_file.read_text(encoding="utf-8") == "retain\n"

    file_destination = tmp_path / "existing_file"
    file_destination.write_text("retain\n", encoding="utf-8")
    with pytest.raises(PublicationError, match="not an ordinary directory"):
        create_evidence_label_bundle(
            sequences_fasta=fasta,
            output_dir=file_destination,
            profile=profile,
            evidence_rules=rules,
            allow_empty_output_dir=True,
        )
    assert file_destination.read_text(encoding="utf-8") == "retain\n"

    symlink_target = tmp_path / "symlink_target"
    symlink_target.mkdir()
    symlink_destination = tmp_path / "symlink_destination"
    symlink_destination.symlink_to(symlink_target, target_is_directory=True)
    with pytest.raises(PublicationError, match="must not be a symbolic link"):
        create_evidence_label_bundle(
            sequences_fasta=fasta,
            output_dir=symlink_destination,
            profile=profile,
            evidence_rules=rules,
            allow_empty_output_dir=True,
        )
    assert symlink_destination.is_symlink()
    assert symlink_target.is_dir()

    with pytest.raises(InputValidationError, match="must be Boolean"):
        create_evidence_label_bundle(
            sequences_fasta=fasta,
            output_dir=tmp_path / "invalid_flag",
            profile=profile,
            evidence_rules=rules,
            allow_empty_output_dir="true",  # type: ignore[arg-type]
        )


def test_built_in_e3_rules_cover_inferable_defaults_and_match_schema() -> None:
    """Every non-catch-all E3 default should have a validated inference rule."""

    root = Path(__file__).parents[1]
    schema = json.loads((root / "configs/evidence_rules.schema.json").read_text())
    source = root / "src/protein_signatures/data/evidence_rules/e3.yaml"
    raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.validate(instance=raw, schema=schema)

    profile = load_profile(source="e3")
    rules = load_evidence_rules(source="e3", profile=profile)
    automatic = {rule.label_id for rule in rules.rules}
    defaults = {
        comparison.target_label_ids[0]
        for comparison in default_profile_comparisons(profile=profile)
    }
    manual_catch_alls = {label for label in defaults if label.endswith(":other_reviewed")}
    assert len(defaults) == 73
    assert len(automatic) == 67
    assert defaults - automatic == manual_catch_alls
    assert len(manual_catch_alls) == 6
    assert "e3:ubiquitin:crl:crl1_scf:f_box" in automatic


def test_cli_uses_named_generic_evidence_arguments() -> None:
    """The public parser should expose generic evidence and explicit verification routes."""

    arguments = build_parser().parse_args(
        (
            "create-evidence-labels",
            "--sequences-fasta",
            "proteins.faa",
            "--output-dir",
            "evidence_bundle",
            "--profile",
            "profile.yaml",
            "--evidence-rules",
            "rules.yaml",
            "--protein-metadata",
            "metadata.tsv",
            "--allow-empty-output-dir",
        )
    )
    assert arguments.command == "create-evidence-labels"
    assert arguments.protein_metadata == Path("metadata.tsv")
    assert arguments.orthofinder_group_type == "HOG"
    assert arguments.allow_empty_output_dir is True
    verifier = build_parser().parse_args(
        ("verify-evidence-labels", "--bundle-dir", "evidence_bundle")
    )
    assert verifier.command == "verify-evidence-labels"


def test_cli_routes_empty_workflow_output_compatibility_flag(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The workflow-only compatibility flag should reach the atomic publisher."""

    captured: dict[str, object] = {}
    marker = tmp_path / "evidence_bundle" / "EVIDENCE_LABELS.json"

    def fake_create_evidence_label_bundle(**kwargs: object) -> Path:
        """Capture named CLI arguments without generating a real bundle."""

        captured.update(kwargs)
        return marker

    monkeypatch.setattr(
        cli_module,
        "create_evidence_label_bundle",
        fake_create_evidence_label_bundle,
    )
    exit_code = cli_module.main(
        [
            "create-evidence-labels",
            "--sequences-fasta",
            str(tmp_path / "proteins.faa"),
            "--output-dir",
            str(marker.parent),
            "--allow-empty-output-dir",
        ]
    )
    assert exit_code == 0
    assert captured["allow_empty_output_dir"] is True
    assert json.loads(capsys.readouterr().out)["marker"] == str(marker)


def test_generic_metadata_requires_complete_strict_boolean_coverage(tmp_path: Path) -> None:
    """Control covariates must not silently accept missing or ambiguous metadata."""

    profile = tmp_path / "profile.yaml"
    rules = tmp_path / "rules.yaml"
    annotations = tmp_path / "annotations.tsv"
    _write_profile(path=profile)
    _write_rules(path=rules)
    fasta, metadata, domains = _write_inputs(root=tmp_path)
    _write_annotations(path=annotations)
    rows = list(
        iter_tsv(
            path=metadata,
            required_fields=("protein_id", "species", "input_candidate"),
        )
    )
    rows[0]["input_candidate"] = "MAYBE"
    write_tsv_atomic(
        path=metadata,
        fieldnames=("protein_id", "species", "input_candidate"),
        records=rows,
    )

    with pytest.raises(InputValidationError, match="must be TRUE or FALSE"):
        create_evidence_label_bundle(
            sequences_fasta=fasta,
            output_dir=tmp_path / "bad_bundle",
            profile=profile,
            evidence_rules=rules,
            protein_metadata=metadata,
            domains=domains,
            external_annotations=annotations,
        )
