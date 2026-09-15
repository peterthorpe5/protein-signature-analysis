"""Command-line interface for validation, analysis and profile inspection."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .automated_test_labels import create_automated_test_labels
from .catalogue import prepare_catalogue
from .e3_workflow_bridge import prepare_e3_workflow_inputs
from .e3_workflow_orchestration import (
    approve_e3_label_review,
    ensure_e3_campaign_marker,
    ensure_e3_preparation_marker,
    stage_e3_label_review,
    verify_e3_label_review,
)
from .errors import ProteinSignatureError
from .evidence_labels import create_evidence_label_bundle, verify_evidence_label_bundle
from .logging_config import configure_logging
from .pipeline import run_campaign, validate_campaign
from .profiles import load_profile
from .publication import verify_completed_result
from .starter import initialise_campaign
from .workflow_markers import publish_validation_marker, publish_verification_marker

LOGGER = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(
        prog="protein-signatures",
        description=(
            "Discover sequence, domain, fold and structural signatures from supplied protein data."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    catalogue_parser = subparsers.add_parser(
        "prepare-catalogue",
        help="Create conservative FASTA and curation starters from a sequence TSV.",
    )
    catalogue_parser.add_argument("--catalogue", required=True, type=Path)
    catalogue_parser.add_argument("--output-dir", required=True, type=Path)
    catalogue_parser.add_argument("--id-column", required=True)
    catalogue_parser.add_argument("--sequence-column", required=True)
    catalogue_parser.add_argument("--name-column", default="")
    catalogue_parser.add_argument("--proposed-category-column", default="")
    catalogue_parser.add_argument("--starter-label-id", default="e3:associated:unknown")
    catalogue_parser.add_argument("--log-level", default="INFO")
    e3_workflow_parser = subparsers.add_parser(
        "prepare-e3-workflow",
        help="Prepare review-gated inputs from a completed E3 end-to-end run.",
    )
    e3_workflow_parser.add_argument("--run-root", required=True, type=Path)
    e3_workflow_parser.add_argument("--output-dir", required=True, type=Path)
    e3_workflow_parser.add_argument("--minimum-mean-plddt", type=float, default=50.0)
    e3_workflow_parser.add_argument("--log-level", default="INFO")
    automated_labels_parser = subparsers.add_parser(
        "create-automated-test-labels",
        help=("Create conspicuously synthetic target/control labels for software testing only."),
    )
    automated_labels_parser.add_argument("--sequences-fasta", required=True, type=Path)
    automated_labels_parser.add_argument("--output-labels", required=True, type=Path)
    automated_labels_parser.add_argument("--marker", required=True, type=Path)
    automated_labels_parser.add_argument("--profile", default="e3")
    automated_labels_parser.add_argument("--target-label", default="ALL")
    automated_labels_parser.add_argument("--template-labels", type=Path)
    automated_labels_parser.add_argument("--structures", type=Path)
    automated_labels_parser.add_argument("--orthofinder-resource", type=Path)
    automated_labels_parser.add_argument("--orthofinder-results", type=Path)
    automated_labels_parser.add_argument(
        "--orthofinder-group-type",
        choices=("HOG", "LEGACY_ORTHOGROUP"),
        default="HOG",
    )
    automated_labels_parser.add_argument("--orthofinder-hierarchy-node", default="N0")
    automated_labels_parser.add_argument("--orthofinder-run-id", default="automated_smoke_test")
    automated_labels_parser.add_argument("--redundancy-clusters", type=Path)
    automated_labels_parser.add_argument("--samples-per-class", type=int, default=20)
    automated_labels_parser.add_argument("--random-seed", type=int, default=1729)
    automated_labels_parser.add_argument("--validation-fraction", type=float, default=0.2)
    automated_labels_parser.add_argument("--log-level", default="INFO")
    evidence_labels_parser = subparsers.add_parser(
        "create-evidence-labels",
        help=("Create provisional evidence-supported labels and outcome-blind matched controls."),
    )
    evidence_labels_parser.add_argument("--sequences-fasta", required=True, type=Path)
    evidence_labels_parser.add_argument("--output-dir", required=True, type=Path)
    evidence_labels_parser.add_argument("--profile", default="e3")
    evidence_labels_parser.add_argument("--evidence-rules", default="e3")
    evidence_labels_parser.add_argument("--protein-metadata", type=Path)
    evidence_labels_parser.add_argument("--review-context", type=Path)
    evidence_labels_parser.add_argument("--domains", type=Path)
    evidence_labels_parser.add_argument("--structures", type=Path)
    evidence_labels_parser.add_argument("--template-labels", type=Path)
    evidence_labels_parser.add_argument("--seed-assignments", type=Path)
    evidence_labels_parser.add_argument("--seed-catalogue", type=Path)
    evidence_labels_parser.add_argument("--external-annotations", type=Path)
    evidence_labels_parser.add_argument("--orthofinder-resource", type=Path)
    evidence_labels_parser.add_argument("--orthofinder-results", type=Path)
    evidence_labels_parser.add_argument(
        "--orthofinder-group-type",
        choices=("HOG", "LEGACY_ORTHOGROUP"),
        default="HOG",
    )
    evidence_labels_parser.add_argument("--orthofinder-hierarchy-node", default="N0")
    evidence_labels_parser.add_argument("--orthofinder-run-id", default="evidence_labelling")
    evidence_labels_parser.add_argument("--redundancy-clusters", type=Path)
    evidence_labels_parser.add_argument("--random-seed", type=int, default=1729)
    evidence_labels_parser.add_argument("--validation-fraction", type=float, default=0.2)
    evidence_labels_parser.add_argument(
        "--allow-empty-output-dir",
        action="store_true",
        help=(
            "Allow an existing empty output directory created by a workflow engine; "
            "non-empty directories remain protected."
        ),
    )
    evidence_labels_parser.add_argument("--log-level", default="INFO")
    verify_evidence_parser = subparsers.add_parser(
        "verify-evidence-labels",
        help="Verify a complete evidence-label bundle and every declared checksum.",
    )
    verify_evidence_parser.add_argument("--bundle-dir", required=True, type=Path)
    verify_evidence_parser.add_argument("--log-level", default="INFO")
    e3_prepare_parser = subparsers.add_parser(
        "workflow-prepare-e3",
        help="Create or verify completed-E3 inputs and publish a workflow marker.",
    )
    e3_prepare_parser.add_argument("--run-root", required=True, type=Path)
    e3_prepare_parser.add_argument("--prepared-dir", required=True, type=Path)
    e3_prepare_parser.add_argument("--minimum-mean-plddt", type=float, default=50.0)
    e3_prepare_parser.add_argument("--marker", required=True, type=Path)
    e3_prepare_parser.add_argument("--log-level", default="INFO")
    e3_stage_review_parser = subparsers.add_parser(
        "workflow-stage-e3-review",
        help="Safely create or adopt the editable E3 label-review authority.",
    )
    e3_stage_review_parser.add_argument("--preparation-marker", required=True, type=Path)
    e3_stage_review_parser.add_argument("--reviewed-labels", required=True, type=Path)
    e3_stage_review_parser.add_argument("--marker", required=True, type=Path)
    e3_stage_review_parser.add_argument("--log-level", default="INFO")
    e3_approve_review_parser = subparsers.add_parser(
        "approve-e3-review",
        help="Validate reviewed labels and bind curator approval to their checksum.",
    )
    e3_approve_review_parser.add_argument("--preparation-marker", required=True, type=Path)
    e3_approve_review_parser.add_argument("--review-marker", required=True, type=Path)
    e3_approve_review_parser.add_argument("--reviewed-labels", required=True, type=Path)
    e3_approve_review_parser.add_argument("--approval-marker", required=True, type=Path)
    e3_approve_review_parser.add_argument("--curator", required=True)
    e3_approve_review_parser.add_argument("--note", default="")
    e3_approve_review_parser.add_argument("--profile", default="e3")
    approval_modes = e3_approve_review_parser.add_mutually_exclusive_group()
    approval_modes.add_argument("--automated-test-mode", action="store_true")
    approval_modes.add_argument("--automated-evidence-mode", action="store_true")
    e3_approve_review_parser.add_argument("--automated-test-marker", type=Path)
    e3_approve_review_parser.add_argument("--evidence-label-marker", type=Path)
    e3_approve_review_parser.add_argument("--log-level", default="INFO")
    e3_verify_review_parser = subparsers.add_parser(
        "workflow-verify-e3-review",
        help="Verify reviewed labels against their checksum-bound approval.",
    )
    e3_verify_review_parser.add_argument("--preparation-marker", required=True, type=Path)
    e3_verify_review_parser.add_argument("--review-marker", required=True, type=Path)
    e3_verify_review_parser.add_argument("--reviewed-labels", required=True, type=Path)
    e3_verify_review_parser.add_argument("--approval-marker", required=True, type=Path)
    e3_verify_review_parser.add_argument("--marker", required=True, type=Path)
    e3_verify_review_parser.add_argument("--profile", default="e3")
    e3_verify_review_parser.add_argument("--log-level", default="INFO")
    e3_initialise_parser = subparsers.add_parser(
        "workflow-initialise-e3",
        help="Create or adopt a reviewed completed-E3 campaign and validate it.",
    )
    e3_initialise_parser.add_argument("--preparation-marker", required=True, type=Path)
    e3_initialise_parser.add_argument("--review-verification-marker", required=True, type=Path)
    e3_initialise_parser.add_argument("--campaign-config", required=True, type=Path)
    e3_initialise_parser.add_argument("--campaign-id", required=True)
    e3_initialise_parser.add_argument("--profile", default="e3")
    e3_initialise_parser.add_argument("--marker", required=True, type=Path)
    e3_initialise_parser.add_argument("--log-level", default="INFO")
    initialise_parser = subparsers.add_parser(
        "initialise",
        help="Create a validated campaign YAML from explicit input authorities.",
    )
    initialise_parser.add_argument("--config", required=True, type=Path)
    initialise_parser.add_argument("--campaign-id", required=True)
    initialise_parser.add_argument("--profile", default="e3")
    initialise_parser.add_argument("--sequences-fasta", required=True, type=Path)
    initialise_parser.add_argument("--label-assignments", required=True, type=Path)
    for option in (
        "label-evidence-marker",
        "label-evidence-audit",
        "control-matching-audit",
        "label-definition-features",
        "class-labelling-summary",
        "unresolved-assignments",
        "features",
        "domains",
        "redundancy-clusters",
        "structures",
        "structure-comparisons",
        "structural-alignment-resource",
        "alphafold-accessions",
        "orthofinder-resource",
        "orthofinder-results",
    ):
        initialise_parser.add_argument(f"--{option}", type=Path)
    initialise_parser.add_argument(
        "--orthofinder-group-type",
        choices=("HOG", "LEGACY_ORTHOGROUP"),
        default="HOG",
    )
    initialise_parser.add_argument("--orthofinder-hierarchy-node", default="N0")
    initialise_parser.add_argument("--orthofinder-run-id", default="")
    initialise_parser.add_argument("--enable-alphafold", action="store_true")
    initialise_parser.add_argument("--enable-foldseek", action="store_true")
    initialise_parser.add_argument("--foldseek-maximum-hits", type=int, default=1000)
    initialise_parser.add_argument("--log-level", default="INFO")
    run_parser = subparsers.add_parser("run-all", help="Run and atomically publish a campaign.")
    run_parser.add_argument("--config", required=True, type=Path)
    run_parser.add_argument("--output-dir", required=True, type=Path)
    run_parser.add_argument("--threads", type=int, default=1)
    run_parser.add_argument("--resume", action="store_true")
    run_parser.add_argument("--log-level", default="INFO")
    validate_parser = subparsers.add_parser(
        "validate", help="Validate local configuration and inputs without downloading models."
    )
    validate_parser.add_argument("--config", required=True, type=Path)
    validate_parser.add_argument("--log-level", default="INFO")
    verify_parser = subparsers.add_parser(
        "verify", help="Verify a completed result manifest and every output checksum."
    )
    verify_parser.add_argument("--resource", required=True, type=Path)
    verify_parser.add_argument("--log-level", default="INFO")
    workflow_validate_parser = subparsers.add_parser(
        "workflow-validate",
        help="Validate a campaign and publish a Snakemake boundary marker.",
    )
    workflow_validate_parser.add_argument("--config", required=True, type=Path)
    workflow_validate_parser.add_argument("--marker", required=True, type=Path)
    workflow_validate_parser.add_argument("--log-level", default="INFO")
    workflow_verify_parser = subparsers.add_parser(
        "workflow-verify",
        help="Verify result and input checksums and publish a Snakemake boundary marker.",
    )
    workflow_verify_parser.add_argument("--resource", required=True, type=Path)
    workflow_verify_parser.add_argument("--marker", required=True, type=Path)
    workflow_verify_parser.add_argument("--log-level", default="INFO")
    profile_parser = subparsers.add_parser(
        "describe-profile", help="Print a built-in or custom profile as JSON."
    )
    profile_parser.add_argument("--profile", default="e3")
    profile_parser.add_argument("--log-level", default="WARNING")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the command-line interface.

    Args:
        argv: Optional argument sequence, excluding program name.

    Returns:
        Process exit code.
    """

    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        configure_logging(level=arguments.log_level)
        if arguments.command == "prepare-catalogue":
            destination = prepare_catalogue(
                catalogue_path=arguments.catalogue,
                output_dir=arguments.output_dir,
                id_column=arguments.id_column,
                sequence_column=arguments.sequence_column,
                name_column=arguments.name_column,
                proposed_category_column=arguments.proposed_category_column,
                starter_label_id=arguments.starter_label_id,
            )
            print(json.dumps({"status": "COMPLETE", "starter_dir": str(destination)}))
        elif arguments.command == "prepare-e3-workflow":
            destination = prepare_e3_workflow_inputs(
                run_root=arguments.run_root,
                output_dir=arguments.output_dir,
                minimum_mean_plddt=arguments.minimum_mean_plddt,
            )
            print(json.dumps({"status": "COMPLETE", "prepared_dir": str(destination)}))
        elif arguments.command == "create-automated-test-labels":
            destination = create_automated_test_labels(
                sequences_fasta=arguments.sequences_fasta,
                output_labels=arguments.output_labels,
                marker_path=arguments.marker,
                profile=arguments.profile,
                target_label=arguments.target_label,
                template_labels=arguments.template_labels,
                structures=arguments.structures,
                orthofinder_resource=arguments.orthofinder_resource,
                orthofinder_results=arguments.orthofinder_results,
                orthofinder_group_type=arguments.orthofinder_group_type,
                orthofinder_hierarchy_node=arguments.orthofinder_hierarchy_node,
                orthofinder_run_id=arguments.orthofinder_run_id,
                redundancy_clusters=arguments.redundancy_clusters,
                samples_per_class=arguments.samples_per_class,
                random_seed=arguments.random_seed,
                validation_fraction=arguments.validation_fraction,
            )
            print(
                json.dumps(
                    {"status": "AUTOMATED_TEST_ONLY", "marker": str(destination)},
                    sort_keys=True,
                )
            )
        elif arguments.command == "create-evidence-labels":
            destination = create_evidence_label_bundle(
                sequences_fasta=arguments.sequences_fasta,
                output_dir=arguments.output_dir,
                profile=arguments.profile,
                evidence_rules=arguments.evidence_rules,
                protein_metadata=arguments.protein_metadata,
                review_context=arguments.review_context,
                domains=arguments.domains,
                structures=arguments.structures,
                template_labels=arguments.template_labels,
                seed_assignments=arguments.seed_assignments,
                seed_catalogue=arguments.seed_catalogue,
                external_annotations=arguments.external_annotations,
                orthofinder_resource=arguments.orthofinder_resource,
                orthofinder_results=arguments.orthofinder_results,
                orthofinder_group_type=arguments.orthofinder_group_type,
                orthofinder_hierarchy_node=arguments.orthofinder_hierarchy_node,
                orthofinder_run_id=arguments.orthofinder_run_id,
                redundancy_clusters=arguments.redundancy_clusters,
                random_seed=arguments.random_seed,
                validation_fraction=arguments.validation_fraction,
                allow_empty_output_dir=arguments.allow_empty_output_dir,
            )
            print(
                json.dumps(
                    {
                        "status": "PROVISIONAL_EVIDENCE_LABELS_COMPLETE",
                        "marker": str(destination),
                    },
                    sort_keys=True,
                )
            )
        elif arguments.command == "verify-evidence-labels":
            document = verify_evidence_label_bundle(bundle_dir=arguments.bundle_dir)
            print(
                json.dumps(
                    {
                        "status": "VALID",
                        "bundle_dir": str(arguments.bundle_dir.expanduser().resolve()),
                        "target_count": document["evidence_supported_target_count"],
                    },
                    sort_keys=True,
                )
            )
        elif arguments.command == "workflow-prepare-e3":
            destination = ensure_e3_preparation_marker(
                run_root=arguments.run_root,
                prepared_dir=arguments.prepared_dir,
                minimum_mean_plddt=arguments.minimum_mean_plddt,
                marker_path=arguments.marker,
            )
            print(json.dumps({"status": "VALID", "marker": str(destination)}, sort_keys=True))
        elif arguments.command == "workflow-stage-e3-review":
            destination = stage_e3_label_review(
                preparation_marker=arguments.preparation_marker,
                reviewed_labels=arguments.reviewed_labels,
                marker_path=arguments.marker,
            )
            print(json.dumps({"status": "VALID", "marker": str(destination)}, sort_keys=True))
        elif arguments.command == "approve-e3-review":
            destination = approve_e3_label_review(
                preparation_marker=arguments.preparation_marker,
                review_marker=arguments.review_marker,
                reviewed_labels=arguments.reviewed_labels,
                approval_marker=arguments.approval_marker,
                curator=arguments.curator,
                note=arguments.note,
                profile=arguments.profile,
                automated_test_mode=arguments.automated_test_mode,
                automated_test_marker=arguments.automated_test_marker,
                automated_evidence_mode=arguments.automated_evidence_mode,
                evidence_label_marker=arguments.evidence_label_marker,
            )
            print(
                json.dumps(
                    {
                        "status": (
                            "AUTOMATED_TEST_ONLY" if arguments.automated_test_mode else "APPROVED"
                        ),
                        "approval_marker": str(destination),
                    },
                    sort_keys=True,
                )
            )
        elif arguments.command == "workflow-verify-e3-review":
            destination = verify_e3_label_review(
                preparation_marker=arguments.preparation_marker,
                review_marker=arguments.review_marker,
                reviewed_labels=arguments.reviewed_labels,
                approval_marker=arguments.approval_marker,
                marker_path=arguments.marker,
                profile=arguments.profile,
            )
            print(json.dumps({"status": "VALID", "marker": str(destination)}, sort_keys=True))
        elif arguments.command == "workflow-initialise-e3":
            destination = ensure_e3_campaign_marker(
                preparation_marker=arguments.preparation_marker,
                review_verification_marker=arguments.review_verification_marker,
                campaign_config=arguments.campaign_config,
                campaign_id=arguments.campaign_id,
                profile=arguments.profile,
                marker_path=arguments.marker,
            )
            print(json.dumps({"status": "VALID", "marker": str(destination)}, sort_keys=True))
        elif arguments.command == "initialise":
            destination = initialise_campaign(
                config_path=arguments.config,
                campaign_id=arguments.campaign_id,
                profile=arguments.profile,
                sequences_fasta=arguments.sequences_fasta,
                label_assignments=arguments.label_assignments,
                label_evidence_marker=arguments.label_evidence_marker,
                label_evidence_audit=arguments.label_evidence_audit,
                control_matching_audit=arguments.control_matching_audit,
                label_definition_features=arguments.label_definition_features,
                class_labelling_summary=arguments.class_labelling_summary,
                unresolved_assignments=arguments.unresolved_assignments,
                features=arguments.features,
                domains=arguments.domains,
                redundancy_clusters=arguments.redundancy_clusters,
                structures=arguments.structures,
                structure_comparisons=arguments.structure_comparisons,
                structural_alignment_resource=arguments.structural_alignment_resource,
                alphafold_accessions=arguments.alphafold_accessions,
                orthofinder_resource=arguments.orthofinder_resource,
                orthofinder_results=arguments.orthofinder_results,
                orthofinder_group_type=arguments.orthofinder_group_type,
                orthofinder_hierarchy_node=arguments.orthofinder_hierarchy_node,
                orthofinder_run_id=arguments.orthofinder_run_id,
                enable_alphafold=arguments.enable_alphafold,
                enable_foldseek=arguments.enable_foldseek,
                foldseek_maximum_hits=arguments.foldseek_maximum_hits,
            )
            print(json.dumps({"status": "CREATED", "config": str(destination)}))
        elif arguments.command == "run-all":
            destination = run_campaign(
                config_path=arguments.config,
                output_dir=arguments.output_dir,
                threads=arguments.threads,
                resume=arguments.resume,
            )
            print(json.dumps({"status": "COMPLETE", "result_dir": str(destination)}))
        elif arguments.command == "validate":
            print(json.dumps(validate_campaign(config_path=arguments.config), sort_keys=True))
        elif arguments.command == "verify":
            verify_completed_result(result_dir=arguments.resource)
            print(
                json.dumps(
                    {
                        "status": "VALID",
                        "result_dir": str(arguments.resource.expanduser().resolve()),
                    },
                    sort_keys=True,
                )
            )
        elif arguments.command == "workflow-validate":
            destination = publish_validation_marker(
                config_path=arguments.config,
                marker_path=arguments.marker,
            )
            print(json.dumps({"status": "VALID", "marker": str(destination)}, sort_keys=True))
        elif arguments.command == "workflow-verify":
            destination = publish_verification_marker(
                result_dir=arguments.resource,
                marker_path=arguments.marker,
            )
            print(json.dumps({"status": "VALID", "marker": str(destination)}, sort_keys=True))
        else:
            profile = load_profile(source=arguments.profile)
            print(
                json.dumps(
                    {
                        "profile_id": profile.profile_id,
                        "profile_version": profile.profile_version,
                        "display_name": profile.display_name,
                        "require_structural_evidence": (profile.require_structural_evidence),
                        "default_comparison": {
                            "target_root_label_id": (profile.default_target_root_label_id),
                            "background_label_id": profile.default_background_label_id,
                            "excluded_label_ids": list(profile.default_excluded_label_ids),
                            "excluded_subtree_label_ids": list(
                                profile.default_excluded_subtree_label_ids
                            ),
                            "description": profile.default_comparison_description,
                        },
                        "labels": [
                            {
                                **vars(label),
                                "aliases": list(label.aliases),
                            }
                            for label in profile.labels
                        ],
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        return 0
    except ProteinSignatureError as error:
        LOGGER.error("%s", error)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        LOGGER.warning("Interrupted by user")
        return 130
