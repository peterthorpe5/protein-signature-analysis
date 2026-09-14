"""Command-line interface for validation, analysis and profile inspection."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .catalogue import prepare_catalogue
from .e3_workflow_bridge import prepare_e3_workflow_inputs
from .errors import ProteinSignatureError
from .logging_config import configure_logging
from .pipeline import run_campaign, validate_campaign
from .profiles import load_profile
from .publication import verify_completed_result
from .starter import initialise_campaign

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
        elif arguments.command == "initialise":
            destination = initialise_campaign(
                config_path=arguments.config,
                campaign_id=arguments.campaign_id,
                profile=arguments.profile,
                sequences_fasta=arguments.sequences_fasta,
                label_assignments=arguments.label_assignments,
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
