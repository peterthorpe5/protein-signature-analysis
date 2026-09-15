"""Create a complete campaign configuration from explicit scientific inputs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

from .config import load_config
from .errors import InputValidationError, ProteinSignatureError, PublicationError
from .io_utils import write_text_atomic
from .validation import validate_identifier, validate_text

LOGGER = logging.getLogger(__name__)


def initialise_campaign(
    *,
    config_path: Path,
    campaign_id: str,
    profile: str,
    sequences_fasta: Path,
    label_assignments: Path,
    label_evidence_marker: Path | None = None,
    label_evidence_audit: Path | None = None,
    control_matching_audit: Path | None = None,
    label_definition_features: Path | None = None,
    class_labelling_summary: Path | None = None,
    unresolved_assignments: Path | None = None,
    features: Path | None = None,
    domains: Path | None = None,
    redundancy_clusters: Path | None = None,
    structures: Path | None = None,
    structure_comparisons: Path | None = None,
    structural_alignment_resource: Path | None = None,
    alphafold_accessions: Path | None = None,
    orthofinder_resource: Path | None = None,
    orthofinder_results: Path | None = None,
    orthofinder_group_type: str = "HOG",
    orthofinder_hierarchy_node: str = "N0",
    orthofinder_run_id: str = "",
    enable_alphafold: bool = False,
    enable_foldseek: bool = False,
    foldseek_maximum_hits: int = 1000,
) -> Path:
    """Write and re-read a ready-to-run campaign YAML file.

    Args:
        config_path: New YAML destination.
        campaign_id: Stable campaign identifier.
        profile: Built-in profile name or profile YAML path.
        sequences_fasta: Protein FASTA authority.
        label_assignments: Reviewed label-assignment TSV authority.
        label_evidence_marker: Optional automated evidence-bundle marker.
        label_evidence_audit: Optional complete label-decision audit TSV.
        control_matching_audit: Optional matched-control audit TSV.
        label_definition_features: Optional circularity-exclusion ledger TSV.
        class_labelling_summary: Optional class-coverage summary TSV.
        unresolved_assignments: Optional abstention and conflict review TSV.
        features: Optional externally derived feature TSV.
        domains: Optional domain-hit and assessment TSV.
        redundancy_clusters: Optional near-redundancy membership TSV.
        structures: Optional structure inventory TSV.
        structure_comparisons: Optional pairwise structure-comparison TSV.
        structural_alignment_resource: Optional completed structural resource.
        alphafold_accessions: Optional protein-to-UniProt accession TSV.
        orthofinder_resource: Preferred completed ``orthofinder-results`` resource.
        orthofinder_results: Raw completed OrthoFinder 2.5.5 or 3 results fallback.
        orthofinder_group_type: ``HOG`` or ``LEGACY_ORTHOGROUP``.
        orthofinder_hierarchy_node: HOG hierarchy node, normally ``N0``.
        orthofinder_run_id: Optional explicit upstream run identifier.
        enable_alphafold: Download requested AlphaFold Database models during analysis.
        enable_foldseek: Run Foldseek over available coordinate models.
        foldseek_maximum_hits: Maximum Foldseek matches retained per query.

    Returns:
        Absolute path to the validated campaign configuration.

    Raises:
        InputValidationError: If supplied authorities are missing or incompatible.
        PublicationError: If the destination already exists.
    """

    destination = Path(config_path).expanduser().resolve()
    if destination.exists():
        raise PublicationError(f"Campaign configuration already exists: {destination}")
    if orthofinder_resource is not None and orthofinder_results is not None:
        raise InputValidationError(
            "Supply either a published orthofinder-results resource or raw OrthoFinder "
            "results, not both."
        )
    if enable_alphafold and alphafold_accessions is None:
        raise InputValidationError(
            "AlphaFold acquisition requires an explicit protein-to-accession TSV."
        )
    if (
        not isinstance(foldseek_maximum_hits, int)
        or isinstance(foldseek_maximum_hits, bool)
        or foldseek_maximum_hits < 1
    ):
        raise InputValidationError("foldseek_maximum_hits must be a positive integer.")
    campaign = validate_identifier(value=campaign_id, field_name="campaign_id")
    profile_source = _profile_source(value=profile)
    file_inputs = {
        "sequences_fasta": _input_file(path=sequences_fasta, field="sequences_fasta"),
        "label_assignments": _input_file(path=label_assignments, field="label_assignments"),
        "label_evidence_marker": _optional_input_file(
            path=label_evidence_marker, field="label_evidence_marker"
        ),
        "label_evidence_audit": _optional_input_file(
            path=label_evidence_audit, field="label_evidence_audit"
        ),
        "control_matching_audit": _optional_input_file(
            path=control_matching_audit, field="control_matching_audit"
        ),
        "label_definition_features": _optional_input_file(
            path=label_definition_features, field="label_definition_features"
        ),
        "class_labelling_summary": _optional_input_file(
            path=class_labelling_summary, field="class_labelling_summary"
        ),
        "unresolved_assignments": _optional_input_file(
            path=unresolved_assignments, field="unresolved_assignments"
        ),
        "features": _optional_input_file(path=features, field="features"),
        "domains": _optional_input_file(path=domains, field="domains"),
        "redundancy_clusters": _optional_input_file(
            path=redundancy_clusters, field="redundancy_clusters"
        ),
        "structures": _optional_input_file(path=structures, field="structures"),
        "structure_comparisons": _optional_input_file(
            path=structure_comparisons, field="structure_comparisons"
        ),
        "alphafold_accessions": _optional_input_file(
            path=alphafold_accessions, field="alphafold_accessions"
        ),
    }
    structural_resource = _optional_input_directory(
        path=structural_alignment_resource,
        field="structural_alignment_resource",
    )
    orthofinder_published = _optional_input_directory(
        path=orthofinder_resource,
        field="orthofinder_resource",
    )
    orthofinder_raw = _optional_input_directory(
        path=orthofinder_results,
        field="orthofinder_results",
    )
    group_type = validate_text(
        value=orthofinder_group_type,
        field_name="orthofinder_group_type",
    ).upper()
    if group_type not in {"HOG", "LEGACY_ORTHOGROUP"}:
        raise InputValidationError("orthofinder_group_type must be HOG or LEGACY_ORTHOGROUP.")
    hierarchy_node = validate_text(
        value=orthofinder_hierarchy_node,
        field_name="orthofinder_hierarchy_node",
        allow_empty=group_type == "LEGACY_ORTHOGROUP",
    )
    if group_type == "LEGACY_ORTHOGROUP" and hierarchy_node:
        raise InputValidationError(
            "orthofinder_hierarchy_node must be blank for LEGACY_ORTHOGROUP."
        )
    document = _campaign_document(
        campaign_id=campaign,
        profile=profile_source,
        file_inputs=file_inputs,
        structural_alignment_resource=structural_resource,
        orthofinder_resource=orthofinder_published,
        orthofinder_results=orthofinder_raw,
        orthofinder_group_type=group_type,
        orthofinder_hierarchy_node=hierarchy_node,
        orthofinder_run_id=orthofinder_run_id,
        enable_alphafold=enable_alphafold,
        enable_foldseek=enable_foldseek,
        foldseek_maximum_hits=foldseek_maximum_hits,
        config_parent=destination.parent,
    )
    content = yaml.safe_dump(document, sort_keys=False, allow_unicode=True)
    write_text_atomic(path=destination, text=content)
    try:
        load_config(path=destination)
    except ProteinSignatureError:
        destination.unlink(missing_ok=True)
        raise
    LOGGER.info("Created validated campaign configuration %s", destination)
    return destination


def _campaign_document(
    *,
    campaign_id: str,
    profile: str,
    file_inputs: dict[str, str | None],
    structural_alignment_resource: str | None,
    orthofinder_resource: str | None,
    orthofinder_results: str | None,
    orthofinder_group_type: str,
    orthofinder_hierarchy_node: str,
    orthofinder_run_id: str,
    enable_alphafold: bool,
    enable_foldseek: bool,
    foldseek_maximum_hits: int,
    config_parent: Path,
) -> dict[str, Any]:
    """Build one serialisable schema-1 campaign document.

    Args:
        campaign_id: Validated campaign identifier.
        profile: Profile source.
        file_inputs: Resolved optional and required file authorities.
        structural_alignment_resource: Optional structural-resource directory.
        orthofinder_resource: Optional published OrthoFinder resource.
        orthofinder_results: Optional raw OrthoFinder result directory.
        orthofinder_group_type: Selected grouping authority.
        orthofinder_hierarchy_node: Selected hierarchy node.
        orthofinder_run_id: Optional run identifier.
        enable_alphafold: Whether AlphaFold acquisition is enabled.
        enable_foldseek: Whether Foldseek is enabled.
        foldseek_maximum_hits: Maximum Foldseek matches retained per query.
        config_parent: Parent used for local cache directories.

    Returns:
        Complete campaign mapping.
    """

    cache_root = config_parent / ".protein_signature_cache"
    return {
        "schema_version": 1,
        "campaign": {"campaign_id": campaign_id, "profile": profile},
        "inputs": {
            **file_inputs,
            "structural_alignment_resource": structural_alignment_resource,
            "orthofinder": {
                "resource_dir": orthofinder_resource,
                "results_dir": orthofinder_results,
                "group_type": orthofinder_group_type,
                "hierarchy_node": orthofinder_hierarchy_node,
                "run_id": validate_text(
                    value=orthofinder_run_id,
                    field_name="orthofinder_run_id",
                    allow_empty=True,
                ),
            },
        },
        "alphafold": {
            "enabled": enable_alphafold,
            "cache_dir": str(cache_root / "alphafold"),
            "timeout_seconds": 60,
            "retries": 3,
            "minimum_mean_plddt": 50,
        },
        "foldseek": {
            "enabled": enable_foldseek,
            "executable": "foldseek",
            "cache_dir": str(cache_root / "foldseek"),
            "e_value_threshold": 0.001,
            "sensitivity": 9.5,
            "maximum_hits": foldseek_maximum_hits,
        },
        "explainable_ml": {
            "enabled": True,
            "minimum_samples_per_class": 10,
            "minimum_groups_per_class": 3,
            "maximum_features": 500,
            "regularisation_strengths": [0.01, 0.1, 1.0],
            "l1_ratio": 0.2,
            "cross_validation_folds": 5,
            "maximum_iterations": 2000,
            "convergence_tolerance": 0.0001,
            "permutation_repeats": 10,
            "maximum_permutation_features": 50,
            "top_local_explanations": 10,
            "maximum_waterfall_plots": 8,
            "plot_cache_dir": str(cache_root / "shap"),
            "exclude_technical_features": True,
        },
        "analysis": {
            "kmer_lengths": [3, 4],
            "minimum_target_proteins": 3,
            "minimum_background_proteins": 3,
            "minimum_feature_proteins": 2,
            "maximum_kmer_features": 250000,
            "validation_fraction": 0.2,
            "fdr_threshold": 0.05,
            "random_seed": 1729,
            "structural_tm_score_threshold": 0.5,
            "structural_minimum_coverage": 0.5,
        },
        "comparisons": "profile_defaults",
    }


def _input_file(*, path: Path, field: str) -> str:
    """Resolve one required non-empty input file.

    Args:
        path: Input path.
        field: Diagnostic field name.

    Returns:
        Absolute path text.

    Raises:
        InputValidationError: If the authority is missing or empty.
    """

    candidate = Path(path).expanduser().resolve()
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"{field} is missing or empty: {candidate}")
    return str(candidate)


def _profile_source(*, value: str) -> str:
    """Resolve a built-in profile name or custom YAML before changing directory.

    Args:
        value: Built-in profile name or user-supplied YAML path.

    Returns:
        Unchanged built-in name or absolute custom-profile path.

    Raises:
        InputValidationError: If a path-like profile is missing or empty.
    """

    source = validate_text(value=value, field_name="profile")
    if "/" not in source and "\\" not in source and not source.endswith((".yaml", ".yml")):
        return source
    candidate = Path(source).expanduser().resolve()
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"Custom profile YAML is missing or empty: {candidate}")
    return str(candidate)


def _optional_input_file(*, path: Path | None, field: str) -> str | None:
    """Resolve an optional input file.

    Args:
        path: Optional path.
        field: Diagnostic field name.

    Returns:
        Absolute path text or ``None``.
    """

    return None if path is None else _input_file(path=path, field=field)


def _optional_input_directory(*, path: Path | None, field: str) -> str | None:
    """Resolve an optional existing input directory.

    Args:
        path: Optional directory.
        field: Diagnostic field name.

    Returns:
        Absolute path text or ``None``.

    Raises:
        InputValidationError: If the supplied directory is absent.
    """

    if path is None:
        return None
    candidate = Path(path).expanduser().resolve()
    if not candidate.is_dir():
        raise InputValidationError(f"{field} directory does not exist: {candidate}")
    return str(candidate)
