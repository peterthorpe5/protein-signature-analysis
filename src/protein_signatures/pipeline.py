"""End-to-end orchestration for a reproducible protein-signature campaign."""

from __future__ import annotations

import logging
import shutil
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Any

from . import __version__
from .alphafold import acquire_alphafold_models, read_alphafold_requests
from .assessment import compose_feature_assessment_universes
from .associations import analyse_feature_associations, summarise_signatures
from .checksums import sha256_json
from .config import config_to_record, load_config
from .errors import ExternalToolError, InputValidationError, PublicationError
from .evidence_labels import (
    CLASS_SUMMARY_FIELDS,
    CONTROL_MATCH_FIELDS,
    EVIDENCE_AUDIT_FIELDS,
    LABEL_DEFINITION_FIELDS,
    UNRESOLVED_FIELDS,
    read_evidence_audit_table,
    verify_evidence_label_bundle,
)
from .explainable_ml import run_explainable_models
from .fasta import read_protein_fasta
from .feature_provenance import (
    ALL_DATA_EXPLORATORY,
    prepare_imported_feature_evidence,
    sequence_cohort_sha256,
)
from .foldseek import foldseek_version, run_foldseek_all_vs_all
from .io_utils import read_json
from .models import (
    AlphaFoldRequest,
    CampaignConfig,
    FeatureRecord,
    FoldEvidenceStatus,
    FoldseekSettings,
    PairwiseStructureComparison,
    ProteinProfile,
    StructureAnalysisEligibility,
    StructureComparisonStatus,
    StructureRecord,
)
from .orthofinder import discover_orthofinder_layout, read_group_memberships
from .orthofinder_resource import (
    discover_orthofinder_resource,
    read_resource_group_context,
    read_resource_memberships,
    resource_input_paths,
)
from .partitions import assign_partitions
from .profiles import (
    default_profile_comparisons,
    expand_positive_memberships,
    load_profile,
    validate_assignment_profile_compatibility,
    validate_comparison_labels,
)
from .publication import publish_result, verify_completed_result, verify_input_authorities
from .redundancy import derive_exact_sequence_clusters, read_redundancy_clusters
from .reporting import build_human_reports
from .sequence_signatures import build_kmer_features
from .structural_resource import import_structural_alignment_resource
from .structural_signatures import (
    derive_structure_assessment_universes,
    derive_structure_features,
)
from .tables import (
    complete_domain_assessments,
    derive_domain_evidence,
    read_domains,
    read_features,
    read_label_assignments,
    read_structure_comparisons,
    read_structures,
)

LOGGER = logging.getLogger(__name__)


def run_campaign(
    *, config_path: Path, output_dir: Path, threads: int = 1, resume: bool = False
) -> Path:
    """Run every enabled analysis stage and atomically publish one result.

    Args:
        config_path: Campaign YAML.
        output_dir: New result directory.
        threads: Positive worker budget recorded for provenance.
        resume: Verify and reuse a completed immutable result.

    Returns:
        Published result directory.

    Raises:
        InputValidationError: If the worker budget or any input is invalid.
    """

    if threads < 1:
        raise InputValidationError("threads must be a positive integer.")
    config = load_config(path=config_path)
    profile = load_profile(source=config.profile_name)
    config = _with_effective_comparisons(config=config, profile=profile)
    validate_comparison_labels(comparisons=config.comparisons, profile=profile)
    run_identity = _run_identity(config=config, profile=profile)
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        if not resume:
            raise PublicationError(
                f"Output already exists and will not be overwritten: {destination}. "
                "Use --resume only for a completed checksum-valid result."
            )
        verify_completed_result(result_dir=destination)
        previous_metadata = read_json(path=destination / "run_metadata.json")
        if (
            not isinstance(previous_metadata, dict)
            or previous_metadata.get("run_identity_sha256") != run_identity
        ):
            raise InputValidationError(
                "Existing result was created from a different configuration, profile "
                "or package version."
            )
        verify_input_authorities(result_dir=destination)
        LOGGER.info("Reused complete checksum-compatible campaign result %s", destination)
        return destination
    data = _prepare_campaign(config=config, profile=profile, threads=threads)
    metadata = {
        "schema_version": 1,
        "package": "protein-signature-analysis",
        "package_version": __version__,
        "run_identity_sha256": run_identity,
        "campaign": config_to_record(config=config),
        "profile": {
            "profile_id": profile.profile_id,
            "profile_version": profile.profile_version,
            "display_name": profile.display_name,
            "require_structural_evidence": profile.require_structural_evidence,
            "default_target_root_label_id": profile.default_target_root_label_id,
            "default_background_label_id": profile.default_background_label_id,
            "default_comparison_description": profile.default_comparison_description,
            "default_excluded_label_ids": list(profile.default_excluded_label_ids),
            "default_excluded_subtree_label_ids": list(profile.default_excluded_subtree_label_ids),
        },
        "threads": threads,
        "counts": {table_name: len(records) for table_name, records in data["tables"].items()},
        "evidence_availability": data["evidence_availability"],
        "profile_structural_evidence": data["profile_structural_evidence"],
        "orthofinder": data["orthofinder"],
        "foldseek": data["foldseek"],
        "imported_structural_alignment": data["imported_structural_alignment"],
        "explainable_ml": data["explainable_ml"],
        "human_reports": data["human_reports"],
        "automated_label_evidence": data["automated_label_evidence"],
        "determinism": {
            "partition_seed": config.analysis.random_seed,
            "stable_sorting": True,
            "atomic_publication": True,
        },
    }
    return publish_result(
        output_dir=output_dir,
        tables=data["tables"],
        metadata=metadata,
        input_paths=data["input_paths"],
        asset_sources=data["asset_sources"],
        resume=resume,
    )


def validate_campaign(*, config_path: Path) -> dict[str, Any]:
    """Validate campaign configuration and local inputs without remote acquisition.

    Args:
        config_path: Campaign YAML.

    Returns:
        JSON-compatible validation summary.
    """

    config = load_config(path=config_path)
    profile = load_profile(source=config.profile_name)
    config = _with_effective_comparisons(config=config, profile=profile)
    validate_comparison_labels(comparisons=config.comparisons, profile=profile)
    sequences = read_protein_fasta(path=config.inputs.sequences_fasta)
    protein_ids = frozenset(item.protein_id for item in sequences)
    assignments = read_label_assignments(
        path=config.inputs.label_assignments,
        protein_ids=protein_ids,
        label_ids=profile.label_ids(),
    )
    validate_assignment_profile_compatibility(assignments=assignments, profile=profile)
    _, label_evidence_metadata = _load_label_evidence_tables(config=config)
    external_feature_records = (
        read_features(path=config.inputs.features, sequences=sequences)
        if config.inputs.features is not None
        else ()
    )
    if config.inputs.domains is not None:
        domain_hits, domain_assessments = read_domains(
            path=config.inputs.domains, sequences=sequences
        )
    else:
        domain_hits, domain_assessments = (), ()
    redundancy_clusters = (
        read_redundancy_clusters(
            path=config.inputs.redundancy_clusters,
            protein_ids=protein_ids,
        )
        if config.inputs.redundancy_clusters is not None
        else ()
    )
    structures = (
        read_structures(path=config.inputs.structures, sequences=sequences)
        if config.inputs.structures is not None
        else ()
    )
    structure_comparisons = (
        read_structure_comparisons(
            path=config.inputs.structure_comparisons,
            protein_ids=protein_ids,
        )
        if config.inputs.structure_comparisons is not None
        else ()
    )
    alphafold_requests = (
        read_alphafold_requests(
            path=config.inputs.alphafold_accessions,
            protein_ids=protein_ids,
        )
        if config.inputs.alphafold_accessions is not None
        else ()
    )
    if config.alphafold.enabled and config.inputs.alphafold_accessions is None:
        raise InputValidationError(
            "alphafold.enabled is true but inputs.alphafold_accessions is not configured."
        )
    foldseek_preflight = (
        _validate_foldseek_preflight(
            settings=config.foldseek,
            structures=structures,
            alphafold_requests=(alphafold_requests if config.alphafold.enabled else ()),
        )
        if config.foldseek.enabled
        else {
            "status": "NOT_SELECTED",
            "tool_version": "",
            "candidate_model_count": 0,
            "candidate_protein_count": 0,
            "maximum_hits": config.foldseek.maximum_hits,
        }
    )
    imported_structural_evidence = (
        import_structural_alignment_resource(
            resource_dir=config.inputs.structural_alignment_resource,
            sequences=sequences,
        )
        if config.inputs.structural_alignment_resource is not None
        else None
    )
    structural_policy = _validate_required_structural_evidence(
        profile=profile,
        structures=structures,
        comparisons=(
            *structure_comparisons,
            *(
                imported_structural_evidence.comparisons
                if imported_structural_evidence is not None
                else ()
            ),
        ),
        imported_feature_count=(
            len(imported_structural_evidence.features)
            if imported_structural_evidence is not None
            else 0
        ),
        pending_alphafold_request_count=(
            len(alphafold_requests) if config.alphafold.enabled and config.foldseek.enabled else 0
        ),
        pending_structural_search_count=(1 if foldseek_preflight["status"] == "VALID" else 0),
    )
    orthofinder_source = None
    if config.inputs.orthofinder_resource is not None:
        orthofinder_source = discover_orthofinder_resource(
            resource_dir=config.inputs.orthofinder_resource
        )
        _validate_requested_resource_run(config=config, resource=orthofinder_source)
        read_resource_memberships(
            resource=orthofinder_source,
            group_type=config.inputs.orthofinder_group_type,
            hierarchy_node=config.inputs.orthofinder_hierarchy_node,
            protein_ids=protein_ids,
        )
    elif config.inputs.orthofinder_results is not None:
        source_run_id = config.inputs.orthofinder_run_id or config.campaign_id
        orthofinder_source = discover_orthofinder_layout(
            results_dir=config.inputs.orthofinder_results
        )
        read_group_memberships(
            layout=orthofinder_source,
            run_id=source_run_id,
            group_type=config.inputs.orthofinder_group_type,
            hierarchy_node=config.inputs.orthofinder_hierarchy_node,
            protein_ids=protein_ids,
        )
    return {
        "status": "VALID",
        "campaign_id": config.campaign_id,
        "profile_id": profile.profile_id,
        "protein_count": len(sequences),
        "assignment_count": len(assignments),
        "external_feature_record_count": len(external_feature_records),
        "domain_hit_count": len(domain_hits),
        "domain_assessment_count": len(domain_assessments),
        "redundancy_membership_count": len(redundancy_clusters),
        "structure_count": len(structures),
        "structure_comparison_count": len(structure_comparisons),
        "alphafold_request_count": len(alphafold_requests),
        "comparison_count": len(config.comparisons),
        "alphafold_enabled": config.alphafold.enabled,
        "foldseek_preflight": foldseek_preflight,
        "explainable_ml_enabled": config.explainable_ml.enabled,
        "orthofinder_version": (
            orthofinder_source.orthofinder_version
            if config.inputs.orthofinder_resource is not None
            else orthofinder_source.version
            if orthofinder_source is not None
            else ""
        ),
        "orthofinder_source_mode": (
            "ORTHOFINDER_RESULTS_RESOURCE"
            if config.inputs.orthofinder_resource is not None
            else orthofinder_source.source_mode
            if orthofinder_source is not None
            else ""
        ),
        "structural_alignment_resource": (
            "VALIDATED"
            if config.inputs.structural_alignment_resource is not None
            else "NOT_SELECTED"
        ),
        "profile_structural_evidence": structural_policy,
        "automated_label_evidence": label_evidence_metadata,
    }


def _validate_required_structural_evidence(
    *,
    profile: ProteinProfile,
    structures: tuple[StructureRecord, ...],
    comparisons: tuple[PairwiseStructureComparison, ...],
    imported_feature_count: int = 0,
    pending_alphafold_request_count: int = 0,
    pending_structural_search_count: int = 0,
    completed_structural_search_count: int = 0,
) -> dict[str, Any]:
    """Enforce a profile's minimum structural-evidence policy.

    Reviewed fold assignments and completed imported comparisons count as
    structural evidence without requiring local coordinate files. During
    read-only validation, a successfully preflighted structural search can be
    pending; the runtime repeats this guard after acquisition and search.

    Args:
        profile: Selected classification profile.
        structures: Parsed supplied and, at runtime, acquired model records.
        comparisons: Parsed imported or generated structural comparisons.
        imported_feature_count: Structural features from a verified resource.
        pending_alphafold_request_count: Requests feeding a validated pending search.
        pending_structural_search_count: Validated but not yet run structural search.
        completed_structural_search_count: Completed searches, including zero-hit runs.

    Returns:
        Auditable policy status and eligible evidence counts.

    Raises:
        InputValidationError: If counts are negative or a required profile has
            neither eligible evidence nor a pending acquisition route.
    """

    counts = (
        imported_feature_count,
        pending_alphafold_request_count,
        pending_structural_search_count,
        completed_structural_search_count,
    )
    if any(count < 0 for count in counts):
        raise InputValidationError("Structural evidence counts cannot be negative.")
    eligible_coordinate_models = sum(
        structure.is_coordinate_analysis_eligible for structure in structures
    )
    fold_eligible_states = {
        StructureAnalysisEligibility.ELIGIBLE,
        StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE,
    }
    eligible_fold_assignments = sum(
        bool(structure.fold_id)
        and structure.fold_evidence_status == FoldEvidenceStatus.ASSESSED_WITH_HIT
        and structure.analysis_eligibility_status in fold_eligible_states
        for structure in structures
    )
    complete_comparisons = sum(
        comparison.comparison_status
        in {
            StructureComparisonStatus.COMPLETE,
            StructureComparisonStatus.SUCCESS,
            StructureComparisonStatus.PASS,
        }
        for comparison in comparisons
    )
    eligible_evidence_items = (
        eligible_fold_assignments
        + complete_comparisons
        + imported_feature_count
        + completed_structural_search_count
    )
    if not profile.require_structural_evidence:
        status = "NOT_REQUIRED"
    elif eligible_evidence_items:
        status = "COMPLETE"
    elif pending_structural_search_count:
        status = "PENDING_STRUCTURAL_SEARCH"
    else:
        raise InputValidationError(
            f"Profile {profile.profile_id!r} requires structural evidence, but no "
            "assessed fold hit, completed structural comparison, verified imported "
            "structural feature or completed structural search was available. "
            "A coordinate model alone is readiness evidence, not an analysis result."
        )
    result = {
        "required": profile.require_structural_evidence,
        "status": status,
        "eligible_coordinate_model_count": eligible_coordinate_models,
        "eligible_fold_assignment_count": eligible_fold_assignments,
        "complete_structure_comparison_count": complete_comparisons,
        "imported_structural_feature_count": imported_feature_count,
        "pending_alphafold_request_count": pending_alphafold_request_count,
        "pending_structural_search_count": pending_structural_search_count,
        "completed_structural_search_count": completed_structural_search_count,
        "eligible_evidence_item_count": eligible_evidence_items,
    }
    LOGGER.info(
        "Validated profile structural policy status=%s eligible_items=%d pending_searches=%d",
        status,
        eligible_evidence_items,
        pending_structural_search_count,
    )
    return result


def _validate_foldseek_preflight(
    *,
    settings: FoldseekSettings,
    structures: tuple[StructureRecord, ...],
    alphafold_requests: tuple[AlphaFoldRequest, ...],
) -> dict[str, Any]:
    """Validate a Foldseek plan without executing a structural search.

    The candidate total is deliberately conservative: it combines supplied,
    locally available and analysis-eligible coordinates with every selected
    AlphaFold request that could yield another model. The runtime repeats the
    exact guard after acquisition using the models that actually qualify.

    Args:
        settings: Validated Foldseek settings.
        structures: Parsed supplied structure records.
        alphafold_requests: Requests selected for remote acquisition later.

    Returns:
        JSON-compatible executable, version and candidate-count evidence.

    Raises:
        ExternalToolError: If the executable or its version is unavailable.
        InputValidationError: If the candidate universe is too small or the
            maximum-hit ceiling could truncate an all-versus-all search.
    """

    executable = shutil.which(settings.executable)
    if executable is None:
        raise ExternalToolError(
            f"Foldseek executable {settings.executable!r} was not found; install it "
            "or disable foldseek."
        )
    tool_version = foldseek_version(executable=executable)
    supplied_candidates = tuple(
        structure for structure in structures if structure.is_coordinate_analysis_eligible
    )
    candidate_count = len(supplied_candidates) + len(alphafold_requests)
    candidate_proteins = {item.protein_id for item in (*supplied_candidates, *alphafold_requests)}
    if len(candidate_proteins) < 2:
        raise InputValidationError(
            "Foldseek is enabled but preflight found fewer than two potential "
            "proteins with analysis-eligible coordinate models."
        )
    if settings.maximum_hits < candidate_count:
        raise InputValidationError(
            "foldseek.maximum_hits must be at least the preflight candidate model "
            f"count ({candidate_count}) so the all-versus-all search cannot be "
            "silently truncated."
        )
    LOGGER.info(
        "Validated Foldseek %s preflight with %d potential models and max-seqs %d",
        tool_version,
        candidate_count,
        settings.maximum_hits,
    )
    return {
        "status": "VALID",
        "executable": executable,
        "tool_version": tool_version,
        "supplied_candidate_model_count": len(supplied_candidates),
        "alphafold_candidate_model_count": len(alphafold_requests),
        "candidate_model_count": candidate_count,
        "candidate_protein_count": len(candidate_proteins),
        "maximum_hits": settings.maximum_hits,
    }


def _prepare_campaign(
    *, config: CampaignConfig, profile: ProteinProfile, threads: int
) -> dict[str, Any]:
    """Load inputs, derive features and assemble canonical result tables.

    Args:
        config: Validated campaign configuration.
        profile: Selected protein classification profile.
        threads: Positive structural worker budget.

    Returns:
        Tables, input authorities and evidence-availability metadata.
    """

    sequences = read_protein_fasta(path=config.inputs.sequences_fasta)
    protein_ids = frozenset(item.protein_id for item in sequences)
    assignments = read_label_assignments(
        path=config.inputs.label_assignments,
        protein_ids=protein_ids,
        label_ids=profile.label_ids(),
    )
    validate_assignment_profile_compatibility(assignments=assignments, profile=profile)
    label_evidence_tables, label_evidence_metadata = _load_label_evidence_tables(config=config)
    exact_redundancy = derive_exact_sequence_clusters(sequences=sequences)
    supplied_redundancy = (
        read_redundancy_clusters(
            path=config.inputs.redundancy_clusters,
            protein_ids=protein_ids,
        )
        if config.inputs.redundancy_clusters is not None
        else ()
    )
    label_memberships = expand_positive_memberships(assignments=assignments, profile=profile)
    external_feature_records = (
        read_features(path=config.inputs.features, sequences=sequences)
        if config.inputs.features is not None
        else ()
    )
    if config.inputs.domains is not None:
        domain_hits, supplied_domain_assessments = read_domains(
            path=config.inputs.domains, sequences=sequences
        )
    else:
        domain_hits, supplied_domain_assessments = (), ()
    domain_assessments = complete_domain_assessments(
        protein_ids=protein_ids, assessments=supplied_domain_assessments
    )
    domain_features, domain_sequences = derive_domain_evidence(
        hits=domain_hits, sequences=sequences
    )
    supplied_structures = (
        read_structures(path=config.inputs.structures, sequences=sequences)
        if config.inputs.structures is not None
        else ()
    )
    requests = (
        read_alphafold_requests(path=config.inputs.alphafold_accessions, protein_ids=protein_ids)
        if config.inputs.alphafold_accessions is not None
        else ()
    )
    acquisitions, acquired_structures = acquire_alphafold_models(
        requests=requests, sequences=sequences, settings=config.alphafold
    )
    structures = _merge_structures(supplied=supplied_structures, acquired=acquired_structures)
    supplied_structure_comparisons = (
        read_structure_comparisons(
            path=config.inputs.structure_comparisons, protein_ids=protein_ids
        )
        if config.inputs.structure_comparisons is not None
        else ()
    )
    imported_structural_evidence = (
        import_structural_alignment_resource(
            resource_dir=config.inputs.structural_alignment_resource,
            sequences=sequences,
        )
        if config.inputs.structural_alignment_resource is not None
        else None
    )
    imported_structure_comparisons = (
        imported_structural_evidence.comparisons if imported_structural_evidence is not None else ()
    )
    foldseek_evidence = None
    generated_structure_comparisons = ()
    if config.foldseek.enabled:
        foldseek_evidence = run_foldseek_all_vs_all(
            structures=structures, settings=config.foldseek, threads=threads
        )
        generated_structure_comparisons = foldseek_evidence.comparisons
    structure_comparisons = _merge_structure_comparisons(
        supplied=(*supplied_structure_comparisons, *imported_structure_comparisons),
        generated=generated_structure_comparisons,
    )
    structural_policy = _validate_required_structural_evidence(
        profile=profile,
        structures=structures,
        comparisons=structure_comparisons,
        imported_feature_count=(
            len(imported_structural_evidence.features)
            if imported_structural_evidence is not None
            else 0
        ),
        completed_structural_search_count=(1 if foldseek_evidence is not None else 0),
    )
    orthofinder_source = None
    memberships = ()
    group_context: tuple[dict[str, Any], ...] = ()
    if config.inputs.orthofinder_resource is not None:
        orthofinder_source = discover_orthofinder_resource(
            resource_dir=config.inputs.orthofinder_resource
        )
        _validate_requested_resource_run(config=config, resource=orthofinder_source)
        memberships = read_resource_memberships(
            resource=orthofinder_source,
            group_type=config.inputs.orthofinder_group_type,
            hierarchy_node=config.inputs.orthofinder_hierarchy_node,
            protein_ids=protein_ids,
        )
        group_context = read_resource_group_context(
            resource=orthofinder_source, memberships=memberships
        )
    elif config.inputs.orthofinder_results is not None:
        orthofinder_source = discover_orthofinder_layout(
            results_dir=config.inputs.orthofinder_results
        )
        memberships = read_group_memberships(
            layout=orthofinder_source,
            run_id=config.inputs.orthofinder_run_id or config.campaign_id,
            group_type=config.inputs.orthofinder_group_type,
            hierarchy_node=config.inputs.orthofinder_hierarchy_node,
            protein_ids=protein_ids,
        )
    partitions = assign_partitions(
        sequences=sequences,
        memberships=memberships,
        validation_fraction=config.analysis.validation_fraction,
        random_seed=config.analysis.random_seed,
        redundancy_memberships=supplied_redundancy,
    )
    discovery_protein_ids = frozenset(
        item.protein_id for item in partitions if item.partition == "DISCOVERY"
    )
    prepared_external = prepare_imported_feature_evidence(
        records=external_feature_records,
        sequences=sequences,
        discovery_protein_ids=discovery_protein_ids,
    )
    structure_features, structure_clusters = derive_structure_features(
        structures=structures,
        comparisons=structure_comparisons,
        tm_score_threshold=config.analysis.structural_tm_score_threshold,
        minimum_coverage=config.analysis.structural_minimum_coverage,
        partitions=partitions,
        derivation_cohort_sha256=sequence_cohort_sha256(
            sequences=sequences,
            protein_ids=discovery_protein_ids,
        ),
    )
    kmer_features = build_kmer_features(
        sequences=sequences,
        discovery_protein_ids=discovery_protein_ids,
        lengths=config.analysis.kmer_lengths,
        minimum_proteins=config.analysis.minimum_feature_proteins,
        maximum_features=config.analysis.maximum_kmer_features,
    )
    features = _merge_features(
        external_features=prepared_external.confirmatory_features,
        domain_features=domain_features,
        structure_features=structure_features,
        kmer_features=kmer_features,
    )
    comparison_universe_members = dict(
        imported_structural_evidence.comparison_universe_members
        if imported_structural_evidence is not None
        else {}
    )
    if foldseek_evidence is not None:
        existing_members = comparison_universe_members.get(foldseek_evidence.comparison_universe_id)
        if (
            existing_members is not None
            and existing_members != foldseek_evidence.assessment_universe
        ):
            raise InputValidationError(
                "Foldseek and imported structural evidence disagree about comparison "
                f"universe {foldseek_evidence.comparison_universe_id!r}."
            )
        comparison_universe_members[foldseek_evidence.comparison_universe_id] = (
            foldseek_evidence.assessment_universe
        )
    structural_assessment_universes = derive_structure_assessment_universes(
        structures=structures,
        structure_clusters=structure_clusters,
        comparison_universe_members=comparison_universe_members,
        alphafold_acquisitions=acquisitions,
    )
    feature_assessment_universes = compose_feature_assessment_universes(
        protein_ids=protein_ids,
        features=features,
        domain_assessments=domain_assessments,
        explicit_assessment_universes=(
            prepared_external.assessment_universes,
            structural_assessment_universes,
        ),
    )
    exploratory_feature_keys = frozenset(
        (feature.feature_type, feature.feature_id)
        for feature in features
        if feature.derivation_scope == ALL_DATA_EXPLORATORY
    )
    associations = analyse_feature_associations(
        comparisons=config.comparisons,
        label_memberships=label_memberships,
        partitions=partitions,
        features=features,
        settings=config.analysis,
        feature_assessment_universes=feature_assessment_universes,
        exploratory_feature_keys=exploratory_feature_keys,
    )
    signatures = summarise_signatures(
        associations=associations,
        fdr_threshold=config.analysis.fdr_threshold,
        exploratory_feature_keys=exploratory_feature_keys,
    )
    explainable_ml = run_explainable_models(
        comparisons=config.comparisons,
        label_memberships=label_memberships,
        partitions=partitions,
        features=features,
        settings=config.explainable_ml,
        random_seed=config.analysis.random_seed,
        feature_assessment_universes=feature_assessment_universes,
        exploratory_feature_keys=exploratory_feature_keys,
    )
    structure_asset_sources, asset_names = _structure_asset_sources(structures=structures)
    tables = {
        "alphafold_acquisitions": tuple(
            _alphafold_record(item=item, asset_names=asset_names) for item in acquisitions
        ),
        "associations": tuple(item.to_record() for item in associations),
        "comparisons": tuple(_comparison_record(item=item) for item in config.comparisons),
        "domain_assessments": tuple(item.to_record() for item in domain_assessments),
        "domain_hits": tuple(item.to_record() for item in domain_hits),
        "domain_sequences": domain_sequences,
        "feature_assessments": tuple(
            item.to_record()
            for item in (
                *prepared_external.assessment_records,
                *(
                    imported_structural_evidence.features
                    if imported_structural_evidence is not None
                    else ()
                ),
            )
        ),
        "features": tuple(item.to_record() for item in features),
        "label_assignments": tuple(item.to_record() for item in assignments),
        **label_evidence_tables,
        "label_memberships": label_memberships,
        "ml_explanations": explainable_ml.explanations,
        "ml_feature_importance": explainable_ml.feature_importance,
        "ml_models": explainable_ml.models,
        "ml_predictions": explainable_ml.predictions,
        "ml_plot_inventory": explainable_ml.plots,
        "orthofinder_memberships": tuple(item.to_record() for item in memberships),
        "orthofinder_group_context": group_context,
        "partitions": tuple(item.to_record() for item in partitions),
        "profile_labels": tuple(_profile_label_record(item=item) for item in profile.labels),
        "proteins": tuple(item.to_record() for item in sequences),
        "redundancy_clusters": tuple(
            item.to_record() for item in (*exact_redundancy, *supplied_redundancy)
        ),
        "signatures": tuple(item.to_record() for item in signatures),
        "imported_structural_group_summaries": (
            imported_structural_evidence.group_summaries
            if imported_structural_evidence is not None
            else ()
        ),
        "structure_clusters": structure_clusters,
        "structure_comparisons": tuple(item.to_record() for item in structure_comparisons),
        "structures": tuple(
            _structure_record(item=item, asset_names=asset_names) for item in structures
        ),
    }
    report_cache = (
        config.explainable_ml.plot_cache_dir.parent
        / "human_reports"
        / _run_identity(config=config, profile=profile)
    )
    human_reports = build_human_reports(
        tables=tables,
        cache_dir=report_cache,
        existing_plot_assets=explainable_ml.plot_assets,
        fdr_threshold=config.analysis.fdr_threshold,
    )
    report_asset_sources = dict(human_reports.assets)
    foldseek_asset_sources = (
        {
            (
                f"assets/foldseek/{foldseek_evidence.cache_key}/foldseek_all_vs_all.tsv"
            ): foldseek_evidence.raw_output_path,
            (
                f"assets/foldseek/{foldseek_evidence.cache_key}/COMPLETED.json"
            ): foldseek_evidence.completion_manifest_path,
        }
        if foldseek_evidence is not None
        else {}
    )
    asset_collections = (
        structure_asset_sources,
        report_asset_sources,
        foldseek_asset_sources,
    )
    all_asset_names = [name for collection in asset_collections for name in collection]
    collisions = {name for name, count in Counter(all_asset_names).items() if count > 1}
    if collisions:
        raise InputValidationError(f"Result asset paths collide: {sorted(collisions)}")
    asset_sources = {
        **structure_asset_sources,
        **report_asset_sources,
        **foldseek_asset_sources,
    }
    input_paths = _input_authorities(
        config=config,
        source=orthofinder_source,
        resource_mode=config.inputs.orthofinder_resource is not None,
    )
    input_paths.extend(
        item.coordinate_path for item in supplied_structures if item.coordinate_path is not None
    )
    if imported_structural_evidence is not None:
        input_paths.extend(imported_structural_evidence.input_paths)
    return {
        "tables": tables,
        "profile_structural_evidence": structural_policy,
        "asset_sources": asset_sources,
        "input_paths": tuple(input_paths),
        "evidence_availability": {
            "domains": "COMPLETE" if config.inputs.domains else "INPUT_UNAVAILABLE",
            "pfam": _pfam_availability(assessments=domain_assessments),
            "structures": "COMPLETE" if structures else "INPUT_UNAVAILABLE",
            "structure_comparisons": ("COMPLETE" if structure_comparisons else "INPUT_UNAVAILABLE"),
            "imported_structural_alignment": (
                "COMPLETE" if imported_structural_evidence else "NOT_SELECTED"
            ),
            "foldseek": "COMPLETE" if foldseek_evidence else "NOT_SELECTED",
            "orthofinder": ("COMPLETE" if orthofinder_source else "INPUT_UNAVAILABLE"),
            "near_redundancy": ("COMPLETE" if supplied_redundancy else "INPUT_UNAVAILABLE"),
            "alphafold_acquisition": ("COMPLETE" if acquisitions else "NOT_SELECTED"),
            "explainable_ml": explainable_ml.status,
            "automated_label_evidence": label_evidence_metadata["status"],
        },
        "orthofinder": (
            orthofinder_source.to_record()
            if orthofinder_source
            else {"status": "INPUT_UNAVAILABLE"}
        ),
        "foldseek": (
            {
                "status": "COMPLETE",
                "tool_version": foldseek_evidence.tool_version,
                "cache_key": foldseek_evidence.cache_key,
                "cache_reused": foldseek_evidence.reused,
                "comparison_universe_id": foldseek_evidence.comparison_universe_id,
                "assessment_universe_size": len(foldseek_evidence.assessment_universe),
                "raw_output_path": (
                    f"assets/foldseek/{foldseek_evidence.cache_key}/foldseek_all_vs_all.tsv"
                ),
                "completion_manifest_path": (
                    f"assets/foldseek/{foldseek_evidence.cache_key}/COMPLETED.json"
                ),
            }
            if foldseek_evidence
            else {"status": "NOT_SELECTED"}
        ),
        "imported_structural_alignment": (
            {
                "status": "COMPLETE",
                "package_version": imported_structural_evidence.package_version,
                "run_digest": imported_structural_evidence.run_digest,
                "comparison_universe_count": len(
                    imported_structural_evidence.comparison_universe_members
                ),
            }
            if imported_structural_evidence
            else {"status": "NOT_SELECTED"}
        ),
        "explainable_ml": {
            "status": explainable_ml.status,
            "implementation_version": explainable_ml.implementation_version,
            "shap_version": explainable_ml.shap_version,
            "matplotlib_version": explainable_ml.matplotlib_version,
            "model_count": len(explainable_ml.models),
            "plot_count": len(explainable_ml.plots),
            "explanation_method": "SHAP_LINEAR_INTERVENTIONAL_LOG_ODDS",
        },
        "human_reports": {
            "status": "COMPLETE",
            "report_file_count": len(human_reports.assets),
            "logical_figure_count": human_reports.figure_count,
            "excel_workbook_count": human_reports.workbook_count,
            "inventory_row_count": len(human_reports.inventory),
        },
        "automated_label_evidence": label_evidence_metadata,
    }


def _merge_structures(
    *, supplied: tuple[StructureRecord, ...], acquired: tuple[StructureRecord, ...]
) -> tuple[StructureRecord, ...]:
    """Merge structure sources while protecting identifier uniqueness.

    Args:
        supplied: User or cluster supplied structures.
        acquired: AlphaFold DB structures acquired by this run.

    Returns:
        Ordered combined structure inventory.

    Raises:
        InputValidationError: If structure identifiers collide.
    """

    rows: dict[str, StructureRecord] = {}
    for structure in (*supplied, *acquired):
        if structure.structure_id in rows:
            raise InputValidationError(
                f"Structure identifier occurs in multiple sources: {structure.structure_id!r}"
            )
        rows[structure.structure_id] = structure
    return tuple(rows[key] for key in sorted(rows))


def _merge_structure_comparisons(
    *,
    supplied: tuple[PairwiseStructureComparison, ...],
    generated: tuple[PairwiseStructureComparison, ...],
) -> tuple[PairwiseStructureComparison, ...]:
    """Merge imported and generated structural comparisons uniquely.

    Args:
        supplied: Imported pairwise evidence.
        generated: Foldseek evidence generated for this campaign.

    Returns:
        Ordered comparison rows.

    Raises:
        InputValidationError: If a source record identifier collides.
    """

    rows: dict[tuple[str, str], PairwiseStructureComparison] = {}
    for comparison in (*supplied, *generated):
        key = (comparison.comparison_tool, comparison.source_record_id)
        if key in rows:
            raise InputValidationError(f"Structural source record occurs more than once: {key!r}")
        rows[key] = comparison
    return tuple(
        rows[key]
        for key in sorted(
            rows,
            key=lambda item: (
                rows[item].protein_a_id,
                rows[item].protein_b_id,
                item[0],
                item[1],
            ),
        )
    )


def _merge_features(
    *,
    external_features: tuple[FeatureRecord, ...],
    domain_features: tuple[FeatureRecord, ...],
    structure_features: tuple[FeatureRecord, ...],
    kmer_features: tuple[FeatureRecord, ...],
) -> tuple[FeatureRecord, ...]:
    """Merge usable features deterministically and remove exact duplicates.

    Args:
        external_features: Imported general features.
        domain_features: Derived domain and architecture features.
        structure_features: Derived fold and structure features.
        kmer_features: Derived amino-acid features.

    Returns:
        Ordered usable feature evidence.
    """

    excluded_states = {"FAILED", "NOT_ASSESSED", "INPUT_UNAVAILABLE", "EXCLUDED"}
    rows: dict[tuple[Any, ...], FeatureRecord] = {}
    for feature in (
        *external_features,
        *domain_features,
        *structure_features,
        *kmer_features,
    ):
        if feature.evidence_status.upper() in excluded_states:
            continue
        key = (
            feature.protein_id,
            feature.feature_type,
            feature.feature_id,
            feature.start,
            feature.end,
            feature.evidence_source,
            feature.evidence_reference,
        )
        rows[key] = feature
    return tuple(
        rows[key]
        for key in sorted(
            rows,
            key=lambda item: (
                str(item[0]),
                str(item[1]),
                str(item[2]),
                int(item[3] or 0),
                int(item[4] or 0),
                str(item[5]),
                str(item[6]),
            ),
        )
    )


def _structure_asset_sources(
    *, structures: tuple[StructureRecord, ...]
) -> tuple[dict[str, Path], dict[str, str]]:
    """Map local coordinates to deterministic portable result asset names.

    Args:
        structures: Combined structure inventory.

    Returns:
        Relative asset-to-source mapping and structure-to-relative-name mapping.
    """

    sources: dict[str, Path] = {}
    names: dict[str, str] = {}
    for structure in structures:
        if structure.coordinate_path is None:
            continue
        source = Path(structure.coordinate_path)
        lower_name = source.name.lower()
        if lower_name.endswith(".cif.gz"):
            suffix = ".cif.gz"
        elif lower_name.endswith(".pdb.gz"):
            suffix = ".pdb.gz"
        elif lower_name.endswith(".cif"):
            suffix = ".cif"
        else:
            suffix = ".pdb"
        relative = f"assets/structures/{structure.coordinate_sha256}{suffix}"
        sources[relative] = source
        names[structure.structure_id] = relative
    return sources, names


def _structure_record(*, item: StructureRecord, asset_names: dict[str, str]) -> dict[str, Any]:
    """Serialise a structure with a portable result-relative coordinate path.

    Args:
        item: Structure record.
        asset_names: Structure identifiers mapped to portable result paths.

    Returns:
        Structure table row.
    """

    row = item.to_record()
    row["coordinate_path"] = asset_names.get(item.structure_id, "")
    return row


def _alphafold_record(*, item: Any, asset_names: dict[str, str]) -> dict[str, Any]:
    """Serialise an AlphaFold outcome with a portable model path.

    Args:
        item: AlphaFold acquisition record.
        asset_names: Structure identifiers mapped to portable result paths.

    Returns:
        AlphaFold acquisition table row.
    """

    row = item.to_record()
    row["coordinate_path"] = asset_names.get(item.structure_id, "")
    return row


def _comparison_record(*, item: Any) -> dict[str, str]:
    """Serialise one comparison definition for tabular publication.

    Args:
        item: Comparison definition.

    Returns:
        Flat comparison row.
    """

    return {
        "comparison_id": item.comparison_id,
        "display_name": item.display_name,
        "target_label_ids": "|".join(item.target_label_ids),
        "background_label_ids": "|".join(item.background_label_ids),
        "description": item.description,
    }


def _profile_label_record(*, item: Any) -> dict[str, Any]:
    """Serialise one profile label for tabular publication.

    Args:
        item: Profile label.

    Returns:
        Flat profile-label row.
    """

    return {
        "label_id": item.label_id,
        "display_name": item.display_name,
        "parent_label_id": item.parent_label_id,
        "level": item.level,
        "system_class": item.system_class,
        "mechanistic_class": item.mechanistic_class,
        "component_role": item.component_role,
        "family": item.family,
        "active_site_expected": item.active_site_expected,
        "active_site_residue": item.active_site_residue,
        "default_analysis": item.default_analysis,
        "default_background_label_id": item.default_background_label_id,
        "assignment_exclusivity_group": item.assignment_exclusivity_group,
        "reviewed_positive_allowed": item.reviewed_positive_allowed,
        "aliases": "|".join(item.aliases),
        "description": item.description,
    }


def _load_label_evidence_tables(
    *, config: CampaignConfig
) -> tuple[dict[str, tuple[dict[str, Any], ...]], dict[str, Any]]:
    """Load and cross-check an optional automated label-evidence bundle.

    Args:
        config: Validated campaign configuration.

    Returns:
        Canonical audit tables and a metadata summary. All audit tables are
        present but empty when automated evidence labelling was not selected.

    Raises:
        InputValidationError: If the bundle is partial or does not match the
            configured label/domain authorities.
    """

    table_contract = {
        "label_evidence_audit": (
            config.inputs.label_evidence_audit,
            EVIDENCE_AUDIT_FIELDS,
        ),
        "control_matching_audit": (
            config.inputs.control_matching_audit,
            CONTROL_MATCH_FIELDS,
        ),
        "label_definition_features": (
            config.inputs.label_definition_features,
            LABEL_DEFINITION_FIELDS,
        ),
        "class_labelling_summary": (
            config.inputs.class_labelling_summary,
            CLASS_SUMMARY_FIELDS,
        ),
        "unresolved_assignments": (
            config.inputs.unresolved_assignments,
            UNRESOLVED_FIELDS,
        ),
    }
    configured_paths = {path for path, _ in table_contract.values() if path is not None}
    marker = config.inputs.label_evidence_marker
    if marker is None and configured_paths:
        raise InputValidationError(
            "Automated label audit tables require inputs.label_evidence_marker."
        )
    if marker is None:
        return (
            {name: () for name in table_contract},
            {"status": "NOT_SELECTED", "interpretation_scope": "NOT_APPLICABLE"},
        )
    missing = [name for name, (path, _) in table_contract.items() if path is None]
    if missing:
        raise InputValidationError(
            f"An automated evidence bundle requires every audit table; missing={missing}."
        )
    marker_path = Path(marker).resolve()
    if marker_path.name != "EVIDENCE_LABELS.json":
        raise InputValidationError(
            "inputs.label_evidence_marker must be named EVIDENCE_LABELS.json."
        )
    document = verify_evidence_label_bundle(bundle_dir=marker_path.parent)
    expected_paths = {
        "label_assignments": config.inputs.label_assignments,
        **{name: path for name, (path, _) in table_contract.items()},
    }
    for field, configured in expected_paths.items():
        if configured is None:
            raise InputValidationError(
                f"Evidence bundle requires a configured {field!r} authority."
            )
        recorded = Path(str(document.get(field, marker_path.parent / f"{field}.tsv"))).resolve()
        if field not in document:
            recorded = (marker_path.parent / f"{field}.tsv").resolve()
        if recorded != Path(configured).resolve():
            raise InputValidationError(
                f"Evidence bundle field {field!r} differs from campaign configuration."
            )
    recorded_domains = document.get("analysis_domains")
    configured_domains = config.inputs.domains
    if recorded_domains is None:
        if configured_domains is not None:
            raise InputValidationError(
                "Campaign domains were configured but the evidence bundle has no "
                "label-definition-safe domain projection."
            )
    elif (
        configured_domains is None
        or Path(str(recorded_domains)).resolve() != Path(configured_domains).resolve()
    ):
        raise InputValidationError(
            "Evidence bundle analysis_domains differs from campaign configuration."
        )
    tables = {
        name: read_evidence_audit_table(path=path, fields=fields)
        for name, (path, fields) in table_contract.items()
    }
    return (
        tables,
        {
            "status": "PROVISIONAL_EVIDENCE_SUPPORTED",
            "interpretation_scope": document["interpretation_scope"],
            "human_review_completed": document["human_review_completed"],
            "ruleset_id": document["ruleset_id"],
            "ruleset_version": document["ruleset_version"],
            "ruleset_sha256": document["ruleset_sha256"],
            "evidence_supported_target_count": document["evidence_supported_target_count"],
            "matched_control_protein_count": document["matched_control_protein_count"],
            "warning": document["warning"],
        },
    )


def _input_authorities(*, config: CampaignConfig, source: Any, resource_mode: bool) -> list[Path]:
    """Collect local scientific authorities for the checksum manifest.

    Args:
        config: Campaign configuration.
        source: Optional raw layout or published-resource identity.
        resource_mode: Whether ``source`` is a published resource.

    Returns:
        Existing input file paths.
    """

    paths = [
        value
        for value in vars(config.inputs).values()
        if isinstance(value, Path) and value.is_file()
    ]
    paths.append(config.config_path)
    profile_path = Path(config.profile_name).expanduser()
    if profile_path.is_file():
        paths.append(profile_path.resolve())
    if source is not None and resource_mode:
        paths.extend(resource_input_paths(resource=source))
    elif source is not None:
        if source.log_path is not None:
            paths.append(source.log_path)
        paths.extend(source.completion_authority_paths)
        if source.orthogroups_path is not None:
            paths.append(source.orthogroups_path)
        paths.extend(source.hog_paths)
        if source.sequence_ids_path is not None:
            paths.append(source.sequence_ids_path)
        if source.species_ids_path is not None:
            paths.append(source.species_ids_path)
    return list(dict.fromkeys(paths))


def _validate_requested_resource_run(*, config: CampaignConfig, resource: Any) -> None:
    """Reject an optional configured run identity that disagrees with a resource.

    Args:
        config: Campaign configuration.
        resource: Validated upstream resource identity.

    Raises:
        InputValidationError: If an explicit run identity differs.
    """

    requested = config.inputs.orthofinder_run_id
    if requested and requested != resource.run_id:
        raise InputValidationError(
            "inputs.orthofinder.run_id does not match the published resource: "
            f"{requested!r} versus {resource.run_id!r}."
        )


def _pfam_availability(*, assessments: tuple[Any, ...]) -> str:
    """Summarise Pfam coverage without conflating absence and non-assessment.

    Args:
        assessments: Complete domain assessment rows.

    Returns:
        Controlled campaign-level coverage state.
    """

    statuses = {
        item.assessment_status.value
        for item in assessments
        if item.domain_authority.casefold() == "pfam"
    }
    assessed = statuses & {"ASSESSED_WITH_HIT", "ASSESSED_NO_HIT"}
    if not assessed:
        return "NOT_ASSESSED"
    if statuses <= {"ASSESSED_WITH_HIT", "ASSESSED_NO_HIT"}:
        return "COMPLETE"
    return "PARTIAL"


def _run_identity(*, config: CampaignConfig, profile: ProteinProfile) -> str:
    """Calculate a stable configuration/profile/package run identity.

    Args:
        config: Validated campaign configuration.
        profile: Validated classification profile.

    Returns:
        SHA-256 identity used to protect resume semantics.
    """

    return sha256_json(
        value={
            "package_version": __version__,
            "configuration": config_to_record(config=config),
            "profile": {
                "profile_id": profile.profile_id,
                "profile_version": profile.profile_version,
                "display_name": profile.display_name,
                "require_structural_evidence": profile.require_structural_evidence,
                "default_target_root_label_id": (profile.default_target_root_label_id),
                "default_background_label_id": profile.default_background_label_id,
                "default_comparison_description": (profile.default_comparison_description),
                "default_excluded_label_ids": list(profile.default_excluded_label_ids),
                "default_excluded_subtree_label_ids": list(
                    profile.default_excluded_subtree_label_ids
                ),
                "labels": [
                    {
                        **vars(label),
                        "aliases": list(label.aliases),
                    }
                    for label in profile.labels
                ],
            },
        }
    )


def _with_effective_comparisons(
    *, config: CampaignConfig, profile: ProteinProfile
) -> CampaignConfig:
    """Resolve profile-default comparisons while preserving explicit campaigns.

    Args:
        config: Parsed campaign configuration.
        profile: Selected classification profile.

    Returns:
        Configuration carrying explicit effective comparisons.
    """

    if config.comparisons:
        return config
    comparisons = default_profile_comparisons(profile=profile)
    LOGGER.info("Expanded profile defaults into %d separate comparisons", len(comparisons))
    return replace(config, comparisons=comparisons)
