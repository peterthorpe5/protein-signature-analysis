"""Auditable evidence-led protein labelling and matched-control selection.

This module is deliberately separate from the synthetic smoke-test generator.
It can create provisional, analysis-eligible assignments only when explicit
evidence rules are satisfied.  It never forces complete classification: weak,
conflicting or structurally unsuitable records remain proposed, ambiguous or
unmapped and cannot enter target or background cohorts.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import logging
import math
import os
import re
import shutil
import tempfile
from collections import Counter, defaultdict
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from . import __version__
from .checksums import sha256_file, sha256_text
from .errors import ConfigurationError, InputValidationError, PublicationError
from .exports import dataframe_to_xlsx_bytes
from .fasta import read_protein_fasta
from .io_utils import iter_tsv, read_json, write_bytes_atomic, write_json_atomic, write_tsv_atomic
from .models import (
    CurationStatus,
    GroupMembership,
    LabelAssignment,
    ProteinProfile,
    SequenceRecord,
    StructureAnalysisEligibility,
)
from .orthofinder import discover_orthofinder_layout, read_group_memberships
from .orthofinder_resource import (
    discover_orthofinder_resource,
    read_resource_memberships,
    resource_input_paths,
)
from .partitions import assign_partitions
from .profiles import default_profile_comparisons, label_ancestors, load_profile
from .redundancy import read_redundancy_clusters
from .tables import (
    DOMAIN_FIELDS,
    LABEL_FIELDS,
    read_domains,
    read_label_assignments,
    read_structures,
)
from .validation import validate_identifier, validate_text

LOGGER = logging.getLogger(__name__)

EVIDENCE_POSITIVE_STATUS = "EVIDENCE_SUPPORTED_POSITIVE"
EVIDENCE_TARGET_STATUS = "AUTOMATED_EVIDENCE_SUPPORTED"
EVIDENCE_CONTROL_STATUS = "AUTOMATED_MATCHED_CONTROL"
EVIDENCE_APPROVER = "AUTOMATED_EVIDENCE_PIPELINE"
EVIDENCE_BUNDLE_STATUS = "PROVISIONAL_EVIDENCE_LABELS_COMPLETE"
EVIDENCE_INTERPRETATION_SCOPE = "PROVISIONAL_HYPOTHESIS_GENERATION"

EXTERNAL_ANNOTATION_FIELDS = (
    "protein_id",
    "label_id",
    "annotation_text",
    "evidence_status",
    "evidence_source",
    "evidence_reference",
    "annotation_scope",
)
REVIEW_CONTEXT_FIELDS = (
    "protein_id",
    "species",
    "input_candidate_states",
    "sequence_length",
    "pfam_accessions",
    "upstream_e3_families",
    "upstream_evidence_roles",
    "annotation_status_details",
    "structure_analysis_eligible",
)
PROTEIN_METADATA_FIELDS = (
    "protein_id",
    "species",
    "input_candidate",
)
SEED_CATALOGUE_FIELDS = (
    "seed_id",
    "seed_category",
    "seed_review_status",
    "associated_seed_categories",
    "associated_seed_review_statuses",
    "protein_sequence",
    "annotation_scope",
    "catalogue_source",
)
EVIDENCE_AUDIT_FIELDS = (
    "protein_id",
    "label_id",
    "rule_id",
    "decision",
    "confidence_tier",
    "score",
    "evidence_group_count",
    "evidence_groups",
    "evidence_items",
    "independence_unit",
    "conflicting_label_ids",
    "reason",
)
CONTROL_MATCH_FIELDS = (
    "background_label_id",
    "target_label_ids",
    "target_protein_id",
    "target_unit_id",
    "control_protein_id",
    "control_unit_id",
    "species_match",
    "structure_eligibility_match",
    "length_log2_difference",
    "domain_count_difference",
    "domain_architecture_jaccard",
    "mean_confidence_difference",
    "match_score",
    "status",
    "reason",
)
LABEL_DEFINITION_FIELDS = (
    "label_id",
    "feature_type",
    "feature_id",
    "feature_name",
    "rule_id",
    "evidence_role",
    "exclusion_scope",
    "reason",
)
CLASS_SUMMARY_FIELDS = (
    "label_id",
    "label_type",
    "direct_positive_protein_count",
    "independent_unit_count",
    "matched_background_label_id",
    "matched_control_protein_count",
    "matched_control_unit_count",
    "status",
)
UNRESOLVED_FIELDS = (
    "protein_id",
    "provisional_label_id",
    "curation_status",
    "best_score",
    "candidate_label_ids",
    "reason",
)

_TRUSTED_EXTERNAL_STATUSES = frozenset({"CURATED", "EXPERIMENTAL", "MANUALLY_REVIEWED", "REVIEWED"})
_TRUE_VALUES = frozenset({"1", "TRUE", "YES"})
_ELIGIBLE_STRUCTURE_STATES = frozenset(
    {
        StructureAnalysisEligibility.ELIGIBLE,
        StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE,
    }
)
_AUDIT_INTEGER_FIELDS = frozenset(
    {
        "evidence_group_count",
        "direct_positive_protein_count",
        "independent_unit_count",
        "matched_control_protein_count",
        "matched_control_unit_count",
        "domain_count_difference",
    }
)
_AUDIT_FLOAT_FIELDS = frozenset(
    {
        "score",
        "length_log2_difference",
        "domain_architecture_jaccard",
        "mean_confidence_difference",
        "match_score",
        "best_score",
    }
)
_AUDIT_BOOLEAN_FIELDS = frozenset({"species_match", "structure_eligibility_match"})


@dataclass(frozen=True)
class EvidenceSettings:
    """Validated scoring, propagation and background-matching policy."""

    annotation_score: float
    specific_annotation_score: float
    pfam_score: float
    supporting_pfam_score: float
    orthology_score: float
    trusted_label_score: float
    minimum_acceptance_score: float
    minimum_evidence_groups: int
    minimum_score_margin: float
    require_target_structure_eligible: bool
    orthology_propagation_enabled: bool
    minimum_orthology_anchor_proteins: int
    control_units_per_target_unit: int
    minimum_control_units_per_background: int
    require_species_match: bool
    require_structure_match: bool
    maximum_log2_length_difference: float
    maximum_domain_count_difference: int
    exclude_input_candidates_from_controls: bool


@dataclass(frozen=True)
class EvidenceRule:
    """One profile-label inference rule compiled from a ruleset YAML."""

    rule_id: str
    label_id: str
    component_role: str
    priority: int
    annotation_patterns: tuple[re.Pattern[str], ...]
    specific_annotation_patterns: tuple[re.Pattern[str], ...]
    exclusion_patterns: tuple[re.Pattern[str], ...]
    required_domain_patterns: tuple[re.Pattern[str], ...]
    supporting_domain_patterns: tuple[re.Pattern[str], ...]
    allow_specific_annotation_only: bool
    propagate_by_orthology: bool


@dataclass(frozen=True)
class EvidenceRuleSet:
    """Complete evidence-labelling ruleset for one profile version."""

    ruleset_id: str
    ruleset_version: str
    profile_id: str
    settings: EvidenceSettings
    rules: tuple[EvidenceRule, ...]
    source_path: Path


@dataclass(frozen=True)
class TextEvidence:
    """One annotation text item with an explicit independence group."""

    text: str
    evidence_group: str
    evidence_source: str
    evidence_reference: str
    reliability: float
    trusted_status: bool


@dataclass(frozen=True)
class ProteinEvidenceContext:
    """Evidence and matching covariates for one authoritative protein."""

    protein_id: str
    sequence_length: int
    sequence_sha256: str
    species: tuple[str, ...]
    domain_ids: tuple[str, ...]
    domain_names: Mapping[str, str]
    domain_count: int
    pfam_assessed: bool
    structure_eligible: bool
    mean_confidence: float | None
    input_candidate: bool
    text_evidence: tuple[TextEvidence, ...]
    direct_labels: tuple[LabelAssignment, ...]
    independence_unit: str


@dataclass(frozen=True)
class EvidenceItem:
    """One scored item contributing to a candidate label."""

    evidence_group: str
    evidence_kind: str
    value: str
    score: float
    evidence_source: str
    evidence_reference: str
    feature_type: str = ""
    feature_id: str = ""
    feature_name: str = ""

    def to_record(self) -> dict[str, Any]:
        """Return a JSON-compatible evidence record.

        Returns:
            Plain mapping suitable for deterministic JSON serialisation.
        """

        return {
            "evidence_group": self.evidence_group,
            "evidence_kind": self.evidence_kind,
            "value": self.value,
            "score": self.score,
            "evidence_source": self.evidence_source,
            "evidence_reference": self.evidence_reference,
            "feature_type": self.feature_type,
            "feature_id": self.feature_id,
            "feature_name": self.feature_name,
        }


@dataclass(frozen=True)
class LabelCandidate:
    """Scored provisional assignment for one protein and profile label."""

    protein_id: str
    label_id: str
    rule_id: str
    component_role: str
    score: float
    evidence_groups: tuple[str, ...]
    evidence_items: tuple[EvidenceItem, ...]
    accepted: bool
    confidence_tier: str
    reason: str
    priority: int
    trusted_direct: bool = False


@dataclass(frozen=True)
class CandidateResolution:
    """Conservative resolution of all candidate labels for one protein."""

    selected: LabelCandidate | None
    candidates: tuple[LabelCandidate, ...]
    decisions: Mapping[tuple[str, str], str]
    status: str
    reason: str
    conflicting_label_ids: tuple[str, ...]


@dataclass(frozen=True)
class GroupEvidence:
    """Loaded homology memberships and their exact source authorities."""

    memberships: tuple[GroupMembership, ...]
    input_paths: tuple[Path, ...]
    source_mode: str
    version: str


def built_in_evidence_rules_path(*, rules_name: str) -> Path:
    """Resolve one packaged evidence-rules YAML by short name.

    Args:
        rules_name: Packaged rule identifier, for example ``e3``.

    Returns:
        Absolute ruleset YAML path.

    Raises:
        ConfigurationError: If the packaged ruleset does not exist.
    """

    name = validate_identifier(value=rules_name, field_name="evidence rules name")
    candidate = files("protein_signatures").joinpath("data", "evidence_rules", f"{name}.yaml")
    if not candidate.is_file():
        raise ConfigurationError(f"Unknown built-in evidence rules: {name!r}")
    return Path(str(candidate)).resolve()


def load_evidence_rules(*, source: str | Path, profile: ProteinProfile) -> EvidenceRuleSet:
    """Load and strictly validate one evidence-labelling ruleset.

    Args:
        source: Built-in short name or custom YAML path.
        profile: Profile whose labels and roles constrain the rules.

    Returns:
        Immutable compiled ruleset.

    Raises:
        ConfigurationError: If YAML is malformed.
        InputValidationError: If fields, values or label references are invalid.
    """

    path = (
        built_in_evidence_rules_path(rules_name=source)
        if isinstance(source, str) and "/" not in source and "\\" not in source
        else Path(source).expanduser().resolve()
    )
    if not path.is_file() or path.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty evidence-rules YAML: {path}")
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        raise ConfigurationError(f"Could not parse evidence-rules YAML {path}: {error}") from error
    if not isinstance(document, Mapping):
        raise InputValidationError("Evidence-rules YAML must contain a mapping.")
    allowed = {
        "schema_version",
        "ruleset_id",
        "ruleset_version",
        "profile_id",
        "settings",
        "label_rules",
    }
    _reject_unknown_fields(value=document, allowed=allowed, field_name="evidence_rules")
    if document.get("schema_version") != 1:
        raise ConfigurationError("Evidence-rules schema_version must be 1.")
    profile_id = validate_identifier(
        value=document.get("profile_id"), field_name="evidence_rules.profile_id"
    )
    if profile_id != profile.profile_id:
        raise InputValidationError(
            f"Evidence rules target profile {profile_id!r}, not {profile.profile_id!r}."
        )
    settings = _parse_evidence_settings(value=document.get("settings"))
    raw_rules = document.get("label_rules")
    if isinstance(raw_rules, (str, bytes)) or not isinstance(raw_rules, Sequence) or not raw_rules:
        raise InputValidationError("evidence_rules.label_rules must be a non-empty list.")
    labels_by_id = {label.label_id: label for label in profile.labels}
    rules = tuple(
        _parse_evidence_rule(value=value, index=index, labels_by_id=labels_by_id)
        for index, value in enumerate(raw_rules)
    )
    rule_ids = [rule.rule_id for rule in rules]
    label_ids = [rule.label_id for rule in rules]
    if len(rule_ids) != len(set(rule_ids)):
        raise InputValidationError("Evidence rule identifiers must be unique.")
    if len(label_ids) != len(set(label_ids)):
        raise InputValidationError("Each profile label may have at most one evidence rule.")
    result = EvidenceRuleSet(
        ruleset_id=validate_identifier(
            value=document.get("ruleset_id"), field_name="evidence_rules.ruleset_id"
        ),
        ruleset_version=validate_text(
            value=document.get("ruleset_version"),
            field_name="evidence_rules.ruleset_version",
        ),
        profile_id=profile_id,
        settings=settings,
        rules=rules,
        source_path=path,
    )
    LOGGER.info(
        "Loaded evidence rules %s version %s with %d label rules",
        result.ruleset_id,
        result.ruleset_version,
        len(result.rules),
    )
    return result


def create_evidence_label_bundle(
    *,
    sequences_fasta: Path,
    output_dir: Path,
    profile: str | Path,
    evidence_rules: str | Path,
    protein_metadata: Path | None = None,
    review_context: Path | None = None,
    domains: Path | None = None,
    structures: Path | None = None,
    template_labels: Path | None = None,
    seed_assignments: Path | None = None,
    seed_catalogue: Path | None = None,
    external_annotations: Path | None = None,
    orthofinder_resource: Path | None = None,
    orthofinder_results: Path | None = None,
    orthofinder_group_type: str = "HOG",
    orthofinder_hierarchy_node: str = "N0",
    orthofinder_run_id: str = "evidence_labelling",
    redundancy_clusters: Path | None = None,
    random_seed: int = 1729,
    validation_fraction: float = 0.2,
) -> Path:
    """Create a provisional evidence-led target and matched-control bundle.

    Accepted target calls require trusted direct labels or the ruleset's minimum
    independent corroboration.  Orthology propagation is one-generation only,
    requires unanimous compatible anchors and the expected domain architecture.
    Backgrounds are selected without using outcome signatures and are matched at
    the independent homology/redundancy-block level.  Weak or conflicting records
    remain non-positive.

    Args:
        sequences_fasta: Authoritative protein FASTA.
        output_dir: New atomically published evidence bundle directory.
        profile: Built-in or custom classification profile.
        evidence_rules: Built-in or custom evidence-rules YAML.
        protein_metadata: Optional generic protein metadata TSV containing
            ``protein_id``, pipe-delimited ``species`` and ``input_candidate``.
        review_context: Optional per-protein predecessor review context TSV.
        domains: Optional explicit Pfam/domain assessment TSV.
        structures: Optional structure inventory used only for eligibility and matching.
        template_labels: Optional complete all-protein starter assignment TSV.
        seed_assignments: Optional trusted, directly curated profile-label TSV.
        seed_catalogue: Optional E3 seed catalogue. Cluster-associated annotations
            remain weak context and are never treated as direct seed labels.
        external_annotations: Optional independent annotation authority TSV.
        orthofinder_resource: Optional published ``orthofinder-results`` resource.
        orthofinder_results: Optional completed OrthoFinder 2.5.5 or 3 results.
        orthofinder_group_type: ``HOG`` or ``LEGACY_ORTHOGROUP``.
        orthofinder_hierarchy_node: HOG hierarchy node, normally ``N0``.
        orthofinder_run_id: Stable identifier for raw OrthoFinder membership.
        redundancy_clusters: Optional near-redundancy membership TSV.
        random_seed: Non-negative deterministic partition seed.
        validation_fraction: Downstream held-out fraction used to freeze blocks.

    Returns:
        Absolute path to ``EVIDENCE_LABELS.json`` inside the published bundle.

    Raises:
        InputValidationError: If evidence or matching contracts are unsatisfied.
        PublicationError: If the destination exists or cannot be published safely.
    """

    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise PublicationError(f"Evidence-label bundle already exists: {destination}")
    if destination == Path(sequences_fasta).expanduser().resolve():
        raise PublicationError("Evidence-label output cannot replace its sequence authority.")
    if isinstance(random_seed, bool) or not isinstance(random_seed, int) or random_seed < 0:
        raise InputValidationError("random_seed must be a non-negative integer.")
    if not 0.0 <= validation_fraction <= 0.9:
        raise InputValidationError("validation_fraction must be between 0.0 and 0.9.")
    if orthofinder_resource is not None and orthofinder_results is not None:
        raise InputValidationError(
            "Supply either orthofinder_resource or orthofinder_results, not both."
        )

    sequences_path = _required_file(path=sequences_fasta, field_name="sequences_fasta")
    sequences = read_protein_fasta(path=sequences_path)
    protein_ids = frozenset(record.protein_id for record in sequences)
    loaded_profile = load_profile(source=profile)
    ruleset = load_evidence_rules(source=evidence_rules, profile=loaded_profile)
    domain_path = _optional_file(path=domains, field_name="domains")
    structure_path = _optional_file(path=structures, field_name="structures")
    review_path = _optional_file(path=review_context, field_name="review_context")
    metadata_path = _optional_file(path=protein_metadata, field_name="protein_metadata")
    template_path = _optional_file(path=template_labels, field_name="template_labels")
    seed_path = _optional_file(path=seed_assignments, field_name="seed_assignments")
    catalogue_path = _optional_file(path=seed_catalogue, field_name="seed_catalogue")
    external_path = _optional_file(path=external_annotations, field_name="external_annotations")

    domains_by_protein, domain_assessed, domain_rows = _load_domain_context(
        path=domain_path,
        sequences=sequences,
    )
    structures_by_protein = _load_structure_context(
        path=structure_path,
        sequences=sequences,
    )
    review_by_protein = _load_review_context(
        path=review_path,
        sequences=sequences,
        domains_by_protein=domains_by_protein,
    )
    generic_metadata = _load_protein_metadata(
        path=metadata_path,
        sequences=sequences,
    )
    context_metadata = _merge_context_metadata(
        review_by_protein=review_by_protein,
        generic_metadata=generic_metadata,
    )
    direct_labels, text_evidence = _load_direct_and_text_evidence(
        protein_ids=protein_ids,
        profile=loaded_profile,
        seed_assignments=seed_path,
        seed_catalogue=catalogue_path,
        external_annotations=external_path,
        sequences=sequences,
        review_by_protein=context_metadata,
    )
    group_evidence = _load_group_evidence(
        protein_ids=protein_ids,
        orthofinder_resource=orthofinder_resource,
        orthofinder_results=orthofinder_results,
        group_type=orthofinder_group_type,
        hierarchy_node=orthofinder_hierarchy_node,
        run_id=orthofinder_run_id,
    )
    supplied_redundancy = (
        read_redundancy_clusters(path=redundancy_clusters, protein_ids=protein_ids)
        if redundancy_clusters is not None
        else ()
    )
    partitions = assign_partitions(
        sequences=sequences,
        memberships=group_evidence.memberships,
        validation_fraction=validation_fraction,
        random_seed=random_seed,
        redundancy_memberships=supplied_redundancy,
    )
    unit_by_protein = {item.protein_id: item.partition_key for item in partitions}
    contexts = _build_protein_contexts(
        sequences=sequences,
        domains_by_protein=domains_by_protein,
        domain_assessed=domain_assessed,
        structures_by_protein=structures_by_protein,
        review_by_protein=context_metadata,
        text_evidence=text_evidence,
        direct_labels=direct_labels,
        unit_by_protein=unit_by_protein,
    )
    candidates = _score_initial_candidates(
        contexts=contexts,
        profile=loaded_profile,
        ruleset=ruleset,
    )
    initial_resolutions = _resolve_candidates(
        contexts=contexts,
        candidates=candidates,
        profile=loaded_profile,
        ruleset=ruleset,
    )
    propagated = _orthology_propagation_candidates(
        contexts=contexts,
        resolutions=initial_resolutions,
        memberships=group_evidence.memberships,
        profile=loaded_profile,
        ruleset=ruleset,
    )
    for protein_id, values in propagated.items():
        candidates.setdefault(protein_id, []).extend(values)
    resolutions = _resolve_candidates(
        contexts=contexts,
        candidates=candidates,
        profile=loaded_profile,
        ruleset=ruleset,
    )
    selected_targets = {
        protein_id: resolution.selected
        for protein_id, resolution in resolutions.items()
        if resolution.selected is not None
    }
    if not selected_targets:
        raise InputValidationError(
            "Evidence rules accepted no target proteins. Supply stronger independent "
            "annotations or revise a versioned ruleset; do not lower thresholds silently."
        )
    control_rows, controls_by_label = _select_matched_controls(
        contexts=contexts,
        resolutions=resolutions,
        selected_targets=selected_targets,
        profile=loaded_profile,
        ruleset=ruleset,
    )
    definitions = _label_definition_rows(
        selected_targets=selected_targets,
        profile=loaded_profile,
    )
    excluded_domain_ids = frozenset(
        row["feature_id"] for row in definitions if row["feature_type"] == "DOMAIN"
    )
    label_rows, unresolved_rows = _build_label_rows(
        sequences=sequences,
        profile=loaded_profile,
        template_labels=template_path,
        resolutions=resolutions,
        selected_targets=selected_targets,
        controls_by_label=controls_by_label,
        ruleset=ruleset,
    )
    audit_rows = _candidate_audit_rows(
        contexts=contexts,
        resolutions=resolutions,
    )
    class_rows = _class_summary_rows(
        profile=loaded_profile,
        selected_targets=selected_targets,
        controls_by_label=controls_by_label,
        contexts=contexts,
    )

    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        tables: dict[str, tuple[dict[str, Any], ...]] = {
            "label_assignments": tuple(label_rows),
            "label_evidence_audit": tuple(audit_rows),
            "control_matching_audit": tuple(control_rows),
            "label_definition_features": tuple(definitions),
            "class_labelling_summary": tuple(class_rows),
            "unresolved_assignments": tuple(unresolved_rows),
        }
        table_fields = {
            "label_assignments": LABEL_FIELDS,
            "label_evidence_audit": EVIDENCE_AUDIT_FIELDS,
            "control_matching_audit": CONTROL_MATCH_FIELDS,
            "label_definition_features": LABEL_DEFINITION_FIELDS,
            "class_labelling_summary": CLASS_SUMMARY_FIELDS,
            "unresolved_assignments": UNRESOLVED_FIELDS,
        }
        for table_name, records in tables.items():
            fields = table_fields[table_name]
            write_tsv_atomic(
                path=staging / f"{table_name}.tsv",
                fieldnames=fields,
                records=records,
            )
            frame = pd.DataFrame.from_records(records, columns=fields)
            write_bytes_atomic(
                path=staging / f"{table_name}.xlsx",
                payload=dataframe_to_xlsx_bytes(
                    frame=frame,
                    title=table_name.replace("_", " ").title(),
                ),
            )
        analysis_domain_path: Path | None = None
        if domain_path is not None:
            filtered_domain_rows = _filter_label_defining_domains(
                rows=domain_rows,
                protein_ids=protein_ids,
                excluded_domain_ids=excluded_domain_ids,
                evidence_reference=(
                    f"source_sha256={sha256_file(path=domain_path)};"
                    f"rules_sha256={sha256_file(path=ruleset.source_path)}"
                ),
            )
            analysis_domain_path = staging / "domains.for_signature_analysis.tsv"
            write_tsv_atomic(
                path=analysis_domain_path,
                fieldnames=DOMAIN_FIELDS,
                records=filtered_domain_rows,
            )
            write_bytes_atomic(
                path=staging / "domains.for_signature_analysis.xlsx",
                payload=dataframe_to_xlsx_bytes(
                    frame=pd.DataFrame.from_records(filtered_domain_rows, columns=DOMAIN_FIELDS),
                    title="Domains For Signature Analysis",
                ),
            )
        input_paths = tuple(
            dict.fromkeys(
                (
                    sequences_path,
                    ruleset.source_path,
                    *(
                        path
                        for path in (
                            review_path,
                            metadata_path,
                            domain_path,
                            structure_path,
                            template_path,
                            seed_path,
                            catalogue_path,
                            external_path,
                            (
                                Path(redundancy_clusters).expanduser().resolve()
                                if redundancy_clusters is not None
                                else None
                            ),
                        )
                        if path is not None
                    ),
                    *group_evidence.input_paths,
                )
            )
        )
        input_inventory = _file_inventory(paths=input_paths)
        output_inventory = _directory_inventory(root=staging)
        marker = staging / "EVIDENCE_LABELS.json"
        write_json_atomic(
            path=marker,
            value={
                "schema_version": 1,
                "status": EVIDENCE_BUNDLE_STATUS,
                "action": "EVIDENCE_LED_LABEL_AND_CONTROL_GENERATION",
                "package_version": __version__,
                "interpretation_scope": EVIDENCE_INTERPRETATION_SCOPE,
                "human_review_completed": False,
                "profile_id": loaded_profile.profile_id,
                "profile_version": loaded_profile.profile_version,
                "ruleset_id": ruleset.ruleset_id,
                "ruleset_version": ruleset.ruleset_version,
                "ruleset_sha256": sha256_file(path=ruleset.source_path),
                "protein_count": len(sequences),
                "evidence_supported_target_count": len(selected_targets),
                "evidence_supported_target_unit_count": len(
                    {contexts[protein_id].independence_unit for protein_id in selected_targets}
                ),
                "matched_control_assignment_count": sum(
                    len(protein_ids_for_label)
                    for protein_ids_for_label in controls_by_label.values()
                ),
                "matched_control_protein_count": len(
                    set().union(*controls_by_label.values()) if controls_by_label else set()
                ),
                "ambiguous_protein_count": sum(
                    resolution.status == "AMBIGUOUS" for resolution in resolutions.values()
                ),
                "proposed_not_accepted_count": sum(
                    resolution.status == "PROPOSED" for resolution in resolutions.values()
                ),
                "label_defining_domain_count": len(excluded_domain_ids),
                "orthofinder_source_mode": group_evidence.source_mode,
                "orthofinder_version": group_evidence.version,
                "orthology_propagation_used": bool(propagated),
                "settings": vars(ruleset.settings),
                "inputs": input_inventory,
                "outputs": output_inventory,
                "label_assignments": str(destination / "label_assignments.tsv"),
                "label_assignments_sha256": sha256_file(path=staging / "label_assignments.tsv"),
                "analysis_domains": (
                    str(destination / analysis_domain_path.name)
                    if analysis_domain_path is not None
                    else None
                ),
                "analysis_domains_sha256": (
                    sha256_file(path=analysis_domain_path)
                    if analysis_domain_path is not None
                    else None
                ),
                "label_evidence_audit": str(destination / "label_evidence_audit.tsv"),
                "control_matching_audit": str(destination / "control_matching_audit.tsv"),
                "label_definition_features": str(destination / "label_definition_features.tsv"),
                "class_labelling_summary": str(destination / "class_labelling_summary.tsv"),
                "unresolved_assignments": str(destination / "unresolved_assignments.tsv"),
                "warning": (
                    "Assignments are evidence-supported automated proposals. They are suitable "
                    "for provisional hypothesis generation only until the evidence summary and "
                    "conflicts receive human scientific review."
                ),
            },
        )
        os.replace(staging, destination)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    published_marker = destination / "EVIDENCE_LABELS.json"
    verify_evidence_label_bundle(bundle_dir=destination)
    LOGGER.info(
        "Published provisional evidence labels targets=%d controls=%d at %s",
        len(selected_targets),
        sum(len(values) for values in controls_by_label.values()),
        destination,
    )
    return published_marker


def verify_evidence_label_bundle(*, bundle_dir: Path) -> Mapping[str, Any]:
    """Verify a published evidence-label bundle and every declared output.

    Args:
        bundle_dir: Published bundle directory.

    Returns:
        Validated marker mapping.

    Raises:
        InputValidationError: If the marker, paths, sizes or checksums disagree.
    """

    root = Path(bundle_dir).expanduser().resolve()
    marker_path = root / "EVIDENCE_LABELS.json"
    document = read_json(path=marker_path)
    if not isinstance(document, Mapping) or document.get("status") != EVIDENCE_BUNDLE_STATUS:
        raise InputValidationError(f"Invalid evidence-label completion marker: {marker_path}")
    if document.get("schema_version") != 1:
        raise InputValidationError("Evidence-label marker schema_version must be 1.")
    if document.get("interpretation_scope") != EVIDENCE_INTERPRETATION_SCOPE:
        raise InputValidationError("Evidence-label marker has an unsafe interpretation scope.")
    if document.get("human_review_completed") is not False:
        raise InputValidationError(
            "Automated evidence-label marker must not claim completed human review."
        )
    for field in (
        "package_version",
        "profile_id",
        "profile_version",
        "ruleset_id",
        "ruleset_version",
        "warning",
    ):
        if not isinstance(document.get(field), str) or not str(document[field]).strip():
            raise InputValidationError(f"Evidence-label marker lacks text field {field!r}.")
    rules_digest = str(document.get("ruleset_sha256") or "")
    if re.fullmatch(r"[0-9a-f]{64}", rules_digest) is None:
        raise InputValidationError("Evidence-label marker has an invalid ruleset checksum.")
    for field in (
        "protein_count",
        "evidence_supported_target_count",
        "evidence_supported_target_unit_count",
        "matched_control_assignment_count",
        "matched_control_protein_count",
        "ambiguous_protein_count",
        "proposed_not_accepted_count",
        "label_defining_domain_count",
    ):
        value = document.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise InputValidationError(
                f"Evidence-label marker field {field!r} must be a non-negative integer."
            )
    inputs = document.get("inputs")
    if not isinstance(inputs, list) or len(inputs) < 2:
        raise InputValidationError("Evidence-label marker has no complete input inventory.")
    declared_inputs: set[Path] = set()
    input_digests: set[str] = set()
    for index, row in enumerate(inputs):
        if not isinstance(row, Mapping):
            raise InputValidationError(f"Malformed evidence input record {index}.")
        candidate = Path(str(row.get("path") or "")).expanduser()
        if not candidate.is_absolute():
            raise InputValidationError(f"Evidence input path is not absolute: {candidate}")
        candidate = candidate.resolve()
        if candidate in declared_inputs:
            raise InputValidationError(f"Evidence input is declared twice: {candidate}")
        declared_inputs.add(candidate)
        if not candidate.is_file():
            raise InputValidationError(f"Evidence input is missing: {candidate}")
        if candidate.stat().st_size != row.get("size_bytes"):
            raise InputValidationError(f"Evidence input size differs: {candidate}")
        digest = sha256_file(path=candidate)
        if digest != row.get("sha256"):
            raise InputValidationError(f"Evidence input checksum differs: {candidate}")
        input_digests.add(digest)
    if rules_digest not in input_digests:
        raise InputValidationError(
            "Evidence-label ruleset checksum is absent from the input inventory."
        )
    outputs = document.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise InputValidationError("Evidence-label marker has no output inventory.")
    declared: set[str] = set()
    for index, row in enumerate(outputs):
        if not isinstance(row, Mapping):
            raise InputValidationError(f"Malformed evidence output record {index}.")
        relative_text = str(row.get("relative_path") or "")
        relative = Path(relative_text)
        if (
            not relative_text
            or relative.is_absolute()
            or ".." in relative.parts
            or relative_text in declared
        ):
            raise InputValidationError(
                f"Unsafe or duplicate evidence output path: {relative_text!r}"
            )
        declared.add(relative_text)
        candidate = (root / relative).resolve()
        if root not in candidate.parents or not candidate.is_file():
            raise InputValidationError(f"Evidence output is missing or unsafe: {candidate}")
        if candidate.stat().st_size != row.get("size_bytes"):
            raise InputValidationError(f"Evidence output size differs: {candidate}")
        if sha256_file(path=candidate) != row.get("sha256"):
            raise InputValidationError(f"Evidence output checksum differs: {candidate}")
    actual = {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path != marker_path
    }
    if actual != declared:
        raise InputValidationError(
            "Evidence-label file inventory differs from its marker: "
            f"undeclared={sorted(actual - declared)}, missing={sorted(declared - actual)}"
        )
    required_outputs = {
        f"{name}.{extension}"
        for name in (
            "label_assignments",
            "label_evidence_audit",
            "control_matching_audit",
            "label_definition_features",
            "class_labelling_summary",
            "unresolved_assignments",
        )
        for extension in ("tsv", "xlsx")
    }
    if not required_outputs <= declared:
        raise InputValidationError(
            "Evidence-label bundle lacks required tabular outputs: "
            f"{sorted(required_outputs - declared)}"
        )
    labels = root / "label_assignments.tsv"
    if document.get("label_assignments") != str(labels):
        raise InputValidationError("Evidence marker label-assignment path differs.")
    if document.get("label_assignments_sha256") != sha256_file(path=labels):
        raise InputValidationError("Evidence marker label-assignment checksum differs.")
    for field in (
        "label_evidence_audit",
        "control_matching_audit",
        "label_definition_features",
        "class_labelling_summary",
        "unresolved_assignments",
    ):
        expected = root / f"{field}.tsv"
        if document.get(field) != str(expected):
            raise InputValidationError(
                f"Evidence marker {field!r} path differs from the bundle authority."
            )
    analysis_domains = document.get("analysis_domains")
    if analysis_domains is not None:
        domains_path = root / "domains.for_signature_analysis.tsv"
        if analysis_domains != str(domains_path):
            raise InputValidationError("Evidence marker analysis-domain path differs.")
        if document.get("analysis_domains_sha256") != sha256_file(path=domains_path):
            raise InputValidationError("Evidence marker analysis-domain checksum differs.")
    elif document.get("analysis_domains_sha256") is not None:
        raise InputValidationError(
            "Evidence marker has an analysis-domain checksum without a domain authority."
        )
    return document


def read_evidence_audit_table(
    *, path: Path | None, fields: Sequence[str]
) -> tuple[dict[str, Any], ...]:
    """Read one optional generated evidence audit table strictly.

    Args:
        path: Optional TSV path.
        fields: Exact required fields for the selected audit relation.

    Returns:
        Ordered schema-ready rows, or an empty tuple when no path is supplied.

    Raises:
        InputValidationError: If fields are invalid or a row is malformed.
    """

    if path is None:
        return ()
    required = tuple(fields)
    if not required or len(required) != len(set(required)):
        raise InputValidationError("Evidence audit fields must be non-empty and unique.")
    rows: list[dict[str, Any]] = []
    for row_number, row in enumerate(
        iter_tsv(path=path, required_fields=required, allow_empty=True),
        start=2,
    ):
        converted: dict[str, Any] = {}
        for field in required:
            value = row[field]
            if field in _AUDIT_INTEGER_FIELDS:
                converted[field] = _audit_integer(
                    value=value,
                    field_name=field,
                    row_number=row_number,
                )
            elif field in _AUDIT_FLOAT_FIELDS:
                converted[field] = _audit_float(
                    value=value,
                    field_name=field,
                    row_number=row_number,
                )
            elif field in _AUDIT_BOOLEAN_FIELDS:
                converted[field] = _audit_boolean(
                    value=value,
                    field_name=field,
                    row_number=row_number,
                )
            else:
                converted[field] = value
        rows.append(converted)
    return tuple(rows)


def _parse_evidence_settings(*, value: Any) -> EvidenceSettings:
    """Validate the ruleset settings mapping."""

    if not isinstance(value, Mapping):
        raise InputValidationError("evidence_rules.settings must be a mapping.")
    fields = {
        "annotation_score",
        "specific_annotation_score",
        "pfam_score",
        "supporting_pfam_score",
        "orthology_score",
        "trusted_label_score",
        "minimum_acceptance_score",
        "minimum_evidence_groups",
        "minimum_score_margin",
        "require_target_structure_eligible",
        "orthology_propagation_enabled",
        "minimum_orthology_anchor_proteins",
        "control_units_per_target_unit",
        "minimum_control_units_per_background",
        "require_species_match",
        "require_structure_match",
        "maximum_log2_length_difference",
        "maximum_domain_count_difference",
        "exclude_input_candidates_from_controls",
    }
    _reject_unknown_fields(value=value, allowed=fields, field_name="evidence_rules.settings")
    numbers = {
        name: _finite_nonnegative_number(value=value.get(name), field_name=name)
        for name in (
            "annotation_score",
            "specific_annotation_score",
            "pfam_score",
            "supporting_pfam_score",
            "orthology_score",
            "trusted_label_score",
            "minimum_acceptance_score",
            "minimum_score_margin",
            "maximum_log2_length_difference",
        )
    }
    integers = {
        name: _positive_integer(value=value.get(name), field_name=name)
        for name in (
            "minimum_evidence_groups",
            "minimum_orthology_anchor_proteins",
            "control_units_per_target_unit",
            "minimum_control_units_per_background",
        )
    }
    maximum_domain_difference = _nonnegative_integer(
        value=value.get("maximum_domain_count_difference"),
        field_name="maximum_domain_count_difference",
    )
    booleans = {
        name: _boolean(value=value.get(name), field_name=name)
        for name in (
            "require_target_structure_eligible",
            "orthology_propagation_enabled",
            "require_species_match",
            "require_structure_match",
            "exclude_input_candidates_from_controls",
        )
    }
    return EvidenceSettings(
        **numbers,
        **integers,
        maximum_domain_count_difference=maximum_domain_difference,
        **booleans,
    )


def _parse_evidence_rule(
    *, value: Any, index: int, labels_by_id: Mapping[str, Any]
) -> EvidenceRule:
    """Validate and compile one label rule."""

    if not isinstance(value, Mapping):
        raise InputValidationError(f"Evidence label rule {index} must be a mapping.")
    fields = {
        "rule_id",
        "label_id",
        "priority",
        "annotation_patterns",
        "specific_annotation_patterns",
        "exclusion_patterns",
        "required_domain_patterns",
        "supporting_domain_patterns",
        "allow_specific_annotation_only",
        "propagate_by_orthology",
    }
    _reject_unknown_fields(
        value=value,
        allowed=fields,
        field_name=f"evidence_rules.label_rules[{index}]",
    )
    label_id = validate_identifier(
        value=value.get("label_id"), field_name=f"label_rules[{index}].label_id"
    )
    label = labels_by_id.get(label_id)
    if label is None:
        raise InputValidationError(f"Evidence rule references unknown profile label: {label_id!r}")
    if not label.reviewed_positive_allowed or label.label_id.startswith("control"):
        raise InputValidationError(
            f"Evidence rule target must permit positive non-control assignments: {label_id!r}"
        )
    annotation_patterns = _compile_patterns(
        value=value.get("annotation_patterns", []),
        field_name=f"label_rules[{index}].annotation_patterns",
    )
    specific_patterns = _compile_patterns(
        value=value.get("specific_annotation_patterns", []),
        field_name=f"label_rules[{index}].specific_annotation_patterns",
    )
    required_domains = _compile_patterns(
        value=value.get("required_domain_patterns", []),
        field_name=f"label_rules[{index}].required_domain_patterns",
    )
    if not annotation_patterns and not specific_patterns:
        raise InputValidationError(
            f"Evidence rule {label_id!r} needs an annotation or specific-annotation pattern."
        )
    return EvidenceRule(
        rule_id=validate_identifier(
            value=value.get("rule_id"), field_name=f"label_rules[{index}].rule_id"
        ),
        label_id=label_id,
        component_role=label.component_role or "UNKNOWN",
        priority=_nonnegative_integer(
            value=value.get("priority", 0), field_name=f"label_rules[{index}].priority"
        ),
        annotation_patterns=annotation_patterns,
        specific_annotation_patterns=specific_patterns,
        exclusion_patterns=_compile_patterns(
            value=value.get("exclusion_patterns", []),
            field_name=f"label_rules[{index}].exclusion_patterns",
        ),
        required_domain_patterns=required_domains,
        supporting_domain_patterns=_compile_patterns(
            value=value.get("supporting_domain_patterns", []),
            field_name=f"label_rules[{index}].supporting_domain_patterns",
        ),
        allow_specific_annotation_only=_boolean(
            value=value.get("allow_specific_annotation_only", False),
            field_name=f"label_rules[{index}].allow_specific_annotation_only",
        ),
        propagate_by_orthology=_boolean(
            value=value.get("propagate_by_orthology", bool(required_domains)),
            field_name=f"label_rules[{index}].propagate_by_orthology",
        ),
    )


def _load_domain_context(
    *, path: Path | None, sequences: tuple[SequenceRecord, ...]
) -> tuple[
    Mapping[str, Mapping[str, str]],
    frozenset[str],
    tuple[dict[str, str], ...],
]:
    """Load domain hits and assessment coverage for evidence scoring."""

    if path is None:
        return {}, frozenset(), ()
    hits, assessments = read_domains(path=path, sequences=sequences)
    domains: dict[str, dict[str, str]] = defaultdict(dict)
    for hit in hits:
        previous = domains[hit.protein_id].setdefault(hit.domain_id, hit.domain_name)
        if previous != hit.domain_name:
            raise InputValidationError(
                f"Conflicting domain names for {(hit.protein_id, hit.domain_id)!r}."
            )
    assessed = frozenset(
        item.protein_id
        for item in assessments
        if item.assessment_status.value in {"ASSESSED_WITH_HIT", "ASSESSED_NO_HIT"}
    )
    rows = tuple(iter_tsv(path=path, required_fields=DOMAIN_FIELDS))
    return domains, assessed, rows


def _load_structure_context(
    *, path: Path | None, sequences: tuple[SequenceRecord, ...]
) -> Mapping[str, tuple[bool, float | None]]:
    """Load structure eligibility and confidence without using folds as labels."""

    if path is None:
        return {}
    records = read_structures(path=path, sequences=sequences)
    result: dict[str, tuple[bool, float | None]] = {}
    for record in records:
        eligible = record.analysis_eligibility_status in _ELIGIBLE_STRUCTURE_STATES
        candidate = (eligible, record.mean_confidence)
        previous = result.get(record.protein_id)
        if previous is None or (
            candidate[0],
            candidate[1] if candidate[1] is not None else -1.0,
        ) > (previous[0], previous[1] if previous[1] is not None else -1.0):
            result[record.protein_id] = candidate
    return result


def _load_review_context(
    *,
    path: Path | None,
    sequences: tuple[SequenceRecord, ...],
    domains_by_protein: Mapping[str, Mapping[str, str]],
) -> Mapping[str, dict[str, Any]]:
    """Load optional predecessor hints and independently verify stable covariates."""

    if path is None:
        return {}
    by_id = {record.protein_id: record for record in sequences}
    result: dict[str, dict[str, Any]] = {}
    for row in iter_tsv(path=path, required_fields=REVIEW_CONTEXT_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        if protein_id not in by_id:
            raise InputValidationError(
                f"Review context references protein absent from FASTA: {protein_id!r}"
            )
        if protein_id in result:
            raise InputValidationError(f"Review context repeats protein: {protein_id!r}")
        try:
            declared_length = int(row["sequence_length"])
        except ValueError as error:
            raise InputValidationError(
                f"Review context has invalid sequence_length for {protein_id!r}."
            ) from error
        if declared_length != by_id[protein_id].sequence_length:
            raise InputValidationError(
                f"Review-context sequence length differs for {protein_id!r}."
            )
        review_domains = frozenset(_split_pipe(value=row["pfam_accessions"]))
        loaded_domains = frozenset(domains_by_protein.get(protein_id, {}))
        if loaded_domains and review_domains != loaded_domains:
            raise InputValidationError(
                f"Review-context Pfam accessions differ from domains.tsv for {protein_id!r}."
            )
        candidate_states = frozenset(
            item.upper() for item in _split_pipe(value=row["input_candidate_states"])
        )
        result[protein_id] = {
            "species": _split_pipe(value=row["species"]),
            "input_candidate": bool(candidate_states & _TRUE_VALUES),
            "upstream_e3_families": row["upstream_e3_families"],
            "upstream_evidence_roles": row["upstream_evidence_roles"],
            "annotation_status_details": row["annotation_status_details"],
            "declared_structure_eligible": row["structure_analysis_eligible"].strip().upper()
            in _TRUE_VALUES,
        }
    missing = frozenset(by_id) - frozenset(result)
    if missing:
        raise InputValidationError(
            f"Review context lacks {len(missing)} FASTA proteins; first={sorted(missing)[:5]}"
        )
    return result


def _load_protein_metadata(
    *, path: Path | None, sequences: tuple[SequenceRecord, ...]
) -> Mapping[str, dict[str, Any]]:
    """Load generic species and candidate-state covariates.

    Args:
        path: Optional generic metadata TSV.
        sequences: Authoritative sequence records.

    Returns:
        Metadata keyed by protein identifier.

    Raises:
        InputValidationError: If coverage, identifiers or Boolean values are invalid.
    """

    if path is None:
        return {}
    protein_ids = frozenset(record.protein_id for record in sequences)
    result: dict[str, dict[str, Any]] = {}
    for row in iter_tsv(path=path, required_fields=PROTEIN_METADATA_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        if protein_id not in protein_ids:
            raise InputValidationError(
                f"Protein metadata references protein absent from FASTA: {protein_id!r}"
            )
        if protein_id in result:
            raise InputValidationError(f"Protein metadata repeats protein: {protein_id!r}")
        candidate_state = row["input_candidate"].strip().upper()
        if candidate_state not in {"TRUE", "FALSE"}:
            raise InputValidationError(
                "protein_metadata.input_candidate must be TRUE or FALSE; "
                f"received {row['input_candidate']!r} for {protein_id!r}."
            )
        result[protein_id] = {
            "species": _split_pipe(value=row["species"]),
            "input_candidate": candidate_state == "TRUE",
        }
    missing = protein_ids - frozenset(result)
    if missing:
        raise InputValidationError(
            f"Protein metadata lacks {len(missing)} FASTA proteins; first={sorted(missing)[:5]}"
        )
    return result


def _merge_context_metadata(
    *,
    review_by_protein: Mapping[str, Mapping[str, Any]],
    generic_metadata: Mapping[str, Mapping[str, Any]],
) -> Mapping[str, dict[str, Any]]:
    """Merge predecessor review context with generic matching metadata.

    Args:
        review_by_protein: Optional predecessor-derived review context.
        generic_metadata: Optional generic metadata authority.

    Returns:
        Merged metadata without silently resolving disagreements.

    Raises:
        InputValidationError: If both authorities disagree on a supplied value.
    """

    result = {protein_id: dict(row) for protein_id, row in review_by_protein.items()}
    for protein_id, metadata in generic_metadata.items():
        existing = result.setdefault(protein_id, {})
        existing_species = tuple(existing.get("species", ()))
        supplied_species = tuple(metadata.get("species", ()))
        if existing_species and supplied_species and existing_species != supplied_species:
            raise InputValidationError(
                f"Generic and predecessor metadata disagree on species for {protein_id!r}."
            )
        if "input_candidate" in existing and (
            bool(existing["input_candidate"]) != bool(metadata.get("input_candidate"))
        ):
            raise InputValidationError(
                f"Generic and predecessor metadata disagree on input_candidate for {protein_id!r}."
            )
        existing["species"] = existing_species or supplied_species
        existing["input_candidate"] = bool(metadata.get("input_candidate"))
    return result


def _load_direct_and_text_evidence(
    *,
    protein_ids: frozenset[str],
    profile: ProteinProfile,
    seed_assignments: Path | None,
    seed_catalogue: Path | None,
    external_annotations: Path | None,
    sequences: tuple[SequenceRecord, ...],
    review_by_protein: Mapping[str, Mapping[str, Any]],
) -> tuple[Mapping[str, tuple[LabelAssignment, ...]], Mapping[str, tuple[TextEvidence, ...]]]:
    """Combine trusted direct labels with independent annotation text evidence."""

    direct: dict[str, list[LabelAssignment]] = defaultdict(list)
    texts: dict[str, list[TextEvidence]] = defaultdict(list)
    for protein_id, row in review_by_protein.items():
        values = " | ".join(
            value
            for value in (
                str(row.get("upstream_e3_families") or ""),
                str(row.get("upstream_evidence_roles") or ""),
                str(row.get("annotation_status_details") or ""),
            )
            if value.strip()
        )
        if values:
            texts[protein_id].append(
                TextEvidence(
                    text=values,
                    evidence_group="UPSTREAM_CONTEXT",
                    evidence_source="Completed predecessor annotation context",
                    evidence_reference="review_context",
                    reliability=1.0,
                    trusted_status=False,
                )
            )
    if seed_assignments is not None:
        assignments = read_label_assignments(
            path=seed_assignments,
            protein_ids=protein_ids,
            label_ids=profile.label_ids(),
        )
        for assignment in assignments:
            if assignment.curation_status == CurationStatus.REVIEWED_POSITIVE:
                direct[assignment.protein_id].append(assignment)
    if seed_catalogue is not None:
        _add_seed_catalogue_evidence(
            path=seed_catalogue,
            sequences=sequences,
            texts=texts,
        )
    if external_annotations is not None:
        for row in iter_tsv(path=external_annotations, required_fields=EXTERNAL_ANNOTATION_FIELDS):
            protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
            if protein_id not in protein_ids:
                raise InputValidationError(
                    f"External annotation references unknown protein: {protein_id!r}"
                )
            source = validate_text(value=row["evidence_source"], field_name="evidence_source")
            reference = validate_text(
                value=row["evidence_reference"],
                field_name="evidence_reference",
                allow_empty=True,
            )
            status = validate_identifier(
                value=row["evidence_status"], field_name="evidence_status"
            ).upper()
            scope = validate_identifier(
                value=row["annotation_scope"], field_name="annotation_scope"
            ).upper()
            label_id = row["label_id"].strip()
            if label_id:
                label_id = validate_identifier(value=label_id, field_name="label_id")
                if label_id not in profile.label_ids():
                    raise InputValidationError(
                        f"External annotation references unknown profile label: {label_id!r}"
                    )
                if status in _TRUSTED_EXTERNAL_STATUSES and scope in {"DIRECT", "EXACT_PROTEIN"}:
                    label = next(item for item in profile.labels if item.label_id == label_id)
                    direct[protein_id].append(
                        LabelAssignment(
                            protein_id=protein_id,
                            label_id=label_id,
                            curation_status=CurationStatus.REVIEWED_POSITIVE,
                            evidence_status=status,
                            evidence_source=source,
                            evidence_reference=reference,
                            component_role=label.component_role or "UNKNOWN",
                            curation_reason=(
                                "Trusted direct label supplied by an external authority."
                            ),
                        )
                    )
            annotation = validate_text(
                value=row["annotation_text"],
                field_name="annotation_text",
                allow_empty=True,
                maximum_length=32_767,
            )
            if annotation:
                group_token = hashlib.sha256(source.encode("utf-8")).hexdigest()[:16]
                texts[protein_id].append(
                    TextEvidence(
                        text=annotation,
                        evidence_group=f"EXTERNAL_{group_token}",
                        evidence_source=source,
                        evidence_reference=reference,
                        reliability=1.0,
                        trusted_status=(
                            status in _TRUSTED_EXTERNAL_STATUSES
                            and scope in {"DIRECT", "EXACT_PROTEIN"}
                        ),
                    )
                )
    return (
        {key: tuple(values) for key, values in direct.items()},
        {key: tuple(values) for key, values in texts.items()},
    )


def _add_seed_catalogue_evidence(
    *,
    path: Path,
    sequences: tuple[SequenceRecord, ...],
    texts: dict[str, list[TextEvidence]],
) -> None:
    """Add direct or cluster-context annotations from the supplied E3 catalogue."""

    sequences_by_id = {record.protein_id: record.sequence for record in sequences}
    catalogue_digest = sha256_file(path=path)
    observed: set[str] = set()
    for row in iter_tsv(path=path, required_fields=SEED_CATALOGUE_FIELDS):
        protein_id = row["seed_id"].strip()
        if not protein_id or protein_id not in sequences_by_id:
            continue
        protein_id = validate_identifier(value=protein_id, field_name="seed_id")
        if protein_id in observed:
            raise InputValidationError(f"Seed catalogue repeats protein identifier: {protein_id!r}")
        observed.add(protein_id)
        catalogue_sequence = "".join(row["protein_sequence"].split()).upper().rstrip("*")
        if catalogue_sequence and catalogue_sequence != sequences_by_id[protein_id]:
            raise InputValidationError(
                f"Seed catalogue sequence differs from FASTA for {protein_id!r}."
            )
        direct_category = row["seed_category"].strip()
        associated = row["associated_seed_categories"].strip()
        scope = row["annotation_scope"].strip().casefold()
        source = row["catalogue_source"].strip() or path.name
        reviewed = row["seed_review_status"].strip().casefold() == "reviewed"
        if direct_category:
            texts[protein_id].append(
                TextEvidence(
                    text=direct_category,
                    evidence_group="CURATED_SEED_CATALOGUE",
                    evidence_source=source,
                    evidence_reference=catalogue_digest,
                    reliability=1.0,
                    trusted_status=reviewed and scope in {"direct", "exact_protein"},
                )
            )
        if associated:
            texts[protein_id].append(
                TextEvidence(
                    text=associated,
                    evidence_group="UPSTREAM_CONTEXT",
                    evidence_source=source,
                    evidence_reference=catalogue_digest,
                    reliability=0.5,
                    trusted_status=False,
                )
            )


def _load_group_evidence(
    *,
    protein_ids: frozenset[str],
    orthofinder_resource: Path | None,
    orthofinder_results: Path | None,
    group_type: str,
    hierarchy_node: str,
    run_id: str,
) -> GroupEvidence:
    """Load optional OrthoFinder memberships and exact authority paths."""

    if orthofinder_resource is not None:
        resource = discover_orthofinder_resource(resource_dir=orthofinder_resource)
        memberships = read_resource_memberships(
            resource=resource,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
            protein_ids=protein_ids,
        )
        return GroupEvidence(
            memberships=memberships,
            input_paths=resource_input_paths(resource=resource),
            source_mode="ORTHOFINDER_RESULTS_RESOURCE",
            version=resource.orthofinder_version,
        )
    if orthofinder_results is not None:
        layout = discover_orthofinder_layout(results_dir=orthofinder_results)
        memberships = read_group_memberships(
            layout=layout,
            run_id=validate_identifier(value=run_id, field_name="orthofinder_run_id"),
            group_type=group_type,
            hierarchy_node=(hierarchy_node if group_type.upper() == "HOG" else ""),
            protein_ids=protein_ids,
        )
        group_path = (
            next(path for path in layout.hog_paths if path.stem == hierarchy_node)
            if group_type.upper() == "HOG"
            else layout.orthogroups_path
        )
        paths = tuple(
            dict.fromkeys(
                path
                for path in (
                    group_path,
                    layout.sequence_ids_path,
                    layout.log_path,
                    *layout.completion_authority_paths,
                )
                if path is not None
            )
        )
        return GroupEvidence(
            memberships=memberships,
            input_paths=paths,
            source_mode=layout.source_mode,
            version=layout.version,
        )
    return GroupEvidence(memberships=(), input_paths=(), source_mode="NOT_SUPPLIED", version="")


def _build_protein_contexts(
    *,
    sequences: tuple[SequenceRecord, ...],
    domains_by_protein: Mapping[str, Mapping[str, str]],
    domain_assessed: frozenset[str],
    structures_by_protein: Mapping[str, tuple[bool, float | None]],
    review_by_protein: Mapping[str, Mapping[str, Any]],
    text_evidence: Mapping[str, tuple[TextEvidence, ...]],
    direct_labels: Mapping[str, tuple[LabelAssignment, ...]],
    unit_by_protein: Mapping[str, str],
) -> Mapping[str, ProteinEvidenceContext]:
    """Assemble immutable per-protein evidence and matching covariates."""

    result: dict[str, ProteinEvidenceContext] = {}
    for sequence in sequences:
        review = review_by_protein.get(sequence.protein_id, {})
        structure_eligible, mean_confidence = structures_by_protein.get(
            sequence.protein_id, (False, None)
        )
        declared_eligible = review.get("declared_structure_eligible")
        if declared_eligible is not None and bool(declared_eligible) != structure_eligible:
            raise InputValidationError(
                f"Review and structure authorities disagree on eligibility for "
                f"{sequence.protein_id!r}."
            )
        domain_map = dict(domains_by_protein.get(sequence.protein_id, {}))
        result[sequence.protein_id] = ProteinEvidenceContext(
            protein_id=sequence.protein_id,
            sequence_length=sequence.sequence_length,
            sequence_sha256=sequence.sequence_sha256,
            species=tuple(review.get("species", ())),
            domain_ids=tuple(sorted(domain_map)),
            domain_names=domain_map,
            domain_count=len(domain_map),
            pfam_assessed=sequence.protein_id in domain_assessed,
            structure_eligible=structure_eligible,
            mean_confidence=mean_confidence,
            input_candidate=bool(review.get("input_candidate", False)),
            text_evidence=tuple(text_evidence.get(sequence.protein_id, ())),
            direct_labels=tuple(direct_labels.get(sequence.protein_id, ())),
            independence_unit=unit_by_protein[sequence.protein_id],
        )
    return result


def _score_initial_candidates(
    *,
    contexts: Mapping[str, ProteinEvidenceContext],
    profile: ProteinProfile,
    ruleset: EvidenceRuleSet,
) -> dict[str, list[LabelCandidate]]:
    """Score trusted labels, annotation text and domain architecture."""

    labels_by_id = {label.label_id: label for label in profile.labels}
    rules_by_label = {rule.label_id: rule for rule in ruleset.rules}
    result: dict[str, list[LabelCandidate]] = defaultdict(list)
    for protein_id, context in contexts.items():
        for assignment in context.direct_labels:
            label = labels_by_id[assignment.label_id]
            if not label.reviewed_positive_allowed or label.label_id.startswith("control"):
                continue
            rule = rules_by_label.get(label.label_id)
            item = EvidenceItem(
                evidence_group="TRUSTED_DIRECT_LABEL",
                evidence_kind="CURATED_LABEL",
                value=label.label_id,
                score=ruleset.settings.trusted_label_score,
                evidence_source=assignment.evidence_source,
                evidence_reference=assignment.evidence_reference,
            )
            structure_ok = _structure_acceptance(
                context=context,
                settings=ruleset.settings,
            )
            result[protein_id].append(
                LabelCandidate(
                    protein_id=protein_id,
                    label_id=label.label_id,
                    rule_id=rule.rule_id if rule is not None else "trusted_direct_label",
                    component_role=label.component_role or assignment.component_role,
                    score=ruleset.settings.trusted_label_score,
                    evidence_groups=("TRUSTED_DIRECT_LABEL",),
                    evidence_items=(item,),
                    accepted=structure_ok,
                    confidence_tier="DIRECT_CURATED",
                    reason=(
                        "Trusted direct curated label."
                        if structure_ok
                        else "Trusted label retained for review but target structure is ineligible."
                    ),
                    priority=rule.priority if rule is not None else 10_000,
                    trusted_direct=True,
                )
            )
        for rule in ruleset.rules:
            candidate = _score_rule(context=context, rule=rule, settings=ruleset.settings)
            if candidate is not None:
                result[protein_id].append(candidate)
    return result


def _score_rule(
    *, context: ProteinEvidenceContext, rule: EvidenceRule, settings: EvidenceSettings
) -> LabelCandidate | None:
    """Score one rule against one protein without forcing a classification."""

    items: list[EvidenceItem] = []
    best_by_group: dict[str, EvidenceItem] = {}
    specific_trusted = False
    for evidence in context.text_evidence:
        if any(pattern.search(evidence.text) for pattern in rule.exclusion_patterns):
            continue
        specific = any(
            pattern.search(evidence.text) for pattern in rule.specific_annotation_patterns
        )
        ordinary = any(pattern.search(evidence.text) for pattern in rule.annotation_patterns)
        if not specific and not ordinary:
            continue
        base_score = settings.specific_annotation_score if specific else settings.annotation_score
        item = EvidenceItem(
            evidence_group=evidence.evidence_group,
            evidence_kind="SPECIFIC_ANNOTATION" if specific else "ANNOTATION",
            value=evidence.text,
            score=base_score * evidence.reliability,
            evidence_source=evidence.evidence_source,
            evidence_reference=evidence.evidence_reference,
        )
        previous = best_by_group.get(evidence.evidence_group)
        if previous is None or item.score > previous.score:
            best_by_group[evidence.evidence_group] = item
        specific_trusted = specific_trusted or (specific and evidence.trusted_status)
    items.extend(best_by_group.values())
    matched_required = _matched_domain_features(
        context=context,
        patterns=rule.required_domain_patterns,
    )
    required_domains_complete = len(matched_required) == len(rule.required_domain_patterns)
    if matched_required:
        items.extend(
            EvidenceItem(
                evidence_group="PFAM",
                evidence_kind="REQUIRED_DOMAIN",
                value=domain_id,
                score=(settings.pfam_score if index == 0 else 0.0),
                evidence_source="Pfam domain assessment",
                evidence_reference="domains.tsv",
                feature_type="DOMAIN",
                feature_id=domain_id,
                feature_name=context.domain_names[domain_id],
            )
            for index, domain_id in enumerate(dict.fromkeys(matched_required))
        )
    matched_supporting = _matched_domain_features(
        context=context,
        patterns=rule.supporting_domain_patterns,
    )
    if matched_supporting:
        items.extend(
            EvidenceItem(
                evidence_group="PFAM",
                evidence_kind="SUPPORTING_DOMAIN",
                value=domain_id,
                score=(settings.supporting_pfam_score if index == 0 else 0.0),
                evidence_source="Pfam domain assessment",
                evidence_reference="domains.tsv",
                feature_type="DOMAIN",
                feature_id=domain_id,
                feature_name=context.domain_names[domain_id],
            )
            for index, domain_id in enumerate(dict.fromkeys(matched_supporting))
        )
    if not items:
        return None
    annotation_items = tuple(
        item for item in items if item.evidence_kind in {"ANNOTATION", "SPECIFIC_ANNOTATION"}
    )
    score = sum(item.score for item in items)
    groups = tuple(sorted({item.evidence_group for item in items if item.score > 0.0}))
    required_ok = not rule.required_domain_patterns or required_domains_complete
    single_specific_ok = (
        rule.allow_specific_annotation_only
        and specific_trusted
        and score >= settings.minimum_acceptance_score
    )
    evidence_ok = (
        bool(annotation_items)
        and required_ok
        and score >= settings.minimum_acceptance_score
        and len(groups) >= settings.minimum_evidence_groups
    ) or single_specific_ok
    structure_ok = _structure_acceptance(context=context, settings=settings)
    accepted = evidence_ok and structure_ok
    if not structure_ok:
        reason = (
            "Target-like evidence found, but the required structure is unavailable or ineligible."
        )
    elif not annotation_items:
        reason = "Domain-only evidence cannot establish protein type."
    elif not required_ok:
        reason = "Annotation matched, but the required domain architecture was not observed."
    elif score < settings.minimum_acceptance_score:
        reason = "Evidence score is below the configured acceptance threshold."
    elif len(groups) < settings.minimum_evidence_groups and not single_specific_ok:
        reason = "Too few independent evidence groups support the assignment."
    else:
        reason = "Evidence and structure policies satisfied."
    confidence = (
        "HIGH"
        if accepted and len(groups) >= settings.minimum_evidence_groups
        else "DIRECT_CURATED"
        if accepted and single_specific_ok
        else "PROPOSED"
    )
    return LabelCandidate(
        protein_id=context.protein_id,
        label_id=rule.label_id,
        rule_id=rule.rule_id,
        component_role=rule.component_role,
        score=score,
        evidence_groups=groups,
        evidence_items=tuple(sorted(items, key=_evidence_item_sort_key)),
        accepted=accepted,
        confidence_tier=confidence,
        reason=reason,
        priority=rule.priority,
    )


def _structure_acceptance(*, context: ProteinEvidenceContext, settings: EvidenceSettings) -> bool:
    """Return whether target structure policy permits analysis membership."""

    return context.structure_eligible if settings.require_target_structure_eligible else True


def _resolve_candidates(
    *,
    contexts: Mapping[str, ProteinEvidenceContext],
    candidates: Mapping[str, Sequence[LabelCandidate]],
    profile: ProteinProfile,
    ruleset: EvidenceRuleSet,
) -> Mapping[str, CandidateResolution]:
    """Resolve ancestry, priority and conflicts without inventing certainty."""

    result: dict[str, CandidateResolution] = {}
    for protein_id in contexts:
        values = tuple(
            sorted(
                candidates.get(protein_id, ()),
                key=lambda item: (
                    not item.accepted,
                    not item.trusted_direct,
                    -item.score,
                    -len(label_ancestors(profile=profile, label_id=item.label_id)),
                    -item.priority,
                    item.label_id,
                    item.rule_id,
                ),
            )
        )
        accepted = tuple(item for item in values if item.accepted)
        decisions: dict[tuple[str, str], str] = {}
        selected: LabelCandidate | None = None
        conflicts: tuple[str, ...] = ()
        if accepted:
            trusted = tuple(item for item in accepted if item.trusted_direct)
            pool = trusted or accepted
            nonredundant = _remove_ancestor_candidates(candidates=pool, profile=profile)
            if len(nonredundant) == 1:
                selected = nonredundant[0]
                status = "ACCEPTED"
                reason = selected.reason
            else:
                ordered = tuple(
                    sorted(
                        nonredundant,
                        key=lambda item: (-item.score, -item.priority, item.label_id),
                    )
                )
                margin = ordered[0].score - ordered[1].score
                if not trusted and margin >= ruleset.settings.minimum_score_margin:
                    selected = ordered[0]
                    status = "ACCEPTED"
                    reason = "Top non-ancestral candidate exceeded the configured conflict margin."
                else:
                    status = "AMBIGUOUS"
                    conflicts = tuple(sorted({item.label_id for item in ordered}))
                    reason = "Conflicting accepted labels were not resolved automatically."
            for item in values:
                key = (item.label_id, item.rule_id)
                if selected is not None and item == selected:
                    decisions[key] = "ACCEPTED"
                elif item.accepted and _is_ancestor_label(
                    ancestor=item.label_id,
                    descendant=(selected.label_id if selected is not None else ""),
                    profile=profile,
                ):
                    decisions[key] = "SUPERSEDED_BY_DESCENDANT"
                elif item.accepted:
                    decisions[key] = "AMBIGUOUS_CONFLICT" if selected is None else "LOWER_PRIORITY"
                else:
                    decisions[key] = "BELOW_ACCEPTANCE_POLICY"
        elif values:
            status = "PROPOSED"
            reason = values[0].reason
            for item in values:
                decisions[(item.label_id, item.rule_id)] = "BELOW_ACCEPTANCE_POLICY"
        else:
            status = "UNMAPPED"
            reason = "No configured target evidence was observed."
        result[protein_id] = CandidateResolution(
            selected=selected,
            candidates=values,
            decisions=decisions,
            status=status,
            reason=reason,
            conflicting_label_ids=conflicts,
        )
    return result


def _orthology_propagation_candidates(
    *,
    contexts: Mapping[str, ProteinEvidenceContext],
    resolutions: Mapping[str, CandidateResolution],
    memberships: tuple[GroupMembership, ...],
    profile: ProteinProfile,
    ruleset: EvidenceRuleSet,
) -> dict[str, list[LabelCandidate]]:
    """Propose one-generation HOG labels from compatible accepted anchors."""

    if not ruleset.settings.orthology_propagation_enabled or not memberships:
        return {}
    rule_by_label = {rule.label_id: rule for rule in ruleset.rules}
    members_by_group: dict[str, list[str]] = defaultdict(list)
    for membership in memberships:
        members_by_group[_membership_group_key(membership=membership)].append(membership.protein_id)
    result: dict[str, list[LabelCandidate]] = defaultdict(list)
    for group_id, protein_ids in sorted(members_by_group.items()):
        anchors = tuple(
            resolution.selected
            for protein_id in protein_ids
            if (resolution := resolutions[protein_id]).selected is not None
        )
        if len(anchors) < ruleset.settings.minimum_orthology_anchor_proteins:
            continue
        consensus = _common_rule_label(
            label_ids=tuple(item.label_id for item in anchors),
            profile=profile,
            rule_label_ids=frozenset(rule_by_label),
        )
        rule = rule_by_label.get(consensus)
        if rule is None or not rule.propagate_by_orthology or not rule.required_domain_patterns:
            continue
        anchor_digest = sha256_text(
            text="\n".join(sorted(f"{item.protein_id}\t{item.label_id}" for item in anchors))
        )
        for protein_id in sorted(set(protein_ids)):
            if resolutions[protein_id].selected is not None:
                continue
            context = contexts[protein_id]
            matched_domains = _matched_domain_features(
                context=context,
                patterns=rule.required_domain_patterns,
            )
            if len(matched_domains) != len(rule.required_domain_patterns):
                continue
            items: list[EvidenceItem] = [
                EvidenceItem(
                    evidence_group="ORTHOFINDER_HOMOLOGY",
                    evidence_kind="UNANIMOUS_HOG_ANCHORS",
                    value=group_id,
                    score=ruleset.settings.orthology_score,
                    evidence_source="Completed OrthoFinder group membership",
                    evidence_reference=f"anchor_sha256={anchor_digest}",
                )
            ]
            items.extend(
                EvidenceItem(
                    evidence_group="PFAM",
                    evidence_kind="ORTHOLOGY_COMPATIBLE_DOMAIN",
                    value=domain_id,
                    score=(ruleset.settings.pfam_score if index == 0 else 0.0),
                    evidence_source="Pfam domain assessment",
                    evidence_reference="domains.tsv",
                    feature_type="DOMAIN",
                    feature_id=domain_id,
                    feature_name=context.domain_names[domain_id],
                )
                for index, domain_id in enumerate(dict.fromkeys(matched_domains))
            )
            score = sum(item.score for item in items)
            structure_ok = _structure_acceptance(
                context=context,
                settings=ruleset.settings,
            )
            accepted = score >= ruleset.settings.minimum_acceptance_score and structure_ok
            result[protein_id].append(
                LabelCandidate(
                    protein_id=protein_id,
                    label_id=rule.label_id,
                    rule_id=f"{rule.rule_id}_hog_propagation",
                    component_role=rule.component_role,
                    score=score,
                    evidence_groups=("ORTHOFINDER_HOMOLOGY", "PFAM"),
                    evidence_items=tuple(items),
                    accepted=accepted,
                    confidence_tier="HIGH" if accepted else "PROPOSED",
                    reason=(
                        "Unanimous HOG anchors plus compatible required domain architecture."
                        if accepted
                        else "HOG/domain evidence found, but target structure is ineligible."
                    ),
                    priority=rule.priority,
                )
            )
    return result


def _select_matched_controls(
    *,
    contexts: Mapping[str, ProteinEvidenceContext],
    resolutions: Mapping[str, CandidateResolution],
    selected_targets: Mapping[str, LabelCandidate],
    profile: ProteinProfile,
    ruleset: EvidenceRuleSet,
) -> tuple[tuple[dict[str, Any], ...], Mapping[str, frozenset[str]]]:
    """Select deterministic background units independently of observed signatures."""

    comparison_rows = default_profile_comparisons(profile=profile)
    targets_by_background: dict[str, set[str]] = defaultdict(set)
    target_labels_by_background: dict[str, set[str]] = defaultdict(set)
    for protein_id, candidate in selected_targets.items():
        ancestors = frozenset(label_ancestors(profile=profile, label_id=candidate.label_id))
        for comparison in comparison_rows:
            target_label = comparison.target_label_ids[0]
            if target_label in ancestors:
                background = comparison.background_label_ids[0]
                targets_by_background[background].add(protein_id)
                target_labels_by_background[background].add(target_label)
    uncertain_ids = frozenset(
        protein_id
        for protein_id, resolution in resolutions.items()
        if resolution.candidates or resolution.selected is not None
    )
    target_units = frozenset(
        contexts[protein_id].independence_unit for protein_id in selected_targets
    )
    clean_candidates = tuple(
        context
        for protein_id, context in contexts.items()
        if protein_id not in uncertain_ids
        and context.independence_unit not in target_units
        and not (
            ruleset.settings.exclude_input_candidates_from_controls and context.input_candidate
        )
    )
    if not clean_candidates:
        raise InputValidationError(
            "No clean control candidates remain after target, ambiguity and homology exclusions."
        )
    rows: list[dict[str, Any]] = []
    controls: dict[str, frozenset[str]] = {}
    for background_label, target_ids in sorted(targets_by_background.items()):
        matches, selected = _match_background_units(
            background_label=background_label,
            target_label_ids=tuple(sorted(target_labels_by_background[background_label])),
            target_ids=frozenset(target_ids),
            contexts=contexts,
            candidates=clean_candidates,
            settings=ruleset.settings,
        )
        matched_units = {contexts[protein_id].independence_unit for protein_id in selected}
        if len(matched_units) < ruleset.settings.minimum_control_units_per_background:
            raise InputValidationError(
                f"Background {background_label!r} has only {len(matched_units)} matched "
                "independent units; revise inputs or a versioned matching policy rather than "
                "using an underpowered or unmatched background."
            )
        rows.extend(matches)
        controls[background_label] = frozenset(selected)
    return tuple(rows), controls


def _match_background_units(
    *,
    background_label: str,
    target_label_ids: tuple[str, ...],
    target_ids: frozenset[str],
    contexts: Mapping[str, ProteinEvidenceContext],
    candidates: tuple[ProteinEvidenceContext, ...],
    settings: EvidenceSettings,
) -> tuple[tuple[dict[str, Any], ...], frozenset[str]]:
    """Greedily match one representative per target and control block."""

    target_representatives = _unit_representatives(
        contexts=tuple(contexts[protein_id] for protein_id in target_ids)
    )
    candidate_representatives = _unit_representatives(contexts=candidates)
    index = _control_bucket_index(contexts=tuple(candidate_representatives.values()))
    used_units: set[str] = set()
    selected: set[str] = set()
    rows: list[dict[str, Any]] = []
    for target_unit, target in sorted(target_representatives.items()):
        for _ in range(settings.control_units_per_target_unit):
            compatible = _nearby_control_candidates(
                target=target,
                index=index,
                used_units=used_units,
                settings=settings,
            )
            scored = [
                (_control_distance(target=target, control=control), control)
                for control in compatible
                if _control_is_compatible(target=target, control=control, settings=settings)
            ]
            if not scored:
                rows.append(
                    _control_match_row(
                        background_label=background_label,
                        target_label_ids=target_label_ids,
                        target=target,
                        control=None,
                        score=None,
                        status="UNMATCHED",
                        reason="No unused candidate satisfied every prespecified caliper.",
                    )
                )
                continue
            score, control = min(
                scored,
                key=lambda item: (item[0], item[1].protein_id),
            )
            used_units.add(control.independence_unit)
            selected.add(control.protein_id)
            rows.append(
                _control_match_row(
                    background_label=background_label,
                    target_label_ids=target_label_ids,
                    target=target,
                    control=control,
                    score=score,
                    status="MATCHED",
                    reason="Nearest unused independent unit within every configured caliper.",
                )
            )
    return tuple(rows), frozenset(selected)


def _control_bucket_index(
    *, contexts: tuple[ProteinEvidenceContext, ...]
) -> Mapping[tuple[str, bool, int], tuple[tuple[int, str, ProteinEvidenceContext], ...]]:
    """Index control candidates by species, structure state and domain count."""

    buckets: dict[
        tuple[str, bool, int],
        list[tuple[int, str, ProteinEvidenceContext]],
    ] = defaultdict(list)
    for context in contexts:
        species_values = context.species or ("",)
        for species in species_values:
            buckets[(species, context.structure_eligible, context.domain_count)].append(
                (context.sequence_length, context.protein_id, context)
            )
    return {
        key: tuple(sorted(values, key=lambda item: (item[0], item[1])))
        for key, values in buckets.items()
    }


def _nearby_control_candidates(
    *,
    target: ProteinEvidenceContext,
    index: Mapping[
        tuple[str, bool, int],
        tuple[tuple[int, str, ProteinEvidenceContext], ...],
    ],
    used_units: Collection[str],
    settings: EvidenceSettings,
    maximum_per_bucket: int = 128,
) -> tuple[ProteinEvidenceContext, ...]:
    """Return bounded nearest-length candidates from compatible matching buckets."""

    species_values = (
        target.species
        if settings.require_species_match
        else tuple(sorted({key[0] for key in index}))
    )
    structures = (target.structure_eligible,) if settings.require_structure_match else (False, True)
    values: dict[str, ProteinEvidenceContext] = {}
    minimum_domains = max(0, target.domain_count - settings.maximum_domain_count_difference)
    maximum_domains = target.domain_count + settings.maximum_domain_count_difference
    for species in species_values or ("",):
        for structure_state in structures:
            for domain_count in range(minimum_domains, maximum_domains + 1):
                bucket = index.get((species, structure_state, domain_count), ())
                if not bucket:
                    continue
                lengths = [item[0] for item in bucket]
                insertion = bisect.bisect_left(lengths, target.sequence_length)
                left = insertion - 1
                right = insertion
                accepted = 0
                while accepted < maximum_per_bucket and (left >= 0 or right < len(bucket)):
                    choose_left = right >= len(bucket) or (
                        left >= 0
                        and abs(bucket[left][0] - target.sequence_length)
                        <= abs(bucket[right][0] - target.sequence_length)
                    )
                    item = bucket[left] if choose_left else bucket[right]
                    left -= int(choose_left)
                    right += int(not choose_left)
                    context = item[2]
                    if context.independence_unit in used_units:
                        continue
                    values.setdefault(context.protein_id, context)
                    accepted += 1
    return tuple(values[key] for key in sorted(values))


def _control_is_compatible(
    *,
    target: ProteinEvidenceContext,
    control: ProteinEvidenceContext,
    settings: EvidenceSettings,
) -> bool:
    """Apply every prespecified background matching caliper."""

    if target.independence_unit == control.independence_unit:
        return False
    if settings.require_species_match and not set(target.species).intersection(control.species):
        return False
    if settings.require_structure_match and target.structure_eligible != control.structure_eligible:
        return False
    length_difference = abs(math.log2((target.sequence_length + 1) / (control.sequence_length + 1)))
    if length_difference > settings.maximum_log2_length_difference:
        return False
    return (
        abs(target.domain_count - control.domain_count) <= settings.maximum_domain_count_difference
    )


def _control_distance(*, target: ProteinEvidenceContext, control: ProteinEvidenceContext) -> float:
    """Calculate an outcome-blind multivariable matching distance."""

    length_difference = abs(math.log2((target.sequence_length + 1) / (control.sequence_length + 1)))
    domain_difference = abs(target.domain_count - control.domain_count)
    target_domains = set(target.domain_ids)
    control_domains = set(control.domain_ids)
    union = target_domains | control_domains
    jaccard = len(target_domains & control_domains) / len(union) if union else 1.0
    if target.mean_confidence is None or control.mean_confidence is None:
        confidence_difference = 0.0
    else:
        confidence_difference = abs(target.mean_confidence - control.mean_confidence) / 100.0
    return (
        2.0 * length_difference
        + 0.25 * domain_difference
        + 0.5 * (1.0 - jaccard)
        + 0.5 * confidence_difference
    )


def _control_match_row(
    *,
    background_label: str,
    target_label_ids: tuple[str, ...],
    target: ProteinEvidenceContext,
    control: ProteinEvidenceContext | None,
    score: float | None,
    status: str,
    reason: str,
) -> dict[str, Any]:
    """Build one complete control-match audit row."""

    if control is None:
        length_difference = None
        domain_difference = None
        jaccard = None
        confidence_difference = None
        species_match = False
        structure_match = False
    else:
        length_difference = abs(
            math.log2((target.sequence_length + 1) / (control.sequence_length + 1))
        )
        domain_difference = abs(target.domain_count - control.domain_count)
        union = set(target.domain_ids) | set(control.domain_ids)
        jaccard = (
            len(set(target.domain_ids) & set(control.domain_ids)) / len(union) if union else 1.0
        )
        confidence_difference = (
            abs(target.mean_confidence - control.mean_confidence)
            if target.mean_confidence is not None and control.mean_confidence is not None
            else None
        )
        species_match = bool(set(target.species).intersection(control.species))
        structure_match = target.structure_eligible == control.structure_eligible
    return {
        "background_label_id": background_label,
        "target_label_ids": "|".join(target_label_ids),
        "target_protein_id": target.protein_id,
        "target_unit_id": target.independence_unit,
        "control_protein_id": control.protein_id if control is not None else "",
        "control_unit_id": control.independence_unit if control is not None else "",
        "species_match": species_match,
        "structure_eligibility_match": structure_match,
        "length_log2_difference": length_difference,
        "domain_count_difference": domain_difference,
        "domain_architecture_jaccard": jaccard,
        "mean_confidence_difference": confidence_difference,
        "match_score": score,
        "status": status,
        "reason": reason,
    }


def _unit_representatives(
    *, contexts: Iterable[ProteinEvidenceContext]
) -> Mapping[str, ProteinEvidenceContext]:
    """Choose one deterministic best-structure representative per independence unit."""

    result: dict[str, ProteinEvidenceContext] = {}
    for context in contexts:
        previous = result.get(context.independence_unit)
        candidate_key = _representative_rank(context=context)
        previous_key = _representative_rank(context=previous) if previous is not None else None
        if previous is None or candidate_key < previous_key:
            result[context.independence_unit] = context
    return result


def _representative_rank(*, context: ProteinEvidenceContext) -> tuple[float, float, int, str]:
    """Return an ascending deterministic rank for a block representative.

    Args:
        context: Candidate protein context.

    Returns:
        Rank favouring eligible structures, high confidence, shorter sequences,
        and finally the lexicographically smallest identifier.
    """

    confidence = context.mean_confidence if context.mean_confidence is not None else -1.0
    return (
        -float(context.structure_eligible),
        -confidence,
        context.sequence_length,
        context.protein_id,
    )


def _build_label_rows(
    *,
    sequences: tuple[SequenceRecord, ...],
    profile: ProteinProfile,
    template_labels: Path | None,
    resolutions: Mapping[str, CandidateResolution],
    selected_targets: Mapping[str, LabelCandidate],
    controls_by_label: Mapping[str, frozenset[str]],
    ruleset: EvidenceRuleSet,
) -> tuple[tuple[dict[str, str], ...], tuple[dict[str, str], ...]]:
    """Create complete assignment coverage and the unresolved review queue."""

    protein_ids = frozenset(record.protein_id for record in sequences)
    base = _base_label_records(
        sequences=sequences,
        profile=profile,
        template_labels=template_labels,
    )
    controls_by_protein: dict[str, list[str]] = defaultdict(list)
    for label_id, values in controls_by_label.items():
        for protein_id in values:
            controls_by_protein[protein_id].append(label_id)
    rows: list[dict[str, str]] = []
    unresolved: list[dict[str, str]] = []
    for protein_id in sorted(protein_ids):
        selected = selected_targets.get(protein_id)
        if selected is not None:
            evidence_payload = json.dumps(
                [item.to_record() for item in selected.evidence_items],
                sort_keys=True,
                separators=(",", ":"),
            )
            rows.append(
                {
                    "protein_id": protein_id,
                    "label_id": selected.label_id,
                    "curation_status": EVIDENCE_POSITIVE_STATUS,
                    "evidence_status": EVIDENCE_TARGET_STATUS,
                    "evidence_source": (
                        f"{ruleset.ruleset_id} evidence rules {ruleset.ruleset_version}"
                    ),
                    "evidence_reference": (
                        "evidence_sha256="
                        + hashlib.sha256(evidence_payload.encode("utf-8")).hexdigest()
                    ),
                    "component_role": selected.component_role,
                    "curation_reason": (
                        f"Automated provisional assignment; {selected.confidence_tier}; "
                        f"score={selected.score:.6g}; {selected.reason}"
                    ),
                }
            )
            continue
        if protein_id in controls_by_protein:
            for label_id in sorted(controls_by_protein[protein_id]):
                rows.append(
                    {
                        "protein_id": protein_id,
                        "label_id": label_id,
                        "curation_status": EVIDENCE_POSITIVE_STATUS,
                        "evidence_status": EVIDENCE_CONTROL_STATUS,
                        "evidence_source": (
                            f"{ruleset.ruleset_id} prespecified covariate matching"
                        ),
                        "evidence_reference": "control_matching_audit.tsv",
                        "component_role": "CONTROL",
                        "curation_reason": (
                            "Outcome-blind matched control selected from a distinct "
                            "homology/redundancy unit."
                        ),
                    }
                )
            continue
        resolution = resolutions[protein_id]
        row = dict(base[protein_id])
        if resolution.status == "AMBIGUOUS":
            row.update(
                {
                    "curation_status": "AMBIGUOUS",
                    "evidence_status": "AUTOMATED_EVIDENCE_CONFLICT",
                    "curation_reason": resolution.reason,
                }
            )
        elif resolution.status == "PROPOSED":
            row.update(
                {
                    "curation_status": "PROPOSED",
                    "evidence_status": "AUTOMATED_EVIDENCE_BELOW_POLICY",
                    "curation_reason": resolution.reason,
                }
            )
        rows.append(row)
        if resolution.status != "UNMAPPED":
            best_score = max((item.score for item in resolution.candidates), default=0.0)
            unresolved.append(
                {
                    "protein_id": protein_id,
                    "provisional_label_id": (
                        resolution.candidates[0].label_id if resolution.candidates else ""
                    ),
                    "curation_status": row["curation_status"],
                    "best_score": format(best_score, ".15g"),
                    "candidate_label_ids": "|".join(
                        sorted({item.label_id for item in resolution.candidates})
                    ),
                    "reason": resolution.reason,
                }
            )
    return tuple(rows), tuple(unresolved)


def _base_label_records(
    *,
    sequences: tuple[SequenceRecord, ...],
    profile: ProteinProfile,
    template_labels: Path | None,
) -> Mapping[str, dict[str, str]]:
    """Load a safe template or build one non-positive row per protein."""

    protein_ids = frozenset(record.protein_id for record in sequences)
    if template_labels is not None:
        assignments = read_label_assignments(
            path=template_labels,
            protein_ids=protein_ids,
            label_ids=profile.label_ids(),
        )
        by_protein: dict[str, list[LabelAssignment]] = defaultdict(list)
        for assignment in assignments:
            if assignment.is_eligible_positive:
                raise InputValidationError(
                    "Evidence-label templates must not contain analysis-positive assignments."
                )
            by_protein[assignment.protein_id].append(assignment)
        if frozenset(by_protein) != protein_ids:
            raise InputValidationError("Evidence-label template must cover every FASTA protein.")
        return {
            protein_id: sorted(
                (assignment.to_record() for assignment in values),
                key=lambda row: row["label_id"],
            )[0]
            for protein_id, values in by_protein.items()
        }
    roots = sorted(label.label_id for label in profile.labels if not label.parent_label_id)
    if len(roots) != 1:
        raise InputValidationError("A template-free profile must have exactly one root label.")
    return {
        protein_id: {
            "protein_id": protein_id,
            "label_id": roots[0],
            "curation_status": "UNMAPPED",
            "evidence_status": "AUTOMATED_EVIDENCE_NOT_MAPPED",
            "evidence_source": "protein-signature-analysis evidence labelling",
            "evidence_reference": "",
            "component_role": "UNKNOWN",
            "curation_reason": "No configured evidence-supported assignment.",
        }
        for protein_id in sorted(protein_ids)
    }


def _candidate_audit_rows(
    *,
    contexts: Mapping[str, ProteinEvidenceContext],
    resolutions: Mapping[str, CandidateResolution],
) -> tuple[dict[str, Any], ...]:
    """Serialise every scored candidate and explicit decision."""

    rows: list[dict[str, Any]] = []
    for protein_id, resolution in resolutions.items():
        for candidate in resolution.candidates:
            rows.append(
                {
                    "protein_id": protein_id,
                    "label_id": candidate.label_id,
                    "rule_id": candidate.rule_id,
                    "decision": resolution.decisions.get(
                        (candidate.label_id, candidate.rule_id), "BELOW_ACCEPTANCE_POLICY"
                    ),
                    "confidence_tier": candidate.confidence_tier,
                    "score": candidate.score,
                    "evidence_group_count": len(candidate.evidence_groups),
                    "evidence_groups": "|".join(candidate.evidence_groups),
                    "evidence_items": json.dumps(
                        [item.to_record() for item in candidate.evidence_items],
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                    "independence_unit": contexts[protein_id].independence_unit,
                    "conflicting_label_ids": "|".join(resolution.conflicting_label_ids),
                    "reason": candidate.reason,
                }
            )
    return tuple(sorted(rows, key=lambda row: (row["protein_id"], row["label_id"], row["rule_id"])))


def _label_definition_rows(
    *,
    selected_targets: Mapping[str, LabelCandidate],
    profile: ProteinProfile,
) -> tuple[dict[str, str], ...]:
    """List domain features used for labels and excluded from confirmatory analysis."""

    rows: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for candidate in selected_targets.values():
        for item in candidate.evidence_items:
            if not item.feature_type or not item.feature_id:
                continue
            for label_id in label_ancestors(profile=profile, label_id=candidate.label_id):
                key = (label_id, item.feature_type, item.feature_id, candidate.rule_id)
                rows[key] = {
                    "label_id": label_id,
                    "feature_type": item.feature_type,
                    "feature_id": item.feature_id,
                    "feature_name": item.feature_name or item.feature_id,
                    "rule_id": candidate.rule_id,
                    "evidence_role": "LABEL_DEFINING",
                    "exclusion_scope": "GLOBAL_CONFIRMATORY_DOMAIN_AND_ML",
                    "reason": (
                        "Feature contributed to automated target assignment and is excluded "
                        "from confirmatory domain association and explainable modelling."
                    ),
                }
    return tuple(rows[key] for key in sorted(rows))


def _class_summary_rows(
    *,
    profile: ProteinProfile,
    selected_targets: Mapping[str, LabelCandidate],
    controls_by_label: Mapping[str, frozenset[str]],
    contexts: Mapping[str, ProteinEvidenceContext],
) -> tuple[dict[str, Any], ...]:
    """Summarise direct target and matched-control class coverage."""

    comparisons = default_profile_comparisons(profile=profile)
    background_by_target = {
        item.target_label_ids[0]: item.background_label_ids[0] for item in comparisons
    }
    target_counts = Counter(item.label_id for item in selected_targets.values())
    target_units: dict[str, set[str]] = defaultdict(set)
    for protein_id, item in selected_targets.items():
        target_units[item.label_id].add(contexts[protein_id].independence_unit)
    rows: list[dict[str, Any]] = []
    for label_id in sorted(target_counts):
        background = background_by_target.get(label_id, "")
        control_ids = controls_by_label.get(background, frozenset())
        rows.append(
            {
                "label_id": label_id,
                "label_type": "TARGET",
                "direct_positive_protein_count": target_counts[label_id],
                "independent_unit_count": len(target_units[label_id]),
                "matched_background_label_id": background,
                "matched_control_protein_count": len(control_ids),
                "matched_control_unit_count": len(
                    {contexts[protein_id].independence_unit for protein_id in control_ids}
                ),
                "status": "PROVISIONAL_EVIDENCE_SUPPORTED",
            }
        )
    for label_id, protein_ids in sorted(controls_by_label.items()):
        rows.append(
            {
                "label_id": label_id,
                "label_type": "BACKGROUND",
                "direct_positive_protein_count": len(protein_ids),
                "independent_unit_count": len(
                    {contexts[protein_id].independence_unit for protein_id in protein_ids}
                ),
                "matched_background_label_id": "",
                "matched_control_protein_count": 0,
                "matched_control_unit_count": 0,
                "status": "PROVISIONAL_MATCHED_CONTROL",
            }
        )
    return tuple(rows)


def _filter_label_defining_domains(
    *,
    rows: tuple[dict[str, str], ...],
    protein_ids: frozenset[str],
    excluded_domain_ids: frozenset[str],
    evidence_reference: str,
) -> tuple[dict[str, str], ...]:
    """Remove label-defining Pfam features while retaining assessment coverage."""

    by_protein: dict[str, list[dict[str, str]]] = defaultdict(list)
    assessed_status: dict[str, str] = {}
    for row in rows:
        protein_id = row["protein_id"]
        status = row["assessment_status"]
        assessed_status[protein_id] = status
        domain_id = row["domain_id"]
        if domain_id and domain_id in excluded_domain_ids:
            continue
        if domain_id:
            by_protein[protein_id].append(dict(row))
        elif status != "ASSESSED_WITH_HIT":
            by_protein[protein_id].append(dict(row))
    output: list[dict[str, str]] = []
    for protein_id in sorted(protein_ids):
        retained = by_protein.get(protein_id, [])
        if retained:
            output.extend(retained)
            continue
        original_status = assessed_status.get(protein_id, "NOT_ASSESSED")
        status = "ASSESSED_NO_HIT" if original_status == "ASSESSED_WITH_HIT" else original_status
        output.append(
            {
                "protein_id": protein_id,
                "domain_authority": "Pfam",
                "assessment_status": status,
                "domain_id": "",
                "domain_name": "",
                "start": "",
                "end": "",
                "score": "",
                "e_value": "",
                "evidence_source": "Label-definition-safe Pfam projection",
                "evidence_reference": evidence_reference,
            }
        )
    return tuple(
        sorted(
            output,
            key=lambda row: (
                row["protein_id"],
                row["domain_id"],
                int(row["start"] or 0),
                int(row["end"] or 0),
            ),
        )
    )


def _matched_domain_features(
    *, context: ProteinEvidenceContext, patterns: tuple[re.Pattern[str], ...]
) -> tuple[str, ...]:
    """Resolve one matching domain for every required regular expression."""

    matches: list[str] = []
    for pattern in patterns:
        candidates = tuple(
            domain_id
            for domain_id in context.domain_ids
            if pattern.search(f"{domain_id} {context.domain_names[domain_id]}")
        )
        if not candidates:
            continue
        matches.append(sorted(candidates)[0])
    return tuple(matches)


def _common_rule_label(
    *, label_ids: tuple[str, ...], profile: ProteinProfile, rule_label_ids: frozenset[str]
) -> str:
    """Return the deepest common rule-bearing ancestor or an empty value."""

    if not label_ids:
        return ""
    common = set(label_ancestors(profile=profile, label_id=label_ids[0]))
    for label_id in label_ids[1:]:
        common.intersection_update(label_ancestors(profile=profile, label_id=label_id))
    eligible = common & set(rule_label_ids)
    if not eligible:
        return ""
    return max(
        eligible,
        key=lambda label_id: len(label_ancestors(profile=profile, label_id=label_id)),
    )


def _remove_ancestor_candidates(
    *, candidates: tuple[LabelCandidate, ...], profile: ProteinProfile
) -> tuple[LabelCandidate, ...]:
    """Discard a broad candidate when an accepted descendant is also present."""

    result = []
    for candidate in candidates:
        if any(
            candidate.label_id != other.label_id
            and _is_ancestor_label(
                ancestor=candidate.label_id,
                descendant=other.label_id,
                profile=profile,
            )
            for other in candidates
        ):
            continue
        result.append(candidate)
    return tuple(result)


def _is_ancestor_label(*, ancestor: str, descendant: str, profile: ProteinProfile) -> bool:
    """Return whether one profile label is an ancestor of another."""

    if not descendant:
        return False
    return ancestor in label_ancestors(profile=profile, label_id=descendant)[1:]


def _membership_group_key(*, membership: GroupMembership) -> str:
    """Create one unambiguous stable OrthoFinder membership key."""

    return "|".join(
        (
            membership.run_id,
            membership.group_type,
            membership.hierarchy_node,
            membership.group_id,
        )
    )


def _evidence_item_sort_key(item: EvidenceItem) -> tuple[str, str, str, str]:
    """Return deterministic ordering fields for evidence items."""

    return (
        item.evidence_group,
        item.evidence_kind,
        item.feature_id,
        item.value,
    )


def _compile_patterns(*, value: Any, field_name: str) -> tuple[re.Pattern[str], ...]:
    """Compile a unique bounded list of case-insensitive regular expressions."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise InputValidationError(f"{field_name} must be a list of regular expressions.")
    patterns: list[re.Pattern[str]] = []
    raw_values: list[str] = []
    for index, raw in enumerate(value):
        text = validate_text(
            value=raw,
            field_name=f"{field_name}[{index}]",
            maximum_length=500,
        )
        if text in raw_values:
            raise InputValidationError(f"{field_name} contains duplicate patterns.")
        try:
            patterns.append(re.compile(text, flags=re.IGNORECASE))
        except re.error as error:
            raise InputValidationError(
                f"Invalid regular expression in {field_name}: {text!r}: {error}"
            ) from error
        raw_values.append(text)
    return tuple(patterns)


def _reject_unknown_fields(
    *, value: Mapping[str, Any], allowed: Collection[str], field_name: str
) -> None:
    """Reject unexpected mapping keys at a public configuration boundary."""

    unknown = sorted(set(value) - set(allowed))
    if unknown:
        raise InputValidationError(f"Unknown fields in {field_name}: {unknown}")


def _finite_nonnegative_number(*, value: Any, field_name: str) -> float:
    """Parse one finite non-negative numeric ruleset value."""

    if isinstance(value, bool):
        raise InputValidationError(f"{field_name} must be a finite non-negative number.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise InputValidationError(f"{field_name} must be a finite non-negative number.") from error
    if not math.isfinite(number) or number < 0.0:
        raise InputValidationError(f"{field_name} must be a finite non-negative number.")
    return number


def _positive_integer(*, value: Any, field_name: str) -> int:
    """Validate one strictly positive integer ruleset value."""

    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise InputValidationError(f"{field_name} must be a positive integer.")
    return value


def _nonnegative_integer(*, value: Any, field_name: str) -> int:
    """Validate one non-negative integer ruleset value."""

    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InputValidationError(f"{field_name} must be a non-negative integer.")
    return value


def _audit_integer(*, value: str, field_name: str, row_number: int) -> int | None:
    """Parse one optional integer cell from a generated audit table.

    Args:
        value: Raw TSV value.
        field_name: Column name for diagnostics.
        row_number: One-based source row number.

    Returns:
        Parsed integer or ``None`` for an empty cell.

    Raises:
        InputValidationError: If the cell is not an integer.
    """

    if not value.strip():
        return None
    try:
        number = int(value)
    except ValueError as error:
        raise InputValidationError(
            f"Invalid integer in {field_name!r} at row {row_number}: {value!r}"
        ) from error
    if number < 0:
        raise InputValidationError(
            f"Negative integer in {field_name!r} at row {row_number}: {value!r}"
        )
    return number


def _audit_float(*, value: str, field_name: str, row_number: int) -> float | None:
    """Parse one optional finite floating-point audit cell.

    Args:
        value: Raw TSV value.
        field_name: Column name for diagnostics.
        row_number: One-based source row number.

    Returns:
        Parsed finite number or ``None`` for an empty cell.

    Raises:
        InputValidationError: If the cell is not finite numeric text.
    """

    if not value.strip():
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise InputValidationError(
            f"Invalid number in {field_name!r} at row {row_number}: {value!r}"
        ) from error
    if not math.isfinite(number):
        raise InputValidationError(
            f"Non-finite number in {field_name!r} at row {row_number}: {value!r}"
        )
    return number


def _audit_boolean(*, value: str, field_name: str, row_number: int) -> bool:
    """Parse one strict Boolean audit cell.

    Args:
        value: Raw TSV value.
        field_name: Column name for diagnostics.
        row_number: One-based source row number.

    Returns:
        Parsed Boolean.

    Raises:
        InputValidationError: If the cell is not TRUE or FALSE.
    """

    normalised = value.strip().upper()
    if normalised not in {"TRUE", "FALSE"}:
        raise InputValidationError(
            f"Invalid Boolean in {field_name!r} at row {row_number}: {value!r}"
        )
    return normalised == "TRUE"


def _boolean(*, value: Any, field_name: str) -> bool:
    """Validate one strict Boolean ruleset value."""

    if not isinstance(value, bool):
        raise InputValidationError(f"{field_name} must be true or false.")
    return value


def _split_pipe(*, value: str) -> tuple[str, ...]:
    """Split a deterministic pipe-delimited review cell."""

    return tuple(sorted({item.strip() for item in str(value).split("|") if item.strip()}))


def _required_file(*, path: Path, field_name: str) -> Path:
    """Resolve one required non-empty input file."""

    candidate = Path(path).expanduser().resolve()
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"{field_name} is missing or empty: {candidate}")
    return candidate


def _optional_file(*, path: Path | None, field_name: str) -> Path | None:
    """Resolve one optional non-empty input file."""

    return None if path is None else _required_file(path=path, field_name=field_name)


def _file_inventory(*, paths: Sequence[Path]) -> tuple[dict[str, Any], ...]:
    """Build an exact checksummed input-file inventory."""

    return tuple(
        {
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in sorted({Path(item).resolve() for item in paths}, key=str)
    )


def _directory_inventory(*, root: Path) -> list[dict[str, Any]]:
    """Inventory every current file below a staging directory."""

    return [
        {
            "relative_path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path=path),
        }
        for path in sorted((item for item in root.rglob("*") if item.is_file()), key=str)
    ]
