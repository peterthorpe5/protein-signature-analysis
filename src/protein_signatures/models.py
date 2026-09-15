"""Immutable models shared across analysis stages."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any


class AnalysisStatus(StrEnum):
    """Controlled assessment and availability states."""

    COMPLETE = "COMPLETE"
    NOT_SELECTED = "NOT_SELECTED"
    INPUT_UNAVAILABLE = "INPUT_UNAVAILABLE"
    INSUFFICIENT_SAMPLE_SIZE = "INSUFFICIENT_SAMPLE_SIZE"
    NO_ELIGIBLE_FEATURES = "NO_ELIGIBLE_FEATURES"
    NO_SIGNIFICANT_SIGNATURE = "NO_SIGNIFICANT_SIGNATURE"
    FAILED = "FAILED"
    EXCLUDED = "EXCLUDED"
    AMBIGUOUS_CLASS = "AMBIGUOUS_CLASS"
    MAPPING_FAILED = "MAPPING_FAILED"


class CurationStatus(StrEnum):
    """Controlled states for protein-to-label assignments."""

    REVIEWED_POSITIVE = "REVIEWED_POSITIVE"
    EVIDENCE_SUPPORTED_POSITIVE = "EVIDENCE_SUPPORTED_POSITIVE"
    REVIEWED_NEGATIVE = "REVIEWED_NEGATIVE"
    REVIEWED_COMPONENT_NOT_CATALYTIC = "REVIEWED_COMPONENT_NOT_CATALYTIC"
    PROPOSED = "PROPOSED"
    AMBIGUOUS = "AMBIGUOUS"
    UNMAPPED = "UNMAPPED"
    EXCLUDED = "EXCLUDED"


class DomainAssessmentStatus(StrEnum):
    """Controlled states for a protein/domain-authority assessment."""

    ASSESSED_WITH_HIT = "ASSESSED_WITH_HIT"
    ASSESSED_NO_HIT = "ASSESSED_NO_HIT"
    NOT_ASSESSED = "NOT_ASSESSED"
    FAILED = "FAILED"


class StructureAnalysisEligibility(StrEnum):
    """Controlled eligibility states for local coordinate-model analysis."""

    ELIGIBLE = "ELIGIBLE"
    NOT_APPLICABLE_EXTERNAL_EVIDENCE = "NOT_APPLICABLE_EXTERNAL_EVIDENCE"
    INELIGIBLE_LOW_CONFIDENCE = "INELIGIBLE_LOW_CONFIDENCE"
    INELIGIBLE_CONFIDENCE_UNAVAILABLE = "INELIGIBLE_CONFIDENCE_UNAVAILABLE"
    INELIGIBLE_SEQUENCE_UNVERIFIED = "INELIGIBLE_SEQUENCE_UNVERIFIED"
    INELIGIBLE_SEQUENCE_MISMATCH = "INELIGIBLE_SEQUENCE_MISMATCH"
    INELIGIBLE_COORDINATE_UNAVAILABLE = "INELIGIBLE_COORDINATE_UNAVAILABLE"
    INELIGIBLE_USER_EXCLUDED = "INELIGIBLE_USER_EXCLUDED"
    INELIGIBLE_VALIDATION_FAILED = "INELIGIBLE_VALIDATION_FAILED"


class StructureCoverageScope(StrEnum):
    """Controlled denominator scopes for structural-alignment coverage."""

    FULL_SEQUENCE = "FULL_SEQUENCE"
    STRUCTURE_MODEL_RESIDUES = "STRUCTURE_MODEL_RESIDUES"
    DOMAIN_OR_CONSTRUCT = "DOMAIN_OR_CONSTRUCT"


class FoldEvidenceStatus(StrEnum):
    """Controlled assessment states for one fold authority."""

    ASSESSED_WITH_HIT = "ASSESSED_WITH_HIT"
    ASSESSED_NO_HIT = "ASSESSED_NO_HIT"
    NOT_ASSESSED = "NOT_ASSESSED"
    FAILED = "FAILED"
    INPUT_UNAVAILABLE = "INPUT_UNAVAILABLE"


class StructureComparisonStatus(StrEnum):
    """Controlled completion states for pairwise structural evidence."""

    COMPLETE = "COMPLETE"
    SUCCESS = "SUCCESS"
    PASS = "PASS"
    FAILED = "FAILED"
    NOT_ASSESSED = "NOT_ASSESSED"
    INPUT_UNAVAILABLE = "INPUT_UNAVAILABLE"
    EXCLUDED = "EXCLUDED"


@dataclass(frozen=True)
class SequenceRecord:
    """One validated protein sequence."""

    protein_id: str
    description: str
    sequence: str
    sequence_length: int
    sequence_sha256: str

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable protein inventory row.

        Returns:
            Dictionary representation of this record.
        """

        return asdict(self)


@dataclass(frozen=True)
class LabelAssignment:
    """One evidence-bearing assignment of a protein to a controlled label."""

    protein_id: str
    label_id: str
    curation_status: CurationStatus
    evidence_status: str
    evidence_source: str
    evidence_reference: str
    component_role: str
    curation_reason: str

    @property
    def is_eligible_positive(self) -> bool:
        """Return whether this assignment can enter positive analysis sets.

        Returns:
            ``True`` for human-reviewed or explicitly evidence-supported
            provisional positive assignments.
        """

        return self.curation_status in {
            CurationStatus.REVIEWED_POSITIVE,
            CurationStatus.EVIDENCE_SUPPORTED_POSITIVE,
        }

    def to_record(self) -> dict[str, str]:
        """Return a serialisable label-assignment row.

        Returns:
            Dictionary representation of this assignment.
        """

        row = asdict(self)
        row["curation_status"] = self.curation_status.value
        return row


@dataclass(frozen=True)
class FeatureRecord:
    """One sequence, domain, fold, pocket or structure-derived feature."""

    protein_id: str
    feature_type: str
    feature_id: str
    feature_name: str
    start: int | None
    end: int | None
    evidence_status: str
    evidence_source: str
    evidence_reference: str
    derivation_scope: str = ""
    feature_definition_sha256: str = ""
    derivation_cohort_sha256: str = ""

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable feature row.

        Returns:
            Dictionary representation of this feature.
        """

        return asdict(self)


@dataclass(frozen=True)
class DomainHit:
    """One coordinate-resolved domain annotation from a named authority."""

    protein_id: str
    domain_authority: str
    domain_id: str
    domain_name: str
    start: int
    end: int
    score: float | None
    e_value: float | None
    evidence_source: str
    evidence_reference: str

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable domain-hit row.

        Returns:
            Dictionary representation of this hit.
        """

        return asdict(self)


@dataclass(frozen=True)
class DomainAssessment:
    """Assessment coverage for one protein and domain authority."""

    protein_id: str
    domain_authority: str
    assessment_status: DomainAssessmentStatus
    hit_count: int
    evidence_source: str
    evidence_reference: str

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable domain-assessment row.

        Returns:
            Dictionary representation with controlled status text.
        """

        row = asdict(self)
        row["assessment_status"] = self.assessment_status.value
        return row


@dataclass(frozen=True)
class StructureRecord:
    """One structure model and its optional fold assignment."""

    protein_id: str
    structure_id: str
    structure_source: str
    structure_version: str
    coordinate_path: Path | None
    coordinate_sha256: str
    availability_status: str
    mean_confidence: float | None
    fold_id: str
    fold_name: str
    fold_authority: str
    fold_authority_version: str
    fold_evidence_reference: str
    fold_evidence_status: FoldEvidenceStatus
    analysis_eligibility_status: StructureAnalysisEligibility
    comparison_universe_ids: tuple[str, ...]

    @property
    def is_coordinate_analysis_eligible(self) -> bool:
        """Return whether this record can enter local coordinate analysis.

        Returns:
            ``True`` only for an explicitly eligible, physically available model.
        """

        return (
            self.analysis_eligibility_status == StructureAnalysisEligibility.ELIGIBLE
            and self.coordinate_path is not None
            and self.availability_status in {"AVAILABLE", "COMPLETE"}
        )

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable structure row.

        Returns:
            Dictionary representation with a portable path string.
        """

        row = asdict(self)
        row["coordinate_path"] = str(self.coordinate_path) if self.coordinate_path else ""
        row["fold_evidence_status"] = self.fold_evidence_status.value
        row["analysis_eligibility_status"] = self.analysis_eligibility_status.value
        row["comparison_universe_ids"] = "|".join(self.comparison_universe_ids)
        return row


@dataclass(frozen=True)
class PairwiseStructureComparison:
    """One validated pairwise structural comparison."""

    protein_a_id: str
    protein_b_id: str
    comparison_tool: str
    comparison_tool_version: str
    tm_score: float | None
    rmsd_angstrom: float | None
    aligned_residue_count: int | None
    coverage_a: float | None
    coverage_b: float | None
    comparison_status: StructureComparisonStatus
    source_record_id: str
    comparison_universe_id: str
    coverage_scope: StructureCoverageScope

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable comparison row.

        Returns:
            Dictionary representation of this comparison.
        """

        row = asdict(self)
        row["comparison_status"] = self.comparison_status.value
        row["coverage_scope"] = self.coverage_scope.value
        return row


@dataclass(frozen=True)
class FoldseekRunEvidence:
    """Completed or reused Foldseek all-versus-all evidence."""

    comparisons: tuple[PairwiseStructureComparison, ...]
    raw_output_path: Path
    completion_manifest_path: Path
    tool_version: str
    cache_key: str
    comparison_universe_id: str
    assessment_universe: frozenset[str]
    reused: bool


@dataclass(frozen=True)
class RedundancyClusterMembership:
    """One protein in an exact or externally supplied redundancy cluster."""

    protein_id: str
    cluster_id: str
    cluster_type: str
    method: str
    method_version: str
    identity_threshold: float
    coverage_threshold: float
    evidence_reference: str

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable redundancy membership.

        Returns:
            Dictionary representation of this membership.
        """

        return asdict(self)


@dataclass(frozen=True)
class StructuralAlignmentImport:
    """Evidence imported from a completed structural-alignment resource."""

    comparisons: tuple[PairwiseStructureComparison, ...]
    features: tuple[FeatureRecord, ...]
    group_summaries: tuple[dict[str, Any], ...]
    input_paths: tuple[Path, ...]
    package_version: str
    run_digest: str
    comparison_universe_members: Mapping[str, frozenset[str]]


@dataclass(frozen=True)
class ComparisonDefinition:
    """One explicit target-versus-background scientific comparison."""

    comparison_id: str
    display_name: str
    target_label_ids: tuple[str, ...]
    background_label_ids: tuple[str, ...]
    description: str


@dataclass(frozen=True)
class AnalysisSettings:
    """Validated settings for signature discovery and validation."""

    kmer_lengths: tuple[int, ...]
    minimum_target_proteins: int
    minimum_background_proteins: int
    minimum_feature_proteins: int
    maximum_kmer_features: int
    validation_fraction: float
    fdr_threshold: float
    random_seed: int
    structural_tm_score_threshold: float
    structural_minimum_coverage: float


@dataclass(frozen=True)
class AlphaFoldSettings:
    """Settings for optional AlphaFold Database model acquisition."""

    enabled: bool
    cache_dir: Path
    timeout_seconds: float
    retries: int
    minimum_mean_plddt: float


@dataclass(frozen=True)
class FoldseekSettings:
    """Settings for optional Foldseek all-versus-all structural analysis."""

    enabled: bool
    executable: str
    cache_dir: Path
    e_value_threshold: float
    sensitivity: float
    maximum_hits: int


@dataclass(frozen=True)
class ExplainableMLSettings:
    """Settings for mandatory group-aware explainable classification."""

    enabled: bool
    minimum_samples_per_class: int
    minimum_groups_per_class: int
    maximum_features: int
    regularisation_strengths: tuple[float, ...]
    l1_ratio: float
    cross_validation_folds: int
    maximum_iterations: int
    convergence_tolerance: float
    permutation_repeats: int
    maximum_permutation_features: int
    top_local_explanations: int
    maximum_waterfall_plots: int
    plot_cache_dir: Path
    exclude_technical_features: bool


@dataclass(frozen=True)
class InputPaths:
    """Resolved input paths for one campaign."""

    sequences_fasta: Path
    label_assignments: Path
    label_evidence_marker: Path | None
    label_evidence_audit: Path | None
    control_matching_audit: Path | None
    label_definition_features: Path | None
    class_labelling_summary: Path | None
    unresolved_assignments: Path | None
    features: Path | None
    domains: Path | None
    redundancy_clusters: Path | None
    structures: Path | None
    structure_comparisons: Path | None
    structural_alignment_resource: Path | None
    alphafold_accessions: Path | None
    orthofinder_resource: Path | None
    orthofinder_results: Path | None
    orthofinder_group_type: str
    orthofinder_hierarchy_node: str
    orthofinder_run_id: str


@dataclass(frozen=True)
class CampaignConfig:
    """Complete validated campaign configuration."""

    schema_version: int
    campaign_id: str
    profile_name: str
    config_path: Path
    inputs: InputPaths
    comparisons: tuple[ComparisonDefinition, ...]
    analysis: AnalysisSettings
    alphafold: AlphaFoldSettings
    foldseek: FoldseekSettings
    explainable_ml: ExplainableMLSettings


@dataclass(frozen=True)
class AlphaFoldRequest:
    """Mapping between a campaign protein and an AlphaFold DB accession."""

    protein_id: str
    uniprot_accession: str


@dataclass(frozen=True)
class AlphaFoldAcquisition:
    """Auditable outcome of one AlphaFold Database acquisition request."""

    protein_id: str
    uniprot_accession: str
    acquisition_status: str
    structure_id: str
    model_version: str
    api_url: str
    coordinate_path: Path | None
    coordinate_sha256: str
    sequence_match: bool | None
    mean_plddt: float | None
    message: str

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable acquisition row.

        Returns:
            Dictionary representation with a portable coordinate path.
        """

        row = asdict(self)
        row["coordinate_path"] = str(self.coordinate_path) if self.coordinate_path else ""
        return row


@dataclass(frozen=True)
class ProfileLabel:
    """One node in a protein-type profile hierarchy."""

    label_id: str
    display_name: str
    parent_label_id: str
    level: str
    system_class: str
    mechanistic_class: str
    component_role: str
    family: str
    active_site_expected: str
    active_site_residue: str
    default_analysis: bool
    default_background_label_id: str
    assignment_exclusivity_group: str
    reviewed_positive_allowed: bool
    aliases: tuple[str, ...]
    description: str


@dataclass(frozen=True)
class ProteinProfile:
    """A validated configurable hierarchy of protein types."""

    profile_id: str
    profile_version: str
    display_name: str
    require_structural_evidence: bool
    labels: tuple[ProfileLabel, ...]
    default_target_root_label_id: str
    default_background_label_id: str
    default_comparison_description: str
    default_excluded_label_ids: tuple[str, ...]
    default_excluded_subtree_label_ids: tuple[str, ...]

    def label_ids(self) -> frozenset[str]:
        """Return every canonical label identifier.

        Returns:
            Immutable set of label identifiers.
        """

        return frozenset(label.label_id for label in self.labels)


@dataclass(frozen=True)
class OrthoFinderLayout:
    """Relevant authorities discovered in completed OrthoFinder output."""

    results_dir: Path
    version: str
    major_version: int
    adapter_name: str
    primary_group_authority: str
    source_mode: str
    log_path: Path | None
    completion_authority_paths: tuple[Path, ...]
    orthogroups_path: Path | None
    hog_paths: tuple[Path, ...]
    species_ids_path: Path | None
    sequence_ids_path: Path | None

    def to_record(self) -> dict[str, Any]:
        """Return a JSON-compatible layout record.

        Returns:
            Dictionary with paths represented as strings.
        """

        return {
            "results_dir": str(self.results_dir),
            "version": self.version,
            "major_version": self.major_version,
            "adapter_name": self.adapter_name,
            "primary_group_authority": self.primary_group_authority,
            "log_path": str(self.log_path or ""),
            "completion_authority_paths": [str(path) for path in self.completion_authority_paths],
            "orthogroups_path": str(self.orthogroups_path or ""),
            "hog_paths": [str(path) for path in self.hog_paths],
            "species_ids_path": str(self.species_ids_path or ""),
            "sequence_ids_path": str(self.sequence_ids_path or ""),
            "source_mode": self.source_mode,
        }


@dataclass(frozen=True)
class OrthoFinderResource:
    """Validated published resource from the companion orthofinder-results package."""

    resource_dir: Path
    manifest_path: Path
    database_path: Path
    run_id: str
    schema_version: int
    package_version: str
    orthofinder_version: str
    adapter_name: str
    primary_group_authority: str
    relations: frozenset[str]
    focus_analysis_available: bool

    def to_record(self) -> dict[str, Any]:
        """Return a JSON-compatible resource description.

        Returns:
            Resource identity, compatibility and relation inventory.
        """

        return {
            "resource_dir": str(self.resource_dir),
            "manifest_path": str(self.manifest_path),
            "database_path": str(self.database_path),
            "run_id": self.run_id,
            "schema_version": self.schema_version,
            "package_version": self.package_version,
            "orthofinder_version": self.orthofinder_version,
            "adapter_name": self.adapter_name,
            "primary_group_authority": self.primary_group_authority,
            "relations": sorted(self.relations),
            "focus_analysis_available": self.focus_analysis_available,
            "source_mode": "ORTHOFINDER_RESULTS_RESOURCE",
        }


@dataclass(frozen=True)
class GroupMembership:
    """One protein membership in a completed OrthoFinder group."""

    run_id: str
    group_type: str
    hierarchy_node: str
    group_id: str
    legacy_orthogroup_id: str
    gene_tree_parent_clade: str
    species_label: str
    protein_id: str

    def to_record(self) -> dict[str, str]:
        """Return a serialisable group-membership row.

        Returns:
            Dictionary representation of this membership.
        """

        return asdict(self)


@dataclass(frozen=True)
class PartitionAssignment:
    """One protein assigned to discovery or held-out validation."""

    protein_id: str
    partition: str
    partition_unit: str
    partition_key: str

    def to_record(self) -> dict[str, str]:
        """Return a serialisable partition row.

        Returns:
            Dictionary representation of this assignment.
        """

        return asdict(self)


@dataclass(frozen=True)
class AssociationResult:
    """Feature enrichment result for one explicit comparison and partition."""

    comparison_id: str
    partition: str
    analysis_unit: str
    feature_type: str
    feature_id: str
    feature_name: str
    target_protein_count: int
    background_protein_count: int
    target_with_feature: int
    background_with_feature: int
    target_unit_count: int
    background_unit_count: int
    target_assessed_unit_count: int
    background_assessed_unit_count: int
    target_unknown_unit_count: int
    background_unknown_unit_count: int
    target_units_with_feature: int
    background_units_with_feature: int
    excluded_mixed_unit_count: int
    target_prevalence: float | None
    background_prevalence: float | None
    target_prevalence_ci_lower: float | None
    target_prevalence_ci_upper: float | None
    background_prevalence_ci_lower: float | None
    background_prevalence_ci_upper: float | None
    prevalence_ci_method: str
    prevalence_difference: float | None
    odds_ratio: float | None
    p_value: float | None
    q_value: float | None
    study_q_value: float | None
    direction: str
    status: AnalysisStatus

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable association row.

        Returns:
            Dictionary representation with controlled status text.
        """

        row = asdict(self)
        row["status"] = self.status.value
        return row


@dataclass(frozen=True)
class SignatureSummary:
    """Discovery and held-out evidence for one candidate signature."""

    comparison_id: str
    feature_type: str
    feature_id: str
    feature_name: str
    discovery_q_value: float | None
    discovery_study_q_value: float | None
    discovery_prevalence_difference: float | None
    validation_q_value: float | None
    validation_study_q_value: float | None
    validation_prevalence_difference: float | None
    evidence_class: str
    status: AnalysisStatus

    def to_record(self) -> dict[str, Any]:
        """Return a serialisable signature-summary row.

        Returns:
            Dictionary representation with a controlled status value.
        """

        row = asdict(self)
        row["status"] = self.status.value
        return row
