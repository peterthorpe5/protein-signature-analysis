"""Canonical Arrow schemas for every published result table."""

from __future__ import annotations

import pyarrow as pa

from .errors import PublicationError


def table_schemas() -> dict[str, pa.Schema]:
    """Return typed schemas keyed by stable result-table name.

    Returns:
        New mapping of table names to immutable Arrow schemas.
    """

    text = pa.string()
    integer = pa.int64()
    number = pa.float64()
    boolean = pa.bool_()
    return {
        "proteins": pa.schema(
            [
                ("protein_id", text),
                ("description", text),
                ("sequence", text),
                ("sequence_length", integer),
                ("sequence_sha256", text),
            ]
        ),
        "redundancy_clusters": pa.schema(
            [
                ("protein_id", text),
                ("cluster_id", text),
                ("cluster_type", text),
                ("method", text),
                ("method_version", text),
                ("identity_threshold", number),
                ("coverage_threshold", number),
                ("evidence_reference", text),
            ]
        ),
        "profile_labels": pa.schema(
            [
                ("label_id", text),
                ("display_name", text),
                ("parent_label_id", text),
                ("level", text),
                ("system_class", text),
                ("mechanistic_class", text),
                ("component_role", text),
                ("family", text),
                ("active_site_expected", text),
                ("active_site_residue", text),
                ("default_analysis", boolean),
                ("default_background_label_id", text),
                ("assignment_exclusivity_group", text),
                ("reviewed_positive_allowed", boolean),
                ("aliases", text),
                ("description", text),
            ]
        ),
        "label_assignments": pa.schema(
            [
                ("protein_id", text),
                ("label_id", text),
                ("curation_status", text),
                ("evidence_status", text),
                ("evidence_source", text),
                ("evidence_reference", text),
                ("component_role", text),
                ("curation_reason", text),
            ]
        ),
        "label_memberships": pa.schema(
            [
                ("protein_id", text),
                ("label_id", text),
                ("membership_source", text),
                ("direct_label_id", text),
            ]
        ),
        "comparisons": pa.schema(
            [
                ("comparison_id", text),
                ("display_name", text),
                ("target_label_ids", text),
                ("background_label_ids", text),
                ("description", text),
            ]
        ),
        "features": pa.schema(
            [
                ("protein_id", text),
                ("feature_type", text),
                ("feature_id", text),
                ("feature_name", text),
                ("start", integer),
                ("end", integer),
                ("evidence_status", text),
                ("evidence_source", text),
                ("evidence_reference", text),
                ("derivation_scope", text),
                ("feature_definition_sha256", text),
                ("derivation_cohort_sha256", text),
            ]
        ),
        "feature_assessments": pa.schema(
            [
                ("protein_id", text),
                ("feature_type", text),
                ("feature_id", text),
                ("feature_name", text),
                ("start", integer),
                ("end", integer),
                ("evidence_status", text),
                ("evidence_source", text),
                ("evidence_reference", text),
                ("derivation_scope", text),
                ("feature_definition_sha256", text),
                ("derivation_cohort_sha256", text),
            ]
        ),
        "domain_hits": pa.schema(
            [
                ("protein_id", text),
                ("domain_authority", text),
                ("domain_id", text),
                ("domain_name", text),
                ("start", integer),
                ("end", integer),
                ("score", number),
                ("e_value", number),
                ("evidence_source", text),
                ("evidence_reference", text),
            ]
        ),
        "domain_assessments": pa.schema(
            [
                ("protein_id", text),
                ("domain_authority", text),
                ("assessment_status", text),
                ("hit_count", integer),
                ("evidence_source", text),
                ("evidence_reference", text),
            ]
        ),
        "domain_sequences": pa.schema(
            [
                ("protein_id", text),
                ("domain_authority", text),
                ("domain_id", text),
                ("start", integer),
                ("end", integer),
                ("domain_sequence", text),
            ]
        ),
        "structures": pa.schema(
            [
                ("protein_id", text),
                ("structure_id", text),
                ("structure_source", text),
                ("structure_version", text),
                ("coordinate_path", text),
                ("coordinate_sha256", text),
                ("availability_status", text),
                ("mean_confidence", number),
                ("fold_id", text),
                ("fold_name", text),
                ("fold_authority", text),
                ("fold_authority_version", text),
                ("fold_evidence_reference", text),
                ("fold_evidence_status", text),
                ("analysis_eligibility_status", text),
                ("comparison_universe_ids", text),
            ]
        ),
        "alphafold_acquisitions": pa.schema(
            [
                ("protein_id", text),
                ("uniprot_accession", text),
                ("acquisition_status", text),
                ("structure_id", text),
                ("model_version", text),
                ("api_url", text),
                ("coordinate_path", text),
                ("coordinate_sha256", text),
                ("sequence_match", boolean),
                ("mean_plddt", number),
                ("message", text),
            ]
        ),
        "structure_comparisons": pa.schema(
            [
                ("protein_a_id", text),
                ("protein_b_id", text),
                ("comparison_tool", text),
                ("comparison_tool_version", text),
                ("tm_score", number),
                ("rmsd_angstrom", number),
                ("aligned_residue_count", integer),
                ("coverage_a", number),
                ("coverage_b", number),
                ("comparison_status", text),
                ("source_record_id", text),
                ("comparison_universe_id", text),
                ("coverage_scope", text),
            ]
        ),
        "structure_clusters": pa.schema(
            [
                ("cluster_id", text),
                ("protein_id", text),
                ("member_count", integer),
                ("reference_member_count", integer),
                ("edge_count", integer),
                ("reference_partition", text),
                ("membership_method", text),
                ("supporting_edge_count", integer),
                ("best_tm_score", number),
                ("tm_score_threshold", number),
                ("minimum_coverage", number),
                ("clustering_method", text),
                ("comparison_universe_id", text),
                ("coverage_scope", text),
                ("comparison_tool", text),
                ("comparison_tool_version", text),
            ]
        ),
        "imported_structural_group_summaries": pa.schema(
            [
                ("cluster_id", text),
                ("primary_group_type", text),
                ("primary_group_id", text),
                ("reference_accession", text),
                ("alignment_tools", text),
                ("alignment_tool_count", integer),
                ("selected_accession_count", integer),
                ("model_available_accession_count", integer),
                ("aligned_accession_count", integer),
                ("supported_accession_count", integer),
                ("position_supported_accession_count", integer),
                ("group_support_fraction", number),
                ("group_position_support_fraction", number),
                ("mean_minimum_tm_score", number),
                ("mean_pocket_overlap_fraction", number),
                ("median_centroid_distance_angstrom", number),
                ("position_alignment_status", text),
                ("alignment_status", text),
                ("interpretation", text),
            ]
        ),
        "orthofinder_memberships": pa.schema(
            [
                ("run_id", text),
                ("group_type", text),
                ("hierarchy_node", text),
                ("group_id", text),
                ("legacy_orthogroup_id", text),
                ("gene_tree_parent_clade", text),
                ("species_label", text),
                ("protein_id", text),
            ]
        ),
        "orthofinder_group_context": pa.schema(
            [
                ("run_id", text),
                ("group_type", text),
                ("hierarchy_node", text),
                ("group_id", text),
                ("legacy_orthogroup_id", text),
                ("gene_tree_parent_clade", text),
                ("member_count", integer),
                ("species_count", integer),
                ("single_copy_species_count", integer),
                ("max_copies_per_species", integer),
                ("mean_copies_per_species", number),
                ("is_singleton", boolean),
                ("distance_method", text),
                ("computation_status", text),
                ("total_member_count", integer),
                ("sampled_member_count", integer),
                ("distance_pair_count", integer),
                ("unresolved_pair_count", integer),
                ("minimum_distance", number),
                ("q25_distance", number),
                ("median_distance", number),
                ("mean_distance", number),
                ("q75_distance", number),
                ("maximum_distance", number),
                ("population_stddev_distance", number),
                ("failure_reason", text),
            ]
        ),
        "partitions": pa.schema(
            [
                ("protein_id", text),
                ("partition", text),
                ("partition_unit", text),
                ("partition_key", text),
            ]
        ),
        "associations": pa.schema(
            [
                ("comparison_id", text),
                ("partition", text),
                ("analysis_unit", text),
                ("feature_type", text),
                ("feature_id", text),
                ("feature_name", text),
                ("target_protein_count", integer),
                ("background_protein_count", integer),
                ("target_with_feature", integer),
                ("background_with_feature", integer),
                ("target_unit_count", integer),
                ("background_unit_count", integer),
                ("target_assessed_unit_count", integer),
                ("background_assessed_unit_count", integer),
                ("target_unknown_unit_count", integer),
                ("background_unknown_unit_count", integer),
                ("target_units_with_feature", integer),
                ("background_units_with_feature", integer),
                ("excluded_mixed_unit_count", integer),
                ("target_prevalence", number),
                ("background_prevalence", number),
                ("target_prevalence_ci_lower", number),
                ("target_prevalence_ci_upper", number),
                ("background_prevalence_ci_lower", number),
                ("background_prevalence_ci_upper", number),
                ("prevalence_ci_method", text),
                ("prevalence_difference", number),
                ("odds_ratio", number),
                ("p_value", number),
                ("q_value", number),
                ("study_q_value", number),
                ("direction", text),
                ("status", text),
            ]
        ),
        "signatures": pa.schema(
            [
                ("comparison_id", text),
                ("feature_type", text),
                ("feature_id", text),
                ("feature_name", text),
                ("discovery_q_value", number),
                ("discovery_study_q_value", number),
                ("discovery_prevalence_difference", number),
                ("validation_q_value", number),
                ("validation_study_q_value", number),
                ("validation_prevalence_difference", number),
                ("evidence_class", text),
                ("status", text),
            ]
        ),
        "ml_models": pa.schema(
            [
                ("comparison_id", text),
                ("model_type", text),
                ("status", text),
                ("status_message", text),
                ("analysis_unit", text),
                ("discovery_target_count", integer),
                ("discovery_background_count", integer),
                ("discovery_group_count", integer),
                ("validation_target_count", integer),
                ("validation_background_count", integer),
                ("validation_group_count", integer),
                ("feature_count", integer),
                ("excluded_technical_feature_rows", integer),
                ("cross_validation_folds", integer),
                ("regularisation_strength", number),
                ("l1_ratio", number),
                ("cv_roc_auc", number),
                ("validation_roc_auc", number),
                ("validation_average_precision", number),
                ("validation_balanced_accuracy", number),
                ("validation_matthews_correlation", number),
                ("validation_brier_score", number),
                ("decision_threshold", number),
                ("explanation_method", text),
                ("shap_background_partition", text),
                ("shap_explained_partition", text),
                ("shap_explained_sample_count", integer),
                ("shap_expected_value_log_odds", number),
            ]
        ),
        "ml_feature_importance": pa.schema(
            [
                ("comparison_id", text),
                ("feature_type", text),
                ("feature_id", text),
                ("feature_name", text),
                ("coefficient_log_odds", number),
                ("odds_multiplier", number),
                ("discovery_target_prevalence", number),
                ("discovery_background_prevalence", number),
                ("discovery_prevalence_difference", number),
                ("validation_permutation_importance_mean", number),
                ("validation_permutation_importance_stddev", number),
                ("mean_absolute_validation_contribution", number),
                ("importance_rank", integer),
            ]
        ),
        "ml_predictions": pa.schema(
            [
                ("comparison_id", text),
                ("protein_id", text),
                ("partition", text),
                ("partition_key", text),
                ("true_class", text),
                ("predicted_probability", number),
                ("predicted_class", text),
                ("correct", boolean),
            ]
        ),
        "ml_explanations": pa.schema(
            [
                ("comparison_id", text),
                ("protein_id", text),
                ("partition", text),
                ("true_class", text),
                ("predicted_probability", number),
                ("feature_type", text),
                ("feature_id", text),
                ("feature_name", text),
                ("feature_value", integer),
                ("discovery_background_mean", number),
                ("coefficient_log_odds", number),
                ("log_odds_contribution", number),
                ("absolute_rank", integer),
                ("explanation_method", text),
            ]
        ),
        "ml_plot_inventory": pa.schema(
            [
                ("comparison_id", text),
                ("plot_type", text),
                ("partition", text),
                ("protein_id", text),
                ("file_format", text),
                ("asset_path", text),
                ("sample_count", integer),
                ("feature_count", integer),
                ("explanation_method", text),
            ]
        ),
    }


def schema_for(*, table_name: str) -> pa.Schema:
    """Resolve one canonical result-table schema.

    Args:
        table_name: Stable result table name.

    Returns:
        Arrow schema.

    Raises:
        PublicationError: If the table name is unknown.
    """

    schemas = table_schemas()
    if table_name not in schemas:
        raise PublicationError(f"No canonical schema is defined for table {table_name!r}.")
    return schemas[table_name]
