"""Shared scientific column definitions for result tables and Excel exports."""

from __future__ import annotations

import logging

from .errors import InputValidationError

LOGGER = logging.getLogger(__name__)

COLUMN_DEFINITIONS = {
    "plot_negative_log10_q_value": (
        "Display-only −log10 discovery q-value, capped at 300 for recorded zeros or values "
        "below 1e-300; the published q-value is retained unchanged."
    ),
    "plot_q_value_clipped": (
        "True when a recorded zero or q-value below 1e-300 reaches the logarithmic plot's "
        "display cap. This is a display flag, not a scientific outcome."
    ),
    "plot_feature_label": (
        "Readable plot label combining a shortened feature name with its recorded stable ID."
    ),
    "protein_id": "Exact FASTA protein identifier; no isoform or accession rewriting is implied.",
    "description": "Description supplied by the relevant input or profile authority.",
    "sequence": "Authoritative amino-acid sequence in one-letter codes; positions start at 1.",
    "sequence_length": "Number of amino-acid residues in the authoritative sequence.",
    "sequence_sha256": "SHA-256 content checksum of the sequence, not a biological score.",
    "start": "First residue of an interval, numbered from 1 and included in the interval.",
    "end": "Last residue of an interval, included; interval length is end minus start plus 1.",
    "comparison_id": "Stable identifier for one declared target-versus-background comparison.",
    "display_name": "Readable name from the class profile or comparison definition.",
    "feature_type": "Evidence family, such as exact k-mers, Pfam domains or structural clusters.",
    "feature_id": "Stable key within an evidence family; read feature_name and its explanation.",
    "feature_name": "Published readable feature name; an opaque cluster key is not a fold name.",
    "feature_explanation": "Biological meaning and interpretation limits of the feature key.",
    "feature_value": "Binary model input: 1 means the feature is present and 0 means absent.",
    "feature_count": "Number of feature definitions in the stated table, model or plot.",
    "feature_definition_sha256": "Checksum binding a feature to its frozen definition.",
    "derivation_cohort_sha256": "Checksum of the protein cohort used to define the feature.",
    "derivation_scope": (
        "Feature provenance: FIXED_EXTERNAL, label-blind DISCOVERY_DERIVED, "
        "audit-only DISCOVERY_SUPERVISED or ALL_DATA_EXPLORATORY."
    ),
    "evidence_source": (
        "Named producer or input authority; importing evidence is not computing it."
    ),
    "evidence_reference": (
        "Traceable release, record, path or checksum reference for the evidence."
    ),
    "evidence_status": "Producer-specific evidence state; missing or unassessed is not a no-hit.",
    "evidence_class": "Discovery/validation decision tier, including local versus study-wide FDR.",
    "status": "Outcome of the stated operation. COMPLETE means software completion within scope.",
    "status_message": "Published explanation of a model status, including exclusions and counts.",
    "analysis_status": "Summary of a comparison outcome; inspect its underlying status rows.",
    "reason": "Recorded reason for a decision, exclusion, failure or unavailable assessment.",
    "message": "Producer's diagnostic message; retain it when investigating a status.",
    "failure_reason": "Recorded reason a calculation failed or could not be assessed.",
    "interpretation": "Producer's interpretation and scientific limitations of this record.",
    "partition": "DISCOVERY defines/selects features; VALIDATION tests frozen candidates.",
    "partition_key": "Connected independence-block key; all linked proteins stay in one split.",
    "partition_unit": "Authority forming a block, such as exact sequence or OrthoFinder links.",
    "analysis_unit": (
        "Inferential unit. Association and held-out metrics use pure independent blocks."
    ),
    "target_protein_count": (
        "Raw target proteins in this partition; not the Fisher-test denominator."
    ),
    "background_protein_count": "Raw background proteins; not the Fisher-test denominator.",
    "target_with_feature": (
        "Raw target proteins with a positive feature; may include excluded blocks."
    ),
    "background_with_feature": "Raw background proteins with the feature; not block prevalence.",
    "target_unit_count": "Pure target independence blocks after mixed-class blocks are excluded.",
    "background_unit_count": "Pure background blocks after mixed-class blocks are excluded.",
    "target_assessed_unit_count": (
        "Target blocks with known feature presence/absence; test denominator."
    ),
    "background_assessed_unit_count": (
        "Background blocks with known presence/absence; test denominator."
    ),
    "target_unknown_unit_count": (
        "Pure target blocks with unknown feature absence; excluded from testing."
    ),
    "background_unknown_unit_count": (
        "Background blocks with unknown absence; excluded from testing."
    ),
    "target_units_with_feature": (
        "Assessed target blocks with at least one relevant positive protein."
    ),
    "background_units_with_feature": (
        "Assessed background blocks with at least one positive protein."
    ),
    "excluded_mixed_unit_count": (
        "Blocks containing both classes, excluded rather than counted twice."
    ),
    "target_prevalence": "Fraction of assessed pure target blocks with the feature, from 0 to 1.",
    "background_prevalence": (
        "Fraction of assessed pure background blocks with the feature, 0 to 1."
    ),
    "target_prevalence_ci_lower": (
        "Lower limit of the target prevalence's 95% Wilson score interval."
    ),
    "target_prevalence_ci_upper": (
        "Upper limit of the target prevalence's 95% Wilson score interval."
    ),
    "background_prevalence_ci_lower": (
        "Lower limit of the background prevalence's 95% Wilson interval."
    ),
    "background_prevalence_ci_upper": (
        "Upper limit of the background prevalence's 95% Wilson interval."
    ),
    "prevalence_ci_method": (
        "Interval method and denominator authority; not an interval for the difference."
    ),
    "prevalence_difference": (
        "Target minus background prevalence, -1 to 1; +0.20 is 20 percentage points."
    ),
    "odds_ratio": "Feature-presence odds in target divided by background odds; 1 is equal odds.",
    "p_value": "Two-sided Fisher exact p-value on assessed pure blocks; not probability of truth.",
    "q_value": (
        "Benjamini-Hochberg adjusted p-value within comparison, partition and feature family."
    ),
    "study_q_value": "BH adjustment across comparisons in the same partition and feature family.",
    "direction": "Sign of target-minus-background prevalence: enrichment or depletion in target.",
    "label_id": "Stable controlled protein-class label in the versioned profile.",
    "label_type": "TARGET or BACKGROUND in an automated labelling summary.",
    "parent_label_id": "Immediate parent in the class hierarchy; inherited counts can overlap.",
    "level": "Profile hierarchy level; it is not a confidence or activity score.",
    "system_class": (
        "Biological system classification, separate from mechanism and component role."
    ),
    "mechanistic_class": "Profile-defined biochemical mechanism; not inferred from a SHAP value.",
    "component_role": (
        "Component function, for example catalytic protein, scaffold or substrate receptor."
    ),
    "family": "Profile-defined protein family within its class hierarchy.",
    "active_site_expected": "Profile expectation of an active site; not a measured residue call.",
    "active_site_residue": "Profile-described expected residue chemistry, not observed activity.",
    "default_analysis": "Whether the profile requests a default comparison for this label.",
    "default_background_label_id": "Profile's default comparison-background label.",
    "assignment_exclusivity_group": (
        "Labels that cannot be assigned together under the profile policy."
    ),
    "reviewed_positive_allowed": (
        "Whether the profile permits a human-reviewed positive for this label."
    ),
    "aliases": "Alternative names declared by the profile; not additional independent evidence.",
    "curation_status": (
        "Review state, including reviewed, provisional, proposed, ambiguous or unmapped."
    ),
    "curation_reason": "Authority's recorded reason for the label's review or proposal state.",
    "direct_label_id": "Original direct label from which this membership was obtained.",
    "membership_source": "DIRECT label membership or INHERITED membership through a parent label.",
    "target_label_ids": "Exact target labels for this comparison or pooled matching stratum.",
    "background_label_ids": "Exact background labels defining this comparison.",
    "background_label_id": "Named pooled background used by the control-matching audit.",
    "matched_background_label_id": "Background pool allocated for this automated target label.",
    "direct_positive_protein_count": (
        "Accepted directly labelled proteins, before comparison filtering."
    ),
    "independent_unit_count": (
        "Distinct independence units in this labelled cohort or pooled summary."
    ),
    "matched_control_protein_count": (
        "Selected control representatives; a pooled count is not a denominator."
    ),
    "matched_control_unit_count": (
        "Distinct selected control units in the named pooled background."
    ),
    "independence_unit": "Connected orthology, exact-sequence and supplied redundancy unit.",
    "rule_id": (
        "Versioned evidence rule responsible for a label proposal or defining-feature exclusion."
    ),
    "decision": "Automated rule decision such as ACCEPTED, PROPOSED or AMBIGUOUS.",
    "confidence_tier": "Rule-support category, not a calibrated probability of protein function.",
    "score": (
        "Producer-specific evidence or domain score; not a p-value or calibrated probability."
    ),
    "best_score": "Highest rule-support score among the proposed labels for this protein.",
    "evidence_group_count": "Independent evidence categories, not the number of database rows.",
    "evidence_groups": "Categories of supporting evidence counted by the labelling rule.",
    "evidence_items": "Detailed source records supporting or contradicting the label decision.",
    "conflicting_label_ids": "Incompatible candidate labels retained for review.",
    "candidate_label_ids": "Labels considered for an unresolved protein; not accepted positives.",
    "provisional_label_id": "Best proposed label in the unresolved review queue.",
    "evidence_role": (
        "Role of a feature in defining a label; defining domains are excluded downstream."
    ),
    "exclusion_scope": "Where label-defining features were removed to reduce direct circularity.",
    "target_protein_id": "Representative protein of a target matching unit.",
    "target_unit_id": "Independence unit for a target control request.",
    "control_protein_id": (
        "Selected representative of a control unit; blank for unmatched requests."
    ),
    "control_unit_id": "Selected control independence unit; not reused within the named pool.",
    "species_match": "Whether target and control share a species label under the matching rule.",
    "structure_eligibility_match": (
        "Whether target/control have matching structural eligibility states."
    ),
    "length_log2_difference": "Absolute log2 target/control length difference; smaller is closer.",
    "domain_count_difference": (
        "Absolute difference in annotated domain counts between representatives."
    ),
    "domain_architecture_jaccard": (
        "Domain-set intersection divided by union, 0 to 1; higher is closer."
    ),
    "mean_confidence_difference": (
        "Absolute difference in representatives' reported mean confidence."
    ),
    "match_score": (
        "Deterministic covariate distance for caliper-passing controls; smaller is closer."
    ),
    "domain_authority": "Domain catalogue or scan authority, usually a pinned Pfam release.",
    "domain_id": "Stable domain accession from its authority, such as a Pfam PF identifier.",
    "domain_name": "Readable domain name from the source authority.",
    "domain_sequence": "Exact authoritative subsequence from start through end, inclusive.",
    "assessment_status": (
        "Hit, assessed no-hit, not assessed or failed; only assessed states define absence."
    ),
    "hit_count": "Number of retained domain hits for the protein and authority.",
    "e_value": (
        "Expected chance-hit count under the source tool's model; not an association q-value."
    ),
    "structure_id": "Stable source/model key; multiple models can exist for one protein.",
    "structure_source": (
        "Authority that produced or supplied coordinates; may contain import provenance."
    ),
    "structure_version": "Version of the coordinate model from its source authority.",
    "coordinate_path": "Coordinate asset path; keep it with the whole portable completed result.",
    "coordinate_sha256": "SHA-256 checksum of the packaged coordinate asset.",
    "availability_status": (
        "Whether coordinates were available; unavailable is not a structural negative."
    ),
    "mean_confidence": (
        "Source-reported average model confidence; AlphaFold pLDDT runs from 0 to 100."
    ),
    "analysis_eligibility_status": (
        "Recorded structural eligibility, including low-confidence exclusion."
    ),
    "fold_id": (
        "Named fold identifier from an explicit authority; blank means no named assignment."
    ),
    "fold_name": "Readable named fold from its authority; distinct from an alignment-cluster key.",
    "fold_authority": (
        "Authority assigning named folds; structural similarity alone supplies no fold name."
    ),
    "fold_authority_version": "Release/version of the named-fold authority.",
    "fold_evidence_reference": "Source record or release supporting the fold assignment.",
    "fold_evidence_status": "Explicit fold-assessment state; NOT_ASSESSED does not mean no fold.",
    "comparison_universe_ids": (
        "Declared structural-search universes used to establish assessed absence."
    ),
    "comparison_universe_id": (
        "Structural-search universe key; retained hits alone are not a denominator."
    ),
    "uniprot_accession": (
        "Exact canonical UniProt accession used for controlled model acquisition."
    ),
    "acquisition_status": (
        "Terminal model-fetch outcome, including unavailable or sequence mismatch."
    ),
    "model_version": "Coordinate model version reported by the acquisition source.",
    "api_url": "Recorded model-source endpoint used for acquisition.",
    "sequence_match": (
        "Whether the downloaded model sequence exactly matched the published protein."
    ),
    "mean_plddt": "Mean AlphaFold per-residue confidence, 0 to 100; not a measure of enrichment.",
    "protein_a_id": "Exact identifier of the first protein in a pairwise structural comparison.",
    "protein_b_id": "Exact identifier of the second protein in a pairwise structural comparison.",
    "comparison_tool": "Actual structural comparison tool supplying the score.",
    "comparison_tool_version": "Recorded version of the structural comparison tool.",
    "tm_score": (
        "Normalised structural similarity, typically 0 to 1; inspect coverage and score scope."
    ),
    "rmsd_angstrom": (
        "Root-mean-square deviation of aligned coordinates in angstroms; lower is closer."
    ),
    "aligned_residue_count": "Number of aligned residues reported by the comparison producer.",
    "coverage_a": "Aligned fraction of protein A, 0 to 1, under the declared coverage scope.",
    "coverage_b": "Aligned fraction of protein B, 0 to 1, under the declared coverage scope.",
    "comparison_status": "Recorded pairwise comparison outcome; failed is not a no-hit.",
    "source_record_id": "Opaque upstream comparison record key retained for traceability.",
    "coverage_scope": (
        "What the search covered, for example a complete search universe or retained hits."
    ),
    "cluster_id": (
        "Stable group key; an SC_ key denotes a whole-model similarity component, not a motif."
    ),
    "cluster_type": "Exact-sequence or supplied near-redundancy grouping method.",
    "method": "Named source method for this derived or imported record.",
    "method_version": "Recorded implementation or authority version for the method.",
    "identity_threshold": "Minimum sequence identity used by the declared redundancy method.",
    "coverage_threshold": (
        "Coverage cut-off of the source redundancy method, not enrichment significance."
    ),
    "member_count": (
        "Number of members in the stated group; structural components include projections."
    ),
    "reference_member_count": (
        "Number of members defining the frozen discovery structural component."
    ),
    "reference_partition": "Partition used to define the reference component; normally DISCOVERY.",
    "edge_count": "Qualifying discovery structural edges; not the number of independent proteins.",
    "membership_method": (
        "DISCOVERY_COMPONENT, VALIDATION_PROJECTION or exploratory all-data component."
    ),
    "supporting_edge_count": "Qualifying discovery links supporting this structural membership.",
    "best_tm_score": "Strongest supporting structural similarity for this cluster membership.",
    "tm_score_threshold": "Minimum TM-score used to form/project the frozen structural component.",
    "minimum_coverage": (
        "Minimum aligned fraction required for both proteins, or plotted bilateral minimum."
    ),
    "clustering_method": (
        "Component construction method; chain connectivity is not all-pairs similarity."
    ),
    "primary_group_type": (
        "Type of the upstream group summarised by the imported structural resource."
    ),
    "primary_group_id": "Opaque upstream group identifier; not a newly discovered motif.",
    "reference_accession": "Upstream accession used as the structural group reference.",
    "alignment_tools": "Named upstream tools contributing to this imported summary.",
    "alignment_tool_count": "Number of distinct alignment tools in the upstream summary.",
    "selected_accession_count": "Accessions selected for the upstream structural group analysis.",
    "model_available_accession_count": (
        "Selected upstream accessions with an available coordinate model."
    ),
    "aligned_accession_count": "Upstream accessions with a completed alignment.",
    "supported_accession_count": (
        "Upstream accessions satisfying the declared structural support rule."
    ),
    "position_supported_accession_count": (
        "Accessions with upstream position support; not local motif q-values."
    ),
    "group_support_fraction": (
        "Upstream structural-support fraction; inspect its source denominator."
    ),
    "group_position_support_fraction": (
        "Upstream positional-support fraction; not campaign residue enrichment."
    ),
    "mean_minimum_tm_score": (
        "Upstream mean of the minimum required TM-score across its alignments."
    ),
    "mean_pocket_overlap_fraction": (
        "Imported pocket-overlap summary; no residue-level discovery is implied."
    ),
    "median_centroid_distance_angstrom": (
        "Median upstream coordinate-centroid distance in angstroms."
    ),
    "position_alignment_status": "Upstream positional-alignment assessment state.",
    "alignment_status": "Upstream overall structural-alignment completion state.",
    "run_id": "Identifier of the exact OrthoFinder run authority.",
    "group_type": "HOG or legacy orthogroup grouping; not interchangeable.",
    "hierarchy_node": "Exact HOG hierarchy node, often N0; absent for legacy orthogroups.",
    "group_id": "Group key together with run_id, group_type and hierarchy_node.",
    "legacy_orthogroup_id": (
        "Legacy orthogroup context; not a replacement for the chosen HOG identity."
    ),
    "gene_tree_parent_clade": "Upstream gene-tree parent-clade context.",
    "species_label": "Species label from the OrthoFinder input authority.",
    "species_count": "Number of represented species in the stated group.",
    "single_copy_species_count": "Represented species having exactly one member in this group.",
    "max_copies_per_species": "Largest per-species member count in the group.",
    "mean_copies_per_species": "Mean copies across the group's represented species.",
    "is_singleton": "Whether the group has one member; not whether every species is single-copy.",
    "distance_method": (
        "Declared upstream sequence-distance method; interpret values in that method."
    ),
    "computation_status": "Completion/availability state of the imported distance summary.",
    "total_member_count": "All members in the source group, before any upstream subsampling.",
    "sampled_member_count": "Source members actually used for the imported distance calculation.",
    "distance_pair_count": "Resolved pairwise distances contributing to the summary.",
    "unresolved_pair_count": "Requested distance pairs with unavailable or failed values.",
    "minimum_distance": "Smallest resolved pairwise distance under the source method.",
    "q25_distance": "25th percentile of resolved source pairwise distances.",
    "median_distance": "50th percentile of resolved source pairwise distances.",
    "mean_distance": "Arithmetic mean of resolved source pairwise distances.",
    "q75_distance": "75th percentile of resolved source pairwise distances.",
    "maximum_distance": "Largest resolved pairwise distance under the source method.",
    "population_stddev_distance": "Population standard deviation of resolved pairwise distances.",
    "model_type": "Published classifier family, usually elastic-net logistic regression.",
    "discovery_target_count": (
        "Target protein rows used by the model's discovery cohort after exclusions."
    ),
    "discovery_background_count": "Background protein rows in the model's discovery cohort.",
    "discovery_group_count": "Distinct discovery independence groups across both model classes.",
    "validation_target_count": (
        "Raw held-out target protein rows; pure-block evaluation may exclude more."
    ),
    "validation_background_count": "Raw held-out background proteins, not the metric denominator.",
    "validation_group_count": (
        "All held-out groups; inspect status_message for evaluated pure blocks."
    ),
    "excluded_technical_feature_rows": (
        "Technical availability/quality rows removed from model inputs."
    ),
    "cross_validation_folds": "Completed group-aware discovery CV folds; each needs both classes.",
    "regularisation_strength": (
        "Selected penalty strength lambda; logistic regression uses C = 1/lambda. "
        "Larger recorded values impose stronger regularisation."
    ),
    "l1_ratio": "Elastic-net L1 fraction, 0 to 1; the remainder is L2 regularisation.",
    "cv_roc_auc": (
        "Mean group-aware discovery CV ROC AUC used to select C; not final held-out evidence."
    ),
    "validation_roc_auc": "Held-out pure-block ROC AUC: 0.5 is chance ranking, 1 perfect ranking.",
    "validation_average_precision": (
        "Held-out average precision (AP), not trapezoidal PR AUC; compare class prevalence."
    ),
    "validation_balanced_accuracy": (
        "Mean sensitivity and specificity at threshold 0.5; 0.5 is chance, 1 perfect."
    ),
    "validation_matthews_correlation": (
        "MCC at threshold 0.5: -1 reversed, 0 no correlation, +1 perfect."
    ),
    "validation_brier_score": (
        "Mean squared probability error, 0 to 1; 0 is perfect and lower is better."
    ),
    "decision_threshold": "Probability cut-off for a TARGET call; this pipeline uses 0.5.",
    "explanation_method": (
        "Published attribution method; SHAP is model explanation, not causation."
    ),
    "shap_background_partition": (
        "Partition supplying the SHAP reference distribution, normally DISCOVERY."
    ),
    "shap_explained_partition": (
        "Partition supplying explained proteins; DISCOVERY is not held-out evidence."
    ),
    "shap_explained_sample_count": "Number of protein predictions receiving SHAP explanations.",
    "shap_expected_value_log_odds": (
        "SHAP reference prediction in log-odds, before feature contributions."
    ),
    "coefficient_log_odds": (
        "Signed logistic coefficient; positive favours TARGET, conditional on other inputs."
    ),
    "odds_multiplier": (
        "Exponential of the coefficient; multiplicative fitted odds effect, not causation."
    ),
    "discovery_background_mean": (
        "Mean binary feature presence across the full pure discovery-block SHAP reference."
    ),
    "validation_permutation_importance_mean": (
        "Mean held-out pure-block balanced-accuracy drop after shuffling this feature "
        "at threshold 0.5; negative can occur."
    ),
    "validation_permutation_importance_stddev": (
        "Population standard deviation of repeated balanced-accuracy drops; "
        "not a confidence interval."
    ),
    "mean_absolute_validation_contribution": (
        "Mean absolute linear-model SHAP log-odds contribution across validation "
        "proteins, using the full pure discovery-block reference; missing without validation."
    ),
    "importance_rank": (
        "Rank by absolute fitted coefficient, with feature-key tie-breaking. "
        "This is not the ordering of the mean-absolute-SHAP bar chart."
    ),
    "true_class": (
        "Published input class (TARGET/BACKGROUND); provisional labels are not experimental truth."
    ),
    "predicted_probability": (
        "Fitted TARGET probability conditional on sampled cohorts; not calibrated biological truth."
    ),
    "predicted_class": (
        "TARGET if fitted probability meets the decision threshold, otherwise BACKGROUND."
    ),
    "correct": (
        "Agreement of predicted class with the published input class, not proven biological "
        "correctness."
    ),
    "log_odds_contribution": "Signed SHAP contribution; positive moves this model towards TARGET.",
    "absolute_rank": "Within-protein rank by absolute SHAP contribution magnitude.",
    "plot_type": "Graphic type, such as SHAP_BEESWARM, SHAP_GLOBAL_BAR or SHAP_WATERFALL.",
    "file_format": "Published storage/download format, such as TSV, Parquet, PNG, SVG or PDF.",
    "asset_path": "Relative path of a verified graphic or report inside the completed result.",
    "sample_count": "Number of samples represented in the stated model graphic or record.",
    "column_name": "Exact machine-readable field name retained in downloads.",
    "declared_type": "Declared or exported data type; distinct from the biological meaning.",
    "definition": "Meaning, units and interpretation limits of a field or glossary term.",
    "category": "Group of related glossary definitions.",
    "term": "Searchable scientific term, status code or exact table field name.",
    "source": "Authority or origin recorded for this exported value.",
    "region_type": "Type of uploaded interval evidence, for example a separately analysed pocket.",
    "region_id": "Stable identifier of the uploaded residue interval or annotation.",
    "relative_path": "Storage path relative to the whole verified completed result.",
    "size_bytes": "Published file size in bytes; part of integrity verification.",
    "sha256": "SHA-256 content checksum, not a scientific score.",
    "rows": "Number of records, not necessarily distinct proteins or independent blocks.",
    "row_count": "Total canonical-table records, rather than the size of its bounded preview.",
    "complete": "Whether this download contains the whole table rather than a bounded summary.",
    "proteins": "Protein count in the stated group; inherited labels can overlap across groups.",
    "models": "Coordinate/model record count; one protein may have several records.",
    "features": "Distinct feature definitions in the stated evidence family and assessment scope.",
    "assignments": (
        "Protein-label rows; multiple labels per protein can make this exceed protein count."
    ),
    "blocks": "Distinct connected independence blocks, not raw protein counts.",
    "signature_count": (
        "Completed signature rows across comparisons; features may occur in several comparisons."
    ),
    "reviewed_proteins": (
        "Published positive label memberships; automated proposals may be included."
    ),
    "labelled_proteins": (
        "Published positive members of this label, including inherited membership."
    ),
    "observed_proteins": (
        "Published positive members of this label; inherited and provisional labels may overlap."
    ),
    "feature_label": "Readable chart label; the feature ID retains its reproducible identity.",
    "has_validation": (
        "Model status is COMPLETE and held-out target proteins are recorded; inspect metrics."
    ),
    "discovery_members": "Proteins defining the frozen discovery structural component.",
    "projected_members": "Validation proteins assigned to the frozen discovery component.",
    "total_members": "All component members, including validation projections.",
    "discovery_edges": "Passing discovery-to-discovery structural edges forming the component.",
    "tm_threshold": "Minimum TM-score required for the structural component.",
    "coverage": "Minimum aligned fraction required for both proteins in the structural component.",
    "enriched_count": (
        "Completed decision-candidate rows with positive discovery difference and local q <= 0.05."
    ),
    "validated_within_count": (
        "Positive discovery candidates validated within comparison, "
        "including study-wide validation."
    ),
    "validated_study_count": (
        "Positive discovery candidates in the validated study-wide decision tier."
    ),
    "complete_count": (
        "Completed signature rows; a count, not the number of completed comparisons."
    ),
    "insufficient_count": "Signature status rows reporting insufficient independent sample size.",
    "no_signature_count": "Rows reporting a completed screen with no significant signature.",
    "target_blocks": "Pure target independence blocks available in this partition.",
    "background_blocks": "Pure background independence blocks available in this partition.",
    "excluded_mixed_blocks": "Mixed target/background blocks excluded from the comparison.",
    "target_units": (
        "Requested target matching units, including matched and unmatched units "
        "in the stated cohort."
    ),
    "matched_target_units": "Target units allocated accepted controls for this comparison.",
    "excluded_unmatched_target_units": (
        "Requested target units excluded from this comparison "
        "because no accepted control was matched."
    ),
    "analysed_target_proteins": (
        "Target proteins retained in the published matched comparison cohort."
    ),
    "excluded_unmatched_target_proteins": (
        "Target proteins excluded from this comparison because their target unit was unmatched."
    ),
    "matched_control_proteins": (
        "Control protein representatives retained in this comparison cohort."
    ),
    "covered_target_units": (
        "Target units having at least one matched control in this background pool."
    ),
    "matched_control_units": (
        "Distinct matched control units in the stated comparison or pooled audit."
    ),
    "target_coverage_fraction": (
        "Target units with at least one control divided by requested target units."
    ),
    "comparison": "Readable name of a target-versus-background comparison.",
    "display_label": "Short plotting label; hover or download for the full stable identifier.",
    "minimum_coverage_bin": (
        "Lower boundary of the 0.05-wide bin for the smaller protein-pair alignment coverage."
    ),
    "tm_score_bin": (
        "Lower boundary of the 0.05-wide TM-score bin, rather than an exact pair score."
    ),
    "comparison_count": "Number of published protein-pair comparisons in this score/coverage bin.",
    "within_comparison_rank": (
        "Rank by ascending discovery q-value within this comparison; effect and feature "
        "identity break ties. The preview applies a separate row limit to each comparison."
    ),
    "label": "Short plotting label; inspect the stable feature ID and comparison for its identity.",
    "alignment_column": "1-based column in the global sequence alignment, including gap columns.",
    "reference_position": "1-based reference-sequence position; missing at an alignment gap.",
    "reference_residue": "Reference amino acid at this alignment column, or a gap marker.",
    "comparison_position": "1-based partner-sequence position; missing at an alignment gap.",
    "comparison_residue": "Partner amino acid at this alignment column, or a gap marker.",
    "reference_enrichment": (
        "Mapped feature q-value colour score on the reference, not a residue test."
    ),
    "comparison_enrichment": (
        "Mapped feature q-value colour score on the partner, not a residue test."
    ),
    "identity": "Whether both aligned residues exist and have the same amino-acid letter.",
    "protein_count": "Distinct proteins represented by positive features in the stated group.",
    "evidence": "Evidence family or optional stage named in the result metadata.",
    "count": (
        "Number of records in the stated group; inspect whether proteins or blocks were counted."
    ),
}


def column_definition(*, column_name: str) -> str:
    """Resolve a shared definition without inventing semantics for unknown fields.

    Args:
        column_name: Exact table column name.

    Returns:
        Definition including scope or an explicit unknown-field explanation.

    Raises:
        InputValidationError: If the name is not non-empty text.
    """
    if not isinstance(column_name, str) or not column_name.strip():
        raise InputValidationError("A help field name must be non-empty text.")
    name = column_name.strip()
    if name in COLUMN_DEFINITIONS:
        return COLUMN_DEFINITIONS[name]
    for prefix, scope in (
        ("discovery_", "Discovery partition"),
        ("validation_", "Held-out validation"),
    ):
        base = name.removeprefix(prefix)
        if name.startswith(prefix) and base in COLUMN_DEFINITIONS:
            return f"{scope}: {COLUMN_DEFINITIONS[base]}"
    LOGGER.debug("No registered scientific definition for field %s", name)
    return (
        "Additional field from the selected producer. Its scientific meaning is not "
        "registered in this version; consult its evidence source rather than inferring a meaning."
    )
