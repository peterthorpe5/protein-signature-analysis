"""Contextual methods, outcome explanations and recorded campaign limits."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from numbers import Real
from typing import Any

from protein_signatures.errors import InputValidationError
from protein_signatures.result_help import COLUMN_DEFINITIONS, column_definition

PAGE_METHODS = {
    "Overview": (
        "Counts describe published records. Completed signature rows may represent the same "
        "feature in several comparisons and can include target depletion. Discovery selects "
        "features; held-out validation tests frozen definitions on separate connected blocks. "
        "The app reads the verified result and does not rerun the scientific analysis."
    ),
    "Signature explorer": (
        "Two-sided Fisher exact tests use feature-assessed pure independence blocks. "
        "A block is positive if a relevant member has the feature; absence is known only "
        "when every relevant member was assessed without a hit. Mixed and unknown blocks "
        "are excluded. Benjamini-Hochberg correction is applied within each comparison, "
        "partition and feature type, then across comparisons within that same partition "
        "and feature type. Study-wide correction does not pool unrelated feature families."
    ),
    "Explainable prediction": (
        "Elastic-net logistic regression is fitted on discovery proteins with inverse-block-size "
        "and class-balanced weights. Group-aware CV selects the penalty. Held-out metrics "
        "use pure independent blocks with features aggregated by presence; the prediction "
        "table and SHAP plots show individual protein predictions. These are different "
        "evaluation units. Mixed blocks are excluded from metrics. SHAP uses a discovery "
        "reference and contributes in model log-odds, not residue-level significance."
    ),
    "Protein & Pfam": (
        "Labels come from the recorded reviewed or provisional authority. DIRECT membership "
        "is the original class; INHERITED membership follows its parents. Pfam/domain "
        "intervals are 1-based inclusive and refer to the authoritative sequence. Positive "
        "features and assessment states are separate: no row is not an assessed absence."
    ),
    "Classes & roles": (
        "Profile labels define class hierarchy, mechanism and component role independently. "
        "Parent and child memberships overlap, so adding label counts double-counts proteins. "
        "Automated accepted assignments remain provisional. Pooled control totals can serve "
        "several comparisons and are not comparison-specific statistical denominators."
    ),
    "Structures & folds": (
        "Qualifying discovery-to-discovery structural edges define frozen connected components. "
        "Validation proteins project onto those components without merging them. A chain of "
        "passing edges does not guarantee every pair is similar. Named folds require an "
        "explicit fold authority. Whole-model components do not supply local motif residues. "
        "Model eligibility retains the input authority's recorded confidence policy."
    ),
    "Model & alignment explorer": (
        "Positive enriched feature memberships are mapped onto their recorded intervals or "
        "every exact overlapping k-mer occurrence in the published sequence. The colour "
        "uses the chosen feature q-value, with the smallest q at overlaps; it is not a "
        "residue-specific test or pLDDT. Coordinates receive colours only after exact full "
        "sequence and numbering checks. The paired view is a global sequence alignment, "
        "not a structural superposition. Uploaded annotations remain a separate evidence layer."
    ),
    "Orthology & partitions": (
        "Connected OrthoFinder, exact-sequence and supplied near-redundancy links form "
        "indivisible partition blocks. A deterministic seed/hash assigns each whole block "
        "to discovery or validation. The HOG identity is run, group type, hierarchy node "
        "and group ID together. An empty optional group-context table does not remove "
        "the memberships or partitions that were actually published."
    ),
    "Canonical data & downloads": (
        "A preview is bounded and is not the complete dataset. Complete TSV/Parquet downloads "
        "retain exact field names and values. A compact workbook summarises a large table; "
        "the filtered exporter is explicitly capped. Newly generated app workbooks use "
        "the same field dictionary as the tooltips and glossary. Previously published "
        "workbooks retain their original definitions and checksums."
    ),
    "Data quality & provenance": (
        "Completion/checksum verification establishes file integrity, not scientific validity. "
        "Review provisional labels, unmatched target blocks, assessed denominators and "
        "label-defining feature exclusions. Defining domains are removed to reduce direct "
        "circularity; correlated domains and phylogeny can still confound the analysis. "
        "Imported pocket assessment records are not new residue-mapped pocket discoveries."
    ),
    "Glossary & help": (
        "This dictionary covers exact canonical field names, displayed summary fields, "
        "status codes and scientific terms. Search is literal and case-insensitive. "
        "Numerical limits below come from this result's metadata, never today's software "
        "defaults. A missing historical setting is shown as Not recorded."
    ),
}

STATUS_DEFINITIONS = {
    "COMPLETE": (
        "The stated computation completed; inspect evidence class, label status and validation."
    ),
    "INSUFFICIENT_SAMPLE_SIZE": (
        "The relevant class has fewer eligible samples than its recorded limit."
    ),
    "INSUFFICIENT_INDEPENDENT_GROUPS": (
        "Too few pure discovery independence groups in either class."
    ),
    "COMPLETE_VALIDATION_UNDERPOWERED": (
        "Model fitted, but held-out class thresholds failed; metrics are withheld."
    ),
    "COMPLETE_NO_VALIDATION": (
        "Model fitted without held-out proteins; SHAP describes discovery predictions."
    ),
    "CROSS_VALIDATION_FAILED": (
        "No usable group-aware fold contained both classes; no classifier was fitted."
    ),
    "NO_ELIGIBLE_FEATURES": (
        "No eligible feature family survived assessment/variation requirements."
    ),
    "NO_SIGNIFICANT_SIGNATURE": (
        "An eligible discovery screen completed without a local-FDR signature."
    ),
    "NOT_SELECTED": (
        "The optional stage or candidate was not selected, rather than biologically absent."
    ),
    "NO_PUBLISHED_SIGNATURE": (
        "No signature outcome was published for this comparison in this resource."
    ),
    "SIGNIFICANT_POSITIVE": (
        "Positive discovery prevalence difference with local q <= 0.05 in the viewer."
    ),
    "DISCOVERY_ENRICHED": (
        "At least one completed positive discovery decision candidate meets local q <= 0.05; "
        "inspect held-out and study-wide validation separately."
    ),
    "COMPLETE_NO_POSITIVE_ENRICHMENT": (
        "Completed rows exist but none pass the viewer's positive-enrichment rule."
    ),
    "PROVISIONAL_EVIDENCE_SUPPORTED": (
        "Automated evidence-supported proposal for hypothesis generation."
    ),
    "REVIEWED_POSITIVE": (
        "Positive assignment explicitly reviewed under the recorded human authority."
    ),
    "ACCEPTED": (
        "Accepted by an automated evidence rule; it remains a provisional biological label."
    ),
    "PROPOSED": "Some rule support, but the acceptance conditions were not met.",
    "AMBIGUOUS": "Incompatible supported labels could not be resolved under the recorded rule.",
    "UNMAPPED": "No configured label evidence supported an assignment.",
    "MATCHED": "A caliper-compatible control unit was allocated to this target matching request.",
    "UNMATCHED": "No eligible control was allocated; not evidence the target lacks a feature.",
    "ASSESSED_WITH_HIT": (
        "The authority completed an assessment and found at least one retained hit."
    ),
    "ASSESSED_NO_HIT": "The authority completed an assessment and found no retained hit.",
    "NOT_ASSESSED": "No completed assessment under this authority; absence is unknown.",
    "FAILED": "Assessment failed; no-hit or absence cannot be inferred.",
    "EXCLUDED": "Evidence was deliberately excluded under the declared analysis policy.",
    "AVAILABLE": "Coordinates were available; this alone does not establish analysis eligibility.",
    "UNAVAILABLE": "The source did not provide usable evidence or coordinates.",
    "ELIGIBLE": "The model passed the input authority's structural eligibility rule.",
    "INELIGIBLE_LOW_CONFIDENCE": (
        "The model failed its source's mean-confidence eligibility threshold."
    ),
    "DISCOVERY_COMPONENT": "Member of the frozen discovery structural component.",
    "VALIDATION_PROJECTION": "Validation protein projected onto a frozen discovery component.",
    "ALL_DATA_COMPONENT": (
        "Helper component derived using all data; exploratory, not held-out confirmation."
    ),
    "FIXED_EXTERNAL": "Feature definition came from a fixed external authority.",
    "DISCOVERY_DERIVED": "Feature definition used discovery data only.",
    "DISCOVERY_SUPERVISED": (
        "Label-aware discovery definition retained for audit, excluded from inference."
    ),
    "ASSESSED_WITH_FEATURE": (
        "Feature assessment completed with retained positive feature evidence."
    ),
    "ASSESSED_NO_FEATURE": (
        "Feature assessment completed with a known absence under its declared scope."
    ),
    "ALL_DATA_EXPLORATORY": (
        "Definition used all data and cannot claim independent held-out confirmation."
    ),
    "IMPORTED_COMPLETE": (
        "The upstream resource completed its assessment; not a new calculation by this package."
    ),
    "DIRECT": "Original assigned label membership.",
    "INHERITED": "Parent-label membership inherited from an original direct label.",
    "UNKNOWN": "The source did not specify this value; not proof of an absent function.",
    "NO_ELIGIBLE_MODELS": (
        "No configured comparison met the requirements for fitting a classifier."
    ),
    "CONVERGENCE_WARNING": (
        "The optimiser reached its iteration limit; inspect coefficients and rerun policy."
    ),
    "INPUT_UNAVAILABLE": "A required input was not supplied; no absence assessment was made.",
    "AMBIGUOUS_CLASS": "The relevant target/background class could not be assigned unambiguously.",
    "MAPPING_FAILED": "Evidence could not be mapped to the authoritative protein identifiers.",
    "EVIDENCE_SUPPORTED_POSITIVE": "Automated accepted positive label; not human-reviewed truth.",
    "REVIEWED_NEGATIVE": "Explicit reviewed negative assignment under its recorded authority.",
    "REVIEWED_COMPONENT_NOT_CATALYTIC": (
        "Reviewed system component without a catalytic-positive assignment."
    ),
    "NOT_APPLICABLE_EXTERNAL_EVIDENCE": (
        "Imported evidence does not require local coordinate eligibility."
    ),
    "INELIGIBLE_CONFIDENCE_UNAVAILABLE": "Required model confidence was not recorded.",
    "INELIGIBLE_SEQUENCE_UNVERIFIED": (
        "The model sequence could not be verified against the campaign sequence."
    ),
    "INELIGIBLE_SEQUENCE_MISMATCH": (
        "The model sequence differs from the authoritative campaign sequence."
    ),
    "INELIGIBLE_COORDINATE_UNAVAILABLE": "Usable local coordinates were unavailable.",
    "INELIGIBLE_USER_EXCLUDED": "The input authority explicitly excluded this model from analysis.",
    "INELIGIBLE_VALIDATION_FAILED": "The model failed the recorded coordinate validation policy.",
    "SUCCESS": "The structural comparison producer reported successful completion.",
    "PASS": (
        "The structural comparison producer reported its passing state; inspect its thresholds."
    ),
    "FULL_SEQUENCE": "Structural coverage denominator is the full authoritative sequence length.",
    "STRUCTURE_MODEL_RESIDUES": "Structural coverage denominator is the residues in the model.",
    "DOMAIN_OR_CONSTRUCT": "Structural coverage denominator is a supplied domain or construct.",
    "AUTOMATED_EVIDENCE_CONFLICT": "Conflicting automated label evidence; no accepted assignment.",
    "AUTOMATED_EVIDENCE_BELOW_POLICY": (
        "Some label evidence exists but fails the acceptance policy."
    ),
    "AUTOMATED_EVIDENCE_NOT_MAPPED": "No configured automated label evidence mapped this protein.",
    "PROVISIONAL_MATCHED_CONTROL": (
        "Automatically selected control pool; not a reviewed biological negative."
    ),
    "NOT_RUN": "This computation has not been run; no result can be interpreted.",
    "SEQUENCE_MISMATCH": (
        "Acquired model sequence differs from the campaign FASTA; model was rejected."
    ),
    "ACQUIRED": (
        "AlphaFold coordinates acquired and passed the recorded sequence/confidence policy."
    ),
    "ACQUIRED_SEQUENCE_UNVERIFIED": (
        "Coordinates acquired but sequence unverified; ineligible for analysis."
    ),
    "ACQUIRED_CONFIDENCE_UNAVAILABLE": (
        "Coordinates acquired without required confidence; ineligible."
    ),
    "ACQUIRED_LOW_CONFIDENCE": (
        "Coordinates acquired below the recorded confidence minimum; ineligible."
    ),
    "MODEL_NOT_AVAILABLE": (
        "No AlphaFold model was returned for this request; not proof of no structure."
    ),
    "NO_VALIDATION_DATA": "No completed eligible held-out test exists for this discovery feature.",
    "DIRECTION_DISCORDANT": "Held-out prevalence difference has the opposite sign to discovery.",
    "VALIDATED_STUDY_WIDE": (
        "Same-direction held-out support passes local and study-wide FDR limits."
    ),
    "VALIDATED_WITHIN_COMPARISON": "Same-direction held-out support passes local FDR only.",
    "DISCOVERY_ONLY": (
        "Discovery passed its decision rule; held-out local significance did not pass."
    ),
}

EVIDENCE_SCOPE_HELP = {
    "DECISION_CANDIDATE": "Discovery passed local and study-wide FDR within its feature family.",
    "EXPLORATORY_LOCAL_ONLY": "Discovery passed local FDR without study-wide support.",
    "EXPLORATORY_ALL_DATA_DERIVATION": (
        "Definition used all data; held-out wording cannot restore independence."
    ),
    "QC_TECHNICAL_NON_BIOLOGICAL": (
        "Technical evidence retained for quality control, not a biological signature."
    ),
}
STATUS_DEFINITIONS.update(
    {
        f"{scope}__{validation}": f"{definition} {STATUS_DEFINITIONS[validation]}"
        for scope, definition in EVIDENCE_SCOPE_HELP.items()
        for validation in (
            "NO_VALIDATION_DATA",
            "DIRECTION_DISCORDANT",
            "VALIDATED_STUDY_WIDE",
            "VALIDATED_WITHIN_COMPARISON",
            "DISCOVERY_ONLY",
        )
    }
)
STATUS_DEFINITIONS.update(
    {
        f"{code}_CONVERGENCE_WARNING": (
            f"{STATUS_DEFINITIONS[code]} {STATUS_DEFINITIONS['CONVERGENCE_WARNING']}"
        )
        for code in ("COMPLETE", "COMPLETE_NO_VALIDATION", "COMPLETE_VALIDATION_UNDERPOWERED")
    }
)

METRIC_HELP = {
    "cv_roc_auc": (
        "ROC is receiver operating characteristic; AUC is area under that curve. "
        "0.5 means chance ranking, 1 perfect ranking, below 0.5 reversed ranking. "
        "CV is group-aware cross-validation inside discovery and selects the penalty. "
        "Use final held-out validation for independent performance; AUC is not accuracy."
    ),
    "validation_roc_auc": (
        "Receiver operating characteristic area under the curve across decision thresholds. "
        "On held-out pure blocks, 0.875 means a random target is ranked above a random "
        "background about 87.5% of the time, counting ties by half. 0.5 is chance; 1 is "
        "perfect ranking. This does not measure calibration or 87.5% correct classifications."
    ),
    "validation_average_precision": (
        "Average precision (AP) summarises precision versus recall, from 0 to 1. "
        "Precision is the target fraction among positive calls; recall is the fraction "
        "of targets recovered. Higher AP is better, relative to target prevalence in "
        "the evaluated pure-block cohort. AP is not trapezoidal PR AUC. No universal "
        "good-value cut-off applies across differently balanced cohorts."
    ),
    "validation_matthews_correlation": (
        "Matthews correlation coefficient (MCC) summarises all four confusion-matrix cells "
        "at the recorded threshold (0.5 in this pipeline). +1 is perfect, 0 no correlation, "
        "-1 perfectly reversed. A value of 0.604 is a positive association between the "
        "calls and supplied labels, not 60.4% accuracy or proof of biological function."
    ),
    "validation_balanced_accuracy": (
        "Mean of target sensitivity and background specificity at threshold 0.5. "
        "0.5 is chance and 1 is perfect; unlike raw accuracy it weights both classes equally."
    ),
    "validation_brier_score": (
        "Mean squared error between predicted probability and binary class on pure blocks. "
        "0 is perfect and lower is better. Compare with a constant-prevalence predictor; "
        "this score alone does not establish calibration on other organisms or proteomes."
    ),
}

EXTRA_TERMS = (
    ("Prediction", "ROC", METRIC_HELP["validation_roc_auc"]),
    ("Prediction", "AUC", METRIC_HELP["cv_roc_auc"]),
    ("Prediction", "AP / average precision / PR AUC", METRIC_HELP["validation_average_precision"]),
    ("Prediction", "MCC", METRIC_HELP["validation_matthews_correlation"]),
    ("Prediction", "CV / cross-validation", METRIC_HELP["cv_roc_auc"]),
    ("Prediction", "Balanced accuracy", METRIC_HELP["validation_balanced_accuracy"]),
    ("Prediction", "Brier score", METRIC_HELP["validation_brier_score"]),
    (
        "Prediction",
        "SHAP",
        "Additive contributions to a fitted prediction, not causal or residue-level tests.",
    ),
    (
        "Prediction",
        "Log-odds",
        "Log(p/(1-p)); 0 means p=0.5. Additive effects here are not probabilities.",
    ),
    (
        "Prediction",
        "Elastic net",
        "Logistic regression regularised by a mixture of L1 and L2 penalties.",
    ),
    (
        "Statistics",
        "Insufficient sample size",
        "See the recorded campaign limits: blocks for association, proteins for ML.",
    ),
    (
        "Statistics",
        "Independent block",
        "Connected orthology/redundancy unit counted once; distinct from raw proteins.",
    ),
    (
        "Statistics",
        "Fisher exact test",
        "Two-sided association test on a feature's assessed target/background block counts.",
    ),
    (
        "Statistics",
        "Wilson interval",
        "95% interval for one class prevalence; not the confidence interval of its difference.",
    ),
    ("Structure", "HOG", "Hierarchical orthogroup at an exact OrthoFinder hierarchy node."),
    (
        "Structure",
        "Pfam",
        "Protein-domain/family authority with explicit assessed-hit and no-hit states.",
    ),
    ("Structure", "RMSD", COLUMN_DEFINITIONS["rmsd_angstrom"]),
    ("Structure", "Cα", "Alpha-carbon atom trace used for the simplified model display."),
    (
        "Structure",
        "mmCIF / PDB",
        "Coordinate-file formats; neither proves predicted models have biological activity.",
    ),
    (
        "Features",
        "k-mer",
        "Exact k-residue amino-acid word. k3:LPD means Leu-Pro-Asp anywhere in a sequence.",
    ),
    (
        "Features",
        "Pfam architecture",
        "Ordered domain identifiers along a protein, not one contiguous sequence motif.",
    ),
    (
        "Structure",
        "SC_ feature ID",
        "Opaque stable key of a whole-model similarity component; not a named fold or local motif.",
    ),
    (
        "Study design",
        "Caliper",
        "Prespecified maximum covariate difference allowed when matching controls.",
    ),
    (
        "Study design",
        "Circularity",
        "Reusing label-defining evidence as independent confirmation of that label.",
    ),
    (
        "Exports",
        "Preview",
        "Bounded visible rows; TSV of this preview is not automatically the complete "
        "canonical table.",
    ),
    ("Features", "Feature ID", COLUMN_DEFINITIONS["feature_id"]),
    (
        "Statistics",
        "Enrichment colour",
        "Display of feature q-value at mapped intervals, not residue-specific significance.",
    ),
    (
        "Exports",
        "TSV",
        "Tab-separated UTF-8 text with exact column names and no spreadsheet formulas.",
    ),
    ("Exports", "Parquet", "Typed columnar storage for complete canonical datasets."),
    (
        "Exports",
        "SHA-256 / checksum",
        "Content fingerprint for integrity verification; not a measure of scientific validity.",
    ),
    ("Prediction", "Precision", "Target fraction among positive calls at a stated threshold."),
    (
        "Prediction",
        "Recall / sensitivity",
        "Fraction of supplied targets correctly called target at a stated threshold.",
    ),
    (
        "Prediction",
        "Specificity",
        "Fraction of supplied backgrounds correctly called background at a stated threshold.",
    ),
    (
        "Prediction",
        "Calibration",
        "Agreement between predicted probabilities and observed class frequencies in a cohort.",
    ),
)

LIMIT_SPECS = (
    (
        "analysis",
        "minimum_target_proteins",
        "Association target minimum",
        "assessed pure blocks per partition",
    ),
    (
        "analysis",
        "minimum_background_proteins",
        "Association background minimum",
        "assessed pure blocks per partition",
    ),
    (
        "analysis",
        "minimum_feature_proteins",
        "Candidate prevalence minimum",
        "positive proteins for feature selection",
    ),
    (
        "explainable_ml",
        "minimum_samples_per_class",
        "Model sample minimum",
        "pure-group proteins in EACH class",
    ),
    (
        "explainable_ml",
        "minimum_groups_per_class",
        "Model group minimum",
        "pure groups in EACH class",
    ),
    (
        "explainable_ml",
        "cross_validation_folds",
        "Requested CV folds",
        "actual completed folds are in ml_models",
    ),
    ("analysis", "fdr_threshold", "Campaign FDR threshold", "association decision rule"),
    (
        "analysis",
        "validation_fraction",
        "Requested validation fraction",
        "blocks; actual class allocation varies",
    ),
    (
        "analysis",
        "structural_tm_score_threshold",
        "Structural TM-score threshold",
        "frozen component construction/projection",
    ),
    (
        "analysis",
        "structural_minimum_coverage",
        "Structural coverage threshold",
        "fraction required for BOTH proteins",
    ),
    (
        "alphafold",
        "minimum_mean_plddt",
        "AlphaFold acquisition confidence",
        "0-100; imported models retain source policy",
    ),
)

PAGE_TERMS = {
    "Overview": ("Independent block", "FDR", "Discovery partition", "Validation partition"),
    "Signature explorer": ("Fisher exact test", "q-value", "Study-wide q-value", "Wilson interval"),
    "Explainable prediction": ("ROC", "AUC", "AP / average precision / PR AUC", "MCC", "SHAP"),
    "Protein & Pfam": ("Pfam", "k-mer", "Pfam architecture", "Coordinate projection"),
    "Classes & roles": ("Target class", "Background class", "Circularity", "Caliper"),
    "Structures & folds": ("Structural cluster", "TM-score", "Bilateral coverage", "pLDDT", "RMSD"),
    "Model & alignment explorer": (
        "k-mer",
        "SC_ feature ID",
        "Cα",
        "Coordinate projection",
        "Sequence alignment",
    ),
    "Orthology & partitions": (
        "HOG",
        "Independent block",
        "Discovery partition",
        "Validation partition",
    ),
    "Canonical data & downloads": ("Preview", "Completed result", "Feature ID"),
    "Data quality & provenance": ("Caliper", "Circularity", "Independent block", "Pocket"),
    "Glossary & help": ("Insufficient sample size", "ROC", "AUC", "MCC"),
}

GRAPH_HELP = {
    "signature_types": (
        "Bars count completed signature rows by evidence family across comparisons. "
        "One feature can occur in several comparisons. Count is not effect size, and "
        "completion does not establish held-out validation or positive enrichment."
    ),
    "feature_coverage": (
        "Bars show proteins represented in positive feature records, coloured by distinct "
        "feature definitions. This is evidence coverage, not absence assessment or enrichment."
    ),
    "signature_scatter": (
        "Each point is a feature. The horizontal axis is target-minus-background block "
        "prevalence; +0.20 means 20 percentage points more common in targets. The vertical "
        "axis can show −log10 of the local discovery q-value, so q=0.01 is 2 and "
        "q=1e-8 is 8. Recorded zeros and values below 1e-300 appear at a clearly stated "
        "display cap of 300, with the original value and clipping flag in hover. "
        "The optional recorded-value axis is reversed so smaller q-values appear higher. "
        "The q=0.05 line is a viewer reference, not a replacement for the campaign's rule. "
        "A near-zero q-value is not a large biological effect. Check held-out and "
        "study-wide evidence."
    ),
    "model_coefficients": (
        "Bars show signed logistic coefficients for individual features in log-odds. "
        "Right of zero favours TARGET, left favours BACKGROUND. They are conditional "
        "fitted effects, not enrichment p-values. Correlated features must be considered together."
    ),
    "validation_individual_shap": (
        "Bars show up to 30 individual features ranked by their published mean absolute "
        "SHAP contribution across validation proteins, in model log-odds. No other-features "
        "aggregate is included. Magnitude does not show target/background direction or "
        "biochemical causation. The ordering differs from coefficient ranks, and these "
        "protein explanations do not use the pure-block metric evaluation unit."
    ),
    "probabilities": (
        "Histograms show individual protein predictions separated by supplied class. "
        "Overlap reveals ambiguity; separation is descriptive. Reported validation metrics "
        "instead evaluate aggregated pure blocks. Probabilities near 1 are model output, "
        "not proof of function, calibration or transfer to another organism."
    ),
    "class_label_coverage": (
        "Bars show positive label memberships, including inherited parent memberships. "
        "Parent and child bars overlap, so their sum is not a unique-protein total. "
        "These counts are descriptive and are not each comparison's assessed denominators."
    ),
    "structure_eligibility_summary": (
        "Bars count coordinate records by recorded eligibility, coloured by fold-assessment "
        "state. Eligible coordinates can still lack named fold assignments. Mean model "
        "confidence is not evidence that a local motif or biochemical function is correct."
    ),
    "structural_alignment_score_coverage": (
        "Each point is a published protein pair: TM-score versus the smaller aligned fraction "
        "of the two proteins. High similarity on a short region may not support whole-model "
        "equivalence. Inspect tool, coverage scope and coordinates; this plot supplies no "
        "aligned residue pairs or residue-specific significance."
    ),
    "enrichment_map": (
        "Colour marks occurrences of positively enriched features along the sequence, "
        "using the selected q-value tier. Blue is weaker significant evidence and red stronger; "
        "white is no mapped significant evidence, not assessed absence. Overlaps use the "
        "smallest feature q-value. This is not a residue-by-residue statistical test or pLDDT."
    ),
    "trace": (
        "The rotatable Cα trace shows the selected model. Colour carries the same mapped "
        "feature evidence as the sequence view. A grey model indicates projection is "
        "unavailable or unsafe after sequence/numbering checks; it does not imply no motif."
    ),
    "alignment": (
        "Columns show a global amino-acid sequence alignment (+2 match, -1 substitution, "
        "-2 gap), with each protein's mapped enrichment. Gaps have no residue. This is "
        "not a 3D superposition; structural scores shown elsewhere are aggregate comparisons."
    ),
    "ranked_association": (
        "The strongest visible positive features are compared by prevalence difference. "
        "The display is bounded and can span several comparisons. q-values and stable IDs "
        "are in hover/downloads. Ranking is prioritisation, not an independent residue test."
    ),
    "SHAP_BEESWARM": (
        "Each point represents one explained protein. Right of zero favours TARGET; left favours "
        "BACKGROUND relative to the discovery reference. Red/blue shows high/low input "
        "feature value, not statistical significance. The grouped other-features row "
        "combines omitted contributions; it is not one motif. Read the explained partition."
    ),
    "SHAP_GLOBAL_BAR": (
        "Bars show mean absolute SHAP log-odds contribution, so they have magnitude but "
        "no direction. The sum of other features can dominate because it combines many "
        "features; it is not a single discovered motif. The adjacent table describes "
        "individual coefficients and permutation importance; its coefficient ranks "
        "are not the SHAP bar ordering."
    ),
    "SHAP_WATERFALL": (
        "Start at the discovery-reference expected log-odds. Red contributions move the "
        "model towards TARGET and blue towards BACKGROUND, reaching this protein's log-odds. "
        "Values are not probabilities or causal effects; correlated features share attribution."
    ),
    "control_matching_coverage": (
        "Bars show the fraction of requested target units with at least one matched control "
        "in each background pool. A low fraction can select a restricted target subset. "
        "The pool may serve several comparisons; use the association table's assessed "
        "pure-block counts for each statistical denominator."
    ),
    "partition_allocation": (
        "Bars count whole connected independence blocks by partition and block authority. "
        "Protein counts can have different proportions because blocks differ in size. "
        "A requested split fraction does not guarantee exact class or protein balance; "
        "inspect each comparison's eligible pure-block counts separately."
    ),
}


def graph_explanation(*, graph_name: str) -> str:
    """Resolve interpretation help for one chart or a safe generic fallback.

    Args:
        graph_name: Stable chart download identity or published SHAP plot type.

    Returns:
        Chart-specific axes, units and interpretation limits.

    Raises:
        InputValidationError: If the name is not non-empty text.
    """
    if not isinstance(graph_name, str) or not graph_name.strip():
        raise InputValidationError("A graph help identity must be non-empty text.")
    if "three_dimensional_model" in graph_name:
        return GRAPH_HELP["trace"]
    if "top_enriched_features" in graph_name:
        return GRAPH_HELP["ranked_association"]
    for token, definition in GRAPH_HELP.items():
        if token in graph_name:
            return definition
    return (
        "Read the labelled axes, units, legend and hover values. This plot describes the "
        "selected published evidence; inspect assessed sample counts, comparison policy "
        "and validation before making a biological claim."
    )


def metric_reading(*, metric_name: str, value: object, baseline: float | None = None) -> str:
    """Interpret a numeric model metric without imposing arbitrary quality grades.

    Args:
        metric_name: Registered model metric field.
        value: Recorded metric or an unavailable value.
        baseline: Target fraction among evaluated pure validation blocks, if known.

    Returns:
        Numeric interpretation with a metric-specific baseline where available.

    Raises:
        InputValidationError: If the metric identity is unknown.
    """
    if metric_name not in METRIC_HELP:
        raise InputValidationError(f"Unknown model metric {metric_name!r}.")
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        return "Not calculated or not recorded; inspect model status and per-class sample limits."
    lower = -1 if metric_name == "validation_matthews_correlation" else 0
    if not lower <= value <= 1:
        return "Outside this metric's valid range; inspect the published model record."
    if metric_name in {"cv_roc_auc", "validation_roc_auc"}:
        scope = "discovery CV" if metric_name == "cv_roc_auc" else "held-out pure-block"
        return (
            f"{value:.3f}: about {100 * value:.1f}% target-versus-background pair ordering "
            f"under {scope} evaluation, with ties counted by half. This is not accuracy."
        )
    if metric_name == "validation_matthews_correlation":
        return f"{value:+.3f} correlation of class calls with supplied labels; 0 is no correlation."
    if metric_name == "validation_balanced_accuracy":
        return f"{100 * value:.1f}% mean of target sensitivity and background specificity."
    known_baseline = (
        not isinstance(baseline, bool)
        and isinstance(baseline, Real)
        and math.isfinite(baseline)
        and 0 <= baseline <= 1
    )
    if metric_name == "validation_average_precision":
        reference = (
            f"{baseline:.3f} target prevalence"
            if known_baseline
            else "unrecorded target prevalence"
        )
        return (
            f"AP {value:.3f}; compare with {reference} in evaluated blocks, "
            "not a universal cut-off."
        )
    reference = f"{baseline * (1 - baseline):.3f}" if known_baseline else "not available"
    return f"Mean squared probability error {value:.3f}; constant-prevalence reference {reference}."


def feature_explanation(*, feature_type: str, feature_id: str) -> str:
    """Explain a stable feature key while retaining its biological claim limits.

    Args:
        feature_type: Published evidence family.
        feature_id: Exact published feature identifier.

    Returns:
        Readable definition independent of class labels or discovered outcomes.

    Raises:
        InputValidationError: If either identifier is not text.
    """
    if not isinstance(feature_type, str) or not isinstance(feature_id, str):
        raise InputValidationError("Feature explanations require text identifiers.")
    if feature_type == "AMINO_ACID_KMER":
        match = re.fullmatch(r"k([1-9][0-9]*):([A-Z]+)", feature_id)
        if match is not None and int(match.group(1)) == len(match.group(2)):
            return (
                f"Exact {match.group(1)}-residue amino-acid word {match.group(2)} "
                "anywhere in the sequence; occurrence is not a functional or structural "
                "motif claim."
            )
        return "Unrecognised exact k-mer key; consult its original feature definition."
    if feature_type == "STRUCTURE_CLUSTER":
        return (
            "Whole-model structural-similarity component; no named fold or local residue interval."
        )
    if feature_type in {"PFAM_ARCHITECTURE", "DOMAIN_ARCHITECTURE"}:
        return (
            "Ordered detected domain identifiers along a protein, not a contiguous "
            "amino-acid motif."
        )
    if feature_type in {"PFAM_DOMAIN", "DOMAIN"}:
        return (
            "Presence of the named domain under its scan authority; consult its residue intervals."
        )
    if feature_type == "FOLD":
        return (
            "Named fold assignment under an explicit authority, not a newly discovered local motif."
        )
    if feature_type == "STRUCTURAL_POCKET":
        return (
            "Imported pocket evidence; assessment alone supplies no enrichment-tested "
            "residue interval."
        )
    if feature_type in {
        "STRUCTURE_AVAILABILITY",
        "STRUCTURE_QUALITY",
        "STRUCTURE_AVAILABLE",
        "EVIDENCE_AVAILABILITY",
    }:
        return "Technical model availability/quality feature, not a biochemical signature."
    return "Feature from the declared producer; consult its source and assessment authority."


def recorded_setting(*, metadata: Mapping[str, Any], section: str, key: str) -> str:
    """Read one numerical limit without replacing missing metadata with defaults.

    Args:
        metadata: Published run metadata.
        section: Campaign configuration section.
        key: Exact numerical setting.

    Returns:
        Recorded value or an explicit unavailable label.
    """
    campaign = metadata.get("campaign", {}) if isinstance(metadata, Mapping) else {}
    settings = campaign.get(section, {}) if isinstance(campaign, Mapping) else {}
    value = settings.get(key) if isinstance(settings, Mapping) else None
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        return "Not recorded"
    if (
        key.startswith("minimum_") and key != "minimum_mean_plddt"
    ) or key == "cross_validation_folds":
        if value < 1 or int(value) != value:
            return "Not recorded"
    elif key == "minimum_mean_plddt":
        if not 0 <= value <= 100:
            return "Not recorded"
    elif not 0 <= value <= 1:
        return "Not recorded"
    return f"{value:g}"


def campaign_limit_rows(*, metadata: Mapping[str, Any]) -> tuple[dict[str, str], ...]:
    """Return all recorded decision limits with their true evaluation units.

    Args:
        metadata: Published run metadata.

    Returns:
        Deterministically ordered limit descriptions for the help panel.
    """
    return tuple(
        {
            "setting": f"{section}.{key}",
            "limit": recorded_setting(metadata=metadata, section=section, key=key),
            "meaning": label,
            "unit": unit,
        }
        for section, key, label, unit in LIMIT_SPECS
    )


def status_explanation(
    *, status: str, metadata: Mapping[str, Any], context: str = "association"
) -> str:
    """Explain an outcome with the appropriate recorded sample requirements.

    Args:
        status: Exact published status code.
        metadata: Published run metadata.
        context: Association or model status context.

    Returns:
        Human-readable outcome and, where relevant, recorded numerical limits.

    Raises:
        InputValidationError: If the status or context is invalid.
    """
    if not isinstance(status, str) or not status.strip() or context not in {"association", "model"}:
        raise InputValidationError("Status help requires a non-empty code and valid context.")
    code = status.strip()
    warning = code.endswith("_CONVERGENCE_WARNING")
    base = code.removesuffix("_CONVERGENCE_WARNING") if warning else code
    text = STATUS_DEFINITIONS.get(
        base, "Producer-specific state; consult its recorded reason and evidence source."
    )
    if base == "INSUFFICIENT_SAMPLE_SIZE" and context == "association":
        target = recorded_setting(
            metadata=metadata, section="analysis", key="minimum_target_proteins"
        )
        background = recorded_setting(
            metadata=metadata, section="analysis", key="minimum_background_proteins"
        )
        text += (
            f" Association requires at least {target} target and {background} background "
            "independent blocks in the tested partition. Feature tests further require "
            "those blocks to be assessed. Mixed blocks and unknown absences do not count. "
            "Despite the configuration names ending in proteins, these are BLOCK limits."
        )
    if (
        base
        in {
            "INSUFFICIENT_SAMPLE_SIZE",
            "INSUFFICIENT_INDEPENDENT_GROUPS",
            "COMPLETE_VALIDATION_UNDERPOWERED",
        }
        and context == "model"
    ):
        samples = recorded_setting(
            metadata=metadata, section="explainable_ml", key="minimum_samples_per_class"
        )
        groups = recorded_setting(
            metadata=metadata, section="explainable_ml", key="minimum_groups_per_class"
        )
        text += (
            f" Each discovery class requires {samples} protein samples AND {groups} pure "
            "independent groups. Held-out metrics require the same two per-class limits "
            "in validation after mixed groups are excluded. Meeting a software minimum "
            "does not establish adequate scientific power. Read the published status_message."
        )
    if warning:
        text += " " + STATUS_DEFINITIONS["CONVERGENCE_WARNING"]
    return text


def glossary_rows(
    *, base_rows: tuple[tuple[str, str, str], ...]
) -> tuple[tuple[str, str, str], ...]:
    """Combine scientific terms with every canonical and displayed field definition.

    Args:
        base_rows: Existing glossary entries to preserve.

    Returns:
        Searchable definitions in deterministic category and term order.
    """
    entries = {(category, term): definition for category, term, definition in base_rows}
    entries.update({(category, term): definition for category, term, definition in EXTRA_TERMS})
    entries.update(
        {("Status codes", code): definition for code, definition in STATUS_DEFINITIONS.items()}
    )
    fields = set(COLUMN_DEFINITIONS)
    for prefix in ("discovery_", "validation_"):
        fields.update(
            prefix + name
            for name in (
                "q_value",
                "study_q_value",
                "prevalence_difference",
                "target_prevalence",
                "background_prevalence",
            )
        )
    entries.update({("Table fields", name): column_definition(column_name=name) for name in fields})
    entries.update(
        {
            ("Campaign settings", f"{section}.{key}"): f"{label}; units: {unit}. "
            "Read the actual value in the result's recorded-limit panel. "
            "Missing values are not defaulted."
            for section, key, label, unit in LIMIT_SPECS
        }
    )
    return tuple(
        (category, term, definition) for (category, term), definition in sorted(entries.items())
    )
