# App critical review and contextual help

Review date: 5 October 2026.

Original review base: `939f53b55353022e9120102af21529fc29cf9075` (`start up speed`).
Compatible patch base: current `main`,
`41f8ffe9ea010c7c9d96e6b2b8bb266596ed2549` (`back up`).

## Scope and assessment

The app has a useful evidence-exploration structure, immutable result verification,
explicit unknown/failed assessment states, independent discovery/validation blocks and
portable data exports. Its main weakness is interpretation: readers can mistake software
completion, provisional labels, feature-level q-values and model predictions for stronger
biological conclusions than the analysis supports. Wide tables and opaque feature IDs
make those mistakes easier.

This review examined all eleven pages, the implementation of their calculations and
exports, and all thirteen supplied screenshots. Twelve screenshots show this protein
signature app; the 09:33:28 screenshot shows the separate OrthoFinder app and was used
only as a presentation reference. Several screenshots already contain local help and
layout changes absent from the original review base. Those changes were subsequently
pushed in `41f8ffe` while the original update still targeted `939f53b`. The compatible
patch reconciles both versions and retains the pushed commit as its base. The original
working checkout, original update ZIP and full recovery ZIP remain preserved.

The full user campaign database and the earlier Word document were not available in this
workspace. Screenshot values are observations, not independently recomputed campaign
results. Automated app verification uses the repository's offline example and focused
DuckDB fixtures. It cannot establish scientific power, external generalisation, performance
on the full campaign or pixel-perfect behaviour on the user's Mac/browser.

## Compatibility with the pushed update

The combined app retains the newer comparison outcome summary, discovery/validation
evidence-tier shortlists, model choices, comparison-specific matched cohorts,
full-pair structural score/coverage density view and aligned residue letter strip.
It adds the recorded-threshold help, field/status glossary and numeric metric
interpretation, alongside the selectable q-value scatter and individual SHAP chart.
Detailed label audits load only after an explicit request, and their preview/export
is capped at 5,000 rows. All upstream app functions remain except `_page_guidance`,
which is replaced by the more complete `_render_page_help`.

The apply script checks the bundle, backs up a clean main checkout and Git history,
pulls using `--ff-only`, prints an explicit reason for a base mismatch, checks the
patch and then applies the overlay with rsync. It does not commit or push.

## Highest-priority findings

| Finding | Consequence | Change or remaining requirement |
| --- | --- | --- |
| “Insufficient sample size” did not reveal the counted unit or saved limits | A reader could compare a protein total with a block threshold or apply current defaults to an old run | Every page shows limits from its own result metadata; table status help distinguishes association blocks from modelling protein/group minimums |
| “Validation PR AUC” displayed `average_precision_score` | The card named a different numerical summary | Renamed to **Validation average precision (AP)**; explains AP versus trapezoidal PR AUC and gives the pure-block target-fraction reference |
| Model metrics, prediction rows and SHAP graphics have different units | Visually separated protein predictions could appear inconsistent with block-level metrics | Explains protein-level predictions/SHAP versus aggregated pure-block evaluation and excludes mixed blocks from reference counts |
| Statistical and technical terms were scattered or unexplained | Readers could not reliably interpret tables without source-code knowledge | Shared field dictionary, table-header help, adjacent expandable field/status definitions, complete searchable glossary and page-specific methods/terms |
| Parent/child class memberships overlap | A sunburst suggests additive parts of a unique total and can double-count members | Replaced by a positive-membership bar chart with parent/child overlap explained; removes the misleading “reviewed” count name |
| All four label-audit tabs queried their full datasets eagerly | Large hidden audits increased rerun cost and memory use | Only the requested selected audit is loaded, explicitly capped at 5,000 rows; complete downloads remain in the canonical-data browser |
| Empty optional OrthoFinder context looked like absent orthology | The page could contradict large membership/group counts | Empty-context message now preserves the authority of published memberships and partitions |
| A high score or tiny q-value could be read as biological certainty | Prioritisation evidence could become an unsupported function claim | Adds precise interpretation, assessed denominators, provisional-label cautions and feature-versus-residue distinctions |

## Page-by-page evaluation

| Page | Critical assessment | Delivered response | Remaining limit |
| --- | --- | --- | --- |
| Overview | Raw signature-table row count includes status placeholders and repeated features across comparisons. Counts are not effect size or independent sample size. | Renames the card **Signature outcome rows**; adds metric tooltips, methods, terms and recorded limits. Existing completed-signature chart remains separate. | Full feature coverage is still an explicit expensive action. No full-campaign timing benchmark was available. |
| Signature explorer | Very small q-values do not establish a large effect. Local and study-wide q-values have different scopes; absence requires assessment. | Retains evidence-tier shortlists across feature families; header help and scientific notation; selectable −log10 q-value chart with a stated display cap of 300, recorded-value hover and a clipping flag for q-values below 1e-300; original q-values remain unchanged. Explains assessed pure-block denominators and multiplicity. | The display cannot recover significance beyond a recorded numerical zero. Missing/non-finite/out-of-range plot coordinates are omitted with an explicit count and retained in the table. |
| Explainable prediction | An alphabetic chooser can open an unfitted model. AP naming, metric units, SHAP remainder bars and coefficient ranks are easy to confuse. | Complete models with held-out target evidence appear first with readable comparison/status labels. Metric tooltips and numerical interpretation cover ROC/AUC/AP/MCC/balanced accuracy/Brier. Adds pure-block counts and a separate interactive chart of the 30 largest individual validation SHAP contributions, ordered independently of coefficients and without a combined remainder. | No metric confidence intervals, calibration curve, external test cohort or taxonomic sensitivity analysis are calculated. The original static plots remain unchanged; absent validation SHAP summaries are reported explicitly. |
| Protein & Pfam | Wide tables bury protein/feature meaning. A missing positive hit can be mistaken for absence. | Promotes readable identity/feature columns, decodes feature keys and explains Pfam assessment, architecture, direct/inherited memberships and inclusive coordinates. | Domain presence and a common k-mer do not establish catalytic mechanism. Review source authority and sequence completeness. |
| Classes & roles | “Reviewed proteins” can include accepted provisional evidence. Overlapping hierarchy counts should not be drawn as additive parts. Pooled controls are not per-comparison denominators. | Uses **labelled proteins** in the summary and overlapping-membership bars; retains parent IDs and role/mechanistic fields in hover/table help. | Automated accepted labels still need biological review. Counts alone do not assess error rates or control representativeness. |
| Structures & folds | Missing fold names, whole-model components and confident coordinates answer different questions. Connected components are not all-pairs similarity. | Explains eligible/ineligible outcomes, named-fold authority, TM-score, bilateral coverage, RMSD and discovery-frozen projection. | Whole-model clusters have no local motif boundaries. High mean confidence can conceal locally uncertain residues. Imported eligibility follows its source's policy. |
| Model & alignment explorer | Coloured k-mer occurrences can look like independently tested residues; grey coordinates can look like no motif. Pairwise sequence alignment can look like structural superposition. | Explains exact overlapping k-mer mapping, selected feature q-value tiers, overlap rule, white/grey states, sequence/numbering checks and sequence-versus-structure alignment. Adds definitions for alignment columns. | Whole-model clusters cannot supply local intervals. A mapped exact word is not a validated binding/catalytic site. Uploaded annotations remain separate provenance-bearing evidence. |
| Orthology & partitions | Composite group identity and unequal protein/block split sizes need explanation. Missing optional context must not negate memberships. A zero near-redundancy count is not proof of biological uniqueness. | Defines composite HOG authority and count units, plots whole-block allocation, displays **Not assessed** when near-redundancy input was unavailable, and corrects the empty optional-context message. | Connected blocks reduce obvious leakage but do not replace phylogenetic comparative modelling. Optional group context can remain unavailable. |
| Canonical data & downloads | A preview export can be confused with a complete dataset; large workbook limits need visibility. | Page methods and every field are defined. Newly generated workbook dictionaries use the same authoritative definitions as the app. | Previously published workbooks retain their original definitions/checksums. Very large canonical tables still require full TSV/Parquet or the explicitly bounded filtered exporter. |
| Data quality & provenance | Checksum success establishes integrity, not scientific validity. Unmatched targets, incomplete assessments and label circularity matter more than an attractive graph. | Adds pooled matched-control coverage, explains its denominator and selection risk, and loads only one bounded audit preview. Defines all visible assessment/curation/outcome fields. | Defining-feature exclusions cannot eliminate correlated evidence or phylogenetic confounding. Matching coverage is descriptive and does not demonstrate exchangeability. |
| Glossary & help | The original small glossary did not cover the scientific tables or status codes. | Merges scientific concepts, every canonical field, displayed summary fields, controlled states, evidence-class combinations and saved-setting names into one literal searchable/downloadable dictionary. Tests fail when a new canonical field/page term lacks an entry. | Unregistered fields from a future external producer receive explicit unknown-meaning help; no semantics are invented. |

## Exact interpretation of “insufficient”

For association testing, `analysis.minimum_target_proteins` and
`analysis.minimum_background_proteins` are misleadingly named: the implementation compares
them with **pure independent block counts**, separately in each partition. A feature test
also needs enough **feature-assessed** blocks in each class. A relevant positive member
makes a block positive; a block is known negative only if all relevant members were
assessed without the feature. Mixed target/background blocks and unknown absences are
excluded. `analysis.minimum_feature_proteins`, in contrast, is a positive-protein feature
vocabulary filter.

For modelling, each discovery class must meet **both**
`explainable_ml.minimum_samples_per_class` (proteins in pure groups) and
`explainable_ml.minimum_groups_per_class` (pure groups). Held-out metrics and permutation
importance require those same two per-class minimums in validation after mixed groups
are excluded. A fitted model with insufficient held-out counts has
`COMPLETE_VALIDATION_UNDERPOWERED`; descriptive protein predictions and SHAP remain
available while held-out metrics are withheld.

The app displays values saved with the result, including non-default limits. Missing,
malformed or non-finite historical settings display **Not recorded**. It never inserts
current defaults into an old run. Meeting these software execution limits is not a power
calculation, an uncertainty estimate or proof that the study has enough independent
examples to support its biological claim.

The campaign's saved FDR threshold is shown separately from the viewer's existing
positive-discovery suggestion and positional-colour rule, which uses `q <= 0.05`.
This update explains that distinction and does not silently change the scientific
decision rule in a completed result.

## Reading prediction metrics and SHAP

| Quantity | Meaning | How to read it |
| --- | --- | --- |
| CV ROC AUC | Mean pure-block ranking score across group-aware discovery CV folds used to select the penalty | 0.5 is chance ranking; 1 is perfect. Selection CV is not the final independent validation result. |
| Validation ROC AUC | Ranking of a randomly chosen target block relative to a background block, counting ties by half | 0.875 corresponds to about 87.5% pair ordering in the evaluated cohort. It does not mean 87.5% accuracy. |
| Average precision (AP) | Recall-weighted precision summary | Compare with the fraction of targets among evaluated pure blocks. No fixed “good” threshold transfers across differently balanced cohorts. It is not trapezoidal PR AUC. |
| MCC | Correlation between thresholded class calls and supplied labels | +1 perfect, 0 no correlation, -1 reversed. 0.604 is not 60.4% accuracy. |
| Balanced accuracy | Mean sensitivity and specificity at probability threshold 0.5 | 0.5 is chance, 1 perfect. Both classes receive equal weight. |
| Brier score | Mean squared probability error | Lower is better; 0 is perfect. The displayed reference for constant prevalence `p` is `p(1-p)`. It is not a complete calibration assessment. |
| Coefficient | Signed fitted effect in log-odds conditional on other inputs | Positive favours TARGET; negative favours BACKGROUND. Correlated features must be interpreted together. |
| `regularisation_strength` | Selected penalty strength lambda | The classifier uses `C = 1/lambda`; larger recorded values mean stronger regularisation. |
| `importance_rank` | Rank by absolute fitted coefficient | It is not a SHAP rank. The displayed coefficient table is limited to 200 features. |
| Permutation importance | Held-out pure-block balanced-accuracy loss when a feature is shuffled | Negative losses can occur. Repeat standard deviation is not a confidence interval. |
| SHAP beeswarm | Individual explained protein contributions relative to the discovery-block reference | Sign indicates target/background direction. Red/blue indicates feature value, not significance. Read the explained partition. |
| SHAP global bar | Mean absolute contribution in model log-odds | Magnitude has no direction. A grouped remainder combines many features and is not one motif. |
| SHAP waterfall | Additive movement from expected log-odds to this protein's log-odds | Contributions are neither probabilities nor causal/physical effects. |

The source uses the full pure discovery-block reference for its linear SHAP masker.
The table's `mean_absolute_validation_contribution` describes validation proteins;
it is absent without validation. Protein-level explanations can differ from aggregated
block-level predictions and metrics even when both are correctly computed.

The official scikit-learn definitions were checked against the implementation:

- [ROC AUC](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html)
- [Average precision](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html)
- [Matthews correlation](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html)

## Specific screenshot observations

- The data-quality screenshot shows **284 of 459 requested substrate-receptor target
  units** with a matched control, approximately **61.9%** coverage; 175 requested units
  lack a matched control in that pool. The displayed 476 control units are pooled
  matching records, not the assessed denominator of each individual feature test.
  This is a material selection concern, not a numerical defect that tooltips can cure.
  Compare matched and unmatched targets and inspect comparison-specific pure-block
  counts before generalising results to the full target family.
- The orthology screenshot contains **121,423 memberships**, **2,412 groups** and
  **2,623 partition blocks** despite empty optional group context. Group-context absence
  must not be reported as absence of orthology. The actual split shown has different
  protein and block proportions; a block-based allocation need not produce an exact
  protein-level fraction in every class.
- The model screenshot shows **41 localisable feature occurrences covering 94 of
  367 residues**. These are mapped occurrences of selected enriched features, not
  41 independent motif discoveries or 94 independently significant residues. Overlap
  and repeated k-mer occurrences make this distinction essential.
- Long `SC_...` IDs in the coefficient and SHAP plots are stable component keys,
  not functional annotations. The help and added feature-explanation column describe
  their whole-model meaning while preserving exact IDs in downloads.
- Very small q-values and near-one probabilities need appropriate display precision.
  The patch uses three significant figures for p/q/E-values and six decimal places
  for protein probabilities. Six decimals can still round values extremely near one;
  downloads preserve the underlying published numerical precision.
- Large blank table areas and horizontally crowded identifiers reduce readability.
  The patch bounds display height by preview row count and places readable feature,
  protein, status and statistical fields before provenance columns. Exact field names
  remain available in heading help and downloads.

## Delivered architecture and verification

The table below records the first update verification against `939f53b`.
The compatible update verification is recorded separately below.

`protein_signatures.result_help` is the shared field dictionary, available without app
dependencies. `protein_signature_app.help_content` holds pure, testable scientific text,
recorded-limit validation, status/feature resolution and numeric metric interpretation.
The app's common table and plot renderers attach help consistently, including canonical
previews and static SHAP images. Export dictionaries use the same shared definitions.

New functions use Google-style docstrings, defensive input handling, named arguments and
repository-compatible formatting. Unknown producer fields are logged at debug level and
explained explicitly. No scientific CSV files, shell-embedded Python, pipeline threshold
changes, result-schema changes or writes to completed results are introduced.

| Check | Result |
| --- | --- |
| Clean unmodified base | **430 tests passed**; total statement-and-branch coverage **89.09%**, below the existing 95% gate |
| First patched full suite on `939f53b` | **548 tests passed**; total statement-and-branch coverage **89.44%**, below the unchanged 95% gate |
| Affected-page/dictionary/export cases included in the full suite | **147 passing cases**, including headless navigation through all eleven pages |
| Both new help modules | **100% statement-and-branch coverage** |
| Dictionary/export and new chart tests in the separate patched clean clone | **111 passed** |
| Update command syntax and documentation links | **Passed**; the Mac/cluster commands were syntax-checked, not executed on the user's systems |
| Compilation, PEP8, Google docstrings, Ruff format/lint, patch check | **Passed** |

The updated suite adds **118 test cases** compared with the clean base. These cover
missing/invalid settings, non-default limits, unknown fields/statuses, every canonical
column and controlled state, numerical metric interpretation, table tooltips, pure-block
reference counts, duplicate/mixed blocks, bounded audit selection, empty coverage,
q-value underflow/invalid coordinates and both axis choices, plus bounded individual
SHAP ranking independent of the coefficient preview.
The complete logs are included under `verification/` in the downloadable bundle.

Python compilation, PEP8 checking at the repository's 100-character line limit,
Google-convention docstring checks, Ruff formatting and Ruff lint passed. The patch was
successfully checked and applied to a separate clean clone of the base commit; the
dictionary/export tests also passed there. The 95% project coverage requirement is
preserved. Its pre-existing shortfall is a release-readiness issue and is not represented
as a passing quality gate. Closing the uncovered legacy workflow branches requires
separate meaningful tests; lowering the threshold would conceal the problem.

## Compatible update verification against `41f8ffe`

| Check | Result |
| --- | --- |
| Full combined suite | **556 tests passed** in 491 seconds; all eleven app pages opened in headless checks |
| Full-suite project coverage | **88.85%**; the unchanged 95% threshold fails, so the full quality gate is not passed |
| Both new help modules | **100% statement-and-branch coverage** |
| Focused final checks | **153 passed**, with two long application checks excluded and covered by the full suite |
| Independent installer-applied checkout | **140 tests passed**; source path explicitly pointed at that checkout |
| Installer on a matching clean main | **Passed**; all 12 overlay files matched, and pre-update source/Git-history backup was verified |
| Installer on a mismatched base | **Refused** with both commit IDs and no overlay copied |
| Installer with existing local changes | **Refused**; the existing untracked file was retained byte for byte |
| Compilation, PEP8, Google docstrings, Ruff format/lint | **Passed** |
| Mac/cluster command blocks | All **12** Bash blocks syntax-checked; not executed on the user's Mac or cluster |

Compared with the first update, the pushed GitHub commit adds substantial viewer
code and further test cases. Coverage totals therefore describe different source
trees; the earlier 89.44% result must not be used for this combined version.
The repository's 95% requirement is retained. No tests failed, and no coverage gate
has been lowered. Closing uncovered workflow and viewer branches remains follow-up
work. The full user campaign was not available for a live Mac/HPC check.

## Follow-up priorities beyond this patch

1. Validate the provisional labels and unmatched-target selection before interpreting
   apparent predictive performance as protein function. Report comparison-specific
   assessed pure-block counts, not pooled control totals.
2. Add uncertainty intervals for held-out metrics and sensitivity analyses by taxonomy,
   matching coverage and assessment availability. Use a new untouched evaluation set
   after iterative redesign influenced by these validation results.
3. Benchmark rerun time and memory on the full campaign. Bounded output does not eliminate
   the database scan/sort needed by some aggregations and ordered previews.
4. The new interactive individual-feature SHAP chart improves labels and ranking without
   regenerating archived plots. A future scientific run can also improve the original
   static figure layout; existing result checksums and plot bytes remain intact.
5. Develop residue-level structural motif/pocket methods as a separate validated stage
   with explicit residue pairs, assessment universes, provenance and multiplicity policy.
   Whole-model connectivity and mapped exact words cannot substitute for that analysis.
