# Output contracts

A result is valid only when `COMPLETED.json` points to the exact checksum of
`manifest.json`. The manifest records every input authority and every published file with
its size and SHA-256. `protein-signature-app` verifies the result before opening it.

## Storage forms

Every canonical table is written twice through bounded 50,000-row batches:

- `tables/<name>.tsv` for manageable portable tables, or
  `tables/<name>.tsv.gz` when the table exceeds one million rows; and
- `tables/<name>.parquet/part-*.parquet` with the fixed Arrow schema for typed
  analysis. Each deterministic part contains at most one million rows.

`protein_signatures.duckdb` contains physical copies of all tables and a
`signature_evidence` convenience view. Null numeric values remain null in Parquet/DuckDB and
blank in TSV.

## Numbered human reports

The pipeline also writes file-first reports that do not require the app:

| Directory | Contents |
|---|---|
| `analysis/00_run_information` | Reading guide and TSV/XLSX report inventory |
| `analysis/01_proteins_and_curation` | Protein, label and comparison tables; curation figures |
| `analysis/02_homology_and_partitions` | Redundancy, OrthoFinder and frozen-partition tables/figures |
| `analysis/03_sequence_and_domains` | Feature, domain and extracted-sequence tables/figures |
| `analysis/04_structures_and_folds` | Models, acquisition, pairwise alignment and cluster tables/figures |
| `analysis/05_association_statistics` | Complete hypothesis ledger, signature summary and comparison figures |
| `analysis/06_explainable_models` | Model, prediction, importance, numeric SHAP and graphical SHAP outputs |
| `analysis/99_final_results` | Copies of principal decision tables and every decision-facing figure |

Each manageable table in a numbered stage is available as UTF-8 TSV and a formatted XLSX
workbook. A table exceeding 500,000 rows retains its complete canonical Parquet and
TSV/TSV.GZ files and receives a compact TSV/XLSX summary index instead. This avoids treating
Excel's worksheet boundary as a scalable scientific-data format. The application creates
bounded filtered feature exports as TSV and XLSX when a user selects at least one filter.

Workbooks have frozen headers, filters, widths, numeric formats and status highlighting;
Text longer than Excel's cell limit is replaced by a pointer and preserved in numbered
`Long_Text` continuation sheets. TSV/Parquet remain the preferred lossless programmatic
authorities.

Each static figure is published as PNG for quick inspection, SVG for editing and vector PDF
for reporting. The same source file may be copied to an owning stage and `99_final_results`;
both result paths are checksum-inventoried. `report_inventory.tsv` is the entry point.

The application's canonical-data browser covers every table listed below. It shows a bounded
DuckDB preview and the exact result-relative Parquet plus TSV/TSV.GZ storage paths. Complete
formatted XLSX is offered for manageable tables; large tables expose their summary workbook
and bounded filtered exporter. The report inventory itself is downloadable in both formats.

## Analytical checkpoint and resume

Before human reporting begins, all canonical tables and non-tabular analysis state are sealed
under `.protein_signature_cache/analysis_checkpoints/<run-identity>/`. The marker binds the
manifest, input authorities, row counts, table schemas, every Parquet part and TSV/TSV.GZ
file by SHA-256. It is published only after every table finishes. `--resume` accepts it only when the
package/configuration/profile identity and all original input checksums still match.

This checkpoint is computational state rather than a completed result. It permits report,
HTML, DuckDB and final-publication recovery without rerunning feature derivation, 73 E3
comparisons or explainable models. Only the final `result/COMPLETED.json` and result manifest
make a portable campaign complete.

## Canonical tables

| Table | Unit and purpose |
|---|---|
| `proteins` | One authoritative sequence and SHA-256 per protein |
| `redundancy_clusters` | Exact and supplied near-redundancy memberships |
| `profile_labels` | Controlled hierarchy plus default-analysis, background, role and exclusivity policies |
| `label_assignments` | Direct curation authority, including ineligible states |
| `label_memberships` | Reviewed-positive direct and inherited memberships |
| `comparisons` | Exact target/background hypotheses |
| `features` | Unified positive feature occurrences used by inference |
| `feature_assessments` | Complete imported feature hit/no-hit/not-assessed/failed ledger with derivation provenance |
| `domain_hits` | Coordinate-resolved Pfam/other domain hits |
| `domain_assessments` | Complete hit/no-hit/not-assessed/failed coverage matrix |
| `domain_sequences` | Extracted domain subsequences with original coordinates |
| `structures` | Portable structure inventory, eligibility, confidence and versioned fold annotation provenance |
| `alphafold_acquisitions` | Success and failure outcomes for every requested model |
| `structure_comparisons` | Imported or generated pairwise alignment evidence |
| `structure_clusters` | Discovery-defined components and validation projections |
| `imported_structural_group_summaries` | Stable subset of predecessor within-group structural summaries |
| `orthofinder_memberships` | Composite-key HOG/orthogroup memberships |
| `orthofinder_group_context` | Imported precursor group and evolutionary-distance context |
| `partitions` | One discovery/validation block assignment per protein |
| `associations` | Complete per-partition contingency tests and explicit non-test statuses |
| `signatures` | Discovery candidates paired with held-out evidence classes |
| `ml_models` | One fitted or explicit non-fitted status per comparison |
| `ml_feature_importance` | Coefficients, prevalence and held-out permutation importance |
| `ml_predictions` | Discovery and validation probabilities and classes |
| `ml_explanations` | Ranked local SHAP contributions for explained samples |
| `ml_plot_inventory` | Every native SHAP PNG/SVG/PDF asset and its scientific scope |

Column order and physical types are defined in `src/protein_signatures/schemas.py` and are
covered by publication tests.

## Association fields

`target_protein_count` and `background_protein_count` describe raw labelled proteins.
`target_unit_count` and `background_unit_count` are all retained pure independence blocks.
The feature-specific inferential denominators are `target_assessed_unit_count` and
`background_assessed_unit_count`; `target_unknown_unit_count` and
`background_unknown_unit_count` expose the blocks excluded because absence was not
assessed. Both protein-level feature counts and block-level feature counts are retained to
make the reduction transparent. Unknown values are never encoded as feature absence.

Prevalence intervals are Wilson 95% intervals over blocks. `p_value` is the two-sided Fisher
exact result. `q_value` is local to comparison/partition/feature type; `study_q_value` spans
comparisons within partition/feature type.

Feature-specific rows that could not be tested because their assessed denominator is too
small retain the feature identity, assessed/unknown counts, a controlled status and null
inferential values. Campaign-level non-test rows use `feature_type=ANALYSIS`.

## Signature evidence classes

Every biological evidence class is the combination of one scope prefix and one validation
suffix.

| Scope prefix | Interpretation |
|---|---|
| `DECISION_CANDIDATE__` | Biological discovery passes study-wide FDR |
| `EXPLORATORY_LOCAL_ONLY__` | Biological discovery passes local but not study-wide FDR |
| `EXPLORATORY_ALL_DATA_DERIVATION__` | Feature definition used all campaign data and cannot be confirmatory |
| `QC_TECHNICAL_NON_BIOLOGICAL__` | Technical availability association; not a biological signature |

| Validation suffix | Interpretation |
|---|---|
| `VALIDATED_STUDY_WIDE` | Same direction and validation q-values pass both local and study-wide FDR |
| `VALIDATED_WITHIN_COMPARISON` | Validation passes local FDR but not study-wide FDR |
| `DISCOVERY_ONLY` | Discovery candidate does not pass validation FDR |
| `DIRECTION_DISCORDANT` | Validation effect reverses direction |
| `NO_VALIDATION_DATA` | No complete matching held-out test exists |

When completed discovery tests produce no local-FDR hit, the summary evidence class is
`NO_DISCOVERY_SIGNATURE`. If no candidate feature was eligible for testing, the status is
`NO_ELIGIBLE_FEATURES` and the evidence class is
`DISCOVERY_NOT_TESTED__NO_ELIGIBLE_FEATURES`; these outcomes are not interchangeable.

The signature table is a prioritised view, not the complete hypothesis ledger. Use
`associations` for all tests. Discovery rows retain the complete eligible hypothesis ledger.
Validation rows cover the frozen locally significant discovery candidates, including
candidate features with zero held-out positives; discovery selection never consults
validation prevalence or labels.

## Model and SHAP fields

`ml_models.status` begins with `COMPLETE` only for a fitted discovery model. Other statuses
explain whether samples, independent groups, variable features or valid cross-validation
folds were missing. A convergence warning is appended without hiding the model.

`shap_background_partition` is always `DISCOVERY`. `shap_explained_partition` is validation
when any held-out samples exist, otherwise discovery. The expected value and contributions
are model log odds from the interventional linear explainer.

Global plot rows have a blank `protein_id`; waterfall rows carry the exact explained protein.
Assets use portable result-relative paths below
`analysis/06_explainable_models/figures/shap/`. PNG is intended for quick viewing, SVG for
editing and PDF for reports. Every form is checksummed and mirrored into the final-results
figure tree.

## Structural cluster fields

`reference_member_count` counts discovery members used to define the component.
`member_count` includes projected validation members. `membership_method` is
`DISCOVERY_COMPONENT`, `VALIDATION_PROJECTION` or `ALL_DATA_COMPONENT` when no partition map
was supplied to the pure helper. Edge count is discovery-definition evidence; each member
also retains supporting-edge count and best TM score.

## Metadata and manifests

`run_metadata.json` records normalised configuration, package/profile versions, evidence
availability, upstream resource identity, Foldseek cache identity, SHAP/scikit-learn/
Matplotlib versions, table counts and deterministic-policy fields. Profile metadata includes
`require_structural_evidence` and the complete default-comparison policy.
It also records the human-report file, workbook, logical-figure and inventory counts.

`manifest.json` records inputs using absolute execution-time paths. Those paths are used only
for computational resume verification; the portable result can be moved and opened without
the original inputs. `COMPLETED.json` is deliberately tiny so a scheduler or downstream job
can check completion without guessing from directory existence.

## Compatibility policy

Adding a nullable table/column requires a documented schema release. Renaming, changing a
type or changing scientific semantics requires a schema-version increment. Consumers should
select columns explicitly and must not infer completion from row count alone.
