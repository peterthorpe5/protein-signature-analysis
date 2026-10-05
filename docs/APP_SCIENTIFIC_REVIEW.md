# Scientific and presentation review of the first app run

Reviewed on 2–3 October 2026 against the thirteen screenshots from 16:19–16:21,
the additional feature-ID crop at 16:34, and the code that queries the completed
`e3_all1972_caliper_complete_20260928` result. This is a review of the displayed
evidence, not an independent re-analysis of every canonical row.

## Findings that affect interpretation

1. **The current run is chiefly hypothesis generation.** The labels are automated
   evidence-supported proposals without completed human review. Only three of 73
   comparisons have positive discovery results in the displayed run. The others
   must be separated into insufficient sample size, complete without positive
   enrichment, and other published states. A missing significant row does not
   demonstrate absence of a class-specific feature.
2. **Control coverage is a limiting factor.** A previous matched-control summary
   reported 284 of 459 substrate-receptor target units with a control, despite
   476 matched control units. A large background pool does not mean every target
   is represented. The app now separates the *published comparison-specific*
   analysis cohorts from the pooled background match audit. Comparisons that share
   a background can have different included target counts. The matching policy
   and any remaining selection bias still require scientific review.
3. **Short k-mers are sequence words, not established structural motifs.**
   `k3:LPD` means a leucine–proline–aspartate word occurs somewhere in a protein.
   Overlapping short words are correlated and may reflect composition, shared
   ancestry or a label-defining domain. Increasing k to five changes the candidate
   vocabulary but does not by itself establish a local 3D motif. The pending k5
   campaign must be analysed separately after completion.
4. **Structural clusters are whole-model groups.** `SC_...` is a stable fingerprint
   of a connected structural-similarity group, not a named fold, local motif or
   pocket. The screenshot includes groups with hundreds to thousands of models at
   TM-score and coverage thresholds of 0.5. Single-linkage connections can join
   diverse members through intermediary models. Membership enrichment should be
   checked against homology, model coverage and matched controls.
5. **The old pairwise scatter was a biased sample.** It selected the 5,000 pairs
   with the highest TM-score and plotted them as if they described the overall
   score/coverage distribution. The page now bins **all** pairs with both scores
   and coverage; its high-scoring individual-row table is explicitly optional.
6. **Positional colour is a projection, not a residue test.** An enriched k-mer is
   first tested for presence per protein. Its exact sequence occurrences inherit
   that feature's adjusted q-value for the 2D and 3D displays. Individual residues
   were not independently tested. White means no mapped significant feature, not
   verified absence; red reflects a smaller adjusted q-value, not higher effect
   size, pLDDT or a binding pocket.
7. **Missing source assessments matter.** The screenshot lists 49,779 Pfam
   assessments as `NOT_ASSESSED`, about 41% of the 121,635 proteins. Those are
   unknown, not Pfam-negative. Imported `STRUCTURAL_POCKET` assessment records
   have `ALL_DATA_EXPLORATORY` scope; they are not residue-level pocket
   enrichment results from this campaign. The app now calls out both limitations.
8. **A verified file is not a validated biological conclusion.** Checksums show
   that the app read the published result unchanged. They do not establish label
   correctness, adequate matching, calibration, independent validation or
   mechanistic interpretation.

## Page-by-page inspection

| Screenshots | Page and displayed problem | Change in this update | Remaining scientific work |
| --- | --- | --- | --- |
| 16:19:10 | Overview began with raw inventory counts and a feature-coverage bar chart, without showing which classes were analysable. “Candidate signatures” included published records of different states. | Added comparison-outcome census, validation counts, a short class chart, provisional-label warning and “signature records” label. Full feature scan remains on request. | Review each comparison's sample sizes and matched controls before class claims. |
| 16:19:25 and 16:34:36 | Signature explorer led with a wide machine-ID table and an undirected thousands-of-points scatter. `feature_name` repeated `k3:LPD` rather than explaining it. | Added biological explanations beside IDs, a per-family shortlist, explicit discovery/validation tier, effect-size bars, scientific-format q-values and expandable full ledger. | Check correlated words, domain circularity, independent held-out effect and local motif specificity. |
| 16:19:30–16:19:36 | Explainable prediction selected an alphabetically early comparison with `INSUFFICIENT_SAMPLE_SIZE`, zero fitted features and missing metrics. The large raw model row filled the page. | Complete held-out models sort first; selection displays readable class/status; raw model record is expandable and a non-fitted reason is stated. Published SHAP graphics gain interpretation notes. | Check whether any held-out model is sufficiently powered and calibrated; SHAP is not causal evidence. |
| 16:19:50 | Protein/Pfam opened on an alphabetically early, unlabelled protein and an empty labels tab. | Suggests labelled proteins with packaged models first; exact ID search remains available. Evidence sources load separately and feature rows carry explanations. | Audit Pfam coverage and whether selected proteins represent independent units. |
| 16:19:59 | The provisional class summary showed many accepted controls but concealed unmatched target coverage. The small sunburst mixed overlapping ancestors and children. The field `reviewed_proteins` was inaccurate for provisional memberships. | Comparison-specific matched and excluded target units come from the published analysis metadata; pooled background totals are shown separately. The hierarchy chart is a bounded label bar view, and the count is `observed_proteins`. | Human review of target labels, ambiguous roles and overlapping hierarchies. |
| 16:20:09 | Structure eligibility was shown, but no named folds existed; cluster hashes and huge components looked like meaningful local structures. The pairwise scatter sampled only top TM-scores. | Added cluster descriptions and chaining caution; replaced scatter with full-distribution score/coverage bins; separated optional top-pair inspection. Imported pocket/group summaries are marked exploratory. | Develop local residue-corresponded structural features and test them with homology-aware controls. |
| 16:20:27–16:20:44 | Cross-class ranking came before the selected protein, used raw IDs and many long q-values. The default discovery projection coloured many overlapping short words. The 3D model occupied a large mostly empty chart; paired heat blocks omitted amino-acid letters. | Model and mapped regions come first, strongest available evidence tier is the default, cross-class and pair panels load on request, Cα chart is compact, and the pair window has aligned residue letters and gaps. | Confirm exact model numbering, then evaluate whether mapped regions have structural locality and class specificity. |
| 16:20:52 | OrthoFinder memberships existed, yet the empty optional detailed context read like a missing orthogroup run. Zero near-redundancy clusters appeared to mean no redundancy, although the input was unavailable. | Explains the optional context separately, shows partition blocks, and displays “Not assessed” when near-redundancy input was unavailable. | Verify group topology and the separation of related proteins between discovery and validation. |
| 16:20:59 | Canonical page showed the first physical protein rows and very long sequences without a scientific selection. | Kept this as a clearly bounded raw-data/export page with page-specific methods. | Treat previews as storage samples, not ranked evidence; use feature and comparison views for conclusions. |
| 16:21:10 | Quality page opened with low-level rejected label decisions, buried matched-control coverage, showed imported pocket assessments without scope, and left AlphaFold acquisition blank. | Comparison-specific cohort counts and pooled background diagnostics appear before audits; audits load one section at a time; Pfam unknowns and pocket scope are explicit; empty acquisition is explained; full block scan is on request. | Resolve evidence conflicts and review label rules, control matching and missing assessments. |

Each page now has expandable purpose and methods/limitations boxes. Each rendered
graph has an expandable interpretation. Machine IDs are retained in exports so
figures can be traced to the canonical tables; their readable meanings are added
to exploratory tables and the glossary.

## Before claiming a transferable structural signature

- Complete the k5 trial and compare discovery and held-out validation on the
  same precisely defined independent units and control policy. Inspect the
  feature-family shortlist rather than counting correlated k-mers as motifs.
- Review automated labels and target/control coverage, especially substrate
  receptors. Report how many targets lack a match and how missing Pfam scans
  affect the class definition.
- For a genuinely structural motif, require a local, reproducible set of
  aligned residue positions or another well-defined 3D feature, then test its
  enrichment and background prevalence on independent homologous groups and
  other organisms. Whole-model connected clusters alone do not provide this.
- Validate the interface on the Mac with the full copied result: measure first
  open, page switches, model selection, annotation upload, optional pair view,
  TSV/Excel/image downloads and all graph explanations. A screenshot cannot
  verify the underlying numerical rows or responsiveness.
