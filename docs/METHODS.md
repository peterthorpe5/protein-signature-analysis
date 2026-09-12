# Methods and statistical interpretation

## Analysis question

Each comparison asks whether a pre-specified protein type carries a categorical sequence,
domain, fold or structural feature at a different rate from an explicit background. The
comparison is defined before feature testing by target and background label IDs. It is not
created by selecting whichever proteins carry a promising feature.

The workflow reports two related quantities:

- **commonness**, the prevalence of a feature among independent target blocks; and
- **association**, the difference between target and background prevalence, together with
  a hypothesis test and false-discovery-rate correction.

A feature can be common but uninformative when it is equally common in controls. A rare
feature can be strongly associated but have limited sensitivity. Both prevalence and effect
direction must therefore accompany every q-value.

## Evidence units

The FASTA protein is the join key, but the default inferential unit is an independence block.
An undirected graph joins proteins that share any supplied OrthoFinder group, exact sequence
cluster or supplied near-redundancy cluster. The connected components of that union graph
are indivisible blocks. Transitive closure matters: if A shares a HOG with B and B is nearly
identical to C, A, B and C remain together even if A and C have no direct relationship.

Each block is assigned once to discovery or validation by a stable hash of its sorted
membership, the configured random seed and the validation fraction. This prevents an exact
duplicate, close homologue or connected paralogue block from contributing to both stages.

For a given comparison, a block containing at least one target and at least one background
protein is ambiguous. It is excluded from both cells of that comparison and reported as an
`excluded_mixed_unit_count`. It is never counted twice or resolved after looking at the
feature.

## Feature construction

All positive evidence is represented by a stable `(feature_type, feature_id)` pair. Presence
is binary at protein level and then binary at independence-block level: ten occurrences in
one related block still contribute one vote.

### Amino-acid sequence

Configured k-mers are enumerated as exact, overlapping uppercase amino-acid strings. A k-mer
is retained only when it occurs in at least `minimum_feature_proteins` proteins assigned to
the frozen discovery partition, subject to the deterministic `maximum_kmer_features` limit
over discovery candidates. The retained vocabulary is then scanned unchanged across every
FASTA protein, including held-out validation proteins. Each feature records the exact
discovery sequence-cohort digest and a deterministic definition digest. The feature states
that a string is present, not how many times it occurs.

This deliberately simple native layer is useful for short exact signals and compositional
checks. It is not equivalent to a gapped motif, position-specific probability matrix,
profile HMM, evolutionary-rate estimate or experimentally verified active site. Calls from
MEME, STREME, FIMO, HMMER, conservation or other reviewed producers can be imported through
`features.tsv` with their exact definition and derivation provenance. `DISCOVERY_DERIVED`
is a submitter attestation that feature definition was label-blind as well as restricted to
the exact discovery sequence cohort; the software verifies the cohort digest but cannot
inspect an external producer's behaviour. Label-aware discovery definitions use
`DISCOVERY_SUPERVISED` and remain audit-only in v0.1 because reusing their selection labels
for a discovery Fisher test would produce post-selection p-values. Features learned from the
full campaign are likewise retained as exploratory evidence but cannot enter confirmatory
held-out association or explainable modelling, and can never be promoted to a decision
candidate. The v0.1 production pipeline keeps these rows audit-ledger-only so they cannot
alter confirmatory feature-family FDR. The low-level association API can label all-data
features in a separately scoped discovery-only exploratory analysis when a complete
assessment universe is available; callers must not merge that exploratory test family into
confirmatory FDR.

### Domains and architectures

Coordinate-resolved domain hits create one domain-presence feature per authority/accession.
The ordered accession series creates a separate architecture feature. Domain subsequences
are extracted using the original one-based inclusive coordinates. Pfam assessment states
remain explicit, so `ASSESSED_NO_HIT`, `NOT_ASSESSED` and `FAILED` cannot collapse into the
same biological absence.

### Folds and structure availability

Explicit `ASSESSED_WITH_HIT` fold IDs from an eligible model or reviewed external-evidence
row become fold features. Only a physically present, explicitly analysis-eligible coordinate
creates a structure-availability feature. Source-specific assessment universes distinguish a
known lack of a qualifying model from missing or failed evidence. These technical availability
features are excluded from the default classifier.

AlphaFold DB acquisition is a model-retrieval adapter, not local structure prediction. A
model must map through an explicit protein-to-UniProt table. The coordinate checksum,
reported model version, sequence agreement when available and mean PDB B-factor/pLDDT are
recorded. Exact sequence verification and the configured mean-pLDDT threshold jointly control
whether the acquired structure is usable; excluded downloads remain auditable and cannot enter
Foldseek or structural feature derivation. This does not prove that a low-confidence region is
disordered. `MODEL_NOT_AVAILABLE`, `SEQUENCE_MISMATCH` and completed acquisition exclusions
are assessed negatives for an analysis-eligible AlphaFold model; `FAILED` and `NOT_SELECTED`
remain unknown.

### Alignment-derived structural clusters

A pairwise structural edge passes only when its status is complete, its TM score is at least
`structural_tm_score_threshold`, and coverage for **both** proteins is at least
`structural_minimum_coverage`. Discovery-to-discovery passing edges define connected
components. Those components are frozen. A validation protein may be projected into a
component through a passing edge to a discovery member, but validation-to-validation edges
cannot create or merge a discovery signature. Connected components are calculated separately
for every comparison-universe, coverage-scope, tool and version tuple. A cluster absence is
assessed only when complete universe membership is declared; sparse retained hits do not
define negatives.

Imported precursor pocket positives are therefore published in the assessment ledger but
are not association-test candidates: the precursor does not prove a pocket-specific negative
universe, and its all-data derivation is exploratory. A future exploratory pocket analysis
must declare its own assessment and multiple-testing scope.

Consequently, `STRUCTURAL_CLUSTER` association tests ask whether membership in a
discovery-defined structural neighbourhood is enriched. They do not assert a formal SCOP,
CATH or ECOD fold unless a separate reviewed fold authority was supplied.

Imported predecessor pocket calls remain exploratory `STRUCTURAL_POCKET` feature-assessment
rows, retain the upstream resource manifest and are namespaced by upstream cluster. Pocket
positions from different predecessor clusters are not treated as the same feature. Global
structural clusters and within-group conserved-pocket claims are intentionally different
evidence families; only the former enters v0.1 association testing.

## Prevalence and confidence intervals

For a feature, let:

- \(a\) be target blocks with the feature;
- \(b\) be target blocks without it;
- \(c\) be background blocks with the feature; and
- \(d\) be background blocks without it.

Target and background prevalence are

\[
\hat p_T = \frac{a}{a+b}, \qquad
\hat p_B = \frac{c}{c+d},
\]

and the reported effect is \(\hat p_T-\hat p_B\). Positive values mean target enrichment;
negative values mean background enrichment.

Each prevalence has a two-sided 95% Wilson score interval with \(z=1.9599639845\):

\[
\frac{\hat p + z^2/(2n) \pm
z\sqrt{\hat p(1-\hat p)/n + z^2/(4n^2)}}{1+z^2/n}.
\]

The interval describes uncertainty in block prevalence. It is not a confidence interval for
the prevalence difference or odds ratio.

## Fisher exact association test

The workflow uses a two-sided Fisher exact test on

| | Feature present | Feature absent |
|---|---:|---:|
| Target blocks | \(a\) | \(b\) |
| Background blocks | \(c\) | \(d\) |

Conditioning on the margins, the probability of the observed upper-left cell is
hypergeometric. The two-sided p-value sums every feasible table whose probability is no
greater than the observed table, with a small numerical tolerance. The odds ratio is
\(ad/(bc)\), reported as infinite when the numerator is non-zero and denominator is zero,
and null when it is not identifiable.

Fisher's test is exact conditional on the table margins; it does not compensate for a poor
control set, unmodelled phylogeny, label error, ascertainment bias or availability bias.

## Two Benjamini–Hochberg corrections

For ordered p-values \(p_{(1)} \le \cdots \le p_{(m)}\), Benjamini–Hochberg starts with
\(m p_{(i)}/i\), caps at one and applies the reverse cumulative minimum. Ties are handled by
stable deterministic ordering.

Two correction scopes are retained:

1. `q_value`: all complete tests in one comparison × partition × feature type.
2. `study_q_value`: all complete tests across comparisons in one partition × feature type.

Discovery applies the configured minimum positive-protein filter and local FDR. Features
passing local discovery FDR with a non-zero effect form the frozen validation hypothesis
family. Validation tests every such feature regardless of held-out prevalence, including
features with no validation positives. Thus an assessed absence is recorded as
non-replication rather than mislabelled `NO_VALIDATION_DATA`, without multiplying the entire
high-cardinality discovery vocabulary across the held-out ledger. Because selection uses
discovery data only, held-out Benjamini–Hochberg correction is applied over the frozen
selected family within each feature type.

Feature types are separate pre-declared evidence families. This prevents a very large k-mer
screen from consuming the entire multiplicity budget of a smaller Pfam or structural-fold
family. It also means splitting or renaming a feature type changes the statistical design.

The study-wide value protects a broad survey of many E3 subclasses more strongly than the
local value. It is not a hierarchical ontology-aware procedure, and correlated parent/child
comparisons remain biologically correlated. The full p/q ledger is therefore published; no
single q-value should be interpreted without its pre-specified scope.

## Discovery and held-out synthesis

A discovery candidate requires:

- a complete discovery association;
- local `q_value <= fdr_threshold`; and
- a non-zero prevalence difference.

Before held-out support is appended, each discovery is assigned an inferential scope:

- `DECISION_CANDIDATE__` only when a biological feature passes study-wide FDR;
- `EXPLORATORY_LOCAL_ONLY__` when it passes local but not study-wide FDR;
- `EXPLORATORY_ALL_DATA_DERIVATION__` when its definition used the full campaign; or
- `QC_TECHNICAL_NON_BIOLOGICAL__` for technical availability evidence.

The matching validation row is then classified without refitting the association, and the
following suffix is appended to that scope:

| Evidence class | Rule |
|---|---|
| `VALIDATED_STUDY_WIDE` | Same direction; validation local and study-wide q-values pass |
| `VALIDATED_WITHIN_COMPARISON` | Same direction; validation local q-value passes |
| `DISCOVERY_ONLY` | Direction is not reversed, but validation local FDR does not pass |
| `DIRECTION_DISCORDANT` | Discovery and validation prevalence differences have opposite signs |
| `NO_VALIDATION_DATA` | Matching complete validation test is unavailable |

If completed discovery tests produce no local-FDR hit, the comparison receives an explicit
`NO_DISCOVERY_SIGNATURE` evidence class with status `NO_SIGNIFICANT_SIGNATURE`. If no feature
was eligible for a discovery test, it instead receives status `NO_ELIGIBLE_FEATURES` and
evidence class `DISCOVERY_NOT_TESTED__NO_ELIGIBLE_FEATURES`. “No significant signature” means
assessed under the configured design; it is not equivalent to not run.

## Explainable classification

The predictive layer is separate corroborating analysis. For each comparison it:

1. takes reviewed target/background proteins and their frozen block identities;
2. removes technical feature types by default;
3. uses only discovery proteins to discard constant features;
4. ranks features by discovery document frequency within each feature type and selects them
   round-robin across sorted feature types, without consulting class labels;
5. selects an elastic-net logistic regularisation strength using shuffled
   `StratifiedGroupKFold`, with blocks kept intact and mean ROC AUC as the score;
6. fits a class-balanced `saga` logistic model to discovery data only; and
7. evaluates untouched validation proteins.

The model uses \(C=1/\lambda\), where the configuration lists \(\lambda\) values. A fixed 0.5
decision threshold produces balanced accuracy and Matthews correlation. Validation ROC AUC,
average precision and Brier score are reported only when both validation classes exist.
Permutation importance is the drop in held-out balanced accuracy after independently
permuting a feature, repeated with the configured seed.

This workflow does not claim an unbiased estimate after arbitrary iterative exploration of
the same validation set. Once validation results influence new labels, features or
thresholds, a new untouched test set is required.

## Mandatory SHAP explanations

Every fitted model uses `shap.LinearExplainer` with an independent masker built only from
the discovery matrix. Validation samples are explained when available; otherwise the
discovery samples are explicitly marked as the explained partition. For a linear model and
the interventional background, the contribution is in additive log-odds space and is
equivalent to the coefficient multiplied by displacement from the background expectation.

Published outputs include the numeric top local contributions, a beeswarm, a global bar
plot and bounded per-protein waterfalls in PNG, SVG and PDF. SHAP answers why this fitted
model produced a prediction. It does not show causation, physical contact, sufficiency or a
universal family mechanism. Correlated domains, k-mers and structural features can divide or
duplicate attribution and should be interpreted together.

## Reproducibility

Input authorities and every result file are SHA-256 inventoried. Configuration paths are
resolved before run identity calculation. Feature ordering, partitions, model randomness,
SHAP jitter and file naming are deterministic. Foldseek caches include coordinate content,
effective parameters and tool version. Publication uses a sibling staging directory and one
same-filesystem atomic rename; a directory without a valid completion marker and manifest
is not complete.
