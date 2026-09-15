# Evidence-led labels and matched controls

## Purpose and status

The evidence-led workflow turns heterogeneous annotation authorities into conservative,
auditable protein-type proposals and prespecified matched backgrounds. It is generic: a
profile defines the biological vocabulary and comparisons, while a separate ruleset defines
what evidence is sufficient for each label.

The resulting assignments are marked `EVIDENCE_SUPPORTED_POSITIVE`, not
`REVIEWED_POSITIVE`. They are suitable for provisional hypothesis generation after an
operator explicitly accepts that scope. They are never represented as human-curated or
experimentally proven classifications.

The synthetic `--automated-test-labels` route is different. It creates random, balanced
software fixtures and prohibits biological interpretation. Never use synthetic smoke output
as a scientific authority.

## Workflow overview

1. Load and checksum the FASTA, profile, rules and optional evidence authorities.
2. Build exact-protein contexts for annotations, domains, structures, species and candidate
   state.
3. Form independence blocks from OrthoFinder groups, supplied redundancy clusters and exact
   sequence identity.
4. Score each configured label independently.
5. Resolve descendants, incompatible alternatives, score margins and explicit abstentions.
6. Optionally propagate one generation within a HOG when accepted anchors are unanimous and
   the recipient has the required domain architecture.
7. Select clean, outcome-blind controls at the independence-block level.
8. remove label-defining domains from the confirmatory domain/ML authority.
9. Publish TSV and formatted XLSX audits plus one checksum-complete marker.
10. Stop for inspection unless provisional evidence is explicitly accepted.
11. Run the standard sequence, domain, fold, structural, statistical and SHAP workflow.

## Required generic authorities

### FASTA

`proteins.faa` is the identifier and sequence authority. Identifiers must be unique, stable
and valid under the package identifier contract. The package does not silently rewrite an
identifier.

### Classification profile

The profile YAML defines:

- canonical target, parent and control labels;
- aliases and display names;
- mechanistic class and component role;
- which labels receive default comparisons;
- the background label for each target; and
- whether the downstream campaign requires structural evidence.

Use `configs/profile.example.yaml` with `configs/evidence_rules.example.yaml` as the paired
generic example, and the built-in `e3` profile as a large production vocabulary. A profile
is not an annotation ruleset.

### Evidence rules

The rules YAML is validated against `configs/evidence_rules.schema.json`. It contains one
rule per automatically inferable target label and versioned global policy settings.

Each label rule can specify:

- ordinary annotation regular expressions;
- highly specific annotation regular expressions;
- exclusion expressions that veto misleading text;
- required Pfam/domain expressions;
- supporting domain expressions;
- whether one trusted highly specific annotation may stand alone; and
- whether compatible HOG propagation is allowed.

All patterns are case-insensitive. Rules are exact scientific authorities: change the
ruleset version whenever a pattern, score, threshold or matching policy changes.

## Optional input contracts

### Generic protein metadata

`protein_metadata.tsv` has exact columns:

| Column | Meaning |
|---|---|
| `protein_id` | FASTA identifier |
| `species` | Species label; multiple values use `|` |
| `input_candidate` | Strict `TRUE` or `FALSE` |

The file must cover every FASTA protein. Species is used only as a control-matching
covariate. `input_candidate` prevents upstream candidates from silently becoming negative
controls when the ruleset requests that exclusion.

### External annotations

`external_annotations.tsv` has exact columns:

| Column | Meaning |
|---|---|
| `protein_id` | FASTA identifier |
| `label_id` | Optional exact profile label |
| `annotation_text` | Optional text searched by rules |
| `evidence_status` | Authority state, such as `CURATED` |
| `evidence_source` | Named database, release or curation authority |
| `evidence_reference` | Record, release, URL or checksum reference |
| `annotation_scope` | `DIRECT`/`EXACT_PROTEIN` or a weaker declared scope |

An exact label becomes trusted direct evidence only when its status is one of `CURATED`,
`EXPERIMENTAL`, `MANUALLY_REVIEWED` or `REVIEWED`, and its scope is `DIRECT` or
`EXACT_PROTEIN`. Text rows from the same named source form one evidence group, so duplicate
records from one database cannot masquerade as independent corroboration.

### Seed assignments and catalogues

`seed_assignments.tsv` uses the standard label-assignment schema. Only rows already marked
`REVIEWED_POSITIVE` become trusted direct anchors.

The E3 seed catalogue may provide exact seed text and weaker cluster-associated context.
Cluster-associated text is deliberately down-weighted and cannot, with Pfam alone, reach the
built-in E3 acceptance threshold. The catalogue never converts a nearby protein into a
reviewed biochemical truth.

### Domains

`domains.tsv` uses the standard explicit assessment schema. `ASSESSED_WITH_HIT`,
`ASSESSED_NO_HIT`, `NOT_ASSESSED` and `FAILED` remain distinct. A domain by itself never
establishes the protein type. Required domains corroborate annotation or HOG evidence.

### Structures

`structures.tsv` uses the standard model/fold inventory. The labeller reads only controlled
analysis eligibility and mean confidence. It does not use an observed fold, structural
cluster, pocket or similarity score to define a target. This separation lets downstream
structural analysis test for association rather than restate the label rule.

### OrthoFinder

Use either a published `orthofinder-results` resource or a completed raw OrthoFinder 2.5.5/3
directory. The package consumes OrthoFinder but never executes it. HOG mode requires a
hierarchy node, normally `N0`; legacy orthogroup mode requires a blank hierarchy node.

### Redundancy clusters

An optional near-redundancy TSV joins proteins into shared analysis blocks. These blocks are
combined transitively with exact-sequence and OrthoFinder links, preventing related copies
from crossing target/control or discovery/validation boundaries.

## Label decision policy

### Independent evidence groups

Scores are accumulated from independent groups, not raw record count. Typical built-in E3
groups are annotation context, Pfam and OrthoFinder homology. Within one annotation source,
only the strongest matching row contributes.

The default E3 scores are:

| Evidence | Score |
|---|---:|
| Ordinary annotation | 2.0 |
| Specific annotation | 4.0 |
| Required Pfam evidence | 2.0 |
| Supporting Pfam evidence | 0.5 |
| Compatible unanimous HOG propagation | 2.0 |
| Trusted exact label | 5.0 |

Normal rule-based acceptance requires a score of at least 3.5, at least two independent
groups, all required domain expressions and a qualifying structure. A trusted direct label
still has to meet the selected structure-eligibility policy for entry into the structural
campaign. These defaults live in the versioned YAML and are not hidden constants.

### Specific annotation exception

A ruleset may explicitly allow one trusted, exact-protein, highly specific annotation to
stand alone. The built-in E3 rules do not currently enable that shortcut. Ordinary annotation
alone is never accepted.

### Hierarchy and conflicts

If an accepted descendant and its accepted ancestor both occur, the descendant wins. Two
non-ancestral alternatives are resolved only when neither is a trusted direct conflict and
the leading score exceeds the configured margin. Otherwise the protein is `AMBIGUOUS` and is
excluded from both target and control cohorts.

### Abstention states

- `ACCEPTED`: policy, evidence and structure requirements passed.
- `PROPOSED`: target-like evidence exists but is below policy or incomplete.
- `AMBIGUOUS`: incompatible accepted alternatives remain unresolved.
- `UNMAPPED`: no configured target evidence was observed.

Proposed and ambiguous rows are retained in the review queue but are not positive analysis
memberships.

## Conservative HOG propagation

Propagation is optional and only one generation deep. A recipient must satisfy all of the
following:

1. belong to the same selected OrthoFinder group as at least the configured number of
   accepted anchor proteins;
2. have anchors whose labels resolve to one compatible rule-bearing label or ancestor;
3. use a rule that explicitly permits propagation;
4. contain every required domain pattern for that rule; and
5. satisfy the target structure policy.

Propagated calls do not recursively become new anchors during the same run. The anchor list
is checksummed in the evidence reference.

## Background selection

### Why a single housekeeping family is unsuitable

A narrow housekeeping family has its own conserved fold, domain architecture, length range,
subcellular biology and evolutionary history. Using it as the universal negative would make
those properties look like target signatures. The package therefore selects controls from a
broad clean pool and matches each target block on prespecified technical and biological
covariates.

### Clean-control exclusions

A protein cannot enter the control pool if:

- it has any accepted, proposed or ambiguous target-rule evidence;
- it belongs to any target homology/redundancy block;
- it is an upstream input candidate when that exclusion is enabled; or
- it fails a configured matching caliper.

These exclusions are outcome-blind: no k-mer, fold enrichment, structural cluster or model
prediction is inspected during control construction.

### Matching unit and variables

One representative is selected deterministically per independence block, favouring an
eligible structure, higher confidence, shorter sequence and then lexical identifier. Each
target block receives the configured number of unused control blocks. Default E3 matching
requires:

- an intersecting species label;
- the same structural-eligibility state;
- absolute log2 length difference no greater than 0.75; and
- domain-count difference no greater than 2.

Among candidates passing every caliper, the closest is selected by a deterministic distance
combining length, domain-count difference, domain-architecture Jaccard dissimilarity and
scaled mean-confidence difference. The audit contains every successful or unmatched request.

No target class is analysed unless its profile-resolved background has the configured minimum
number of independent control blocks.

## Circularity protection

Pfam features that contributed to any accepted automated label are listed in
`label_definition_features.tsv` with scope
`GLOBAL_CONFIRMATORY_DOMAIN_AND_ML`. Those accessions are removed from the domain file passed
to downstream association testing and SHAP modelling. Assessment coverage is retained with
explicit no-hit rows after removal.

This is intentionally global within a campaign. A domain used to label one E3 class cannot
be advertised as a novel confirmatory signature for another class in that same automated
campaign. The unfiltered source remains checksummed in the evidence bundle for audit.

This protection does not make the whole study causal. Annotation, orthology and control
selection can still correlate with phylogeny, taxon sampling and data availability.

## Evidence-bundle outputs

Every table is emitted as machine-oriented TSV and formatted XLSX:

| File stem | Purpose |
|---|---|
| `label_assignments` | Complete downstream assignment authority |
| `label_evidence_audit` | Every scored candidate, item, decision and conflict |
| `control_matching_audit` | Target/control block matches and distances |
| `label_definition_features` | Features excluded to prevent circular inference |
| `class_labelling_summary` | Per-class target/control coverage |
| `unresolved_assignments` | Proposed and ambiguous review queue |
| `domains.for_signature_analysis` | Circularity-safe domain authority, when supplied |

`EVIDENCE_LABELS.json` records package/profile/rules versions, settings, counts, input
checksums and every output checksum. Any modification invalidates verification.

## Generic execution

First create and inspect only the evidence bundle:

```bash
./run_evidence_signature_workflow.sh \
  --work-dir /absolute/persistent/campaign \
  --campaign-id example_family_20260915 \
  --sequences-fasta /absolute/inputs/proteins.faa \
  --profile /absolute/inputs/profile.yaml \
  --evidence-rules /absolute/inputs/evidence_rules.yaml \
  --protein-metadata /absolute/inputs/protein_metadata.tsv \
  --domains /absolute/inputs/domains.tsv \
  --structures /absolute/inputs/structures.tsv \
  --external-annotations /absolute/inputs/external_annotations.tsv \
  --orthofinder-results /absolute/orthofinder/Results \
  --threads 24
```

After audit, rerun the same command with:

```text
--accept-provisional-evidence-labels
```

For an outer Slurm allocation, also add:

```text
--submit-slurm --slurm-account barton --slurm-partition barton \
--slurm-memory 128G --slurm-time 2-00:00:00
```

The outer job uses the local Snakemake profile inside its allocation. Logs and workflow state
remain under the campaign work directory. The immutable result is `result/`.

## Completed E3 precursor execution

Use `run_completed_e3_workflow.sh --evidence-led-labels`. Without explicit provisional
acceptance, its DAG stops after the evidence bundle. With acceptance, it runs preparation,
evidence proposal, approval-marker generation, campaign initialisation, analysis and final
verification unattended.

The E3 rules cover 67 inferable default targets, including F-box. Six catch-all labels ending
in `other_reviewed` intentionally have no inference rule; only trusted direct assignments can
populate them. The profile still defines all 73 comparisons.

## Review priorities

Before treating the result as more than provisional, review:

1. the highest-population and highest-impact target classes;
2. every ambiguous record and near-threshold proposal;
3. target classes with few independent units;
4. unmatched target rows or controls near a caliper;
5. whether taxon and structure-availability balance is scientifically plausible;
6. every HOG-propagated call and its anchors;
7. the complete label-defining-feature exclusion ledger; and
8. sensitivity to alternate ruleset and matching-policy versions.

Human review should create a new, separately versioned assignment authority. It must not edit
a published bundle in place.
