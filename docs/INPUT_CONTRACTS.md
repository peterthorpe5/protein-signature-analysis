# Input contracts

All scientific tables are UTF-8, tab-delimited and have one header row. Fields may be blank
only where stated. Identifiers must be non-empty, whitespace-free controlled tokens. Paths
in campaign YAML are resolved relative to the YAML file.

The parser, rather than this prose, is the executable authority. Example files are under
`examples/minimal_e3/`.

## Campaign YAML

Schema version 1 has these top-level keys:

| Key | Requirement | Purpose |
|---|---|---|
| `schema_version` | Required; exactly `1` | Configuration contract |
| `campaign` | Required | Stable campaign ID and profile source |
| `inputs` | Required | Explicit file/resource authorities |
| `analysis` | Optional mapping | Feature, partition, FDR and structural thresholds |
| `alphafold` | Optional mapping | AlphaFold DB retrieval policy |
| `foldseek` | Optional mapping | All-versus-all structural search policy |
| `explainable_ml` | Optional mapping, always enabled | Group-aware model and SHAP policy |
| `comparisons` | `profile_defaults` or a non-overlapping list | Target/background hypotheses |

Use [the complete example](../configs/campaign.example.yaml) and
[JSON schema](../configs/schema.json). `protein-signatures validate` is the required
pre-flight check for the configuration and every configured local TSV: FASTA labels,
generic features, domain assessments, redundancy, structures, pairwise comparisons and
AlphaFold request mappings. Selected completed-resource adapters are also verified. When
Foldseek is enabled, pre-flight resolves the executable, records its version and proves that
the candidate count is compatible with `maximum_hits`; it never runs a Foldseek search or
downloads an AlphaFold model.

For production work, use `start_from_inputs.sh --initialise-only` to create and validate a
new `campaign.yaml`, inspect and edit that file to pre-specify all thresholds and comparisons,
then run the reviewed configuration with `start_from_inputs.sh --work-dir ... --resume`.
Once the YAML exists it is the sole configuration authority; config-defining command-line
options are rejected on the resumed invocation.

## Protein FASTA

- Plain text or gzip-compressed.
- The first whitespace-delimited header token is the exact `protein_id`.
- Identifiers must be unique.
- Sequences are upper-cased; a single terminal `*` is removed.
- Internal stops and unsupported residue symbols fail validation.
- No identifier or isoform normalisation is performed.

## `label_assignments.tsv`

Required header:

```text
protein_id	label_id	curation_status	evidence_status	evidence_source	evidence_reference	component_role	curation_reason
```

`label_id` must exist in the chosen profile. Accepted `curation_status` values are:

- `REVIEWED_POSITIVE`
- `REVIEWED_NEGATIVE`
- `REVIEWED_COMPONENT_NOT_CATALYTIC`
- `PROPOSED`
- `AMBIGUOUS`
- `UNMAPPED`
- `EXCLUDED`

Only `REVIEWED_POSITIVE` is expanded into analysis membership. A background control must
therefore be positively reviewed as a member of its background label, not merely marked
`REVIEWED_NEGATIVE` against a target label. `(protein_id, label_id)` must be unique.

## `features.tsv`

This optional generic bridge accepts versioned feature calls from tools not run natively.

```text
protein_id	feature_type	feature_id	feature_name	start	end	evidence_status	evidence_source	evidence_reference	derivation_scope	feature_definition_sha256	derivation_cohort_sha256
```

`start` and `end` are either both blank or one-based inclusive coordinates within the
protein. An exact protein/type/feature/coordinate tuple is unique. This is an assessment
ledger, so a missing protein/feature row is unknown rather than an observed absence.

| `evidence_status` | Inference meaning |
|---|---|
| `ASSESSED_WITH_FEATURE` | Assessed positive; coordinates are permitted |
| `ASSESSED_NO_FEATURE` | Assessed negative; coordinates must be blank |
| `NOT_ASSESSED` | Unknown; coordinates must be blank |
| `FAILED` | Unknown because assessment failed; coordinates must be blank |
| `EXCLUDED` | Deliberately excluded and unknown; coordinates must be blank |

Every definition must also declare one leakage policy:

| `derivation_scope` | Digest contract | Confirmatory use |
|---|---|---|
| `FIXED_EXTERNAL` | Fixed independently before this campaign; cohort digest blank | Eligible |
| `DISCOVERY_DERIVED` | Cohort digest equals the exact frozen discovery cohort; submitter attests derivation was label-blind | Eligible after unchanged validation projection |
| `DISCOVERY_SUPERVISED` | Cohort digest equals the exact frozen discovery cohort; definition used target/background labels | Audit only in v0.1; excluded from association and modelling |
| `ALL_DATA_EXPLORATORY` | Cohort digest equals the full campaign cohort | Audit only; excluded from confirmatory inference |

`feature_definition_sha256` is always the 64-character lower-case SHA-256 of the immutable
motif, HMM, rule or other definition artifact. It and all derivation fields must agree across
rows sharing `(feature_type, feature_id)`. The cohort digest is the SHA-256 of canonical
UTF-8 rows `protein_id<TAB>sequence_sha256<NEWLINE>`, sorted by exact protein ID. This binds
a learned definition to the precise sequences allowed to teach it.

The cohort digest proves which sequences were declared; it cannot inspect an external
producer and therefore cannot prove that labels were unused. Selecting `DISCOVERY_DERIVED`
is an explicit submitter attestation that motif/profile/feature selection was label-blind.
Any definition selected using discovery target/background labels must instead declare
`DISCOVERY_SUPERVISED`. That honest audit-only scope is reserved for future validation-only
inference and cannot contribute discovery p-values, FDR, signatures or ML in v0.1.

All imported rows are published in `feature_assessments`; only assessed-positive,
confirmatory-eligible rows enter canonical `features`. `ASSESSED_WITH_FEATURE` and
`ASSESSED_NO_FEATURE` rows define each feature's denominator. Missing, not-assessed, failed
and excluded rows remain unknown and are never silently converted to zero. See
`examples/features.template.tsv`.

Recommended stable types include `MEME_MOTIF`, `STREME_MOTIF`, `HMM_PROFILE`,
`CONSERVED_RESIDUE`, `DISORDER_REGION` and `EXPERIMENTAL_SITE`. Record the producer and exact
version/reference; do not encode tool versions into an unstable display name.

## `domains.tsv`

```text
protein_id	domain_authority	assessment_status	domain_id	domain_name	start	end	score	e_value	evidence_source	evidence_reference
```

One file contains both hit and sentinel rows:

| `assessment_status` | Domain payload |
|---|---|
| `ASSESSED_WITH_HIT` | ID, name and one-based inclusive coordinates required |
| `ASSESSED_NO_HIT` | Domain payload must be blank |
| `NOT_ASSESSED` | Domain payload must be blank |
| `FAILED` | Domain payload must be blank |

Every protein is completed to at least one Pfam assessment row. Missing protein/Pfam pairs
become `NOT_ASSESSED`; they never become no-hit calls. Multiple hit rows may share the same
protein/authority status and provenance but must have unique coordinates.

Pfam hits generate both `PFAM_DOMAIN` presence and a coordinate-ordered
`PFAM_ARCHITECTURE`. Domain subsequences are extracted into the result for downstream use.

## `redundancy_clusters.tsv`

```text
protein_id	cluster_id	method	method_version	identity_threshold	coverage_threshold	evidence_reference
```

Each protein can occur in at most one supplied near-redundancy cluster. Thresholds are
fractions from zero to one and must be populated. The package also creates an exact
sequence-SHA-256 cluster for every protein.

## `structures.tsv`

```text
protein_id	structure_id	structure_source	structure_version	coordinate_path	coordinate_sha256	availability_status	mean_confidence	fold_id	fold_name	fold_authority	fold_authority_version	fold_evidence_reference	fold_evidence_status	analysis_eligibility_status	comparison_universe_ids
```

- `structure_id` is globally unique within the campaign.
- `coordinate_path` is optional but must point to a non-empty local PDB/mmCIF, optionally
  gzip-compressed, when populated.
- A local coordinate digest is calculated. A supplied digest must match it.
- `mean_confidence` is optional and bounded from 0 to 100.
- `availability_status` is controlled. `AVAILABLE` or `COMPLETE` requires a coordinate;
  every other state requires it to be absent.
- `analysis_eligibility_status` is required. `ELIGIBLE` requires an available coordinate.
  Confidence, sequence or validation exclusions use the matching `INELIGIBLE_*` state.
  `NOT_APPLICABLE_EXTERNAL_EVIDENCE` is reserved for reviewed fold/comparison evidence whose
  coordinates are not transferred into this campaign.
- `fold_evidence_status` is one of `ASSESSED_WITH_HIT`, `ASSESSED_NO_HIT`, `NOT_ASSESSED`,
  `FAILED` or `INPUT_UNAVAILABLE`. Every hit, assessed no-hit or failed assessment requires
  `fold_authority`, its exact `fold_authority_version`, and a reproducible
  `fold_evidence_reference`; a hit additionally requires `fold_id`. Unassessed or unavailable
  rows leave all fold provenance and payload fields blank. One authority release is allowed
  per campaign, preventing version-mixed denominators. Hit rows for the same authority/fold
  definition must use the same method or definition reference; protein/model identity remains
  available separately in `structure_id`.
- `comparison_universe_ids` is blank or a pipe-separated list of complete comparison
  universes containing this model. Declare every member, including models with no retained
  passing edge, so absence can be assessed safely.

Coordinate files are copied to content-addressed paths in the result.

## `structure_comparisons.tsv`

```text
protein_a_id	protein_b_id	comparison_tool	comparison_tool_version	tm_score	rmsd_angstrom	aligned_residue_count	coverage_a	coverage_b	comparison_status	source_record_id	comparison_universe_id	coverage_scope
```

Protein pairs must be distinct. `comparison_universe_id` names one complete search/alignment
campaign. `coverage_scope` is exactly `FULL_SEQUENCE`, `STRUCTURE_MODEL_RESIDUES` or
`DOMAIN_OR_CONSTRUCT` and states the denominator used for both coverage values. Successful
status `COMPLETE`, `SUCCESS` or `PASS` requires a TM score and bilateral coverage; unsuccessful
states must carry no numerical metrics. Components never bridge universes, scopes, tools or
tool versions.

## `alphafold_accessions.tsv`

```text
protein_id	uniprot_accession
```

Both columns are required and unique by protein. Acquisition is performed only when
`alphafold.enabled: true`. The adapter uses the official AlphaFold DB prediction endpoint,
accepts only HTTPS downloads from EBI, caches the returned PDB and records:

- API URL and model version;
- coordinate SHA-256;
- full-sequence agreement when the API supplies a sequence; and
- mean model pLDDT inferred from PDB B-factors when present.

A 404 is `MODEL_NOT_AVAILABLE`; other exhausted failures are `FAILED`. Downloaded models are
retained for audit, but only exact sequence agreement plus mean pLDDT at or above
`minimum_mean_plddt` produces `analysis_eligibility_status=ELIGIBLE`. Low-confidence,
confidence-unavailable and sequence-unverified coordinates never enter Foldseek or create a
structure-availability feature.
`MODEL_NOT_AVAILABLE`, `SEQUENCE_MISMATCH` and completed quality/verification exclusions enter
the AlphaFold assessment denominator as known non-eligible outcomes. `FAILED` and
`NOT_SELECTED` remain unknown rather than biological or structural absence.

## Completed OrthoFinder input

Choose exactly one route.

### Preferred: published `orthofinder-results` resource

Set `inputs.orthofinder.resource_dir`. Supported resource schema versions are 3 and 4. The
adapter requires a complete `run_manifest.json`, verifies every declared output, inspects
the physical DuckDB and requires a single consistent run identity. FASTA identifiers are
matched through the precursor sequence authority.

### Fallback: raw OrthoFinder output

Set `inputs.orthofinder.results_dir`. The directory must contain a `Log.txt` identifying
exactly OrthoFinder 2.5.5 or an OrthoFinder 3.x release, together with either
`Orthogroups.tsv` or hierarchical `N*.tsv` tables. The adapter supports `HOG` at an explicit
hierarchy node or `LEGACY_ORTHOGROUP`. `SequenceIDs.txt` is honoured where present.

Raw results are accepted only when `Log.txt` contains the exact
`OrthoFinder run completed` record written by OrthoFinder's official entry point after
result production. Version, completion, selected authorities and campaign memberships are
all checked during `protein-signatures validate`. No OrthoFinder command is constructed or
run by this package. Retain the scheduler record as additional operational provenance.

## Completed predecessor structural resource

`inputs.structural_alignment_resource` may point to the structural result itself or the
supported parent end-to-end layout. Exactly one resource root must resolve. Its completion
manifest and declared checksums are verified before three Parquet authorities are read:

- `structural_alignments.parquet`
- `pocket_comparisons.parquet`
- `structural_alignment_summary.parquet`

Only identifiers matching the campaign FASTA are imported. Within-group pocket evidence and
cross-class structural clustering remain separately named in the outputs.

## Explicit comparisons

Each comparison has:

```yaml
- comparison_id: kinase_vs_matched_non_kinase
  display_name: Kinases versus matched non-kinases
  target_label_ids: [protein:kinase]
  background_label_ids: [control:matched_non_kinase]
  description: Prespecified primary contrast.
```

Target and background label lists must be non-empty and disjoint. A protein inheriting into
both classes is a hard error. Matching is a responsibility of the supplied authority; the
package does not invent controls after seeing outcomes.
