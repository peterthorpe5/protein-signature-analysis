# Structural precursor boundary

## Scope and ownership

This document defines the compatibility boundary between `protein-signature-analysis` and
`peterthorpe5/E3_project_draft/e3_structural_alignment` v0.6.0. The precursor remains the
authority for its E3-specific model selection, within-group global US-align/TM-align runs,
pocket comparisons and group interpretations. This package is the authority for combining
that published evidence with a campaign FASTA, other structural routes, frozen
discovery/validation partitions, structural-feature derivation, association tests, models
and the portable result bundle.

The boundary is a completed, checksum-inventoried resource, not a Python import and not an
intermediate work directory. `protein-signature-analysis` does not invoke the precursor,
resume it, infer its missing outputs or recalculate its pocket analysis.

## Accepted resource layouts

Set `inputs.structural_alignment_resource` to either the completed structural result itself
or an end-to-end parent with exactly one supported child:

```text
RESOURCE/
├── provenance/
│   └── run_manifest.json
└── tables/
    ├── structural_alignments.parquet
    ├── pocket_comparisons.parquet
    └── structural_alignment_summary.parquet
```

```text
END_TO_END/09b_structural_alignment/structural_alignment/...
END_TO_END/structural_alignment/...
```

Resolution fails when no candidate or more than one candidate contains
`provenance/run_manifest.json`. The manifest must be a JSON object with:

- `status` equal to lowercase `complete`;
- `run_digest` equal to a 64-character lowercase hexadecimal SHA-256 value;
- `outputs`, a non-empty list of unique, safe, resource-relative paths, each with the exact
  non-negative `size_bytes` and a 64-character lowercase hexadecimal `sha256`; and
- `package_version`, which is retained as provenance (`unknown` when absent).

Every file declared in `outputs` is checked for path safety, existence, byte size and SHA-256
before a Parquet table is read. All three required tables must be declared in that inventory
and must exist as non-empty files; an unmanifested table is rejected even when present. The
adapter records `package_version` but does not use it as a version gate, so compatibility is
established by the manifest and column contracts below.

Configure the import as:

```yaml
inputs:
  structural_alignment_resource: /absolute/path/completed_structural_resource
```

`protein-signatures validate` reads and validates this resource. It does not run either the
precursor or Foldseek.

## Required Parquet contracts

Extra columns are allowed and ignored. Required columns are selected by name.

### `structural_alignments.parquet`

| Required column | Imported meaning and validation |
|---|---|
| `cluster_id` | Non-empty upstream group/cluster provenance. |
| `reference_accession` | Exact campaign FASTA `protein_id`; non-matching rows are not imported. |
| `mobile_accession` | Exact campaign FASTA `protein_id`; must differ from the reference. |
| `alignment_tool` | Non-empty upstream tool name, normally `US-align` or `TM-align`. |
| `status` | Case-normalised member of `COMPLETE`, `SUCCESS`, `PASS`, `FAILED`, `NOT_ASSESSED`, `INPUT_UNAVAILABLE`, `EXCLUDED`. |
| `tool_version` | Non-empty upstream tool version. |
| `aligned_length` | Optional positive integer. Required for successful rows because it supplies both coverage numerators. |
| `rmsd_angstrom` | Optional finite, non-negative number. |
| `minimum_tm_score` | Optional finite number in `[0, 1]`; required for successful rows. |

For `COMPLETE`, `SUCCESS` or `PASS`, TM score and both calculated coverages must be present.
For every other status, `aligned_length`, `rmsd_angstrom` and `minimum_tm_score` must all be
null or blank. The imported coverage values are
`min(1, aligned_length / campaign_FASTA_length)` for each accession and therefore have
`coverage_scope=FULL_SEQUENCE`. RMSD may remain null on a successful row.

The source identity is constructed as
`cluster_id|alignment_tool|reference_accession|mobile_accession`. An exact duplicate of that
four-part row identity is rejected. If the table contains rows but none match the campaign
FASTA, import fails.

### `pocket_comparisons.parquet`

| Required column | Imported meaning and validation |
|---|---|
| `cluster_id` | Required non-empty provenance for an assessed, campaign-matching row. |
| `reference_accession` | Exact campaign FASTA `protein_id`; non-matching rows are ignored. |
| `mobile_accession` | Exact campaign FASTA `protein_id`; non-matching rows are ignored. |
| `alignment_tool` | Required non-empty provenance for an assessed, campaign-matching row. |
| `status` | Only `ASSESSED`, `COMPLETE` and `PASS`, case-insensitively, can emit features. Other values are skipped. |
| `same_pocket_position_supported` | Explicit Boolean: Boolean value or `true`/`false`, `1`/`0`, `yes`/`no`. |
| `pocket_structure_conserved` | Explicit Boolean using the same accepted values. |

A true support flag creates a positive `STRUCTURAL_POCKET` feature for both proteins. The
two local feature names are `SAME_3D_POCKET_POSITION` and `CONSERVED_3D_POCKET`. They are
prefixed with an `ES3A_...` namespace derived from the upstream `run_digest` and
`cluster_id`, preventing collisions between precursor runs or groups. These are imported as
`ALL_DATA_EXPLORATORY` with `evidence_source=e3_structural_alignment`; they are not promoted
to discovery-derived confirmatory definitions.

The table supplies positive supported-pair evidence. This adapter does not construct a
complete pocket no-hit assessment universe from negative or skipped rows.

### `structural_alignment_summary.parquet`

The stable imported subset is:

| Kind | Required columns |
|---|---|
| Identity/provenance text | `cluster_id`, `primary_group_type`, `primary_group_id`, `reference_accession`, `alignment_tools`, `position_alignment_status`, `alignment_status`, `interpretation` |
| Non-negative integers | `alignment_tool_count`, `selected_accession_count`, `model_available_accession_count`, `aligned_accession_count`, `supported_accession_count`, `position_supported_accession_count` |
| Optional finite numbers | `group_support_fraction`, `group_position_support_fraction`, `mean_minimum_tm_score`, `mean_pocket_overlap_fraction`, `median_centroid_distance_angstrom` |

`(cluster_id, primary_group_type, primary_group_id)` must be non-empty and unique. The rows
are republished, in this stable subset, as `imported_structural_group_summaries` for audit
and presentation. They are not recomputed by this package.

For each imported alignment cluster, `selected_accession_count` is also a completeness
cross-check: the number of distinct campaign-matching accessions observed as alignment
endpoints must equal the declared count. A missing summary or unequal count is fatal because
the package could not justify absence calls for that comparison universe. This check proves
the membership set; it does not independently prove that every possible pair is present.

## What is reused and what is deliberately local

| Reused from the v0.6.0 resource | Deliberately owned or duplicated locally |
|---|---|
| Completed global US-align/TM-align rows, including tool/version, status, aligned length, RMSD and minimum TM score | A narrow read-only compatibility adapter, so this package has no runtime dependency on `E3_project_draft` |
| Assessed same-position and conserved-pocket support calls | Translation into canonical `structure_comparisons` and exploratory `STRUCTURAL_POCKET` features |
| Stable group summary values and interpretations | FASTA-based bilateral coverage calculation and exact campaign identifier filtering |
| Upstream `run_digest`, `package_version`, manifest and declared file hashes | Namespacing, assessment-universe checks, discovery-safe cluster construction, statistical analysis and result publication |

Canonical copies in a `protein-signature-analysis` result are intentional publication and
provenance duplication. They do not transfer scientific ownership or imply that the
predecessor analysis was rerun. Conversely, precursor-specific intermediate tables,
coordinates, logs, scripts and unsupported columns are not adopted merely because they are
present in the resource.

## Predecessor alignments versus Foldseek

The two routes can coexist, but they answer different questions and are never silently
substituted for one another.

| Property | Imported precursor alignment | Campaign Foldseek search |
|---|---|---|
| Scope | Pairs within an upstream `cluster_id` selected by the E3 precursor | Every explicitly eligible local model against the campaign model collection |
| Tool evidence | Published US-align/TM-align result | Locally executed `foldseek easy-search --alignment-type 1` |
| TM score | Published `minimum_tm_score` | Conservative minimum of query- and target-normalised TM scores |
| Coverage denominator | Authoritative full FASTA sequence length | Foldseek query/target model residue length |
| Coverage scope | `FULL_SEQUENCE` | `STRUCTURE_MODEL_RESIDUES` |
| RMSD | Preserved when published | Not emitted by the current Foldseek parser |
| Universe identity | Hash of precursor `run_digest` plus upstream `cluster_id` | Hash-keyed set of all eligible Foldseek input proteins |
| Reuse | The resource itself is immutable imported evidence | A checksum-valid cache keyed by coordinates, parameters and exact Foldseek version |

Foldseek requires at least two distinct proteins with explicitly eligible, available local
coordinate models. `maximum_hits` must be at least the number of eligible models so a query
cannot be silently truncated by the hit cap; the configured E-value still filters retained
hits. Imported and generated rows are merged only when `(comparison_tool,
source_record_id)` remains unique.

## Cluster and assessment-universe semantics

Three identifiers must not be conflated:

1. The precursor `cluster_id` names an upstream E3 group and defines the imported alignment
   assessment universe.
2. `comparison_universe_id` is a namespaced identity for the set in which qualifying
   presence and absence can be interpreted. Precursor universes are derived separately for
   each `(run_digest, cluster_id)`; Foldseek has one cache-keyed eligible-model universe.
3. The downstream `SC_...` `structure_clusters.cluster_id` names a thresholded connected
   component derived by this package. It is not the precursor cluster ID and is not a formal
   SCOP, CATH or ECOD fold assignment.

Passing edges require a successful status, inclusive minimum TM score, and inclusive
coverage threshold for both proteins. Edges are grouped by comparison universe, coverage
scope, tool and tool version before connected components are calculated. US-align,
TM-align and Foldseek therefore do not fuse into one component merely because they compare
the same proteins.

With campaign partitions, only discovery-to-discovery passing edges define a component.
A validation protein may be projected into a frozen component through a passing edge to a
discovery member; validation-to-validation edges cannot create or merge components. Because
the clustering method is connected-components single linkage, not every pair of members is
required to pass the thresholds.

Only a declared complete comparison universe supports negative `STRUCTURE_CLUSTER`
assessment. Observed positive endpoints alone are insufficient. The imported precursor
route supplies its checked cluster-local membership sets; Foldseek supplies its complete
eligible-query set. A supplied generic comparison table must separately declare universe
membership through `structures.tsv` if absence is to be interpreted.

## AlphaFold model boundary

Importing the precursor resource imports alignment and pocket evidence only. It does not
copy precursor coordinate files into `assets/structures`, create `structures` rows, or make
those models available for a new Foldseek run.

Models needed for campaign-native analysis must come through one of the package's model
routes:

- list existing models in `structures.tsv`, with exact source, version, coordinate path,
  checksum, confidence and eligibility provenance. Any supplied fold annotation also needs
  its classification-authority release and reproducible evidence reference; or
- provide an exact `protein_id` to UniProt accession mapping and enable AlphaFold DB
  acquisition.

AlphaFold DB acquisition downloads released models from the official EBI service. It does
not run AlphaFold prediction. The package records success and failure, model version,
coordinate checksum, sequence agreement and mean PDB B-factor/pLDDT; the configured mean
pLDDT threshold determines local-analysis eligibility. A precursor-generated or other local
prediction must be described by its actual predictor and version, not relabelled as an
AlphaFold DB acquisition.

## Safe archive and transfer examples

Archive the smallest completed precursor resource containing its manifest, required tables
and any assets needed for independent inspection. These commands do not modify the source:

```bash
SOURCE_DIR=/absolute/persistent/path/top_200_completed_result
ARCHIVE_DIR=/absolute/persistent/path
ARCHIVE_NAME=top_200_completed_result.tar.gz
ARCHIVE="${ARCHIVE_DIR}/${ARCHIVE_NAME}"

test -d "${SOURCE_DIR}"
find "${SOURCE_DIR}" -maxdepth 4 -type f -print | sort | sed -n '1,300p'
tar -C "$(dirname "${SOURCE_DIR}")" \
  -czf "${ARCHIVE}" \
  "$(basename "${SOURCE_DIR}")"
(cd "${ARCHIVE_DIR}" && sha256sum "${ARCHIVE_NAME}" > "${ARCHIVE_NAME}.sha256")
tar -tzf "${ARCHIVE}" | sed -n '1,300p'
```

On macOS, replace the checksum line with:

```bash
(cd "${ARCHIVE_DIR}" && shasum -a 256 "${ARCHIVE_NAME}" > "${ARCHIVE_NAME}.sha256")
```

Transfer the archive and checksum together into a new staging directory, verify before
extraction, and preserve the completed resource unchanged:

```bash
INCOMING=/absolute/persistent/path/incoming/e3_structural_alignment_v0.6.0
mkdir "${INCOMING}"
rsync -a /transfer/source/top_200_completed_result.tar.gz \
  /transfer/source/top_200_completed_result.tar.gz.sha256 \
  "${INCOMING}/"
cd "${INCOMING}"
sha256sum --check top_200_completed_result.tar.gz.sha256
tar -tzf top_200_completed_result.tar.gz | sed -n '1,300p'
tar -xzf top_200_completed_result.tar.gz
```

Then point `inputs.structural_alignment_resource` at the extracted result and run
`protein-signatures validate --config /absolute/path/campaign.yaml`. Before transfer, inspect
the archive for credentials, private scheduler configuration and unrelated data. Do not use
`rsync --delete` for a published authority.

## Current repository limitation

The top-200/full-1000 real-run inputs are not present in this repository. The repository
contains a small synthetic offline example and unit-test fixtures for the v0.6.0 adapter,
but not the real top-200 precursor result, its archive, or the full-1000 campaign FASTA,
labels, models, completed precursor resource or other production inputs. Consequently the
checked-in tests validate the contract and failure handling, not a scientific reproduction
of either real run. A real top-200 or full-1000 run must be transferred as an external,
checksum-preserved authority and reviewed before execution.

For operational hand-off and scheduler staging, continue with the
[structural cluster runbook](CLUSTER_RUNBOOK.md). See [input contracts](INPUT_CONTRACTS.md)
for all campaign inputs, [methods](METHODS.md) for the downstream statistical semantics, and
[scientific boundaries](SCIENTIFIC_BOUNDARIES.md) before interpreting structural evidence.
