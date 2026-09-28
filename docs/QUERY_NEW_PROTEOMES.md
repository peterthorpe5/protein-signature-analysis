# Query a new proteome or protein group

The reference campaign defines a protein class through a profile, explicit
target/background comparisons and discovery/validation evidence. A query applies
the **frozen, positively enriched, held-out validated** signatures to new
proteins. The same commands work with an E3 subclass, a custom profile or another
protein family. Query output is candidate evidence for review, not a class
probability or a verified functional annotation.

## Reference and query authorities

First finish and verify a reference campaign. Supply its canonical
`tables/signatures.tsv` and `tables/comparisons.tsv` (or their `.tsv.gz` forms).
These files fix the feature IDs and exact target/background label definitions.
By default, only `DECISION_CANDIDATE__VALIDATED_STUDY_WIDE` features with a
positive held-out prevalence difference are selected. Add
`--include-within-comparison` to admit the weaker, locally validated tier.
Features that were not validated do not become query markers.

The query FASTA is the complete protein identifier authority. Exact amino-acid
k-mers are scanned directly against every sequence with their reference IDs;
each is explicitly present or assessed absent. Supply one or more other feature
assessment TSVs with these required columns:

| Column | Meaning |
|---|---|
| `protein_id` | Exact FASTA identifier |
| `feature_type` | Reference signature feature type |
| `feature_id` | Frozen reference feature ID |
| `evidence_status` | `ASSESSED_WITH_FEATURE`, `ASSESSED_NO_FEATURE`, `NOT_ASSESSED`, `FAILED` or `EXCLUDED` |

Additional columns are allowed. A missing row means **unknown**, never assessed
absence. A no-feature row is appropriate only after a successful, complete
assessment of that protein and definition. Duplicate identical assessments are
accepted; conflicting calls stop publication. Supplied k-mer calls must agree
with the FASTA scan. Rows for valid feature IDs outside the selected signatures
are ignored and counted in the log.

```bash
protein-signatures query-signatures \
  --signatures /reference/result/tables/signatures.tsv \
  --comparisons /reference/result/tables/comparisons.tsv \
  --query-fasta /new_organism/proteins.faa \
  --query-features /new_organism/domain_assessments.tsv \
  --query-features /new_organism/other_feature_assessments.tsv \
  --comparison-id selected_class_vs_controls \
  --output-dir /new_organism/signature_query
```

`--query-features` is optional for an exact-k-mer-only query. Other feature types
stay unknown without their own assessment rows. The command creates a new
directory atomically; an existing directory is never replaced. Repeat
`--comparison-id` to query several exact IDs from `comparisons.tsv`; omit it
to query every class in the reference campaign.

## Query orthogroups

Supply a two-column `protein_id`, `unit_id` TSV and add
`--query-units /new_organism/orthogroup_members.tsv`. Every FASTA protein must
belong to at least one unit; overlapping groups are allowed. For one feature,
a group is present if any member is positive, assessed absent only when **all**
members were assessed negative, and otherwise unknown. A group hit lists the
positive member IDs. Group evidence cannot by itself assign every member the
same molecular function.

## Project frozen structural neighbourhoods

If new predicted structures and query-to-reference structural alignments are
available, project only onto **discovery members** of the reference campaign's
`tables/structure_clusters.tsv`:

```bash
protein-signatures project-structural-signatures \
  --reference-clusters /reference/result/tables/structure_clusters.tsv \
  --query-comparisons /new_organism/query_to_reference_alignments.tsv \
  --query-fasta /new_organism/proteins.faa \
  --output-dir /new_organism/structural_projection

protein-signatures query-signatures \
  --signatures /reference/result/tables/signatures.tsv \
  --comparisons /reference/result/tables/comparisons.tsv \
  --query-fasta /new_organism/proteins.faa \
  --query-features /new_organism/structural_projection/projected_structural_features.tsv \
  --output-dir /new_organism/signature_query_with_structure
```

The alignment input is a TSV with `protein_a_id`, `protein_b_id`,
`comparison_tool`, `comparison_tool_version`, `tm_score`, `coverage_a`,
`coverage_b`, `comparison_status`, `comparison_universe_id` and
`coverage_scope`. Each row must have exactly one query endpoint; query IDs must
not overlap the reference IDs. Passing pairs must meet the frozen TM-score and
**both** coverage thresholds with the same tool, version and coverage scope.
The command takes supplied alignments; it does not run structure prediction or
Foldseek. It emits positive memberships only. An omitted hit remains unknown
because a retained-hit file cannot establish an exhaustive negative result.
Structural clusters are connected neighbourhoods; the command does not locate
a local residue motif or prove every member is pairwise similar.

## Read the result

`candidate_summary.tsv` contains one row per query protein or group and
reference comparison, including counts of present, assessed absent and unknown
features. `candidate_feature_hits.tsv` contains the positive reference markers,
their held-out q-value and, for groups, their positive members.
`QUERY_COMPLETED.json` records input and output SHA-256 hashes, selected tier and
row counts. Structural projection has its own checksum-bound marker.

Presence of one or more associated features prioritises a candidate for review.
It is not a calibrated posterior, a claim of enzyme activity or an independent
validation on the new organism. Differences in proteome assembly, taxon, domain
scanner, structure quality and alignment completeness can change assessment
coverage. Confirm the candidates with reviewed orthology, domain architecture,
local structural geometry and independent biological evidence.
