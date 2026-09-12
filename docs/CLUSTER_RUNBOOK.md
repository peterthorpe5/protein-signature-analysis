# Structural cluster runbook

## Purpose

This runbook covers the hand-off from an existing cluster structural campaign into the
standalone signature workflow. It does not rerun OrthoFinder. It can either import the
checksum-complete predecessor E3 structural result or consume a generic structure inventory
and pairwise comparison table.

Do not start from an arbitrary intermediate directory merely because files are present. A
source is complete only when its declared completion manifest and outputs validate.

## Route A: predecessor structural resource

The current adapter accepts either the exact structural result root or a parent containing
one of these unambiguous layouts:

```text
SOURCE/
├── provenance/run_manifest.json
└── tables/
    ├── structural_alignments.parquet
    ├── pocket_comparisons.parquet
    └── structural_alignment_summary.parquet
```

or:

```text
SOURCE/09b_structural_alignment/structural_alignment/...
SOURCE/structural_alignment/...
```

The adapter verifies every manifest-declared output before reading the three Parquet
authorities. It imports:

- global US-align/TM-align records as pairwise structural comparisons;
- supported conserved/same-position pocket calls as `STRUCTURAL_POCKET` features; and
- group-level pocket/alignment summaries for audit and visualisation.

Only accessions matching the campaign FASTA are retained. Coverage is calculated against
those authoritative sequence lengths. The upstream package version and run digest remain in
run metadata.

Configure:

```yaml
inputs:
  structural_alignment_resource: /persistent/path/completed_structural_run
```

This route complements an all-versus-all Foldseek screen. It is not duplicated or silently
recalculated.

## Route B: generic models and pairwise comparisons

Supply:

- `structures.tsv`, with one structure/model identity, provenance, confidence, controlled
  analysis eligibility, optional fold call, optional coordinate path and pipe-separated
  complete comparison-universe memberships per row; and
- `structure_comparisons.tsv`, with complete TM score, RMSD, aligned length and bilateral
  coverage when already calculated, plus a universe ID and explicit coverage-denominator
  scope.

Enable Foldseek if usable coordinate paths are present and a new all-model search is
required. Every usable model is queried against the campaign collection. Pre-flight and
runtime require `maximum_hits` (1,000 by default) to be at least the candidate-model count,
so the result cap cannot silently truncate a query; the configured E-value still filters the
published edge graph. Foldseek execution is cached by coordinate content, parameters and
exact tool version. It is launched as an argument vector, never through an interpolated shell
command.

Do not import the same comparison record through two routes: `(comparison_tool,
source_record_id)` must be unique.

## AlphaFold models

For missing models, create an exact mapping:

```text
protein_id	uniprot_accession
```

and enable AlphaFold DB retrieval. This downloads released models from the official EBI
service; it does not run AlphaFold inference. A missing model, failed transfer, sequence
mismatch and low-confidence model remain separate outcomes. Downloaded coordinates are kept
for provenance, but only exact sequence-verified models meeting the configured mean-pLDDT
threshold enter Foldseek.

For models produced by another local predictor or cluster pipeline, list them in
`structures.tsv` with the exact predictor/version and coordinate checksum rather than
mislabeling them as AlphaFold DB models.

## Package a top-200 or completed 1,000-cluster run for adapter inspection

Set `SOURCE_DIR` to the smallest completed result root that contains its manifest, tables
and required assets. The command is read-only with respect to the source:

```bash
SOURCE_DIR=/absolute/persistent/path/top_200_completed_result
ARCHIVE=/absolute/persistent/path/top_200_completed_result.tar.gz

test -d "${SOURCE_DIR}"
find "${SOURCE_DIR}" -maxdepth 4 -type f -print | sort | sed -n '1,300p'
tar -C "$(dirname "${SOURCE_DIR}")" \
  -czf "${ARCHIVE}" \
  "$(basename "${SOURCE_DIR}")"
sha256sum "${ARCHIVE}" > "${ARCHIVE}.sha256"
tar -tzf "${ARCHIVE}" | sed -n '1,300p'
```

On macOS use:

```bash
shasum -a 256 "${ARCHIVE}" > "${ARCHIVE}.sha256"
```

Before upload, confirm the archive contains no credentials, private scheduler configuration
or unrelated data. Upload both the archive and checksum. If the top-200 layout differs from
the contract above, preserve it exactly; the adapter should be extended from evidence rather
than filenames being guessed.

## Scheduler pattern

Heavy temporary work belongs under scheduler-provided `TMPDIR`; persistent outputs belong
under an explicit campaign directory. Never hard-code `/tmp` and never use `rsync --delete`.
This pattern assumes `campaign.yaml` was first created with `start_from_inputs.sh
--initialise-only`, scientifically reviewed, validated and frozen with its referenced inputs.

```bash
#!/usr/bin/env bash
set -Eeuo pipefail

: "${TMPDIR:?The scheduler must provide TMPDIR}"
PERSISTENT_RUN=/persistent/path/protein_signatures/e3_1000_2026_09
JOB_TOKEN="${SLURM_JOB_ID:-manual}"
JOB_WORK_DIR="${TMPDIR}/protein_signatures_${JOB_TOKEN}"
FINAL_RESULT="${PERSISTENT_RUN}/result"
TRANSFER_RESULT="${PERSISTENT_RUN}/.result.transfer_${JOB_TOKEN}"

mkdir -p "${JOB_WORK_DIR}/inputs" "${PERSISTENT_RUN}"
test ! -e "${FINAL_RESULT}"
test ! -e "${TRANSFER_RESULT}"
rsync -a /persistent/path/frozen_inputs/ "${JOB_WORK_DIR}/inputs/"

./run_protein_signature_analysis.sh \
  --config "${JOB_WORK_DIR}/inputs/campaign.yaml" \
  --output-dir "${JOB_WORK_DIR}/result" \
  --threads "${SLURM_CPUS_PER_TASK:-1}"

conda run --name protein_signature_analysis \
  protein-signatures verify --resource "${JOB_WORK_DIR}/result"

mkdir "${TRANSFER_RESULT}"
rsync -a "${JOB_WORK_DIR}/result/" "${TRANSFER_RESULT}/"
conda run --name protein_signature_analysis \
  protein-signatures verify --resource "${TRANSFER_RESULT}"
mv "${TRANSFER_RESULT}" "${FINAL_RESULT}"
conda run --name protein_signature_analysis \
  protein-signatures verify --resource "${FINAL_RESULT}"
```

The nonexistent-destination checks prevent stale undeclared files from surviving an earlier
copy. The transfer is verified in a sibling directory and renamed on the persistent
filesystem only after it is complete. The final `verify` command must succeed before the
scheduler job is treated as complete.

## Pre-flight checklist

1. Freeze the FASTA and label-assignment authority.
2. Confirm all structure `protein_id` values match the first exact FASTA header token.
3. Confirm coordinate files are non-empty and checksums match when supplied.
4. Confirm upstream structural and OrthoFinder manifests are complete.
5. Select structural TM-score and bilateral-coverage thresholds before viewing enrichment.
6. Select target/background comparisons before viewing structural clusters.
7. Run `protein-signatures validate` on the cluster login node.
8. Submit to a new output directory; never overwrite a partial or published result.

## What to inspect after completion

Start with files; the app is optional:

```text
result/analysis/
├── 00_run_information/
├── 01_proteins_and_curation/
├── 02_homology_and_partitions/
├── 03_sequence_and_domains/
├── 04_structures_and_folds/
├── 05_association_statistics/
├── 06_explainable_models/
└── 99_final_results/
```

Every canonical table has a TSV and formatted Excel copy in its owning stage. Key decision
tables and figures are mirrored into `99_final_results`. Figures are PNG for quick review,
SVG for editing and PDF for publication. `report_inventory.tsv` provides a complete index.

Before biological interpretation:

- check `COMPLETED.json` and run `protein-signatures verify`;
- inspect structure and Pfam availability plots for class imbalance;
- inspect `excluded_mixed_unit_count` and discovery/validation block sizes;
- distinguish structure clusters from reviewed fold assignments;
- prefer study-wide q-values when surveying many subclasses; and
- read SHAP beside held-out metrics and association results.
