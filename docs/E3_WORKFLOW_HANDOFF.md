# Completed E3 workflow hand-off

## Purpose

This guide starts an independent protein-signature campaign from a completed
`E3_project_draft/e3_end_to_end_workflow` result. The predecessor remains immutable and
authoritative for its own selection, domain retrieval, AlphaFold assets, within-group
US-align/TM-align analysis and pocket analysis. This package performs the separate
cross-class signature analysis.

The bridge does not run OrthoFinder, AlphaFold inference or the predecessor workflow. It
also never converts an upstream family annotation into a reviewed class label.

## Source map

Set `--run-root` to the end-to-end run directory itself. The adapter uses these authorities:

| Predecessor source | Signature use |
|---|---|
| `11_app_ready/stage_manifest.json` | Formal completed-run gate |
| `05_orthology/orthology/tables/candidate_group_member_sequences.parquet` | Exact sequence, accession and HOG context |
| `06_domains/tables/domain_hits.parquet` | Coordinate-resolved Pfam hits and E3 review hints |
| `06_domains/tables/domain_summary.parquet` | Explicit assessed-hit, assessed-no-hit and unavailable Pfam state |
| `09_ligandability/tables/reused_asset_manifest.parquet` | Checksum-bound AlphaFold coordinate inventory for new Foldseek analysis |
| `09_ligandability/tables/reused_model_quality.parquet` | Mean pLDDT used to gate coordinates before Foldseek |
| `09b_structural_alignment/structural_alignment/` | Completed within-group US-align/TM-align and pocket evidence |
| `04_orthofinder/Results/` | Published OrthoFinder membership and identifier tables |
| `04_orthofinder/stage_manifest.json` | Completed Stage 04 gate and checksum inventory |
| `04_orthofinder/orthofinder_authority.tsv` | Reviewed archive, version and decision authority |
| `04_orthofinder/orthofinder_reuse_validation.tsv` | Exact five-file extraction-validation ledger |

Every consumed Stage 05, 06 and 09 Parquet file must be listed exactly once in its complete
stage manifest, with matching byte size and SHA-256. Stage 09b supports both published
contracts: a standalone component manifest with an `outputs` inventory, and the end-to-end
aggregate manifest with a `datasets` inventory. For the aggregate, every declared scientific
Parquet checksum is cross-checked against
`09b_structural_alignment/stage_manifest.json`, including the aggregate run manifest itself.
Input-file and coordinate-file checksums are verified again by the signature campaign.

The reused Stage 04 intentionally has no OrthoFinder `Log.txt`. The three adjacent workflow
authorities replace that raw-run completion evidence: they must agree on the exact required
files, sizes and SHA-256 digests. The signature adapter reads the existing files in place;
do not copy them into the new work directory and do not add or link a log file. Direct raw
OrthoFinder 2.5.5/3 inputs outside this wrapper continue to require OrthoFinder's official
completed-run marker.

`10_integrated_resource/final_results` is deliberately absent from this map. Those files
rank predecessor candidates; using their ranks as signature labels would leak the outcome
of the earlier prioritisation into this analysis.

## Why the workflow has phases

A completed structural run proves that models and pairwise evidence exist. It does not prove
that every HOG member is a particular E3 subclass, and it does not create matched non-E3
controls. Automatically propagating a seed category to every orthologue would make the
signature analysis circular.

The completed-E3 Snakemake DAG therefore enforces this order:

1. `ensure_e3_preparation` creates or verifies the immutable review bundle;
2. `stage_e3_label_review` safely creates or adopts the separate editable label authority;
3. a curator assigns exact E3 class/component labels and matched-control labels;
4. `approve` validates those assignments and binds curator approval to their exact SHA-256;
5. `verify_e3_label_review` fails if either the approval chain or reviewed file changed;
6. `initialise_e3_campaign` creates or adopts and validates `campaign.yaml`;
7. `run_e3_campaign` revalidates the campaign and executes the atomic analysis; and
8. `verify_e3_campaign` independently checks every result and input checksum.

Only step 3 is manual. Once the checksum-bound approval exists, `--phase all` executes steps
5–8 unattended. Calling `--phase all` before approval stops cleanly at step 2 rather than
guessing labels.

## Phase 1: prepare inputs

```bash
RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"
SIGNATURE_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_v0_1_0_20260914"

./run_completed_e3_workflow.sh \
  --phase prepare \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_WORK}" \
  --campaign-id "e3_all1972_signatures_20260914" \
  --minimum-mean-plddt 50 \
  --log-level INFO \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 32G \
  --slurm-time 04:00:00 \
  --threads 18 \
  --slurm-dry-run
```

Repeat without `--slurm-dry-run` after inspecting the exact `sbatch` command. This phase runs
the E3 DAG through `review_ready`. The real all-1972 run recorded about 4 GiB maximum RSS;
32 GiB leaves conservative headroom. An existing prepared bundle is checksum-verified and
reused rather than rebuilt. Follow the returned job identifier and retained logs:

```bash
squeue --job JOB_ID
tail -f "${SIGNATURE_WORK}/slurm_logs/protein_signature_prepare_JOB_ID.out"
```

The new `prepared_inputs/` directory contains:

| File | Meaning |
|---|---|
| `PREPARED.json` | Completion marker, source identities, counts and checksums |
| `proteins.faa` | Unique exact accession-bearing Stage 05 sequences |
| `domains.tsv` | Pfam hits plus explicit no-hit/not-assessed sentinel rows |
| `structures.tsv` | Sequence-matched, checksum-verified Stage 09 coordinates, pLDDT and eligibility |
| `label_assignments.REVIEW_REQUIRED.tsv` | All-`UNMAPPED` safe starter; never use directly for analysis |
| `e3_label_curation_review.tsv` | Exact source-cluster JSON arrays, HOG, species, Pfam, upstream-family/role hints and empty decision columns |
| `alphafold_accessions.MISSING_MODELS_REVIEW_REQUIRED.tsv` | Canonical accessions lacking a reused model, for optional bounded follow-up |
| `source_inventory.tsv` | Exact files and SHA-256 values used during preparation |

The bridge merges repeated accessions only when their sequences agree exactly. Conflicting
sequence, Pfam, model-quality, path, byte-count or coordinate-checksum records fail closed.
Stage 05 `cluster_id` values are opaque provenance rather than signature-package identifiers;
therefore composite DeepClust values such as
`Arabidopsis_thaliana@@sp|B3H578|PHD1_ARATH` are retained exactly. Because `|` is valid inside
these values, the `cluster_ids` review column is a deterministic JSON array rather than an
ambiguous delimiter-joined string. HOG and orthogroup identifiers remain strictly validated.
Coordinates below the declared mean-pLDDT threshold, or without a model-quality value, remain
in the audit table but receive an explicit ineligible state and never enter Foldseek. Preparation
requires at least two eligible models and records both the available and eligible counts in
`PREPARED.json`.

Snakemake state is kept separately under `workflow_state/e3/`. At this checkpoint the key
markers are `01_preparation/PREPARED_VERIFIED.json` and
`02_label_review/REVIEW_READY.json`. The editable
`reviewed_label_assignments.tsv` is deliberately not a Snakemake output: failed-job cleanup
can never remove curator work.

## Phase 2: curate and initialise

The `stage_e3_label_review` rule makes the separate copy on a new campaign. If, as in the
current all-1972 run, it was already copied manually, the rule adopts it and records that
decision without changing a byte. Add or duplicate rows when a protein has more than one
reviewed hierarchical label. Only `REVIEWED_POSITIVE` assignments enter analysis membership.
A control protein must be reviewed positively into the relevant `control:...` label; marking
it `REVIEWED_NEGATIVE` against a target does not create background membership.

```bash
# Curate ${SIGNATURE_WORK}/reviewed_label_assignments.tsv.

./run_completed_e3_workflow.sh \
  --phase approve \
  --work-dir "${SIGNATURE_WORK}" \
  --curator "Peter Thorpe" \
  --review-note "Reviewed E3 subclasses, component roles and controls"
```

Approval refuses the byte-identical all-`UNMAPPED` template, incomplete FASTA coverage,
unknown labels, incompatible component roles, unreviewed evidence states, assignments without
at least one reviewed-positive analysis target, and any populated target lacking its
profile-resolved matched background. The marker
`workflow_state/e3/02_label_review/REVIEW_APPROVED.json` is immutable. To change approved
labels, use a new versioned work directory and issue a new approval; do not edit an existing
approval marker.

To stop after campaign creation for a deliberate configuration review:

```bash
./run_completed_e3_workflow.sh \
  --phase initialise \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_WORK}" \
  --campaign-id "e3_all1972_signatures_20260914" \
  --threads 18
```

The generated `campaign.yaml` is now the sole configuration authority. Initialisation does
not run the analysis. Inspect at least:

- all enabled target/background comparisons and their reviewed counts;
- `analysis.structural_tm_score_threshold` and `structural_minimum_coverage`;
- HOG-grouped discovery/validation partition settings;
- k-mer and minimum-sample settings;
- `foldseek.maximum_hits`, automatically set to the number of prepared eligible models so
  an all-model query cannot be silently truncated; and
- Foldseek cache and final output paths on persistent storage.

The Stage 09b import and new Foldseek route are complementary. Stage 09b contributes precise
within-predecessor-group global alignment and pocket evidence. Foldseek searches every
eligible prepared model against the entire prepared model collection, which is the route
that can discover folds shared across different HOGs or E3 classes.

## Phase 3: run start-to-finish through Snakemake and Slurm

After approval, a single submitted allocation can run the completed-E3 Snakemake DAG from
review verification through final result verification. It creates `campaign.yaml` if absent,
or adopts it only if all critical input authorities are unchanged:

```bash
./run_completed_e3_workflow.sh \
  --phase all \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_WORK}" \
  --campaign-id "e3_all1972_signatures_20260914" \
  --threads 24 \
  --log-level INFO \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 128G \
  --slurm-time 2-00:00:00
```

After `--phase initialise`, the resulting YAML is also a normal generic campaign and can be
run with the generic logout-safe Snakemake controller and Slurm executor:

```bash
./submit_protein_signature_workflow_slurm.sh \
  --config "${SIGNATURE_WORK}/campaign.yaml" \
  --work-dir "${SIGNATURE_WORK}" \
  --account barton \
  --partition barton \
  --threads 24 \
  --memory-mb 128000 \
  --runtime-minutes 2880 \
  --max-jobs 10 \
  --dry-run
```

Remove `--dry-run` after review. The first command uses one Slurm allocation for the complete
E3 adapter DAG. The second submits a small durable controller and lets the generic DAG request
each rule allocation. Do not launch both for the same campaign.

Run into a new result directory. A pre-existing partial directory is not overwritten. The
`--resume` option only verifies and reuses a complete immutable result with the same run
identity and unchanged input authorities; it does not continue a partial result.

Foldseek can be the dominant work unit. Before submission, record the prepared structure
count from `PREPARED.json`, confirm storage capacity for its cache and choose CPU, memory and
wall time from the site's top-200/top-1,000 benchmark experience.
The full execution and recovery contract is in [Snakemake and Slurm](SNAKEMAKE_SLURM.md).

## Phase 4: verify and inspect

```bash
./run_completed_e3_workflow.sh \
  --phase verify \
  --work-dir "${SIGNATURE_WORK}"

./run_protein_signature_app.sh \
  --resource "${SIGNATURE_WORK}/result"
```

File-first results are under numbered directories in `result/analysis/`, with key tables and
figures mirrored into `99_final_results`. Canonical tables are supplied as TSV and formatted
XLSX; figures are supplied as PNG, SVG and PDF. The read-only app uses the same verified
Parquet/DuckDB authorities.

## Interpretation boundary

Association tests use explicit assessed universes, Fisher exact tests and both local and
study-wide Benjamini-Hochberg FDR. SHAP is always generated for fitted group-aware models.
Neither a significant Pfam/fold/pocket association nor a high model score proves E3
function, substrate recognition or ligandability. Check held-out HOG performance, class and
control counts, structure/Pfam missingness and correlated parent/child hypotheses before
biological interpretation.
