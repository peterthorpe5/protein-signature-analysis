# Bounded 3/4/5-mer E3 trial

The optional `analysis.kmer_vocabulary_policy: prevalence_ranked` limits the
number of tested k-mer definitions while keeping a separate limit on all raw
discovery candidates. It distributes the tested-feature budget across the
configured lengths. The mode uses discovery sequences without class labels and
projects its frozen vocabulary to validation proteins. `strict` remains the
default for other campaigns.

The 5-mer trial is exploratory. An omitted candidate cannot be called absent or
unenriched; longer exact strings do not by themselves locate a structural motif.
Inspect the selected counts by length in the analysis log and use held-out
evidence before interpreting a candidate.

## Campaign setup

Use the supplied trial `campaign.yaml` at
`$NEW_WORK/campaign.yaml`. Its `inputs` paths deliberately name new
`prepared_inputs/` and `evidence_label_bundle/` directories. The completed-E3
workflow generates those authorities and their path-bound checksum markers.
Copying the old evidence bundle into the new directory would invalidate its
marker. Validation of the new YAML succeeds only *after* those new authorities
have been created.

The source E3 run, its coordinate models, the structural resource and the
OrthoFinder results stay under `$RUN_ROOT`. The old signature campaign is no
longer needed after the new result is independently verified. Its Foldseek cache
can be copied to the new cache location without changing the content-keyed
search identity.

```bash
cd /gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/protein-signature-analysis
git pull --ff-only origin main

RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"
BASE="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs"
OLD_WORK="$BASE/e3_all1972_caliper_complete_20260928"
NEW_WORK="$BASE/e3_all1972_k5_trial_20261002"

if [[ -d "$OLD_WORK/.protein_signature_cache/foldseek" ]]; then
  mkdir -p "$NEW_WORK/.protein_signature_cache/foldseek"
  rsync -a "$OLD_WORK/.protein_signature_cache/foldseek/" \
    "$NEW_WORK/.protein_signature_cache/foldseek/"
fi

# The YAML must be in place at $NEW_WORK/campaign.yaml first.
./run_completed_e3_workflow.sh \
  --phase all \
  --run-root "$RUN_ROOT" \
  --work-dir "$NEW_WORK" \
  --campaign-id e3_all1972_k5_trial_20261002 \
  --evidence-led-labels \
  --evidence-rules e3 \
  --minimum-mean-plddt 50 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 256G \
  --slurm-time 2-00:00:00 \
  --threads 24
```

The first submission stops after it generates the new evidence bundle. When
that job finishes, check the marker and control summary:

```bash
conda run --no-capture-output --name protein_signature_analysis \
  protein-signatures verify-evidence-labels \
  --bundle-dir "$NEW_WORK/evidence_label_bundle"

cat "$NEW_WORK/evidence_label_bundle/control_match_summary.tsv"
```

If the evidence proposal is suitable for this provisional trial, submit the
same command with `--accept-provisional-evidence-labels` added. That submission
adopts the supplied YAML, validates the newly generated input authorities and
runs the analysis. Do not use the generic launcher with the old path-bound
evidence bundle for this setup.

When the job finishes, verify the complete new result:

```bash
conda run --no-capture-output --name protein_signature_analysis \
  protein-signatures verify --resource "$NEW_WORK/result"
```

After `VALID` and after checking that no other work needs the old campaign,
the old HPC work directory may be removed. Keep `$RUN_ROOT`: it contains
external input authorities still used by the new campaign.
