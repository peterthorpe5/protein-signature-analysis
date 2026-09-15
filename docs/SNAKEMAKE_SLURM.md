# Snakemake and Slurm execution

## Scope

The packaged workflow is protein-agnostic. Its required authority is a reviewed
`campaign.yaml`, normally created by `start_from_inputs.sh`. That YAML can refer to any
supported protein FASTA, label profile, Pfam/domain ledger, external feature table, structure
inventory, AlphaFold request table, structural comparison resource and OrthoFinder 2.5.5 or
3.x result. The workflow consumes OrthoFinder output and never runs OrthoFinder.

`workflow/Snakefile` is the generic three-boundary DAG. The optional
`workflow/E3Snakefile` extends it for one completed-predecessor layout with preparation,
human-review staging, checksum-bound approval verification and campaign initialisation. It
still produces the same generic `campaign.yaml` and never runs OrthoFinder.

## Why the analysis rule is transactional

The Python engine calculates interdependent sequence, domain, structural, association,
held-out validation and SHAP tables in one controlled transaction. It publishes `result/`
only after every table, figure, Parquet file, DuckDB table and checksum is complete. Splitting
those internals into independently publishable Snakemake rules would weaken that atomicity and
duplicate scientific logic.

Snakemake therefore owns three clear boundaries:

1. `validate_campaign` checks the configuration and current input authorities and writes a
   checksum-bound validation marker;
2. `run_campaign` performs the complete atomic analysis; and
3. `verify_campaign` independently verifies the immutable result and every original input
   checksum before writing the final workflow marker.

The launcher forces `verify_campaign` on every non-dry invocation. A pre-existing completed
result is reusable only with `--resume`, matching configuration/run identity and unchanged
input authorities. A partial result is never resumed or overwritten.

Only files below `workflow_state/` are declared as Snakemake outputs. Files inside the
immutable `result/` are deliberately not declared, so Snakemake's failed-job cleanup can
never remove a scientific manifest or completion marker. The analysis rule writes its own
external marker only after the Python publisher and an immediate checksum verification both
succeed.

## Persistent layout

For `--work-dir /data/campaign_001`, the recommended layout is:

| Path | Purpose |
|---|---|
| `campaign.yaml` | Frozen scientific configuration |
| `result/` | Immutable, manifested scientific result |
| `workflow_state/01_validation/VALIDATED.json` | Validation boundary marker |
| `workflow_state/02_analysis/RESULT_VERIFIED.json` | External marker for the atomic result transaction |
| `workflow_state/03_verification/VERIFIED.json` | Final result-and-input verification marker |
| `workflow_state/logs/` | Rule-level application logs |
| `workflow_state/slurm/` | Controller stdout and stderr |
| `workflow_state/.snakemake/` | Snakemake metadata and Slurm executor logs |

All paths supplied to Slurm should be absolute and reside on persistent shared storage.

## Local or single-allocation execution

Use the local profile on a workstation or from inside a Slurm allocation that already has
enough memory for the whole analysis:

```bash
./run_protein_signature_analysis.sh \
  --config /data/campaign_001/campaign.yaml \
  --output-dir /data/campaign_001/result \
  --workflow-state-dir /data/campaign_001/workflow_state \
  --profile local \
  --threads 24 \
  --memory-mb 64000 \
  --runtime-minutes 1440
```

`memory-mb` and `runtime-minutes` are recorded as Snakemake resources. The enclosing local
or Slurm allocation remains the actual resource limit.

## Durable Slurm controller

The Dundee E3 production convention uses account `barton`, partition `barton`, Snakemake 9
and the Slurm executor plugin. The generic submitter adopts those defaults but exposes them
as flags:

```bash
./submit_protein_signature_workflow_slurm.sh \
  --config /data/campaign_001/campaign.yaml \
  --work-dir /data/campaign_001 \
  --account barton \
  --partition barton \
  --threads 24 \
  --memory-mb 64000 \
  --runtime-minutes 1440 \
  --max-jobs 10 \
  --controller-memory 4G \
  --controller-time 3-00:00:00 \
  --dry-run
```

Review the shell-escaped command, then remove `--dry-run`. The login-node process returns
after `sbatch` accepts the controller. The controller holds an exclusive `flock` lock for the
workflow state and uses the packaged `profiles/slurm/config.v8+.yaml` profile to submit rule
jobs. Do not submit a second controller for the same state directory.

Monitor the returned controller identifier with:

```bash
squeue --job JOB_ID
tail -f /data/campaign_001/workflow_state/slurm/controller_JOB_ID.out
```

If the controller ends unexpectedly, inspect its stderr, the rule logs and Snakemake's Slurm
logs before resubmission. Correct the external cause, then repeat the submit command with
`--resume`. A valid completed result is verified and reused; an incomplete directory still
fails closed and must be investigated rather than deleted automatically.

## Completed-E3 start-to-finish DAG

The predecessor adapter has one deliberate human checkpoint. First submit through the
`review_ready` target:

```bash
./run_completed_e3_workflow.sh \
  --phase prepare \
  --run-root /absolute/path/completed_e3_run \
  --work-dir /absolute/path/signature_campaign \
  --campaign-id e3_signatures_20260914 \
  --minimum-mean-plddt 50 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 32G \
  --slurm-time 04:00:00 \
  --threads 18 \
  --slurm-dry-run
```

Remove only `--slurm-dry-run` to submit. The worker command omits `--submit-slurm`, preventing
recursive submission, and confirms that the allocation contains at least the requested
thread count. Standard output and error are retained under `WORK_DIR/slurm_logs/` using the
Slurm job identifier. The DAG creates or verifies `prepared_inputs/`, then safely creates or
adopts `reviewed_label_assignments.tsv`; it never overwrites that human-owned file.

Preparation streams projected Parquet columns in bounded batches and streams FASTA output;
it no longer materialises every full upstream Parquet row twice. It still retains one exact
deduplicated sequence and its merged context per protein. The observed all-1972 job used
about 4 GiB peak RSS, so the documented 32 GiB request includes substantial headroom.

After curating that TSV, issue a checksum-bound approval from the login node:

```bash
./run_completed_e3_workflow.sh \
  --phase approve \
  --work-dir /absolute/path/signature_campaign \
  --curator "Curator name"
```

Then submit the same campaign through the full DAG:

```bash
./run_completed_e3_workflow.sh \
  --phase all \
  --run-root /absolute/path/completed_e3_run \
  --work-dir /absolute/path/signature_campaign \
  --campaign-id e3_signatures_20260914 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 128G \
  --slurm-time 2-00:00:00 \
  --threads 24
```

The full target verifies the approval checksum, creates or adopts and validates
`campaign.yaml`, revalidates it immediately before compute, atomically publishes the result
and independently verifies result plus original-input checksums. If `--phase all` is called
without approval it runs only through `review_ready`, reports the pause and exits cleanly.

## Recovery and safety rules

- Use `--dry-run` before the first controller or preparation submission.
- Use a new work directory for a scientifically changed campaign.
- Use `--resume` only for an unchanged configuration and completed checksum-valid result.
- Use `--unlock` only after confirming no controller or Snakemake process owns the workflow.
- Never edit files under a completed `result/`; any extra or changed file invalidates it.
- Preserve the workflow state, result manifest, scheduler logs and frozen configuration with
  the analysis record.
