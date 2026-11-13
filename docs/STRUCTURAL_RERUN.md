# Compare the original structures with LoL-align

[Reifenrath et al. (2026), Nature Methods](https://doi.org/10.1038/s41592-026-03249-w)
introduced LoL-align for local **pairwise 3D structure alignment**. It is implemented
as alignment type 3 in [Foldseek](https://github.com/steineggerlab/foldseek#alignment-mode).
MAFFT alignments and FastTree trees still use amino-acid sequences. The original
statistical analysis, reference signatures, models and structural evidence stay intact.

## Install and verify Foldseek

On Linux x86-64 with AVX2, use an existing LoL-capable Foldseek, or install the
tested revision without sudo:

```bash
bash scripts/install_foldseek_lolalign.sh \
    --prefix /data/pete/software/foldseek_lolalign_39b02b9
```

Use the Sabrina software root `/home/ubuntu/data/pete/software` when installing there.
The installer checks revision `39b02b9fd8fea53a26a13ce638bfda616a0e2003` and the
archive SHA-256 before installing. The upstream binary URL is rolling. If its bytes
have changed, supply the exact archived build with `--archive`, or install another
verified LoL-capable revision. The installer does not replace an existing prefix.
The tested revision's [source is available upstream](https://github.com/steineggerlab/foldseek/tree/39b02b9fd8fea53a26a13ce638bfda616a0e2003).

Set this in the existing signature workbench TOML on Stefan:

```toml
foldseek = "/data/pete/software/foldseek_lolalign_39b02b9/bin/foldseek"
```

The worker checks that `foldseek easy-search -h` advertises
`--lolalign-multidomain` before comparing anything. An older Foldseek binary can
still serve the existing signature pipeline; it cannot run the new LoL route.

## Try one pair first

Under **Protein research workbench → Research tools → LoL-align structures**,
select two single-chain PDB or mmCIF models. You can upload models or reuse
completed ESMFold predictions. Choose the scoring mode, queue the job and inspect
the finished pair. The comparison publishes the aligned residue correspondence,
geometric metrics, separate LoL scores, raw Foldseek TSV and mobile-to-reference
transform. The trace overlay applies that transform. It does not infer pocket
conservation or add features to a completed signature analysis.

## Rerun the exact pairs from an earlier structural run

Use the **original E3 structural producer's** `tables/structural_alignments.parquet`
from the completed Stage 09b or standalone structural run. The reduced tables in
the signature result do not contain the original model paths and SHA-256 values.
The original coordinate files are required. If copied to Stefan, use the
recorded path prefix as `--source-model-root` and the current equivalent folder
as `--models-root`. The paths below are examples to replace with your actual
table, coordinate folder and group.

```bash
PS_BIN="/data/pete/software/miniforge3/envs/protein_signature_analysis/bin"
PS_PROFILE="/data/pete/protein_signature_research/config/workbench.toml"
PS_PRIOR_PAIRS="/path/to/original/structural_alignment/tables/structural_alignments.parquet"
PS_MODEL_ROOT="/data/pete/path/to/copied_original_models"
PS_OLD_MODEL_ROOT="/path/recorded/in/original/model/paths"
PS_GROUP="REPLACE_WITH_EXACT_CLUSTER_ID"
PS_RERUN="/data/pete/protein_signature_structural_reruns/lolalign_pilot_20261111"

"$PS_BIN/protein-signature-workbench" --config "$PS_PROFILE" \
    --action plan-structures --pair-table "$PS_PRIOR_PAIRS" \
    --cluster-id "$PS_GROUP" --max-pairs 100

"$PS_BIN/protein-signature-workbench" --config "$PS_PROFILE" \
    --action compare-structures --pair-table "$PS_PRIOR_PAIRS" \
    --models-root "$PS_MODEL_ROOT" --source-model-root "$PS_OLD_MODEL_ROOT" \
    --cluster-id "$PS_GROUP" --max-pairs 100 --output "$PS_RERUN"
```

If model paths have not moved, omit `--source-model-root`. For an extended run,
use a new output directory, repeat `--cluster-id` for selected groups, or omit
that flag for all prior US-align/TM-align pairs. Raise `--max-pairs` only after
reviewing the plan. Counts include separate rows for the original methods.
Repeated identical ordered model pairs reuse one new comparison. The batch
runs one pair at a time and does not duplicate original coordinate files in
the report directory. Run a long batch in an independent tmux session on the
server; it is a CLI command, separate from the app worker.

The batch refuses an existing output directory and checks that both models
match the original hashes before each comparison. Missing or changed models
are `NOT_ASSESSED`. A successful tool call with no correspondence is
`NO_ALIGNMENT_FOUND`; failed commands are `FAILED`. Those states have blank
unmeasured scores, never zero scores. `COMPARISON.json` binds the prior table,
report, executable revision and status counts. Upload
`structural_method_comparison.tsv` on the LoL-align page to filter by previous
method or new status. Keep the pair directories and provenance with the TSV.

## Read the comparison

| Field | Meaning |
| --- | --- |
| `previous_*` | Original exact model pair, method, status and recorded metrics. |
| `lol_normalised_score` | LoL score from the Foldseek `evalue` slot; **not a statistical E-value**. |
| `lol_homology_probability` | Tool-estimated homology probability; not a function or pocket probability. |
| `lol_query_tm_score`, `lol_target_tm_score` | Geometric scores on LoL's local correspondence; no global TM-align was run for these new rows. |
| `lol_query_model_coverage`, `lol_target_model_coverage` | Aligned residue pairs divided by each coordinate model length; not full FASTA coverage. |

The publication benchmarks support evaluating LoL-align, but do not establish
that it improves your specific E3 models. Start with models that span short
and multi-domain proteins, disagreements, and missing evidence. Examine local
coverage as well as scores before deciding whether to change any analysis
threshold. Pocket overlap or signature re-estimation requires a separately
planned analysis with the matching residue mapping and controls.

The provider's original E3 LoL alignment workflow remains available in
[E3_project_draft](https://github.com/peterthorpe5/E3_project_draft/blob/main/e3_end_to_end_workflow/docs/LOLALIGN.md).
This workbench update does not change or run that workflow.
