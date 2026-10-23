# Protein research workbench on de.NBI

Version 0.3.0 adds a generic family-expansion workbench alongside the completed-result
viewer. E3 is an example profile, not a requirement. Your Mac supplies the browser and
SSH tunnel. Sabrina reads its installed DeepClust resources and runs the analyses.

## Install the update

Keep running production jobs on their existing code/environment until they finish.
The commands below update your checkout only when you deliberately run them.

On your Mac, unpack the supplied archive in Downloads. It contains a clean source
tree without a `.git` directory. Set the destination to your existing checkout:

```bash
PS_UPDATE="$HOME/Downloads/protein_signature_analysis_v0_3_0/protein-signature-analysis"
PS_REPO="$HOME/path/to/protein-signature-analysis"

git -C "$PS_REPO" status --short
rsync -av --dry-run \
  --exclude '.git/' --exclude '.venv/' --exclude '__pycache__/' \
  "$PS_UPDATE/" "$PS_REPO/"
PS_BACKUP="$HOME/Downloads/protein_signature_backup_$(date +%Y%m%d_%H%M%S)"
rsync -av --backup --backup-dir="$PS_BACKUP" \
  --exclude '.git/' --exclude '.venv/' --exclude '__pycache__/' \
  "$PS_UPDATE/" "$PS_REPO/"

cd "$PS_REPO"
conda activate protein_signature_analysis
python -m pip install -e '.[app,dev]'
python -m pytest tests/workbench -q
python -m ruff check src tests
python -m ruff format --check src tests
git diff --stat
git add src tests configs docs pyproject.toml environment.yml environment.denbi.yml \
  launch_protein_signature_denbi.sh README.md CHANGELOG.md CITATION.cff RELEASE_QA.md \
  MANIFEST.in run_tests.sh .github/workflows/ci.yml
git diff --cached --stat
git commit -m "Add de.NBI DeepClust protein research workbench"
git push origin HEAD
```

`rsync` deliberately has no `--delete`; it preserves other local files. Check local
modifications before copying because changed files with the same names are updated.
The archive records its exact GitHub base commit in `BASE_REVISION.txt`.

On Sabrina, use the large data disk and the relocated Miniforge installation:

```bash
ssh denbi

PS_CONDA="/home/ubuntu/data/pete/software/miniforge3/bin/conda"
PS_SERVER_REPO="/home/ubuntu/data/pete/software/protein-signature-analysis"

mkdir -p /home/ubuntu/data/pete/software
git clone https://github.com/peterthorpe5/protein-signature-analysis.git "$PS_SERVER_REPO"
cd "$PS_SERVER_REPO"
"$PS_CONDA" env create --file environment.denbi.yml
"$PS_CONDA" run --no-capture-output --name protein_signature_analysis \
  python -m pip install -e '.[app,dev]'
```

If that checkout already exists, replace the clone command with
`git -C "$PS_SERVER_REPO" pull --ff-only`. If the environment already exists,
use `"$PS_CONDA" env update --name protein_signature_analysis --file environment.denbi.yml`
instead of creating it. No `--prune` is used.

Alternatively, after applying the update on your Mac, transfer source directly:

```bash
rsync -av --exclude '.git/' --exclude '.venv/' --exclude '__pycache__/' \
  "$PS_REPO/" denbi:/home/ubuntu/data/pete/software/protein-signature-analysis/
```

## Reuse the installed E3 database paths

The E3 launcher specifies the existing profile at
`/home/ubuntu/data/protein_research/config/research.toml`. The individual DeepClust
source paths are read from that file, rather than guessed from the server layout.

On Sabrina:

```bash
PS_BIN="/home/ubuntu/data/pete/software/miniforge3/envs/protein_signature_analysis/bin"
PS_PROFILE="/home/ubuntu/data/protein_signature_research/config/workbench.toml"

"$PS_BIN/protein-signature-workbench" \
  --config "$PS_PROFILE" --action configure-denbi \
  --source-config /home/ubuntu/data/protein_research/config/research.toml \
  --workspace /home/ubuntu/data/protein_signature_research

"$PS_BIN/protein-signature-workbench" --config "$PS_PROFILE" --action doctor
```

This imports supported paths, executable settings and limits into an independent profile
and workspace. It does not modify the E3 configuration or either completed result.
An existing destination profile is protected. Subsequent changes should be made by
editing the existing TOML deliberately, not repeatedly running configuration import.
After changing a profile/source, submit new jobs: old queued requests fail with an
explicit configuration-change message rather than switching data silently.

Optional executables configured as absolute paths in the E3 environment may remain
valid there. The new environment also installs DIAMOND, MAFFT, FastTree and HMMER;
set their fields to `diamond`, `mafft`, `FastTree` and `hmmscan` to use these copies.

The application supports a mounted official release with:

| Purpose | Official file | Supported columns |
| --- | --- | --- |
| Locate complete families | `clust_index_RowGroup.parquet` | `CLUSTER`, `NMEMBER`, `RowGroup` |
| Extract proteins | `joined_with_index_RowGroupFinal.parquet` | `f0` ID, `f1` sequence, `f2` cluster |
| Map seed accessions | `SeqIdMapClustId.parquet` | `SEQID`, `CLUSTERID` |

Other explicit aliases used by the E3 workbench are supported. The official extraction
index omits clusters below three members: a mapped unindexed cluster remains visible
but cannot be retrieved through this indexed route.

## Open from your Mac

```bash
cd "$PS_REPO"
bash launch_protein_signature_denbi.sh
```

Default Mac URL: `http://127.0.0.1:8510`; remote app: `127.0.0.1:8505`.
These ports and the `protein-signatures` tmux session are separate from the E3 launcher.
Open **Protein research workbench** when a reference result is also loaded. Without
a reference result, the guided workbench opens directly.

```bash
bash launch_protein_signature_denbi.sh --status
bash launch_protein_signature_denbi.sh --no-browser
bash launch_protein_signature_denbi.sh --stop-tunnel
bash launch_protein_signature_denbi.sh --restart-app
```

The launcher starts/reuses an app and one independent worker under tmux. It forwards
only loopback ports and retains server logs under
`/home/ubuntu/data/protein_signature_research/logs`. Worker jobs survive browser and
SSH disconnection. Re-running the launcher reconnects. Stopping a tunnel does not stop
the app or worker. Only the owned GUI is restarted by `--restart-app`.

## The guided analysis

1. **Choose seeds.** Upload or paste an unaligned protein FASTA, or select exact IDs,
   local HOGs, profile labels or feature members from a completed reference result.
   The complete seed selection must fit the limit: no hidden top-N sampling occurs.
2. **Find families.** Map identifiers for existing membership, or search the global
   representatives using DIAMOND. Inspect mapping status, hit statistics and family size.
3. **Retrieve sequences.** Select exact global centroid IDs. Complete extraction uses
   only indexed physical row groups; it rejects missing IDs, duplicates and count mismatches.
4. **Screen members.** Re-align each member against seeds. Identity, bilateral coverage
   and E-value filters all apply. Exact-sequence deduplication is optional and audited.
   Optional taxonomy filtering uses configured reviewed metadata and exact NCBI taxon IDs;
   descendant lineages and taxonomy inferred from names are not assumed.
5. **Investigate.** Align a bounded exact set with MAFFT, inspect conservation and optionally
   infer a WAG+gamma FastTree. Local Pfam scans use `hmmscan --cut_ga`. A completed reference
   result enables frozen-signature evidence checks; positive Pfam evidence can be supplied.
6. **Prepare a new campaign.** Combine seeds and candidates, keep one exact sequence with
   seed priority and retain every original ID in the provenance ledger. Review the label
   template, select controls and recompute homology blocks before running new statistics.

All sequence outputs are FASTA; tables are TSV. Jobs retain their requests, parameters,
input checksums, source references, command arguments, tool version probes, logs, output
manifest and status.
Global source identity uses paths, sizes and modification times; it is explicitly **not**
a full cryptographic checksum of the very large DeepClust release. Record the actual
release in the profile and keep its source files immutable.

## Representative search database

The full sequence Parquet and mapping index do not themselves provide a DIAMOND search
database. If the E3 profile already has `centroid_database`, it is imported. Otherwise set
that field to an installed `.dmnd` built from the same release's exact global centroids.

If only a global representative FASTA is available, set `centroid_fasta`, then use
**Start here → Queue representative database build**. On completion, set
`centroid_database` to that job's `representatives.dmnd` path. This is an explicitly
requested, potentially large operation. There is no automatic full-database extraction
or download. A missing representative database disables only sequence-based searches;
exact member mapping and complete indexed retrieval remain available.

Search defaults are 30% identity, 50% coverage of each protein, E-value ≤ 1e-5 and up
to 50 representative targets per seed, using DIAMOND's very-sensitive mode. These are
editable starting points, not universal biological cut-offs. Searches are heuristic
and bounded: the app does not promise every related protein across the tree of life.

## Optional reference result and Pfam

Set `resource` in the workbench profile to a complete copied signature result. Keep
`COMPLETED.json`, `manifest.json`, DuckDB, Parquet and all declared assets together;
the application verifies the immutable result. The workbench workspace must remain
outside that result tree. A compact review archive without canonical assets is insufficient.

For local domains, install the desired Pfam release, run `hmmpress` against Pfam-A.hmm
and set `pfam_hmm`. The workbench does not automatically download a domain database.
Frozen-signature checks assess exact k-mers natively and can use positive Pfam evidence.
Unassessed domains, architectures and structures remain unknown. A no-hit in a different
Pfam release is not silently converted into absence of a frozen signature.

## Bounds and recovery

The GUI previews at most 1000 table rows and offers browser downloads up to 100 MiB.
Larger files have a direct `rsync` instruction. Uploads are limited to 20 MiB; larger
families flow between server jobs without repeated browser uploads. Default retrieval
limit: 2,000,000 proteins. Default alignment limits: 1000 proteins and 5,000,000 cells.
An oversized alignment fails explicitly; select exact IDs and try a smaller set.

`memory_gb` is a DuckDB memory bound. DIAMOND uses a conservative block size derived
from it, not a hard operating-system memory limit. MAFFT/FastTree run only after the
alignment input bounds pass. Use server/container/cgroup limits for an absolute process cap.

Only one worker executes jobs per workspace. COMPLETE outputs are checksum-verified
before downstream reuse or download. A failed task retains partial output and diagnostics
but offers no completed artifacts. Waiting jobs can be cancelled when a worker is idle.
An interrupted RUNNING task is never automatically retried: inspect its recorded PID,
server logs and process state, then submit a new request. Do not edit an active request.
Upgrading code should be done after the worker is idle and deliberately restarted.

A RUNNING record blocks further claims after interruption. Use the job's explicit
**Mark a stopped, interrupted job as failed** button, or the command below. Recovery
requires the same host and verifies that both the recorded worker and any saved
external-tool process group have stopped. It retains all files and never retries a job.

```bash
"$PS_BIN/protein-signature-workbench" --config "$PS_PROFILE" \
  --action recover --job-id job_REPLACE_WITH_THE_EXACT_JOB_ID
```

After a later code update, once Jobs shows no RUNNING or QUEUED work, stop the idle
worker window and reconnect with the Mac launcher to load the new worker code:

```bash
"$PS_BIN/tmux" kill-window -t protein-signatures:=worker
```

## Scientific interpretation

DeepClust expands sequence coverage. It does not automatically supply reviewed functions,
independent validation, structures or negative controls. Homologues and paralogues may
share motifs while differing in function. Related sequences must not leak between
discovery and validation; recompute connected homology/independence blocks across the
combined collection. More proteins may still represent one independent family.

The handover's limitations, including the earlier provisional labels and control-matching
concern, are not corrected by adding sequences. Existing completed results remain the
reference; a revised study requires an explicitly new campaign with reviewed design.

## Focused verification

```bash
python -m pytest tests/workbench -q
python -m pytest tests/test_help_content.py -q
python -m ruff check src tests
python -m ruff format --check src tests
bash -n launch_protein_signature_denbi.sh
```

Upstream contracts:
[DeepClust data retrieval](https://github.com/drostlab/deepclust_dataretrieval),
[DIAMOND](https://github.com/bbuchfink/diamond), and the
[E3 research workbench](https://github.com/peterthorpe5/E3_project_draft/tree/main/e3_python_app).
