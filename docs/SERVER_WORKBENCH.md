# Run the protein signature workbench on Stefan or Sabrina

Version 0.5.0 uses the same application and worker on both hosts. Each host has its own
job workspace, configuration and SSH tunnel. Mounts and completed signature
results are inputs; all new jobs use a separate writable directory.

## Copy the update from your Mac

Unpack the supplied ZIP into Downloads. You can upgrade directly from 0.3.0. The archive identifies the exact base
commit and includes a patch. To update an existing checkout, first check your
current changes. The file list copies only the changes for this release.

```bash
PS_RELEASE="$HOME/Downloads/protein_signature_analysis_v0_5_0"
PS_REPO="$HOME/github_repos/protein-signature-analysis"
PS_BACKUP="$HOME/Downloads/protein_signature_backup_$(date +%Y%m%d_%H%M%S)"

git -C "$PS_REPO" status --short
rsync -avh --dry-run --files-from="$PS_RELEASE/CHANGED_PATHS.txt" \
    "$PS_RELEASE/protein-signature-analysis/" "$PS_REPO/"
rsync -avh --backup --backup-dir="$PS_BACKUP" \
    --files-from="$PS_RELEASE/CHANGED_PATHS.txt" \
    "$PS_RELEASE/protein-signature-analysis/" "$PS_REPO/"

cd "$PS_REPO"
conda activate protein_signature_analysis
python -m pip install -e '.[app,dev]'
python -m pytest tests/workbench -q
git diff --check
git add --pathspec-from-file="$PS_RELEASE/CHANGED_PATHS.txt"
git diff --cached --stat
git commit -m "Add HOG workflows and audited DeepClust research"
git push origin HEAD
```

The second rsync backs up files it replaces. If the checkout has diverged from
`BASE_REVISION.txt`, review local changes first or use the archive's patch
with `git apply --check` before applying it. Use **one** update route.

## Install on Stefan

Before replacing the installed code, let the current worker finish running and
queued work. From your Mac, copy the source to Stefan using the SSH alias you already use.
Set `PS_SSH` to your working `pete@host` if it differs from `stefan`.

```bash
PS_SSH=stefan
PS_SERVER_REPO="/data/pete/software/protein-signature-analysis"
ssh "$PS_SSH" "mkdir -p '$PS_SERVER_REPO'"
rsync -avh --exclude='.git/' --exclude='.venv/' --exclude='__pycache__/' \
    --exclude='*.egg-info/' --exclude='.pytest_cache/' --exclude='.coverage' \
    "$PS_REPO/" "$PS_SSH:$PS_SERVER_REPO/"
ssh "$PS_SSH"
```

On Stefan, use your existing Conda installation. This example uses the prefix
expected by the Mac launcher. An existing different prefix can be passed
with `--environment` when launching.

```bash
PS_SERVER_REPO="/data/pete/software/protein-signature-analysis"
PS_ENV="/data/pete/software/miniforge3/envs/protein_signature_analysis"
cd "$PS_SERVER_REPO"
if [[ -d "$PS_ENV/conda-meta" ]]; then
    conda env update --prefix "$PS_ENV" --file environment.denbi.yml
else
    conda env create --prefix "$PS_ENV" --file environment.denbi.yml
fi
```

The environment installs the editable app, DIAMOND 2.2.1 or newer, MAFFT, FastTree, HMMER
and tmux. The LoL-enabled Foldseek build is installed separately:
[structural comparison instructions](STRUCTURAL_RERUN.md).

## Attach the installed DeepClust data

On Stefan, the mount helper gives you `/data/pete/ken_deepclust_data`. Check
that it is available and that the actual files have arrived. The data link
should lead into the read-only SSHFS mount.

```bash
PS_DEEPCLUST="/data/pete/ken_deepclust_data"
find -L "$PS_DEEPCLUST" -maxdepth 3 -type f \( \
    -name 'clust_index_RowGroup.parquet' -o \
    -name 'joined_with_index_RowGroupFinal.parquet' -o \
    -name 'SeqIdMapClustId.parquet' -o -name 'centroids.dedup.dmnd' \) -print
```

Set `PS_DEEPCLUST` to the directory containing all three official Parquet
files. This command creates a **new** profile and references the mounted files;
it does not copy or build the database.

```bash
PS_BIN="$PS_ENV/bin"
PS_WORKSPACE="/data/pete/protein_signature_research"
PS_PROFILE="$PS_WORKSPACE/config/workbench.toml"

"$PS_BIN/protein-signature-workbench" --config "$PS_PROFILE" \
    --action configure-mount --deepclust-root "$PS_DEEPCLUST" \
    --workspace "$PS_WORKSPACE" --release "Ken's installed DeepClust release"

"$PS_BIN/protein-signature-workbench" --config "$PS_PROFILE" --action doctor
```

If the files are in separate directories, start with the
[example TOML](../configs/workbench.example.toml), set the three exact absolute
paths, and run doctor. On either host, an existing profile can be edited in
place after jobs finish. Do not rerun `configure-mount` over an existing file.

If the matching representative DIAMOND database is already mounted, set
`centroid_database` to its absolute `.dmnd` file, or supply
`--centroid-database` when creating the profile. The complete index and sequence
Parquet files support exact membership and indexed retrieval. The
representative database separately enables sequence-similarity search.
Searches return bounded candidate lists, with no promise to find every relative.

For Sabrina, retain the existing environment at
`/home/ubuntu/data/pete/software/miniforge3/envs/protein_signature_analysis`
and the existing profile at
`/home/ubuntu/data/protein_signature_research/config/workbench.toml`.
Use the same env update command with that prefix, and keep the existing
DeepClust paths. If a profile is needed from scratch, the
[Sabrina reference setup](DENBI_WORKBENCH.md#reuse-the-installed-e3-database-paths)
shows how to import it from the E3 research configuration.

## Optional completed result

To retain your all1972 reference view, copy the complete result once if it is
not on Stefan already. Run on your Mac:

```bash
PS_CAMPAIGN="e3_all1972_caliper_complete_20260928"
PS_RESULT_MAC="$HOME/protein_signature_results/$PS_CAMPAIGN/result"
PS_RESULT_STEFAN="/data/pete/protein_signature_results/$PS_CAMPAIGN/result"
ssh "$PS_SSH" "mkdir -p '$PS_RESULT_STEFAN'"
rsync -avh --progress --partial --exclude='.DS_Store' --exclude='._*' \
    "$PS_RESULT_MAC/" "$PS_SSH:$PS_RESULT_STEFAN/"
```

Set `resource` in the server TOML to that result directory. Copy
`COMPLETED.json`, the manifest and every declared asset together. The result
and the workbench must be separate directories. Leave `resource` empty
to open the workbench without a reference result.

## Launch from your Mac

```bash
cd "$PS_REPO"
bash launch_protein_signature_denbi.sh --server stefan
bash launch_protein_signature_denbi.sh --server stefan --status
```

Stefan opens at `http://127.0.0.1:8511`. Sabrina uses
`bash launch_protein_signature_denbi.sh --server sabrina` and port **8510**.
Both use a server loopback port of **8505** and their own local socket.
Use `--ssh-host`, `--environment` or `--config` to override defaults.

The `protein-signatures` tmux session runs the GUI and worker in separate
windows. Closing the browser or tunnel leaves server jobs running. Re-run
the same launch command to reconnect. `--stop-tunnel` disconnects the Mac.
`--restart-app` replaces only the GUI. After updating the package, stop the
**idle worker** on that server and restart the GUI from your Mac:

```bash
# On the chosen server, with PS_BIN set as above:
"$PS_BIN/tmux" kill-window -t protein-signatures:=worker

# On the Mac, using the same server and overrides as the original launch:
bash launch_protein_signature_denbi.sh --server stefan --restart-app
```

Completed old jobs remain viewable. Requests submitted under the old code
must finish before restarting; version-bound requests do not silently use
new code or settings.

## Research tools

Choose exact seed proteins, search the installed full release or map membership,
retrieve selected complete families and screen them against the original seeds.
Use **Research tools** or the direct links on completed jobs. The
[family research guide](FAMILY_RESEARCH.md) covers the new analysis and SSH options:

| Tool | Input and output |
| --- | --- |
| JBrowse MSA | Start MAFFT from a released HOG or open a saved alignment. Inspect exact identity counts, Newick and Auspice JSON, and download an offline viewer. |
| Audited subclustering | Stricter DIAMOND DeepClust plus exact centroid audits, complete input comparisons, species context, network and GraphML. |
| Nightingale | Published HOG members, completed Pfam/InterProScan results, or matching FASTA and portable annotation TSV. The exact ID is shown and positional predictions are tied to sequence SHA-256. |
| ESMFold & models | Small exact selection, opt-in API prediction, portable FASTA/request ZIP, and inspection of an imported model against its original FASTA. |
| Predict annotations | Local InterProScan, opt-in EBI service, or local eggNOG. Explicit per-protein assessments; eggNOG transfers do not invent residue positions. |
| LoL-align structures | Pairwise local correspondence, separate scores and rotatable overlay; previous/new method report filters. |
| Import DeepClust results | Bounded portable Parquet export, exact row/cluster selection and validated FASTA. Imported rows are not claimed as complete global membership. |

To enable the ESMFold API, set `esm_endpoint` to a compatible operator-approved
HTTPS service. For the Atlas plain-sequence route use
`https://api.esmatlas.com/foldSequence/v1/pdb/` and `esm_format = "sequence"`.
For a deployed compatible endpoint accepting a JSON `inputs` string, use
`esm_format = "json_inputs"`. Set `esm_token_env` only if the service needs a
bearer token and provide that environment variable to the **worker**. Keep the
token out of profiles, downloads and commands saved in shell history. The
default length limit is 400 residues, with no truncation. API availability and
model revision depend on the provider. No sequence is submitted by opening
a viewer or preparing a portable request.

For annotations, configure `interproscan` and its installed data. For the
optional EBI service, enter a per-job contact in the app or set an operator
default with `interpro_email`. Set `emapper` and
`eggnog_data` for a local eggNOG database. Pfam remains an independent local
HMMER route. A predicted hit does not establish function, and a no-hit is not
evidence that the protein lacks function.

The bundled JBrowse and Nightingale assets, sample alignment and display
adaptations come from the MIT-licensed
[E3 v0.23.1 reference](https://github.com/peterthorpe5/E3_project_draft/tree/245ffa9bdf068a5f91aa5658451ca1f02c549d2e/e3_python_app).
Their checksums and third-party notices ship with the Python package.
