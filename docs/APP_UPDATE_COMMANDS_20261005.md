# Apply the app update and reopen a completed result

5 October 2026. Repository: `peterthorpe5/protein-signature-analysis`.
Patch base: `41f8ffe9ea010c7c9d96e6b2b8bb266596ed2549`.

This update adds help and improves presentation of existing results. A new cluster
analysis is **not required**. If your completed result is already on the Mac, follow
steps 1 and 2, then step 5. Steps 3 and 4 cover updating the cluster checkout,
verifying a completed result and downloading it when needed.

The bundle contains `APPLY_UPDATE.sh`, `overlay/`, `changes.patch`, `BASE_COMMIT.txt`, this guide,
`REVIEW.md` and `verification/`. It does not depend on a `README_APPLY.txt` file.
No data files, credentials or Git metadata are copied by the overlay.

## 1. Mac: extract and apply the compatible update

This ZIP targets commit `41f8ffe9ea010c7c9d96e6b2b8bb266596ed2549`
(`back up`, pushed at 17:01 UK time on 5 October). It combines the code already
pushed to GitHub with the remaining help and chart update. Keep the original ZIP
and recovery ZIP, but use this compatible ZIP for installation.

The earlier `git apply --check` failed before applying its patch. The later
chained command stopped when its base-commit comparison was false; it did not
reach rsync. This new script prints the expected and actual commits if they differ.

Download and extract `protein_signature_app_compatible_update_20261005.zip`
into Downloads. Finder may do the extraction automatically. For Terminal extraction:

```bash
unzip -n "$HOME/Downloads/protein_signature_app_compatible_update_20261005.zip" \
  -d "$HOME/Downloads"
```

Then run this block. It explicitly replaces the old `SIGNATURE_UPDATE` variable.

```bash
SIGNATURE_UPDATE="$HOME/Downloads/protein_signature_app_compatible_update_20261005"
SIGNATURE_MAC_REPO="/Users/PThorpe001/github_repos/protein-signature-analysis"

bash "$SIGNATURE_UPDATE/APPLY_UPDATE.sh" --repo "$SIGNATURE_MAC_REPO"
```

The script checks the bundle checksums, requires a clean `main` checkout, saves a
local Git-history bundle and a copy of working files beside your checkout, then
pulls `main` using `--ff-only`. Ignored data and environments are left in place.
It checks the exact commit and patch before running the rsync dry run and copy.
The backup path is printed. The working-file copy includes tracked files and
untracked files visible to Git; the script stops if such local changes are present.

The included script contains the requested rsync commands:

```bash
git apply --check "$SIGNATURE_UPDATE/changes.patch"
rsync -avc --dry-run "$SIGNATURE_UPDATE/overlay/" "$SIGNATURE_MAC_REPO/"
rsync -avc "$SIGNATURE_UPDATE/overlay/" "$SIGNATURE_MAC_REPO/"
```

These lines document what the script already runs. Continue to step 2 after its
`APPLIED` message. If it stops, retain the printed output and backup; the error
identifies the reason. The script does not reset commits, delete files, stash work,
commit changes or push. There is no `--delete` rsync option.

## 2. Mac: install, check, commit and push

Use the existing `protein_signature_analysis` Conda environment. The update adds
no dependencies. An editable installation ensures the launcher uses this checkout.

```bash
cd "/Users/PThorpe001/github_repos/protein-signature-analysis"
conda run --no-capture-output --name protein_signature_analysis \
  python -m pip install --no-deps --editable '.[app,dev]'

conda run --no-capture-output --name protein_signature_analysis \
  python -m pytest -q \
  tests/test_help_content.py tests/test_exports.py tests/test_app_ui.py

conda run --no-capture-output --name protein_signature_analysis \
  python -m pycodestyle src tests --max-line-length=100
conda run --no-capture-output --name protein_signature_analysis \
  python -m pydocstyle --convention=google --add-ignore=D105,D107,D202 \
  src/protein_signatures src/protein_signature_app
conda run --no-capture-output --name protein_signature_analysis \
  python -m ruff format --check src tests
conda run --no-capture-output --name protein_signature_analysis \
  python -m ruff check src tests
git diff --check
git diff --stat
```

If development tools are missing, install the declared development extras in that
environment: `python -m pip install --editable '.[app,dev]'` through `conda run`.
The optional full repository gate remains:

```bash
conda run --no-capture-output --name protein_signature_analysis bash run_tests.sh
```

The original review base already fell below the 95% coverage gate. The combined
update passes all 556 tests, but reports 88.85% project coverage, so that gate
still fails. The threshold is retained. `REVIEW.md` and the verification logs
record the result; a coverage-gate failure is separate from failing tests.

After inspecting the changes and results, commit the explicit update files and
push to your existing `main` branch:

```bash
git add README.md \
  docs/APP_CRITICAL_REVIEW_20261005.md docs/APP_UPDATE_COMMANDS_20261005.md \
  docs/METHODS.md docs/TEST_TRACEABILITY.tsv \
  src/protein_signature_app/app.py src/protein_signature_app/help_content.py \
  src/protein_signatures/exports.py src/protein_signatures/result_help.py \
  tests/test_app_ui.py tests/test_exports.py tests/test_help_content.py
git commit -m "Add complete app help and readable evidence charts"
git push origin main
```

## 3. Cluster: pull, install and verify the result

Run these commands in your existing HPC SSH session after pushing the Mac commit.
The last handover recorded k=3,4,5 trial job `434616` as running; its current state
has not been checked here. Inspect it before changing a checkout/environment used
by a running job. These checks do not submit another job.

```bash
squeue --jobs 434616
sacct --jobs 434616 --format=JobID,JobName,State,ExitCode,Elapsed --parsable2
```

If that job is still running or pending, let it finish before updating its shared
checkout or environment. The default result below is the completed caliper
campaign shown in your screenshots. To use a completed and verified k5 trial,
change `SIGNATURE_CAMPAIGN` to `e3_all1972_k5_trial_20261002`.

```bash
SIGNATURE_HPC_REPO="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/protein-signature-analysis"
SIGNATURE_HPC_BASE="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs"
SIGNATURE_CAMPAIGN="e3_all1972_caliper_complete_20260928"
SIGNATURE_HPC_RESULT="$SIGNATURE_HPC_BASE/$SIGNATURE_CAMPAIGN/result"

cd "$SIGNATURE_HPC_REPO"
git status --short
git pull --ff-only origin main
git log -1 --oneline
conda run --no-capture-output --name protein_signature_analysis \
  python -m pip install --no-deps --editable '.[app,dev]'
conda run --no-capture-output --name protein_signature_analysis \
  protein-signatures verify --resource "$SIGNATURE_HPC_RESULT"
```

Continue with a complete result only when verification reports `VALID`. A pending
or failed trial is not a published result. No new `sbatch` command is needed for
this app update, and existing published workbooks retain their recorded definitions
and checksums. Workbooks newly exported in the updated app use the new definitions.

If you deliberately need to resume an incomplete campaign, use its original
workflow and inspect its existing configuration and logs. The separate
[k-mer trial guide](https://github.com/peterthorpe5/protein-signature-analysis/blob/41f8ffe9ea010c7c9d96e6b2b8bb266596ed2549/docs/KMER_LENGTH_TRIAL.md) describes that scientific workflow;
it is not a prerequisite for installing this app update. New submissions for this
setup should use `--slurm-account barton --slurm-partition general`.

## 4. Mac: copy a verified complete result from the cluster

Run this on the **Mac**, not inside the cluster SSH session. Substitute the HPC
login hostname or SSH alias you normally use for `YOUR_HPC_LOGIN_HOST`; that
externally accessible host is not recorded in the available handover.
The campaign value must match the result you verified in step 3.

```bash
SIGNATURE_HPC_LOGIN="pthorpe001@YOUR_HPC_LOGIN_HOST"
SIGNATURE_HPC_BASE="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs"
SIGNATURE_CAMPAIGN="e3_all1972_caliper_complete_20260928"
SIGNATURE_HPC_RESULT="$SIGNATURE_HPC_BASE/$SIGNATURE_CAMPAIGN/result"
SIGNATURE_MAC_RESULT="$HOME/protein_signature_results/$SIGNATURE_CAMPAIGN/result"

mkdir -p "$SIGNATURE_MAC_RESULT"
rsync -avh --partial --progress \
  "$SIGNATURE_HPC_LOGIN:$SIGNATURE_HPC_RESULT/" "$SIGNATURE_MAC_RESULT/"

conda run --no-capture-output --name protein_signature_analysis \
  protein-signatures verify --resource "$SIGNATURE_MAC_RESULT"
```

The complete result directory is transferred, including its tables, metadata,
figures and published model assets. `--delete` is not used. Failed transfers can
be resumed with the same command. The app verifies the copied result again.

## 5. Mac: launch the updated app

Stop the existing protein-signature app with Control-C in its Terminal. Use port
8504 as shown in your supplied screenshots; the wrapper also accepts another
free port. The direct launcher below logs to a dated file and the Terminal.

For the completed result copied using step 4:

```bash
SIGNATURE_CAMPAIGN="e3_all1972_caliper_complete_20260928"
SIGNATURE_MAC_RESULT="$HOME/protein_signature_results/$SIGNATURE_CAMPAIGN/result"
SIGNATURE_MAC_LOG_DIR="$HOME/protein_signature_app_logs"

mkdir -p "$SIGNATURE_MAC_LOG_DIR"
cd "/Users/PThorpe001/github_repos/protein-signature-analysis"
conda run --no-capture-output --name protein_signature_analysis \
  protein-signature-app --resource "$SIGNATURE_MAC_RESULT" \
  --port 8504 --address localhost \
  2>&1 | tee "$SIGNATURE_MAC_LOG_DIR/app_$(date -u +%Y%m%dT%H%M%SZ).log"
```

In another Mac Terminal, open the app:

```bash
open "http://localhost:8504"
```

If your existing verified result lives elsewhere, set `SIGNATURE_MAC_RESULT` to
the exact directory you already pass to `--resource`; no fresh transfer is needed.

The repository wrapper provides the same launch without the extra log capture:

```bash
./run_protein_signature_app.sh --resource "$SIGNATURE_MAC_RESULT" \
  --conda-environment protein_signature_analysis --port 8504 --address localhost
```

Open each page's question-mark panels, hover over table headings, inspect the
recorded limits, and use the searchable glossary. In Signature explorer, change
the q-value axis. In Explainable prediction, compare the individual validation
SHAP chart with the preserved published graphics. The completed campaign's
calculations and statistical decisions are unchanged by the display update.
