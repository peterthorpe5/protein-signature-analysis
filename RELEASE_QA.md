# Release QA: version 0.1.0

QA date: 12 September 2026  
Candidate: `protein-signature-analysis` 0.1.0  
Local platform: Linux x86_64, Python 3.12.14

## Release decision

The candidate is ready for Mac/HPC handover. Functional, numerical, file-format,
application, build and installation checks passed. The real top-200 and full top-1,000
structural campaign outputs were not supplied to this repository QA, so reproducing those
scientific results remains deliberately pending.

One independent release confirmation is still required after transfer to the Mac: run the
complete quality gate on the exact extracted tree and allow GitHub Actions to repeat it after
the first push. This is recorded explicitly because the final QA changes added 14 focused
regression cases after the most recent complete coverage run.

## Code-quality and test gates

| Gate | Result |
|---|---|
| Python compilation | PASS |
| Ruff formatting and lint | PASS |
| PEP 8 (`pycodestyle`, 100-character line limit) | PASS |
| Google-style docstrings (`pydocstyle`) | PASS |
| Bash syntax, including the macOS Bash 3.2 compatibility path | PASS |
| Complete pre-release baseline | 237 passed, 3 third-party warnings, 299.80 seconds |
| Baseline statement/branch coverage | 95.42% (threshold 95.00%) |
| Current collected suite | 251 cases |
| Current focused OrthoFinder adapter suite | 20 passed |
| Current shell/repository contract suite | 9 passed |
| Canonical app browser and pagination regression cases | 6 passed |
| External-file DuckDB view rejection | PASS for app and OrthoFinder resource adapters |

The three warnings come from SHAP's use of Matplotlib colour-map methods marked for future
deprecation. They are outside this package, do not change results and are retained rather
than hidden.

## Fresh offline end-to-end campaign

A new output and new SHAP/report caches were used; no earlier completed result was reused.
The packaged E3/F-box example validated and completed in 10 seconds.

| Check | Result |
|---|---|
| Input proteins | 24 |
| E3 profile | version 1.1.0; 118 labels; structural evidence required |
| Pfam/domain evidence | COMPLETE; 12 hits and 24 explicit assessments |
| Fold/structural evidence | COMPLETE; 8 fold assignments and 7 pairwise comparisons |
| Canonical tables | 26 |
| TSV = Parquet = DuckDB | PASS for all tables and all 3,629 rows |
| Manifested inputs/outputs | 6 / 237; all size and SHA-256 checks passed |
| Physical result files | 239; no links or undeclared files |
| Numbered analysis files | 183, indexed by 181 inventory rows |
| Owning/final table pairs | 26 / 10, each as TSV and formatted XLSX |
| Owning/final logical figures | 18 / 18 |
| Figure files | 36 PDF, 36 PNG and 36 SVG; all parsed successfully |
| Mandatory SHAP output | 18 assets representing 6 logical plots; validation partition |
| App canonical downloads | complete TSV/XLSX pair for every one of 26 tables |
| Result verification | PASS |
| Completed-result resume | PASS; checksum-compatible result reused |

All 37 Excel workbooks passed ZIP/package integrity and formatting checks. Values in 71,189
data cells matched their TSV authorities, including typed numbers and missing values. Five
workbooks were intentionally header-only. No worksheet cell formula nodes were present.
Representative small, wide and report-inventory workbooks also opened and converted to PDF
with LibreOffice. The generated SHAP beeswarm and converted inventory were inspected
visually.

## E3 seed-catalogue conversion

The supplied 1,000-row `e3_seed_catalogue.tsv` was converted twice into independent starter
directories. The five outputs were byte-identical between runs.

| Check | Result |
|---|---|
| Source SHA-256 | `10945b7cae2212f3bc425bd80a37481c96f3194742724bdabe41c0148e73db52` |
| Source dimensions | 1,000 unique identifiers; 32 columns |
| Sequence content | 409,431 amino acids; lengths 52-3,758 |
| Proposed-category strings | 127 distinct values retained for review |
| Conservative assignment state | all `UNMAPPED`; none training eligible |
| Promotion state | all `NOT_PROMOTED_REVIEW_REQUIRED` |
| AlphaFold starter mappings | 1,000 syntactic accession candidates |
| Duplicate-run reproducibility | PASS, byte for byte |

The proposed categories are deliberately not treated as reviewed E3 truth. A curator must
map each starter row to the controlled profile before association or model training.

## Independent numerical check

The native two-sided Fisher exact implementation was compared with SciPy for all 2,400
non-empty 2x2 tables whose cells ranged from zero through six. The largest absolute p-value
difference was `9.658940314238862e-15` (ours `0.9999999999999903`, SciPy `1.0`), comfortably
below `1e-12`.

## Distribution and installation checks

The source distribution and universal wheel built offline with setuptools 84.0.0 and wheel
0.48.0. Archive paths were traversal-safe; no links, devices, caches, Git internals, virtual
environment files or common secret/key filenames were present. All four shell launchers
retained executable mode. Every wheel `RECORD` size and SHA-256 entry verified.

The wheel and source distribution were installed separately into clean virtual environments
without resolving dependencies. In both installations:

- `protein-signatures --version` returned `0.1.0`;
- the imported module came from the new environment, not the working checkout;
- the bundled E3 profile contained 118 labels and generated 73 default comparisons;
- `protein-signatures describe-profile --profile e3` succeeded; and
- `protein-signature-app --help` succeeded.

## Required independent Mac confirmation

Run this from the extracted repository on the Mac before creating the release tag:

```bash
cd /Users/PThorpe001/github_repos/protein-signature-analysis
conda env create --file environment.yml
conda activate protein_signature_analysis
python -m pip install --no-deps --editable '.[app,dev]'
./run_tests.sh | tee release_qa_macos.log
```

If the dedicated environment already exists, replace `conda env create` with:

```bash
conda env update --file environment.yml --name protein_signature_analysis --prune
```

The expected outcome is 251 passed tests, only the three SHAP/Matplotlib warnings above and
at least 95.00% total coverage. Do not tag the release if the count, coverage threshold or
any test differs; retain the log and report the complete output.

## Scientific and operational boundaries

- The actual top-200 precursor output and full top-1,000 structural run were not available.
  Adapter contracts are tested, but empirical reproduction awaits those completed files.
- Native sequence discovery in 0.1.0 is exact overlapping amino-acid k-mers. MEME, STREME,
  FIMO, HMMER, conservation and disorder evidence can be imported but are not invoked.
- A formally named fold association requires supplied/imported SCOP, CATH, ECOD or another
  versioned fold authority. Native Foldseek clusters are structural neighbourhoods, not
  automatic formal fold classifications.
- The local container did not contain Foldseek or Conda and did not make live AlphaFold DB
  requests. Command construction, parsing, caching, acquisition states and launch behaviour
  are covered with fixtures/mocks; live cluster and provider smoke tests remain external.
- OrthoFinder 2.5.5 and 3.x raw layouts and the preferred published-resource contract are
  tested with synthetic fixtures. This package consumes OrthoFinder results and never runs
  OrthoFinder.
- Fisher/BH results are observational associations. Separate local feature-family and
  study-wide FDR values do not remove phylogenetic, curation or control-selection bias.
- Complete-case, group-aware elastic-net models can lose power when evidence is missing.
  SHAP explains the fitted model in its declared partition; it is not a mechanism or causal
  explanation.
- Dynamic Plotly-to-PDF downloads require Chrome/Chromium through Kaleido. Pipeline PDFs and
  native SHAP PDFs do not. The app limits one complete downloadable asset to 512 MiB.
- XLSX exports disable formula and URL conversion. TSV remains lossless machine-oriented text
  and may begin with spreadsheet-active characters; use the hardened XLSX for manual opening.
- Inputs should be immutable and in access-controlled local storage for the duration of a
  run. Checksums detect ordinary pre/post-run changes but are not a transactional snapshot of
  multiple independently changing files.
- Do not launch concurrent writers against the same result or catalogue destination. Atomic
  staging protects normal completion, but version 0.1.0 does not claim protection from a
  hostile same-account process bypassing the launch protocol.
- Dependency ranges support reproducible provenance records but are not a fully hash-locked
  software supply-chain environment. Preserve `run_metadata.json` and the manifest with each
  scientific result.

## Persistence checkpoint

Before final QA, a clean source checkpoint was written to persistent storage as
`protein-signature-analysis-QA-checkpoint-20260912.tar.gz` with SHA-256
`31604dcc8ce451bebe1fa7ce57fe10a66d652f4e98f0f2751760df248bfa4798`. The final release
archive is a separate, newer artefact and should be preferred after its checksum verifies.
