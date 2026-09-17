# Changelog

## Unreleased

- Represent universally assessed amino-acid k-mer denominators with one shared
  immutable campaign-protein universe instead of materialising one full set per
  feature. Validation now checks large positive-membership collections without
  copying them, preserving identical association denominators while preventing
  the all-species campaign's 20.46-billion-membership allocation.
- Ignore only recognised macOS `._*` AppleDouble sidecars and `.DS_Store` files
  when verifying an evidence-label bundle, log every ignored path and continue
  to reject every other undeclared file.
- Derive the completed-E3 Snakemake `e3_memory_mb` resource from the submitted
  `--slurm-memory` request; direct execution can set the same recorded resource
  explicitly with `--memory-mb`.
- Replace structural-cluster membership's repeated full-universe scans with
  bounded edge indexes while preserving discovery components, held-out
  projection, support counts, best scores, cluster identities and deterministic
  output. Add progress logs for the formerly silent post-Foldseek stage.
- Add allocation-time Slurm scratch discovery with explicit, `SLURM_TMPDIR`,
  `TMPDIR`, node `/tmp` and persistent-work-directory fallbacks; validate real
  writability and free space, export `TMPDIR`/`TMP`/`TEMP`, record provenance,
  clean successful job scratch and retain failed-job scratch for diagnosis.
- Make both evidence-led Snakemake DAGs compatible with atomic bundle
  publication by explicitly allowing removal of the empty output directory that
  Snakemake creates before the rule command. Existing files, non-empty
  directories and symbolic links remain protected.
- Import the structural precursor's explicit `REFERENCE` diagonal rows as
  assessment-universe membership rather than invalid pairwise comparisons. The
  rows remain checksum-bound, are validated against their group summaries and
  sentinel metrics, are counted in run metadata, and never enter association or
  structural-clustering evidence.
- Cross-check imported structural universes against `aligned_accession_count`
  rather than `selected_accession_count`, because selected proteins without
  usable coordinates were never members of the structural comparison universe.
- Reconcile raw OrthoFinder group members with campaign FASTA accessions using
  exact identifiers, controlled UniProt accession/entry aliases and the
  documented single `sample@@identifier` qualifier. Ambiguous mappings fail
  closed rather than selecting one identifier silently.
- Add regression coverage for the completed-E3 prepared-accession versus raw
  OrthoFinder identifier contract used by automated smoke and evidence-led runs.

## 0.1.0 - 2026-09-11

- Establish the standalone, protein-agnostic analysis engine and read-only
  Protein Signature Explorer.
- Add an extensible E3 ubiquitin-ligase profile covering mechanistic classes,
  complex-component roles, major subclasses and explicit unresolved states.
- Add a checksum-complete published-resource adapter plus raw-layout adapters tested with
  exactly OrthoFinder 2.5.5 and OrthoFinder 3.x; raw fallback requires OrthoFinder's official
  completion-log marker and the package never launches OrthoFinder.
- Add amino-acid k-mer, domain, fold, structure-cluster and pocket-feature
  association analyses with discovery and held-out group partitions. K-mer vocabularies are
  defined only on discovery sequences and projected unchanged onto validation proteins.
- Add machine-checkable generic-feature definition/cohort provenance, explicit per-feature
  assessment universes and an audit ledger; supervised-discovery and all-data-derived
  definitions have honest controlled audit-only scopes.
- Add typed TSV/Parquet publication, a physical DuckDB, checksums, validation
  records and an immutable run manifest.
- Add deterministic Foldseek command construction and result parsing for
  optional fold-library or structure-cluster evidence, with executable/version and
  candidate-universe pre-flight checks that do not run a search; `maximum_hits` cannot
  truncate all-versus-all campaigns.
- Separate retained structure provenance from controlled analysis eligibility; low-confidence,
  confidence-unavailable and sequence-unverified AlphaFold downloads cannot enter Foldseek or
  structure features.
- Add controlled fold/comparison states, strict payload consistency, explicit comparison
  universes and coverage scopes, cluster-local predecessor pockets, and complete structural
  assessment-universe ledgers.
- Add mandatory group-aware elastic-net modelling, held-out evaluation, numeric
  SHAP explanations, and native SHAP PNG/SVG/PDF graphics for every fitted model.
- Add numbered file-first analysis reports with TSV/formatted-XLSX table pairs,
  PNG/SVG/PDF figures, a report inventory and `99_final_results` mirror.
- Add a canonical-data app page with complete, manifested TSV/XLSX downloads for all 26
  result tables, alongside PDF downloads for every app plot.
- Lock app and published OrthoFinder-resource DuckDB connections against external file,
  extension and persistent-secret access.
- Add a Python 3.12 macOS launcher smoke job and tracked secret/transient-filename gate.
- Add a checksum-verifying `protein-signatures verify` command.
- Add a review-before-run shell launcher with safe existing-config semantics and Conda
  synchronisation, a synthetic E3 example and comprehensive tests.
- Add a checksum-gated bridge and phased launcher for completed E3 end-to-end workflow runs,
  converting exact Stage 05 sequences, Stage 06 Pfam assessments and Stage 09 AlphaFold assets
  while importing Stage 09b structural evidence and refusing automatic label promotion; Stage 09
  mean pLDDT now gates model eligibility without hiding low-confidence or unassessed coordinates.
- Allow campaign initialisation to set Foldseek's complete-search hit limit explicitly.
- Accept the completed end-to-end Stage 09b aggregate `datasets` manifest while retaining
  strict cross-checks against the outer stage checksum inventory; standalone structural
  component manifests remain supported.
- Accept the predecessor's intentionally log-less reused Stage 04 OrthoFinder publication
  only when its complete stage manifest, reviewed-archive authority, exact five-file
  validation ledger and current file checksums all agree; direct raw results still require
  OrthoFinder's official completion marker.
- Preserve Stage 05 composite DeepClust cluster identifiers as bounded, TSV-safe provenance
  and serialise their review-cell collections as deterministic JSON arrays; canonical HOG,
  orthogroup and protein identifiers retain their stricter validation contracts.
- Stream only projected predecessor Parquet columns in bounded DuckDB batches and publish the
  prepared FASTA incrementally, avoiding the former whole-table and whole-FASTA duplicate
  allocations during completed-E3 preparation.
- Add validated direct Slurm submission for completed-E3 prepare, initialise, full-DAG and
  run phases, non-recursive workers, dry-run output and persistent job-specific stdout/stderr
  paths; accept scheduler allocations that safely exceed the requested CPU count.
- Add a generic Snakemake 9 workflow with explicit validation, atomic analysis and independent
  result/input verification rules; provide local and Slurm-executor profiles plus a durable,
  lock-protected Slurm controller patterned after the E3 production workflow.
- Add a completed-E3 Snakemake DAG spanning preparation, no-clobber review staging,
  checksum-bound curator approval verification, campaign initialisation, atomic analysis and
  independent final verification. The all-`UNMAPPED` template cannot be approved, and every
  populated analysis target must have its profile-resolved matched background.
- Add an isolated, fully automated software-smoke route that deterministically exercises all
  profile-default comparisons without human label review. It uses structure-eligible proteins,
  keeps whole HOG/exact-sequence/redundancy blocks exclusive to one label, checksum-binds
  synthetic labels to test-only approval and prohibits scientific interpretation.
- Make the synthetic-label generator profile-agnostic for arbitrary FASTA, custom profile,
  structure, OrthoFinder 2.5.5/3, published OrthoFinder-resource and redundancy inputs; the E3
  default covers all 73 comparisons (including F-box) through terminal targets and all 14
  profile-defined matched controls.
- Add `gawk`, the strict `nodefaults` channel, Snakemake and its Slurm executor plugin to the
  Conda environment; retain Kaleido as a pip-installed dependency for channel compatibility.
- Add generic, versioned evidence-led production labelling with trusted exact annotations,
  Pfam corroboration, conservative one-generation OrthoFinder propagation, explicit
  abstention and complete decision audits.
- Add outcome-blind control matching at joined homology/redundancy-block level by species,
  structure eligibility, length and domain complexity; exclude all target-like and upstream
  candidate proteins from the clean-control pool.
- Prevent label leakage by recording every label-defining Pfam feature and removing it from
  the confirmatory domain and SHAP input authority while retaining explicit assessment
  coverage.
- Add generic and completed-E3 Snakemake routes that stop after the evidence bundle by
  default, or run the complete provisional campaign only after an explicit acceptance flag.
- Publish label, matching, exclusion, class-coverage and unresolved-record audits as TSV,
  formatted XLSX, canonical Parquet/DuckDB tables, static figures and app downloads.
