# Changelog

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
