# Workbench verification

Run `python -m pytest tests/workbench -q` from the repository root.

Fixtures use tiny real Parquet files matching the official DeepClust schema, including a
family crossing physical row-group boundaries, identical sequences with different IDs,
unindexed clusters, reviewed taxonomy and ambiguous identifier mappings.

Tests exercise real DuckDB queries, streaming input, durable requests, hashes, CLI actions,
headless Streamlit navigation, phylograms, frozen-signature scoring and source protection.
DIAMOND, MAFFT, FastTree and HMMER outputs are mocked at the executable boundary; the
local SSH launcher uses fake SSH/curl commands and never connects to Sabrina. Process
execution, timeout and termination tests also use short local subprocesses.

The global production DeepClust release and live de.NBI services are not available in
these tests. Use `--action doctor` and start with one small seed family after installation.
