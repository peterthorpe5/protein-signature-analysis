# Contributing

## Scientific changes first

Open an issue before changing a label hierarchy, feature meaning, FDR family, partition rule,
structural threshold interpretation or canonical schema. These are scientific-interface
changes even when the Python diff is small. State the biological question, expected
authority, missing-data behaviour and compatibility impact.

Do not silently broaden E3 labels, reclassify complex components as catalytic proteins or
promote proposed catalogue records to reviewed positives.

## Development setup

```bash
conda env create --file environment.yml
conda activate protein_signature_analysis
python -m pip install --no-deps --editable '.[app,dev]'
./run_tests.sh
```

The supported Python versions are 3.11 and 3.12. Foldseek is supplied by the Conda
environment. Plotly PDF downloads require Chrome/Chromium; `plotly_get_chrome` can install a
compatible runtime when policy permits.

## Code requirements

- PEP 8 with a 100-character line limit.
- Google-style docstrings for every module, class and function.
- Type hints and named arguments for workflow APIs.
- Deterministic ordering and explicit random seeds.
- Structured logging at stage and work-unit boundaries.
- Read-only authorities, checksum-bound caches and atomic output publication.
- Defensive validation with typed contextual errors.
- UTF-8 TSV rather than CSV for scientific exchange.
- UK English in user-facing text.
- No secrets, credentials, private paths or unlicensed reference data.

Every new function needs direct unit coverage of its normal and relevant failure states.
Update [`docs/TEST_TRACEABILITY.tsv`](docs/TEST_TRACEABILITY.tsv) when responsibility moves.
The repository gate compiles code, runs Ruff, pycodestyle and pydocstyle, executes unit,
integration and Streamlit tests with branch measurement, and requires at least 95% coverage.

## External tools and evidence

Construct commands as argument arrays and never through a shell string. Record executable
version, effective parameters, input checksums and output completion. Parser tests remain
mandatory even when a binary smoke test cannot run in CI.

An imported feature must name its true producer and version/reference. Do not imply that
this package calculated MEME, HMMER, conservation, disorder, fold or pocket evidence when it
only imported it.

## Pull requests

Keep changes focused and include:

1. the problem and scientific interpretation;
2. affected contracts and compatibility decision;
3. tests for success, malformed input, missing evidence and determinism;
4. documentation and example updates; and
5. the full `./run_tests.sh` result.

Never add a generated campaign result or large coordinate/model file to Git. Use tiny
synthetic fixtures with clear provenance.
