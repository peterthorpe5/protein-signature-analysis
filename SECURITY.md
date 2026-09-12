# Security policy

## Supported version

Security fixes are applied to the latest released minor version. Scientific results remain
immutable; rerun a campaign with the fixed release rather than editing a published bundle.

## Reporting a vulnerability

Use GitHub's private **Security → Report a vulnerability** workflow for
`peterthorpe5/protein-signature-analysis`. Do not open a public issue containing an exploit,
credential, private dataset path or embargoed scientific input.

Include the affected version, operating system, minimal reproduction, impact and whether the
problem involves path traversal, archive handling, untrusted TSV/YAML, DuckDB queries,
external commands, network retrieval or the Streamlit app.

## Security model

- Result and upstream-resource paths are resolved and checked against path escape.
- Streamlit and published OrthoFinder-resource queries use locked, read-only DuckDB
  connections with external file access, extension loading and persistent secrets disabled.
- Foldseek is invoked without a shell.
- AlphaFold model downloads require HTTPS and an allowlisted EBI host.
- Inputs and published outputs are checksum-inventoried.
- Excel formula and URL conversion is disabled for untrusted cell text.
- Existing outputs are not overwritten; publication uses staging and an atomic rename.

These controls do not make arbitrary third-party input safe to trust. Run unreviewed
campaigns with ordinary user privileges, inspect manifests before opening external files and
do not place credentials in configuration, TSV, FASTA descriptions or uploaded archives.
Keep scientific inputs immutable and access-controlled for the complete run, and do not start
concurrent writers against the same result or catalogue destination. Open formatted XLSX
exports for manual spreadsheet review; raw TSV intentionally preserves leading characters
exactly and is not a spreadsheet-sanitised representation.
