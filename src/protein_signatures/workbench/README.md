# Protein research workbench

This package reads mounted DeepClust resources and writes new, independent analysis jobs.
It does not modify completed protein-signature campaigns or the global database.

- `config.py`: strict operator profiles, source identity, E3 path import and diagnostics.
- `catalogue.py`: lazy exact mapping, official column aliases and complete row-group retrieval.
- `analysis.py`: streaming FASTA indexing, explicit similarity controls, bounded tools and alignment.
- `jobs.py`: durable requests, exclusive workers, input/output checksums and interrupted-job recovery.
- `operations.py`: search, member screening, alignment, local Pfam, frozen signatures and curation starters.
- `seeds.py`: bounded seed export from canonical protein, HOG, label or feature selections.
- `cli.py`: named-option operator actions and worker entry point.

See `docs/DENBI_WORKBENCH.md` at the repository root for setup and scientific interpretation.
The Streamlit pages are implemented in `protein_signature_app/workbench.py` and
`protein_signature_app/workbench_viewers.py`. Focused tests live in `tests/workbench`.
