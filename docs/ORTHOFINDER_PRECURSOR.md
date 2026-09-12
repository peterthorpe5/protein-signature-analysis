# OrthoFinder precursor boundary

## Ownership

This package consumes OrthoFinder evidence but never runs OrthoFinder. The standalone
[`orthofinder-results`](https://github.com/peterthorpe5/orthofinder-results) project remains
the preferred authority for generic parsing, taxonomy reconciliation, species/gene trees,
selection coverage and evolutionary-distance summaries.

`protein-signature-analysis` imports only the group relationships and context required to:

- keep homologous proteins in one discovery/validation block;
- expose exact group provenance beside a signature; and
- support group-aware model cross-validation.

It does not duplicate the OrthoFinder app or reinterpret its taxonomy tree.

## Preferred route: published resource

Configure:

```yaml
inputs:
  orthofinder:
    resource_dir: /absolute/path/completed_orthofinder_resource
    results_dir: null
    group_type: HOG
    hierarchy_node: N0
    run_id: ""
```

The adapter accepts `orthofinder-results` resource schema 3 or 4. It requires:

- `run_manifest.json` with `status: complete`;
- every declared output to exist with the stated size and SHA-256;
- `duckdb/orthofinder_results.duckdb`;
- one consistent `run_id` in manifest and database; and
- the required relations `distance_statistics`, `group_statistics`, `hog_memberships`,
  `legacy_orthogroup_memberships`, `resource_metadata` and `sequences`.

When `inputs.orthofinder.run_id` is non-empty it must equal the resource identity. The
adapter imports either one exact HOG node or the legacy orthogroup authority, plus available
group size/species/copy-number and evolutionary-distance context. It matches through the
precursor sequence authority and fails if the campaign FASTA cannot be mapped.

This route is preferred because the precursor has already reconciled the source layout and
published a checksum-complete authority.

## Fallback route: raw OrthoFinder output

Configure:

```yaml
inputs:
  orthofinder:
    resource_dir: null
    results_dir: /absolute/path/OrthoFinder/Results_campaign
    group_type: HOG
    hierarchy_node: N0
    run_id: orthofinder_campaign_001
```

The read-only adapter supports exactly OrthoFinder 2.5.5 and OrthoFinder 3.x layouts. It
recursively locates the shallowest `Log.txt`, extracts the recorded version, and requires
the exact `OrthoFinder run completed` marker emitted by the official OrthoFinder entry
point ([2.5.5 source](https://github.com/davidemms/OrthoFinder/blob/2.5.5/scripts_of/__main__.py#L1829),
[3.x source](https://github.com/davidemms/OrthoFinder/blob/e3e8d59b8162be5ba2b9e4a0a10f321e813a69b7/scripts_of/__main__.py#L1422)).
It can read:

- `Phylogenetic_Hierarchical_Orthogroups/N*.tsv`-style hierarchical tables using a selected
  `HOG` node; or
- `Orthogroups.tsv` using `LEGACY_ORTHOGROUP` and a blank hierarchy node.

`SequenceIDs.txt` is used when present to map OrthoFinder internal identifiers back to the
first exact FASTA-header token. `SpeciesIDs.txt` is recorded as provenance when present.
The parser accepts known heading spellings, preserves species-column labels and rejects a
campaign protein mapped to more than one group at the selected authority.

Raw input publishes memberships but not the richer precursor group-context relation. The
absence of that context is explicit; it is not back-filled from filenames or guessed
taxonomy.

This fallback validates the official completion marker and proves that the selected group
authority maps to the campaign FASTA. Independently retain the scheduler state and full
`Log.txt` as operational provenance. Prefer the published-resource route when available
because it additionally verifies a completed manifest, declared checksums and DuckDB schema.

## Composite identity

Never identify a HOG using `group_id` alone. The stable authority is:

```text
(run_id, group_type, hierarchy_node, group_id)
```

`legacy_orthogroup_id` and `gene_tree_parent_clade` are contextual fields. They are not
silently promoted to the selected group identity.

## Selecting HOG versus legacy orthogroup

- Prefer a defined HOG hierarchy node when the scientific question and precursor authority
  use hierarchical orthogroups.
- Use `LEGACY_ORTHOGROUP` only deliberately, with `hierarchy_node: ""`.
- Keep the choice fixed across discovery and validation.
- Changing group type, node or run ID changes the independence graph and therefore creates a
  scientifically different campaign.

## What is duplicated locally

The standalone package necessarily duplicates a small compatibility adapter so it can run
without importing another repository at runtime. It does **not** duplicate OrthoFinder
execution or the precursor's generic analytics. The adapter validates a narrow versioned
contract, converts memberships into this package's canonical schema and binds the exact
upstream files into the campaign manifest.

## Pre-flight checks

```bash
protein-signatures validate --config /path/campaign.yaml
protein-signatures run-all \
  --config /path/campaign.yaml \
  --output-dir /path/new_result \
  --threads 24
```

Validation parses all configured local evidence, reads the selected OrthoFinder authority
and requires the official raw-log completion marker. It does not acquire AlphaFold models or
run a Foldseek search; when Foldseek is enabled, it only resolves the executable, records its
version and validates the candidate/hit-limit plan. Retain independent scheduler records as
additional operational provenance before submitting the expensive campaign.

## Safe archive for inspection

If an unfamiliar precursor layout needs adapter work, package the exact completed resource
without modifying it:

```bash
SOURCE_DIR=/absolute/path/completed_orthofinder_resource
ARCHIVE=/absolute/path/orthofinder_resource_for_review.tar.gz

test -d "${SOURCE_DIR}"
tar -C "$(dirname "${SOURCE_DIR}")" \
  -czf "${ARCHIVE}" \
  "$(basename "${SOURCE_DIR}")"
sha256sum "${ARCHIVE}" > "${ARCHIVE}.sha256"
tar -tzf "${ARCHIVE}" | sed -n '1,200p'
```

On macOS, replace `sha256sum` with `shasum -a 256`. Do not package private credentials,
scheduler secrets or unrelated sequence authorities.
