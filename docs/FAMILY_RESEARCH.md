# HOG and DeepClust family research

The workbench reads the mounted DeepClust release and the optional completed
signature result in place. Every new analysis has its own job directory.
Existing alignments, scores, structures and completed job manifests remain
available after upgrading. This release includes the
[LoL-align comparison and structural rerun tools](STRUCTURAL_RERUN.md).

1. From **Orthology & partitions**, search for a released HOG and choose
   **Align this HOG** or **Browse HOG features**. You can also open either
   tool from the workbench. A HOG selection identifies its OrthoFinder run,
   hierarchy node and group ID.
2. Choose all available released sequences or list an exact member subset.
   The UI counts members with missing sequences. Download the FASTA and
   member context, or run MAFFT. An optional exact reference controls the
   alignment statistics. It does not force MAFFT to preserve a reference
   alignment.
3. Use the HOG as seeds to search the installed global representatives, or
   map exact accessions. Retrieve the complete indexed global clusters that
   you select and screen their members against the seeds. Global DeepClust
   clusters and released HOGs have separate identifiers and provenance.
4. Send the selected HOG, retrieved family, screened candidates, or a saved
   alignment's exact unaligned input into **Audited subclustering**,
   **Predict annotations**, **ESMFold & models**, or MAFFT. A subcluster
   handoff uses its complete audited membership, even when a chart shows
   fewer nodes.

## Alignment statistics and trees

The JBrowse page shows a saved MAFFT alignment and its exact unaligned input.
Column and reference tables define canonical coverage, gaps, Shannon entropy
and the raw BLOSUM62 sum over comparable positions. Protein identity counts
only positions with a canonical amino acid in both rows. A pair with no
comparable residues has undefined identity. Download the pairwise TSV for
counts, or explicitly select rows for the bounded heatmap.

Use **Infer FastTree from this saved alignment** when a matrix is already
complete. FastTree reads that exact saved matrix without repeating MAFFT.
Download the Newick, exact tip mapping and Auspice v2 JSON. Branch lengths
describe divergence; the display root does not establish an ancestor.
A separately uploaded reviewed species Newick can be pruned to exact names
or reviewed taxon IDs. The resulting zero-length links for multiple proteins
from one species describe display grouping only.

## Audited subclusters and taxonomy

Install DIAMOND **2.2.1 or newer**. The job builds a temporary database for
the selected family, runs `deepclust`, then runs `realign` on the native
assignments. It retains the native TSV, exact alignments, tool ID map, final
membership TSV/Parquet, thresholds and logs. The audit derives identity
from identical residues and aligned length, checks coverage in both proteins
and verifies the requested E-value and a bit score greater than 20. Missing
or failing assignments stay in the family as singletons. Edges describe
passing centroid-to-member relationships. No all-pairs claim follows from
the network.

Choose a completed subcluster to compare its full membership with the whole
input: lengths, native and final group counts, taxonomy and exact evidence.
GraphML contains the displayed nodes and edges. Taxon Jaccard and exclusive
overlap counts use reviewed NCBI IDs from the complete selected groups.

`family_context.tsv` binds each metadata row to an exact protein ID and
normalised sequence SHA-256. Published species labels are separate from
reviewed names and IDs. Set `taxonomy` to a reviewed per-sequence mapping,
and optionally set `species_taxonomy` to an E3-compatible reviewed TSV such
as `e3_python_app/src/e3app/resources/species_taxonomy.tsv` from your
installed E3 checkout. Accepted species mapping columns are
`source_species_name`, `canonical_species_name`, `taxon_id`, or the
E3 review-template equivalents `workflow_species_label`,
`accepted_species_name`, `ncbi_taxon_id`. A supplied `mapping_status`
must be `REVIEWED` for a row to count. Missing or conflicting reviewed IDs
remain unknown. For older inputs, explicitly choose the relevant released
HOG to match IDs and sequences against the current result. The old job's
files stay as they were.

## Run jobs on Stefan, Sabrina or a configured SSH worker

Choose **Queue for worker** for the regular server worker, or **Run on app
host** to start just this job in the app's environment. Both share the
workspace lock. If it is busy, select the waiting request in **Jobs &
results** and start it once the active job finishes. You can also use the
displayed terminal command:

```bash
protein-signature-workbench --config "$PS_PROFILE" --action worker \
    --once --job-id "$PS_JOB_ID"
```

For audited subclustering on another server, install this same app version
and a supported DIAMOND there. Set these keys in the app host's profile:

```toml
subcluster_ssh_host = "stefan"
subcluster_ssh_executable = "/data/pete/software/miniforge3/envs/protein_signature_analysis/bin/protein-signature-workbench"
subcluster_ssh_profile = "/data/pete/protein_signature_research/config/subcluster-worker.toml"
```

Use an SSH alias that already works **from the app host** and has an
established host key. Use a separate writable workspace in the remote
profile, plus an absolute `diamond` path when Conda is not on the
non-interactive PATH. The app sends only the selected family FASTA and
thresholds after explicit approval in the UI, then checks the returned
inventory, input identity and artifact hashes. It does not copy the full
global release for this step.

EBI InterProScan accepts a private contact for an individual request, or the
profile's default. Hosted ESMFold and EBI each require consent before sending
sequences to the configured service. Viewing results and preparing a
portable request do not contact those services.

**Jobs & results → Workspace storage and scratch retention** shows the
workspace paths and previews eligible scratch. Apply cleanup after reviewing
the list. Running jobs, job requests, saved outputs and mounted data stay
untouched. Temporary files from older versions that lack an ownership marker
are outside this cleanup.

The family-analysis implementation is adapted from the MIT-licensed
[E3 Python app v0.23.10](https://github.com/peterthorpe5/E3_project_draft/tree/9900d41c25d3a613ec441a23003c63425ff85932/e3_python_app).
