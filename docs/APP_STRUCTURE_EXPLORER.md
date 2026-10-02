# Model and alignment explorer

The completed result already packages the database, figures, and any available coordinate
files under `result/assets/structures/`. Keep the **whole** `result/` directory together: the
application checks `COMPLETED.json`, `manifest.json`, and every declared file size and checksum on
opening. A lone DuckDB file is insufficient.
During an open session, file names, types, sizes and change times are checked on each interaction;
full checksum verification is repeated when that inventory changes. The explicit
`protein-signatures verify --resource /path/to/result` command always checks every byte.
Finder's `.DS_Store` and `._*` AppleDouble sidecars may appear anywhere under the copied
result without blocking the app. All declared outputs still need valid checksums, and other
unexpected files remain an error.

```bash
# On the Mac, from the repository after copying the completed result:
./run_protein_signature_app.sh --resource /absolute/path/to/completed/result
```

For the current HPC test, copy the result directory recursively after it has completed:

```bash
# On the Mac: replace HPC_LOGIN and MAC_DESTINATION with your own paths.
rsync -avP HPC_LOGIN:/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_caliper_complete_20260928/result/ MAC_DESTINATION/result/
```

Navigate to **Model & alignment explorer**. Select a protein and target-versus-background
comparison. The ranked mapped-region table, a 2D sequence heat map, and a rotatable 3D Cα
backbone show positional evidence. Recorded feature intervals and exact k-mer occurrences
require a positive target-minus-background prevalence difference and q-value ≤ 0.05 at the
chosen evidence tier before they colour a residue. Select discovery, held-out within-comparison validation,
or held-out study-wide validation. The stronger tiers require the corresponding validated
evidence class and q-value.
The comparison chooser shows completed positive significant signatures first and flags
comparisons with none. The suggested protein list contains target-class members with packaged
coordinates, ranked by model confidence; an exact protein ID box can inspect any other member.
The **All comparison outcomes** expander lists every configured comparison and its
published status. **Comparison for protein and model mapping** controls the positional
colours and alignment. The **Most significant association results** table has a separate
multi-select, initially including every comparison with positive significant signatures
(up to ten). Its balanced preview shows 20 rows per comparison by default, ranked by
discovery q-value then target prevalence difference; increase the per-comparison limit
when needed. The table and cross-class plot can be downloaded. The **Signature explorer**
provides the complete feature table for an individual comparison, including rows that
do not pass the positive significance filter.
Published k-mer memberships do not carry residue coordinates, so the viewer finds each exact
overlapping occurrence in the chosen protein sequence. These are **sequence occurrences of an
enriched feature**, not individually measured residue-level q-values. Whole-model structural
clusters still have no positional interval to project.
The positional view starts with the 50 most significant positive feature memberships for the
chosen evidence tier; the control can show 10 to 1,000. This prevents thousands of overlapping
short k-mers from making an entire sequence uniformly red. The ranked signature table above
remains separately available for broader inspection.
The **Structures & folds** page separately charts model eligibility and reports when named
fold assignments are absent; alignment-derived structural clusters remain browsable without
implying that they are curated fold names.
At overlapping positions, the smallest q-value wins. White means **no significant mapped
evidence**, not an assessed absence; blue marks weaker significant evidence, and red stronger
evidence. The colour scale saturates at q ≤ 10⁻⁸. The view identifies positional intervals
of enriched features; it does not assign residue-level p-values to individual residues.

The in-app renderer displays PDB and mmCIF Cα traces. A model gets positional colours only
when one model chain has
exactly the published full-length amino-acid sequence and residues numbered 1 through N
without insertion codes. Partial models, mismatches, and renumbered models are shown in
grey to prevent a false residue assignment. Coordinate assets stay in the verified bundle.
When no local model exists, a button can request the current AlphaFold DB PDB by a canonical
UniProt accession. The downloaded sequence and coordinate numbering must match exactly; the
model remains in the Streamlit session and can be downloaded separately. External UniProt,
InterPro, AlphaFold DB and Mol* links open only for canonical accessions. The Foldseek
search link lets you upload the downloaded PDB or mmCIF to search other structures. The
optional model fetch needs internet access from the machine running Streamlit; the rest
of the app uses local data.

The **Paired sequence alignment** panel selects a published structural comparison, shows its
TM-score and bilateral coverage, and computes an exploratory global **sequence** alignment
(+2 match, −1 substitution, −2 gap). Its alignment columns are coloured from the same
positional enrichment evidence for each protein, and the full aligned FASTA and alignment
table are downloadable. The result currently publishes **aggregate** structural comparisons
and whole-model clusters, not residue-pair alignments, so the app never labels these columns
as structurally superposed. The external EMERALD and RCSB pairwise links allow a further
alignment analysis; a given computed structure might be absent from RCSB.

## Optional pocket and residue-level annotations

The present structural-motif run does not publish pockets. If a separate analysis provides
residue intervals with association statistics, upload a UTF-8 TSV using **Optional residue
annotations**. The template button in the app writes the exact header:

```tsv
comparison_id	protein_id	region_type	region_id	start	end	q_value	prevalence_difference	evidence_source
e3_vs_controls	P12345	POCKET	pocket_1	37	42	0.004	0.22	pocket_run_2026_10
```

Use exact published comparison and protein IDs; positions are 1-based inclusive. Supply
your own, appropriately corrected q-values, positive target-minus-background prevalence
differences, and a traceable evidence source. Select **Uploaded residue annotations** to
colour these in a separate layer. Uploaded q-values are **never combined** with the
published signature q-values. Files are limited to 10 MiB and 10,000 rows. The upload is
held in memory and does not alter the checksum-verified result; the displayed annotations
and plots have downloads.

The sidebar **How to use this page** menu and **Glossary & help** page explain evidence
states, structural terms, colour mapping and interpretation. Every visible table offers
an immediate TSV download and a **Prepare formatted Excel workbook** button; prepared
workbooks have frozen headings, filters, types, readable widths and a column dictionary.
Interactive plots export PNG, PDF and self-contained HTML on demand. The
published SHAP images can be downloaded in their original formats.
