"""Generic glossary and contextual help for the result viewer."""

from __future__ import annotations

GLOSSARY = (
    (
        "Study design",
        "Target class",
        ("The protein class or subclass being tested for enrichment."),
    ),
    (
        "Study design",
        "Background class",
        (
            "The comparison proteins used to estimate feature enrichment; check matching "
            "and coverage."
        ),
    ),
    (
        "Study design",
        "Discovery partition",
        ("Proteins used to discover and select features."),
    ),
    (
        "Study design",
        "Validation partition",
        ("Held-out proteins used to test whether selected features reproduce."),
    ),
    (
        "Statistics",
        "Prevalence difference",
        (
            "Target prevalence minus background prevalence; positive values indicate "
            "enrichment in the target."
        ),
    ),
    (
        "Statistics",
        "q-value",
        (
            "A p-value adjusted for multiple tests within a specified family; smaller "
            "values give stronger statistical evidence."
        ),
    ),
    (
        "Statistics",
        "Study-wide q-value",
        ("Adjustment spanning configured comparisons within the same evidence family."),
    ),
    (
        "Statistics",
        "FDR",
        (
            "False discovery rate: the expected fraction of false positives among called "
            "discoveries under the testing procedure."
        ),
    ),
    (
        "Statistics",
        "Evidence class",
        (
            "Summary of discovery, held-out validation and multiplicity support; inspect "
            "the underlying counts."
        ),
    ),
    (
        "Structure",
        "Structural cluster",
        (
            "Whole-model membership based on pairwise structural similarity and coverage. "
            "It has no residue boundaries."
        ),
    ),
    (
        "Structure",
        "TM-score",
        (
            "A normalised structural similarity score; use it with alignment coverage and "
            "the comparison universe."
        ),
    ),
    (
        "Structure",
        "Bilateral coverage",
        ("The aligned fractions of both proteins in a structural comparison."),
    ),
    (
        "Structure",
        "pLDDT",
        (
            "AlphaFold per-residue prediction confidence, not experimental evidence of "
            "function or a statistical enrichment score."
        ),
    ),
    (
        "Structure",
        "Coordinate projection",
        (
            "A recorded 1-based sequence interval mapped to the same positions in an exact "
            "full-length PDB or mmCIF model. Published k-mers can also be projected onto "
            "their exact sequence occurrences."
        ),
    ),
    (
        "Structure",
        "Residue-level annotation",
        (
            "Optional user-supplied interval evidence, including future pocket residues, "
            "with its own q-value and provenance."
        ),
    ),
    (
        "Structure",
        "Pocket",
        (
            "A residue-level cavity or binding-site annotation. The current "
            "structural-cluster output does not supply pocket coordinates."
        ),
    ),
    (
        "Alignment",
        "Sequence alignment",
        (
            "A global alignment of amino-acid strings. It does not locate structurally "
            "superposed residues."
        ),
    ),
    (
        "Alignment",
        "Structural alignment",
        (
            "A comparison of 3D coordinates. This result publishes scores and coverage, not "
            "aligned residue pairs."
        ),
    ),
    (
        "Alignment",
        "Gap",
        ("An inserted alignment column with no residue at that position in one sequence."),
    ),
    (
        "Exports",
        "Completed result",
        (
            "A checksum-verified, immutable bundle containing the database, reports and "
            "coordinate assets."
        ),
    ),
    (
        "Exports",
        "Portable model",
        (
            "A coordinate file inside result/assets/structures; keep it with the complete "
            "result bundle."
        ),
    ),
)

PAGE_HELP = {
    "Overview": (
        "Start with campaign scope, data coverage and the published evidence inventory. "
        "Select another page for detailed tables and downloadable plots."
    ),
    "Signature explorer": (
        "Choose a comparison and feature family. Positive prevalence difference means "
        "target enrichment. Discovery and validation q-values have different roles; "
        "examine both and the evidence class."
    ),
    "Explainable prediction": (
        "Review held-out model metrics and SHAP contributions. Predictive importance is "
        "supporting evidence, not proof that a motif causes a biochemical function."
    ),
    "Protein & Pfam": (
        "Choose a protein to inspect its labels, explicit Pfam sequence intervals, "
        "positive features and assessment states. An unassessed feature is not an "
        "absence."
    ),
    "Classes & roles": (
        "Inspect configured class hierarchies and roles. The same application supports "
        "any published protein profile, not just E3 enzymes."
    ),
    "Structures & folds": (
        "Whole-model folds and clusters are aggregate structural evidence. TM-score, "
        "bilateral coverage and method provenance describe similarity but do not "
        "specify aligned residue positions."
    ),
    "Model & alignment explorer": (
        "Comparisons with significant results appear first. Suggested target proteins have "
        "packaged models; enter an exact ID to inspect any other protein. Colour recorded "
        "intervals and exact enriched k-mer occurrences by target "
        "enrichment (blue to red); white means no significant mapped enrichment. Choose "
        "a structural pair for an exploratory sequence alignment. Pocket intervals may "
        "be supplied separately with a provenance-labelled TSV."
    ),
    "Orthology & partitions": (
        "Check orthogroup membership and which connected sequence or homology blocks "
        "were assigned to discovery or validation."
    ),
    "Canonical data & downloads": (
        "Browse bounded previews and download the complete verified datasets. Excel "
        "workbooks can be prepared on demand and have filters, fixed headings and a "
        "column dictionary."
    ),
    "Data quality & provenance": (
        "Inspect input coverage, exclusions, statuses, evidence provenance and "
        "checksums before drawing scientific conclusions."
    ),
    "Glossary & help": (
        "Search terms by name or definition. Export this dictionary as TSV or formatted "
        "Excel for offline use."
    ),
}
