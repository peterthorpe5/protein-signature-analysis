"""Generic glossary and contextual help for the result viewer."""

from __future__ import annotations

GLOSSARY = (
    (
        "Features",
        "Feature ID",
        (
            "A stable key for a tested sequence word, domain, architecture or structural "
            "group. Keep it for reproducibility; read the feature explanation for its "
            "biological meaning."
        ),
    ),
    (
        "Features",
        "k-mer (k3:LPD)",
        (
            "The exact amino-acid word LPD (leucine–proline–aspartate), found anywhere "
            "in a protein. k3 means three residues; k5 means five. Enrichment of a "
            "short word does not alone establish a functional motif."
        ),
    ),
    (
        "Features",
        "Pfam architecture",
        (
            "The ordered IDs of detected Pfam domains along a protein. It is an "
            "architecture, not a single contiguous amino-acid motif."
        ),
    ),
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
            "An SC_ ID is a stable fingerprint of its definition; it is not a named fold "
            "and has no residue boundaries."
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
        "Start with campaign scope, published counts and completed signatures. Exact "
        "feature coverage across every source row is available on request and can take "
        "several minutes for a large campaign. Select another page for detailed tables "
        "and downloadable plots."
    ),
    "Signature explorer": (
        "Choose a comparison, feature family and evidence tier. Start with the short "
        "list of candidates, then expand the complete ledger and inspect counts. "
        "Discovery and held-out validation answer different questions."
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
        "packaged models; enter an exact ID to inspect any other protein. The ranked "
        "association table can show several significant classes together when requested, "
        "independently "
        "of the comparison used to colour the selected model. All comparison outcomes "
        "are downloadable, including insufficient-sample statuses. Colour recorded "
        "intervals and exact enriched k-mer occurrences by target "
        "enrichment (blue to red); white means no significant mapped enrichment. Choose "
        "a structural pair for an exploratory sequence alignment. These optional panels "
        "load only when selected. Pocket intervals may "
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


PAGE_METHODS = {
    "Overview": (
        "**Method.** Published counts describe the immutable campaign. The comparison "
        "census counts positive features at discovery q ≤ 0.05 and then separately "
        "counts features that passed held-out within-comparison and study-wide checks. "
        "No-signature and insufficient-sample outcomes remain separate. The optional "
        "full feature-coverage scan reads the large membership ledger.\n\n"
        "**Limit.** A count of short sequence words is not a count of independent "
        "structural motifs. Automated target labels need human review."
    ),
    "Signature explorer": (
        "**Method.** Each published feature is compared between target and background "
        "protein or independent-unit cohorts. Positive target-minus-background "
        "prevalence and adjusted q ≤ 0.05 propose a discovery hit. The default "
        "decision-candidate view also requires discovery study-wide correction. Held-out "
        "validation is a separate tier; study-wide adjustment covers configured "
        "comparisons in an evidence family. The shortlist takes up to five entries "
        "per feature family, sorted by the selected q-value then effect size.\n\n"
        "**Limit.** Correlated k-mers can produce many similar rows. A small q-value "
        "alone does not imply a large effect or biochemical mechanism. Review controls, "
        "counts, evidence class and validation."
    ),
    "Explainable prediction": (
        "**Method.** A model is fitted on discovery partitions and evaluated on held-out "
        "proteins. Coefficients, permutation importance and SHAP explain model "
        "predictions at different levels. Missing metrics indicate that a model or "
        "validation estimate was unavailable.\n\n"
        "**Limit.** Predictive performance and SHAP attribution do not establish "
        "causation or validate a residue-level motif."
    ),
    "Protein & Pfam": (
        "**Method.** Published identifiers connect the full-length protein sequence "
        "to labels, coordinate-resolved Pfam hits and positive feature memberships. "
        "Domain assessments record whether a scan was performed even when no hit "
        "was found. Evidence sections load separately.\n\n"
        "**Limit.** A missing positive feature row is not an assessed absence; check "
        "the assessment and derivation state before interpreting it."
    ),
    "Classes & roles": (
        "**Method.** Profile labels encode hierarchy, mechanistic class, system "
        "class and component role separately. Observed counts come from published "
        "memberships; the automated-evidence table shows how provisional labels "
        "were assigned. Published comparison-specific cohorts report matched and "
        "excluded target units separately from the shared background pools.\n\n"
        "**Limit.** Nested labels can overlap, so chart areas and class counts should "
        "not be added together as distinct proteins. Provisional labels are hypotheses."
    ),
    "Structures & folds": (
        "**Method.** Coordinate inventory records source, eligibility and model "
        "confidence. Pairwise structural comparisons provide TM-score and coverage. "
        "Thresholded alignment relationships define whole-model connected clusters; "
        "named folds appear only if supplied by a fold authority.\n\n"
        "**Limit.** A cluster is not a local three-dimensional motif. Aggregate "
        "similarity scores do not give residue pairs or cavity positions."
    ),
    "Model & alignment explorer": (
        "**Method.** Positive published decision-candidate memberships with a selected q-value "
        "tier are placed onto their recorded intervals or exact k-mer occurrences "
        "in the selected sequence. Overlaps take the smallest q-value, with colour "
        "saturating at 10⁻⁸. A full-length sequence and residue-number match is "
        "required for 3D projection. An optional pair view computes a new global "
        "amino-acid alignment; uploaded annotations are handled separately.\n\n"
        "**Limit.** The colour is the significance of an enriched protein-level "
        "feature projected onto its occurrences, not a residue-level p-value, "
        "pLDDT or causal proof. Structural clusters have no motif coordinates."
    ),
    "Orthology & partitions": (
        "**Method.** Published orthogroup memberships give homologous context. "
        "Connected homology and redundancy blocks are assigned together to "
        "discovery or validation to limit leakage.\n\n"
        "**Limit.** An orthogroup is a grouping hypothesis, not proof that all "
        "members share one biochemical role."
    ),
    "Canonical data & downloads": (
        "**Method.** The app reads checksum-verified published tables without "
        "altering them. Previews are bounded; complete storage and published reports "
        "remain available. TSV is immediate, while formatted Excel is prepared "
        "when requested.\n\n"
        "**Limit.** A bounded preview or filtered export is not the complete "
        "feature-membership table."
    ),
    "Data quality & provenance": (
        "**Method.** Published evidence audits distinguish provisional labels, "
        "comparison-specific cohorts, outcome-blind pooled matched controls, "
        "excluded label-defining features, "
        "abstentions, assessment states and inference blocks. The full block scan "
        "is run on request.\n\n"
        "**Limit.** Low matched-control coverage or automated target labels can "
        "weaken biological interpretation even if a feature has a small q-value."
    ),
    "Glossary & help": (
        "**Method.** This searchable dictionary explains the app's terms and exports "
        "as TSV or formatted Excel. Each analysis page also describes its methods "
        "and limitations in context.\n\n"
        "**Limit.** Definitions provide orientation; use the campaign's full "
        "configuration and published provenance for reproducible methods."
    ),
}
