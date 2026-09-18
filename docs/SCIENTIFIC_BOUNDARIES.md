# Scientific boundaries and claim language

## What a result can support

The package can support statements such as:

- “PF00646 was more prevalent in reviewed F-box-receptor blocks than in the specified
  controls under this campaign.”
- “A discovery-defined Foldseek neighbourhood was common in this reviewed class and showed
  the same enrichment direction in held-out homology blocks.”
- “This fitted classifier relied strongly on the following domain, sequence and structural
  features for this prediction.”

It cannot by itself support statements such as:

- “this protein is biochemically an E3”;
- “this residue is catalytic”;
- “this pocket is ligandable or druggable”;
- “absence of a hit proves absence of a domain or biological function”; or
- “the model explanation reveals the causal mechanism.”

Experimental validation and expert review remain necessary.

## Labels are an authority, not an inferred outcome

Only `REVIEWED_POSITIVE` assignments define class membership. Proposed, ambiguous,
associated, pseudoenzyme and non-catalytic-component records remain in their explicit
states. A new catalogue or publication can propose labels, but it must not silently replace
the reviewed seed authority.

The same evidence should not be used circularly both to assign the target label and then be
presented as independent confirmation of that label. When a Pfam domain is part of the
curation rule, its enrichment remains descriptive and should be labelled as expected rather
than independent validation.

## E3 is heterogeneous

There is no supported assumption that every E3-related protein shares one universal short
motif or fold. RING/U-box proteins transfer ubiquitin without an obligatory E3–ubiquitin
thioester; HECT and RBR/RCR mechanisms use catalytic cysteines; Cullin–RING systems separate
scaffold, adaptor, receptor and RING catalytic roles. F-box is a CRL1 substrate-receptor
component, whereas U-box is a RING-like catalytic mechanism. Pooling these without an
explicit exploratory comparison would erase the relevant biology.

A PHD or generic zinc-binding/RING-like annotation is not sufficient proof of E3 activity.
Ambiguity must remain visible.

## Domain evidence

`ASSESSED_NO_HIT` means a declared authority successfully assessed the protein and found no
qualifying hit. `NOT_ASSESSED` means no such conclusion can be drawn. `FAILED` records an
attempt that did not produce valid evidence. These states are never interchangeable.

Pfam and architecture associations can be excellent family signatures while still being
neither necessary nor sufficient for function. Coordinate errors, fragments and database
release changes can alter calls.

## Structure evidence

AlphaFold DB supplies predicted models. Model availability, global mean pLDDT and structural
similarity do not establish experimental conformation, interaction state, pocket geometry or
activity. Low pLDDT is not automatically disorder, and a high global score does not validate
every local residue. PAE and residue-level confidence should accompany local geometric
claims when supplied by an upstream authority. The workflow retains excluded downloads for
provenance but admits only exact sequence-verified models meeting the configured global
confidence threshold to local coordinate analysis.

TM score and coverage threshold a pairwise relationship; a connected component can contain
members that do not all pass pairwise thresholds with one another. It is a structural
neighbourhood, not a proof that all members share a formally assigned fold. Reviewed
SCOP/CATH/ECOD-style fold calls stay separate.

Coverage fractions from full sequences, coordinate-model residues and isolated constructs are
not interchangeable. Their scope is explicit, and graph components cannot bridge comparison
universes, scopes, tools or versions. An absent retained hit is `UNKNOWN` unless a complete
search-universe membership proves that the protein was assessed.

The predecessor E3 structural workflow studies within-group alignments and conserved pocket
position. The generic Foldseek layer screens across campaign structures. Their outputs may
coexist, but neither is silently substituted for the other. Predecessor pocket feature IDs are
cluster-specific; pocket position is not assumed comparable across unrelated alignments.

## Controls, phylogeny and generalisation

Association quality depends on control construction. Controls should be selected before
observing the candidate signature and, where possible, matched for taxonomic distribution,
length, domain complexity, subcellular context, sequence completeness and structure
availability. Shuffled sequences are useful for composition artefacts but are not adequate
as the only biological control.

Connected HOG/redundancy blocks reduce obvious pseudoreplication; they do not implement a
full phylogenetic comparative model. Very deep shared ancestry, sampling imbalance and
species-specific annotation practice can still confound an association.

The validation split is held out from discovery/model fitting, but it is not permanently
independent after a scientist uses it to redesign the workflow. A subsequent confirmatory
study needs a new untouched set.

## Multiple testing

Local and study-wide Benjamini–Hochberg q-values control defined test families, not the
probability that an individual biological claim is true. The profile contains correlated
parent, child and component comparisons; the study-wide correction spans them but is not an
ontology-aware hierarchical FDR model. Report the exact comparison, feature family, effect,
sample/block counts and both q-values.

## Explainable AI

The classifier is intentionally interpretable and group-aware, but a high ROC AUC can still
reflect dataset construction, hidden availability differences or taxonomic shortcuts.
Technical availability features are excluded by default, and label-blind feature selection
reduces one source of optimistic bias. These safeguards do not make the model causal.

SHAP values are local attributions in the fitted model's log-odds space. They are not binding
energies, residue contacts, evolutionary constraints or experimental effect sizes. Strongly
correlated features can share attribution unpredictably. Read SHAP beside the association
ledger and held-out metrics, never as a replacement for them.

## Native and imported methods

Version 0.2.0 natively implements exact k-mer presence, supplied domains/architectures,
reviewed fold calls, AlphaFold DB acquisition, Foldseek structural comparisons,
discovery-frozen structural clusters, block-level association and mandatory explainable
modelling. Generic `features.tsv` can carry results from MEME/STREME/FIMO, HMMER,
conservation, disorder or experimental curation.

Generic feature definitions must be independently prespecified or learned label-blind from
the frozen discovery sequence cohort. Definition and cohort SHA-256 values make definition
identity and declared cohort membership machine-checkable, but they cannot prove whether an
external producer inspected class labels. `DISCOVERY_DERIVED` is therefore an input
attestation of label-blind derivation. Label-aware discovery definitions must declare
`DISCOVERY_SUPERVISED`; they remain audit-only until validation-only inference is
implemented. All-data-derived generic definitions are explicitly exploratory and are also
excluded from confirmatory association and modelling. Per-feature assessed universes keep a
missing or failed call unknown rather than treating it as biological absence.

Import support does not mean the package ran or validated the biological assumptions of the
producer. Every imported row must name its source/version/reference. The package must not be
cited as the native producer of those calls.
