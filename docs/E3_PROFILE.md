# Bundled E3 profile

## Purpose and scope

The built-in `e3` profile is a controlled vocabulary and default comparison policy for core
plant and human E3 systems. It is not a universal catalogue of every bacterial, viral or
lineage-specific E3 mechanism, and it is not an automatic E3 annotator. Protein-to-label
assignments remain an evidence-bearing user authority.

The profile has 118 labels and version `1.1.0`. `comparisons: profile_defaults` creates 73
explicit analyses of protein-level mechanistic, family or component-role strata. Parent and
child hypotheses overlap biologically, so they are separate tests rather than statistically
independent hypotheses. Heterogeneous multiprotein-system headings and unresolved catch-all
groups remain browsable but have `default_analysis: false`.

Structural evidence is mandatory for the bundled E3 profile. An E3 campaign must yield an
assessed fold hit, completed structural comparison, feature from a verified structural
resource, or a completed structural search (including an assessed zero-hit search). An eligible
coordinate model is readiness/QC evidence but does not by itself satisfy this result policy.
Read-only validation may report a pending structural search only after Foldseek pre-flight has
succeeded; an AlphaFold request without that enabled search is not structural result evidence.
The runtime repeats the guard after acquisition/search. Custom profiles default to
`require_structural_evidence: false` unless their author opts in.

## Orthogonal fields

Each label can carry:

| Field | Meaning |
|---|---|
| `system_class` | Modifier/system context, such as ubiquitin E3 or control |
| `mechanistic_class` | Transfer mechanism or multiprotein system, such as RING, HECT, RBR or CRL |
| `component_role` | Catalytic E3, RING component, scaffold, adaptor, receptor or other role |
| `family` | Supported family-level name |
| `active_site_expected` | `YES`, `NO` or `UNKNOWN` for an intrinsic catalytic site |
| `active_site_residue` | Named expected catalytic/recruitment feature when defined |
| `default_analysis` | Whether the label is a specific automatic target stratum |
| `default_background_label_id` | Optional background inherited by descendants; the nearest declaration wins before the profile fallback |
| `assignment_exclusivity_group` | Axis on which conflicting reviewed-positive labels are rejected |
| `reviewed_positive_allowed` | Whether direct `REVIEWED_POSITIVE` assignment is meaningful |

These fields must not be collapsed into one label. In particular, an F-box protein is a
CRL1/SCF substrate receptor with an SKP1-binding domain; it is not independently catalytic.
A U-box protein is a RING-like catalytic E3 class. HECT, RBR and RCR proteins use catalytic
cysteines through different mechanisms.

## Default target strata

The 73 default target labels are grouped as follows.

| Group | Default targets |
|---|---|
| RING | RING parent; H2, HC parent, HCa, HCb, C2, D, G, S/T and v subclasses; TRIM-like; MARCH; other reviewed |
| U-box | U-box parent; PUB/ARM; CHIP-like; UFD2/E4-like; other reviewed |
| HECT | HECT parent; plant UPL; NEDD4; HERC; UBE3; HECTD; HUWE; other reviewed |
| RBR/RCR/RNF213 | RBR parent; Ariadne; Parkin; HOIP; HHARI; other RBR; RCR parent; MYCBP2-like; RNF213 RZ-finger |
| Cullin–RING components | RBX and Elongin B-C; CUL1/SKP1/F-box and six F-box architectures; CUL2/VHL-box; CUL3/BTB families; CUL4/DDB1/DCAF; CUL5/SOCS-box; CUL7/FBXW8; CUL9 |
| APC/C components | APC2 scaffold, APC11 RING component and coactivator receptor |
| SUMO mechanisms | SIZ/PIAS, NSE2/MMS21, RanBP2 and ZNF451 |
| NEDD8/UFM1 roles | RBX and DCN1 NEDD8 roles; UFL1, UFBP1 and CDK5RAP3 UFM1-complex roles |

F-box receptor subclasses are Kelch, LRR, WD40, F-box-associated-domain, Tubby and other
reviewed F-box. BTB receptor subclasses include BTB/BACK/Kelch, MATH/BTB and other reviewed
BTB. Mechanistically coherent parent and child analyses are explicit, so their hypotheses are
correlated and the study-wide FDR column should be preferred for broad profile surveying.

CRL1–CRL7 and APC/C system headings are not automatic targets because inherited membership
would pool unrelated scaffolds, adaptors, receptors and catalytic components. The SUMO,
NEDD8 and UFM1 headings are likewise descriptive; their specific mechanism/role children are
tested. Broad atypical-ubiquitin-E3 and other-UBL catch-alls remain available for curation but
are not tested automatically.

Inspect the live vocabulary, label flags and profile metadata rather than copying IDs from a
document:

```bash
protein-signatures describe-profile --profile e3 > e3_profile.json
```

`describe-profile` does not generate campaign comparisons. Those are resolved when a campaign
is validated or run; `protein-signatures validate --config CAMPAIGN.yaml` reports their count.
The canonical profile source is `src/protein_signatures/data/profiles/e3.yaml`.

## Prespecified reference policies

Every bundled default uses a mechanism- or component-role-specific reference. A background
declared on a parent target is inherited by its children; a declaration on the target itself
takes precedence. The profile-level `control:prespecified_non_e3_reference` is an explicit
fallback for custom extensions, but none of the 73 shipped defaults resolves to it.

| Reference stratum | Default comparisons |
|---|---:|
| RING mechanism | 13 |
| U-box mechanism | 5 |
| HECT mechanism | 8 |
| RBR mechanism | 6 |
| RCR mechanism | 2 |
| RNF213 mechanism | 1 |
| RING component | 2 |
| Scaffold | 7 |
| Adaptor | 3 |
| Substrate receptor | 16 |
| CUL9 | 1 |
| SUMO E3 | 4 |
| NEDD8 E3 roles | 2 |
| UFM1 E3 roles | 3 |

These labels are control-set contracts, not automatically inferred negatives. Production
users must assign `REVIEWED_POSITIVE` membership to the exact reference label used by each
comparison, select those proteins before inspecting signatures, and record how taxonomic
composition, length, domain complexity, cellular context and structure availability were
matched. The names do not establish that matching by themselves.

For confirmatory work, replace a default with an explicit campaign comparison whenever the
target needs a narrower control. In particular, a subtype-versus-reference default can find
features common in that subtype, while an explicit sibling comparison is required to ask what
distinguishes it from another subtype. Alternative within-superfamily, housekeeping,
NLR/resistance and shuffled controls remain available but are never substituted automatically.
Shuffled sequences are technical controls and cannot replace a biological reference set.

## Labels excluded from default testing

The roots `e3`, `e3:ubiquitin` and `e3:ubl` are excluded as pooled targets. Direct
`REVIEWED_POSITIVE` assignments are also disallowed on roots and heterogeneous system
headings; proteins are assigned to a specific protein-level child and inherit system
membership. The entire `e3:associated` subtree is excluded and contains:

- reviewed system components that are not independently catalytic;
- putative E3 records;
- pseudoenzymes;
- ambiguous PHD/RING-like proteins; and
- unresolved E3-related records.

These records remain visible for audit but cannot contaminate positive strata because their
labels set `reviewed_positive_allowed: false`.

## Assignment behaviour

Only `REVIEWED_POSITIVE` assignments enter target/background membership. Here “positive”
means positive membership in the named protein type, not necessarily an independently
catalytic protein. For example, a reviewed F-box receptor uses `REVIEWED_POSITIVE` with
`component_role: SUBSTRATE_RECEPTOR`; the role is checked against the profile.

A reviewed positive assigned to a leaf inherits membership in each ancestor. Direct and
inherited memberships are published with `membership_source` and `direct_label_id`. Direct
membership wins if the same protein/label relationship is also inherited through another
reviewed row. A `REVIEWED_POSITIVE` row cannot carry a pending, failed, ambiguous, unknown or
otherwise unreviewed `evidence_status`.

Target and background labels in one comparison must be disjoint. If a protein inherits into
both, the campaign stops rather than guessing which label should win. Conflicting primary
ubiquitin-E3 mechanisms such as RING and U-box are also rejected through the profile
exclusivity group. Proposed and ambiguous rows remain auditable but do not enter training.

## Seed catalogue workflow

Use `prepare-catalogue` to create FASTA, label-review and AlphaFold-accession templates. The
preparer maps no protein automatically: every proposed assignment is `UNMAPPED` until a
reviewer selects a canonical profile label and changes its curation status deliberately. Long
cluster-associated protein-name fields from the supplied 1,000-row catalogue are preserved.

```bash
protein-signatures prepare-catalogue \
  --catalogue e3_seed_catalogue.tsv \
  --output-dir e3_seed_starter \
  --id-column seed_id \
  --sequence-column protein_sequence \
  --name-column associated_seed_protein_names \
  --proposed-category-column associated_seed_categories
```

Catalogue terms such as “associated”, “PHD”, “putative”, “BTB” or “component” must not be
promoted to catalytic positives solely by string matching. Cluster-associated annotations
that mention both RING and U-box require review rather than dual positive assignment.

## Adding or changing subclasses

1. Copy the YAML profile and change `profile_id` and `profile_version`.
2. Add a globally unique label with one known parent.
3. Keep mechanism, system and role fields scientifically distinct.
4. Add unique aliases only when they do not collide case-insensitively.
5. Set `default_analysis: true` only for a coherent protein-level stratum.
6. Set a target or ancestor background when the profile fallback is not defensible.
7. Add exclusivity and reviewed-positive policies where applicable.
8. Add reviewed assignments and tests before analysing it.
9. Treat the profile change as a new campaign identity; never mutate a published result.

For non-E3 work, start from [`configs/profile.example.yaml`](../configs/profile.example.yaml)
and validate it against [`configs/profile.schema.json`](../configs/profile.schema.json).
