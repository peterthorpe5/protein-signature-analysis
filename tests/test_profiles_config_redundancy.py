"""Tests for profiles, configuration, catalogue preparation and partitions."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

from protein_signatures.catalogue import (
    _is_uniprot_accession,
    _read_catalogue_rows,
    _wrap_sequence,
    prepare_catalogue,
)
from protein_signatures.config import (
    _parse_alphafold,
    _parse_analysis,
    _parse_comparison,
    _parse_foldseek,
    _parse_inputs,
    _resolve_optional_directory,
    _resolve_optional_path,
    _resolve_profile_source,
    _validate_comparisons,
    config_to_record,
    load_config,
)
from protein_signatures.errors import ConfigurationError, InputValidationError, PublicationError
from protein_signatures.models import (
    ComparisonDefinition,
    CurationStatus,
    GroupMembership,
    LabelAssignment,
    RedundancyClusterMembership,
    SequenceRecord,
)
from protein_signatures.partitions import _partition_for_key, assign_partitions
from protein_signatures.profiles import (
    _validate_hierarchy,
    built_in_profile_path,
    default_profile_comparisons,
    expand_positive_memberships,
    label_ancestors,
    load_profile,
    resolve_label,
    validate_assignment_profile_compatibility,
    validate_comparison_labels,
)
from protein_signatures.redundancy import (
    derive_exact_sequence_clusters,
    read_redundancy_clusters,
)
from protein_signatures.starter import (
    _input_file,
    _optional_input_directory,
    _optional_input_file,
    initialise_campaign,
)


def test_e3_profile_defaults_are_complete_and_conservative() -> None:
    """E3 defaults should cover reviewed strata without promoting ambiguous records."""

    profile = load_profile(source="e3")
    assert built_in_profile_path(profile_name="e3").is_file()
    assert resolve_label(profile=profile, term="F-box") == ("e3:ubiquitin:crl:crl1_scf:f_box")
    comparisons = default_profile_comparisons(profile=profile)
    targets = {item.target_label_ids[0] for item in comparisons}
    assert profile.profile_version == "1.1.0"
    assert profile.require_structural_evidence is True
    assert len(profile.labels) == 118
    assert len(comparisons) == 73
    assert "e3:ubiquitin:ring" in targets
    assert "e3:ubiquitin:u_box" in targets
    assert "e3:ubiquitin:hect" in targets
    assert "e3:ubiquitin:rbr" in targets
    assert "e3:ubiquitin:rnf213_rz" in targets
    assert "e3:ubiquitin:crl:crl1_scf:f_box" in targets
    assert "e3:ubl:sumo:siz_pias" in targets
    assert "e3:ubl:nedd8:dcn1" in targets
    assert "e3:ubl:ufm1:ufl1" in targets
    assert not any(target.startswith("e3:associated") for target in targets)
    assert "e3:ubiquitin:crl" not in targets
    assert "e3:ubiquitin:crl:crl1_scf" not in targets
    assert "e3:ubiquitin:apc_c" not in targets
    assert "e3:ubiquitin:atypical_reviewed" not in targets
    assert "e3:ubl:sumo" not in targets
    assert "e3:ubl:other_reviewed" not in targets
    background_counts = Counter(item.background_label_ids[0] for item in comparisons)
    assert background_counts == {
        "control:matched_adaptor_reference": 3,
        "control:matched_cul9_reference": 1,
        "control:matched_hect_mechanism_reference": 8,
        "control:matched_nedd8_e3_reference": 2,
        "control:matched_rbr_mechanism_reference": 6,
        "control:matched_rcr_mechanism_reference": 2,
        "control:matched_ring_component_reference": 2,
        "control:matched_ring_mechanism_reference": 13,
        "control:matched_rnf213_mechanism_reference": 1,
        "control:matched_scaffold_reference": 7,
        "control:matched_substrate_receptor_reference": 16,
        "control:matched_sumo_e3_reference": 4,
        "control:matched_u_box_mechanism_reference": 5,
        "control:matched_ufm1_e3_reference": 3,
    }
    assert "control:prespecified_non_e3_reference" not in background_counts
    assert all("Prespecified comparison" in item.description for item in comparisons)
    assert all("analysed separately" in item.description for item in comparisons)
    labels = {item.label_id: item for item in profile.labels}
    assert (
        "histidines at metal-ligand positions four and five"
        in labels["e3:ubiquitin:ring:ring_h2"].description
    )
    assert "two residues" in labels["e3:ubiquitin:ring:ring_hc:ring_hca"].description
    assert "three or four residues" in labels["e3:ubiquitin:ring:ring_hc:ring_hcb"].description
    hect = next(item for item in profile.labels if item.label_id == "e3:ubiquitin:hect")
    assert hect.active_site_expected == "YES"
    assert hect.active_site_residue == "CATALYTIC_CYSTEINE"
    assert label_ancestors(profile=profile, label_id=hect.label_id)[-1] == "protein"
    with pytest.raises(InputValidationError, match="Unknown"):
        resolve_label(profile=profile, term="not a class")
    with pytest.raises(InputValidationError, match="Unknown"):
        label_ancestors(profile=profile, label_id="missing")
    with pytest.raises(ConfigurationError, match="Unknown built-in"):
        built_in_profile_path(profile_name="missing")


def test_profile_membership_expansion_uses_reviewed_positive_only() -> None:
    """Inherited memberships should retain the most direct reviewed authority."""

    profile = load_profile(source="e3")
    child = "e3:ubiquitin:crl:crl1_scf:f_box:kelch"
    parent = "e3:ubiquitin:crl:crl1_scf:f_box"
    assignments = (
        _assignment(protein_id="p1", label_id=child, status=CurationStatus.REVIEWED_POSITIVE),
        _assignment(protein_id="p1", label_id=parent, status=CurationStatus.REVIEWED_POSITIVE),
        _assignment(protein_id="p2", label_id=child, status=CurationStatus.PROPOSED),
    )
    rows = expand_positive_memberships(assignments=assignments, profile=profile)
    assert not any(row["protein_id"] == "p2" for row in rows)
    direct_parent = next(
        row for row in rows if row["protein_id"] == "p1" and row["label_id"] == parent
    )
    assert direct_parent["membership_source"] == "DIRECT"
    valid = ComparisonDefinition("ok", "OK", (parent,), ("control:matched_non_e3",), "")
    validate_comparison_labels(comparisons=(valid,), profile=profile)
    with pytest.raises(InputValidationError, match="unknown labels"):
        validate_comparison_labels(
            comparisons=(replace(valid, target_label_ids=("missing",)),),
            profile=profile,
        )


def test_profile_assignment_compatibility_rejects_roles_states_and_mechanism_conflicts() -> None:
    """Reviewed positives must obey profile roles, states and exclusive mechanisms."""

    profile = load_profile(source="e3")
    f_box = "e3:ubiquitin:crl:crl1_scf:f_box"
    ring = "e3:ubiquitin:ring:ring_h2"
    u_box = "e3:ubiquitin:u_box"
    valid = replace(
        _assignment(
            protein_id="p1",
            label_id=f_box,
            status=CurationStatus.REVIEWED_POSITIVE,
        ),
        component_role="SUBSTRATE_RECEPTOR",
    )
    validate_assignment_profile_compatibility(assignments=(valid,), profile=profile)
    proposed_wrong_role = replace(valid, curation_status=CurationStatus.PROPOSED)
    validate_assignment_profile_compatibility(
        assignments=(proposed_wrong_role,),
        profile=profile,
    )
    with pytest.raises(InputValidationError, match="expected 'SUBSTRATE_RECEPTOR'"):
        validate_assignment_profile_compatibility(
            assignments=(replace(valid, component_role="CATALYTIC_E3"),),
            profile=profile,
        )
    associated = replace(
        valid,
        label_id="e3:associated:putative",
        component_role="PUTATIVE_CATALYTIC_E3",
    )
    with pytest.raises(InputValidationError, match="does not permit"):
        validate_assignment_profile_compatibility(assignments=(associated,), profile=profile)
    with pytest.raises(InputValidationError, match="incompatible evidence_status"):
        validate_assignment_profile_compatibility(
            assignments=(replace(valid, evidence_status="UNREVIEWED"),),
            profile=profile,
        )
    with pytest.raises(InputValidationError, match="unknown profile label"):
        validate_assignment_profile_compatibility(
            assignments=(replace(valid, label_id="missing"),),
            profile=profile,
        )
    incompatible = (
        replace(valid, label_id=ring, component_role="CATALYTIC_E3"),
        replace(valid, label_id=u_box, component_role="CATALYTIC_E3"),
    )
    with pytest.raises(InputValidationError, match="exclusivity group"):
        validate_assignment_profile_compatibility(assignments=incompatible, profile=profile)


def test_profile_hierarchy_validation_rejects_each_ambiguous_shape() -> None:
    """Empty, duplicate, multi-root, broken, cyclic and alias-conflicting profiles fail."""

    profile = load_profile(source="e3")
    with pytest.raises(InputValidationError, match="at least one"):
        _validate_hierarchy(profile=replace(profile, labels=()))
    with pytest.raises(InputValidationError, match="Duplicate"):
        _validate_hierarchy(profile=replace(profile, labels=(profile.labels[0], profile.labels[0])))
    second_root = replace(profile.labels[1], parent_label_id="")
    with pytest.raises(InputValidationError, match="exactly one root"):
        _validate_hierarchy(profile=replace(profile, labels=(profile.labels[0], second_root)))
    broken = replace(profile.labels[1], parent_label_id="missing")
    with pytest.raises(InputValidationError, match="unknown parent"):
        _validate_hierarchy(profile=replace(profile, labels=(profile.labels[0], broken)))
    alias = replace(profile.labels[1], aliases=(profile.labels[0].label_id,))
    with pytest.raises(InputValidationError, match="shared"):
        _validate_hierarchy(profile=replace(profile, labels=(profile.labels[0], alias)))
    first = replace(profile.labels[1], parent_label_id=profile.labels[2].label_id)
    second = replace(profile.labels[2], parent_label_id=profile.labels[1].label_id)
    with pytest.raises(InputValidationError, match="cycle"):
        _validate_hierarchy(profile=replace(profile, labels=(profile.labels[0], first, second)))
    with pytest.raises(InputValidationError, match="policy"):
        _validate_hierarchy(profile=replace(profile, default_background_label_id="missing"))
    self_background = replace(
        profile.labels[1],
        default_background_label_id=profile.labels[1].label_id,
    )
    with pytest.raises(InputValidationError, match="itself"):
        _validate_hierarchy(profile=replace(profile, labels=(profile.labels[0], self_background)))
    disallowed_default = replace(
        profile.labels[1],
        default_analysis=True,
        reviewed_positive_allowed=False,
    )
    with pytest.raises(InputValidationError, match="must permit"):
        _validate_hierarchy(
            profile=replace(profile, labels=(profile.labels[0], disallowed_default))
        )


def test_profile_loading_and_default_policy_errors(tmp_path: Path) -> None:
    """Missing YAML, malformed YAML and incomplete default policies fail clearly."""

    with pytest.raises(InputValidationError, match="Missing or empty"):
        load_profile(source=tmp_path / "missing.yaml")
    malformed = tmp_path / "malformed.yaml"
    malformed.write_text("labels: [", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Invalid profile YAML"):
        load_profile(source=malformed)
    profile_typo = tmp_path / "profile_typo.yaml"
    profile_document = yaml.safe_load(
        built_in_profile_path(profile_name="e3").read_text(encoding="utf-8")
    )
    profile_document["unexpected"] = True
    profile_typo.write_text(yaml.safe_dump(profile_document), encoding="utf-8")
    with pytest.raises(InputValidationError, match="unknown fields"):
        load_profile(source=profile_typo)
    profile = load_profile(source="e3")
    with pytest.raises(InputValidationError, match="profile_defaults requires"):
        default_profile_comparisons(
            profile=replace(
                profile,
                default_target_root_label_id="",
                default_background_label_id="",
            )
        )
    with pytest.raises(InputValidationError, match="unknown labels"):
        default_profile_comparisons(profile=replace(profile, default_background_label_id="missing"))


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("require_structural_evidence", "yes", "require_structural_evidence"),
        ("default_analysis", "yes", "default_analysis"),
        ("reviewed_positive_allowed", "no", "reviewed_positive_allowed"),
    ),
)
def test_profile_boolean_policies_are_strict(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    """Profile policy booleans must be YAML booleans rather than truthy text."""

    document = yaml.safe_load(built_in_profile_path(profile_name="e3").read_text(encoding="utf-8"))
    if field == "require_structural_evidence":
        document[field] = value
    else:
        document["labels"][0][field] = value
    path = tmp_path / f"invalid_{field}.yaml"
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    with pytest.raises(InputValidationError, match=message):
        load_profile(source=path)


def test_configuration_defaults_and_serialisation(example_dir: Path, tmp_path: Path) -> None:
    """Profile defaults and path resolution should serialise deterministically."""

    document = _example_document(example_dir=example_dir)
    document["comparisons"] = "profile_defaults"
    path = _write_config(path=tmp_path / "campaign.yaml", document=document)
    config = load_config(path=path)
    assert config.comparisons == ()
    record = config_to_record(config=config)
    assert record["analysis"]["kmer_lengths"] == [3, 4]
    assert record["campaign"]["profile"] == "e3"
    assert Path(record["inputs"]["sequences_fasta"]).is_absolute()
    assert _resolve_profile_source(value="e3", base=tmp_path) == "e3"
    custom = _resolve_profile_source(value="custom.yaml", base=tmp_path)
    assert custom == str((tmp_path / "custom.yaml").resolve())


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda row: row.update(schema_version=2), "Unsupported schema_version"),
        (lambda row: row.update(comparisons=[]), "non-empty"),
        (lambda row: row["analysis"].update(kmer_lengths=[]), "kmer_lengths"),
        (lambda row: row["analysis"].update(kmer_lengths=[3, 3]), "duplicates"),
        (lambda row: row["analysis"].update(kmer_lengths=[13]), "kmer_lengths"),
        (lambda row: row["alphafold"].update(enabled="yes"), "enabled"),
        (lambda row: row["foldseek"].update(enabled="yes"), "enabled"),
        (lambda row: row["foldseek"].update(executable="relative/path"), "executable"),
        (lambda row: row.update(analaysis={}), "unknown fields"),
    ],
)
def test_configuration_rejects_invalid_documents(
    example_dir: Path,
    tmp_path: Path,
    mutation: object,
    message: str,
) -> None:
    """Invalid schema, settings and executable values should fail closed."""

    document = _example_document(example_dir=example_dir)
    mutation(document)  # type: ignore[operator]
    path = _write_config(path=tmp_path / "invalid.yaml", document=document)
    with pytest.raises((ConfigurationError, InputValidationError), match=message):
        load_config(path=path)


def test_configuration_parsers_cover_bounds_and_orthofinder_modes(
    example_dir: Path, tmp_path: Path
) -> None:
    """Private parsers should preserve valid alternatives and reject contradictions."""

    inputs = _example_document(example_dir=example_dir)["inputs"]
    inputs["orthofinder"] = {
        "resource_dir": str(tmp_path),
        "results_dir": str(tmp_path),
    }
    with pytest.raises(ConfigurationError, match="mutually exclusive"):
        _parse_inputs(value=inputs, base=tmp_path)
    inputs["orthofinder"] = {"group_type": "other"}
    with pytest.raises(ConfigurationError, match="group_type"):
        _parse_inputs(value=inputs, base=tmp_path)
    inputs["orthofinder"] = {
        "group_type": "LEGACY_ORTHOGROUP",
        "hierarchy_node": "N0",
    }
    with pytest.raises(ConfigurationError, match="must be empty"):
        _parse_inputs(value=inputs, base=tmp_path)
    inputs["orthofinder"] = {
        "group_type": "LEGACY_ORTHOGROUP",
        "hierarchy_node": "",
        "run_id": "run_1",
    }
    parsed = _parse_inputs(value=inputs, base=tmp_path)
    assert parsed.orthofinder_run_id == "run_1"
    assert parsed.orthofinder_hierarchy_node == ""
    absolute_cache = tmp_path / "cache"
    assert _parse_alphafold(value={"cache_dir": str(absolute_cache)}, base=tmp_path).cache_dir == (
        absolute_cache
    )
    assert _parse_foldseek(value={"cache_dir": str(absolute_cache)}, base=tmp_path).cache_dir == (
        absolute_cache
    )
    for value in ({"kmer_lengths": [0]}, {"validation_fraction": 1.0}):
        with pytest.raises((ConfigurationError, InputValidationError)):
            _parse_analysis(value=value)


def test_comparison_and_optional_path_validation(tmp_path: Path) -> None:
    """Comparison set logic and optional path handling should reject ambiguous input."""

    valid = _parse_comparison(
        value={
            "comparison_id": "a",
            "target_label_ids": ["target"],
            "background_label_ids": ["background"],
        },
        index=0,
    )
    assert valid.display_name == "a"
    for value in (
        {"comparison_id": "a", "target_label_ids": [], "background_label_ids": ["b"]},
        {
            "comparison_id": "a",
            "target_label_ids": ["same"],
            "background_label_ids": ["same"],
        },
    ):
        with pytest.raises(ConfigurationError):
            _parse_comparison(value=value, index=0)
    with pytest.raises(ConfigurationError, match="unique"):
        _validate_comparisons(comparisons=(valid, valid))
    assert _resolve_optional_path(value=None, base=tmp_path) is None
    assert _resolve_optional_directory(value="", base=tmp_path) is None
    with pytest.raises(InputValidationError, match="missing or empty"):
        _resolve_optional_path(value="missing.tsv", base=tmp_path)
    with pytest.raises(InputValidationError, match="does not exist"):
        _resolve_optional_directory(value="missing", base=tmp_path)


def test_initialise_campaign_writes_a_valid_generic_configuration(
    example_dir: Path, tmp_path: Path
) -> None:
    """The starter command should make a usable config without hand-written YAML."""

    config_path = tmp_path / "campaign.yaml"
    created = initialise_campaign(
        config_path=config_path,
        campaign_id="starter_1",
        profile="e3",
        sequences_fasta=example_dir / "proteins.faa",
        label_assignments=example_dir / "label_assignments.tsv",
        domains=example_dir / "domains.tsv",
        structures=example_dir / "structures.tsv",
        structure_comparisons=example_dir / "structure_comparisons.tsv",
    )
    assert created == config_path
    assert load_config(path=created).campaign_id == "starter_1"
    with pytest.raises(PublicationError, match="already exists"):
        initialise_campaign(
            config_path=config_path,
            campaign_id="starter_1",
            profile="e3",
            sequences_fasta=example_dir / "proteins.faa",
            label_assignments=example_dir / "label_assignments.tsv",
        )


def test_initialise_campaign_rejects_incompatible_sources(
    example_dir: Path, tmp_path: Path
) -> None:
    """Starter validation should reject ambiguous sources and unsafe structural options."""

    common = {
        "config_path": tmp_path / "campaign.yaml",
        "campaign_id": "starter",
        "profile": "e3",
        "sequences_fasta": example_dir / "proteins.faa",
        "label_assignments": example_dir / "label_assignments.tsv",
    }
    with pytest.raises(InputValidationError, match="not both"):
        initialise_campaign(
            **common,
            orthofinder_resource=tmp_path,
            orthofinder_results=tmp_path,
        )
    with pytest.raises(InputValidationError, match="requires"):
        initialise_campaign(**common, enable_alphafold=True)
    with pytest.raises(InputValidationError, match="group_type"):
        initialise_campaign(**common, orthofinder_group_type="OTHER")
    with pytest.raises(InputValidationError, match="must be blank"):
        initialise_campaign(
            **common,
            orthofinder_group_type="LEGACY_ORTHOGROUP",
            orthofinder_hierarchy_node="N0",
        )
    assert _optional_input_file(path=None, field="x") is None
    assert _optional_input_directory(path=None, field="x") is None
    with pytest.raises(InputValidationError, match="missing or empty"):
        _input_file(path=tmp_path / "missing", field="x")
    with pytest.raises(InputValidationError, match="does not exist"):
        _optional_input_directory(path=tmp_path / "missing", field="x")


def test_catalogue_preparation_never_promotes_unreviewed_categories(tmp_path: Path) -> None:
    """A sequence catalogue should yield review templates with zero training positives."""

    source = tmp_path / "catalogue.tsv"
    long_name = "Long associated annotation " * 120
    source.write_text(
        "identifier\tsequence\tname\tcategory\n"
        f"Q9SA03\tACDEFG\t{long_name}\tF-box\n"
        "local_2\tACDFFG\tProtein two\tRing finger\n",
        encoding="utf-8",
    )
    output = prepare_catalogue(
        catalogue_path=source,
        output_dir=tmp_path / "prepared",
        id_column="identifier",
        sequence_column="sequence",
        name_column="name",
        proposed_category_column="category",
    )
    marker = json.loads((output / "PREPARED.json").read_text(encoding="utf-8"))
    assert marker["protein_count"] == 2
    assert marker["training_eligible_assignment_count"] == 0
    assert long_name.strip() in (output / "proteins.faa").read_text(encoding="utf-8")
    assignments = (output / "label_assignments.REVIEW_REQUIRED.tsv").read_text(encoding="utf-8")
    assert assignments.count("\tUNMAPPED\t") == 2
    accessions = (output / "alphafold_accessions.REVIEW_REQUIRED.tsv").read_text(encoding="utf-8")
    assert "Q9SA03" in accessions and "local_2" not in accessions
    with pytest.raises(PublicationError, match="already exists"):
        prepare_catalogue(
            catalogue_path=source,
            output_dir=output,
            id_column="identifier",
            sequence_column="sequence",
        )
    assert _wrap_sequence(sequence="A" * 81).splitlines() == ["A" * 80, "A"]
    with pytest.raises(InputValidationError):
        _wrap_sequence(sequence="", width=80)
    assert _is_uniprot_accession(value="A0A060D0U3") is True
    assert _is_uniprot_accession(value="local_2") is False


def test_catalogue_parser_rejects_missing_and_duplicate_values(tmp_path: Path) -> None:
    """Catalogue parsing should reject blank sequences and all duplicate identities."""

    with pytest.raises(InputValidationError, match="distinct"):
        prepare_catalogue(
            catalogue_path=tmp_path / "missing",
            output_dir=tmp_path / "out",
            id_column="id",
            sequence_column="id",
        )
    with pytest.raises(InputValidationError, match="required"):
        prepare_catalogue(
            catalogue_path=tmp_path / "missing",
            output_dir=tmp_path / "out",
            id_column="",
            sequence_column="sequence",
        )
    for name, rows, message in (
        ("blank", "p1\t\n", "no sequence"),
        ("duplicate", "p1\tACD\np1\tACD\n", "duplicate identifier"),
        ("conflict", "p1\tACD\np1\tACE\n", "conflicting sequences"),
    ):
        path = tmp_path / f"{name}.tsv"
        path.write_text("id\tsequence\n" + rows, encoding="utf-8")
        with pytest.raises(InputValidationError, match=message):
            _read_catalogue_rows(
                path=path,
                id_column="id",
                sequence_column="sequence",
                name_column="",
                proposed_category_column="",
            )
    invalid = tmp_path / "invalid.tsv"
    invalid.write_text("id\tsequence\np1\tAC1\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="unsupported residue"):
        prepare_catalogue(
            catalogue_path=invalid,
            output_dir=tmp_path / "invalid_out",
            id_column="id",
            sequence_column="sequence",
        )


def test_redundancy_authorities_and_connected_partitions(tmp_path: Path) -> None:
    """Exact, HOG and near-redundancy links should form indivisible partition blocks."""

    sequences = _sequences()
    exact = derive_exact_sequence_clusters(sequences=sequences)
    assert len(exact) == 4
    assert exact[0].cluster_type == "EXACT_SEQUENCE"
    table = tmp_path / "redundancy.tsv"
    table.write_text(
        "protein_id\tcluster_id\tmethod\tmethod_version\tidentity_threshold\t"
        "coverage_threshold\tevidence_reference\n"
        "p2\tNR1\tMMSEQS2\t15.6\t0.9\t0.8\tRUN:1\n"
        "p3\tNR1\tMMSEQS2\t15.6\t0.9\t0.8\tRUN:1\n",
        encoding="utf-8",
    )
    near = read_redundancy_clusters(path=table, protein_ids=frozenset({"p1", "p2", "p3", "p4"}))
    memberships = (
        _membership(protein_id="p1", group_id="HOG1"),
        _membership(protein_id="p2", group_id="HOG2"),
    )
    partitions = assign_partitions(
        sequences=sequences,
        memberships=memberships,
        redundancy_memberships=near,
        validation_fraction=0.5,
        random_seed=7,
    )
    by_id = {item.protein_id: item for item in partitions}
    assert by_id["p1"].partition_key == by_id["p2"].partition_key
    assert by_id["p2"].partition_key == by_id["p3"].partition_key
    assert by_id["p1"].partition == by_id["p3"].partition
    assert by_id["p1"].partition_unit == "REDUNDANCY_BLOCK"
    assert by_id["p4"].partition_unit == "EXACT_SEQUENCE"
    assert all(
        item.partition == "DISCOVERY"
        for item in assign_partitions(
            sequences=sequences,
            memberships=(),
            validation_fraction=0.0,
            random_seed=0,
        )
    )
    assert _partition_for_key(key="x", validation_fraction=1.0, random_seed=0) == "VALIDATION"


def test_redundancy_and_partition_validation_errors(tmp_path: Path) -> None:
    """Malformed cluster rows and contradictory memberships should fail closed."""

    sequences = _sequences()
    protein_ids = frozenset(item.protein_id for item in sequences)
    header = (
        "protein_id\tcluster_id\tmethod\tmethod_version\tidentity_threshold\t"
        "coverage_threshold\tevidence_reference\n"
    )
    for name, rows, message in (
        ("unknown", "missing\tC1\tM\t1\t0.9\t0.8\tR\n", "absent"),
        ("blank", "p1\tC1\tM\t1\t\t0.8\tR\n", "must be populated"),
        (
            "duplicate",
            "p1\tC1\tM\t1\t0.9\t0.8\tR\np1\tC1\tM\t1\t0.9\t0.8\tR\n",
            "more than once",
        ),
    ):
        path = tmp_path / f"{name}.tsv"
        path.write_text(header + rows, encoding="utf-8")
        with pytest.raises(InputValidationError, match=message):
            read_redundancy_clusters(path=path, protein_ids=protein_ids)
    with pytest.raises(InputValidationError, match="validation_fraction"):
        assign_partitions(
            sequences=sequences,
            memberships=(),
            validation_fraction=0.95,
            random_seed=1,
        )
    with pytest.raises(InputValidationError, match="random_seed"):
        assign_partitions(
            sequences=sequences,
            memberships=(),
            validation_fraction=0.2,
            random_seed=-1,
        )
    with pytest.raises(InputValidationError, match="absent"):
        assign_partitions(
            sequences=sequences,
            memberships=(_membership(protein_id="missing", group_id="HOG1"),),
            validation_fraction=0.2,
            random_seed=1,
        )
    with pytest.raises(InputValidationError, match="multiple partition"):
        assign_partitions(
            sequences=sequences,
            memberships=(
                _membership(protein_id="p1", group_id="HOG1"),
                _membership(protein_id="p1", group_id="HOG2"),
            ),
            validation_fraction=0.2,
            random_seed=1,
        )
    unknown = _near(protein_id="missing", cluster_id="C1")
    with pytest.raises(InputValidationError, match="absent"):
        assign_partitions(
            sequences=sequences,
            memberships=(),
            redundancy_memberships=(unknown,),
            validation_fraction=0.2,
            random_seed=1,
        )
    with pytest.raises(InputValidationError, match="multiple redundancy"):
        assign_partitions(
            sequences=sequences,
            memberships=(),
            redundancy_memberships=(
                _near(protein_id="p1", cluster_id="C1"),
                _near(protein_id="p1", cluster_id="C2"),
            ),
            validation_fraction=0.2,
            random_seed=1,
        )


def _example_document(*, example_dir: Path) -> dict[str, object]:
    """Load the example configuration with absolute input paths."""

    document = yaml.safe_load((example_dir / "campaign.yaml").read_text(encoding="utf-8"))
    for field in (
        "sequences_fasta",
        "label_assignments",
        "domains",
        "structures",
        "structure_comparisons",
    ):
        document["inputs"][field] = str((example_dir / document["inputs"][field]).resolve())
    return document


def _write_config(*, path: Path, document: dict[str, object]) -> Path:
    """Write one test YAML document."""

    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def _assignment(*, protein_id: str, label_id: str, status: CurationStatus) -> LabelAssignment:
    """Build one compact curation assignment."""

    return LabelAssignment(protein_id, label_id, status, "REVIEWED", "test", "", "ROLE", "test")


def _sequences() -> tuple[SequenceRecord, ...]:
    """Return records containing one exact duplicate pair."""

    return (
        SequenceRecord("p1", "", "AAAA", 4, "a" * 64),
        SequenceRecord("p2", "", "AAAA", 4, "a" * 64),
        SequenceRecord("p3", "", "CCCC", 4, "c" * 64),
        SequenceRecord("p4", "", "DDDD", 4, "d" * 64),
    )


def _membership(*, protein_id: str, group_id: str) -> GroupMembership:
    """Build one HOG membership."""

    return GroupMembership("run", "HOG", "N0", group_id, "OG", "N0", "species", protein_id)


def _near(*, protein_id: str, cluster_id: str) -> RedundancyClusterMembership:
    """Build one near-redundancy membership."""

    return RedundancyClusterMembership(
        protein_id,
        cluster_id,
        "NEAR_REDUNDANCY",
        "TEST",
        "1",
        0.9,
        0.8,
        "test",
    )
