"""Focused orchestration tests for optional campaign branches and helpers."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import duckdb
import pytest
import yaml

import protein_signatures.pipeline as pipeline_module
from protein_signatures.checksums import sha256_file
from protein_signatures.config import load_config
from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.models import (
    AlphaFoldAcquisition,
    AlphaFoldRequest,
    DomainAssessment,
    DomainAssessmentStatus,
    FeatureRecord,
    FoldEvidenceStatus,
    FoldseekRunEvidence,
    FoldseekSettings,
    PairwiseStructureComparison,
    StructuralAlignmentImport,
    StructureAnalysisEligibility,
    StructureComparisonStatus,
    StructureCoverageScope,
    StructureRecord,
)
from protein_signatures.pipeline import (
    _alphafold_record,
    _input_authorities,
    _merge_features,
    _merge_structure_comparisons,
    _merge_structures,
    _pfam_availability,
    _prepare_campaign,
    _structure_asset_sources,
    _structure_record,
    _validate_foldseek_preflight,
    _validate_requested_resource_run,
    _validate_required_structural_evidence,
    _with_effective_comparisons,
    run_campaign,
    validate_campaign,
)
from protein_signatures.profiles import load_profile


def test_validate_campaign_handles_alphafold_selection_and_missing_mapping(
    example_dir: Path,
    tmp_path: Path,
) -> None:
    """Read-only validation should check selected AlphaFold mappings without downloading."""

    document = _campaign_document(example_dir=example_dir)
    document["alphafold"]["enabled"] = True
    missing_path = _write_document(tmp_path / "missing.yaml", document)
    with pytest.raises(InputValidationError, match="alphafold.enabled"):
        validate_campaign(config_path=missing_path)
    accessions = tmp_path / "alphafold.tsv"
    accessions.write_text("protein_id\tuniprot_accession\nfbox01\tQ9SA03\n", encoding="utf-8")
    document["inputs"]["alphafold_accessions"] = str(accessions)
    valid_path = _write_document(tmp_path / "valid.yaml", document)
    assert validate_campaign(config_path=valid_path)["alphafold_enabled"] is True


def test_e3_profile_requires_eligible_structural_evidence(
    example_dir: Path,
    tmp_path: Path,
) -> None:
    """E3 validation should fail without evidence but allow pending AlphaFold work."""

    profile = load_profile(source="e3")
    with pytest.raises(InputValidationError, match="requires structural evidence"):
        _validate_required_structural_evidence(
            profile=profile,
            structures=(),
            comparisons=(),
        )
    pending = _validate_required_structural_evidence(
        profile=profile,
        structures=(),
        comparisons=(),
        pending_alphafold_request_count=1,
        pending_structural_search_count=1,
    )
    assert pending["status"] == "PENDING_STRUCTURAL_SEARCH"
    with pytest.raises(InputValidationError, match="requires structural evidence"):
        _validate_required_structural_evidence(
            profile=profile,
            structures=(),
            comparisons=(),
            pending_alphafold_request_count=1,
        )
    coordinate = tmp_path / "ready.pdb"
    coordinate.write_text("ATOM\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="coordinate model alone"):
        _validate_required_structural_evidence(
            profile=profile,
            structures=(
                _structure(
                    "p1",
                    "ready",
                    coordinate=coordinate,
                    digest=sha256_file(path=coordinate),
                ),
            ),
            comparisons=(),
        )
    completed_zero_hit_search = _validate_required_structural_evidence(
        profile=profile,
        structures=(),
        comparisons=(),
        completed_structural_search_count=1,
    )
    assert completed_zero_hit_search["status"] == "COMPLETE"
    assert (
        _validate_required_structural_evidence(
            profile=replace(profile, require_structural_evidence=False),
            structures=(),
            comparisons=(),
        )["status"]
        == "NOT_REQUIRED"
    )

    document = _campaign_document(example_dir=example_dir)
    document["inputs"]["structures"] = None
    document["inputs"]["structure_comparisons"] = None
    path = _write_document(tmp_path / "no_structural_evidence.yaml", document)
    with pytest.raises(InputValidationError, match="requires structural evidence"):
        validate_campaign(config_path=path)


def test_validate_campaign_parses_every_optional_evidence_table(
    example_dir: Path,
    tmp_path: Path,
) -> None:
    """Preflight should reject invalid rows in optional feature and comparison TSVs."""

    document = _campaign_document(example_dir=example_dir)
    features = tmp_path / "features.tsv"
    features.write_text(
        "protein_id\tfeature_type\tfeature_id\tfeature_name\tstart\tend\t"
        "evidence_status\tevidence_source\tevidence_reference\tderivation_scope\t"
        "feature_definition_sha256\tderivation_cohort_sha256\n"
        "fbox01\tMOTIF\tM1\tMotif one\t1\t3\tASSESSED_WITH_FEATURE\ttest\t"
        f"row1\tFIXED_EXTERNAL\t{'0' * 64}\t\n",
        encoding="utf-8",
    )
    document["inputs"]["features"] = str(features)
    valid_path = _write_document(tmp_path / "all_optional.yaml", document)
    summary = validate_campaign(config_path=valid_path)
    assert summary["external_feature_record_count"] == 1
    assert summary["domain_hit_count"] == 12
    assert summary["structure_count"] == 8
    assert summary["structure_comparison_count"] == 7

    comparison_path = Path(str(document["inputs"]["structure_comparisons"]))
    invalid_comparisons = tmp_path / "invalid_comparisons.tsv"
    invalid_comparisons.write_text(
        comparison_path.read_text(encoding="utf-8").replace("fbox01", "unknown", 1),
        encoding="utf-8",
    )
    document["inputs"]["structure_comparisons"] = str(invalid_comparisons)
    invalid_path = _write_document(tmp_path / "invalid_optional.yaml", document)
    with pytest.raises(InputValidationError, match="absent from FASTA"):
        validate_campaign(config_path=invalid_path)


def test_validate_campaign_preflights_enabled_foldseek_without_search(
    example_dir: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Campaign validation should inspect Foldseek but never execute easy-search."""

    document = _campaign_document(example_dir=example_dir)
    accessions = tmp_path / "alphafold.tsv"
    accessions.write_text(
        "protein_id\tuniprot_accession\nfbox01\tQ9SA03\nfbox02\tQ9SA04\n",
        encoding="utf-8",
    )
    document["inputs"]["alphafold_accessions"] = str(accessions)
    document["alphafold"]["enabled"] = True
    document["foldseek"]["enabled"] = True
    monkeypatch.setattr(pipeline_module.shutil, "which", lambda _value: "/usr/bin/foldseek")
    monkeypatch.setattr(pipeline_module, "foldseek_version", lambda **_kwargs: "10-test")
    summary = validate_campaign(
        config_path=_write_document(tmp_path / "foldseek_campaign.yaml", document)
    )
    assert summary["foldseek_preflight"]["candidate_model_count"] == 2
    assert summary["foldseek_preflight"]["tool_version"] == "10-test"


def test_foldseek_preflight_checks_tool_candidates_and_full_hit_ceiling(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Read-only Foldseek preflight should fail closed before a structural search."""

    settings = FoldseekSettings(
        enabled=True,
        executable="foldseek",
        cache_dir=tmp_path / "cache",
        e_value_threshold=0.001,
        sensitivity=9.5,
        maximum_hits=2,
    )
    coordinate = tmp_path / "model.pdb"
    coordinate.write_text("ATOM\n", encoding="utf-8")
    structure = _structure(
        "p1",
        "s1",
        coordinate=coordinate,
        digest=sha256_file(path=coordinate),
    )
    monkeypatch.setattr(pipeline_module.shutil, "which", lambda _value: None)
    with pytest.raises(ExternalToolError, match="was not found"):
        _validate_foldseek_preflight(
            settings=settings,
            structures=(structure,),
            alphafold_requests=(AlphaFoldRequest("p2", "Q2"),),
        )

    monkeypatch.setattr(pipeline_module.shutil, "which", lambda _value: "/usr/bin/foldseek")
    monkeypatch.setattr(pipeline_module, "foldseek_version", lambda **_kwargs: "10-test")
    with pytest.raises(InputValidationError, match="fewer than two"):
        _validate_foldseek_preflight(
            settings=settings,
            structures=(structure,),
            alphafold_requests=(),
        )
    with pytest.raises(InputValidationError, match="fewer than two"):
        _validate_foldseek_preflight(
            settings=settings,
            structures=(structure,),
            alphafold_requests=(AlphaFoldRequest("p1", "Q1"),),
        )
    with pytest.raises(InputValidationError, match="maximum_hits"):
        _validate_foldseek_preflight(
            settings=settings,
            structures=(structure,),
            alphafold_requests=(
                AlphaFoldRequest("p2", "Q2"),
                AlphaFoldRequest("p3", "Q3"),
            ),
        )
    summary = _validate_foldseek_preflight(
        settings=replace(settings, maximum_hits=3),
        structures=(structure,),
        alphafold_requests=(
            AlphaFoldRequest("p2", "Q2"),
            AlphaFoldRequest("p3", "Q3"),
        ),
    )
    assert summary == {
        "status": "VALID",
        "executable": "/usr/bin/foldseek",
        "tool_version": "10-test",
        "supplied_candidate_model_count": 1,
        "alphafold_candidate_model_count": 2,
        "candidate_model_count": 3,
        "candidate_protein_count": 3,
        "maximum_hits": 3,
    }


def test_prepare_campaign_covers_acquisition_foldseek_and_structural_import(
    example_dir: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Selected optional structural stages should merge evidence and bind provenance."""

    document = _campaign_document(example_dir=example_dir)
    accessions = tmp_path / "alphafold.tsv"
    accessions.write_text("protein_id\tuniprot_accession\nfbox05\tQ9SA03\n", encoding="utf-8")
    structural_resource = tmp_path / "structural_resource"
    structural_resource.mkdir()
    document["inputs"]["alphafold_accessions"] = str(accessions)
    document["inputs"]["structural_alignment_resource"] = str(structural_resource)
    document["alphafold"]["enabled"] = True
    document["foldseek"]["enabled"] = True
    config = load_config(path=_write_document(tmp_path / "campaign.yaml", document))
    coordinate = tmp_path / "af.pdb"
    coordinate.write_text("ATOM\n", encoding="utf-8")
    digest = sha256_file(path=coordinate)
    acquired = _structure(
        protein_id="fbox05",
        structure_id="AF-Q9SA03-F1",
        coordinate=coordinate,
        digest=digest,
    )
    acquisition = AlphaFoldAcquisition(
        "fbox05",
        "Q9SA03",
        "ACQUIRED",
        "AF-Q9SA03-F1",
        "v4",
        "https://example.invalid",
        coordinate,
        digest,
        True,
        90.0,
        "",
    )
    imported_comparison = _comparison("US-align", "imported")
    imported = StructuralAlignmentImport(
        comparisons=(imported_comparison,),
        features=(
            FeatureRecord(
                "fbox01",
                "STRUCTURAL_POCKET",
                "POCKET",
                "Pocket",
                None,
                None,
                "DERIVED",
                "upstream",
                "row",
            ),
        ),
        group_summaries=(),
        input_paths=(tmp_path / "import_manifest.json",),
        package_version="1.0",
        run_digest="a" * 64,
        comparison_universe_members={"upstream": frozenset({"fbox01", "ctrl01"})},
    )
    raw_output = tmp_path / "foldseek.tsv"
    manifest = tmp_path / "foldseek.json"
    raw_output.write_text("raw", encoding="utf-8")
    manifest.write_text("{}", encoding="utf-8")
    foldseek = FoldseekRunEvidence(
        comparisons=(_comparison("Foldseek", "generated"),),
        raw_output_path=raw_output,
        completion_manifest_path=manifest,
        tool_version="1",
        cache_key="key",
        comparison_universe_id="test_campaign",
        assessment_universe=frozenset({"fbox01", "ctrl01"}),
        reused=False,
    )
    monkeypatch.setattr(
        pipeline_module,
        "acquire_alphafold_models",
        lambda **_kwargs: ((acquisition,), (acquired,)),
    )
    monkeypatch.setattr(
        pipeline_module,
        "import_structural_alignment_resource",
        lambda **_kwargs: imported,
    )
    monkeypatch.setattr(
        pipeline_module,
        "run_foldseek_all_vs_all",
        lambda **_kwargs: foldseek,
    )
    data = _prepare_campaign(config=config, profile=load_profile(source="e3"), threads=2)
    assert data["evidence_availability"]["alphafold_acquisition"] == "COMPLETE"
    assert data["evidence_availability"]["foldseek"] == "COMPLETE"
    assert data["evidence_availability"]["imported_structural_alignment"] == "COMPLETE"
    assert raw_output not in data["input_paths"]
    assert coordinate not in data["input_paths"]
    assert data["asset_sources"]["assets/foldseek/key/foldseek_all_vs_all.tsv"] == raw_output
    assert data["asset_sources"]["assets/foldseek/key/COMPLETED.json"] == manifest
    assert data["foldseek"]["raw_output_path"] == ("assets/foldseek/key/foldseek_all_vs_all.tsv")


def test_explainable_ml_publishes_through_end_to_end_pipeline(
    example_dir: Path,
    tmp_path: Path,
) -> None:
    """An enabled model should publish typed rows into the physical DuckDB result."""

    document = _campaign_document(example_dir=example_dir)
    document["explainable_ml"]["enabled"] = True
    config_path = _write_document(tmp_path / "ml_campaign.yaml", document)
    result = run_campaign(config_path=config_path, output_dir=tmp_path / "result")
    with duckdb.connect(str(result / "protein_signatures.duckdb"), read_only=True) as connection:
        model = connection.execute("SELECT status, validation_roc_auc FROM ml_models").fetchone()
        assert model is not None
        assert str(model[0]).startswith("COMPLETE")
        assert 0.0 <= float(model[1]) <= 1.0
        assert connection.execute("SELECT count(*) FROM ml_predictions").fetchone()[0] == 24


def test_pipeline_merge_helpers_reject_collisions_and_filter_failures(tmp_path: Path) -> None:
    """Structure and feature merges should be unique, ordered and status aware."""

    first = _structure("p1", "s1")
    second = _structure("p2", "s2")
    assert _merge_structures(supplied=(second,), acquired=(first,)) == (first, second)
    with pytest.raises(InputValidationError, match="multiple sources"):
        _merge_structures(supplied=(first,), acquired=(first,))
    comparison = _comparison("Foldseek", "row")
    assert _merge_structure_comparisons(supplied=(comparison,), generated=()) == (comparison,)
    with pytest.raises(InputValidationError, match="more than once"):
        _merge_structure_comparisons(supplied=(comparison,), generated=(comparison,))
    retained = FeatureRecord("p1", "TYPE", "x", "X", None, None, "DERIVED", "test", "x")
    excluded = replace(retained, feature_id="bad", evidence_status="FAILED")
    assert _merge_features(
        external_features=(retained, retained, excluded),
        domain_features=(),
        structure_features=(),
        kmer_features=(),
    ) == (retained,)


def test_structure_asset_suffixes_and_portable_records(tmp_path: Path) -> None:
    """Coordinate assets should use digest names while preserving recognised suffixes."""

    structures = []
    for index, suffix in enumerate((".cif.gz", ".pdb.gz", ".cif", ".ent"), start=1):
        path = tmp_path / f"model_{index}{suffix}"
        path.write_text("ATOM\n", encoding="utf-8")
        structures.append(
            _structure(
                protein_id=f"p{index}",
                structure_id=f"s{index}",
                coordinate=path,
                digest=sha256_file(path=path),
            )
        )
    structures.append(_structure("p5", "s5"))
    sources, names = _structure_asset_sources(structures=tuple(structures))
    assert len(sources) == 4
    assert names["s1"].endswith(".cif.gz")
    assert names["s2"].endswith(".pdb.gz")
    assert names["s3"].endswith(".cif")
    assert names["s4"].endswith(".pdb")
    assert (
        _structure_record(item=structures[0], asset_names=names)["coordinate_path"] == names["s1"]
    )
    acquisition = AlphaFoldAcquisition(
        "p1", "Q9SA03", "COMPLETE", "s1", "v4", "url", None, "", True, 80.0, ""
    )
    assert _alphafold_record(item=acquisition, asset_names={})["coordinate_path"] == ""


def test_input_authorities_resource_and_raw_modes(example_dir: Path) -> None:
    """Input provenance should include the exact selected upstream authority files."""

    config = load_config(path=example_dir / "campaign.yaml")
    resource_file = example_dir / "proteins.faa"
    resource = SimpleNamespace()
    original = pipeline_module.resource_input_paths
    pipeline_module.resource_input_paths = lambda **_kwargs: (resource_file,)
    try:
        paths = _input_authorities(config=config, source=resource, resource_mode=True)
    finally:
        pipeline_module.resource_input_paths = original
    assert resource_file in paths
    raw = SimpleNamespace(
        log_path=example_dir / "campaign.yaml",
        orthogroups_path=example_dir / "domains.tsv",
        hog_paths=(example_dir / "label_assignments.tsv",),
        sequence_ids_path=example_dir / "structures.tsv",
        species_ids_path=example_dir / "structure_comparisons.tsv",
    )
    raw_paths = _input_authorities(config=config, source=raw, resource_mode=False)
    assert raw.log_path in raw_paths
    assert raw.species_ids_path in raw_paths


def test_resource_identity_pfam_states_and_profile_defaults(example_dir: Path) -> None:
    """Resource run selection and Pfam coverage should expose every controlled state."""

    config = load_config(path=example_dir / "campaign.yaml")
    _validate_requested_resource_run(config=config, resource=SimpleNamespace(run_id="run"))
    configured = replace(
        config,
        inputs=replace(config.inputs, orthofinder_run_id="expected"),
    )
    with pytest.raises(InputValidationError, match="does not match"):
        _validate_requested_resource_run(
            config=configured,
            resource=SimpleNamespace(run_id="observed"),
        )
    assert _pfam_availability(assessments=()) == "NOT_ASSESSED"
    complete = (
        DomainAssessment("p1", "Pfam", DomainAssessmentStatus.ASSESSED_WITH_HIT, 1, "test", ""),
    )
    assert _pfam_availability(assessments=complete) == "COMPLETE"
    partial = (
        *complete,
        DomainAssessment("p2", "Pfam", DomainAssessmentStatus.FAILED, 0, "test", ""),
    )
    assert _pfam_availability(assessments=partial) == "PARTIAL"
    profile = load_profile(source="e3")
    defaults = _with_effective_comparisons(config=replace(config, comparisons=()), profile=profile)
    assert len(defaults.comparisons) == 73
    assert _with_effective_comparisons(config=config, profile=profile) is config


def _campaign_document(*, example_dir: Path) -> dict[str, object]:
    """Load the packaged campaign and absolutise every selected file input."""

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


def _write_document(path: Path, document: dict[str, object]) -> Path:
    """Write one YAML campaign document for a test."""

    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def _structure(
    protein_id: str,
    structure_id: str,
    *,
    coordinate: Path | None = None,
    digest: str = "",
) -> StructureRecord:
    """Create one compact structure record."""

    return StructureRecord(
        protein_id,
        structure_id,
        "TEST",
        "1",
        coordinate,
        digest,
        "COMPLETE" if coordinate is not None else "INPUT_UNAVAILABLE",
        80.0,
        "",
        "",
        "",
        "",
        "",
        FoldEvidenceStatus.NOT_ASSESSED,
        (
            StructureAnalysisEligibility.ELIGIBLE
            if coordinate is not None
            else StructureAnalysisEligibility.NOT_APPLICABLE_EXTERNAL_EVIDENCE
        ),
        (),
    )


def _comparison(tool: str, record_id: str) -> PairwiseStructureComparison:
    """Create one compact complete structure comparison."""

    return PairwiseStructureComparison(
        "fbox01",
        "ctrl01",
        tool,
        "1",
        0.8,
        1.0,
        80,
        0.8,
        0.8,
        StructureComparisonStatus.COMPLETE,
        record_id,
        "test_campaign",
        StructureCoverageScope.FULL_SEQUENCE,
    )
