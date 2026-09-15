"""Campaign configuration loading, normalisation and validation."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

from .errors import ConfigurationError, InputValidationError
from .models import (
    AlphaFoldSettings,
    AnalysisSettings,
    CampaignConfig,
    ComparisonDefinition,
    ExplainableMLSettings,
    FoldseekSettings,
    InputPaths,
)
from .validation import (
    parse_optional_float,
    parse_optional_integer,
    reject_unknown_fields,
    require_mapping,
    require_sequence,
    validate_identifier,
    validate_text,
)

LOGGER = logging.getLogger(__name__)


def load_config(*, path: Path) -> CampaignConfig:
    """Load and validate one YAML campaign configuration.

    Args:
        path: Campaign YAML path.

    Returns:
        Validated, path-resolved configuration.

    Raises:
        ConfigurationError: If YAML cannot be parsed.
        InputValidationError: If the configuration contract is violated.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty campaign configuration: {source}")
    try:
        with source.open(mode="r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle)
    except yaml.YAMLError as error:
        raise ConfigurationError(f"Invalid YAML in {source}: {error}") from error
    document = require_mapping(value=raw, field_name="configuration")
    reject_unknown_fields(
        value=document,
        allowed=frozenset(
            {
                "schema_version",
                "campaign",
                "inputs",
                "analysis",
                "alphafold",
                "foldseek",
                "explainable_ml",
                "comparisons",
            }
        ),
        field_name="configuration",
    )
    schema_version = parse_optional_integer(
        value=document.get("schema_version"), field_name="schema_version", minimum=1
    )
    if schema_version != 1:
        raise ConfigurationError(f"Unsupported schema_version: {schema_version!r}; expected 1.")
    campaign = require_mapping(value=document.get("campaign"), field_name="campaign")
    reject_unknown_fields(
        value=campaign,
        allowed=frozenset({"campaign_id", "profile"}),
        field_name="campaign",
    )
    inputs = _parse_inputs(value=document.get("inputs"), base=source.parent)
    analysis = _parse_analysis(value=document.get("analysis", {}))
    alphafold = _parse_alphafold(value=document.get("alphafold", {}), base=source.parent)
    foldseek = _parse_foldseek(value=document.get("foldseek", {}), base=source.parent)
    explainable_ml = _parse_explainable_ml(
        value=document.get("explainable_ml", {}), base=source.parent
    )
    comparisons_value = document.get("comparisons", "profile_defaults")
    if comparisons_value == "profile_defaults":
        comparisons = ()
    else:
        comparison_values = require_sequence(value=comparisons_value, field_name="comparisons")
        if not comparison_values:
            raise ConfigurationError("comparisons must be 'profile_defaults' or a non-empty list.")
        comparisons = tuple(
            _parse_comparison(value=value, index=index)
            for index, value in enumerate(comparison_values)
        )
    _validate_comparisons(comparisons=comparisons)
    config = CampaignConfig(
        schema_version=schema_version,
        campaign_id=validate_identifier(
            value=campaign.get("campaign_id"), field_name="campaign.campaign_id"
        ),
        profile_name=_resolve_profile_source(
            value=campaign.get("profile", "e3"), base=source.parent
        ),
        config_path=source,
        inputs=inputs,
        comparisons=comparisons,
        analysis=analysis,
        alphafold=alphafold,
        foldseek=foldseek,
        explainable_ml=explainable_ml,
    )
    LOGGER.info(
        "Loaded campaign %s with %d explicit comparisons",
        config.campaign_id,
        len(config.comparisons),
    )
    return config


def config_to_record(*, config: CampaignConfig) -> dict[str, Any]:
    """Return a deterministic JSON-compatible configuration record.

    Args:
        config: Validated campaign configuration.

    Returns:
        Nested plain mapping with absolute input paths.
    """

    inputs = {
        key: (str(value) if isinstance(value, Path) else value)
        for key, value in vars(config.inputs).items()
    }
    comparisons = [
        {
            "comparison_id": item.comparison_id,
            "display_name": item.display_name,
            "target_label_ids": list(item.target_label_ids),
            "background_label_ids": list(item.background_label_ids),
            "description": item.description,
        }
        for item in config.comparisons
    ]
    return {
        "schema_version": config.schema_version,
        "campaign": {
            "campaign_id": config.campaign_id,
            "profile": config.profile_name,
        },
        "inputs": inputs,
        "analysis": {
            **vars(config.analysis),
            "kmer_lengths": list(config.analysis.kmer_lengths),
        },
        "alphafold": {
            **vars(config.alphafold),
            "cache_dir": str(config.alphafold.cache_dir),
        },
        "foldseek": {
            **vars(config.foldseek),
            "cache_dir": str(config.foldseek.cache_dir),
        },
        "explainable_ml": {
            **vars(config.explainable_ml),
            "plot_cache_dir": str(config.explainable_ml.plot_cache_dir),
            "regularisation_strengths": list(config.explainable_ml.regularisation_strengths),
        },
        "comparisons": comparisons,
    }


def _parse_inputs(*, value: Any, base: Path) -> InputPaths:
    """Parse and resolve the input section.

    Args:
        value: Raw input mapping.
        base: Configuration directory for relative paths.

    Returns:
        Resolved input paths.
    """

    row = require_mapping(value=value, field_name="inputs")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {
                "sequences_fasta",
                "label_assignments",
                "label_evidence_marker",
                "label_evidence_audit",
                "control_matching_audit",
                "label_definition_features",
                "class_labelling_summary",
                "unresolved_assignments",
                "features",
                "domains",
                "redundancy_clusters",
                "structures",
                "structure_comparisons",
                "structural_alignment_resource",
                "alphafold_accessions",
                "orthofinder",
            }
        ),
        field_name="inputs",
    )
    orthofinder = require_mapping(value=row.get("orthofinder", {}), field_name="inputs.orthofinder")
    reject_unknown_fields(
        value=orthofinder,
        allowed=frozenset(
            {"resource_dir", "results_dir", "group_type", "hierarchy_node", "run_id"}
        ),
        field_name="inputs.orthofinder",
    )
    resource_dir = _resolve_optional_directory(value=orthofinder.get("resource_dir"), base=base)
    results_dir = _resolve_optional_directory(value=orthofinder.get("results_dir"), base=base)
    if resource_dir is not None and results_dir is not None:
        raise ConfigurationError(
            "inputs.orthofinder.resource_dir and results_dir are mutually exclusive."
        )
    group_type = validate_text(
        value=orthofinder.get("group_type", "HOG"),
        field_name="inputs.orthofinder.group_type",
    ).upper()
    if group_type not in {"HOG", "LEGACY_ORTHOGROUP"}:
        raise ConfigurationError("inputs.orthofinder.group_type must be HOG or LEGACY_ORTHOGROUP.")
    hierarchy_node = validate_text(
        value=orthofinder.get("hierarchy_node", "N0" if group_type == "HOG" else ""),
        field_name="inputs.orthofinder.hierarchy_node",
        allow_empty=group_type == "LEGACY_ORTHOGROUP",
    )
    if group_type == "HOG":
        hierarchy_node = validate_identifier(
            value=hierarchy_node,
            field_name="inputs.orthofinder.hierarchy_node",
        )
    elif hierarchy_node:
        raise ConfigurationError(
            "inputs.orthofinder.hierarchy_node must be empty for LEGACY_ORTHOGROUP."
        )
    run_id = validate_text(
        value=orthofinder.get("run_id", ""),
        field_name="inputs.orthofinder.run_id",
        allow_empty=True,
    )
    if run_id:
        run_id = validate_identifier(
            value=run_id,
            field_name="inputs.orthofinder.run_id",
        )
    return InputPaths(
        sequences_fasta=_resolve_required_path(
            value=row.get("sequences_fasta"), base=base, field_name="inputs.sequences_fasta"
        ),
        label_assignments=_resolve_required_path(
            value=row.get("label_assignments"),
            base=base,
            field_name="inputs.label_assignments",
        ),
        label_evidence_marker=_resolve_optional_path(
            value=row.get("label_evidence_marker"), base=base
        ),
        label_evidence_audit=_resolve_optional_path(
            value=row.get("label_evidence_audit"), base=base
        ),
        control_matching_audit=_resolve_optional_path(
            value=row.get("control_matching_audit"), base=base
        ),
        label_definition_features=_resolve_optional_path(
            value=row.get("label_definition_features"), base=base
        ),
        class_labelling_summary=_resolve_optional_path(
            value=row.get("class_labelling_summary"), base=base
        ),
        unresolved_assignments=_resolve_optional_path(
            value=row.get("unresolved_assignments"), base=base
        ),
        features=_resolve_optional_path(value=row.get("features"), base=base),
        domains=_resolve_optional_path(value=row.get("domains"), base=base),
        redundancy_clusters=_resolve_optional_path(value=row.get("redundancy_clusters"), base=base),
        structures=_resolve_optional_path(value=row.get("structures"), base=base),
        structure_comparisons=_resolve_optional_path(
            value=row.get("structure_comparisons"), base=base
        ),
        structural_alignment_resource=_resolve_optional_directory(
            value=row.get("structural_alignment_resource"), base=base
        ),
        alphafold_accessions=_resolve_optional_path(
            value=row.get("alphafold_accessions"), base=base
        ),
        orthofinder_resource=resource_dir,
        orthofinder_results=results_dir,
        orthofinder_group_type=group_type,
        orthofinder_hierarchy_node=hierarchy_node,
        orthofinder_run_id=run_id,
    )


def _parse_analysis(*, value: Any) -> AnalysisSettings:
    """Parse bounded analysis settings.

    Args:
        value: Raw settings mapping.

    Returns:
        Validated analysis settings.
    """

    row = require_mapping(value=value, field_name="analysis")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {
                "kmer_lengths",
                "minimum_target_proteins",
                "minimum_background_proteins",
                "minimum_feature_proteins",
                "maximum_kmer_features",
                "validation_fraction",
                "fdr_threshold",
                "random_seed",
                "structural_tm_score_threshold",
                "structural_minimum_coverage",
            }
        ),
        field_name="analysis",
    )
    lengths = tuple(
        parse_optional_integer(value=item, field_name="analysis.kmer_lengths", minimum=1)
        for item in require_sequence(
            value=row.get("kmer_lengths", [3, 4]), field_name="analysis.kmer_lengths"
        )
    )
    if not lengths or any(length is None or length > 12 for length in lengths):
        raise ConfigurationError("analysis.kmer_lengths must contain integers from 1 to 12.")
    integer_lengths = tuple(int(length) for length in lengths)
    if len(integer_lengths) != len(set(integer_lengths)):
        raise ConfigurationError("analysis.kmer_lengths must not contain duplicates.")
    validation_fraction = parse_optional_float(
        value=row.get("validation_fraction", 0.20),
        field_name="analysis.validation_fraction",
        minimum=0.0,
        maximum=0.9,
    )
    fdr_threshold = parse_optional_float(
        value=row.get("fdr_threshold", 0.05),
        field_name="analysis.fdr_threshold",
        minimum=0.0,
        maximum=1.0,
    )
    return AnalysisSettings(
        kmer_lengths=integer_lengths,
        minimum_target_proteins=int(
            parse_optional_integer(
                value=row.get("minimum_target_proteins", 3),
                field_name="analysis.minimum_target_proteins",
                minimum=1,
            )
        ),
        minimum_background_proteins=int(
            parse_optional_integer(
                value=row.get("minimum_background_proteins", 3),
                field_name="analysis.minimum_background_proteins",
                minimum=1,
            )
        ),
        minimum_feature_proteins=int(
            parse_optional_integer(
                value=row.get("minimum_feature_proteins", 2),
                field_name="analysis.minimum_feature_proteins",
                minimum=1,
            )
        ),
        maximum_kmer_features=int(
            parse_optional_integer(
                value=row.get("maximum_kmer_features", 250_000),
                field_name="analysis.maximum_kmer_features",
                minimum=1,
            )
        ),
        validation_fraction=float(validation_fraction),
        fdr_threshold=float(fdr_threshold),
        random_seed=int(
            parse_optional_integer(
                value=row.get("random_seed", 1729),
                field_name="analysis.random_seed",
                minimum=0,
            )
        ),
        structural_tm_score_threshold=float(
            parse_optional_float(
                value=row.get("structural_tm_score_threshold", 0.50),
                field_name="analysis.structural_tm_score_threshold",
                minimum=0.0,
                maximum=1.0,
            )
        ),
        structural_minimum_coverage=float(
            parse_optional_float(
                value=row.get("structural_minimum_coverage", 0.50),
                field_name="analysis.structural_minimum_coverage",
                minimum=0.0,
                maximum=1.0,
            )
        ),
    )


def _parse_alphafold(*, value: Any, base: Path) -> AlphaFoldSettings:
    """Parse AlphaFold Database acquisition settings.

    Args:
        value: Raw AlphaFold mapping.
        base: Configuration directory for relative cache paths.

    Returns:
        Validated acquisition settings.
    """

    row = require_mapping(value=value, field_name="alphafold")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {"enabled", "cache_dir", "timeout_seconds", "retries", "minimum_mean_plddt"}
        ),
        field_name="alphafold",
    )
    enabled_value = row.get("enabled", False)
    if not isinstance(enabled_value, bool):
        raise ConfigurationError("alphafold.enabled must be true or false.")
    raw_cache = validate_text(
        value=row.get("cache_dir", ".protein_signature_cache/alphafold"),
        field_name="alphafold.cache_dir",
    )
    cache = Path(raw_cache).expanduser()
    if not cache.is_absolute():
        cache = base / cache
    return AlphaFoldSettings(
        enabled=enabled_value,
        cache_dir=cache.resolve(),
        timeout_seconds=float(
            parse_optional_float(
                value=row.get("timeout_seconds", 60.0),
                field_name="alphafold.timeout_seconds",
                minimum=1.0,
                maximum=600.0,
            )
        ),
        retries=int(
            parse_optional_integer(
                value=row.get("retries", 3), field_name="alphafold.retries", minimum=0
            )
        ),
        minimum_mean_plddt=float(
            parse_optional_float(
                value=row.get("minimum_mean_plddt", 50.0),
                field_name="alphafold.minimum_mean_plddt",
                minimum=0.0,
                maximum=100.0,
            )
        ),
    )


def _parse_foldseek(*, value: Any, base: Path) -> FoldseekSettings:
    """Parse optional Foldseek structural-analysis settings.

    Args:
        value: Raw Foldseek mapping.
        base: Configuration directory for relative cache paths.

    Returns:
        Validated Foldseek settings.
    """

    row = require_mapping(value=value, field_name="foldseek")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {
                "enabled",
                "executable",
                "cache_dir",
                "e_value_threshold",
                "sensitivity",
                "maximum_hits",
            }
        ),
        field_name="foldseek",
    )
    enabled_value = row.get("enabled", False)
    if not isinstance(enabled_value, bool):
        raise ConfigurationError("foldseek.enabled must be true or false.")
    executable = validate_text(
        value=row.get("executable", "foldseek"), field_name="foldseek.executable"
    )
    if Path(executable).name != executable and not Path(executable).is_absolute():
        raise ConfigurationError("foldseek.executable must be a command name or an absolute path.")
    raw_cache = validate_text(
        value=row.get("cache_dir", ".protein_signature_cache/foldseek"),
        field_name="foldseek.cache_dir",
    )
    cache = Path(raw_cache).expanduser()
    if not cache.is_absolute():
        cache = base / cache
    return FoldseekSettings(
        enabled=enabled_value,
        executable=executable,
        cache_dir=cache.resolve(),
        e_value_threshold=float(
            parse_optional_float(
                value=row.get("e_value_threshold", 0.001),
                field_name="foldseek.e_value_threshold",
                minimum=0.0,
            )
        ),
        sensitivity=float(
            parse_optional_float(
                value=row.get("sensitivity", 9.5),
                field_name="foldseek.sensitivity",
                minimum=1.0,
                maximum=10.0,
            )
        ),
        maximum_hits=int(
            parse_optional_integer(
                value=row.get("maximum_hits", 1000),
                field_name="foldseek.maximum_hits",
                minimum=1,
            )
        ),
    )


def _parse_explainable_ml(*, value: Any, base: Path) -> ExplainableMLSettings:
    """Parse mandatory, leakage-safe explainable-classification settings.

    Args:
        value: Raw explainable-ML mapping.
        base: Configuration directory for relative plot-cache paths.

    Returns:
        Validated modelling settings.

    Raises:
        ConfigurationError: If Boolean flags or regularisation values are invalid.
    """

    row = require_mapping(value=value, field_name="explainable_ml")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {
                "enabled",
                "minimum_samples_per_class",
                "minimum_groups_per_class",
                "maximum_features",
                "regularisation_strengths",
                "l1_ratio",
                "cross_validation_folds",
                "maximum_iterations",
                "convergence_tolerance",
                "permutation_repeats",
                "maximum_permutation_features",
                "top_local_explanations",
                "maximum_waterfall_plots",
                "plot_cache_dir",
                "exclude_technical_features",
            }
        ),
        field_name="explainable_ml",
    )
    enabled = row.get("enabled", True)
    exclude_technical = row.get("exclude_technical_features", True)
    if not isinstance(enabled, bool):
        raise ConfigurationError("explainable_ml.enabled must be true or false.")
    if not enabled:
        raise ConfigurationError(
            "explainable_ml is a mandatory analysis stage; enabled must be true."
        )
    if not isinstance(exclude_technical, bool):
        raise ConfigurationError("explainable_ml.exclude_technical_features must be true or false.")
    strengths = tuple(
        float(
            parse_optional_float(
                value=item,
                field_name="explainable_ml.regularisation_strengths",
                minimum=1e-12,
            )
        )
        for item in require_sequence(
            value=row.get("regularisation_strengths", [0.01, 0.1, 1.0]),
            field_name="explainable_ml.regularisation_strengths",
        )
    )
    if not strengths or len(strengths) != len(set(strengths)):
        raise ConfigurationError(
            "explainable_ml.regularisation_strengths must be a non-empty unique list."
        )
    raw_plot_cache = validate_text(
        value=row.get("plot_cache_dir", ".protein_signature_cache/shap"),
        field_name="explainable_ml.plot_cache_dir",
    )
    plot_cache = Path(raw_plot_cache).expanduser()
    if not plot_cache.is_absolute():
        plot_cache = base / plot_cache
    return ExplainableMLSettings(
        enabled=enabled,
        minimum_samples_per_class=int(
            parse_optional_integer(
                value=row.get("minimum_samples_per_class", 10),
                field_name="explainable_ml.minimum_samples_per_class",
                minimum=2,
            )
        ),
        minimum_groups_per_class=int(
            parse_optional_integer(
                value=row.get("minimum_groups_per_class", 3),
                field_name="explainable_ml.minimum_groups_per_class",
                minimum=2,
            )
        ),
        maximum_features=int(
            parse_optional_integer(
                value=row.get("maximum_features", 500),
                field_name="explainable_ml.maximum_features",
                minimum=1,
            )
        ),
        regularisation_strengths=tuple(sorted(strengths)),
        l1_ratio=float(
            parse_optional_float(
                value=row.get("l1_ratio", 0.2),
                field_name="explainable_ml.l1_ratio",
                minimum=0.0,
                maximum=1.0,
            )
        ),
        cross_validation_folds=int(
            parse_optional_integer(
                value=row.get("cross_validation_folds", 5),
                field_name="explainable_ml.cross_validation_folds",
                minimum=2,
            )
        ),
        maximum_iterations=int(
            parse_optional_integer(
                value=row.get("maximum_iterations", 2000),
                field_name="explainable_ml.maximum_iterations",
                minimum=50,
            )
        ),
        convergence_tolerance=float(
            parse_optional_float(
                value=row.get("convergence_tolerance", 1e-4),
                field_name="explainable_ml.convergence_tolerance",
                minimum=1e-12,
                maximum=0.1,
            )
        ),
        permutation_repeats=int(
            parse_optional_integer(
                value=row.get("permutation_repeats", 10),
                field_name="explainable_ml.permutation_repeats",
                minimum=1,
            )
        ),
        maximum_permutation_features=int(
            parse_optional_integer(
                value=row.get("maximum_permutation_features", 50),
                field_name="explainable_ml.maximum_permutation_features",
                minimum=1,
            )
        ),
        top_local_explanations=int(
            parse_optional_integer(
                value=row.get("top_local_explanations", 10),
                field_name="explainable_ml.top_local_explanations",
                minimum=1,
            )
        ),
        maximum_waterfall_plots=int(
            parse_optional_integer(
                value=row.get("maximum_waterfall_plots", 8),
                field_name="explainable_ml.maximum_waterfall_plots",
                minimum=0,
            )
        ),
        plot_cache_dir=plot_cache.resolve(),
        exclude_technical_features=exclude_technical,
    )


def _parse_comparison(*, value: Any, index: int) -> ComparisonDefinition:
    """Parse one explicit target/background definition.

    Args:
        value: Raw comparison mapping.
        index: Zero-based index for diagnostics.

    Returns:
        Validated comparison.
    """

    row = require_mapping(value=value, field_name=f"comparisons[{index}]")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {
                "comparison_id",
                "display_name",
                "target_label_ids",
                "background_label_ids",
                "description",
            }
        ),
        field_name=f"comparisons[{index}]",
    )
    target = tuple(
        validate_identifier(value=item, field_name=f"comparisons[{index}].target_label_ids")
        for item in require_sequence(
            value=row.get("target_label_ids"),
            field_name=f"comparisons[{index}].target_label_ids",
        )
    )
    background = tuple(
        validate_identifier(value=item, field_name=f"comparisons[{index}].background_label_ids")
        for item in require_sequence(
            value=row.get("background_label_ids"),
            field_name=f"comparisons[{index}].background_label_ids",
        )
    )
    if not target or not background:
        raise ConfigurationError(
            f"comparisons[{index}] requires target_label_ids and background_label_ids."
        )
    if set(target) & set(background):
        raise ConfigurationError(f"comparisons[{index}] has overlapping target/background labels.")
    comparison_id = validate_identifier(
        value=row.get("comparison_id"), field_name=f"comparisons[{index}].comparison_id"
    )
    return ComparisonDefinition(
        comparison_id=comparison_id,
        display_name=validate_text(
            value=row.get("display_name", comparison_id),
            field_name=f"comparisons[{index}].display_name",
        ),
        target_label_ids=target,
        background_label_ids=background,
        description=validate_text(
            value=row.get("description", ""),
            field_name=f"comparisons[{index}].description",
            allow_empty=True,
        ),
    )


def _validate_comparisons(*, comparisons: tuple[ComparisonDefinition, ...]) -> None:
    """Reject missing or duplicate comparison identifiers.

    Args:
        comparisons: Parsed comparisons.

    Raises:
        ConfigurationError: If comparison identifiers are invalid.
    """

    identifiers = [item.comparison_id for item in comparisons]
    if len(identifiers) != len(set(identifiers)):
        raise ConfigurationError("comparison_id values must be unique.")


def _resolve_profile_source(*, value: Any, base: Path) -> str:
    """Resolve a built-in profile name or configuration-relative YAML path.

    Args:
        value: Raw profile source.
        base: Configuration directory.

    Returns:
        Built-in short name or absolute profile path text.
    """

    source = validate_text(value=value, field_name="campaign.profile")
    if "/" not in source and "\\" not in source and not source.endswith((".yaml", ".yml")):
        return source
    candidate = Path(source).expanduser()
    if not candidate.is_absolute():
        candidate = base / candidate
    return str(candidate.resolve())


def _resolve_required_path(*, value: Any, base: Path, field_name: str) -> Path:
    """Resolve and validate a required file path.

    Args:
        value: Raw path.
        base: Base directory for relative paths.
        field_name: Field name for diagnostics.

    Returns:
        Existing non-empty file.
    """

    candidate = _resolve_path(value=value, base=base, field_name=field_name)
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"{field_name} is missing or empty: {candidate}")
    return candidate


def _resolve_optional_path(*, value: Any, base: Path) -> Path | None:
    """Resolve an optional existing non-empty file.

    Args:
        value: Raw path or blank value.
        base: Base directory for relative paths.

    Returns:
        Existing path or ``None``.
    """

    if value is None or str(value).strip() == "":
        return None
    candidate = _resolve_path(value=value, base=base, field_name="optional input")
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise InputValidationError(f"Optional input is missing or empty: {candidate}")
    return candidate


def _resolve_optional_directory(*, value: Any, base: Path) -> Path | None:
    """Resolve an optional existing directory.

    Args:
        value: Raw directory or blank value.
        base: Base directory for relative paths.

    Returns:
        Existing directory or ``None``.
    """

    if value is None or str(value).strip() == "":
        return None
    candidate = _resolve_path(value=value, base=base, field_name="optional directory")
    if not candidate.is_dir():
        raise InputValidationError(f"Optional input directory does not exist: {candidate}")
    return candidate


def _resolve_path(*, value: Any, base: Path, field_name: str) -> Path:
    """Resolve an absolute or configuration-relative path.

    Args:
        value: Raw path.
        base: Base directory for relative paths.
        field_name: Field name for diagnostics.

    Returns:
        Resolved path without requiring existence.
    """

    text = validate_text(value=value, field_name=field_name)
    candidate = Path(text).expanduser()
    if not candidate.is_absolute():
        candidate = base / candidate
    return candidate.resolve()
