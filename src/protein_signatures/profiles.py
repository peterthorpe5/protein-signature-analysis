"""Load and validate data-driven protein-type classification profiles."""

from __future__ import annotations

import logging
from collections import defaultdict
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml

from .errors import ConfigurationError, InputValidationError
from .models import ComparisonDefinition, LabelAssignment, ProfileLabel, ProteinProfile
from .validation import (
    reject_unknown_fields,
    require_mapping,
    require_sequence,
    validate_identifier,
    validate_text,
)

LOGGER = logging.getLogger(__name__)
_NON_REVIEWED_POSITIVE_EVIDENCE_STATES = frozenset(
    {
        "AMBIGUOUS",
        "EXCLUDED",
        "FAILED",
        "NOT_ASSESSED",
        "NOT_REVIEWED",
        "PENDING",
        "PROPOSED",
        "UNKNOWN",
        "UNMAPPED",
        "UNREVIEWED",
    }
)


def built_in_profile_path(*, profile_name: str) -> Path:
    """Resolve one packaged profile by its short name.

    Args:
        profile_name: Packaged profile identifier, for example ``e3``.

    Returns:
        Filesystem path to the packaged YAML profile.

    Raises:
        ConfigurationError: If the named profile is unavailable.
    """

    name = validate_identifier(value=profile_name, field_name="profile name")
    candidate = files("protein_signatures").joinpath("data", "profiles", f"{name}.yaml")
    if not candidate.is_file():
        raise ConfigurationError(f"Unknown built-in protein profile: {name!r}")
    return Path(str(candidate)).resolve()


def load_profile(*, source: str | Path) -> ProteinProfile:
    """Load a built-in or user-supplied protein classification hierarchy.

    Args:
        source: Built-in short name or YAML path.

    Returns:
        Validated immutable profile.

    Raises:
        ConfigurationError: If YAML parsing fails.
        InputValidationError: If profile content is invalid.
    """

    path = (
        built_in_profile_path(profile_name=source)
        if isinstance(source, str) and "/" not in source and "\\" not in source
        else Path(source).expanduser().resolve()
    )
    if not path.is_file() or path.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty profile YAML: {path}")
    try:
        with path.open(mode="r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle)
    except yaml.YAMLError as error:
        raise ConfigurationError(f"Invalid profile YAML in {path}: {error}") from error
    document = require_mapping(value=raw, field_name="profile")
    reject_unknown_fields(
        value=document,
        allowed=frozenset(
            {
                "profile_id",
                "profile_version",
                "display_name",
                "require_structural_evidence",
                "default_comparison",
                "labels",
            }
        ),
        field_name="profile",
    )
    labels_raw = require_sequence(value=document.get("labels"), field_name="profile.labels")
    labels = tuple(_parse_label(value=value, index=index) for index, value in enumerate(labels_raw))
    default_comparison = require_mapping(
        value=document.get("default_comparison", {}),
        field_name="profile.default_comparison",
    )
    reject_unknown_fields(
        value=default_comparison,
        allowed=frozenset(
            {
                "target_root_label_id",
                "background_label_id",
                "excluded_label_ids",
                "excluded_subtree_label_ids",
                "description",
            }
        ),
        field_name="profile.default_comparison",
    )
    require_structural_evidence = document.get("require_structural_evidence", False)
    if not isinstance(require_structural_evidence, bool):
        raise InputValidationError("profile.require_structural_evidence must be true or false.")
    profile = ProteinProfile(
        profile_id=validate_identifier(
            value=document.get("profile_id"), field_name="profile.profile_id"
        ),
        profile_version=validate_text(
            value=document.get("profile_version"), field_name="profile.profile_version"
        ),
        display_name=validate_text(
            value=document.get("display_name"), field_name="profile.display_name"
        ),
        require_structural_evidence=require_structural_evidence,
        labels=labels,
        default_target_root_label_id=validate_text(
            value=default_comparison.get("target_root_label_id", ""),
            field_name="profile.default_comparison.target_root_label_id",
            allow_empty=True,
        ),
        default_background_label_id=validate_text(
            value=default_comparison.get("background_label_id", ""),
            field_name="profile.default_comparison.background_label_id",
            allow_empty=True,
        ),
        default_comparison_description=validate_text(
            value=default_comparison.get(
                "description",
                "Profile-defined target-versus-background comparison.",
            ),
            field_name="profile.default_comparison.description",
        ),
        default_excluded_label_ids=tuple(
            validate_identifier(
                value=item,
                field_name="profile.default_comparison.excluded_label_ids",
            )
            for item in require_sequence(
                value=default_comparison.get("excluded_label_ids", []),
                field_name="profile.default_comparison.excluded_label_ids",
            )
        ),
        default_excluded_subtree_label_ids=tuple(
            validate_identifier(
                value=item,
                field_name="profile.default_comparison.excluded_subtree_label_ids",
            )
            for item in require_sequence(
                value=default_comparison.get("excluded_subtree_label_ids", []),
                field_name="profile.default_comparison.excluded_subtree_label_ids",
            )
        ),
    )
    _validate_hierarchy(profile=profile)
    LOGGER.info(
        "Loaded protein profile %s version %s with %d labels",
        profile.profile_id,
        profile.profile_version,
        len(profile.labels),
    )
    return profile


def resolve_label(*, profile: ProteinProfile, term: str) -> str:
    """Resolve a canonical profile label or case-insensitive alias.

    Args:
        profile: Validated classification profile.
        term: Canonical label, display name or alias.

    Returns:
        Canonical label identifier.

    Raises:
        InputValidationError: If the term is unknown.
    """

    query = validate_text(value=term, field_name="label term").casefold()
    for label in profile.labels:
        if query in {
            label.label_id.casefold(),
            label.display_name.casefold(),
            *(alias.casefold() for alias in label.aliases),
        }:
            return label.label_id
    raise InputValidationError(f"Unknown label or alias in profile: {term!r}")


def label_ancestors(*, profile: ProteinProfile, label_id: str) -> tuple[str, ...]:
    """Return a label followed by each ancestor up to the root.

    Args:
        profile: Validated classification profile.
        label_id: Canonical label identifier.

    Returns:
        Ordered label and ancestor identifiers.
    """

    by_id = {label.label_id: label for label in profile.labels}
    if label_id not in by_id:
        raise InputValidationError(f"Unknown profile label: {label_id!r}")
    result: list[str] = []
    current = by_id[label_id]
    while True:
        result.append(current.label_id)
        if not current.parent_label_id:
            return tuple(result)
        current = by_id[current.parent_label_id]


def expand_positive_memberships(
    *, assignments: tuple[LabelAssignment, ...], profile: ProteinProfile
) -> tuple[dict[str, str], ...]:
    """Expand reviewed positive labels to their inherited ancestors.

    Args:
        assignments: Evidence-bearing direct assignments.
        profile: Validated classification profile.

    Returns:
        Unique protein/label membership rows with direct/inherited provenance.
    """

    rows: dict[tuple[str, str], dict[str, str]] = {}
    for assignment in assignments:
        if not assignment.is_eligible_positive:
            continue
        for index, label_id in enumerate(
            label_ancestors(profile=profile, label_id=assignment.label_id)
        ):
            key = (assignment.protein_id, label_id)
            source = "DIRECT" if index == 0 else "INHERITED"
            candidate = {
                "protein_id": assignment.protein_id,
                "label_id": label_id,
                "membership_source": source,
                "direct_label_id": assignment.label_id,
            }
            existing = rows.get(key)
            if existing is None or (
                existing["membership_source"] == "INHERITED" and source == "DIRECT"
            ):
                rows[key] = candidate
    return tuple(rows[key] for key in sorted(rows))


def validate_comparison_labels(
    *, comparisons: tuple[ComparisonDefinition, ...], profile: ProteinProfile
) -> None:
    """Require every comparison label to exist in the selected profile.

    Args:
        comparisons: Explicit scientific comparisons.
        profile: Selected classification profile.

    Raises:
        InputValidationError: If a comparison references an unknown label.
    """

    known = profile.label_ids()
    for comparison in comparisons:
        unknown = (set(comparison.target_label_ids) | set(comparison.background_label_ids)) - known
        if unknown:
            raise InputValidationError(
                f"Comparison {comparison.comparison_id!r} references unknown labels: "
                f"{sorted(unknown)}"
            )


def validate_assignment_profile_compatibility(
    *, assignments: tuple[LabelAssignment, ...], profile: ProteinProfile
) -> None:
    """Validate positive-assignment roles and mutually exclusive profile axes.

    Args:
        assignments: Direct evidence-bearing protein-to-label assignments.
        profile: Selected classification profile.

    Raises:
        InputValidationError: If an analysis-positive assignment is disallowed, has the
            wrong component role, or occupies conflicting exclusive strata.
    """

    by_id = {label.label_id: label for label in profile.labels}
    group_owners: dict[tuple[str, str], set[str]] = defaultdict(set)
    for assignment in assignments:
        if assignment.label_id not in by_id:
            raise InputValidationError(
                f"Assignment references unknown profile label: {assignment.label_id!r}"
            )
        if not assignment.is_eligible_positive:
            continue
        label = by_id[assignment.label_id]
        if not label.reviewed_positive_allowed:
            raise InputValidationError(
                f"Label {label.label_id!r} does not permit analysis-positive assignments."
            )
        if assignment.evidence_status.strip().upper() in (_NON_REVIEWED_POSITIVE_EVIDENCE_STATES):
            raise InputValidationError(
                f"Protein {assignment.protein_id!r} has analysis-positive curation "
                f"but incompatible evidence_status {assignment.evidence_status!r}."
            )
        expected_role = label.component_role
        if (
            expected_role
            and expected_role != "UNKNOWN"
            and assignment.component_role != expected_role
        ):
            raise InputValidationError(
                f"Protein {assignment.protein_id!r} has component_role "
                f"{assignment.component_role!r} for {label.label_id!r}; expected "
                f"{expected_role!r}."
            )
        for ancestor_id in label_ancestors(profile=profile, label_id=label.label_id):
            ancestor = by_id[ancestor_id]
            if ancestor.assignment_exclusivity_group:
                group_owners[(assignment.protein_id, ancestor.assignment_exclusivity_group)].add(
                    ancestor.label_id
                )
    for (protein_id, group), owners in sorted(group_owners.items()):
        if len(owners) > 1:
            raise InputValidationError(
                f"Protein {protein_id!r} has conflicting reviewed-positive labels in "
                f"exclusivity group {group!r}: {sorted(owners)}"
            )


def default_profile_comparisons(*, profile: ProteinProfile) -> tuple[ComparisonDefinition, ...]:
    """Create separate comparisons from a profile's declared default policy.

    Args:
        profile: Selected classification profile.

    Returns:
        Target-label versus background-label comparison definitions.

    Raises:
        InputValidationError: If the profile does not declare a valid policy.
    """

    label_ids = profile.label_ids()
    root = profile.default_target_root_label_id
    fallback_background = profile.default_background_label_id
    if not root or not fallback_background:
        raise InputValidationError(
            "profile_defaults requires default_comparison.target_root_label_id and "
            "background_label_id in the selected profile."
        )
    policy_labels = {
        root,
        fallback_background,
        *profile.default_excluded_label_ids,
        *profile.default_excluded_subtree_label_ids,
        *(
            label.default_background_label_id
            for label in profile.labels
            if label.default_background_label_id
        ),
    }
    unknown = policy_labels - label_ids
    if unknown:
        raise InputValidationError(
            f"The profile default-comparison policy references unknown labels: {sorted(unknown)}"
        )
    excluded = set(profile.default_excluded_label_ids) | {root, fallback_background}

    def in_excluded_subtree(label_id: str) -> bool:
        """Return whether a label descends from an excluded policy root."""

        ancestors = label_ancestors(profile=profile, label_id=label_id)
        return bool(set(ancestors) & set(profile.default_excluded_subtree_label_ids))

    labels = [
        label
        for label in profile.labels
        if root in label_ancestors(profile=profile, label_id=label.label_id)
        and label.default_analysis
        and label.label_id not in excluded
        and not in_excluded_subtree(label.label_id)
    ]
    by_id = {label.label_id: label for label in profile.labels}
    comparisons: list[ComparisonDefinition] = []
    for label in labels:
        label_background = _resolved_default_background(
            profile=profile,
            label_id=label.label_id,
        )
        comparisons.append(
            ComparisonDefinition(
                comparison_id=(
                    f"{label.label_id.replace(':', '__').replace('/', '_')}_vs_"
                    f"{label_background.replace(':', '__').replace('/', '_')}"
                ),
                display_name=(
                    f"{label.display_name} versus {by_id[label_background].display_name}"
                ),
                target_label_ids=(label.label_id,),
                background_label_ids=(label_background,),
                description=(
                    f"{profile.default_comparison_description} "
                    "The target stratum is analysed separately; heterogeneous "
                    "multiprotein-system headings are not automatic targets."
                ),
            )
        )
    return tuple(comparisons)


def _resolved_default_background(*, profile: ProteinProfile, label_id: str) -> str:
    """Resolve the nearest label-specific background, then the profile fallback.

    Args:
        profile: Selected classification profile.
        label_id: Default target label.

    Returns:
        Canonical background label identifier.
    """

    by_id = {label.label_id: label for label in profile.labels}
    for candidate_id in label_ancestors(profile=profile, label_id=label_id):
        candidate = by_id[candidate_id].default_background_label_id
        if candidate:
            return candidate
    return profile.default_background_label_id


def _parse_label(*, value: Any, index: int) -> ProfileLabel:
    """Parse one profile-label mapping.

    Args:
        value: Raw YAML value.
        index: Zero-based label index for diagnostics.

    Returns:
        Validated profile label.
    """

    row = require_mapping(value=value, field_name=f"profile.labels[{index}]")
    reject_unknown_fields(
        value=row,
        allowed=frozenset(
            {
                "label_id",
                "display_name",
                "parent_label_id",
                "level",
                "system_class",
                "mechanistic_class",
                "component_role",
                "family",
                "active_site_expected",
                "active_site_residue",
                "default_analysis",
                "default_background_label_id",
                "assignment_exclusivity_group",
                "reviewed_positive_allowed",
                "aliases",
                "description",
            }
        ),
        field_name=f"profile.labels[{index}]",
    )
    aliases = tuple(
        validate_text(value=alias, field_name=f"profile.labels[{index}].aliases")
        for alias in require_sequence(
            value=row.get("aliases", []), field_name=f"profile.labels[{index}].aliases"
        )
    )
    default_analysis = row.get("default_analysis", False)
    if not isinstance(default_analysis, bool):
        raise InputValidationError(
            f"profile.labels[{index}].default_analysis must be true or false."
        )
    reviewed_positive_allowed = row.get("reviewed_positive_allowed", True)
    if not isinstance(reviewed_positive_allowed, bool):
        raise InputValidationError(
            f"profile.labels[{index}].reviewed_positive_allowed must be true or false."
        )
    return ProfileLabel(
        label_id=validate_identifier(
            value=row.get("label_id"), field_name=f"profile.labels[{index}].label_id"
        ),
        display_name=validate_text(
            value=row.get("display_name"), field_name=f"profile.labels[{index}].display_name"
        ),
        parent_label_id=validate_text(
            value=row.get("parent_label_id", ""),
            field_name=f"profile.labels[{index}].parent_label_id",
            allow_empty=True,
        ),
        level=validate_identifier(
            value=row.get("level"), field_name=f"profile.labels[{index}].level"
        ),
        system_class=validate_text(
            value=row.get("system_class", ""),
            field_name=f"profile.labels[{index}].system_class",
            allow_empty=True,
        ),
        mechanistic_class=validate_text(
            value=row.get("mechanistic_class", ""),
            field_name=f"profile.labels[{index}].mechanistic_class",
            allow_empty=True,
        ),
        component_role=validate_text(
            value=row.get("component_role", ""),
            field_name=f"profile.labels[{index}].component_role",
            allow_empty=True,
        ),
        family=validate_text(
            value=row.get("family", ""),
            field_name=f"profile.labels[{index}].family",
            allow_empty=True,
        ),
        active_site_expected=validate_text(
            value=row.get("active_site_expected", "UNKNOWN"),
            field_name=f"profile.labels[{index}].active_site_expected",
        ).upper(),
        active_site_residue=validate_text(
            value=row.get("active_site_residue", ""),
            field_name=f"profile.labels[{index}].active_site_residue",
            allow_empty=True,
        ),
        default_analysis=default_analysis,
        default_background_label_id=validate_text(
            value=row.get("default_background_label_id", ""),
            field_name=f"profile.labels[{index}].default_background_label_id",
            allow_empty=True,
        ),
        assignment_exclusivity_group=validate_text(
            value=row.get("assignment_exclusivity_group", ""),
            field_name=f"profile.labels[{index}].assignment_exclusivity_group",
            allow_empty=True,
        ),
        reviewed_positive_allowed=reviewed_positive_allowed,
        aliases=aliases,
        description=validate_text(
            value=row.get("description"),
            field_name=f"profile.labels[{index}].description",
        ),
    )


def _validate_hierarchy(*, profile: ProteinProfile) -> None:
    """Validate uniqueness, parent links, cycles and aliases.

    Args:
        profile: Candidate profile.

    Raises:
        InputValidationError: If hierarchy invariants fail.
    """

    if not profile.labels:
        raise InputValidationError("A protein profile must contain at least one label.")
    by_id: dict[str, ProfileLabel] = {}
    aliases: dict[str, str] = {}
    for label in profile.labels:
        if label.default_analysis and not label.reviewed_positive_allowed:
            raise InputValidationError(
                f"Default-analysis label {label.label_id!r} must permit "
                "REVIEWED_POSITIVE assignments."
            )
        if label.label_id in by_id:
            raise InputValidationError(f"Duplicate profile label: {label.label_id!r}")
        by_id[label.label_id] = label
        for term in (label.label_id, label.display_name, *label.aliases):
            key = term.casefold()
            previous = aliases.get(key)
            if previous is not None and previous != label.label_id:
                raise InputValidationError(
                    f"Profile alias {term!r} is shared by {previous!r} and {label.label_id!r}."
                )
            aliases[key] = label.label_id
        if label.default_background_label_id == label.label_id:
            raise InputValidationError(
                f"Label {label.label_id!r} cannot use itself as its default background."
            )
    roots = [label for label in profile.labels if not label.parent_label_id]
    if len(roots) != 1:
        raise InputValidationError(
            f"A profile must contain exactly one root label; observed {len(roots)}."
        )
    for label in profile.labels:
        if label.parent_label_id and label.parent_label_id not in by_id:
            raise InputValidationError(
                f"Label {label.label_id!r} has unknown parent {label.parent_label_id!r}."
            )
        visited: set[str] = set()
        current = label
        while current.parent_label_id:
            if current.label_id in visited:
                raise InputValidationError(
                    f"Profile hierarchy contains a cycle at {current.label_id!r}."
                )
            visited.add(current.label_id)
            current = by_id[current.parent_label_id]
    policy_labels = {
        profile.default_target_root_label_id,
        profile.default_background_label_id,
        *profile.default_excluded_label_ids,
        *profile.default_excluded_subtree_label_ids,
        *(
            label.default_background_label_id
            for label in profile.labels
            if label.default_background_label_id
        ),
    } - {""}
    unknown_policy_labels = policy_labels - set(by_id)
    if unknown_policy_labels:
        raise InputValidationError(
            "Profile default-comparison policy references unknown labels: "
            f"{sorted(unknown_policy_labels)}"
        )
    for label in profile.labels:
        if not label.default_background_label_id:
            continue
        background_ancestors = label_ancestors(
            profile=profile,
            label_id=label.default_background_label_id,
        )
        if label.label_id in background_ancestors:
            raise InputValidationError(
                f"Label {label.label_id!r} cannot use one of its descendants as its "
                "default background."
            )
    for label in profile.labels:
        if not label.default_analysis:
            continue
        resolved_background = _resolved_default_background(
            profile=profile,
            label_id=label.label_id,
        )
        if not resolved_background:
            raise InputValidationError(
                f"Default-analysis label {label.label_id!r} has no background policy."
            )
        target_ancestors = label_ancestors(profile=profile, label_id=label.label_id)
        background_ancestors = label_ancestors(
            profile=profile,
            label_id=resolved_background,
        )
        if resolved_background in target_ancestors or label.label_id in background_ancestors:
            raise InputValidationError(
                f"Default-analysis label {label.label_id!r} and resolved background "
                f"{resolved_background!r} overlap in the profile hierarchy."
            )
