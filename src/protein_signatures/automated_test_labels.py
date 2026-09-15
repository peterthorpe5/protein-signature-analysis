"""Create deterministic synthetic labels for end-to-end software smoke tests."""

from __future__ import annotations

import hashlib
import logging
from collections import Counter
from pathlib import Path
from typing import Any

from . import __version__
from .checksums import sha256_file
from .errors import InputValidationError, PublicationError
from .fasta import read_protein_fasta
from .io_utils import iter_tsv, write_json_atomic, write_tsv_atomic
from .models import (
    ComparisonDefinition,
    GroupMembership,
    PartitionAssignment,
    ProteinProfile,
    SequenceRecord,
)
from .orthofinder import discover_orthofinder_layout, read_group_memberships
from .orthofinder_resource import (
    discover_orthofinder_resource,
    read_resource_memberships,
)
from .partitions import assign_partitions
from .profiles import (
    default_profile_comparisons,
    expand_positive_memberships,
    label_ancestors,
    load_profile,
    resolve_label,
    validate_assignment_profile_compatibility,
)
from .redundancy import read_redundancy_clusters
from .tables import LABEL_FIELDS, read_label_assignments
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)
AUTOMATED_TEST_TOKEN = "AUTOMATED_TEST_ONLY"
AUTOMATED_TEST_EVIDENCE_STATUS = "SYNTHETIC_TEST_ONLY"
AUTOMATED_TEST_APPROVER = "AUTOMATED_TEST_MODE"
_MINIMUM_TOTAL_SAMPLES_PER_CLASS = 20
_MINIMUM_DISCOVERY_SAMPLES_PER_CLASS = 10
_MINIMUM_DISCOVERY_GROUPS_PER_CLASS = 3
_MINIMUM_VALIDATION_GROUPS_PER_CLASS = 1


def create_automated_test_labels(
    *,
    sequences_fasta: Path,
    output_labels: Path,
    marker_path: Path,
    profile: str | Path,
    target_label: str = "ALL",
    template_labels: Path | None = None,
    structures: Path | None = None,
    orthofinder_resource: Path | None = None,
    orthofinder_results: Path | None = None,
    orthofinder_group_type: str = "HOG",
    orthofinder_hierarchy_node: str = "N0",
    orthofinder_run_id: str = "automated_smoke_test",
    redundancy_clusters: Path | None = None,
    samples_per_class: int = 20,
    random_seed: int = 1729,
    validation_fraction: float = 0.2,
) -> Path:
    """Generate isolated synthetic target and control memberships for a smoke test.

    The generated labels test software execution only. They are selected without
    biological evidence and must never be interpreted as protein classification.
    Whole leakage-safe partition units are kept exclusive to one synthetic class.

    Args:
        sequences_fasta: Authoritative protein FASTA for any supported campaign.
        output_labels: New test-only label TSV, or a staged byte-identical template.
        marker_path: New JSON audit marker for the synthetic assignment operation.
        profile: Built-in or custom protein-type profile.
        target_label: Default-analysis label or alias, or ``ALL`` to exercise every
            profile-defined default comparison.
        template_labels: Optional one-row-per-protein label template to preserve for
            unselected proteins.
        structures: Optional structure inventory. When supplied, both cohorts are
            selected only from coordinate-analysis-eligible proteins.
        orthofinder_resource: Optional published ``orthofinder-results`` resource.
        orthofinder_results: Optional raw completed OrthoFinder 2.5.5 or 3 output.
        orthofinder_group_type: ``HOG`` or ``LEGACY_ORTHOGROUP``.
        orthofinder_hierarchy_node: HOG hierarchy node, normally ``N0``.
        orthofinder_run_id: Stable identifier for raw OrthoFinder memberships.
        redundancy_clusters: Optional near-redundancy membership TSV.
        samples_per_class: Total synthetic proteins selected for each class.
        random_seed: Non-negative deterministic selection and partition seed.
        validation_fraction: Partition fraction used by the downstream campaign.

    Returns:
        Absolute path to the published test-only audit marker.

    Raises:
        InputValidationError: If inputs cannot form independent test cohorts.
        PublicationError: If a human-owned or existing output would be overwritten.
    """

    destination = Path(output_labels).expanduser().resolve()
    marker = Path(marker_path).expanduser().resolve()
    _require_test_only_destination(path=destination, field_name="output_labels")
    _require_test_only_destination(path=marker, field_name="marker_path")
    if destination == marker:
        raise InputValidationError("Test labels and their audit marker must be different files.")
    if marker.exists():
        raise PublicationError(f"Automated test-label marker already exists: {marker}")
    if (
        isinstance(samples_per_class, bool)
        or not isinstance(samples_per_class, int)
        or samples_per_class < _MINIMUM_TOTAL_SAMPLES_PER_CLASS
    ):
        raise InputValidationError(
            f"samples_per_class must be an integer of at least {_MINIMUM_TOTAL_SAMPLES_PER_CLASS}."
        )
    if isinstance(random_seed, bool) or not isinstance(random_seed, int) or random_seed < 0:
        raise InputValidationError("random_seed must be a non-negative integer.")
    if not 0.0 < validation_fraction < 1.0:
        raise InputValidationError("validation_fraction must be greater than 0 and less than 1.")
    if orthofinder_resource is not None and orthofinder_results is not None:
        raise InputValidationError(
            "Supply either orthofinder_resource or orthofinder_results, not both."
        )

    sequences = read_protein_fasta(path=sequences_fasta)
    protein_ids = frozenset(record.protein_id for record in sequences)
    loaded_profile = load_profile(source=profile)
    comparisons = _resolve_test_comparisons(
        profile=loaded_profile,
        target_label=target_label,
    )
    labels_by_id = {label.label_id: label for label in loaded_profile.labels}
    direct_target_ids = _direct_test_target_labels(
        profile=loaded_profile,
        comparisons=comparisons,
    )
    background_ids = tuple(
        dict.fromkeys(comparison.background_label_ids[0] for comparison in comparisons)
    )
    overlap = frozenset(direct_target_ids) & frozenset(background_ids)
    if overlap:
        raise InputValidationError(
            f"Automated test target and background labels must be disjoint: {sorted(overlap)}"
        )
    cohort_roles = {
        **{label_id: "TARGET" for label_id in direct_target_ids},
        **{label_id: "BACKGROUND" for label_id in background_ids},
    }

    template_path = (
        Path(template_labels).expanduser().resolve() if template_labels is not None else None
    )
    base_records = _base_label_records(
        sequences=sequences,
        profile=loaded_profile,
        template_path=template_path,
    )
    if template_path is not None and destination == template_path:
        raise PublicationError("Automated test labels must not overwrite their source template.")
    if destination.exists():
        if template_path is None or sha256_file(path=destination) != sha256_file(
            path=template_path
        ):
            raise PublicationError(
                "Refusing to overwrite labels that are not the byte-identical staged template: "
                f"{destination}"
            )

    eligible_ids = (
        _eligible_structure_proteins(path=structures, protein_ids=protein_ids)
        if structures is not None
        else protein_ids
    )
    memberships = _load_group_memberships(
        protein_ids=protein_ids,
        orthofinder_resource=orthofinder_resource,
        orthofinder_results=orthofinder_results,
        group_type=orthofinder_group_type,
        hierarchy_node=orthofinder_hierarchy_node,
        run_id=orthofinder_run_id,
    )
    supplied_redundancy = (
        read_redundancy_clusters(path=redundancy_clusters, protein_ids=protein_ids)
        if redundancy_clusters is not None
        else ()
    )
    partitions = assign_partitions(
        sequences=sequences,
        memberships=memberships,
        validation_fraction=validation_fraction,
        random_seed=random_seed,
        redundancy_memberships=supplied_redundancy,
    )
    selections = _select_test_cohorts(
        partitions=partitions,
        eligible_ids=eligible_ids,
        cohort_roles=cohort_roles,
        samples_per_class=samples_per_class,
        validation_fraction=validation_fraction,
        random_seed=random_seed,
    )

    selected_by_protein = {
        row["protein_id"]: row for cohort_rows in selections.values() for row in cohort_rows
    }
    records: list[dict[str, str]] = []
    for protein_id in sorted(base_records):
        record = dict(base_records[protein_id])
        selection = selected_by_protein.get(protein_id)
        if selection is not None:
            selected_label_id = selection["label_id"]
            selected_label = labels_by_id[selected_label_id]
            partition_digest = hashlib.sha256(
                f"{selection['partition_unit']}:{selection['partition_key']}".encode("utf-8")
            ).hexdigest()
            record.update(
                {
                    "label_id": selected_label_id,
                    "curation_status": "REVIEWED_POSITIVE",
                    "evidence_status": AUTOMATED_TEST_EVIDENCE_STATUS,
                    "evidence_source": ("protein-signature-analysis automated software smoke test"),
                    "evidence_reference": (
                        f"{AUTOMATED_TEST_TOKEN}:{selection['partition']}:{partition_digest}"
                    ),
                    "component_role": selected_label.component_role or "UNKNOWN",
                    "curation_reason": (
                        "Synthetic deterministic membership for end-to-end software testing; "
                        "not biological evidence and not scientifically interpretable."
                    ),
                }
            )
        records.append(record)

    write_tsv_atomic(path=destination, fieldnames=LABEL_FIELDS, records=records)
    generated = read_label_assignments(
        path=destination,
        protein_ids=protein_ids,
        label_ids=loaded_profile.label_ids(),
    )
    validate_assignment_profile_compatibility(
        assignments=generated,
        profile=loaded_profile,
    )
    positive_counts = Counter(
        assignment.label_id for assignment in generated if assignment.is_eligible_positive
    )
    expected_counts = dict.fromkeys(cohort_roles, samples_per_class)
    if {label_id: positive_counts[label_id] for label_id in expected_counts} != expected_counts:
        raise PublicationError("Published automated test-label cohort counts changed unexpectedly.")
    expanded_counts = Counter(
        row["label_id"]
        for row in expand_positive_memberships(
            assignments=generated,
            profile=loaded_profile,
        )
    )
    comparison_cohort_counts = [
        {
            "comparison_id": comparison.comparison_id,
            "target_label_id": comparison.target_label_ids[0],
            "target_protein_count": expanded_counts[comparison.target_label_ids[0]],
            "background_label_id": comparison.background_label_ids[0],
            "background_protein_count": expanded_counts[comparison.background_label_ids[0]],
        }
        for comparison in comparisons
    ]
    if any(
        row["target_protein_count"] < samples_per_class
        or row["background_protein_count"] < samples_per_class
        for row in comparison_cohort_counts
    ):
        raise PublicationError(
            "Automated direct and inherited memberships do not cover every comparison."
        )

    marker_value: dict[str, Any] = {
        "schema_version": 1,
        "status": "AUTOMATED_TEST_ONLY",
        "action": "AUTOMATED_TEST_LABEL_GENERATION",
        "package_version": __version__,
        "warning": (
            "Synthetic labels for software testing only; biological or scientific "
            "interpretation is prohibited."
        ),
        "scientific_interpretation_allowed": False,
        "profile_id": loaded_profile.profile_id,
        "profile_version": loaded_profile.profile_version,
        "comparison_count": len(comparisons),
        "comparisons": comparison_cohort_counts,
        "direct_target_label_ids": list(direct_target_ids),
        "background_label_ids": list(background_ids),
        "cohort_count": len(cohort_roles),
        "samples_per_class": samples_per_class,
        "random_seed": random_seed,
        "validation_fraction": validation_fraction,
        "structure_eligible_only": structures is not None,
        "sequence_fasta": str(Path(sequences_fasta).expanduser().resolve()),
        "sequence_fasta_sha256": sha256_file(path=Path(sequences_fasta).expanduser().resolve()),
        "template_labels": str(template_path) if template_path is not None else None,
        "template_labels_sha256": (
            sha256_file(path=template_path) if template_path is not None else None
        ),
        "output_labels": str(destination),
        "output_labels_sha256": sha256_file(path=destination),
        "orthofinder_grouping_used": bool(memberships),
        "near_redundancy_used": bool(supplied_redundancy),
        "selected_cohorts": selections,
    }
    write_json_atomic(path=marker, value=marker_value)
    LOGGER.warning(
        "Published AUTOMATED TEST ONLY labels comparisons=%d cohorts=%d "
        "samples_per_class=%d at %s; these labels are not scientific evidence",
        len(comparisons),
        len(cohort_roles),
        samples_per_class,
        destination,
    )
    return marker


def _resolve_test_comparisons(
    *, profile: ProteinProfile, target_label: str
) -> tuple[ComparisonDefinition, ...]:
    """Resolve all defaults or one single-target, single-background comparison."""

    comparisons = default_profile_comparisons(profile=profile)
    if not comparisons:
        raise InputValidationError("The selected profile has no default testable comparisons.")
    if target_label.strip().upper() == "ALL":
        invalid = tuple(
            comparison
            for comparison in comparisons
            if len(comparison.target_label_ids) != 1 or len(comparison.background_label_ids) != 1
        )
        if invalid:
            raise InputValidationError(
                "Automated all-class testing requires single-target, single-background "
                "profile-default comparisons."
            )
        return comparisons
    target_id = resolve_label(profile=profile, term=target_label)
    matches = tuple(
        comparison
        for comparison in comparisons
        if comparison.target_label_ids == (target_id,) and len(comparison.background_label_ids) == 1
    )
    if len(matches) != 1:
        available = sorted(item.target_label_ids[0] for item in comparisons)
        raise InputValidationError(
            f"Automated test target {target_id!r} is not one uniquely testable profile "
            f"default; choose one of {available}."
        )
    return matches


def _direct_test_target_labels(
    *, profile: ProteinProfile, comparisons: tuple[ComparisonDefinition, ...]
) -> tuple[str, ...]:
    """Return terminal tested targets whose membership covers tested ancestors."""

    targets = tuple(comparison.target_label_ids[0] for comparison in comparisons)
    target_set = frozenset(targets)
    terminal = tuple(
        target
        for target in targets
        if not any(
            target != candidate
            and target in label_ancestors(profile=profile, label_id=candidate)[1:]
            for candidate in target_set
        )
    )
    if not terminal:
        raise InputValidationError("Automated test comparisons have no terminal target labels.")
    return terminal


def _base_label_records(
    *,
    sequences: tuple[SequenceRecord, ...],
    profile: ProteinProfile,
    template_path: Path | None,
) -> dict[str, dict[str, str]]:
    """Load a one-row template or create generic unmapped rows for every protein."""

    protein_ids = frozenset(record.protein_id for record in sequences)
    if template_path is None:
        return {
            protein_id: {
                "protein_id": protein_id,
                "label_id": profile.default_target_root_label_id,
                "curation_status": "UNMAPPED",
                "evidence_status": "AUTOMATED_TEST_NOT_SELECTED",
                "evidence_source": "protein-signature-analysis automated software smoke test",
                "evidence_reference": AUTOMATED_TEST_TOKEN,
                "component_role": "UNKNOWN",
                "curation_reason": (
                    "Not selected for the synthetic software-test target/control cohorts."
                ),
            }
            for protein_id in sorted(protein_ids)
        }
    assignments = read_label_assignments(
        path=template_path,
        protein_ids=protein_ids,
        label_ids=profile.label_ids(),
    )
    counts = Counter(assignment.protein_id for assignment in assignments)
    if frozenset(counts) != protein_ids or any(count != 1 for count in counts.values()):
        raise InputValidationError(
            "Automated test-label templates require exactly one row for every FASTA protein."
        )
    return {assignment.protein_id: assignment.to_record() for assignment in assignments}


def _eligible_structure_proteins(*, path: Path, protein_ids: frozenset[str]) -> frozenset[str]:
    """Return proteins explicitly eligible for coordinate-based analysis."""

    eligible: set[str] = set()
    for row in iter_tsv(
        path=path,
        required_fields=(
            "protein_id",
            "availability_status",
            "coordinate_path",
            "analysis_eligibility_status",
        ),
    ):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        if protein_id not in protein_ids:
            raise InputValidationError(
                f"Structure inventory references protein absent from FASTA: {protein_id!r}"
            )
        if (
            row["analysis_eligibility_status"].strip().upper() == "ELIGIBLE"
            and row["availability_status"].strip().upper() in {"AVAILABLE", "COMPLETE"}
            and row["coordinate_path"].strip()
        ):
            eligible.add(protein_id)
    if not eligible:
        raise InputValidationError(
            "Structure-restricted automated test labels found no eligible coordinate models."
        )
    return frozenset(eligible)


def _load_group_memberships(
    *,
    protein_ids: frozenset[str],
    orthofinder_resource: Path | None,
    orthofinder_results: Path | None,
    group_type: str,
    hierarchy_node: str,
    run_id: str,
) -> tuple[GroupMembership, ...]:
    """Load optional OrthoFinder grouping from either supported authority."""

    authority = group_type.strip().upper()
    node = hierarchy_node if authority == "HOG" else ""
    if orthofinder_resource is not None:
        resource = discover_orthofinder_resource(resource_dir=orthofinder_resource)
        return read_resource_memberships(
            resource=resource,
            group_type=authority,
            hierarchy_node=node,
            protein_ids=protein_ids,
        )
    if orthofinder_results is not None:
        layout = discover_orthofinder_layout(results_dir=orthofinder_results)
        return read_group_memberships(
            layout=layout,
            run_id=validate_identifier(value=run_id, field_name="orthofinder_run_id"),
            group_type=authority,
            hierarchy_node=node,
            protein_ids=protein_ids,
        )
    return ()


def _select_test_cohorts(
    *,
    partitions: tuple[PartitionAssignment, ...],
    eligible_ids: frozenset[str],
    cohort_roles: dict[str, str],
    samples_per_class: int,
    validation_fraction: float,
    random_seed: int,
) -> dict[str, list[dict[str, str]]]:
    """Select disjoint label cohorts without splitting homology units.

    Every eligible partition unit is assigned to at most one label. Discovery
    cohorts span at least three independent units, while held-out validation
    cohorts span at least one. The latter deliberately permits small smoke-test
    validation cohorts: the route exercises training and SHAP generation, but
    does not pretend that synthetic held-out performance is scientific evidence.
    """

    if not cohort_roles:
        raise InputValidationError("At least one automated test cohort is required.")
    invalid_roles = {
        label_id: role
        for label_id, role in cohort_roles.items()
        if role not in {"TARGET", "BACKGROUND"}
    }
    if invalid_roles:
        raise InputValidationError(f"Automated test cohorts contain invalid roles: {invalid_roles}")
    validation_count = max(1, round(samples_per_class * validation_fraction))
    discovery_count = samples_per_class - validation_count
    if discovery_count < _MINIMUM_DISCOVERY_SAMPLES_PER_CLASS:
        raise InputValidationError(
            "samples_per_class and validation_fraction must leave at least "
            f"{_MINIMUM_DISCOVERY_SAMPLES_PER_CLASS} synthetic discovery samples "
            "per class."
        )

    proteins_by_unit: dict[tuple[str, str, str], list[str]] = {}
    for assignment in partitions:
        if assignment.protein_id not in eligible_ids:
            continue
        if assignment.partition not in {"DISCOVERY", "VALIDATION"}:
            raise InputValidationError(
                "Automated test partition assignments must be DISCOVERY or VALIDATION; "
                f"received {assignment.partition!r}."
            )
        unit = (
            assignment.partition,
            assignment.partition_unit,
            assignment.partition_key,
        )
        proteins_by_unit.setdefault(unit, []).append(assignment.protein_id)

    units_by_partition: dict[str, list[dict[str, Any]]] = {
        "DISCOVERY": [],
        "VALIDATION": [],
    }
    for (partition, unit_type, unit_key), protein_ids in proteins_by_unit.items():
        ordered_proteins = sorted(
            protein_ids,
            key=lambda protein_id: (
                _selection_digest(
                    random_seed=random_seed,
                    values=(partition, unit_type, unit_key, protein_id),
                ),
                protein_id,
            ),
        )
        units_by_partition[partition].append(
            {
                "partition_unit": unit_type,
                "partition_key": unit_key,
                "protein_ids": ordered_proteins,
            }
        )
    for partition, units in units_by_partition.items():
        units.sort(
            key=lambda unit: (
                -len(unit["protein_ids"]),
                _selection_digest(
                    random_seed=random_seed,
                    values=(
                        partition,
                        str(unit["partition_unit"]),
                        str(unit["partition_key"]),
                    ),
                ),
            )
        )

    requirements = {
        "DISCOVERY": (
            discovery_count,
            _MINIMUM_DISCOVERY_GROUPS_PER_CLASS,
        ),
        "VALIDATION": (
            validation_count,
            _MINIMUM_VALIDATION_GROUPS_PER_CLASS,
        ),
    }
    allocated: dict[str, dict[str, list[dict[str, Any]]]] = {
        label_id: {"DISCOVERY": [], "VALIDATION": []} for label_id in cohort_roles
    }
    cohort_order = sorted(
        cohort_roles,
        key=lambda label_id: (
            _selection_digest(random_seed=random_seed, values=("COHORT", label_id)),
            label_id,
        ),
    )
    for partition, (sample_count, minimum_group_count) in requirements.items():
        _allocate_partition_units(
            partition=partition,
            available_units=units_by_partition[partition],
            allocated=allocated,
            cohort_order=cohort_order,
            sample_count=sample_count,
            minimum_group_count=minimum_group_count,
        )

    cohorts: dict[str, list[dict[str, str]]] = {label_id: [] for label_id in cohort_roles}
    for label_id in cohort_order:
        for partition, (sample_count, _) in requirements.items():
            selected = _sample_allocated_units(
                units=allocated[label_id][partition],
                sample_count=sample_count,
            )
            cohorts[label_id].extend(
                {
                    "protein_id": protein_id,
                    "label_id": label_id,
                    "cohort_role": cohort_roles[label_id],
                    "partition": partition,
                    "partition_unit": str(unit["partition_unit"]),
                    "partition_key": str(unit["partition_key"]),
                }
                for unit, protein_id in selected
            )
    for values in cohorts.values():
        values.sort(
            key=lambda row: (
                row["partition"],
                row["partition_unit"],
                row["partition_key"],
                row["protein_id"],
            )
        )
    return cohorts


def _allocate_partition_units(
    *,
    partition: str,
    available_units: list[dict[str, Any]],
    allocated: dict[str, dict[str, list[dict[str, Any]]]],
    cohort_order: list[str],
    sample_count: int,
    minimum_group_count: int,
) -> None:
    """Distribute whole units fairly until every label meets one requirement."""

    capacities = {label_id: 0 for label_id in cohort_order}
    cohort_rank = {label_id: index for index, label_id in enumerate(cohort_order)}
    for unit in available_units:
        incomplete = [
            label_id
            for label_id in cohort_order
            if len(allocated[label_id][partition]) < minimum_group_count
            or capacities[label_id] < sample_count
        ]
        if not incomplete:
            break
        label_id = min(
            incomplete,
            key=lambda item: (
                len(allocated[item][partition]) >= minimum_group_count,
                capacities[item] / sample_count,
                len(allocated[item][partition]) / minimum_group_count,
                cohort_rank[item],
            ),
        )
        allocated[label_id][partition].append(unit)
        capacities[label_id] += len(unit["protein_ids"])

    unsatisfied = {
        label_id: {
            "required_samples": sample_count,
            "available_samples": capacities[label_id],
            "required_independent_units": minimum_group_count,
            "available_independent_units": len(allocated[label_id][partition]),
        }
        for label_id in cohort_order
        if capacities[label_id] < sample_count
        or len(allocated[label_id][partition]) < minimum_group_count
    }
    if unsatisfied:
        raise InputValidationError(
            "Insufficient disjoint eligible partition units for automated test cohorts "
            f"in {partition}: eligible_units={len(available_units)}, "
            f"cohorts={len(cohort_order)}, unsatisfied={unsatisfied}."
        )


def _sample_allocated_units(
    *,
    units: list[dict[str, Any]],
    sample_count: int,
) -> list[tuple[dict[str, Any], str]]:
    """Select an exact protein count while retaining every allocated unit."""

    selected = [(unit, str(unit["protein_ids"][0])) for unit in units]
    remaining = [
        (unit, str(protein_id)) for unit in units for protein_id in unit["protein_ids"][1:]
    ]
    selected.extend(remaining[: sample_count - len(selected)])
    if len(selected) != sample_count:
        raise PublicationError(
            "Automated cohort allocation did not retain its requested sample count."
        )
    return selected


def _selection_digest(*, random_seed: int, values: tuple[str, ...]) -> str:
    """Return a stable seeded ordering digest for one selection tuple."""

    payload = ":".join((str(random_seed), *values))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _require_test_only_destination(*, path: Path, field_name: str) -> None:
    """Require conspicuous test-only naming for every synthetic output authority."""

    if AUTOMATED_TEST_TOKEN not in path.name.upper():
        raise InputValidationError(
            f"{field_name} filename must contain {AUTOMATED_TEST_TOKEN!r}: {path}"
        )
