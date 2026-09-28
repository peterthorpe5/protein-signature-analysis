"""Apply frozen, validated class signatures to new proteins or groups."""

from __future__ import annotations

import logging
import math
import os
import shutil
import tempfile
from collections import defaultdict
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .checksums import sha256_file
from .errors import InputValidationError, PublicationError
from .fasta import read_protein_fasta, sequence_kmers
from .feature_provenance import FEATURE_ASSESSMENT_STATUSES
from .io_utils import iter_tsv, write_json_atomic, write_tsv_atomic
from .models import SequenceRecord
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)

SUMMARY_FIELDS = (
    "unit_id",
    "unit_type",
    "member_count",
    "comparison_id",
    "target_label_ids",
    "background_label_ids",
    "reference_feature_count",
    "present_feature_count",
    "assessed_absent_feature_count",
    "unknown_feature_count",
    "present_feature_types",
    "evidence_state",
)
HIT_FIELDS = (
    "unit_id",
    "comparison_id",
    "feature_type",
    "feature_id",
    "feature_name",
    "positive_protein_ids",
    "validation_q_value",
    "validation_fdr_scope",
    "validation_prevalence_difference",
    "evidence_class",
)
STRUCTURAL_PROJECTION_FIELDS = (
    "protein_id",
    "feature_type",
    "feature_id",
    "evidence_status",
    "supporting_reference_count",
    "best_tm_score",
    "comparison_tool",
    "comparison_tool_version",
    "coverage_scope",
    "query_comparison_universe_id",
)


@dataclass(frozen=True)
class _FrozenFeature:
    """One externally validated reference signature feature."""

    feature_type: str
    feature_id: str
    feature_name: str
    q_value: float
    fdr_scope: str
    prevalence_difference: float
    evidence_class: str


@dataclass(frozen=True)
class _ClusterReference:
    """One discovery-frozen structural cluster definition."""

    cluster_id: str
    comparison_tool: str
    comparison_tool_version: str
    coverage_scope: str
    tm_score_threshold: float
    minimum_coverage: float


def project_structural_signatures(
    *,
    reference_clusters_path: Path,
    query_comparisons_path: Path,
    query_fasta: Path,
    output_dir: Path,
) -> Path:
    """Project new structural hits onto frozen discovery cluster members.

    Only passing positive evidence is emitted. A missing hit stays unknown,
    because a retained-hit file cannot prove an exhaustive negative search.
    The output's first four columns can be used by ``query-signatures``.

    Args:
        reference_clusters_path: Canonical discovery cluster membership TSV.
        query_comparisons_path: Query-to-reference pairwise alignment TSV.
        query_fasta: Exact query protein identifier authority.
        output_dir: New directory for the immutable projection result.

    Returns:
        Directory with positive feature rows and checksum-bound provenance.

    Raises:
        InputValidationError: If a pairwise record or frozen definition conflicts.
        PublicationError: If the output directory exists or publication fails.
    """

    query_ids = frozenset(item.protein_id for item in read_protein_fasta(path=query_fasta))
    references = _load_cluster_references(path=reference_clusters_path)
    overlapping = query_ids & references.keys()
    if overlapping:
        raise InputValidationError(
            f"Query and reference identifiers overlap: {sorted(overlapping)[:10]}"
        )
    support: dict[tuple[str, str], tuple[set[str], float, _ClusterReference, str]] = {}
    required = (
        "protein_a_id",
        "protein_b_id",
        "comparison_tool",
        "comparison_tool_version",
        "tm_score",
        "coverage_a",
        "coverage_b",
        "comparison_status",
        "comparison_universe_id",
        "coverage_scope",
    )
    for row in iter_tsv(path=query_comparisons_path, required_fields=required, allow_empty=True):
        a_query = row["protein_a_id"] in query_ids
        b_query = row["protein_b_id"] in query_ids
        if a_query == b_query:
            raise InputValidationError("Each query alignment must have exactly one query endpoint.")
        query_id = row["protein_a_id"] if a_query else row["protein_b_id"]
        reference_id = row["protein_b_id"] if a_query else row["protein_a_id"]
        if reference_id not in references:
            continue
        if row["comparison_status"] not in {"COMPLETE", "SUCCESS", "PASS"}:
            continue
        try:
            tm_score = float(row["tm_score"])
            query_coverage = float(row["coverage_a"] if a_query else row["coverage_b"])
            reference_coverage = float(row["coverage_b"] if a_query else row["coverage_a"])
        except ValueError as error:
            raise InputValidationError(
                "Complete query alignments require numeric scores."
            ) from error
        if any(
            not math.isfinite(value) or not 0 <= value <= 1
            for value in (tm_score, query_coverage, reference_coverage)
        ):
            raise InputValidationError("Structural scores and coverages must lie between 0 and 1.")
        universe_id = validate_identifier(
            value=row["comparison_universe_id"], field_name="comparison_universe_id"
        )
        compatible_definitions = tuple(
            definition
            for definition in references[reference_id]
            if (row["comparison_tool"], row["comparison_tool_version"], row["coverage_scope"])
            == (
                definition.comparison_tool,
                definition.comparison_tool_version,
                definition.coverage_scope,
            )
        )
        if not compatible_definitions:
            raise InputValidationError(
                f"Query comparison is incompatible with frozen reference {reference_id!r}."
            )
        for definition in compatible_definitions:
            if (
                tm_score < definition.tm_score_threshold
                or query_coverage < definition.minimum_coverage
                or reference_coverage < definition.minimum_coverage
            ):
                continue
            key = (query_id, definition.cluster_id)
            previous = support.get(key)
            if previous is None:
                support[key] = ({reference_id}, tm_score, definition, universe_id)
            else:
                references_seen, previous_score, previous_definition, previous_universe = previous
                if universe_id != previous_universe or definition != previous_definition:
                    raise InputValidationError(
                        f"Conflicting structural projection provenance for {key!r}."
                    )
                references_seen.add(reference_id)
                support[key] = (
                    references_seen,
                    max(previous_score, tm_score),
                    definition,
                    universe_id,
                )
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise PublicationError(f"Projection output directory already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        table = staging / "projected_structural_features.tsv"
        count = write_tsv_atomic(
            path=table,
            fieldnames=STRUCTURAL_PROJECTION_FIELDS,
            records=(
                {
                    "protein_id": query_id,
                    "feature_type": "STRUCTURE_CLUSTER",
                    "feature_id": cluster_id,
                    "evidence_status": "ASSESSED_WITH_FEATURE",
                    "supporting_reference_count": len(members),
                    "best_tm_score": score,
                    "comparison_tool": definition.comparison_tool,
                    "comparison_tool_version": definition.comparison_tool_version,
                    "coverage_scope": definition.coverage_scope,
                    "query_comparison_universe_id": universe_id,
                }
                for (query_id, cluster_id), (members, score, definition, universe_id) in sorted(
                    support.items()
                )
            ),
        )
        write_json_atomic(
            path=staging / "PROJECTION_COMPLETED.json",
            value={
                "schema_version": 1,
                "status": "POSITIVE_STRUCTURAL_PROJECTIONS_COMPLETE",
                "package_version": __version__,
                "interpretation": "NO_HIT_IS_UNKNOWN_NOT_ASSESSED_NEGATIVE",
                "positive_projection_count": count,
                "reference_clusters_sha256": sha256_file(path=reference_clusters_path),
                "query_comparisons_sha256": sha256_file(path=query_comparisons_path),
                "query_fasta_sha256": sha256_file(path=query_fasta),
                "output_sha256": sha256_file(path=table),
            },
        )
        os.replace(staging, destination)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    LOGGER.info("Projected %d positive structural cluster memberships at %s", count, destination)
    return destination


def _load_cluster_references(*, path: Path) -> dict[str, tuple[_ClusterReference, ...]]:
    """Load only the original discovery members of frozen structural clusters.

    Args:
        path: Canonical structure_clusters TSV or TSV.GZ.

    Returns:
        Reference protein identifiers mapped to compatible cluster definitions.
    """

    required = (
        "cluster_id",
        "protein_id",
        "reference_partition",
        "membership_method",
        "tm_score_threshold",
        "minimum_coverage",
        "coverage_scope",
        "comparison_tool",
        "comparison_tool_version",
    )
    by_protein: dict[str, set[_ClusterReference]] = defaultdict(set)
    by_cluster: dict[str, _ClusterReference] = {}
    for row in iter_tsv(path=path, required_fields=required):
        if row["reference_partition"] != "DISCOVERY":
            continue
        if row["membership_method"] != "DISCOVERY_COMPONENT":
            continue
        try:
            tm_score = float(row["tm_score_threshold"])
            coverage = float(row["minimum_coverage"])
        except ValueError as error:
            raise InputValidationError("Frozen cluster thresholds must be numeric.") from error
        if any(not math.isfinite(value) or not 0 <= value <= 1 for value in (tm_score, coverage)):
            raise InputValidationError("Frozen cluster thresholds must lie between 0 and 1.")
        cluster_id = validate_identifier(value=row["cluster_id"], field_name="cluster_id")
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        definition = _ClusterReference(
            cluster_id=cluster_id,
            comparison_tool=row["comparison_tool"],
            comparison_tool_version=row["comparison_tool_version"],
            coverage_scope=row["coverage_scope"],
            tm_score_threshold=tm_score,
            minimum_coverage=coverage,
        )
        previous = by_cluster.setdefault(cluster_id, definition)
        if previous != definition:
            raise InputValidationError(f"Conflicting frozen structural cluster: {cluster_id!r}")
        if definition in by_protein[protein_id]:
            raise InputValidationError(f"Duplicate reference cluster membership: {protein_id!r}")
        by_protein[protein_id].add(definition)
    if not by_protein:
        raise InputValidationError(
            "No discovery reference members were found in structural clusters."
        )
    return {
        protein_id: tuple(sorted(definitions, key=lambda item: item.cluster_id))
        for protein_id, definitions in by_protein.items()
    }


def query_signature_evidence(
    *,
    signatures_path: Path,
    comparisons_path: Path,
    query_fasta: Path,
    query_features_path: Path | None,
    output_dir: Path,
    query_units_path: Path | None = None,
    additional_query_features_paths: tuple[Path, ...] = (),
    include_within_comparison: bool = False,
    comparison_ids: tuple[str, ...] = (),
) -> Path:
    """Audit frozen signature evidence in a new protein collection.

    This is candidate prioritisation, not a trained classifier. A missing or
    failed feature assessment is never interpreted as a negative result. When
    groups are supplied, presence means at least one member is positive and
    absence requires every member to be explicitly assessed negative.

    Args:
        signatures_path: Immutable canonical signatures TSV or TSV.GZ.
        comparisons_path: Reference comparison definitions TSV or TSV.GZ.
        query_fasta: Authoritative query protein FASTA.
        query_features_path: Optional feature assessments using frozen reference IDs.
            Exact amino-acid k-mers are assessed natively from query FASTA.
        output_dir: New directory for the atomic query result.
        query_units_path: Optional protein-to-unit TSV; may contain overlapping
            orthogroups, but every query protein must appear in at least one.
        additional_query_features_paths: Further feature assessment tables,
            such as positive structural projections alongside sequence calls.
        include_within_comparison: Include within-comparison validated features
            as well as study-wide validated features.
        comparison_ids: Optional exact reference comparison IDs to query. Empty
            selects every comparison from the frozen reference authority.

    Returns:
        Published result directory containing summary, positive hits and marker.

    Raises:
        InputValidationError: If feature states or reference mappings conflict.
        PublicationError: If the destination exists or publishing fails.
    """

    if not isinstance(include_within_comparison, bool):
        raise InputValidationError("include_within_comparison must be Boolean.")
    proteins = read_protein_fasta(path=query_fasta)
    protein_ids = frozenset(item.protein_id for item in proteins)
    all_comparisons = _load_comparisons(path=comparisons_path)
    comparisons = all_comparisons
    if comparison_ids:
        if len(comparison_ids) != len(set(comparison_ids)):
            raise InputValidationError("Query comparison IDs must be unique.")
        unknown = set(comparison_ids) - all_comparisons.keys()
        if unknown:
            raise InputValidationError(f"Unknown query comparison IDs: {sorted(unknown)}")
        comparisons = {
            comparison_id: all_comparisons[comparison_id]
            for comparison_id in sorted(comparison_ids)
        }
    frozen = _load_frozen_signatures(
        path=signatures_path,
        comparison_ids=frozenset(all_comparisons),
        include_within_comparison=include_within_comparison,
    )
    frozen = {
        comparison_id: features
        for comparison_id, features in frozen.items()
        if comparison_id in comparisons
    }
    feature_keys = frozenset(
        (feature.feature_type, feature.feature_id)
        for features in frozen.values()
        for feature in features
    )
    feature_paths = tuple(
        path for path in (query_features_path, *additional_query_features_paths) if path is not None
    )
    native_kmers = _scan_reference_kmers(proteins=proteins, feature_keys=feature_keys)
    states = _load_query_features(
        paths=feature_paths,
        protein_ids=protein_ids,
        feature_keys=feature_keys,
    )
    for (protein_id, feature_type, feature_id), status in states.items():
        if feature_type == "AMINO_ACID_KMER":
            native_status = (
                "ASSESSED_WITH_FEATURE"
                if feature_id in native_kmers[protein_id]
                else "ASSESSED_NO_FEATURE"
            )
            if status != native_status:
                raise InputValidationError(
                    f"Supplied k-mer assessment disagrees with query FASTA: "
                    f"{protein_id!r}, {feature_id!r}."
                )
    units = _load_query_units(path=query_units_path, protein_ids=protein_ids)
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise PublicationError(f"Query output directory already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        summary_count = write_tsv_atomic(
            path=staging / "candidate_summary.tsv",
            fieldnames=SUMMARY_FIELDS,
            records=(
                summary
                for summary, _ in _iter_query_assessments(
                    comparisons=comparisons,
                    frozen=frozen,
                    units=units,
                    states=states,
                    native_kmers=native_kmers,
                    unit_type="GROUP" if query_units_path is not None else "PROTEIN",
                )
            ),
        )
        hit_count = write_tsv_atomic(
            path=staging / "candidate_feature_hits.tsv",
            fieldnames=HIT_FIELDS,
            records=(
                hit
                for _, hits in _iter_query_assessments(
                    comparisons=comparisons,
                    frozen=frozen,
                    units=units,
                    states=states,
                    native_kmers=native_kmers,
                    unit_type="GROUP" if query_units_path is not None else "PROTEIN",
                )
                for hit in hits
            ),
        )
        inputs = {
            "signatures_sha256": sha256_file(path=signatures_path),
            "comparisons_sha256": sha256_file(path=comparisons_path),
            "query_fasta_sha256": sha256_file(path=query_fasta),
            "query_feature_tables": [
                {"path": str(path.expanduser().resolve()), "sha256": sha256_file(path=path)}
                for path in feature_paths
            ],
            "query_units_sha256": (
                sha256_file(path=query_units_path) if query_units_path is not None else None
            ),
        }
        write_json_atomic(
            path=staging / "QUERY_COMPLETED.json",
            value={
                "schema_version": 1,
                "status": "QUERY_EVIDENCE_COMPLETE",
                "package_version": __version__,
                "interpretation": "UNREVIEWED_CANDIDATE_EVIDENCE_NOT_CLASS_PROBABILITY",
                "unit_type": "GROUP" if query_units_path is not None else "PROTEIN",
                "unit_count": len(units),
                "comparison_count": len(comparisons),
                "selected_comparison_ids": sorted(comparisons),
                "summary_row_count": summary_count,
                "positive_hit_row_count": hit_count,
                "include_within_comparison": include_within_comparison,
                "exact_kmers_assessed_from_query_fasta": True,
                "inputs": inputs,
                "outputs": {
                    name: sha256_file(path=staging / name)
                    for name in ("candidate_summary.tsv", "candidate_feature_hits.tsv")
                },
            },
        )
        os.replace(staging, destination)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    LOGGER.info(
        "Published candidate evidence units=%d comparisons=%d hits=%d at %s",
        len(units),
        len(comparisons),
        hit_count,
        destination,
    )
    return destination


def _load_comparisons(*, path: Path) -> dict[str, tuple[str, str]]:
    """Load one target/background mapping for every reference comparison.

    Args:
        path: Canonical comparisons TSV.

    Returns:
        Comparison IDs mapped to target and background label ID strings.
    """

    result: dict[str, tuple[str, str]] = {}
    for row in iter_tsv(
        path=path,
        required_fields=("comparison_id", "target_label_ids", "background_label_ids"),
    ):
        comparison_id = validate_identifier(value=row["comparison_id"], field_name="comparison_id")
        target = row["target_label_ids"].strip()
        background = row["background_label_ids"].strip()
        if not target or not background or comparison_id in result:
            raise InputValidationError(
                f"Missing or duplicate comparison definition: {comparison_id}"
            )
        result[comparison_id] = (target, background)
    return result


def _load_frozen_signatures(
    *,
    path: Path,
    comparison_ids: frozenset[str],
    include_within_comparison: bool,
) -> dict[str, tuple[_FrozenFeature, ...]]:
    """Select positive, direction-consistent validated reference features.

    Args:
        path: Canonical signatures table.
        comparison_ids: IDs declared in the supplied comparisons authority.
        include_within_comparison: Whether to admit local validation evidence.

    Returns:
        Frozen reference features grouped by comparison.
    """

    accepted = {"DECISION_CANDIDATE__VALIDATED_STUDY_WIDE"}
    if include_within_comparison:
        accepted.add("DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON")
    features: dict[str, list[_FrozenFeature]] = defaultdict(list)
    seen: set[tuple[str, str, str]] = set()
    required = (
        "comparison_id",
        "feature_type",
        "feature_id",
        "feature_name",
        "validation_q_value",
        "validation_study_q_value",
        "validation_prevalence_difference",
        "evidence_class",
        "status",
    )
    for row in iter_tsv(path=path, required_fields=required):
        comparison_id = row["comparison_id"]
        if comparison_id not in comparison_ids:
            raise InputValidationError(f"Signature has unknown comparison: {comparison_id!r}")
        if row["evidence_class"] not in accepted or row["status"] != "COMPLETE":
            continue
        try:
            q_value = float(
                row["validation_study_q_value"]
                if row["evidence_class"].endswith("STUDY_WIDE")
                else row["validation_q_value"]
            )
            difference = float(row["validation_prevalence_difference"])
        except ValueError as error:
            raise InputValidationError(
                "Validated signatures require numeric held-out evidence."
            ) from error
        if not math.isfinite(q_value) or not 0 <= q_value <= 1 or not math.isfinite(difference):
            raise InputValidationError("Validated signature statistics are outside valid bounds.")
        if difference <= 0:
            continue
        feature_type = validate_identifier(value=row["feature_type"], field_name="feature_type")
        feature_id = validate_identifier(value=row["feature_id"], field_name="feature_id")
        key = (comparison_id, feature_type, feature_id)
        if key in seen:
            raise InputValidationError(f"Duplicate validated signature: {key!r}")
        seen.add(key)
        features[comparison_id].append(
            _FrozenFeature(
                feature_type=feature_type,
                feature_id=feature_id,
                feature_name=row["feature_name"],
                q_value=q_value,
                fdr_scope=(
                    "STUDY_WIDE"
                    if row["evidence_class"].endswith("STUDY_WIDE")
                    else "WITHIN_COMPARISON"
                ),
                prevalence_difference=difference,
                evidence_class=row["evidence_class"],
            )
        )
    return {
        comparison_id: tuple(sorted(rows, key=lambda row: (row.feature_type, row.feature_id)))
        for comparison_id, rows in features.items()
    }


def _load_query_features(
    *,
    paths: tuple[Path, ...],
    protein_ids: frozenset[str],
    feature_keys: frozenset[tuple[str, str]],
) -> dict[tuple[str, str, str], str]:
    """Read positive, assessed-negative and unknown query evidence strictly.

    Args:
        paths: Four-column query features TSVs or richer canonical feature tables.
        protein_ids: FASTA identifiers that may be assessed.
        feature_keys: Frozen reference definitions eligible for querying.

    Returns:
        State keyed by protein and reference feature.
    """

    result: dict[tuple[str, str, str], str] = {}
    ignored = 0
    for path in paths:
        for row in iter_tsv(
            path=path,
            required_fields=("protein_id", "feature_type", "feature_id", "evidence_status"),
            allow_empty=True,
        ):
            protein_id = row["protein_id"]
            if protein_id not in protein_ids:
                raise InputValidationError(
                    f"Query feature references an unknown protein: {protein_id!r}"
                )
            key = (row["feature_type"], row["feature_id"])
            status = row["evidence_status"].strip().upper()
            if status not in FEATURE_ASSESSMENT_STATUSES:
                raise InputValidationError(f"Unknown query feature assessment status: {status!r}")
            if key not in feature_keys:
                ignored += 1
                continue
            lookup = (protein_id, *key)
            previous = result.setdefault(lookup, status)
            if previous != status:
                raise InputValidationError(f"Conflicting query feature assessments: {lookup!r}")
    LOGGER.info(
        "Read %d frozen query assessments; ignored %d other feature rows",
        len(result),
        ignored,
    )
    return result


def _scan_reference_kmers(
    *, proteins: tuple[SequenceRecord, ...], feature_keys: frozenset[tuple[str, str]]
) -> dict[str, frozenset[str]]:
    """Assess frozen exact k-mers directly against every query sequence.

    Args:
        proteins: Validated query FASTA records.
        feature_keys: Validated reference feature type and ID pairs.

    Returns:
        Present reference k-mer IDs by query protein. Missing IDs are assessed absent.

    Raises:
        InputValidationError: If a frozen k-mer ID does not match its native definition.
    """

    reference_by_length: dict[int, set[str]] = defaultdict(set)
    for feature_type, feature_id in feature_keys:
        if feature_type != "AMINO_ACID_KMER":
            continue
        prefix, separator, residues = feature_id.partition(":")
        if not separator or not prefix.startswith("k") or not prefix[1:].isdigit():
            raise InputValidationError(f"Invalid frozen amino-acid k-mer ID: {feature_id!r}")
        length = int(prefix[1:])
        if not 1 <= length <= 12 or len(residues) != length or residues != residues.upper():
            raise InputValidationError(f"Invalid frozen amino-acid k-mer ID: {feature_id!r}")
        reference_by_length[length].add(feature_id)
    result: dict[str, frozenset[str]] = {}
    for protein in proteins:
        present: set[str] = set()
        for length, reference in reference_by_length.items():
            present.update(
                f"k{length}:{residues}"
                for residues in sequence_kmers(sequence=protein.sequence, length=length)
                if f"k{length}:{residues}" in reference
            )
        result[protein.protein_id] = frozenset(present)
    LOGGER.info(
        "Assessed %d exact k-mer definitions across %d query proteins",
        sum(map(len, reference_by_length.values())),
        len(proteins),
    )
    return result


def _load_query_units(
    *, path: Path | None, protein_ids: frozenset[str]
) -> dict[str, tuple[str, ...]]:
    """Map every query protein to itself or to one or more supplied groups.

    Args:
        path: Optional protein-to-unit membership TSV.
        protein_ids: Complete authoritative FASTA collection.

    Returns:
        Sorted member identifiers for each query unit.
    """

    if path is None:
        return {protein_id: (protein_id,) for protein_id in sorted(protein_ids)}
    result: dict[str, set[str]] = defaultdict(set)
    assigned: set[str] = set()
    for row in iter_tsv(path=path, required_fields=("protein_id", "unit_id")):
        protein_id = row["protein_id"]
        if protein_id not in protein_ids:
            raise InputValidationError(f"Query group has an unknown protein: {protein_id!r}")
        unit_id = validate_identifier(value=row["unit_id"], field_name="unit_id")
        if protein_id in result[unit_id]:
            raise InputValidationError(
                f"Duplicate protein and unit membership: {protein_id}, {unit_id}"
            )
        result[unit_id].add(protein_id)
        assigned.add(protein_id)
    missing = protein_ids - assigned
    if missing:
        raise InputValidationError(
            f"Query group table omits FASTA proteins: {sorted(missing)[:10]}"
        )
    return {unit_id: tuple(sorted(members)) for unit_id, members in sorted(result.items())}


def _iter_query_assessments(
    *,
    comparisons: Mapping[str, tuple[str, str]],
    frozen: Mapping[str, tuple[_FrozenFeature, ...]],
    units: Mapping[str, tuple[str, ...]],
    states: Mapping[tuple[str, str, str], str],
    native_kmers: Mapping[str, frozenset[str]],
    unit_type: str,
) -> Iterator[tuple[dict[str, object], tuple[dict[str, object], ...]]]:
    """Assess each unit and class using bounded per-unit feature rows.

    Args:
        comparisons: Frozen target/background label mapping.
        frozen: Validated positive reference features per comparison.
        units: Query protein or group members.
        states: Explicit feature assessments.
        native_kmers: FASTA-assessed positive exact k-mers by protein.
        unit_type: PROTEIN or GROUP output value.

    Yields:
        A summary and its positive feature evidence rows.
    """

    for unit_id, members in sorted(units.items()):
        for comparison_id, (targets, backgrounds) in sorted(comparisons.items()):
            present = 0
            absent = 0
            unknown = 0
            types: set[str] = set()
            hits: list[dict[str, object]] = []
            reference = frozen.get(comparison_id, ())
            for feature in reference:
                calls = tuple(
                    (
                        "ASSESSED_WITH_FEATURE"
                        if feature.feature_id in native_kmers[protein_id]
                        else "ASSESSED_NO_FEATURE"
                    )
                    if feature.feature_type == "AMINO_ACID_KMER"
                    else states.get((protein_id, feature.feature_type, feature.feature_id))
                    for protein_id in members
                )
                positive = tuple(
                    protein_id
                    for protein_id, call in zip(members, calls, strict=True)
                    if call == "ASSESSED_WITH_FEATURE"
                )
                if positive:
                    present += 1
                    types.add(feature.feature_type)
                    hits.append(
                        {
                            "unit_id": unit_id,
                            "comparison_id": comparison_id,
                            "feature_type": feature.feature_type,
                            "feature_id": feature.feature_id,
                            "feature_name": feature.feature_name,
                            "positive_protein_ids": "|".join(positive),
                            "validation_q_value": feature.q_value,
                            "validation_fdr_scope": feature.fdr_scope,
                            "validation_prevalence_difference": feature.prevalence_difference,
                            "evidence_class": feature.evidence_class,
                        }
                    )
                elif all(call == "ASSESSED_NO_FEATURE" for call in calls):
                    absent += 1
                else:
                    unknown += 1
            evidence_state = (
                "NO_VALIDATED_REFERENCE_FEATURES"
                if not reference
                else "FEATURE_EVIDENCE_PRESENT"
                if present
                else "UNKNOWN_ASSESSMENT"
                if unknown
                else "ASSESSED_NO_FEATURE"
            )
            yield (
                {
                    "unit_id": unit_id,
                    "unit_type": unit_type,
                    "member_count": len(members),
                    "comparison_id": comparison_id,
                    "target_label_ids": targets,
                    "background_label_ids": backgrounds,
                    "reference_feature_count": len(reference),
                    "present_feature_count": present,
                    "assessed_absent_feature_count": absent,
                    "unknown_feature_count": unknown,
                    "present_feature_types": "|".join(sorted(types)),
                    "evidence_state": evidence_state,
                },
                tuple(hits),
            )
