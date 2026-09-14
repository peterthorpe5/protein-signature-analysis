"""Prepare generic signature inputs from a completed E3 end-to-end workflow.

The bridge is deliberately read-only with respect to the predecessor run.  It
copies no scientific conclusions and never promotes upstream family hints into
reviewed labels.  Instead, it publishes exact sequences, Pfam assessment rows,
available coordinate models and a curation worksheet for an independent
protein-signature campaign.
"""

from __future__ import annotations

import hashlib
import logging
import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb

from .catalogue import _is_uniprot_accession, _output_inventory, _wrap_sequence
from .checksums import sha256_file
from .errors import InputValidationError, PublicationError
from .fasta import read_protein_fasta
from .io_utils import read_json, write_json_atomic, write_text_atomic, write_tsv_atomic
from .orthofinder import discover_orthofinder_layout
from .structural_resource import (
    resolve_structural_resource_root,
    verify_structural_resource_outputs,
)
from .validation import validate_identifier, validate_text

LOGGER = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}")
_SEQUENCE_CANDIDATES = (
    "orthology/tables/candidate_group_member_sequences.parquet",
    "tables/candidate_group_member_sequences.parquet",
)
_DOMAIN_HIT_CANDIDATES = ("tables/domain_hits.parquet",)
_DOMAIN_SUMMARY_CANDIDATES = ("tables/domain_summary.parquet",)
_ASSET_CANDIDATES = ("tables/reused_asset_manifest.parquet",)
_MODEL_QUALITY_CANDIDATES = ("tables/reused_model_quality.parquet",)
_COORDINATE_SUFFIXES = (".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz")


@dataclass(frozen=True)
class E3WorkflowPaths:
    """Resolved, checksum-verified predecessor authorities."""

    run_root: Path
    final_manifest: Path
    sequence_table: Path
    domain_hits_table: Path
    domain_summary_table: Path
    asset_manifest_table: Path
    model_quality_table: Path
    structural_resource: Path
    structural_stage_manifest: Path | None
    orthofinder_results: Path
    orthofinder_version: str
    orthofinder_source_mode: str
    orthofinder_log_path: Path | None
    orthofinder_completion_authority_paths: tuple[Path, ...]


@dataclass(frozen=True)
class SequencePreparation:
    """Unique sequences and merged upstream membership context."""

    sequences: Mapping[str, str]
    audit_records: tuple[dict[str, str], ...]
    skipped_unmapped_rows: int


@dataclass(frozen=True)
class DomainPreparation:
    """Pfam input rows and non-authoritative E3 curation hints."""

    domain_records: tuple[dict[str, str], ...]
    audit_by_protein: Mapping[str, dict[str, str]]


@dataclass(frozen=True)
class StructurePreparation:
    """Available coordinate inventory derived from the asset manifest."""

    structure_records: tuple[dict[str, str], ...]
    protein_ids: frozenset[str]
    eligible_protein_ids: frozenset[str]
    skipped_non_coordinate_rows: int
    skipped_unmatched_rows: int


def resolve_e3_workflow_paths(*, run_root: Path) -> E3WorkflowPaths:
    """Resolve and verify the completed predecessor inputs used by the bridge.

    Args:
        run_root: Completed E3 end-to-end workflow directory.

    Returns:
        Absolute paths for the supported predecessor authorities.

    Raises:
        InputValidationError: If completion, layout or checksums are invalid.
    """

    root = Path(run_root).expanduser().resolve()
    if not root.is_dir():
        raise InputValidationError(f"E3 workflow run directory does not exist: {root}")
    final_manifest = root / "11_app_ready" / "stage_manifest.json"
    _require_complete_manifest(path=final_manifest)

    stage05 = root / "05_orthology"
    sequence_table = _resolve_first_file(
        root=stage05,
        candidates=_SEQUENCE_CANDIDATES,
        label="candidate group-member sequences",
    )
    _verify_manifested_files(stage_root=stage05, paths=(sequence_table,))

    stage06 = root / "06_domains"
    domain_hits = _resolve_first_file(
        root=stage06,
        candidates=_DOMAIN_HIT_CANDIDATES,
        label="domain hits",
    )
    domain_summary = _resolve_first_file(
        root=stage06,
        candidates=_DOMAIN_SUMMARY_CANDIDATES,
        label="domain summary",
    )
    _verify_manifested_files(stage_root=stage06, paths=(domain_hits, domain_summary))

    stage09 = root / "09_ligandability"
    assets = _resolve_first_file(
        root=stage09,
        candidates=_ASSET_CANDIDATES,
        label="reused structure-asset manifest",
    )
    model_quality = _resolve_first_file(
        root=stage09,
        candidates=_MODEL_QUALITY_CANDIDATES,
        label="reused structure model-quality table",
    )
    _verify_manifested_files(stage_root=stage09, paths=(assets, model_quality))

    structural_resource = resolve_structural_resource_root(path=root)
    structural_manifest = read_json(path=structural_resource / "provenance" / "run_manifest.json")
    if not isinstance(structural_manifest, Mapping):
        raise InputValidationError("Structural run manifest must contain an object.")
    verify_structural_resource_outputs(
        resource_dir=structural_resource,
        manifest=structural_manifest,
    )
    structural_stage_manifest = None
    if "datasets" in structural_manifest:
        structural_stage_manifest = structural_resource.parent / "stage_manifest.json"

    orthofinder_results = root / "04_orthofinder" / "Results"
    if not orthofinder_results.is_dir():
        raise InputValidationError(
            "Completed E3 workflow lacks the expected OrthoFinder Results directory: "
            f"{orthofinder_results}"
        )
    orthofinder_layout = discover_orthofinder_layout(results_dir=orthofinder_results)
    return E3WorkflowPaths(
        run_root=root,
        final_manifest=final_manifest,
        sequence_table=sequence_table,
        domain_hits_table=domain_hits,
        domain_summary_table=domain_summary,
        asset_manifest_table=assets,
        model_quality_table=model_quality,
        structural_resource=structural_resource,
        structural_stage_manifest=structural_stage_manifest,
        orthofinder_results=orthofinder_results,
        orthofinder_version=orthofinder_layout.version,
        orthofinder_source_mode=orthofinder_layout.source_mode,
        orthofinder_log_path=orthofinder_layout.log_path,
        orthofinder_completion_authority_paths=(orthofinder_layout.completion_authority_paths),
    )


def prepare_e3_workflow_inputs(
    *,
    run_root: Path,
    output_dir: Path,
    minimum_mean_plddt: float = 50.0,
) -> Path:
    """Publish a conservative input-review bundle from a completed E3 run.

    The generated label table contains only ``UNMAPPED`` assignments.  A curator
    must copy and review it before campaign initialisation; upstream family and
    role values are retained solely as hints in a separate audit worksheet.

    Args:
        run_root: Completed E3 end-to-end workflow directory.
        output_dir: New directory for prepared input authorities.
        minimum_mean_plddt: Inclusive Foldseek eligibility threshold on the
            predecessor model-quality scale from zero to 100.

    Returns:
        Atomically published absolute output directory.

    Raises:
        InputValidationError: If predecessor data are incomplete or inconsistent.
        PublicationError: If publication cannot complete atomically.
    """

    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise PublicationError(f"E3 workflow input bundle already exists: {destination}")
    if isinstance(minimum_mean_plddt, bool):
        raise InputValidationError("minimum_mean_plddt must be a finite number from 0 to 100.")
    try:
        confidence_threshold = float(minimum_mean_plddt)
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            "minimum_mean_plddt must be a finite number from 0 to 100."
        ) from error
    if not math.isfinite(confidence_threshold) or not 0.0 <= confidence_threshold <= 100.0:
        raise InputValidationError("minimum_mean_plddt must be a finite number from 0 to 100.")
    paths = resolve_e3_workflow_paths(run_root=run_root)
    sequence_rows = _read_parquet_records(
        path=paths.sequence_table,
        required=(
            "parsed_accession",
            "protein_sequence",
            "sequence_length",
            "sequence_sha256",
        ),
    )
    sequences = _prepare_sequences(rows=sequence_rows)
    if not sequences.sequences:
        raise InputValidationError("No exact accession-bearing protein sequences were available.")

    domain_reference = (
        f"domain_hits_sha256={sha256_file(path=paths.domain_hits_table)};"
        f"domain_summary_sha256={sha256_file(path=paths.domain_summary_table)}"
    )
    domains = _prepare_domains(
        hit_rows=_read_parquet_records(
            path=paths.domain_hits_table,
            required=(
                "member_accession",
                "source_database",
                "entry_accession",
                "entry_name",
                "location_start",
                "location_end",
                "score",
            ),
            allow_empty=True,
        ),
        summary_rows=_read_parquet_records(
            path=paths.domain_summary_table,
            required=(
                "member_accession",
                "annotation_availability_status",
                "annotation_status_detail",
                "pfam_hit_count",
                "interpro_version",
            ),
        ),
        sequences=sequences.sequences,
        evidence_reference=domain_reference,
    )
    structures = _prepare_structures(
        rows=_read_parquet_records(
            path=paths.asset_manifest_table,
            required=("accession", "path", "sha256"),
        ),
        quality_rows=_read_parquet_records(
            path=paths.model_quality_table,
            required=("accession", "mean_plddt"),
            allow_empty=True,
        ),
        sequence_ids=frozenset(sequences.sequences),
        run_root=paths.run_root,
        asset_manifest_path=paths.asset_manifest_table,
        minimum_mean_plddt=confidence_threshold,
    )
    if len(structures.eligible_protein_ids) < 2:
        raise InputValidationError(
            "At least two sequence-matched, confidence-eligible coordinate models are required "
            f"for Foldseek; found {len(structures.eligible_protein_ids)} at mean pLDDT >= "
            f"{confidence_threshold:g}."
        )
    return _publish_bundle(
        destination=destination,
        paths=paths,
        sequences=sequences,
        domains=domains,
        structures=structures,
        minimum_mean_plddt=confidence_threshold,
    )


def _require_complete_manifest(*, path: Path) -> Mapping[str, Any]:
    """Read one manifest and require a lowercase complete status.

    Args:
        path: Manifest JSON path.

    Returns:
        Decoded manifest mapping.

    Raises:
        InputValidationError: If the document or completion state is invalid.
    """

    value = read_json(path=path)
    if not isinstance(value, Mapping):
        raise InputValidationError(f"Stage manifest must contain an object: {path}")
    if value.get("status") != "complete":
        raise InputValidationError(f"Stage manifest is not marked complete: {path}")
    return value


def _resolve_first_file(*, root: Path, candidates: Sequence[str], label: str) -> Path:
    """Resolve the first non-empty supported path beneath a stage root.

    Args:
        root: Stage directory.
        candidates: Relative paths in preference order.
        label: Human-readable authority name.

    Returns:
        Absolute selected file path.

    Raises:
        InputValidationError: If no candidate is a non-empty file.
    """

    for relative in candidates:
        candidate = (root / relative).resolve()
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise InputValidationError(
        f"Could not find {label} below {root}; expected one of: " + "; ".join(candidates)
    )


def _verify_manifested_files(*, stage_root: Path, paths: Sequence[Path]) -> None:
    """Require consumed stage files to match their checksum inventory.

    Args:
        stage_root: Stage directory containing ``stage_manifest.json``.
        paths: Consumed files that must occur exactly once in the inventory.

    Raises:
        InputValidationError: If an inventory entry, size or checksum is invalid.
    """

    manifest_path = stage_root / "stage_manifest.json"
    manifest = _require_complete_manifest(path=manifest_path)
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise InputValidationError(f"Stage manifest has no output inventory: {manifest_path}")
    for source in paths:
        matches: list[Mapping[str, Any]] = []
        for row in outputs:
            if not isinstance(row, Mapping):
                raise InputValidationError(f"Malformed output record in {manifest_path}")
            raw_path = row.get("path", row.get("relative_path"))
            relative = Path(str(raw_path or ""))
            if not relative.parts or relative.is_absolute() or ".." in relative.parts:
                raise InputValidationError(f"Unsafe output path in {manifest_path}: {raw_path!r}")
            candidates = {(stage_root / relative).resolve()}
            if relative.parts[0] == stage_root.name:
                candidates.add((stage_root.parent / relative).resolve())
            if source.resolve() in candidates:
                matches.append(row)
        if len(matches) != 1:
            raise InputValidationError(
                f"Expected exactly one checksum entry for {source}; found {len(matches)}."
            )
        record = matches[0]
        try:
            expected_size = int(record.get("size_bytes"))
        except (TypeError, ValueError) as error:
            raise InputValidationError(
                f"Invalid size_bytes for manifested output: {source}"
            ) from error
        expected_digest = str(record.get("sha256", ""))
        if expected_size < 0 or _DIGEST.fullmatch(expected_digest) is None:
            raise InputValidationError(f"Invalid checksum record for manifested output: {source}")
        if source.stat().st_size != expected_size:
            raise InputValidationError(f"Manifested output size mismatch: {source}")
        if sha256_file(path=source) != expected_digest:
            raise InputValidationError(f"Manifested output checksum mismatch: {source}")


def _read_parquet_records(
    *, path: Path, required: Sequence[str], allow_empty: bool = False
) -> tuple[dict[str, Any], ...]:
    """Read a Parquet authority after validating its named column contract.

    Args:
        path: Existing Parquet file.
        required: Required unique columns.
        allow_empty: Whether a zero-row table is valid.

    Returns:
        Ordered raw row dictionaries.

    Raises:
        InputValidationError: If the file, schema or rows are invalid.
    """

    source = Path(path).expanduser().resolve()
    fields = tuple(required)
    if len(fields) != len(set(fields)) or not fields:
        raise InputValidationError("Required Parquet fields must be non-empty and unique.")
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty Parquet authority: {source}")
    try:
        with duckdb.connect(":memory:") as connection:
            cursor = connection.execute("SELECT * FROM read_parquet(?)", [str(source)])
            columns = tuple(str(item[0]) for item in cursor.description)
            missing = sorted(set(fields).difference(columns))
            if missing:
                raise InputValidationError(
                    f"{source.name} is missing required columns: {', '.join(missing)}"
                )
            rows = tuple(dict(zip(columns, values, strict=True)) for values in cursor.fetchall())
    except duckdb.Error as error:
        raise InputValidationError(f"Could not read Parquet authority {source}: {error}") from error
    if not rows and not allow_empty:
        raise InputValidationError(f"Parquet authority contains no records: {source}")
    return rows


def _prepare_sequences(*, rows: Sequence[Mapping[str, Any]]) -> SequencePreparation:
    """Validate and deduplicate accession-bearing upstream sequences.

    Args:
        rows: Candidate group-member sequence rows.

    Returns:
        Unique sequences, merged group context and skipped-row count.

    Raises:
        InputValidationError: If sequence identity or checksum records conflict.
    """

    sequences: dict[str, str] = {}
    contexts: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    skipped = 0
    for row_number, row in enumerate(rows, start=1):
        raw_id = str(row.get("parsed_accession") or "").strip()
        if not raw_id:
            skipped += 1
            continue
        protein_id = validate_identifier(
            value=raw_id,
            field_name=f"parsed_accession row {row_number}",
        )
        sequence = str(row.get("protein_sequence") or "").strip().upper()
        if not sequence or any(character.isspace() for character in sequence):
            raise InputValidationError(
                f"Protein {protein_id!r} has an empty or whitespace-bearing sequence."
            )
        try:
            declared_length = int(row.get("sequence_length"))
        except (TypeError, ValueError) as error:
            raise InputValidationError(
                f"Protein {protein_id!r} has an invalid sequence_length."
            ) from error
        declared_digest = str(row.get("sequence_sha256") or "").strip().lower()
        try:
            encoded_sequence = sequence.encode("ascii", errors="strict")
        except UnicodeEncodeError as error:
            raise InputValidationError(
                f"Protein {protein_id!r} contains a non-ASCII residue symbol."
            ) from error
        observed_digest = hashlib.sha256(encoded_sequence).hexdigest()
        if declared_length != len(sequence) or declared_digest != observed_digest:
            raise InputValidationError(
                f"Protein {protein_id!r} has inconsistent length or SHA-256 metadata."
            )
        previous = sequences.setdefault(protein_id, sequence)
        if previous != sequence:
            raise InputValidationError(f"Conflicting sequences for accession {protein_id!r}.")
        for source_field, destination_field in (
            ("cluster_id", "cluster_ids"),
            ("group_id", "group_ids"),
            ("orthogroup_id", "orthogroup_ids"),
            ("species", "species"),
        ):
            value = str(row.get(source_field) or "").strip()
            if value:
                contexts[protein_id][destination_field].add(
                    validate_identifier(value=value, field_name=source_field)
                )
        candidate_raw = row.get("is_input_candidate")
        candidate_value = str(candidate_raw if candidate_raw is not None else "").strip().casefold()
        if candidate_value in {"1", "true", "yes"}:
            contexts[protein_id]["input_candidate_states"].add("TRUE")
        elif candidate_value in {"0", "false", "no"}:
            contexts[protein_id]["input_candidate_states"].add("FALSE")
        else:
            contexts[protein_id]["input_candidate_states"].add("UNKNOWN")
    audit = tuple(
        {
            "protein_id": protein_id,
            "cluster_ids": "|".join(sorted(contexts[protein_id]["cluster_ids"])),
            "group_ids": "|".join(sorted(contexts[protein_id]["group_ids"])),
            "orthogroup_ids": "|".join(sorted(contexts[protein_id]["orthogroup_ids"])),
            "species": "|".join(sorted(contexts[protein_id]["species"])),
            "input_candidate_states": "|".join(
                sorted(contexts[protein_id]["input_candidate_states"])
            ),
            "sequence_length": str(len(sequences[protein_id])),
            "sequence_sha256": hashlib.sha256(sequences[protein_id].encode("ascii")).hexdigest(),
        }
        for protein_id in sorted(sequences)
    )
    return SequencePreparation(
        sequences=sequences,
        audit_records=audit,
        skipped_unmapped_rows=skipped,
    )


def _prepare_domains(
    *,
    hit_rows: Sequence[Mapping[str, Any]],
    summary_rows: Sequence[Mapping[str, Any]],
    sequences: Mapping[str, str],
    evidence_reference: str,
) -> DomainPreparation:
    """Convert upstream InterPro/Pfam rows into the explicit domain ledger.

    Args:
        hit_rows: Upstream domain-hit rows.
        summary_rows: Upstream per-protein assessment rows.
        sequences: Exact protein sequences keyed by campaign identifier.
        evidence_reference: Checksum-bound predecessor reference.

    Returns:
        Canonical domain rows and curation hints.

    Raises:
        InputValidationError: If coordinates, counts or duplicate hits conflict.
    """

    source = "Completed E3 workflow InterPro/Pfam authority"
    reference = validate_text(
        value=evidence_reference,
        field_name="domain evidence_reference",
    )
    summary: dict[str, dict[str, Any]] = defaultdict(
        lambda: {"available": False, "declared_hits": 0, "details": set(), "versions": set()}
    )
    for row in summary_rows:
        protein_id = str(row.get("member_accession") or "").strip()
        if not protein_id or protein_id not in sequences:
            continue
        status = str(row.get("annotation_availability_status") or "").strip().upper()
        if status == "AVAILABLE":
            summary[protein_id]["available"] = True
        try:
            count = int(row.get("pfam_hit_count") or 0)
        except (TypeError, ValueError) as error:
            raise InputValidationError(
                f"Invalid pfam_hit_count for protein {protein_id!r}."
            ) from error
        if count < 0:
            raise InputValidationError(f"Negative pfam_hit_count for protein {protein_id!r}.")
        summary[protein_id]["declared_hits"] = max(int(summary[protein_id]["declared_hits"]), count)
        for field, destination in (
            ("annotation_status_detail", "details"),
            ("interpro_version", "versions"),
        ):
            value = str(row.get(field) or "").strip()
            if value:
                summary[protein_id][destination].add(validate_text(value=value, field_name=field))

    hits: dict[tuple[str, str, int, int], dict[str, str]] = {}
    hints: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in hit_rows:
        protein_id = str(row.get("member_accession") or "").strip()
        if not protein_id or protein_id not in sequences:
            continue
        for source_field, destination in (
            ("e3_family", "families"),
            ("evidence_role", "roles"),
        ):
            value = str(row.get(source_field) or "").strip()
            if value:
                hints[protein_id][destination].add(
                    validate_text(value=value, field_name=source_field)
                )
        if str(row.get("source_database") or "").strip().casefold() != "pfam":
            continue
        domain_id = validate_identifier(
            value=row.get("entry_accession"),
            field_name="Pfam entry_accession",
        )
        try:
            start = int(row.get("location_start"))
            end = int(row.get("location_end"))
        except (TypeError, ValueError) as error:
            raise InputValidationError(
                f"Invalid Pfam coordinates for {(protein_id, domain_id)!r}."
            ) from error
        if start < 1 or end < start or end > len(sequences[protein_id]):
            raise InputValidationError(
                f"Pfam coordinates exceed sequence bounds for {(protein_id, domain_id)!r}."
            )
        score = ""
        raw_score = row.get("score")
        if raw_score is not None and str(raw_score).strip():
            try:
                parsed_score = float(raw_score)
            except (TypeError, ValueError) as error:
                raise InputValidationError(
                    f"Invalid Pfam score for {(protein_id, domain_id)!r}."
                ) from error
            if not math.isfinite(parsed_score):
                raise InputValidationError(
                    f"Non-finite Pfam score for {(protein_id, domain_id)!r}."
                )
            score = format(parsed_score, ".15g")
        name = str(row.get("entry_name") or "").strip() or domain_id
        record = {
            "protein_id": protein_id,
            "domain_authority": "Pfam",
            "assessment_status": "ASSESSED_WITH_HIT",
            "domain_id": domain_id,
            "domain_name": validate_text(value=name, field_name="Pfam entry_name"),
            "start": str(start),
            "end": str(end),
            "score": score,
            "e_value": "",
            "evidence_source": source,
            "evidence_reference": reference,
        }
        key = (protein_id, domain_id, start, end)
        previous = hits.setdefault(key, record)
        if previous != record:
            raise InputValidationError(f"Conflicting duplicate Pfam hit: {key!r}")
        hints[protein_id]["pfam_ids"].add(domain_id)

    hits_by_protein: dict[str, list[dict[str, str]]] = defaultdict(list)
    for key, record in sorted(hits.items()):
        hits_by_protein[key[0]].append(record)
    records: list[dict[str, str]] = []
    audit: dict[str, dict[str, str]] = {}
    hit_proteins = {key[0] for key in hits}
    for protein_id in sorted(sequences):
        if protein_id in hit_proteins:
            records.extend(hits_by_protein[protein_id])
            assessment = "ASSESSED_WITH_HIT"
        elif bool(summary[protein_id]["available"]):
            if int(summary[protein_id]["declared_hits"]) > 0:
                raise InputValidationError(
                    f"Pfam summary declares hits but no hit rows exist for {protein_id!r}."
                )
            assessment = "ASSESSED_NO_HIT"
            records.append(
                {
                    "protein_id": protein_id,
                    "domain_authority": "Pfam",
                    "assessment_status": assessment,
                    "domain_id": "",
                    "domain_name": "",
                    "start": "",
                    "end": "",
                    "score": "",
                    "e_value": "",
                    "evidence_source": source,
                    "evidence_reference": reference,
                }
            )
        else:
            assessment = "NOT_ASSESSED"
            records.append(
                {
                    "protein_id": protein_id,
                    "domain_authority": "Pfam",
                    "assessment_status": assessment,
                    "domain_id": "",
                    "domain_name": "",
                    "start": "",
                    "end": "",
                    "score": "",
                    "e_value": "",
                    "evidence_source": source,
                    "evidence_reference": reference,
                }
            )
        audit[protein_id] = {
            "pfam_assessment_status": assessment,
            "pfam_accessions": "|".join(sorted(hints[protein_id]["pfam_ids"])),
            "upstream_e3_families": "|".join(sorted(hints[protein_id]["families"])),
            "upstream_evidence_roles": "|".join(sorted(hints[protein_id]["roles"])),
            "annotation_status_details": "|".join(sorted(summary[protein_id]["details"])),
            "interpro_versions": "|".join(sorted(summary[protein_id]["versions"])),
        }
    return DomainPreparation(
        domain_records=tuple(
            sorted(
                records,
                key=lambda item: (
                    item["protein_id"],
                    item["domain_id"],
                    int(item["start"] or 0),
                    int(item["end"] or 0),
                ),
            )
        ),
        audit_by_protein=audit,
    )


def _prepare_structures(
    *,
    rows: Sequence[Mapping[str, Any]],
    quality_rows: Sequence[Mapping[str, Any]],
    sequence_ids: frozenset[str],
    run_root: Path,
    asset_manifest_path: Path,
    minimum_mean_plddt: float,
) -> StructurePreparation:
    """Resolve and checksum upstream coordinate assets for Foldseek input.

    Args:
        rows: Upstream structure-asset manifest rows.
        quality_rows: Upstream per-accession model-confidence rows.
        sequence_ids: Exact campaign FASTA identifiers.
        run_root: Completed predecessor root.
        asset_manifest_path: Source asset-manifest path.
        minimum_mean_plddt: Inclusive eligibility threshold from zero to 100.

    Returns:
        Canonical structure rows plus exclusion counts.

    Raises:
        InputValidationError: If paths, confidence, checksums or identities conflict.
    """

    if isinstance(minimum_mean_plddt, bool):
        raise InputValidationError("minimum_mean_plddt must be a finite number from 0 to 100.")
    try:
        confidence_threshold = float(minimum_mean_plddt)
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            "minimum_mean_plddt must be a finite number from 0 to 100."
        ) from error
    if not math.isfinite(confidence_threshold) or not 0.0 <= confidence_threshold <= 100.0:
        raise InputValidationError("minimum_mean_plddt must be a finite number from 0 to 100.")
    quality_by_protein: dict[str, float | None] = {}
    for row_number, row in enumerate(quality_rows, start=1):
        protein_id = validate_identifier(
            value=row.get("accession"),
            field_name=f"model-quality accession row {row_number}",
        )
        if protein_id not in sequence_ids:
            continue
        raw_confidence = row.get("mean_plddt")
        if raw_confidence is None or not str(raw_confidence).strip():
            confidence = None
        else:
            if isinstance(raw_confidence, bool):
                raise InputValidationError(f"Invalid mean_plddt for protein {protein_id!r}.")
            try:
                confidence = float(raw_confidence)
            except (TypeError, ValueError) as error:
                raise InputValidationError(
                    f"Invalid mean_plddt for protein {protein_id!r}."
                ) from error
            if not math.isfinite(confidence) or not 0.0 <= confidence <= 100.0:
                raise InputValidationError(
                    f"mean_plddt must be finite and between 0 and 100 for {protein_id!r}."
                )
        if protein_id in quality_by_protein and quality_by_protein[protein_id] != confidence:
            raise InputValidationError(
                f"Conflicting model-quality rows exist for protein {protein_id!r}."
            )
        quality_by_protein[protein_id] = confidence

    selected: dict[str, dict[str, str]] = {}
    skipped_non_coordinate = 0
    skipped_unmatched = 0
    for row_number, row in enumerate(rows, start=1):
        raw_path = str(row.get("path") or "").strip()
        if not raw_path or not raw_path.casefold().endswith(_COORDINATE_SUFFIXES):
            skipped_non_coordinate += 1
            continue
        protein_id = validate_identifier(
            value=row.get("accession"),
            field_name=f"asset accession row {row_number}",
        )
        if protein_id not in sequence_ids:
            skipped_unmatched += 1
            continue
        coordinate = _resolve_asset_path(
            value=raw_path,
            run_root=run_root,
            asset_manifest_path=asset_manifest_path,
        )
        digest = str(row.get("sha256") or "").strip().lower()
        if _DIGEST.fullmatch(digest) is None:
            raise InputValidationError(f"Invalid coordinate SHA-256 for {protein_id!r}.")
        if digest != sha256_file(path=coordinate):
            raise InputValidationError(f"Coordinate checksum mismatch for {protein_id!r}.")
        raw_bytes = row.get("bytes")
        if raw_bytes is not None and str(raw_bytes).strip():
            try:
                declared_bytes = int(raw_bytes)
            except (TypeError, ValueError) as error:
                raise InputValidationError(
                    f"Invalid coordinate byte count for {protein_id!r}."
                ) from error
            if declared_bytes != coordinate.stat().st_size:
                raise InputValidationError(f"Coordinate byte-count mismatch for {protein_id!r}.")
        structure_key = hashlib.sha256(f"{protein_id}\t{digest}".encode("utf-8")).hexdigest()[:24]
        action = validate_text(
            value=str(row.get("action") or "").strip(),
            field_name="asset action",
            allow_empty=True,
        )
        structure_source = "Completed E3 workflow structure asset"
        if action:
            structure_source = f"{structure_source}; publication_action={action}"
        confidence = quality_by_protein.get(protein_id)
        if confidence is None:
            eligibility = "INELIGIBLE_CONFIDENCE_UNAVAILABLE"
            mean_confidence = ""
        elif confidence < confidence_threshold:
            eligibility = "INELIGIBLE_LOW_CONFIDENCE"
            mean_confidence = format(confidence, ".15g")
        else:
            eligibility = "ELIGIBLE"
            mean_confidence = format(confidence, ".15g")
        record = {
            "protein_id": protein_id,
            "structure_id": f"E3WF_{structure_key}",
            "structure_source": structure_source,
            "structure_version": "",
            "coordinate_path": str(coordinate),
            "coordinate_sha256": digest,
            "availability_status": "AVAILABLE",
            "mean_confidence": mean_confidence,
            "fold_id": "",
            "fold_name": "",
            "fold_authority": "",
            "fold_authority_version": "",
            "fold_evidence_reference": "",
            "fold_evidence_status": "NOT_ASSESSED",
            "analysis_eligibility_status": eligibility,
            "comparison_universe_ids": "",
        }
        previous = selected.setdefault(protein_id, record)
        if previous["coordinate_sha256"] != digest:
            raise InputValidationError(
                f"Multiple different coordinate assets exist for protein {protein_id!r}."
            )
        if record["coordinate_path"] < previous["coordinate_path"]:
            selected[protein_id] = record
    records = tuple(selected[protein_id] for protein_id in sorted(selected))
    eligible = frozenset(
        protein_id
        for protein_id, record in selected.items()
        if record["analysis_eligibility_status"] == "ELIGIBLE"
    )
    LOGGER.info(
        "Prepared %d structures (%d Foldseek eligible at mean pLDDT >= %g); skipped "
        "non-coordinate=%d unmatched=%d",
        len(records),
        len(eligible),
        confidence_threshold,
        skipped_non_coordinate,
        skipped_unmatched,
    )
    return StructurePreparation(
        structure_records=records,
        protein_ids=frozenset(selected),
        eligible_protein_ids=eligible,
        skipped_non_coordinate_rows=skipped_non_coordinate,
        skipped_unmatched_rows=skipped_unmatched,
    )


def _resolve_asset_path(*, value: str, run_root: Path, asset_manifest_path: Path) -> Path:
    """Resolve one absolute or known-relative upstream asset path.

    Args:
        value: Recorded path value.
        run_root: Completed predecessor root.
        asset_manifest_path: Asset manifest used to supply the value.

    Returns:
        Unique existing non-empty coordinate path.

    Raises:
        InputValidationError: If no unique file can be resolved.
    """

    raw = Path(value).expanduser()
    if raw.is_absolute():
        candidates = (raw.resolve(),)
    else:
        stage_root = asset_manifest_path.resolve().parents[1]
        candidates = tuple(
            dict.fromkeys(
                candidate.resolve()
                for candidate in (
                    asset_manifest_path.resolve().parent / raw,
                    stage_root / raw,
                    run_root.resolve() / raw,
                    run_root.resolve() / "09_ligandability" / raw,
                )
            )
        )
    matches = tuple(
        candidate
        for candidate in candidates
        if candidate.is_file() and candidate.stat().st_size > 0
    )
    if len(matches) != 1:
        raise InputValidationError(
            f"Expected one existing coordinate for asset path {value!r}; found {len(matches)}."
        )
    return matches[0]


def _publish_bundle(
    *,
    destination: Path,
    paths: E3WorkflowPaths,
    sequences: SequencePreparation,
    domains: DomainPreparation,
    structures: StructurePreparation,
    minimum_mean_plddt: float,
) -> Path:
    """Atomically write the prepared review bundle and completion marker.

    Args:
        destination: New final bundle directory.
        paths: Verified predecessor authorities.
        sequences: Prepared unique sequences and context.
        domains: Prepared Pfam ledger and hints.
        structures: Prepared coordinate inventory.
        minimum_mean_plddt: Confidence threshold used for Foldseek eligibility.

    Returns:
        Published absolute bundle directory.

    Raises:
        PublicationError: If output validation or publication fails.
    """

    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging.", dir=destination.parent))
    try:
        fasta_text = "".join(
            f">{protein_id}\n{_wrap_sequence(sequence=sequences.sequences[protein_id])}"
            for protein_id in sorted(sequences.sequences)
        )
        fasta_path = staging / "proteins.faa"
        write_text_atomic(path=fasta_path, text=fasta_text)
        parsed = read_protein_fasta(path=fasta_path)
        if len(parsed) != len(sequences.sequences):
            raise PublicationError("Prepared FASTA record count changed during validation.")

        write_tsv_atomic(
            path=staging / "domains.tsv",
            fieldnames=(
                "protein_id",
                "domain_authority",
                "assessment_status",
                "domain_id",
                "domain_name",
                "start",
                "end",
                "score",
                "e_value",
                "evidence_source",
                "evidence_reference",
            ),
            records=domains.domain_records,
        )
        write_tsv_atomic(
            path=staging / "structures.tsv",
            fieldnames=(
                "protein_id",
                "structure_id",
                "structure_source",
                "structure_version",
                "coordinate_path",
                "coordinate_sha256",
                "availability_status",
                "mean_confidence",
                "fold_id",
                "fold_name",
                "fold_authority",
                "fold_authority_version",
                "fold_evidence_reference",
                "fold_evidence_status",
                "analysis_eligibility_status",
                "comparison_universe_ids",
            ),
            records=structures.structure_records,
        )
        sequence_authority_digest = sha256_file(path=paths.sequence_table)
        label_records = tuple(
            {
                "protein_id": protein_id,
                "label_id": "e3:associated:unknown",
                "curation_status": "UNMAPPED",
                "evidence_status": "UNREVIEWED_PREDECESSOR_CONTEXT",
                "evidence_source": "Completed E3 workflow review bridge",
                "evidence_reference": sequence_authority_digest,
                "component_role": "UNKNOWN",
                "curation_reason": "Manual class, role and control review required.",
            }
            for protein_id in sorted(sequences.sequences)
        )
        write_tsv_atomic(
            path=staging / "label_assignments.REVIEW_REQUIRED.tsv",
            fieldnames=(
                "protein_id",
                "label_id",
                "curation_status",
                "evidence_status",
                "evidence_source",
                "evidence_reference",
                "component_role",
                "curation_reason",
            ),
            records=label_records,
        )
        sequence_context = {record["protein_id"]: record for record in sequences.audit_records}
        curation_records = tuple(
            {
                **sequence_context[protein_id],
                **domains.audit_by_protein[protein_id],
                "structure_available": (
                    "TRUE" if protein_id in structures.protein_ids else "FALSE"
                ),
                "structure_analysis_eligible": (
                    "TRUE" if protein_id in structures.eligible_protein_ids else "FALSE"
                ),
                "profile_label_id_to_assign": "",
                "curation_status_to_assign": "",
                "component_role_to_assign": "",
                "evidence_reference_to_assign": "",
                "curation_note": "",
            }
            for protein_id in sorted(sequences.sequences)
        )
        write_tsv_atomic(
            path=staging / "e3_label_curation_review.tsv",
            fieldnames=tuple(curation_records[0]),
            records=curation_records,
        )
        write_tsv_atomic(
            path=staging / "alphafold_accessions.MISSING_MODELS_REVIEW_REQUIRED.tsv",
            fieldnames=("protein_id", "uniprot_accession"),
            records=(
                {"protein_id": protein_id, "uniprot_accession": protein_id}
                for protein_id in sorted(sequences.sequences)
                if protein_id not in structures.protein_ids
                and _is_uniprot_accession(value=protein_id)
            ),
        )
        source_records = _source_inventory(paths=paths)
        write_tsv_atomic(
            path=staging / "source_inventory.tsv",
            fieldnames=("authority", "path", "size_bytes", "sha256"),
            records=source_records,
        )
        outputs = _output_inventory(root=staging)
        write_json_atomic(
            path=staging / "PREPARED.json",
            value={
                "schema_version": 1,
                "status": "COMPLETE",
                "source_run_root": str(paths.run_root),
                "structural_alignment_resource": str(paths.structural_resource),
                "orthofinder_results": str(paths.orthofinder_results),
                "orthofinder_version": paths.orthofinder_version,
                "orthofinder_source_mode": paths.orthofinder_source_mode,
                "protein_count": len(sequences.sequences),
                "pfam_record_count": len(domains.domain_records),
                "structure_count": len(structures.structure_records),
                "foldseek_eligible_structure_count": len(structures.eligible_protein_ids),
                "minimum_mean_plddt": minimum_mean_plddt,
                "training_eligible_assignment_count": 0,
                "skipped_unmapped_sequence_rows": sequences.skipped_unmapped_rows,
                "skipped_non_coordinate_asset_rows": (structures.skipped_non_coordinate_rows),
                "skipped_unmatched_asset_rows": structures.skipped_unmatched_rows,
                "next_action": "CURATE_LABEL_ASSIGNMENTS_AND_CONTROLS",
                "outputs": outputs,
            },
        )
        os.replace(staging, destination)
    except (InputValidationError, PublicationError):
        shutil.rmtree(staging, ignore_errors=True)
        raise
    except (OSError, UnicodeError) as error:
        shutil.rmtree(staging, ignore_errors=True)
        raise PublicationError(f"Could not publish E3 workflow input bundle: {error}") from error
    LOGGER.info(
        "Published E3 workflow input review bundle with %d proteins, %d structures and %d "
        "Foldseek-eligible structures at %s",
        len(sequences.sequences),
        len(structures.structure_records),
        len(structures.eligible_protein_ids),
        destination,
    )
    return destination


def _source_inventory(*, paths: E3WorkflowPaths) -> tuple[dict[str, str], ...]:
    """Build the checksum-bound list of predecessor files used by preparation.

    Args:
        paths: Resolved predecessor authorities.

    Returns:
        Deterministically ordered source inventory rows.
    """

    authorities = {
        "final_stage_manifest": paths.final_manifest,
        "sequence_table": paths.sequence_table,
        "domain_hits_table": paths.domain_hits_table,
        "domain_summary_table": paths.domain_summary_table,
        "asset_manifest_table": paths.asset_manifest_table,
        "model_quality_table": paths.model_quality_table,
        "structural_run_manifest": (paths.structural_resource / "provenance" / "run_manifest.json"),
    }
    if paths.structural_stage_manifest is not None:
        authorities["structural_stage_manifest"] = paths.structural_stage_manifest
    if paths.orthofinder_log_path is not None:
        authorities["orthofinder_log"] = paths.orthofinder_log_path
    completion_names = {
        "stage_manifest.json": "orthofinder_stage_manifest",
        "orthofinder_authority.tsv": "orthofinder_authority",
        "orthofinder_reuse_validation.tsv": "orthofinder_reuse_validation",
    }
    for path in paths.orthofinder_completion_authority_paths:
        authority = completion_names.get(path.name)
        if authority is None or authority in authorities:
            raise InputValidationError(f"Unexpected OrthoFinder completion authority path: {path}")
        authorities[authority] = path
    return tuple(
        {
            "authority": authority,
            "path": str(path),
            "size_bytes": str(path.stat().st_size),
            "sha256": sha256_file(path=path),
        }
        for authority, path in sorted(authorities.items())
    )
