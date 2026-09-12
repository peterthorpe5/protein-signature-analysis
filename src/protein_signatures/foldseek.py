"""Safe, cached Foldseek all-versus-all structural-alignment adapter."""

from __future__ import annotations

import csv
import logging
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .checksums import sha256_file, sha256_json
from .errors import ExternalToolError, InputValidationError
from .io_utils import open_text, read_json, write_json_atomic
from .models import (
    FoldseekRunEvidence,
    FoldseekSettings,
    PairwiseStructureComparison,
    StructureComparisonStatus,
    StructureCoverageScope,
    StructureRecord,
)
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)
FOLDSEEK_COLUMNS = (
    "query",
    "target",
    "alnlen",
    "qlen",
    "tlen",
    "evalue",
    "bits",
    "alntmscore",
    "qtmscore",
    "ttmscore",
)


def build_foldseek_command(
    *,
    executable: str,
    query_dir: Path,
    output_path: Path,
    temporary_dir: Path,
    settings: FoldseekSettings,
    threads: int,
) -> tuple[str, ...]:
    """Build an argument-safe Foldseek all-versus-all command.

    Args:
        executable: Resolved Foldseek executable.
        query_dir: Directory of uniquely named coordinate models.
        output_path: Raw tabular result path.
        temporary_dir: Scheduler-local or cache-local temporary directory.
        settings: Search thresholds and sensitivity.
        threads: Positive worker count.

    Returns:
        Command argument tuple that never invokes a shell.

    Raises:
        InputValidationError: If the worker count is invalid.
    """

    if threads < 1:
        raise InputValidationError("Foldseek threads must be positive.")
    return (
        executable,
        "easy-search",
        str(query_dir),
        str(query_dir),
        str(output_path),
        str(temporary_dir),
        "--alignment-type",
        "1",
        "-e",
        str(settings.e_value_threshold),
        "-s",
        str(settings.sensitivity),
        "--max-seqs",
        str(settings.maximum_hits),
        "--threads",
        str(threads),
        "--format-output",
        ",".join(FOLDSEEK_COLUMNS),
    )


def run_foldseek_all_vs_all(
    *,
    structures: tuple[StructureRecord, ...],
    settings: FoldseekSettings,
    threads: int,
) -> FoldseekRunEvidence:
    """Run or reuse a checksum-keyed Foldseek all-versus-all comparison.

    Args:
        structures: Structure inventory with local coordinate files.
        settings: Foldseek executable and cache settings.
        threads: Positive worker count.

    Returns:
        Parsed pairwise comparisons plus cache provenance.

    Raises:
        ExternalToolError: If Foldseek is unavailable or execution fails.
        InputValidationError: If fewer than two usable models are available.
    """

    usable = tuple(
        structure for structure in structures if structure.is_coordinate_analysis_eligible
    )
    excluded_count = len(structures) - len(usable)
    if excluded_count:
        LOGGER.info(
            "Excluded %d of %d structure records from Foldseek by analysis eligibility",
            excluded_count,
            len(structures),
        )
    usable_protein_ids = frozenset(structure.protein_id for structure in usable)
    if len(usable) < 2 or len(usable_protein_ids) < 2:
        raise InputValidationError(
            "Foldseek all-versus-all analysis requires explicitly eligible coordinate "
            "models for at least two distinct proteins."
        )
    if settings.maximum_hits < len(usable):
        raise InputValidationError(
            "foldseek.maximum_hits must be at least the number of eligible models "
            f"({len(usable)}) for an all-versus-all campaign; received "
            f"{settings.maximum_hits}."
        )
    executable = _resolve_executable(value=settings.executable)
    version = foldseek_version(executable=executable)
    cache_key = _cache_key(structures=usable, settings=settings, tool_version=version)
    comparison_universe_id = f"FOLDSEEK_{cache_key}"
    run_dir = settings.cache_dir / cache_key
    raw_output = run_dir / "foldseek_all_vs_all.tsv"
    completion = run_dir / "COMPLETED.json"
    query_map = {_query_id(structure_id=item.structure_id): item.protein_id for item in usable}
    if completion.is_file() and raw_output.is_file():
        marker = read_json(path=completion)
        if (
            isinstance(marker, dict)
            and marker.get("cache_key") == cache_key
            and marker.get("raw_output_sha256") == sha256_file(path=raw_output)
        ):
            comparisons = parse_foldseek_all_vs_all(
                path=raw_output,
                query_to_protein=query_map,
                tool_version=version,
                comparison_universe_id=comparison_universe_id,
            )
            LOGGER.info("Reused checksum-valid Foldseek cache unit %s", cache_key)
            return FoldseekRunEvidence(
                comparisons=comparisons,
                raw_output_path=raw_output,
                completion_manifest_path=completion,
                tool_version=version,
                cache_key=cache_key,
                comparison_universe_id=comparison_universe_id,
                assessment_universe=usable_protein_ids,
                reused=True,
            )
    settings.cache_dir.mkdir(parents=True, exist_ok=True)
    lock_path = settings.cache_dir / f".{cache_key}.lock"
    lock_descriptor = _acquire_lock(path=lock_path)
    staging = Path(tempfile.mkdtemp(prefix=f".{cache_key}.staging.", dir=settings.cache_dir))
    try:
        query_dir = staging / "queries"
        query_dir.mkdir()
        _link_queries(structures=usable, query_dir=query_dir)
        raw_staging = staging / raw_output.name
        temporary_base = Path(os.environ.get("TMPDIR", str(settings.cache_dir))).resolve()
        temporary_base.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f"protein_signatures_foldseek_{cache_key[:8]}_", dir=temporary_base
        ) as temporary_name:
            command = build_foldseek_command(
                executable=executable,
                query_dir=query_dir,
                output_path=raw_staging,
                temporary_dir=Path(temporary_name),
                settings=settings,
                threads=threads,
            )
            LOGGER.info("Running Foldseek cache unit %s with %d models", cache_key, len(usable))
            process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.returncode != 0:
            diagnostic = (process.stderr or process.stdout).strip()[-2000:]
            raise ExternalToolError(f"Foldseek exited with code {process.returncode}: {diagnostic}")
        if not raw_staging.is_file() or raw_staging.stat().st_size == 0:
            raise ExternalToolError("Foldseek completed without a non-empty result table.")
        comparisons = parse_foldseek_all_vs_all(
            path=raw_staging,
            query_to_protein=query_map,
            tool_version=version,
            comparison_universe_id=comparison_universe_id,
        )
        write_json_atomic(
            path=staging / completion.name,
            value={
                "status": "COMPLETE",
                "cache_key": cache_key,
                "raw_output_sha256": sha256_file(path=raw_staging),
                "tool": "Foldseek",
                "tool_version": version,
                "comparison_universe_id": comparison_universe_id,
                "command": list(command),
            },
        )
        if run_dir.exists():
            shutil.rmtree(run_dir)
        os.replace(staging, run_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        os.close(lock_descriptor)
        lock_path.unlink(missing_ok=True)
    return FoldseekRunEvidence(
        comparisons=comparisons,
        raw_output_path=raw_output,
        completion_manifest_path=completion,
        tool_version=version,
        cache_key=cache_key,
        comparison_universe_id=comparison_universe_id,
        assessment_universe=usable_protein_ids,
        reused=False,
    )


def parse_foldseek_all_vs_all(
    *,
    path: Path,
    query_to_protein: dict[str, str],
    tool_version: str,
    comparison_universe_id: str,
) -> tuple[PairwiseStructureComparison, ...]:
    """Parse and symmetrise one declared Foldseek tabular format.

    Args:
        path: Raw headerless Foldseek result.
        query_to_protein: Query identifiers mapped to campaign proteins.
        tool_version: Recorded Foldseek version.
        comparison_universe_id: Stable identifier for this complete search universe.

    Returns:
        Best conservative comparison per unordered protein pair.

    Raises:
        InputValidationError: If rows or identifiers violate the adapter contract.
    """

    universe_id = validate_identifier(
        value=comparison_universe_id,
        field_name="comparison_universe_id",
    )
    best: dict[tuple[str, str], PairwiseStructureComparison] = {}
    with open_text(path=path) as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=FOLDSEEK_COLUMNS)
        for row_number, row in enumerate(reader, start=1):
            if None in row or any(value is None for value in row.values()):
                raise InputValidationError(
                    f"Malformed Foldseek field count at row {row_number}: {path}"
                )
            query = _resolve_query(value=row["query"], query_to_protein=query_to_protein)
            target = _resolve_query(value=row["target"], query_to_protein=query_to_protein)
            if query == target:
                continue
            aligned = _parse_integer(value=row["alnlen"], field="alnlen", row_number=row_number)
            query_length = _parse_integer(value=row["qlen"], field="qlen", row_number=row_number)
            target_length = _parse_integer(value=row["tlen"], field="tlen", row_number=row_number)
            query_tm = _parse_unit_float(
                value=row["qtmscore"], field="qtmscore", row_number=row_number
            )
            target_tm = _parse_unit_float(
                value=row["ttmscore"], field="ttmscore", row_number=row_number
            )
            pair = tuple(sorted((query, target)))
            conservative_tm = min(query_tm, target_tm)
            record_id = f"foldseek:{pair[0]}:{pair[1]}"
            candidate = PairwiseStructureComparison(
                protein_a_id=pair[0],
                protein_b_id=pair[1],
                comparison_tool="Foldseek",
                comparison_tool_version=tool_version,
                tm_score=conservative_tm,
                rmsd_angstrom=None,
                aligned_residue_count=aligned,
                coverage_a=min(
                    1.0,
                    aligned / query_length if query == pair[0] else aligned / target_length,
                ),
                coverage_b=min(
                    1.0,
                    aligned / target_length if target == pair[1] else aligned / query_length,
                ),
                comparison_status=StructureComparisonStatus.COMPLETE,
                source_record_id=record_id,
                comparison_universe_id=universe_id,
                coverage_scope=StructureCoverageScope.STRUCTURE_MODEL_RESIDUES,
            )
            previous = best.get(pair)
            if previous is None or conservative_tm > float(previous.tm_score or -1.0):
                best[pair] = candidate
    return tuple(best[key] for key in sorted(best))


def foldseek_version(*, executable: str) -> str:
    """Read a bounded Foldseek version string.

    Args:
        executable: Resolved executable path.

    Returns:
        Non-empty version text.

    Raises:
        ExternalToolError: If the version command fails.
    """

    process = subprocess.run((executable, "version"), capture_output=True, text=True, check=False)
    version = (process.stdout or process.stderr).strip().splitlines()
    if process.returncode != 0 or not version:
        raise ExternalToolError(
            f"Could not determine Foldseek version (exit {process.returncode})."
        )
    return version[0][:200]


def _cache_key(
    *,
    structures: tuple[StructureRecord, ...],
    settings: FoldseekSettings,
    tool_version: str,
) -> str:
    """Calculate the structural work-unit cache key.

    Args:
        structures: Usable coordinate records.
        settings: Foldseek search settings.
        tool_version: Resolved external-tool version.

    Returns:
        SHA-256 cache key.
    """

    return sha256_json(
        value={
            "structures": [
                {
                    "structure_id": item.structure_id,
                    "coordinate_sha256": item.coordinate_sha256
                    or sha256_file(path=Path(item.coordinate_path)),
                }
                for item in sorted(structures, key=lambda record: record.structure_id)
            ],
            "tool_version": tool_version,
            "e_value_threshold": settings.e_value_threshold,
            "sensitivity": settings.sensitivity,
            "maximum_hits": settings.maximum_hits,
            "format": FOLDSEEK_COLUMNS,
        }
    )


def _query_id(*, structure_id: str) -> str:
    """Create a filesystem-safe opaque query identifier.

    Args:
        structure_id: Stable structure identifier.

    Returns:
        Opaque Foldseek query identifier.
    """

    return f"PSQ_{sha256_json(value=structure_id)[:20]}"


def _link_queries(*, structures: tuple[StructureRecord, ...], query_dir: Path) -> None:
    """Create uniquely named read-only symlinks for Foldseek input.

    Args:
        structures: Usable coordinate records.
        query_dir: Empty staging query directory.
    """

    for structure in structures:
        source = Path(structure.coordinate_path).resolve()
        suffix = ".cif" if ".cif" in source.name.lower() else ".pdb"
        destination = query_dir / f"{_query_id(structure_id=structure.structure_id)}{suffix}"
        destination.symlink_to(source)


def _resolve_query(*, value: str, query_to_protein: dict[str, str]) -> str:
    """Resolve a Foldseek chain/file identifier to a campaign protein.

    Args:
        value: Raw Foldseek identifier.
        query_to_protein: Opaque query-to-protein mapping.

    Returns:
        Campaign protein identifier.

    Raises:
        InputValidationError: If no unambiguous query prefix matches.
    """

    identifier = value.strip()
    matches = [
        protein_id
        for query_id, protein_id in query_to_protein.items()
        if identifier == query_id
        or identifier.startswith(f"{query_id}.")
        or identifier.startswith(f"{query_id}_")
    ]
    if len(matches) != 1:
        raise InputValidationError(f"Unknown or ambiguous Foldseek query identifier: {value!r}")
    return matches[0]


def _parse_integer(*, value: str, field: str, row_number: int) -> int:
    """Parse one positive Foldseek integer field.

    Args:
        value: Raw text.
        field: Field name.
        row_number: Source row number.

    Returns:
        Positive integer.
    """

    try:
        parsed = int(value)
    except ValueError as error:
        raise InputValidationError(
            f"Foldseek {field} must be an integer at row {row_number}: {value!r}"
        ) from error
    if parsed < 1:
        raise InputValidationError(
            f"Foldseek {field} must be positive at row {row_number}: {parsed}"
        )
    return parsed


def _parse_unit_float(*, value: str, field: str, row_number: int) -> float:
    """Parse one finite Foldseek score from zero to one.

    Args:
        value: Raw text.
        field: Field name.
        row_number: Source row number.

    Returns:
        Valid score.
    """

    try:
        parsed = float(value)
    except ValueError as error:
        raise InputValidationError(
            f"Foldseek {field} must be numeric at row {row_number}: {value!r}"
        ) from error
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise InputValidationError(
            f"Foldseek {field} must be from zero to one at row {row_number}: {parsed}"
        )
    return parsed


def _resolve_executable(*, value: str) -> str:
    """Resolve a configured Foldseek executable without invoking a shell.

    Args:
        value: Command name or absolute executable path.

    Returns:
        Executable path.

    Raises:
        ExternalToolError: If no executable is available.
    """

    resolved = shutil.which(value)
    if resolved is None:
        raise ExternalToolError(
            f"Foldseek executable {value!r} was not found; install it or disable foldseek."
        )
    return resolved


def _acquire_lock(*, path: Path) -> int:
    """Acquire an exclusive structural cache lock.

    Args:
        path: Lock-file path.

    Returns:
        Open lock descriptor.

    Raises:
        ExternalToolError: If another process owns the unit.
    """

    try:
        return os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise ExternalToolError(f"Foldseek cache unit is already running: {path}") from error
