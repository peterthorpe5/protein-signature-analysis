"""Operator configuration for immutable DeepClust sources and a separate workspace."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import shutil
import tomllib
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from protein_signatures.errors import InputValidationError

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True, kw_only=True)
class WorkbenchConfig:
    """Trusted server paths, executable names and explicit resource limits."""

    workspace: str
    index: str = ""
    sequences: str = ""
    sequence_map: str = ""
    identifier_bridge: str = ""
    taxonomy: str = ""
    centroid_database: str = ""
    centroid_fasta: str = ""
    resource: str = ""
    release: str = "Operator-configured DeepClust release"
    identifier_prefixes: tuple[str, ...] = ("UNIPROT_", "NR_")
    diamond: str = "diamond"
    mafft: str = "mafft"
    fasttree: str = "FastTree"
    hmmscan: str = "hmmscan"
    pfam_hmm: str = ""
    threads: int = 4
    memory_gb: int = 4
    max_sequences: int = 2_000_000
    max_seed_sequences: int = 10_000
    max_alignment_sequences: int = 1_000
    max_alignment_cells: int = 5_000_000
    command_timeout: int = 86_400


def bounded_number(
    *, value: object, name: str, minimum: float, maximum: float, integer: bool = False
) -> int | float:
    """Validate a finite scientific parameter within inclusive bounds.

    Args:
        value: Requested value.
        name: Field shown in an error.
        minimum: Inclusive lower bound.
        maximum: Inclusive upper bound.
        integer: Require an integral value.

    Returns:
        Validated number.

    Raises:
        InputValidationError: If the value is invalid.
    """
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise InputValidationError(f"{name} must be numeric.") from error
    if isinstance(value, bool) or not math.isfinite(number) or not minimum <= number <= maximum:
        raise InputValidationError(f"{name} must be between {minimum} and {maximum}.")
    if integer and not number.is_integer():
        raise InputValidationError(f"{name} must be an integer.")
    return int(number) if integer else number


def load_config(*, path: Path) -> WorkbenchConfig:
    """Load a strict TOML profile with paths relative to the profile directory.

    Args:
        path: Operator-owned profile, supplied at application launch.

    Returns:
        Validated, absolute-path configuration.

    Raises:
        InputValidationError: If the profile cannot be read or is malformed.
    """
    source = Path(path).expanduser().resolve()
    try:
        with source.open(mode="rb") as handle:
            values = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise InputValidationError(f"Could not read workbench profile: {source}") from error
    allowed = {field.name for field in fields(WorkbenchConfig)}
    unknown = set(values) - allowed
    if unknown or not values.get("workspace"):
        raise InputValidationError(f"Profile needs workspace; unknown fields: {sorted(unknown)}")
    numbers = {
        "threads": 256,
        "memory_gb": 4096,
        "max_sequences": 100_000_000,
        "max_seed_sequences": 100_000,
        "max_alignment_sequences": 10_000,
        "max_alignment_cells": 100_000_000,
        "command_timeout": 604_800,
    }
    for name, maximum in numbers.items():
        if name in values:
            values[name] = bounded_number(
                value=values[name], name=name, minimum=1, maximum=maximum, integer=True
            )
    prefixes = values.get("identifier_prefixes", ["UNIPROT_", "NR_"])
    if not isinstance(prefixes, list) or any(
        not isinstance(value, str) or not value or any(c.isspace() for c in value)
        for value in prefixes
    ):
        raise InputValidationError("identifier_prefixes must contain non-empty namespace strings.")
    values["identifier_prefixes"] = tuple(prefixes)
    for name in allowed - numbers.keys() - {"identifier_prefixes"}:
        if name in values and (
            not isinstance(values[name], str) or any(ord(c) < 32 for c in values[name])
        ):
            raise InputValidationError(f"{name} must be a plain text string.")
    for name in (
        "workspace",
        "index",
        "sequences",
        "sequence_map",
        "identifier_bridge",
        "taxonomy",
        "centroid_database",
        "centroid_fasta",
        "resource",
        "pfam_hmm",
    ):
        if values.get(name):
            candidate = Path(values[name]).expanduser()
            values[name] = str((source.parent / candidate).resolve())
    config = WorkbenchConfig(**values)
    root = Path(config.workspace)
    for name in ("index", "sequences", "sequence_map", "centroid_fasta"):
        if getattr(config, name) and Path(getattr(config, name)).is_relative_to(root):
            raise InputValidationError(
                "Keep DeepClust source files outside the writable workspace."
            )
    if config.resource and (
        root.is_relative_to(Path(config.resource)) or Path(config.resource).is_relative_to(root)
    ):
        raise InputValidationError("The workbench and completed result must be separate trees.")
    LOGGER.info("Loaded workbench profile %s", source)
    return config


def config_identity(*, config: WorkbenchConfig) -> str:
    """Fingerprint settings and source metadata without hashing the global database.

    Args:
        config: Validated profile.

    Returns:
        SHA-256 of settings plus size and modification time of configured inputs.
    """
    document = {"settings": asdict(config), "sources": {}}
    for name in (
        "index",
        "sequences",
        "sequence_map",
        "identifier_bridge",
        "taxonomy",
        "centroid_database",
        "centroid_fasta",
        "pfam_hmm",
    ):
        path = Path(getattr(config, name)) if getattr(config, name) else None
        if path is not None and path.is_file():
            stat = path.stat()
            document["sources"][name] = {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    return hashlib.sha256(json.dumps(obj=document, sort_keys=True).encode()).hexdigest()


def configure_denbi(*, source: Path, destination: Path, workspace: Path) -> Path:
    """Import installed E3 resource paths into a separate signature-workbench profile.

    Args:
        source: Existing E3 research.toml on Sabrina.
        destination: New profile; existing files are protected.
        workspace: Separate writable signature-analysis directory on the data disk.

    Returns:
        Validated profile path.

    Raises:
        InputValidationError: If paths overlap, input is invalid or destination exists.
    """
    try:
        with source.open(mode="rb") as handle:
            previous = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise InputValidationError(f"Could not read the existing E3 profile: {source}") from error
    allowed = {field.name for field in fields(WorkbenchConfig)} - {"workspace", "resource"}
    values = {key: value for key, value in previous.items() if key in allowed}
    # Resolve inherited relative resource paths against their original profile.
    path_names = (
        "index",
        "sequences",
        "sequence_map",
        "identifier_bridge",
        "taxonomy",
        "centroid_database",
        "centroid_fasta",
        "pfam_hmm",
    )
    for name in path_names:
        if values.get(name):
            values[name] = str(
                (source.resolve().parent / Path(values[name]).expanduser()).resolve()
            )
    if previous.get("workspace") and workspace.resolve() == Path(previous["workspace"]).resolve():
        raise InputValidationError("Use a separate workspace from the E3 application.")
    values["workspace"] = str(workspace.expanduser().resolve())
    values.setdefault("max_sequences", 2_000_000)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with destination.open(mode="x", encoding="utf-8") as handle:
            handle.write("# Imported E3 paths; independent signature-workbench jobs.\n")
            for key, value in sorted(values.items()):
                handle.write(f"{key} = {json.dumps(obj=value)}\n")
        load_config(path=destination)
    except FileExistsError as error:
        raise InputValidationError(f"Profile already exists: {destination}") from error
    LOGGER.info("Created separate de.NBI profile %s", destination)
    return destination.resolve()


def diagnostics(*, config: WorkbenchConfig) -> list[dict[str, str]]:
    """Inspect capabilities using file metadata; no global scan or analysis is run.

    Args:
        config: Operator profile.

    Returns:
        Parallel capability, status and detail records.
    """
    rows = []
    for name in ("index", "sequences", "sequence_map", "centroid_database", "taxonomy", "pfam_hmm"):
        value = getattr(config, name)
        ready = bool(value) and Path(value).is_file() and os.access(value, os.R_OK)
        rows.append(
            {
                "capability": name,
                "status": "READY" if ready else "UNAVAILABLE",
                "detail": value or "Not configured",
            }
        )
    for name in ("diamond", "mafft", "fasttree", "hmmscan"):
        resolved = shutil.which(cmd=getattr(config, name))
        rows.append(
            {
                "capability": name,
                "status": "READY" if resolved else "UNAVAILABLE",
                "detail": resolved or getattr(config, name),
            }
        )
    return rows
