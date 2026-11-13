"""Operator configuration for immutable DeepClust sources and a separate workspace."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import shutil
import tomllib
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from urllib.parse import urlsplit

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
    species_taxonomy: str = ""
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
    foldseek: str = "foldseek"
    interproscan: str = "interproscan.sh"
    emapper: str = "emapper.py"
    eggnog_data: str = ""
    interpro_email: str = ""
    esm_endpoint: str = ""
    esm_format: str = "sequence"
    esm_token_env: str = ""
    esm_max_length: int = 400
    max_structure_residues: int = 3000
    threads: int = 4
    memory_gb: int = 4
    max_sequences: int = 2_000_000
    max_seed_sequences: int = 10_000
    max_alignment_sequences: int = 1_000
    max_alignment_cells: int = 5_000_000
    max_subcluster_sequences: int = 50_000
    max_subcluster_residues: int = 20_000_000
    subcluster_ssh_host: str = ""
    subcluster_ssh_executable: str = ""
    subcluster_ssh_profile: str = ""
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
        "max_subcluster_sequences": 500_000,
        "max_subcluster_residues": 1_000_000_000,
        "command_timeout": 604_800,
        "esm_max_length": 1024,
        "max_structure_residues": 10000,
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
        "species_taxonomy",
        "centroid_database",
        "centroid_fasta",
        "resource",
        "pfam_hmm",
        "eggnog_data",
    ):
        if values.get(name):
            candidate = Path(values[name]).expanduser()
            values[name] = str((source.parent / candidate).resolve())
    config = WorkbenchConfig(**values)
    if config.esm_format not in {"sequence", "json_inputs"}:
        raise InputValidationError("esm_format must be sequence or json_inputs.")
    if config.esm_endpoint:
        https_endpoint(value=config.esm_endpoint)
    if config.esm_token_env and (
        not config.esm_token_env.isascii() or not config.esm_token_env.isidentifier()
    ):
        raise InputValidationError("esm_token_env must be an environment variable name.")
    ssh_fields = (
        config.subcluster_ssh_host,
        config.subcluster_ssh_executable,
        config.subcluster_ssh_profile,
    )
    if any(ssh_fields):
        validate_subcluster_ssh(config=config)
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


def validate_subcluster_ssh(*, config: WorkbenchConfig) -> None:
    """Require a trusted SSH alias and two absolute shell-safe remote paths."""
    if not re.fullmatch(
        r"(?:[A-Za-z0-9_][A-Za-z0-9_.-]*@)?[A-Za-z0-9][A-Za-z0-9_.-]*",
        config.subcluster_ssh_host,
    ) or not all(
        re.fullmatch(r"/[A-Za-z0-9_./-]+", value)
        for value in (
            config.subcluster_ssh_executable,
            config.subcluster_ssh_profile,
        )
    ):
        raise InputValidationError(
            "Set an SSH host alias plus absolute executable and profile paths without spaces."
        )


def https_endpoint(*, value: str) -> str:
    """Validate an operator-owned HTTPS endpoint without inline credentials.

    Args:
        value: Configured provider URL.

    Returns:
        Validated URL.

    Raises:
        InputValidationError: If the URL contains credentials, a query or a fragment.
    """
    try:
        endpoint = urlsplit(value)
        _ = endpoint.port
    except ValueError as error:
        raise InputValidationError("Invalid protein provider URL.") from error
    if (
        endpoint.scheme != "https"
        or not endpoint.hostname
        or endpoint.username is not None
        or endpoint.password is not None
        or endpoint.query
        or endpoint.fragment
        or any(c.isspace() or ord(c) < 32 for c in value)
    ):
        raise InputValidationError("Protein providers require HTTPS without inline credentials.")
    return value


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
        "species_taxonomy",
        "centroid_database",
        "centroid_fasta",
        "pfam_hmm",
        "eggnog_data",
    ):
        path = Path(getattr(config, name)) if getattr(config, name) else None
        if path is not None and (path.is_file() or path.is_dir()):
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
    for old, new in (
        ("email", "interpro_email"),
        ("endpoint_max_length", "esm_max_length"),
    ):
        if old in previous:
            values.setdefault(new, previous[old])
    if values.get("esm_endpoint"):
        values.setdefault("esm_format", "json_inputs")
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
        "eggnog_data",
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


def configure_mount(
    *,
    root: Path,
    destination: Path,
    workspace: Path,
    release: str,
    resource: Path | None = None,
    centroid_database: Path | None = None,
) -> Path:
    """Configure exact files from an installed DeepClust data directory.

    Args:
        root: Mounted release directory containing the three official Parquet files.
        destination: New TOML profile; never overwritten.
        workspace: Independent writable job root.
        release: Operator-supplied database release name.
        resource: Optional immutable signature result.
        centroid_database: Explicit matching representative DIAMOND database, if present.

    Returns:
        Absolute path of the newly written profile.

    Raises:
        InputValidationError: If required files are unavailable or paths overlap.
    """
    root, workspace, destination = (
        path.expanduser().resolve() for path in (root, workspace, destination)
    )
    if not release.strip() or any(ord(c) < 32 for c in release):
        raise InputValidationError("Supply a plain-text DeepClust release label.")
    if workspace.is_relative_to(root) or destination.is_relative_to(root):
        raise InputValidationError("Keep the profile and workspace outside the read-only mount.")
    values = {"workspace": str(workspace), "release": release}
    for name, filename in (
        ("index", "clust_index_RowGroup.parquet"),
        ("sequences", "joined_with_index_RowGroupFinal.parquet"),
        ("sequence_map", "SeqIdMapClustId.parquet"),
    ):
        path = root / filename
        if not path.is_file() or not os.access(path, os.R_OK):
            raise InputValidationError(f"Missing {path}; specify the directory containing it.")
        if path.resolve().is_relative_to(workspace):
            raise InputValidationError("DeepClust inputs must be outside the writable workspace.")
        values[name] = str(path.resolve())
    if resource is not None:
        resource = resource.expanduser().resolve()
        if (
            not resource.is_dir()
            or workspace.is_relative_to(resource)
            or (resource.is_relative_to(workspace) or destination.is_relative_to(resource))
        ):
            raise InputValidationError("Use a separate existing completed result directory.")
        values["resource"] = str(resource)
    if centroid_database is not None:
        centroid_database = centroid_database.expanduser().resolve()
        if not centroid_database.is_file() or centroid_database.suffix != ".dmnd":
            raise InputValidationError("Supply an existing representative .dmnd database.")
        values["centroid_database"] = str(centroid_database)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with destination.open(mode="x", encoding="utf-8") as handle:
            for key, value in values.items():
                handle.write(f"{key} = {json.dumps(obj=value, ensure_ascii=False)}\n")
    except FileExistsError as error:
        raise InputValidationError(f"Profile already exists: {destination}") from error
    LOGGER.info("Created mounted DeepClust profile %s", destination)
    return destination


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
    for name in ("diamond", "mafft", "fasttree", "hmmscan", "foldseek", "interproscan", "emapper"):
        resolved = shutil.which(cmd=getattr(config, name))
        rows.append(
            {
                "capability": name,
                "status": "READY" if resolved else "UNAVAILABLE",
                "detail": resolved or getattr(config, name),
            }
        )
    rows.append(
        {
            "capability": "esmfold_api",
            "status": "CONFIGURED" if config.esm_endpoint else "UNAVAILABLE",
            "detail": config.esm_endpoint or "No prediction endpoint configured",
        }
    )
    return rows
