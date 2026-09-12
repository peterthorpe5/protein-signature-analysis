"""Content checksums and deterministic serialisation helpers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .errors import InputValidationError


def sha256_file(*, path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Calculate the SHA-256 digest of a readable regular file.

    Args:
        path: Input file.
        chunk_size: Number of bytes read per iteration.

    Returns:
        Lower-case hexadecimal digest.

    Raises:
        InputValidationError: If the path or chunk size is invalid.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise InputValidationError(f"Checksum input is not a regular file: {source}")
    if chunk_size < 1:
        raise InputValidationError(f"chunk_size must be positive: {chunk_size}")
    digest = hashlib.sha256()
    with source.open(mode="rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(*, text: str) -> str:
    """Calculate a SHA-256 digest for UTF-8 text.

    Args:
        text: Source text.

    Returns:
        Lower-case hexadecimal digest.
    """

    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def stable_json(*, value: Any) -> str:
    """Serialise JSON deterministically for manifests and cache keys.

    Args:
        value: JSON-compatible value.

    Returns:
        Canonical compact JSON text.
    """

    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_json(*, value: Any) -> str:
    """Calculate a checksum over canonical JSON.

    Args:
        value: JSON-compatible value.

    Returns:
        Lower-case hexadecimal digest.
    """

    return sha256_text(text=stable_json(value=value))
