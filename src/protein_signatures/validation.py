"""Reusable defensive validators for configuration and tabular input."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from typing import Any

from .errors import InputValidationError

_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:+/-]{0,254}$")


def validate_identifier(*, value: Any, field_name: str) -> str:
    """Validate a stable machine-readable identifier.

    Args:
        value: Candidate identifier.
        field_name: Human-readable field name for diagnostics.

    Returns:
        Stripped identifier.

    Raises:
        InputValidationError: If the identifier is empty or unsafe.
    """

    text = str(value).strip() if value is not None else ""
    if not _IDENTIFIER.fullmatch(text):
        raise InputValidationError(
            f"{field_name} must match {_IDENTIFIER.pattern!r}; received {text!r}."
        )
    return text


def validate_text(
    *,
    value: Any,
    field_name: str,
    allow_empty: bool = False,
    maximum_length: int = 10_000,
) -> str:
    """Validate bounded human-readable text without control characters.

    Args:
        value: Candidate text.
        field_name: Human-readable field name for diagnostics.
        allow_empty: Whether an empty value is accepted.
        maximum_length: Maximum accepted character count.

    Returns:
        Stripped text.

    Raises:
        InputValidationError: If the text violates the contract.
    """

    text = str(value).strip() if value is not None else ""
    if not text and not allow_empty:
        raise InputValidationError(f"{field_name} must not be empty.")
    if len(text) > maximum_length:
        raise InputValidationError(f"{field_name} exceeds the {maximum_length}-character limit.")
    if any(ord(character) < 32 and character not in "\t\n\r" for character in text):
        raise InputValidationError(f"{field_name} contains a control character.")
    if "\t" in text or "\n" in text or "\r" in text:
        raise InputValidationError(f"{field_name} contains a TSV-breaking character.")
    return text


def parse_optional_integer(
    *, value: Any, field_name: str, minimum: int | None = None
) -> int | None:
    """Parse an optional integer with a lower bound.

    Args:
        value: Candidate number or blank value.
        field_name: Field name for diagnostics.
        minimum: Optional inclusive lower bound.

    Returns:
        Parsed integer or ``None`` for blank input.

    Raises:
        InputValidationError: If parsing or range validation fails.
    """

    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = int(str(value).strip())
    except ValueError as error:
        raise InputValidationError(f"{field_name} must be an integer: {value!r}") from error
    if minimum is not None and parsed < minimum:
        raise InputValidationError(f"{field_name} must be at least {minimum}: {parsed}")
    return parsed


def parse_optional_float(
    *,
    value: Any,
    field_name: str,
    minimum: float | None = None,
    maximum: float | None = None,
) -> float | None:
    """Parse an optional finite floating-point value with bounds.

    Args:
        value: Candidate number or blank value.
        field_name: Field name for diagnostics.
        minimum: Optional inclusive lower bound.
        maximum: Optional inclusive upper bound.

    Returns:
        Parsed float or ``None`` for blank input.

    Raises:
        InputValidationError: If parsing, finiteness or bounds fail.
    """

    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = float(str(value).strip())
    except ValueError as error:
        raise InputValidationError(f"{field_name} must be numeric: {value!r}") from error
    if not math.isfinite(parsed):
        raise InputValidationError(f"{field_name} must be finite: {value!r}")
    if minimum is not None and parsed < minimum:
        raise InputValidationError(f"{field_name} must be at least {minimum}: {parsed}")
    if maximum is not None and parsed > maximum:
        raise InputValidationError(f"{field_name} must be at most {maximum}: {parsed}")
    return parsed


def require_mapping(*, value: Any, field_name: str) -> dict[str, Any]:
    """Return one configuration mapping or fail with context.

    Args:
        value: Candidate mapping.
        field_name: Configuration field name.

    Returns:
        Plain dictionary copy.

    Raises:
        InputValidationError: If the value is not a mapping.
    """

    if not isinstance(value, dict):
        raise InputValidationError(f"{field_name} must be a mapping.")
    return dict(value)


def require_sequence(*, value: Any, field_name: str) -> list[Any]:
    """Return one configuration sequence while rejecting text values.

    Args:
        value: Candidate list or tuple.
        field_name: Configuration field name.

    Returns:
        Plain list copy.

    Raises:
        InputValidationError: If the value is not a list or tuple.
    """

    if not isinstance(value, (list, tuple)):
        raise InputValidationError(f"{field_name} must be a list.")
    return list(value)


def reject_unknown_fields(
    *, value: Mapping[str, Any], allowed: frozenset[str], field_name: str
) -> None:
    """Reject misspelled or unsupported keys in a configuration mapping.

    Args:
        value: Already validated mapping.
        allowed: Exact accepted field names.
        field_name: Mapping name for diagnostics.

    Raises:
        InputValidationError: If one or more keys are not supported strings.
    """

    unknown = sorted(str(key) for key in value if not isinstance(key, str) or key not in allowed)
    if unknown:
        raise InputValidationError(f"{field_name} contains unknown fields: {unknown}.")
