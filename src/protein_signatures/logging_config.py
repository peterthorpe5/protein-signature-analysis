"""Consistent console and file logging configuration."""

from __future__ import annotations

import logging
from pathlib import Path

from .errors import ConfigurationError

_FORMAT = "%(asctime)s %(levelname)s %(name)s %(message)s"


def configure_logging(*, level: str = "INFO", log_path: Path | None = None) -> None:
    """Configure deterministic package logging.

    Args:
        level: Standard Python logging level name.
        log_path: Optional persistent log destination.

    Raises:
        ConfigurationError: If the level is unknown.
    """

    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ConfigurationError(f"Unsupported log level: {level!r}")
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_path is not None:
        destination = Path(log_path).expanduser().resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(destination, encoding="utf-8"))
    logging.basicConfig(
        level=numeric_level,
        format=_FORMAT,
        handlers=handlers,
        force=True,
    )
