"""Console launcher for the Streamlit protein-signature application."""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

from protein_signatures.errors import ProteinSignatureError

from .backend import app_database_path


def build_parser() -> argparse.ArgumentParser:
    """Build the application launcher parser.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(prog="protein-signature-app")
    parser.add_argument("--resource", required=True, type=Path)
    parser.add_argument("--port", default=8501, type=int)
    parser.add_argument("--address", default="localhost")
    return parser


def build_streamlit_command(*, database: Path, port: int, address: str) -> tuple[str, ...]:
    """Build the argument-safe Streamlit launch command.

    Args:
        database: Verified DuckDB path.
        port: TCP port from 1 to 65535.
        address: Streamlit bind address.

    Returns:
        Subprocess argument tuple.

    Raises:
        ValueError: If the port or address is invalid.
    """

    if not 1 <= port <= 65535:
        raise ValueError("port must be between 1 and 65535")
    if not address or any(character.isspace() for character in address):
        raise ValueError("address must be non-empty and contain no whitespace")
    app_path = Path(__file__).with_name("app.py").resolve()
    return (
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_path),
        "--server.port",
        str(port),
        "--server.address",
        address,
        "--",
        "--resource",
        str(database.parent),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Check required paths and launch the verified Streamlit viewer.

    Args:
        argv: Optional launcher arguments.

    Returns:
        Child process exit code or validation failure code.
    """

    arguments = build_parser().parse_args(argv)
    try:
        database = app_database_path(resource=arguments.resource)
        command = build_streamlit_command(
            database=database, port=arguments.port, address=arguments.address
        )
        return subprocess.run(command, check=False).returncode
    except (ProteinSignatureError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
