"""Console launcher for the Streamlit protein-signature application."""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

from protein_signatures.errors import ProteinSignatureError
from protein_signatures.workbench.config import load_config

from .backend import app_database_path


def build_parser() -> argparse.ArgumentParser:
    """Build the application launcher parser.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(prog="protein-signature-app")
    parser.add_argument("--resource", type=Path)
    parser.add_argument("--workbench-config", type=Path)
    parser.add_argument("--port", default=8501, type=int)
    parser.add_argument("--address", default="localhost")
    return parser


def build_streamlit_command(
    *, database: Path | None, port: int, address: str, workbench_config: Path | None = None
) -> tuple[str, ...]:
    """Build the argument-safe Streamlit launch command.

    Args:
        database: Verified DuckDB path.
        port: TCP port from 1 to 65535.
        address: Streamlit bind address.
        workbench_config: Optional trusted server profile for independent analysis jobs.

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
    command = (
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
    )
    if database is not None:
        command += ("--resource", str(database.parent))
    if workbench_config is not None:
        command += ("--workbench-config", str(workbench_config.expanduser().resolve()))
    return command


def main(argv: Sequence[str] | None = None) -> int:
    """Check required paths and launch the verified Streamlit viewer.

    Args:
        argv: Optional launcher arguments.

    Returns:
        Child process exit code or validation failure code.
    """

    arguments = build_parser().parse_args(argv)
    try:
        if arguments.resource is None and arguments.workbench_config is None:
            raise ValueError("Supply --resource, --workbench-config, or both.")
        if arguments.workbench_config is not None:
            profile = load_config(path=arguments.workbench_config)
            if arguments.resource is None and profile.resource:
                arguments.resource = Path(profile.resource)
        database = app_database_path(resource=arguments.resource) if arguments.resource else None
        command = build_streamlit_command(
            database=database,
            port=arguments.port,
            address=arguments.address,
            workbench_config=arguments.workbench_config,
        )
        return subprocess.run(command, check=False).returncode
    except (ProteinSignatureError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
