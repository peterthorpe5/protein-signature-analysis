"""Named-option interface for durable server-side workbench jobs."""

from __future__ import annotations

import argparse
import json
import logging
import signal
import sys
import time
from collections.abc import Sequence
from pathlib import Path

from protein_signatures.errors import ProteinSignatureError

from .config import configure_denbi, diagnostics, load_config
from .jobs import OPERATIONS, cancel_job, list_jobs, recover_interrupted_job, run_next, submit_job


def build_parser() -> argparse.ArgumentParser:
    """Build a parser whose actions and scientific inputs are all named options.

    Returns:
        Workbench command-line parser.
    """
    parser = argparse.ArgumentParser(prog="protein-signature-workbench")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--action",
        required=True,
        choices=("configure-denbi", "doctor", "submit", "worker", "status", "cancel", "recover"),
    )
    parser.add_argument(
        "--source-config",
        type=Path,
        default=Path("/home/ubuntu/data/protein_research/config/research.toml"),
    )
    parser.add_argument(
        "--workspace", type=Path, default=Path("/home/ubuntu/data/protein_signature_research")
    )
    parser.add_argument("--operation", choices=OPERATIONS)
    parser.add_argument(
        "--parameters", type=Path, help="JSON object containing named scientific controls."
    )
    parser.add_argument("--family-fasta", type=Path)
    parser.add_argument("--seed-fasta", type=Path)
    parser.add_argument("--family-job")
    parser.add_argument("--family-artifact")
    parser.add_argument("--seed-job")
    parser.add_argument("--seed-artifact")
    parser.add_argument("--features-job")
    parser.add_argument("--features-artifact", default="query_features.tsv")
    parser.add_argument("--job-id")
    parser.add_argument(
        "--once", action="store_true", help="Run at most one waiting job and return."
    )
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO"
    )
    return parser


def main(*, argv: Sequence[str] | None = None) -> int:
    """Run diagnostics, queue work or execute the persistent single-workspace worker.

    Args:
        argv: Optional explicit arguments, excluding executable name.

    Returns:
        Zero on success or two on expected validation failure.
    """
    arguments = build_parser().parse_args(args=argv)
    logging.basicConfig(level=arguments.log_level, format="%(asctime)s %(levelname)s %(message)s")
    try:
        if arguments.action == "configure-denbi":
            path = configure_denbi(
                source=arguments.source_config,
                destination=arguments.config,
                workspace=arguments.workspace,
            )
            print(path)
            return 0
        config = load_config(path=arguments.config)
        if arguments.action == "doctor":
            print("capability\tstatus\tdetail")
            for row in diagnostics(config=config):
                print("\t".join(row.values()))
        elif arguments.action == "status":
            print("job_id\toperation\tstatus\tcreated_at\terror")
            for job in list_jobs(config=config):
                print(
                    "\t".join(
                        str(job.get(name, "")).replace("\n", " ").replace("\t", " ")
                        for name in ("job_id", "operation", "status", "created_at", "error")
                    )
                )
        elif arguments.action == "cancel":
            cancel_job(config=config, job_id=arguments.job_id)
        elif arguments.action == "recover":
            recover_interrupted_job(config=config, job_id=arguments.job_id)
        elif arguments.action == "submit":
            parameters = (
                json.loads(s=arguments.parameters.read_text(encoding="utf-8"))
                if arguments.parameters
                else {}
            )
            uploads, sources = {}, {}
            for name, path, job, artifact in (
                ("family", arguments.family_fasta, arguments.family_job, arguments.family_artifact),
                ("seeds", arguments.seed_fasta, arguments.seed_job, arguments.seed_artifact),
            ):
                if path:
                    if path.stat().st_size > 20 * 1024 * 1024:
                        raise ValueError(
                            (
                                "CLI uploads are limited to 20 MiB; reuse a completed job for "
                                "larger "
                                "input."
                            )
                        )
                    uploads[name] = path.read_bytes()
                if job or artifact:
                    if not job or not artifact:
                        raise ValueError(
                            "Both source-job and source-artifact options are required."
                        )
                    sources[name] = {"job_id": job, "name": artifact}
            if arguments.features_job:
                sources["features"] = {
                    "job_id": arguments.features_job,
                    "name": arguments.features_artifact,
                }
            print(
                submit_job(
                    config=config,
                    operation=arguments.operation,
                    parameters=parameters,
                    uploads=uploads,
                    sources=sources,
                )
            )
        elif arguments.action == "worker":
            previous_handler = signal.signal(signal.SIGTERM, signal.default_int_handler)
            try:
                while True:
                    claimed = run_next(config=config)
                    if arguments.once:
                        if claimed:
                            print(claimed)
                        return 0
                    if claimed is None:
                        time.sleep(2)
            finally:
                signal.signal(signal.SIGTERM, previous_handler)
        return 0
    except (ProteinSignatureError, OSError, ValueError, TypeError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        logging.getLogger(__name__).info(
            "Worker interrupted; running jobs require explicit inspection."
        )
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
