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

from .config import configure_denbi, configure_mount, diagnostics, load_config
from .jobs import (
    MAX_UPLOAD_BYTES,
    OPERATIONS,
    cancel_job,
    list_jobs,
    recover_interrupted_job,
    run_next,
    submit_job,
)
from .structures import compare_previous, old_pairs


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
        choices=(
            "configure-denbi",
            "configure-mount",
            "doctor",
            "submit",
            "worker",
            "status",
            "cancel",
            "recover",
            "workspace",
            "plan-structures",
            "compare-structures",
            "subcluster-stream",
        ),
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
    parser.add_argument("--alignment-fasta", type=Path)
    parser.add_argument("--alignment-job")
    parser.add_argument("--alignment-artifact", default="alignment.fasta")
    parser.add_argument("--context-tsv", type=Path)
    parser.add_argument("--ebi-contact-file", type=Path)
    parser.add_argument("--subcluster-id")
    parser.add_argument("--seed-fasta", type=Path)
    parser.add_argument("--query-model", type=Path)
    parser.add_argument("--target-model", type=Path)
    parser.add_argument("--deepclust-export", type=Path)
    parser.add_argument("--deepclust-root", type=Path)
    parser.add_argument("--centroid-database", type=Path)
    parser.add_argument("--resource", type=Path)
    parser.add_argument("--release")
    parser.add_argument("--pair-table", type=Path)
    parser.add_argument("--models-root", type=Path)
    parser.add_argument("--source-model-root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cluster-id", action="append", default=[])
    parser.add_argument("--max-pairs", type=int, default=100)
    parser.add_argument("--single-domain", action="store_true")
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
        if arguments.action == "configure-mount":
            if not arguments.deepclust_root or not arguments.release:
                raise ValueError(
                    "Supply --deepclust-root and --release for the installed database."
                )
            path = configure_mount(
                root=arguments.deepclust_root,
                destination=arguments.config,
                workspace=arguments.workspace,
                release=arguments.release,
                resource=arguments.resource,
                centroid_database=arguments.centroid_database,
            )
            print(path)
            return 0
        config = load_config(path=arguments.config)
        if arguments.action == "workspace":
            print(config.workspace)
        elif arguments.action == "doctor":
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
        elif arguments.action in {"plan-structures", "compare-structures"}:
            if not arguments.pair_table:
                raise ValueError(
                    "Supply the original structural_alignments.parquet as --pair-table."
                )
            if arguments.action == "plan-structures":
                rows = old_pairs(
                    table=arguments.pair_table,
                    clusters=arguments.cluster_id,
                    maximum=arguments.max_pairs,
                )
                print("metric\tcount")
                print(f"original_rows\t{len(rows)}")
                print(f"groups\t{len({row['cluster_id'] for row in rows})}")
                print(
                    "Model paths and SHA-256 will be checked before each comparison.",
                    file=sys.stderr,
                )
            else:
                if not arguments.models_root or not arguments.output:
                    raise ValueError("Supply --models-root and a new --output directory.")
                report = compare_previous(
                    config=config,
                    table=arguments.pair_table,
                    models_root=arguments.models_root,
                    source_model_root=arguments.source_model_root
                    or arguments.models_root.resolve(),
                    output=arguments.output,
                    clusters=arguments.cluster_id,
                    maximum=arguments.max_pairs,
                    multidomain=not arguments.single_domain,
                )
                print(report)
        elif arguments.action == "submit":
            parameters = (
                json.loads(s=arguments.parameters.read_text(encoding="utf-8"))
                if arguments.parameters
                else {}
            )
            uploads, sources = {}, {}
            for name, path in (
                ("query_model", arguments.query_model),
                ("target_model", arguments.target_model),
                ("deepclust_export", arguments.deepclust_export),
                ("alignment", arguments.alignment_fasta),
                ("context", arguments.context_tsv),
                ("ebi_email", arguments.ebi_contact_file),
            ):
                if path:
                    if path.stat().st_size > MAX_UPLOAD_BYTES:
                        raise ValueError("CLI uploads are limited to 20 MiB per input.")
                    uploads[name] = path.read_bytes()
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
            if arguments.alignment_job:
                sources["alignment"] = {
                    "job_id": arguments.alignment_job,
                    "name": arguments.alignment_artifact,
                }
            if arguments.subcluster_id:
                sources["subcluster_membership"] = {
                    "job_id": arguments.family_job,
                    "name": "membership.parquet",
                    "subcluster_id": arguments.subcluster_id,
                }
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
        elif arguments.action == "subcluster-stream":
            from .subcluster_remote import write_stream_response

            previous_handler = signal.signal(signal.SIGTERM, signal.default_int_handler)
            previous_hup = signal.signal(signal.SIGHUP, signal.default_int_handler)
            try:
                write_stream_response(
                    config=config,
                    source=sys.stdin.buffer,
                    destination=sys.stdout.buffer,
                )
            finally:
                signal.signal(signal.SIGTERM, previous_handler)
                signal.signal(signal.SIGHUP, previous_hup)
        elif arguments.action == "worker":
            if arguments.job_id and not arguments.once:
                raise ValueError("Use --once with an exact --job-id.")
            previous_handler = signal.signal(signal.SIGTERM, signal.default_int_handler)
            try:
                while True:
                    claimed = run_next(config=config, job_id=arguments.job_id)
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
