"""Start a specific durable request using the app host's installed interpreter."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path
from threading import Thread

from protein_signatures.errors import InputValidationError

from .config import WorkbenchConfig, load_config
from .jobs import job_directory, read_job

LOGGER = logging.getLogger(__name__)


def start_local_job(*, config: WorkbenchConfig, config_path: Path, job_id: str) -> int:
    """Launch an exact queued job in a detached child and return the child PID.

    Args:
        config: Profile used to submit the job.
        config_path: Trusted, launch-time TOML path on the app host.
        job_id: Selected waiting request. The same worker lock applies.

    Returns:
        Child process ID. The job remains queued if the workspace is busy.
    """
    path = config_path.resolve(strict=True)
    if load_config(path=path) != config:
        raise InputValidationError("The launch profile changed; reload the workbench first.")
    if read_job(config=config, job_id=job_id)["status"] != "QUEUED":
        raise InputValidationError("Only a queued job can be started on the app host.")
    log = job_directory(config=config, job_id=job_id) / "app_worker.log"
    descriptor = os.open(log, os.O_WRONLY | os.O_CREAT | os.O_APPEND | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, mode="ab") as handle:
        child = subprocess.Popen(
            args=[
                sys.executable,
                "-m",
                "protein_signatures.workbench.cli",
                "--config",
                str(path),
                "--action",
                "worker",
                "--once",
                "--job-id",
                job_id,
            ],
            stdin=subprocess.DEVNULL,
            stdout=handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    Thread(target=child.wait, daemon=True, name=f"reap-{job_id}").start()
    LOGGER.info("Started app-host worker %d for %s", child.pid, job_id)
    return child.pid
