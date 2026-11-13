"""Transfer a bounded family to an operator-configured SSH subclustering worker."""

from __future__ import annotations

import fcntl
import json
import logging
import re
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO
from zipfile import ZIP_DEFLATED, BadZipFile, ZipFile

from protein_signatures import __version__
from protein_signatures.checksums import sha256_file
from protein_signatures.errors import InputValidationError

from .analysis import iter_fasta, run_command
from .config import WorkbenchConfig, validate_subcluster_ssh
from .storage import temporary_root
from .subclustering import subcluster_controls, subcluster_local

LOGGER = logging.getLogger(__name__)
REQUEST_LIMIT = 64 * 1024 * 1024
RESPONSE_LIMIT = 192 * 1024 * 1024
EXPANDED_LIMIT = 512 * 1024 * 1024
BLOCK = 1024 * 1024
MANIFEST = "subcluster_response.json"
PROTOCOL = "protein-signatures-subcluster-v1"
REQUIRED = frozenset(
    {
        "sequences.fasta",
        "subcluster_id_map.tsv",
        "native_membership.tsv",
        "exact_alignments.tsv",
        "membership.parquet",
        "membership.tsv",
        "subclusters.tsv",
        "subcluster_settings.json",
        "diamond_version.txt",
    }
)
OPTIONAL = frozenset({"commands.jsonl", "tool_stdout.log", "tool_stderr.log"})


def _write_request(*, source: Path, request: Path, controls: dict) -> dict:
    """Write finite metadata and exact FASTA bytes for a request-bound response."""
    size = source.stat().st_size
    if not 1 <= size <= REQUEST_LIMIT:
        raise InputValidationError("SSH family FASTA must contain 1 byte to 64 MiB.")
    record = {
        "protocol": PROTOCOL,
        "package_version": __version__,
        "input_bytes": size,
        "input_sha256": sha256_file(path=source),
        "controls": controls,
    }
    with request.open(mode="xb") as stream, source.open(mode="rb") as input_stream:
        stream.write(json.dumps(obj=record, sort_keys=True).encode() + b"\n")
        shutil.copyfileobj(input_stream, stream, length=BLOCK)
    return record


def _read_request(*, stream: BinaryIO, destination: Path) -> dict:
    """Read exactly the authorised input length and verify its content hash."""
    raw = stream.readline(8193)
    if len(raw) > 8192 or not raw.endswith(b"\n"):
        raise InputValidationError("Invalid SSH subcluster request header.")
    try:
        record = json.loads(s=raw)
    except (ValueError, UnicodeError) as error:
        raise InputValidationError("Invalid SSH subcluster request header.") from error
    if (
        not isinstance(record, dict)
        or set(record) != {"protocol", "package_version", "input_bytes", "input_sha256", "controls"}
        or record["protocol"] != PROTOCOL
        or record["package_version"] != __version__
        or type(record["input_bytes"]) is not int
        or not 1 <= record["input_bytes"] <= REQUEST_LIMIT
        or not isinstance(record["input_sha256"], str)
        or not re.fullmatch(r"[a-f0-9]{64}", record["input_sha256"])
        or not isinstance(record["controls"], dict)
        or set(record["controls"]) != {"identity", "mutual_coverage", "evalue"}
    ):
        raise InputValidationError("SSH subcluster request is invalid or uses another app version.")
    record["controls"] = subcluster_controls(parameters=record["controls"])
    remaining = record["input_bytes"]
    with destination.open(mode="xb") as output:
        while remaining:
            block = stream.read(min(remaining, BLOCK))
            if not block:
                raise InputValidationError("SSH subcluster FASTA was truncated.")
            output.write(block)
            remaining -= len(block)
    if stream.read(1) or sha256_file(path=destination) != record["input_sha256"]:
        raise InputValidationError("SSH subcluster input failed its checksum or length.")
    return record


def _write_response(*, output: Path, archive: Path, request: dict) -> None:
    """Build an allowlisted results ZIP with exact sizes and hashes."""
    files = {}
    for name in sorted(REQUIRED | OPTIONAL):
        path = output / name
        if path.is_file():
            files[name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path=path)}
    if not REQUIRED <= files.keys() or sum(row["bytes"] for row in files.values()) > EXPANDED_LIMIT:
        raise InputValidationError(
            "SSH subcluster response is incomplete or exceeds the size limit."
        )
    with ZipFile(file=archive, mode="x", compression=ZIP_DEFLATED) as bundle:
        for name in files:
            bundle.write(filename=output / name, arcname=name)
        bundle.writestr(
            zinfo_or_arcname=MANIFEST,
            data=json.dumps(
                obj={
                    "request": request,
                    "files": files,
                },
                sort_keys=True,
            ),
        )
    if archive.stat().st_size > RESPONSE_LIMIT:
        raise InputValidationError("SSH subcluster archive exceeds the response limit.")


def write_stream_response(
    *,
    config: WorkbenchConfig,
    source: BinaryIO,
    destination: BinaryIO,
) -> None:
    """Serve one authenticated SSH request under the remote workspace's worker lock."""
    workspace = Path(config.workspace)
    workspace.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (workspace / "worker.lock").open(mode="a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise InputValidationError("The configured remote workspace is busy.") from error
        with TemporaryDirectory(
            dir=temporary_root(config=config), prefix="ssh_subcluster_"
        ) as root:
            path = Path(root)
            fasta, output, archive = path / "input.fasta", path / "results", path / "response.zip"
            request = _read_request(stream=source, destination=fasta)
            output.mkdir(mode=0o700)
            subcluster_local(
                config=config, source=fasta, output=output, controls=request["controls"]
            )
            _write_response(output=output, archive=archive, request=request)
            with archive.open(mode="rb") as payload:
                shutil.copyfileobj(payload, destination, length=BLOCK)
            destination.flush()


def _read_response(*, archive: Path, output: Path, request: dict) -> None:
    """Check inventory, request, sizes and hashes before using any remote result."""
    try:
        with ZipFile(file=archive) as bundle:
            members = bundle.infolist()
            names = {item.filename for item in members}
            if (
                len(names) != len(members)
                or not (REQUIRED | {MANIFEST}) <= names <= (REQUIRED | OPTIONAL | {MANIFEST})
                or any(item.flag_bits & 1 for item in members)
                or sum(item.file_size for item in members) > EXPANDED_LIMIT + 32_768
                or bundle.getinfo(MANIFEST).file_size > 32_768
            ):
                raise InputValidationError("SSH subcluster response inventory is invalid.")
            manifest = json.loads(s=bundle.read(name=MANIFEST))
            if not isinstance(manifest, dict) or manifest.get("request") != request:
                raise InputValidationError("SSH subcluster response belongs to another request.")
            files = manifest.get("files")
            if not isinstance(files, dict) or set(files) != names - {MANIFEST}:
                raise InputValidationError("SSH subcluster response manifest is incomplete.")
            for name, record in files.items():
                info = bundle.getinfo(name)
                if (
                    not isinstance(record, dict)
                    or record.get("bytes") != info.file_size
                    or not isinstance(record.get("sha256"), str)
                ):
                    raise InputValidationError(
                        "SSH subcluster response size or checksum is invalid."
                    )
                path = output / (f"remote_{name}" if name in OPTIONAL else name)
                with bundle.open(name) as payload, path.open(mode="xb") as handle:
                    shutil.copyfileobj(payload, handle, length=BLOCK)
                if sha256_file(path=path) != record["sha256"]:
                    raise InputValidationError("SSH subcluster response failed its checksum.")
    except (BadZipFile, ValueError, RuntimeError, NotImplementedError) as error:
        raise InputValidationError("Could not read the verified SSH subcluster archive.") from error


def subcluster_remote(
    *,
    config: WorkbenchConfig,
    source: Path,
    output: Path,
    controls: dict,
) -> None:
    """Send a finite request to a trusted host and verify the exact returned proteins."""
    validate_subcluster_ssh(config=config)
    if not 1 <= source.stat().st_size <= REQUEST_LIMIT:
        raise InputValidationError("SSH family FASTA must contain 1 byte to 64 MiB.")
    expected = {}
    for identifier, sequence in iter_fasta(path=source):
        if identifier in expected or len(expected) >= config.max_subcluster_sequences:
            raise InputValidationError("SSH family contains duplicate IDs or exceeds the limit.")
        expected[identifier] = sequence
    if sum(map(len, expected.values())) > config.max_subcluster_residues:
        raise InputValidationError("SSH family exceeds the configured residue limit.")
    with TemporaryDirectory(
        dir=temporary_root(config=config), prefix="subcluster_request_"
    ) as root:
        path = Path(root)
        request, archive = path / "request.bin", path / "response.zip"
        record = _write_request(source=source, request=request, controls=controls)
        run_command(
            config=config,
            output=output,
            stdin_path=request,
            stdout_path=archive,
            stdout_limit=RESPONSE_LIMIT,
            arguments=[
                "ssh",
                "-T",
                "-o",
                "BatchMode=yes",
                "-o",
                "ClearAllForwardings=yes",
                "-o",
                "PermitLocalCommand=no",
                "-o",
                "ConnectTimeout=15",
                "-o",
                "StrictHostKeyChecking=yes",
                config.subcluster_ssh_host,
                config.subcluster_ssh_executable,
                "--config",
                config.subcluster_ssh_profile,
                "--action",
                "subcluster-stream",
            ],
        )
        _read_response(archive=archive, output=output, request=record)
    returned = list(iter_fasta(path=output / "sequences.fasta"))
    if len(returned) != len(expected) or dict(returned) != expected:
        raise InputValidationError(
            "SSH response changed the submitted protein identifiers or sequences."
        )
    LOGGER.info("Verified SSH subcluster response for %d proteins", len(expected))
