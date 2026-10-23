"""Streaming sequence preparation and bounded, audited external-tool execution."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import os
import signal
import subprocess
from collections import Counter
from collections.abc import Iterator, Sequence
from pathlib import Path

import pyarrow as pa

from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.io_utils import write_json_atomic

from .catalogue import identifiers, normalise_sequence, open_catalogue
from .config import WorkbenchConfig, bounded_number

LOGGER = logging.getLogger(__name__)
HIT_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)


def iter_fasta(*, path: Path) -> Iterator[tuple[str, str]]:
    """Yield unaligned proteins without retaining an in-memory identifier catalogue.

    Args:
        path: UTF-8 protein FASTA; uniqueness is checked by index_fasta.

    Yields:
        Exact ID and normalised protein sequence.

    Raises:
        InputValidationError: If FASTA syntax or sequences are invalid.
    """
    current, parts = None, []
    with path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                if current is not None:
                    yield current, normalise_sequence(value="".join(parts))
                header = value[1:].split()
                if not header:
                    raise InputValidationError("FASTA contains an empty header.")
                current = identifiers(values=[header[0]])[0]
                parts = []
            elif current is None:
                raise InputValidationError("Sequence text precedes the first FASTA header.")
            else:
                parts.append(value)
    if current is None:
        raise InputValidationError("FASTA contains no protein records.")
    yield current, normalise_sequence(value="".join(parts))


def index_fasta(*, connection: object, path: Path, name: str, maximum: int) -> int:
    """Index FASTA in bounded Arrow batches with SQL uniqueness checks.

    Args:
        connection: Private DuckDB connection with disk spill.
        path: Source FASTA.
        name: Internal fixed table name, family or seeds.
        maximum: Maximum complete protein count; no silent truncation.

    Returns:
        Complete protein count.

    Raises:
        InputValidationError: If IDs repeat, proteins are invalid or limits are exceeded.
    """
    if name not in {"family", "seeds"}:
        raise InputValidationError("Unsupported internal sequence table.")
    connection.execute(
        f"CREATE TABLE {name} (sequence_id VARCHAR,sequence VARCHAR,sequence_sha256 VARCHAR)"
    )
    batch, count = [], 0
    for identifier, protein in iter_fasta(path=path):
        count += 1
        if count > maximum or len(protein) > 100_000:
            raise InputValidationError(
                "Complete protein input exceeds configured count or length bounds."
            )
        batch.append(
            {
                "sequence_id": identifier,
                "sequence": protein,
                "sequence_sha256": hashlib.sha256(protein.encode()).hexdigest(),
            }
        )
        if len(batch) == 4096:
            connection.register(
                view_name="input_batch", python_object=pa.Table.from_pylist(mapping=batch)
            )
            connection.execute(f"INSERT INTO {name} SELECT * FROM input_batch")
            connection.unregister(view_name="input_batch")
            batch.clear()
    if batch:
        connection.register(
            view_name="input_batch", python_object=pa.Table.from_pylist(mapping=batch)
        )
        connection.execute(f"INSERT INTO {name} SELECT * FROM input_batch")
        connection.unregister(view_name="input_batch")
    if connection.execute(
        f"SELECT 1 FROM {name} GROUP BY sequence_id HAVING count(*)>1 LIMIT 1"
    ).fetchone():
        raise InputValidationError("FASTA contains duplicate identifiers.")
    return count


def write_fasta(*, connection: object, query: str, path: Path) -> int:
    """Stream an internal two-column query into a FASTA file.

    Args:
        connection: Private analytical connection.
        query: Trusted sequence-ID and sequence SELECT.
        path: New destination.

    Returns:
        Written protein count.
    """
    cursor = connection.execute(query)
    count = 0
    with path.open(mode="x", encoding="utf-8") as handle:
        while batch := cursor.fetchmany(size=4096):
            for identifier, sequence in batch:
                handle.write(f">{identifier}\n{sequence}\n")
                count += 1
    return count


def run_command(
    *,
    config: WorkbenchConfig,
    arguments: Sequence[str],
    output: Path,
    stdout_name: str | None = None,
) -> None:
    """Run an argument vector, retaining logs and terminating its process group on timeout.

    Args:
        config: Command timeout and workspace limits.
        arguments: Explicit executable and arguments; no shell interpolation.
        output: Owned output directory.
        stdout_name: Optional output receiving stdout, for MAFFT and FastTree.

    Raises:
        ExternalToolError: If a command cannot run, fails or times out.
    """
    LOGGER.info("Running %s", arguments[0])
    with (output / "commands.jsonl").open(mode="a", encoding="utf-8") as audit:
        audit.write(
            json.dumps(
                obj={
                    "argv": list(arguments),
                    "tool_version": executable_version(config=config, executable=arguments[0]),
                }
            )
            + "\n"
        )
    environment = dict(os.environ)
    environment["TMPDIR"] = str(Path(config.workspace) / "temporary")
    with (
        (output / (stdout_name or "tool_stdout.log")).open(mode="ab") as stdout,
        (output / "tool_stderr.log").open(mode="ab") as stderr,
    ):
        try:
            process = subprocess.Popen(
                args=list(arguments),
                stdout=stdout,
                stderr=stderr,
                env=environment,
                start_new_session=True,
            )
        except OSError as error:
            raise ExternalToolError(f"Executable unavailable: {arguments[0]}") from error
        marker = output / "active_process.json"
        try:
            write_json_atomic(
                path=marker, value={"process_group": process.pid, "argv": list(arguments)}
            )
            status = process.wait(timeout=config.command_timeout)
        except subprocess.TimeoutExpired as error:
            terminate_process(process=process)
            raise ExternalToolError(f"Command timed out: {arguments[0]}") from error
        except BaseException:
            # An interrupted worker must not leave its separate tool group running.
            terminate_process(process=process)
            raise
        finally:
            marker.unlink(missing_ok=True)
    if status != 0:
        raise ExternalToolError(
            f"{arguments[0]} exited with status {status}; inspect tool_stderr.log."
        )


def terminate_process(*, process: subprocess.Popen) -> None:
    """Stop an owned external process group and reap the child with a bounded grace period.

    Args:
        process: Child started in its own process session.
    """
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait()
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def executable_version(*, config: WorkbenchConfig, executable: str) -> dict:
    """Record a bounded version probe for a configured analytical executable.

    Args:
        config: Trusted executable names and timeout.
        executable: Configured tool being invoked.

    Returns:
        Probe arguments, exit status and up to 2000 output characters, or explicit unavailability.
    """
    flags = {
        config.diamond: ["version"],
        config.mafft: ["--version"],
        config.fasttree: ["-help"],
        config.hmmscan: ["-h"],
    }
    if executable not in flags:
        return {"status": "NOT_CONFIGURED_ANALYTICAL_TOOL"}
    arguments = [executable, *flags[executable]]
    try:
        result = subprocess.run(
            args=arguments,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="utf-8",
            errors="replace",
            check=False,
            timeout=min(5, config.command_timeout),
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"argv": arguments, "status": "UNAVAILABLE", "detail": str(error)[:500]}
    return {"argv": arguments, "status": result.returncode, "output": result.stdout[:2000]}


def thresholds(*, parameters: dict) -> tuple[float, float, float]:
    """Validate independent percent-identity, bilateral-coverage and E-value thresholds.

    Args:
        parameters: User scientific controls.

    Returns:
        Identity percent, coverage percent and maximum E-value.
    """
    identity = bounded_number(
        value=parameters.get("identity", 30), name="Identity", minimum=0, maximum=100
    )
    coverage = bounded_number(
        value=parameters.get("coverage", 50), name="Bilateral coverage", minimum=0, maximum=100
    )
    evalue = bounded_number(
        value=parameters.get("evalue", 1e-5), name="E-value", minimum=1e-300, maximum=1
    )
    return float(identity), float(coverage), float(evalue)


def search_sequences(
    *,
    config: WorkbenchConfig,
    query: Path,
    database: str,
    parameters: dict,
    output: Path,
    maximum_hits: int,
) -> Path:
    """Run DIAMOND with explicit controls and validate every retained alignment row.

    Args:
        config: Installed DIAMOND and worker bounds.
        query: Normalised query FASTA.
        database: Local DIAMOND database.
        parameters: Identity, coverage and E-value thresholds.
        output: Owned job output.
        maximum_hits: Maximum targets per query; this search is explicitly bounded.

    Returns:
        Headered, validated alignment TSV with both measured coverages.
    """
    identity, coverage, evalue = thresholds(parameters=parameters)
    raw = output / "hits.raw.tsv"
    run_command(
        config=config,
        output=output,
        arguments=[
            config.diamond,
            "blastp",
            "--query",
            str(query),
            "--db",
            database,
            "--out",
            str(raw),
            "--outfmt",
            "6",
            *HIT_FIELDS,
            "--threads",
            str(config.threads),
            "--tmpdir",
            str(Path(config.workspace) / "temporary"),
            "--block-size",
            str(min(0.4, config.memory_gb / 12)),
            "--index-chunks",
            "4",
            "--max-target-seqs",
            str(maximum_hits),
            "--evalue",
            str(evalue),
            "--id",
            str(identity),
            "--query-cover",
            str(coverage),
            "--subject-cover",
            str(coverage),
            "--very-sensitive",
            "--no-parse-seqids",
        ],
    )
    result = output / "hits.tsv"
    with (
        raw.open(encoding="utf-8") as source,
        result.open(mode="x", encoding="utf-8", newline="") as destination,
    ):
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow([*HIT_FIELDS, "query_coverage_percent", "subject_coverage_percent"])
        for fields in csv.reader(source, delimiter="\t"):
            if len(fields) != len(HIT_FIELDS):
                raise InputValidationError("DIAMOND returned an invalid alignment record.")
            row = dict(zip(HIT_FIELDS, fields, strict=True))
            try:
                numbers = {name: float(row[name]) for name in HIT_FIELDS[2:]}
                qcov = 100 * (abs(numbers["qend"] - numbers["qstart"]) + 1) / numbers["qlen"]
                scov = 100 * (abs(numbers["send"] - numbers["sstart"]) + 1) / numbers["slen"]
            except (ValueError, ZeroDivisionError) as error:
                raise InputValidationError("DIAMOND returned invalid numeric evidence.") from error
            if (
                any(not math.isfinite(value) for value in numbers.values())
                or not 0 <= numbers["pident"] <= 100
                or not 0 <= numbers["evalue"]
                or not 0 < qcov <= 100
                or not 0 < scov <= 100
                or not 1 <= numbers["qstart"] <= numbers["qend"] <= numbers["qlen"]
                or not 1 <= numbers["sstart"] <= numbers["send"] <= numbers["slen"]
                or any(
                    not numbers[name].is_integer()
                    for name in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send")
                )
            ):
                raise InputValidationError("DIAMOND returned non-finite or out-of-range evidence.")
            identifiers(values=[row["qseqid"], row["sseqid"]])
            if (
                numbers["pident"] >= identity
                and qcov >= coverage
                and scov >= coverage
                and numbers["evalue"] <= evalue
            ):
                writer.writerow([*fields, qcov, scov])
    return result


def align_family(*, config: WorkbenchConfig, source: Path, parameters: dict, output: Path) -> None:
    """Align a bounded explicit selection and optionally infer a WAG+gamma tree.

    Args:
        config: MAFFT, FastTree and alignment resource limits.
        source: Unaligned source proteins.
        parameters: Optional exact selected IDs and tree Boolean.
        output: Owned job directory.

    Raises:
        InputValidationError: If selection, alignment or output contracts fail.
    """
    with open_catalogue(config=config, with_index=False) as connection:
        index_fasta(connection=connection, path=source, name="family", maximum=config.max_sequences)
        selected = parameters.get("selected_ids", [])
        if selected:
            values = identifiers(values=selected, maximum=config.max_alignment_sequences)
            connection.execute(
                "CREATE TABLE selected AS SELECT unnest(?) AS sequence_id", [list(values)]
            )
            if connection.execute(
                "SELECT sequence_id FROM selected EXCEPT SELECT sequence_id FROM family"
            ).fetchone():
                raise InputValidationError("An alignment selection ID is absent from the input.")
            connection.execute(
                (
                    "CREATE VIEW alignment_input AS SELECT f.* FROM family f JOIN selected "
                    "USING(sequence_id)"
                )
            )
        else:
            connection.execute("CREATE VIEW alignment_input AS SELECT * FROM family")
        count, max_length = connection.execute(
            "SELECT count(*),max(length(sequence)) FROM alignment_input"
        ).fetchone()
        if (
            count < 2
            or count > config.max_alignment_sequences
            or count * max_length > config.max_alignment_cells
        ):
            raise InputValidationError(
                (
                    "Alignment exceeds configured sequence/cell limits or has fewer than "
                    "two sequences. Select a smaller explicit subset."
                )
            )
        # Short internal names prevent Newick/MAFFT identifier truncation or punctuation loss.
        rows = connection.execute(
            "SELECT sequence_id,sequence FROM alignment_input ORDER BY sequence_id"
        ).fetchall()
    mapping = {f"seq_{i:06d}": identifier for i, (identifier, _) in enumerate(rows)}
    input_path = output / "tool_input.fasta"
    input_path.write_text(
        data="".join(f">seq_{i:06d}\n{sequence}\n" for i, (_, sequence) in enumerate(rows)),
        encoding="utf-8",
    )
    with (output / "tip_labels.tsv").open(mode="x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tree_tip", "sequence_id"])
        writer.writerows(mapping.items())
    run_command(
        config=config,
        output=output,
        stdout_name="alignment.internal.fasta",
        arguments=[
            config.mafft,
            "--thread",
            str(config.threads),
            "--auto",
            str(input_path),
        ],
    )
    alignment, parts, current = {}, [], None
    with (output / "alignment.internal.fasta").open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if current is not None:
                    if current in alignment:
                        raise InputValidationError("MAFFT duplicated an identifier.")
                    alignment[current] = "".join(parts).upper()
                current, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
        if current is not None:
            if current in alignment:
                raise InputValidationError("MAFFT duplicated an identifier.")
            alignment[current] = "".join(parts).upper()
    lengths = {len(value) for value in alignment.values()}
    if set(alignment) != set(mapping) or len(lengths) != 1 or not next(iter(lengths), 0):
        raise InputValidationError("MAFFT returned missing IDs or an invalid alignment matrix.")
    width = next(iter(lengths))
    if count * width > config.max_alignment_cells:
        raise InputValidationError("Aligned matrix exceeds max_alignment_cells.")
    for i, (_, sequence) in enumerate(rows):
        if alignment[f"seq_{i:06d}"].replace("-", "") != sequence:
            raise InputValidationError("MAFFT altered an input protein sequence.")
    (output / "alignment.fasta").write_text(
        data="".join(f">{mapping[key]}\n{sequence}\n" for key, sequence in alignment.items()),
        encoding="utf-8",
    )
    (output / "sequences.fasta").write_text(
        data="".join(f">{key}\n{sequence}\n" for key, sequence in rows), encoding="utf-8"
    )
    with (output / "conservation.tsv").open(mode="x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "alignment_position",
                "consensus",
                "consensus_fraction_non_gap",
                "gap_fraction",
                "shannon_entropy_bits",
            ]
        )
        for position in range(width):
            residues = Counter(sequence[position] for sequence in alignment.values())
            gaps = residues.pop("-", 0)
            denominator = sum(residues.values())
            consensus, frequency = (
                sorted(residues.items(), key=lambda item: (-item[1], item[0]))[0]
                if residues
                else ("-", 0)
            )
            entropy = (
                -sum(
                    (value / denominator) * math.log2(value / denominator)
                    for value in residues.values()
                )
                if denominator
                else 0
            )
            writer.writerow(
                [
                    position + 1,
                    consensus,
                    frequency / denominator if denominator else 0,
                    gaps / count,
                    entropy,
                ]
            )
    tree = parameters.get("tree", False)
    if not isinstance(tree, bool):
        raise InputValidationError("tree must be a Boolean.")
    if tree:
        if count < 3:
            raise InputValidationError("Tree inference requires at least three sequences.")
        run_command(
            config=config,
            output=output,
            stdout_name="tree.nwk",
            arguments=[
                config.fasttree,
                "-wag",
                "-gamma",
                str(output / "alignment.internal.fasta"),
            ],
        )
        if not (output / "tree.nwk").read_text(encoding="utf-8").strip().endswith(";"):
            raise InputValidationError("FastTree did not return a complete Newick tree.")
