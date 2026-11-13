"""Conservative DIAMOND DeepClust partitions audited by exact realignment."""

from __future__ import annotations

import csv
import hashlib
import logging
import math
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

import pyarrow as pa
import pyarrow.parquet as pq

from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic, write_tsv_atomic

from .analysis import index_fasta, run_command, write_fasta
from .catalogue import open_catalogue
from .config import WorkbenchConfig, bounded_number
from .storage import temporary_root

LOGGER = logging.getLogger(__name__)
REALIGN_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "length",
    "evalue",
    "bitscore",
    "nident",
)
PASS = "EXACT_THRESHOLDS_PASS"
CENTROID = "CENTROID_SELF"
MEMBERSHIP_FIELDS = (
    "subcluster_id",
    "audited_representative_id",
    "sequence_id",
    "native_representative_id",
    "alignment_status",
    "sequence_length",
    "pident",
    "representative_coverage",
    "member_coverage",
    "evalue",
    "bitscore",
    "alignment_length",
    "identical_count",
)
NATIVE_HEADERS = {("centroid", "member"), ("representative", "member"), ("cluster", "member")}


def subcluster_controls(*, parameters: dict) -> dict[str, float]:
    """Validate identity, mutual coverage and E-value cutoffs in user-facing units."""
    return {
        "identity": float(
            bounded_number(
                value=parameters.get("identity", 70),
                name="Identity percent",
                minimum=1,
                maximum=100,
            )
        ),
        "mutual_coverage": float(
            bounded_number(
                value=parameters.get("mutual_coverage", 80),
                name="Mutual coverage percent",
                minimum=1,
                maximum=100,
            )
        ),
        "evalue": float(
            bounded_number(
                value=parameters.get("evalue", 1e-10),
                name="E-value",
                minimum=1e-300,
                maximum=1,
            )
        ),
    }


def native_partition(*, path: Path, aliases: dict[str, tuple[str, int]]) -> dict[str, str]:
    """Read a complete, disjoint native assignment with a self row for each centroid."""
    mapping = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        if tuple(next(reader, [])) not in NATIVE_HEADERS:
            raise InputValidationError("Unrecognised DIAMOND DeepClust membership header.")
        for row in reader:
            if len(row) != 2:
                raise InputValidationError("Malformed native DIAMOND membership row.")
            centroid, member = row
            if centroid not in aliases or member not in aliases or member in mapping:
                raise InputValidationError("Native membership contains an unknown or duplicate ID.")
            mapping[member] = centroid
    if set(mapping) != set(aliases) or any(mapping.get(rep) != rep for rep in mapping.values()):
        raise InputValidationError("Native membership must assign every input to a valid centroid.")
    return mapping


def exact_evidence(*, row: list[str], query_length: int, subject_length: int) -> dict | None:
    """Validate a realignment and derive exact identity and both sequence coverages."""
    try:
        values = {
            field: float(value)
            for field, value in zip(
                REALIGN_FIELDS[2:],
                row[2:],
                strict=True,
            )
        }
    except ValueError:
        return None
    if not all(math.isfinite(value) for value in values.values()):
        return None
    integers = ("qlen", "slen", "qstart", "qend", "sstart", "send", "length", "nident")
    qspan, sspan = values["qend"] - values["qstart"] + 1, values["send"] - values["sstart"] + 1
    if (
        any(not values[name].is_integer() for name in integers)
        or values["qlen"] != query_length
        or values["slen"] != subject_length
        or not 1 <= values["qstart"] <= values["qend"] <= query_length
        or not 1 <= values["sstart"] <= values["send"] <= subject_length
        or not max(qspan, sspan) <= values["length"] <= qspan + sspan
        or not 0 <= values["nident"] <= min(qspan, sspan)
        or not 0 <= values["pident"] <= 100
        or values["evalue"] < 0
    ):
        return None
    return {
        "pident": 100 * values["nident"] / values["length"],
        "representative_coverage": 100 * qspan / query_length,
        "member_coverage": 100 * sspan / subject_length,
        "evalue": values["evalue"],
        "bitscore": values["bitscore"],
        "alignment_length": int(values["length"]),
        "identical_count": int(values["nident"]),
    }


def audit_partition(
    *,
    native: dict[str, str],
    aliases: dict[str, tuple[str, int]],
    exact: Path,
    controls: dict,
) -> list[dict]:
    """Keep only passing centroid assignments; rescue missing and failed edges as singletons."""
    best, invalid = {}, set()
    with exact.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        for index, row in enumerate(reader, start=1):
            if index > len(aliases) * 20 or len(row) != len(REALIGN_FIELDS):
                raise InputValidationError("Unexpected size or column count in DIAMOND realign.")
            first, second = row[:2]
            if second not in native or native[second] != first:
                raise InputValidationError("Realignment contains a pair outside native membership.")
            record = exact_evidence(
                row=row,
                query_length=aliases[first][1],
                subject_length=aliases[second][1],
            )
            if record is None:
                invalid.add(second)
                continue
            order = record["bitscore"], -record["evalue"], record["pident"]
            if second not in best or order > best[second][0]:
                best[second] = (order, record)
    rows = []
    for member, centroid in native.items():
        identifier, length = aliases[member]
        native_rep = aliases[centroid][0]
        evidence = best[member][1] if member in best else {}
        if member == centroid:
            status = CENTROID
        elif not evidence:
            status = "INVALID_REALIGNMENT" if member in invalid else "NO_REALIGNMENT"
        elif (
            evidence["pident"] >= controls["identity"]
            and evidence["representative_coverage"] >= controls["mutual_coverage"]
            and evidence["member_coverage"] >= controls["mutual_coverage"]
            and evidence["evalue"] <= controls["evalue"]
            and evidence["bitscore"] > 20
        ):
            status = PASS
        else:
            status = "EXACT_THRESHOLDS_FAIL"
        rep = native_rep if status in {CENTROID, PASS} else identifier
        rows.append(
            {
                "subcluster_id": "sub_" + hashlib.sha256(rep.encode()).hexdigest()[:24],
                "audited_representative_id": rep,
                "sequence_id": identifier,
                "native_representative_id": native_rep,
                "alignment_status": status,
                "sequence_length": length,
                **{key: evidence.get(key) for key in MEMBERSHIP_FIELDS[6:]},
            }
        )
    return sorted(rows, key=lambda row: (row["subcluster_id"], row["sequence_id"]))


def write_partition(*, rows: list[dict], output: Path, controls: dict) -> None:
    """Publish audited membership, partition sizes and threshold provenance."""
    types = [
        *[pa.string()] * 5,
        pa.int64(),
        *[pa.float64()] * 5,
        pa.int64(),
        pa.int64(),
    ]
    schema = pa.schema(fields=list(zip(MEMBERSHIP_FIELDS, types, strict=True)))
    pq.write_table(
        table=pa.Table.from_pylist(mapping=rows, schema=schema),
        where=output / "membership.parquet",
        compression="zstd",
    )
    write_tsv_atomic(path=output / "membership.tsv", fieldnames=MEMBERSHIP_FIELDS, records=rows)
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["subcluster_id"]].append(row)
    summaries = []
    for group, members in grouped.items():
        passed = [row for row in members if row["alignment_status"] == PASS]
        summaries.append(
            {
                "subcluster_id": group,
                "representative_id": members[0]["audited_representative_id"],
                "member_count": len(members),
                "rescued_singleton_count": sum(
                    row["alignment_status"] not in {CENTROID, PASS} for row in members
                ),
                "minimum_passing_identity": min((row["pident"] for row in passed), default=None),
                "minimum_passing_representative_coverage": min(
                    (row["representative_coverage"] for row in passed), default=None
                ),
                "minimum_passing_member_coverage": min(
                    (row["member_coverage"] for row in passed), default=None
                ),
            }
        )
    summaries.sort(key=lambda row: (-row["member_count"], row["subcluster_id"]))
    write_tsv_atomic(
        path=output / "subclusters.tsv",
        fieldnames=tuple(summaries[0]),
        records=summaries,
    )
    write_json_atomic(
        path=output / "subcluster_settings.json",
        value={
            "input_sequences": len(rows),
            "native_subclusters": len({row["native_representative_id"] for row in rows}),
            "audited_subclusters": len(grouped),
            "alignment_status_counts": dict(Counter(row["alignment_status"] for row in rows)),
            "exact_identity_percent": controls["identity"],
            "mutual_coverage_percent": controls["mutual_coverage"],
            "maximum_evalue": controls["evalue"],
            "minimum_bitscore_exclusive": 20,
            "interpretation": "Failed or absent centroid links yield singletons.",
        },
    )


def subcluster_local(
    *,
    config: WorkbenchConfig,
    source: Path,
    output: Path,
    controls: dict,
) -> None:
    """Run DIAMOND 2.2.1+ against a temporary family database and audit all assignments."""
    run_command(
        config=config,
        output=output,
        arguments=[config.diamond, "version"],
        stdout_name="diamond_version.txt",
    )
    version = (output / "diamond_version.txt").read_text(encoding="utf-8")
    found = re.search(r"(?im)^\s*(?:diamond(?:\s+version)?\s+)?(\d+)\.(\d+)\.(\d+)\b", version)
    if not found or tuple(map(int, found.groups())) < (2, 2, 1):
        raise InputValidationError(
            "Audited DeepClust subclustering requires DIAMOND 2.2.1 or newer."
        )
    with TemporaryDirectory(dir=temporary_root(config=config), prefix="subcluster_") as scratch:
        root = Path(scratch)
        tool_input = root / "proteins.fasta"
        with open_catalogue(config=config, with_index=False) as connection:
            index_fasta(
                connection=connection,
                path=source,
                name="family",
                maximum=config.max_subcluster_sequences,
            )
            residues = connection.execute("SELECT sum(length(sequence)) FROM family").fetchone()[0]
            if residues > config.max_subcluster_residues:
                raise InputValidationError("Subcluster input exceeds the configured residue limit.")
            connection.execute(
                "CREATE TABLE aliases AS SELECT 'p' || lpad(CAST(row_number() "
                "OVER (ORDER BY sequence_id) AS VARCHAR),9,'0') AS tool_id,* FROM family"
            )
            write_fasta(
                connection=connection,
                query="SELECT tool_id,sequence FROM aliases ORDER BY tool_id",
                path=tool_input,
            )
            write_fasta(
                connection=connection,
                query="SELECT sequence_id,sequence FROM aliases ORDER BY tool_id",
                path=output / "sequences.fasta",
            )
            aliases = {
                key: (identifier, length)
                for key, identifier, length in connection.execute(
                    "SELECT tool_id,sequence_id,length(sequence) FROM aliases ORDER BY tool_id"
                ).fetchall()
            }
        write_tsv_atomic(
            path=output / "subcluster_id_map.tsv",
            fieldnames=("tool_id", "sequence_id", "sequence_length"),
            records=(
                {"tool_id": key, "sequence_id": record[0], "sequence_length": record[1]}
                for key, record in aliases.items()
            ),
        )
        database = root / "proteins"
        run_command(
            config=config,
            output=output,
            arguments=[
                config.diamond,
                "makedb",
                "--in",
                str(tool_input),
                "--db",
                str(database),
                "--threads",
                str(config.threads),
            ],
        )
        common = [
            "--db",
            str(database.with_suffix(".dmnd")),
            "--threads",
            str(config.threads),
            "--memory-limit",
            f"{config.memory_gb}G",
            "--tmpdir",
            scratch,
        ]
        native_path, exact = output / "native_membership.tsv", output / "exact_alignments.tsv"
        run_command(
            config=config,
            output=output,
            arguments=[
                config.diamond,
                "deepclust",
                *common,
                "--id",
                str(controls["identity"]),
                "--mutual-cover",
                str(controls["mutual_coverage"]),
                "--evalue",
                str(controls["evalue"]),
                "--out",
                str(native_path),
                "--header",
            ],
        )
        native = native_partition(path=native_path, aliases=aliases)
        pairs = root / "realign_membership.tsv"
        with (
            native_path.open(encoding="utf-8") as original,
            pairs.open(
                mode="x",
                encoding="utf-8",
            ) as destination,
        ):
            original.readline()
            shutil.copyfileobj(original, destination)
        run_command(
            config=config,
            output=output,
            arguments=[
                config.diamond,
                "realign",
                *common,
                "--clusters",
                str(pairs),
                "--out",
                str(exact),
                "--outfmt",
                "6",
                *REALIGN_FIELDS,
                "--comp-based-stats",
                "0",
            ],
        )
        rows = audit_partition(native=native, aliases=aliases, exact=exact, controls=controls)
    write_partition(rows=rows, output=output, controls=controls)
    LOGGER.info("Audited %d DIAMOND DeepClust members", len(rows))


def subcluster(*, config: WorkbenchConfig, source: Path, output: Path, parameters: dict) -> None:
    """Run one audited family partition locally or on an explicitly configured SSH host."""
    controls = subcluster_controls(parameters=parameters)
    transport = parameters.get("transport", "local")
    if transport == "local":
        subcluster_local(config=config, source=source, output=output, controls=controls)
    elif transport == "ssh" and parameters.get("remote_allowed") is True:
        from .subcluster_remote import subcluster_remote

        subcluster_remote(config=config, source=source, output=output, controls=controls)
    else:
        raise InputValidationError(
            "Remote subclustering requires explicit sequence transfer consent."
        )
