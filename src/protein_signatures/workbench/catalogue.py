"""Bounded DuckDB queries and row-group extraction of mounted DeepClust data."""

from __future__ import annotations

import csv
import hashlib
import logging
import re
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from protein_signatures.errors import InputValidationError

from .config import WorkbenchConfig

LOGGER = logging.getLogger(__name__)
RESIDUES = frozenset("ABCDEFGHIKLMNOPQRSTUVWXYZ")


def identifiers(*, values: Sequence[str], maximum: int = 10_000) -> tuple[str, ...]:
    """Validate exact IDs without case folding or removing accession versions.

    Args:
        values: Input identifiers.
        maximum: Maximum input count.

    Returns:
        Unique identifiers in input order.

    Raises:
        InputValidationError: If the collection or an ID is malformed.
    """
    if isinstance(values, (str, bytes)) or not 1 <= len(values) <= maximum:
        raise InputValidationError(f"Supply between 1 and {maximum} exact identifiers.")
    if any(
        not isinstance(value, str)
        or not value
        or len(value) > 2000
        or any(c.isspace() or ord(c) < 32 or ord(c) == 127 for c in value)
        for value in values
    ):
        raise InputValidationError("Identifiers must be non-empty and contain no whitespace.")
    return tuple(dict.fromkeys(values))


def column(*, columns: Sequence[str], choices: Sequence[str], optional: bool = False) -> str:
    """Resolve explicit source aliases and safely quote the actual column.

    Args:
        columns: Physical names.
        choices: Supported aliases in preference order.
        optional: Return SQL NULL when absent.

    Returns:
        Quoted column or NULL.

    Raises:
        InputValidationError: If no required alias exists.
    """
    name = next((item for item in choices if item in columns), None)
    if name is None:
        if optional:
            return "NULL"
        raise InputValidationError(f"DeepClust source needs one of: {', '.join(choices)}")
    return '"' + name.replace('"', '""') + '"'


def source_view(*, connection: object, path: str, name: str) -> list[str]:
    """Open a local Parquet or TSV lazily without collecting its records.

    Args:
        connection: Private DuckDB connection.
        path: Trusted operator source file.
        name: Internal fixed view name.

    Returns:
        Physical column names.

    Raises:
        InputValidationError: If the source is unavailable.
    """
    if not path or not Path(path).is_file():
        raise InputValidationError(f"Configure a readable local source for {name}: {path}")
    relation = (
        connection.read_csv(path_or_buffer=path, delimiter="\t", header=True, all_varchar=True)
        if Path(path).suffix == ".tsv"
        else connection.read_parquet(path_or_buffer=path)
    )
    relation.create_view(view_name=name)
    return relation.columns


@contextmanager
def open_catalogue(*, config: WorkbenchConfig, with_index: bool = True) -> Iterator[object]:
    """Open bounded local analytical views with spill files on the workspace disk.

    Args:
        config: Trusted server configuration.
        with_index: Expose the cluster index.

    Yields:
        Private DuckDB connection; source files remain unchanged.
    """
    spill = Path(config.workspace) / "temporary"
    spill.mkdir(parents=True, exist_ok=True, mode=0o700)
    connection = duckdb.connect(
        config={
            "threads": config.threads,
            "memory_limit": f"{config.memory_gb}GB",
            "temp_directory": str(spill),
            "preserve_insertion_order": False,
            "autoinstall_known_extensions": False,
            "autoload_known_extensions": False,
        }
    )
    try:
        if with_index:
            columns = source_view(connection=connection, path=config.index, name="index_source")
            key = column(columns=columns, choices=("CLUSTER", "cluster_id", "ClusterID", "f2"))
            size = column(columns=columns, choices=("NMEMBER", "member_count", "size"))
            groups = column(columns=columns, choices=("RowGroup", "row_groups"))
            connection.execute(
                f"CREATE VIEW clusters AS SELECT CAST({key} AS VARCHAR) AS cluster_id, "
                f"TRY_CAST({size} AS BIGINT) AS member_count, "
                f"CAST({groups} AS VARCHAR) AS row_groups FROM index_source"
            )
        yield connection
    except duckdb.Error as error:
        LOGGER.exception("DeepClust query failed")
        raise InputValidationError(f"DeepClust query failed: {error}") from error
    finally:
        connection.close()


def resolve_ids(*, config: WorkbenchConfig, values: Sequence[str]) -> pd.DataFrame:
    """Resolve exact membership, preserving unmatched, ambiguous and unindexed records.

    Args:
        config: Mounted index, mapping and optional identifier bridge.
        values: Exact protein IDs; HOG IDs are not global protein IDs.

    Returns:
        Many-to-many mapping and explicit per-input status records.
    """
    requested = identifiers(values=values, maximum=config.max_seed_sequences)
    aliases = []
    for value in requested:
        parts = value.split("|")
        base = [value, parts[1]] if len(parts) >= 3 and parts[0] in {"sp", "tr"} else [value]
        for alias in dict.fromkeys(
            [*base, *(p + v for p in config.identifier_prefixes for v in base)]
        ):
            aliases.append({"input_id": value, "alias": alias, "mapping_evidence": "EXACT_ALIAS"})
    with open_catalogue(config=config) as connection:
        columns = source_view(connection=connection, path=config.sequence_map, name="map_source")
        seq = column(columns=columns, choices=("SEQID", "sequence_id", "SeqID", "SeqId", "f0"))
        key = column(
            columns=columns, choices=("CLUSTERID", "cluster_id", "CLUSTER", "ClustID", "f2", "f1")
        )
        connection.register(view_name="aliases", python_object=pd.DataFrame(data=aliases))
        bridge = ""
        if config.identifier_bridge:
            source_view(connection=connection, path=config.identifier_bridge, name="bridge")
            connection.execute("SELECT input_id,sequence_id FROM bridge LIMIT 0")
            bridge = (
                " UNION SELECT a.input_id,CAST(b.sequence_id AS VARCHAR),'DECLARED_BRIDGE' "
                "FROM aliases a JOIN bridge b ON a.alias=b.input_id"
            )
        result = connection.execute(
            "WITH expanded AS (SELECT * FROM aliases" + bridge + "), matched AS ("
            f"SELECT DISTINCT a.input_id,a.alias,a.mapping_evidence,CAST(m.{seq} AS VARCHAR) "
            f"AS sequence_id,CAST(m.{key} AS VARCHAR) AS cluster_id "
            f"FROM expanded a JOIN map_source m ON a.alias=CAST(m.{seq} AS VARCHAR)), "
            "inputs AS (SELECT DISTINCT input_id FROM aliases) "
            "SELECT i.input_id,m.alias,m.sequence_id,m.cluster_id,c.member_count, "
            "m.mapping_evidence,CASE WHEN m.cluster_id IS NULL THEN 'UNMATCHED' "
            "WHEN c.cluster_id IS NULL THEN 'MAPPED_UNINDEXED' ELSE 'MAPPED' END AS status "
            "FROM inputs i LEFT JOIN matched m USING(input_id) "
            "LEFT JOIN clusters c USING(cluster_id) ORDER BY input_id,cluster_id,sequence_id "
            "LIMIT 1000001"
        ).fetchdf()
    if len(result) > 1_000_000:
        raise InputValidationError(
            "Mapping is too ambiguous; split this input into smaller batches."
        )
    return result


def row_groups(*, value: object, maximum: int) -> tuple[int, ...]:
    """Parse official scalar or multi-group index references.

    Args:
        value: Published row-group reference.
        maximum: Physical row-group count.

    Returns:
        Sorted unique offsets.

    Raises:
        InputValidationError: If a reference is invalid.
    """
    tokens = re.split(pattern=r"[,;\s\[\]]+", string=str(value))
    tokens = [item for item in tokens if item]
    if not tokens or any(not item.isdigit() for item in tokens):
        raise InputValidationError("Missing or malformed DeepClust row-group reference.")
    groups = tuple(sorted({int(item) for item in tokens}))
    if any(item >= maximum for item in groups):
        raise InputValidationError("DeepClust row-group reference exceeds the sequence file.")
    return groups


def normalise_sequence(*, value: object) -> str:
    """Normalise protein whitespace and terminal stops, rejecting corrupt sequences.

    Args:
        value: Source sequence.

    Returns:
        Ungapped upper-case amino-acid sequence.

    Raises:
        InputValidationError: If the sequence is empty, gapped or invalid.
    """
    if not isinstance(value, str):
        raise InputValidationError("Protein sequence is missing or is not text.")
    sequence = "".join(value.split()).upper().removesuffix("*")
    if not sequence or set(sequence) - RESIDUES:
        raise InputValidationError("Protein contains gaps, internal stops or invalid residues.")
    return sequence


def retrieve_clusters(*, config: WorkbenchConfig, values: Sequence[str], output: Path) -> None:
    """Stream complete indexed families; verify every cluster count before publication.

    Args:
        config: Mounted DeepClust sources and extraction bounds.
        values: Exact global centroids, at most 1000.
        output: Empty job output directory.

    Raises:
        InputValidationError: If IDs, counts, sequences or limits are invalid.
    """
    selected = identifiers(values=values, maximum=1000)
    with open_catalogue(config=config) as connection:
        index = connection.execute(
            "SELECT * FROM clusters WHERE cluster_id IN (SELECT unnest(?))", [list(selected)]
        ).fetchdf()
        if (
            set(index.cluster_id) != set(selected)
            or index.cluster_id.duplicated().any()
            or index.member_count.isna().any()
            or (index.member_count < 1).any()
        ):
            raise InputValidationError(
                "Cluster missing, duplicated or invalid in the row-group index."
            )
        expected = dict(zip(index.cluster_id, index.member_count, strict=True))
        if sum(expected.values()) > config.max_sequences:
            raise InputValidationError(
                "Complete clusters exceed max_sequences; choose fewer clusters."
            )
        if not config.sequences or not Path(config.sequences).is_file():
            raise InputValidationError("Configure the mounted DeepClust sequence Parquet.")
        with pq.ParquetFile(source=config.sequences) as parquet:
            columns = parquet.schema_arrow.names
            key = column(columns=columns, choices=("f0", "sequence_id")).strip('"')
            seq = column(columns=columns, choices=("f1", "protein_sequence", "sequence")).strip('"')
            cluster = column(columns=columns, choices=("f2", "cluster_id")).strip('"')
            groups = sorted(
                {
                    g
                    for value in index.row_groups
                    for g in row_groups(value=value, maximum=parquet.num_row_groups)
                }
            )
            counts = dict.fromkeys(selected, 0)
            with (
                (output / "family.fasta").open(mode="x", encoding="utf-8") as fasta,
                (output / "members.tsv").open(mode="x", encoding="utf-8", newline="") as table,
            ):
                writer = csv.writer(table, delimiter="\t", lineterminator="\n")
                writer.writerow(["sequence_id", "cluster_id", "sequence_length", "sequence_sha256"])
                for batch in parquet.iter_batches(
                    batch_size=4096, row_groups=groups, columns=[key, seq, cluster]
                ):
                    keep = batch.filter(
                        pc.is_in(batch.column(cluster), value_set=pa.array(selected))
                    )
                    for record in keep.to_pylist():
                        identifier = identifiers(values=[record[key]])[0]
                        protein = normalise_sequence(value=record[seq])
                        parent = record[cluster]
                        counts[parent] += 1
                        if counts[parent] > expected[parent]:
                            raise InputValidationError(
                                "Retrieved more records than the published count."
                            )
                        fasta.write(f">{identifier}\n{protein}\n")
                        writer.writerow(
                            [
                                identifier,
                                parent,
                                len(protein),
                                hashlib.sha256(protein.encode()).hexdigest(),
                            ]
                        )
            if counts != expected:
                raise InputValidationError(
                    "Retrieved counts disagree with the published cluster index."
                )
        source_view(connection=connection, path=str(output / "members.tsv"), name="members")
        duplicate = connection.execute(
            "SELECT sequence_id FROM members GROUP BY sequence_id HAVING count(*)>1 LIMIT 1"
        ).fetchone()
        if duplicate:
            raise InputValidationError(f"Duplicate global sequence identifier: {duplicate[0]}")
        summary = connection.execute(
            "SELECT cluster_id,count(*) AS retrieved_members,"
            "min(CAST(sequence_length AS BIGINT)) AS minimum_length,"
            "max(CAST(sequence_length AS BIGINT)) AS maximum_length FROM members "
            "GROUP BY cluster_id ORDER BY cluster_id"
        ).fetchdf()
        summary["published_members"] = summary.cluster_id.map(arg=expected)
        summary["physical_row_groups_read"] = len(groups)
        summary.to_csv(path_or_buf=output / "clusters.tsv", sep="\t", index=False)
        if config.taxonomy:
            source_view(connection=connection, path=config.taxonomy, name="taxonomy_source")
            # Conflicting/unreviewed metadata stays unclassified; names are never inferred.
            columns = connection.execute("SELECT * FROM taxonomy_source LIMIT 0").description
            reviewed = (
                "WHERE mapping_status='REVIEWED'"
                if "mapping_status" in {c[0] for c in columns}
                else ""
            )
            connection.execute(
                "CREATE VIEW taxonomy AS SELECT sequence_id,min(species) AS species,"
                "min(taxon_id) AS taxon_id FROM taxonomy_source "
                + reviewed
                + " GROUP BY sequence_id HAVING count(DISTINCT taxon_id)=1"
            )
            taxa = connection.sql(
                query=(
                    "SELECT m.cluster_id,coalesce(t.species,'Unclassified') AS species,"
                    "t.taxon_id,count(*) AS sequence_count FROM members m LEFT JOIN taxonomy t "
                    "USING(sequence_id) GROUP BY m.cluster_id,t.species,t.taxon_id "
                    "ORDER BY m.cluster_id,sequence_count DESC"
                )
            )
            taxa.write_csv(file_name=str(output / "taxonomy_summary.tsv"), sep="\t", header=True)
    LOGGER.info("Retrieved %d complete clusters, %d sequences", len(selected), sum(counts.values()))
