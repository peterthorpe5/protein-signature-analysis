"""Small official-schema DeepClust fixtures with physical row-group boundaries."""

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from protein_signatures.workbench.config import WorkbenchConfig


@pytest.fixture(name="workbench_config")
def fixture_workbench_config(tmp_path: Path) -> WorkbenchConfig:
    """Provide real Parquet files containing two indexed families and one omitted record."""
    source = tmp_path / "release"
    source.mkdir()
    pq.write_table(
        table=pa.table(
            data={
                "f0": ["UNIPROT_P1", "B", "C", "NR_D.1", "E", "F"],
                "f1": ["ACDEFGHIK*", "ACDEFGHIK", "ACDEFGHIL", "MNPQRSTVW", "MNPQRSTVY", "AAAAAAA"],
                "f2": ["UNIPROT_P1"] * 3 + ["NR_D.1"] * 2 + ["F"],
                "f3": list(range(6)),
            }
        ),
        where=source / "sequences.parquet",
        row_group_size=2,
    )
    pq.write_table(
        table=pa.table(
            data={
                "CLUSTER": ["UNIPROT_P1", "NR_D.1"],
                "NMEMBER": [3, 2],
                "RowStart": [0, 3],
                "RowGroup": ["0,1", "1,2"],
            }
        ),
        where=source / "index.parquet",
    )
    pq.write_table(
        table=pa.table(
            data={
                "SEQID": ["UNIPROT_P1", "B", "C", "NR_D.1", "E", "F", "multi", "multi"],
                "CLUSTERID": ["UNIPROT_P1"] * 3 + ["NR_D.1"] * 2 + ["F", "UNIPROT_P1", "NR_D.1"],
            }
        ),
        where=source / "map.parquet",
    )
    (source / "taxonomy.tsv").write_text(
        data=(
            "sequence_id\tspecies\ttaxon_id\tmapping_status\n"
            "UNIPROT_P1\tPlant\t3702\tREVIEWED\nB\tPlant\t3702\tREVIEWED\n"
            "C\tHuman\t9606\tREVIEWED\nNR_D.1\tUnknown\t1\tUNREVIEWED\n"
        ),
        encoding="utf-8",
    )
    (source / "bridge.tsv").write_text(data="input_id\tsequence_id\nlocal_b\tB\n", encoding="utf-8")
    return WorkbenchConfig(
        workspace=str(tmp_path / "workspace"),
        index=str(source / "index.parquet"),
        sequences=str(source / "sequences.parquet"),
        sequence_map=str(source / "map.parquet"),
        taxonomy=str(source / "taxonomy.tsv"),
        identifier_bridge=str(source / "bridge.tsv"),
        threads=1,
        memory_gb=1,
        max_alignment_cells=1000,
    )


@pytest.fixture(name="seed_bytes")
def fixture_seed_bytes() -> bytes:
    """Return one query using a local identifier absent from the global map."""
    return b">local_seed\nACDEFGHIK\n"


@pytest.fixture(name="config_path")
def fixture_config_path(workbench_config: WorkbenchConfig, tmp_path: Path) -> Path:
    """Write the operator fixture in the same format used on Sabrina."""
    import json
    from dataclasses import asdict

    path = tmp_path / "workbench.toml"
    path.write_text(
        data="\n".join(
            f"{key} = {json.dumps(obj=value)}" for key, value in asdict(workbench_config).items()
        ),
        encoding="utf-8",
    )
    return path


@pytest.fixture(name="fake_tools")
def fixture_fake_tools(monkeypatch: pytest.MonkeyPatch):
    """Replace only external binaries, leaving real Parquet, queue and SQL paths active."""
    from protein_signatures.workbench import analysis, operations

    calls = []

    def fake_command(*, config, arguments, output, stdout_name=None):
        """Supply deterministic outputs following each external tool's actual contract."""
        calls.append(list(arguments))
        if "makedb" in arguments:
            database = Path(arguments[arguments.index("--db") + 1])
            database.with_suffix(".dmnd").write_bytes(data=b"fixture database")
        elif "blastp" in arguments:
            query = Path(arguments[arguments.index("--query") + 1])
            destination = Path(arguments[arguments.index("--out") + 1])
            database = arguments[arguments.index("--db") + 1]
            lines = []
            for identifier, sequence in analysis.iter_fasta(path=query):
                subject = "local_seed" if "seed_database" in database else "UNIPROT_P1"
                lines.append(
                    f"{identifier}\t{subject}\t90\t9\t{len(sequence)}\t9\t1\t9\t1\t9\t1e-20\t100\n"
                )
            destination.write_text(data="".join(lines), encoding="utf-8")
        elif arguments[0] == config.mafft:
            (output / stdout_name).write_bytes(data=Path(arguments[-1]).read_bytes())
        elif arguments[0] == config.fasttree:
            (output / stdout_name).write_text(
                data="(seq_000000:0.1,seq_000001:0.2,seq_000002:0.3);\n", encoding="utf-8"
            )
        elif arguments[0] == config.hmmscan:
            destination = Path(arguments[arguments.index("--domtblout") + 1])
            destination.write_text(
                data=(
                    "# mock HMMER domtblout\nFbox PF00646.1 40 B - 9 1e-8 40 0 1 1 1e-8 1e-8"
                    " 40 0 1 9 1 9 1 9 0.99 domain\n"
                ),
                encoding="utf-8",
            )
        else:
            raise AssertionError(arguments)

    monkeypatch.setattr(analysis, "run_command", fake_command)
    monkeypatch.setattr(operations, "run_command", fake_command)
    return calls
