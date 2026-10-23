"""Failure injection at tool boundaries and real multi-batch sequence processing."""

import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench import analysis, operations
from protein_signatures.workbench.analysis import (
    align_family,
    executable_version,
    index_fasta,
    run_command,
    search_sequences,
)
from protein_signatures.workbench.catalogue import open_catalogue


def test_fasta_batches_count_limits_and_no_silent_sampling(workbench_config, tmp_path):
    """The 4096-record boundary remains complete without collecting all sequences in Python."""
    path = tmp_path / "many.fa"
    path.write_text(data="\n" + "".join(f">p{i}\nACD\n" for i in range(4096)), encoding="utf-8")
    with open_catalogue(config=workbench_config, with_index=False) as connection:
        assert index_fasta(connection=connection, path=path, name="family", maximum=4096) == 4096
        assert connection.execute("SELECT count(*) FROM family").fetchone()[0] == 4096
        with pytest.raises(InputValidationError):
            index_fasta(connection=connection, path=path, name="other", maximum=4096)
    with open_catalogue(config=workbench_config, with_index=False) as connection:
        with pytest.raises(InputValidationError, match="exceeds"):
            index_fasta(connection=connection, path=path, name="family", maximum=10)


def test_tool_version_provenance_is_bounded_and_unavailability_is_explicit(
    workbench_config, monkeypatch
):
    """Scientific job provenance includes actual installed-tool output without blocking analyses."""
    result = executable_version(
        config=replace(workbench_config, mafft=sys.executable), executable=sys.executable
    )
    assert result["status"] == 0 and "Python" in result["output"]
    assert (
        executable_version(config=workbench_config, executable="other")["status"]
        == "NOT_CONFIGURED_ANALYTICAL_TOOL"
    )
    monkeypatch.setattr(
        analysis.subprocess,
        "run",
        lambda **_kwargs: (_ for _ in ()).throw(subprocess.TimeoutExpired(cmd="tool", timeout=1)),
    )
    assert (
        executable_version(config=workbench_config, executable=workbench_config.diamond)["status"]
        == "UNAVAILABLE"
    )


@pytest.mark.parametrize(
    "raw",
    [
        "malformed\n",
        "a\tb\t90\t3\t0\t3\t1\t3\t1\t3\t1e-20\t100\n",
        "a\tb\tnan\t3\t3\t3\t1\t3\t1\t3\t1e-20\t100\n",
        "a\tb\t101\t3\t3\t3\t1\t3\t1\t3\t1e-20\t100\n",
        "a\tb\t90\t3\t3\t3\t1\t3\t1\t3\t-1\t100\n",
        "a\tb\t90\t3\t3\t3\t1\t4\t1\t3\t1e-20\t100\n",
        "a\tb\t90\t3\t3\t3\t1\t3\t1\t4\t1e-20\t100\n",
        "a\tb\t90\t3\t3.5\t3\t1\t3\t1\t3\t1e-20\t100\n",
    ],
)
def test_malformed_diamond_evidence_is_rejected(workbench_config, tmp_path, monkeypatch, raw):
    """Malformed tool output may not be promoted into biological similarity evidence."""
    monkeypatch.setattr(
        analysis,
        "run_command",
        lambda **_kwargs: (tmp_path / "hits.raw.tsv").write_text(data=raw, encoding="utf-8"),
    )
    with pytest.raises(InputValidationError):
        search_sequences(
            config=workbench_config,
            query=tmp_path / "query",
            database="db",
            parameters={},
            output=tmp_path,
            maximum_hits=1,
        )


@pytest.mark.parametrize(
    "identity,qlen,slen,evalue",
    [(10, 3, 3, "1e-20"), (90, 10, 3, "1e-20"), (90, 3, 10, "1e-20"), (90, 3, 3, "1")],
)
def test_all_member_thresholds_are_rechecked(
    workbench_config, tmp_path, monkeypatch, identity, qlen, slen, evalue
):
    """Each failed threshold drops the hit even if an external tool unexpectedly returns it."""
    raw = f"a\tb\t{identity}\t3\t{qlen}\t{slen}\t1\t3\t1\t3\t{evalue}\t100\n"
    monkeypatch.setattr(
        analysis,
        "run_command",
        lambda **_kwargs: (tmp_path / "hits.raw.tsv").write_text(data=raw, encoding="utf-8"),
    )
    result = search_sequences(
        config=workbench_config,
        query=tmp_path / "query",
        database="db",
        parameters={},
        output=tmp_path,
        maximum_hits=1,
    )
    assert pd.read_csv(result, sep="\t").empty


@pytest.mark.parametrize(
    "matrix",
    [
        "",
        ">seq_000000\nACD\n",
        ">seq_000000\nACD\n>seq_000001\nACDE\n",
        ">seq_000000\nACF\n>seq_000001\nACE\n",
        ">seq_000000\nACD\n>seq_000000\nACD\n>seq_000001\nACE\n",
        ">seq_000000\nACD\n>seq_000001\nACE\n>seq_000000\nACD\n",
    ],
)
def test_invalid_mafft_output_is_never_published(workbench_config, tmp_path, monkeypatch, matrix):
    """An alignment must preserve each exact input ID and every non-gap residue."""
    source = tmp_path / "family.fa"
    source.write_text(data=">a\nACD\n>b\nACE\n", encoding="utf-8")
    monkeypatch.setattr(
        analysis,
        "run_command",
        lambda **kwargs: (tmp_path / kwargs["stdout_name"]).write_text(
            data=matrix, encoding="utf-8"
        ),
    )
    with pytest.raises(InputValidationError):
        align_family(config=workbench_config, source=source, parameters={}, output=tmp_path)


def test_expanded_alignment_matrix_limit_and_all_gap_conservation(
    workbench_config, tmp_path, monkeypatch
):
    """The output matrix, not just unaligned input size, is bounded and gaps stay explicit."""
    source = tmp_path / "family.fa"
    source.write_text(data=">a\nACD\n>b\nACE\n", encoding="utf-8")
    matrix = "\n>seq_000000\nACD-\n>seq_000001\nACE-\n"
    monkeypatch.setattr(
        analysis,
        "run_command",
        lambda **kwargs: (Path(kwargs["output"]) / kwargs["stdout_name"]).write_text(
            data=matrix, encoding="utf-8"
        ),
    )
    output = tmp_path / "bounded"
    output.mkdir()
    with pytest.raises(InputValidationError, match="matrix exceeds"):
        align_family(
            config=replace(workbench_config, max_alignment_cells=6),
            source=source,
            parameters={},
            output=output,
        )
    output = tmp_path / "complete"
    output.mkdir()
    align_family(
        config=workbench_config,
        source=source,
        parameters={"selected_ids": ["a", "b"]},
        output=output,
    )
    table = pd.read_csv(output / "conservation.tsv", sep="\t")
    assert table.loc[3, "gap_fraction"] == 1 and table.loc[3, "consensus"] == "-"


def test_run_command_success_and_interruption_cleanup(workbench_config, tmp_path, monkeypatch):
    """A successful command is logged and interrupted tools are terminated and unmarked."""
    run_command(config=workbench_config, arguments=[sys.executable, "--version"], output=tmp_path)
    assert not (tmp_path / "active_process.json").exists()
    process = SimpleNamespace(pid=2_000_000_000)
    process.wait = lambda **_kwargs: (_ for _ in ()).throw(KeyboardInterrupt())
    monkeypatch.setattr(analysis.subprocess, "Popen", lambda **_kwargs: process)
    stopped = []
    monkeypatch.setattr(
        analysis, "terminate_process", lambda **kwargs: stopped.append(kwargs["process"])
    )
    with pytest.raises(KeyboardInterrupt):
        run_command(config=workbench_config, arguments=["mock"], output=tmp_path)
    assert stopped == [process] and not (tmp_path / "active_process.json").exists()


@pytest.mark.parametrize(
    "parameters", [{"deduplicate": "yes"}, {"taxon_ids": "3702"}, {"taxon_ids": [3702]}]
)
def test_screen_controls_require_explicit_reviewed_metadata(workbench_config, tmp_path, parameters):
    """Type coercion cannot silently enable filters or turn unknown taxonomy into evidence."""
    with pytest.raises(InputValidationError):
        operations.screen_family(
            config=replace(workbench_config, taxonomy=""),
            source=tmp_path / "source",
            seeds=tmp_path / "seeds",
            parameters=parameters,
            output=tmp_path,
        )


def test_screen_without_deduplication_and_no_hit_audit(
    workbench_config, tmp_path, fake_tools, monkeypatch
):
    """All input proteins retain audit rows, including sequences with no retained hit."""
    source = tmp_path / "family.fa"
    seeds = tmp_path / "seeds.fa"
    source.write_text(data=">a\nACDEFGHIK\n>b\nACDEFGHIK\n", encoding="utf-8")
    seeds.write_text(data=">local_seed\nACDEFGHIK\n", encoding="utf-8")
    output = tmp_path / "all"
    output.mkdir()
    operations.screen_family(
        config=workbench_config,
        source=source,
        seeds=seeds,
        parameters={"deduplicate": False},
        output=output,
    )
    assert (output / "candidates.fasta").read_text().count(">") == 2

    def no_hits(**kwargs):
        """Supply a valid header-only search result from a successful no-hit search."""
        path = kwargs["output"] / "empty_hits.tsv"
        path.write_text(
            data="\t".join(
                [*analysis.HIT_FIELDS, "query_coverage_percent", "subject_coverage_percent"]
            )
            + "\n",
            encoding="utf-8",
        )
        return path

    monkeypatch.setattr(operations, "search_sequences", no_hits)
    output = tmp_path / "none"
    output.mkdir()
    operations.screen_family(
        config=workbench_config, source=source, seeds=seeds, parameters={}, output=output
    )
    audit = pd.read_csv(output / "candidate_audit.tsv", sep="\t")
    assert set(audit.selection_status) == {"NO_RETAINED_HIT"}
