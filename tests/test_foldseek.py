"""Tests for safe Foldseek commands, parsing, caching and failure handling."""

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import protein_signatures.foldseek as fs
from protein_signatures.checksums import sha256_file
from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.models import (
    FoldEvidenceStatus,
    FoldseekSettings,
    StructureAnalysisEligibility,
    StructureCoverageScope,
    StructureRecord,
)


def _line(
    query: str,
    target: str,
    *,
    aligned: int | str = 8,
    qlen: int = 10,
    tlen: int = 10,
    qtm: float | str = 0.8,
    ttm: float = 0.7,
) -> str:
    """Return one row in the exact declared Foldseek output format."""

    values = (
        query,
        target,
        aligned,
        qlen,
        tlen,
        1e-10,
        100,
        0.75,
        qtm,
        ttm,
    )
    return "\t".join(str(value) for value in values) + "\n"


def test_build_command_and_version(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Foldseek invocation should be shell-free, explicit and versioned."""

    settings = _settings(cache=tmp_path / "cache")
    command = fs.build_foldseek_command(
        executable="/tools/foldseek",
        query_dir=tmp_path / "queries",
        output_path=tmp_path / "out.tsv",
        temporary_dir=tmp_path / "work",
        settings=settings,
        threads=8,
    )
    assert command[:2] == ("/tools/foldseek", "easy-search")
    assert "--format-output" in command
    assert command[-1] == ",".join(fs.FOLDSEEK_COLUMNS)
    with pytest.raises(InputValidationError, match="threads"):
        fs.build_foldseek_command(
            executable="foldseek",
            query_dir=tmp_path,
            output_path=tmp_path / "x",
            temporary_dir=tmp_path,
            settings=settings,
            threads=0,
        )
    monkeypatch.setattr(
        fs.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(
            returncode=0, stdout="Foldseek 10.0\nextra\n", stderr=""
        ),
    )
    assert fs.foldseek_version(executable="foldseek") == "Foldseek 10.0"
    monkeypatch.setattr(
        fs.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=1, stdout="", stderr="bad"),
    )
    with pytest.raises(ExternalToolError, match="determine"):
        fs.foldseek_version(executable="foldseek")


def test_parser_symmetrises_and_keeps_best_hit(tmp_path: Path) -> None:
    """Directional hits should become one conservative best unordered comparison."""

    path = tmp_path / "foldseek.tsv"
    path.write_text(
        _line("q1", "q1", aligned=100, qlen=100, tlen=100, qtm=1.0, ttm=1.0)
        + _line("q1.pdb_A", "q2_1", aligned=80, qlen=100, tlen=160, qtm=0.8, ttm=0.6)
        + _line("q2", "q1", aligned=90, qlen=160, tlen=100, qtm=0.9, ttm=0.85),
        encoding="utf-8",
    )
    rows = fs.parse_foldseek_all_vs_all(
        path=path,
        query_to_protein={"q1": "p1", "q2": "p2"},
        tool_version="v1",
        comparison_universe_id="campaign_1",
    )
    assert len(rows) == 1
    row = rows[0]
    assert row.protein_a_id == "p1"
    assert row.protein_b_id == "p2"
    assert row.tm_score == pytest.approx(0.85)
    assert row.aligned_residue_count == 90
    assert row.coverage_a == pytest.approx(0.9)
    assert row.coverage_b == pytest.approx(90 / 160)
    assert row.comparison_universe_id == "campaign_1"
    assert row.coverage_scope == StructureCoverageScope.STRUCTURE_MODEL_RESIDUES


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("unknown\tq2\t1\t1\t1\t0\t1\t1\t1\t1\n", "query identifier"),
        (_line("q1", "q2", aligned="bad"), "must be an integer"),
        (_line("q1", "q2", aligned=0), "must be positive"),
        (_line("q1", "q2", qtm="bad"), "must be numeric"),
        (_line("q1", "q2", qtm=2), "zero to one"),
        ("q1\tq2\t1\n", "field count"),
    ],
)
def test_parser_rejects_malformed_rows(tmp_path: Path, line: str, message: str) -> None:
    """Malformed identifiers, numbers and row widths should fail with context."""

    path = tmp_path / "bad.tsv"
    path.write_text(line, encoding="utf-8")
    with pytest.raises(InputValidationError, match=message):
        fs.parse_foldseek_all_vs_all(
            path=path,
            query_to_protein={"q1": "p1", "q2": "p2"},
            tool_version="v1",
            comparison_universe_id="campaign_1",
        )


def test_run_foldseek_creates_and_reuses_checksum_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A successful all-versus-all unit should be atomically cached and reused."""

    structures = _structures(root=tmp_path)
    settings = _settings(cache=tmp_path / "cache")
    monkeypatch.setattr(fs, "_resolve_executable", lambda **_kwargs: "/bin/foldseek")
    monkeypatch.setattr(fs, "foldseek_version", lambda **_kwargs: "foldseek-v1")
    calls = {"easy": 0}

    def fake_run(command: tuple[str, ...], **_kwargs: object) -> SimpleNamespace:
        calls["easy"] += 1
        query_ids = sorted(path.stem for path in Path(command[2]).iterdir())
        Path(command[4]).write_text(
            _line(query_ids[0], query_ids[1], aligned=4, qlen=4, tlen=4),
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(fs.subprocess, "run", fake_run)
    first = fs.run_foldseek_all_vs_all(structures=structures, settings=settings, threads=2)
    assert first.reused is False
    assert len(first.comparisons) == 1
    assert first.raw_output_path.is_file()
    marker = json.loads(first.completion_manifest_path.read_text(encoding="utf-8"))
    assert marker["raw_output_sha256"] == sha256_file(path=first.raw_output_path)
    assert marker["comparison_universe_id"].startswith("FOLDSEEK_")
    second = fs.run_foldseek_all_vs_all(structures=structures, settings=settings, threads=2)
    assert second.reused is True
    assert second.cache_key == first.cache_key
    assert calls["easy"] == 1


def test_run_foldseek_reports_tool_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Insufficient models, failed commands and empty output must remain distinct failures."""

    structures = _structures(root=tmp_path)
    settings = _settings(cache=tmp_path / "cache")
    with pytest.raises(InputValidationError, match="at least two"):
        fs.run_foldseek_all_vs_all(structures=structures[:1], settings=settings, threads=1)
    monkeypatch.setattr(fs, "_resolve_executable", lambda **_kwargs: "/bin/foldseek")
    monkeypatch.setattr(fs, "foldseek_version", lambda **_kwargs: "foldseek-v1")
    monkeypatch.setattr(
        fs.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=3, stdout="", stderr="tool failed"),
    )
    with pytest.raises(ExternalToolError, match="tool failed"):
        fs.run_foldseek_all_vs_all(structures=structures, settings=settings, threads=1)
    monkeypatch.setattr(
        fs.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )
    with pytest.raises(ExternalToolError, match="non-empty"):
        fs.run_foldseek_all_vs_all(
            structures=structures,
            settings=_settings(cache=tmp_path / "cache2"),
            threads=1,
        )


def test_foldseek_excludes_ineligible_models_and_prevents_truncation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Only eligible models may run and max-seqs must cover the full model universe."""

    structures = _structures(root=tmp_path, count=3)
    ineligible = replace(
        structures[2],
        analysis_eligibility_status=StructureAnalysisEligibility.INELIGIBLE_LOW_CONFIDENCE,
    )
    monkeypatch.setattr(fs, "_resolve_executable", lambda **_kwargs: "/bin/foldseek")
    monkeypatch.setattr(fs, "foldseek_version", lambda **_kwargs: "foldseek-v1")
    observed = {"query_count": 0}

    def fake_run(command: tuple[str, ...], **_kwargs: object) -> SimpleNamespace:
        query_ids = sorted(path.stem for path in Path(command[2]).iterdir())
        observed["query_count"] = len(query_ids)
        Path(command[4]).write_text(
            _line(query_ids[0], query_ids[1], aligned=4, qlen=4, tlen=4),
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(fs.subprocess, "run", fake_run)
    evidence = fs.run_foldseek_all_vs_all(
        structures=(*structures[:2], ineligible),
        settings=_settings(cache=tmp_path / "eligible-cache", maximum_hits=2),
        threads=1,
    )
    assert len(evidence.comparisons) == 1
    assert observed["query_count"] == 2
    with pytest.raises(InputValidationError, match="must be at least"):
        fs.run_foldseek_all_vs_all(
            structures=structures,
            settings=_settings(cache=tmp_path / "truncated-cache", maximum_hits=2),
            threads=1,
        )


def test_cache_and_query_helpers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Cache identities, query links, executable lookup and locks should be deterministic."""

    structures = _structures(root=tmp_path)
    settings = _settings(cache=tmp_path / "cache")
    key = fs._cache_key(structures=structures, settings=settings, tool_version="v1")
    assert len(key) == 64
    assert key == fs._cache_key(
        structures=tuple(reversed(structures)), settings=settings, tool_version="v1"
    )
    query_id = fs._query_id(structure_id="unsafe/structure")
    assert query_id.startswith("PSQ_") and "/" not in query_id
    query_dir = tmp_path / "queries"
    query_dir.mkdir()
    fs._link_queries(structures=structures, query_dir=query_dir)
    assert len(tuple(query_dir.iterdir())) == 2
    mapping = {query_id: "protein"}
    assert fs._resolve_query(value=f"{query_id}.pdb_A", query_to_protein=mapping) == "protein"
    with pytest.raises(InputValidationError):
        fs._resolve_query(value="unknown", query_to_protein=mapping)
    monkeypatch.setattr(fs.shutil, "which", lambda value: f"/tools/{value}")
    assert fs._resolve_executable(value="foldseek") == "/tools/foldseek"
    monkeypatch.setattr(fs.shutil, "which", lambda _value: None)
    with pytest.raises(ExternalToolError, match="not found"):
        fs._resolve_executable(value="foldseek")
    lock = tmp_path / "unit.lock"
    descriptor = fs._acquire_lock(path=lock)
    try:
        with pytest.raises(ExternalToolError, match="already running"):
            fs._acquire_lock(path=lock)
    finally:
        os.close(descriptor)
        lock.unlink()


def _settings(*, cache: Path, maximum_hits: int = 100) -> FoldseekSettings:
    """Return bounded test Foldseek settings."""

    return FoldseekSettings(
        enabled=True,
        executable="foldseek",
        cache_dir=cache,
        e_value_threshold=0.001,
        sensitivity=9.5,
        maximum_hits=maximum_hits,
    )


def _structures(*, root: Path, count: int = 2) -> tuple[StructureRecord, ...]:
    """Create checksum-bound eligible coordinate records."""

    rows: list[StructureRecord] = []
    for index in range(1, count + 1):
        suffix = ".cif" if index % 2 == 0 else ".pdb"
        path = root / f"model{index}{suffix}"
        path.write_text(f"MODEL {index}\n", encoding="utf-8")
        rows.append(
            StructureRecord(
                protein_id=f"p{index}",
                structure_id=f"structure{index}",
                structure_source="test",
                structure_version="1",
                coordinate_path=path,
                coordinate_sha256=sha256_file(path=path),
                availability_status="AVAILABLE",
                mean_confidence=90.0,
                fold_id="",
                fold_name="",
                fold_authority="",
                fold_authority_version="",
                fold_evidence_reference="",
                fold_evidence_status=FoldEvidenceStatus.NOT_ASSESSED,
                analysis_eligibility_status=StructureAnalysisEligibility.ELIGIBLE,
                comparison_universe_ids=(),
            )
        )
    return tuple(rows)
