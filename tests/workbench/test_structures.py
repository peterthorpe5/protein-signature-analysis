"""Check local LoL-align score semantics and preserve original structural inputs."""

import hashlib
import json
import os
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from Bio.PDB import PDBIO, PDBParser

from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.workbench import structures
from protein_signatures.workbench.cli import main
from protein_signatures.workbench.jobs import artifact_path, read_job, run_next, submit_job


@pytest.fixture(name="monomer")
def fixture_monomer():
    """Provide a tiny PDB fixture; no external sequence or model service is needed."""
    return (
        b"".join(
            (
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"{i * 2.0:8.3f}{1.0:8.3f}{2.0:8.3f}  1.00 85.00           C\n"
            ).encode()
            for i in range(1, 5)
        )
        + b"END\n"
    )


VALID_LOL = (
    "query\ttarget\t4\t4\t1\t4\t1\t4\tAAAA\tAAAA\t0.9\t0.8\t0.4\t-0.5\t0.6\t"
    "1,0,0,0,1,0,0,0,1\t0,0,0\n"
)


def test_parse_local_model_and_score_semantics(monomer):
    """Keep the LoL score separate from probability, coverage and geometric scores."""
    suffix, sequence, residues = structures.model_residues(content=monomer, maximum=4)
    row, transform = structures.parse_lolalign(text=VALID_LOL, query=sequence, target=sequence)
    assert suffix == ".pdb" and sequence == "AAAA"
    assert residues[2]["position"] == 3 and residues[2]["bfactor"] == 85
    assert row["lol_normalised_score"] == -0.5
    assert row["lol_homology_probability"] == 0.6
    assert row["query_tm_score"] == 0.9 and row["target_tm_score"] == 0.8
    assert row["aligned_residue_pairs"] == 4
    assert transform["direction"] == "target_to_query"
    assert structures.parse_lolalign(text="", query=sequence, target=sequence)[0]["status"] == (
        "NO_ALIGNMENT_FOUND"
    )


@pytest.mark.parametrize(
    "text",
    [
        VALID_LOL + VALID_LOL,
        VALID_LOL.replace("query", "wrong"),
        "short",
        VALID_LOL.replace("\t0.9\t", "\tnan\t"),
        VALID_LOL.replace("\t0.6\t", "\t1.1\t"),
        VALID_LOL.replace("\tAAAA\tAAAA\t", "\tA-AA\tAAAA\t"),
        VALID_LOL.replace("1,0,0,0,1,0,0,0,1", "1,0,0,0,1,0,0,0,-1"),
        VALID_LOL.replace("\t0,0,0", "\t0,0"),
        VALID_LOL.replace("\t0.4\t", "\t-1\t"),
    ],
)
def test_inconsistent_tool_rows_fail(text):
    """Reject scientifically inconsistent, malformed and ambiguous results."""
    with pytest.raises(ExternalToolError):
        structures.parse_lolalign(text=text, query="AAAA", target="AAAA")


def test_model_rejects_other_chains_and_limits(monomer, tmp_path, monkeypatch):
    """Do not silently choose a chain or truncate a large coordinate file."""
    with pytest.raises(InputValidationError, match="one selected"):
        structures.model_residues(
            content=monomer.replace(b"END\n", b"") + monomer.replace(b"ALA A", b"ALA B"),
            maximum=10,
        )
    with pytest.raises(InputValidationError, match="limit"):
        structures.model_residues(content=monomer, maximum=3)
    with pytest.raises(InputValidationError):
        structures.model_residues(content=b"not a structure", maximum=4)
    model = tmp_path / "model"
    model.write_bytes(monomer)
    assert structures.read_model(path=model) == monomer
    monkeypatch.setattr(structures, "MODEL_BYTES", 12)
    with pytest.raises(InputValidationError, match="Coordinate files"):
        structures.read_model(path=model)


def test_capability_uses_installed_binary(monkeypatch):
    """Recognise the advertised LoL flag, not a guessed Foldseek version floor."""

    def probe(*, args, **_kwargs):
        return SimpleNamespace(
            stdout="tested\n" if args[-1] == "version" else "--lolalign-multidomain",
            stderr="",
        )

    monkeypatch.setattr(structures.subprocess, "run", probe)
    assert structures.lolalign_version(executable="foldseek") == "tested"
    monkeypatch.setattr(
        structures.subprocess, "run", lambda **_: SimpleNamespace(stdout="old", stderr="")
    )
    with pytest.raises(ExternalToolError, match="does not support"):
        structures.lolalign_version(executable="foldseek")
    monkeypatch.undo()
    with pytest.raises(ExternalToolError, match="Cannot inspect"):
        structures.lolalign_version(executable="/no/such/foldseek")


def test_queued_structure_comparison_keeps_input_hashes(
    workbench_config,
    monomer,
    monkeypatch,
):
    """Run through the durable queue and keep correspondence and transform artifacts."""
    monkeypatch.setattr(structures, "lolalign_version", lambda **_: "tested")

    def binary(*, config, arguments, output):
        assert arguments[arguments.index("--alignment-type") + 1] == "3"
        assert arguments[arguments.index("--exhaustive-search") + 1] == "1"
        assert arguments[arguments.index("--lolalign-multidomain") + 1] == "0"
        (output / "foldseek_raw.tsv").write_text(VALID_LOL)

    monkeypatch.setattr(structures, "run_command", binary)
    job = submit_job(
        config=workbench_config,
        operation="structure_pair",
        parameters={"multidomain": False},
        uploads={"query_model": monomer, "target_model": monomer},
    )
    assert run_next(config=workbench_config) == job
    state = read_job(config=workbench_config, job_id=job)
    assert state["status"] == "COMPLETE", state
    path = artifact_path(config=workbench_config, job_id=job, name="structure_provenance.json")
    assert (
        json.loads(path.read_text())["models"]["query"]["sha256"]
        == hashlib.sha256(monomer).hexdigest()
    )
    assert "structural_correspondence.fasta" in state["artifacts"]
    assert "target_to_query.json" in state["artifacts"]


def test_batch_keeps_prior_rows_and_audits_modified_models(
    workbench_config,
    config_path,
    monomer,
    tmp_path,
    monkeypatch,
    capsys,
):
    """Relocate original models by prefix and compare only checksum-matching inputs."""
    models = tmp_path / "models"
    models.mkdir()
    digest = hashlib.sha256(monomer).hexdigest()
    (models / "a.pdb").write_bytes(monomer)
    (models / "b.pdb").write_bytes(monomer)
    (models / "c.pdb").write_bytes(monomer + b"REMARK changed\n")
    original = tmp_path / "original"
    original.mkdir()
    table = original / "structural_alignments.parquet"
    base = dict(
        cluster_id="OG1",
        reference_accession="a",
        mobile_accession="b",
        alignment_tool="TM-align",
        status="COMPLETE",
        reference_model_path="/previous/models/a.pdb",
        mobile_model_path="/previous/models/b.pdb",
        reference_model_sha256=digest,
        mobile_model_sha256=digest,
        aligned_length=4,
        rmsd_angstrom=0.1,
        minimum_tm_score=0.5,
    )
    rows = [
        base,
        {**base, "alignment_tool": "US-align"},
        {**base, "mobile_accession": "c", "mobile_model_path": "/previous/models/c.pdb"},
        {**base, "alignment_tool": "LoL-align"},
        {**base, "status": "REFERENCE"},
    ]
    pq.write_table(pa.Table.from_pylist(rows), table)
    before = table.read_bytes()
    calls = []
    executable = tmp_path / "foldseek"
    executable.write_text("not executed")
    config = replace(workbench_config, foldseek=str(executable))
    monkeypatch.setattr(structures, "lolalign_version", lambda **_: "tested")
    monkeypatch.setattr(
        structures, "compare_pair", lambda **kw: calls.append(kw) or {"status": "ALIGNED"}
    )
    assert (
        main(
            argv=[
                "--config",
                str(config_path),
                "--action",
                "plan-structures",
                "--pair-table",
                str(table),
                "--cluster-id",
                "OG1",
            ]
        )
        == 0
    )
    assert "original_rows\t3" in capsys.readouterr().out
    with pytest.raises(InputValidationError, match="exceeds"):
        structures.old_pairs(table=table, clusters=[], maximum=2)
    output = tmp_path / "new"
    kwargs = dict(
        config=config,
        table=table,
        models_root=models,
        source_model_root=Path("/previous/models"),
        output=output,
        clusters=["OG1"],
        maximum=3,
    )
    report = structures.compare_previous(**kwargs)
    assert len(calls) == 1
    assert "SHA-256 differs" in report.read_text()
    assert json.loads((output / "COMPARISON.json").read_text())["status_counts"] == {
        "ALIGNED": 2,
        "NOT_ASSESSED": 1,
    }
    assert table.read_bytes() == before
    with pytest.raises(InputValidationError, match="new output"):
        structures.compare_previous(**kwargs)
    kwargs["output"] = original / "do-not-write"
    with pytest.raises(InputValidationError, match="new output"):
        structures.compare_previous(**kwargs)


@pytest.mark.skipif(
    not (os.getenv("FOLDSEEK_EXE") and os.getenv("FOLDSEEK_REFERENCE_MODEL")),
    reason="Requires a LoL-enabled binary and a real single-chain reference model",
)
def test_real_foldseek_recovers_a_rigid_transform(tmp_path):
    """Confirm alignment and transform direction with a rotated known monomer."""
    original = Path(os.environ["FOLDSEEK_REFERENCE_MODEL"])
    reference = PDBParser(QUIET=True).get_structure("reference", str(original))
    for atom in reference.get_atoms():
        atom.transform(np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]), np.array([12, -5, 8]))
    mobile = tmp_path / "mobile.pdb"
    writer = PDBIO()
    writer.set_structure(reference)
    writer.save(str(mobile))
    output = tmp_path / "comparison"
    output.mkdir()
    config = structures.WorkbenchConfig(
        workspace=str(tmp_path), foldseek=os.environ["FOLDSEEK_EXE"], command_timeout=120
    )
    result = structures.compare_pair(
        config=config, query=original.read_bytes(), target=mobile.read_bytes(), output=output
    )
    assert result["status"] == "ALIGNED"
    assert result["query_model_coverage"] > 0.99 and result["rmsd_angstrom"] < 0.05
    transform = json.loads((output / "target_to_query.json").read_text())
    axes = ("x", "y", "z")
    _, _, original_records = structures.model_residues(content=original.read_bytes(), maximum=3000)
    _, _, mobile_records = structures.model_residues(content=mobile.read_bytes(), maximum=3000)
    expected = np.asarray([[row[axis] for axis in axes] for row in original_records])
    observed = np.asarray([[row[axis] for axis in axes] for row in mobile_records])
    recovered = observed @ np.asarray(transform["rotation"]).T + np.asarray(
        transform["translation"]
    )
    assert np.max(np.linalg.norm(expected - recovered, axis=1)) < 0.1
