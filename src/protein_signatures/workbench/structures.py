"""Local LoL-align comparisons and read-only reruns of E3 structural pairs.

Foldseek score and transform conventions follow the MIT-licensed E3 implementation
(Peter Thorpe, E3_project_draft, e3structalign/lolalign.py).
"""

from __future__ import annotations

import hashlib
import io
import logging
import math
import re
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyarrow.parquet as pq
from Bio.PDB import FastMMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.SeqUtils import seq1

from protein_signatures.checksums import sha256_file
from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.io_utils import write_json_atomic, write_tsv_atomic

from .analysis import run_command
from .config import WorkbenchConfig, bounded_number

LOGGER = logging.getLogger(__name__)
MODEL_BYTES = 20 * 1024 * 1024
LOL_FIELDS = (
    "query",
    "target",
    "qlen",
    "tlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "qaln",
    "taln",
    "qtmscore",
    "ttmscore",
    "rmsd",
    "evalue",
    "prob",
    "u",
    "t",
)
OLD_FIELDS = (
    "cluster_id",
    "reference_accession",
    "mobile_accession",
    "alignment_tool",
    "status",
    "reference_model_path",
    "mobile_model_path",
    "reference_model_sha256",
    "mobile_model_sha256",
    "aligned_length",
    "rmsd_angstrom",
    "minimum_tm_score",
)
METRIC_FIELDS = (
    "status",
    "alignment_scope",
    "aligned_residue_pairs",
    "alignment_columns",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "query_model_coverage",
    "target_model_coverage",
    "sequence_identity",
    "query_tm_score",
    "target_tm_score",
    "rmsd_angstrom",
    "lol_normalised_score",
    "lol_homology_probability",
)


def read_model(*, path: Path) -> bytes:
    """Read a bounded coordinate snapshot without loading oversized artifacts."""
    with path.open(mode="rb") as handle:
        content = handle.read(MODEL_BYTES + 1)
    if not 0 < len(content) <= MODEL_BYTES:
        raise InputValidationError("Coordinate files must contain 1 byte to 20 MiB.")
    return content


def model_residues(*, content: bytes, maximum: int) -> tuple[str, str, list[dict]]:
    """Parse one monomer model and return its coordinate sequence and C-alpha records.

    Args:
        content: Bounded PDB or mmCIF coordinate bytes.
        maximum: Maximum C-alpha residues; no truncation or chain selection.

    Returns:
        File suffix, coordinate sequence and residue records.

    Raises:
        InputValidationError: If the model or coordinates cannot be interpreted unambiguously.
    """
    if not 0 < len(content) <= MODEL_BYTES:
        raise InputValidationError("Each coordinate model must be at most 20 MiB.")
    try:
        text = content.decode("utf-8")
        cif = text.lstrip().startswith("data_")
        parser = FastMMCIFParser(QUIET=True) if cif else PDBParser(QUIET=True)
        structure = parser.get_structure("model", io.StringIO(text))
        models = list(structure.get_models())
        if len(models) != 1:
            raise InputValidationError("Provide one coordinate model, with one protein chain.")
        chains = [
            (
                chain.id,
                [
                    residue
                    for residue in chain
                    if "CA" in residue and (residue.id[0] == " " or residue.resname == "MSE")
                ],
            )
            for chain in models[0]
        ]
        chains = [(name, residues) for name, residues in chains if residues]
        if len(chains) != 1:
            raise InputValidationError("Provide one selected protein chain with C-alpha atoms.")
        name, residues = chains[0]
        if not 1 <= len(residues) <= maximum:
            raise InputValidationError(f"Coordinate model exceeds the {maximum} residue limit.")
        records = []
        for index, residue in enumerate(residues, start=1):
            atom = residue["CA"]
            x, y, z = (float(value) for value in atom.coord)
            bfactor = float(atom.bfactor)
            if not all(math.isfinite(value) for value in (x, y, z, bfactor)):
                raise InputValidationError("Coordinate model contains non-finite values.")
            records.append(
                {
                    "position": index,
                    "chain": name,
                    "model_residue": f"{residue.id[1]}{residue.id[2].strip()}",
                    "residue": seq1(residue.resname, custom_map={"MSE": "M"}),
                    "x": x,
                    "y": y,
                    "z": z,
                    "bfactor": bfactor,
                }
            )
    except (
        UnicodeError,
        ValueError,
        KeyError,
        IndexError,
        TypeError,
        PDBConstructionException,
        StopIteration,
    ) as error:
        raise InputValidationError("Could not parse the supplied coordinate model.") from error
    return ".cif" if cif else ".pdb", "".join(row["residue"] for row in records), records


def lolalign_version(*, executable: str) -> str:
    """Check the actual installed Foldseek binary for LoL-align capability.

    Args:
        executable: Operator-configured binary, or its name on PATH.

    Returns:
        Exact version or revision reported by the tool.

    Raises:
        ExternalToolError: If an older or missing binary cannot run LoL-align.
    """
    try:
        version = subprocess.run(
            args=[executable, "version"], capture_output=True, text=True, timeout=15, check=True
        ).stdout.strip()
        help_result = subprocess.run(
            args=[executable, "easy-search", "-h"],
            capture_output=True,
            text=True,
            timeout=15,
            check=True,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise ExternalToolError("Cannot inspect the configured Foldseek executable.") from error
    if (
        not version
        or len(version.splitlines()) != 1
        or "--lolalign-multidomain" not in help_result.stdout + help_result.stderr
    ):
        raise ExternalToolError(
            "This Foldseek build does not support LoL-align; see docs/STRUCTURAL_RERUN.md."
        )
    return version


def parse_lolalign(*, text: str, query: str, target: str) -> tuple[dict, dict]:
    """Validate one LoL-align correspondence, geometry and local scores.

    Args:
        text: Exact headerless output for the requested pair.
        query: Full coordinate sequence of the reference model.
        target: Full coordinate sequence of the mobile model.

    Returns:
        Result row and validated target-to-query rigid transform, or explicit no hit.

    Raises:
        ExternalToolError: If the tool reports incompatible or ambiguous evidence.
    """
    rows = [line for line in text.splitlines() if line.strip()]
    if not rows:
        return {"status": "NO_ALIGNMENT_FOUND", "alignment_scope": "LOCAL"}, {}
    if len(rows) != 1 or len(rows[0].split("\t")) != len(LOL_FIELDS):
        raise ExternalToolError("Expected one LoL-align row in the requested format.")
    row = dict(zip(LOL_FIELDS, rows[0].split("\t"), strict=True))
    try:
        qstart, qend, tstart, tend = (
            int(row[name]) for name in ("qstart", "qend", "tstart", "tend")
        )
        qa, ta = row["qaln"].upper(), row["taln"].upper()
        if (
            Path(row["query"]).stem != "query"
            or Path(row["target"]).stem != "target"
            or int(row["qlen"]) != len(query)
            or int(row["tlen"]) != len(target)
            or not 1 <= qstart <= qend <= len(query)
            or not 1 <= tstart <= tend <= len(target)
            or len(qa) != len(ta)
            or qa.replace("-", "") != query[qstart - 1 : qend]
            or ta.replace("-", "") != target[tstart - 1 : tend]
        ):
            raise ValueError("Model identifiers, lengths or coordinates differ")
        aligned = [(q, t) for q, t in zip(qa, ta, strict=True) if q != "-" and t != "-"]
        if not aligned or any(q == t == "-" for q, t in zip(qa, ta, strict=True)):
            raise ValueError("Invalid alignment gap pattern")
        qt, tt, rmsd, score, probability = (
            float(row[key]) for key in ("qtmscore", "ttmscore", "rmsd", "evalue", "prob")
        )
        rotation = np.array([float(value) for value in row["u"].split(",")]).reshape(3, 3)
        translation = np.array([float(value) for value in row["t"].split(",")])
        if (
            not all(math.isfinite(value) for value in (qt, tt, rmsd, score, probability))
            or not all(0 <= value <= 1 for value in (qt, tt, probability))
            or rmsd < 0
            or score > 1.00001
            or translation.shape != (3,)
            or not np.isfinite(translation).all()
            or not np.allclose(rotation @ rotation.T, np.eye(3), atol=0.005)
            or not math.isclose(float(np.linalg.det(rotation)), 1.0, abs_tol=0.01)
        ):
            raise ValueError("Invalid scores or rigid transform")
    except (ValueError, OverflowError) as error:
        raise ExternalToolError(f"Invalid LoL-align correspondence: {error}.") from error
    count = len(aligned)
    return {
        "status": "ALIGNED",
        "alignment_scope": "LOCAL",
        "aligned_residue_pairs": count,
        "alignment_columns": len(qa),
        "query_start": qstart,
        "query_end": qend,
        "target_start": tstart,
        "target_end": tend,
        "query_model_coverage": count / len(query),
        "target_model_coverage": count / len(target),
        "sequence_identity": sum(q == t for q, t in aligned) / count,
        "query_tm_score": qt,
        "target_tm_score": tt,
        "rmsd_angstrom": rmsd,
        "lol_normalised_score": score,
        "lol_homology_probability": probability,
    }, {
        "direction": "target_to_query",
        "convention": "query_xyz = rotation @ target_xyz + translation",
        "rotation": rotation.tolist(),
        "translation": translation.tolist(),
    }


def compare_pair(
    *,
    config: WorkbenchConfig,
    query: bytes,
    target: bytes,
    output: Path,
    multidomain: bool = True,
    version: str = "",
    retain_models: bool = True,
) -> dict:
    """Run an exhaustive local comparison using a new output directory.

    Args:
        config: Worker limits and Foldseek executable.
        query: Reference coordinates, in PDB or mmCIF.
        target: Mobile coordinates, in PDB or mmCIF.
        output: Empty owned destination.
        multidomain: True for multi-domain LoL scoring, false for single-domain.
        version: Already checked version for an operator-controlled batch.
        retain_models: Retain coordinate snapshots for interactive pair viewing.

    Returns:
        Per-pair scientific evidence and status.

    Raises:
        InputValidationError: For unsupported inputs or options.
        ExternalToolError: For an unavailable tool or inconsistent result.
    """
    if not isinstance(multidomain, bool):
        raise InputValidationError("multidomain must be true or false.")
    checked_version = version or lolalign_version(executable=config.foldseek)
    paths, sequences, provenance = {}, {}, {}
    raw = output / "foldseek_raw.tsv"
    with TemporaryDirectory(prefix="lolalign_", dir=output) as temporary:
        model_dir = output if retain_models else Path(temporary)
        for label, content in (("query", query), ("target", target)):
            suffix, sequence, _ = model_residues(
                content=content, maximum=config.max_structure_residues
            )
            path = model_dir / f"{label}{suffix}"
            with path.open(mode="xb") as handle:
                handle.write(content)
            paths[label], sequences[label] = path, sequence
            provenance[label] = {
                "sha256": hashlib.sha256(content).hexdigest(),
                "residues": len(sequence),
            }
        run_command(
            config=config,
            output=output,
            arguments=[
                config.foldseek,
                "easy-search",
                str(paths["query"]),
                str(paths["target"]),
                str(raw),
                temporary,
                "--alignment-type",
                "3",
                "--lolalign-multidomain",
                str(int(multidomain)),
                "--exhaustive-search",
                "1",
                "--threads",
                "1",
                "--max-seqs",
                "1",
                "-e",
                "1e30",
                "-a",
                "1",
                "--format-output",
                ",".join(LOL_FIELDS),
            ],
        )
    if not raw.is_file() or raw.stat().st_size > 4 * MODEL_BYTES:
        raise ExternalToolError("Missing or oversized Foldseek pairwise result.")
    text = raw.read_text(encoding="utf-8")
    metrics, transform = parse_lolalign(
        text=text, query=sequences["query"], target=sequences["target"]
    )
    if transform:
        write_json_atomic(path=output / "target_to_query.json", value=transform)
        raw_row = dict(zip(LOL_FIELDS, text.strip().split("\t"), strict=True))
        (output / "structural_correspondence.fasta").write_text(
            data=f">query\n{raw_row['qaln']}\n>target\n{raw_row['taln']}\n", encoding="utf-8"
        )
    write_tsv_atomic(
        path=output / "lolalign.tsv",
        fieldnames=METRIC_FIELDS,
        records=[{name: metrics.get(name) for name in METRIC_FIELDS}],
    )
    write_json_atomic(
        path=output / "structure_provenance.json",
        value={
            "tool": "LoL-align via Foldseek",
            "version": checked_version,
            "models": provenance,
            "multidomain": multidomain,
            "scope": "LOCAL",
            "score_note": "evalue slot contains normalised LoL score, not a statistical E-value",
            "tm_note": "Geometric TM metrics for the LoL residue correspondence",
        },
    )
    LOGGER.info("LoL-align pair status: %s", metrics["status"])
    return metrics


def old_pairs(*, table: Path, clusters: list[str], maximum: int) -> list[dict]:
    """Read a bounded exact selection of original E3 US-align/TM-align rows.

    Args:
        table: Original structural_alignments.parquet.
        clusters: Exact groups to include; empty means every group.
        maximum: Maximum selected comparison rows; no silent sampling.

    Returns:
        Prior comparison rows with model paths, hashes and original metrics.
    """
    maximum = int(
        bounded_number(value=maximum, name="max-pairs", minimum=1, maximum=100_000, integer=True)
    )
    selected = []
    with pq.ParquetFile(
        source=table, thrift_string_size_limit=16_000_000, thrift_container_size_limit=1_000_000
    ) as reader:
        if not set(OLD_FIELDS).issubset(reader.schema_arrow.names):
            raise InputValidationError("Select the original E3 structural_alignments.parquet.")
        for batch in reader.iter_batches(batch_size=2048, columns=OLD_FIELDS):
            for row in batch.to_pylist():
                if str(row["alignment_tool"]).casefold() not in {"us-align", "tm-align"}:
                    continue
                if row["status"] == "REFERENCE" or (clusters and row["cluster_id"] not in clusters):
                    continue
                if row["reference_accession"] == row["mobile_accession"]:
                    raise InputValidationError("Unexpected self-pair in the prior alignment table.")
                selected.append(row)
                if len(selected) > maximum:
                    raise InputValidationError(
                        "Pair count exceeds --max-pairs; select exact groups."
                    )
    if not selected:
        raise InputValidationError("No original US-align or TM-align pairs match the selection.")
    return selected


def original_model(*, path: str, digest: str, old_root: Path, root: Path) -> bytes:
    """Read the exact previously assessed model from a trusted relocated tree.

    Args:
        path: Prior model path recorded in the structural table.
        digest: Required original SHA-256.
        old_root: Prefix used by the previous producer.
        root: Existing model root on the current server.

    Returns:
        Original coordinate bytes.

    Raises:
        InputValidationError: If the model is missing, unsafe or no longer identical.
    """
    if (
        not isinstance(path, str)
        or not path
        or not isinstance(digest, str)
        or not re.fullmatch(r"[a-f0-9]{64}", digest)
    ):
        raise InputValidationError("Previous model path or SHA-256 is missing.")
    try:
        relative = Path(path).relative_to(old_root)
        source = (root / relative).resolve()
        if not source.is_relative_to(root) or not source.is_file():
            raise ValueError("Model is missing or outside the declared tree")
        data = read_model(path=source)
    except (OSError, ValueError) as error:
        raise InputValidationError(
            "Previous coordinate model is missing or outside its root."
        ) from error
    if hashlib.sha256(data).hexdigest() != digest:
        raise InputValidationError("Model size or SHA-256 differs from the previous run.")
    return data


def compare_previous(
    *,
    config: WorkbenchConfig,
    table: Path,
    models_root: Path,
    source_model_root: Path,
    output: Path,
    clusters: list[str],
    maximum: int,
    multidomain: bool = True,
) -> Path:
    """Compare prior structural pairs and retain both methods in a new audit table.

    Args:
        config: Trusted execution limits and completed reference location.
        table: Read-only original alignment table.
        models_root: Current location of the original coordinate tree.
        source_model_root: Absolute prefix stored in the prior model paths.
        output: New directory outside both original input trees.
        clusters: Exact cluster IDs, or empty for all.
        maximum: Explicit maximum prior rows to assess.
        multidomain: LoL-align scoring mode.

    Returns:
        Path to the separately published comparison report.
    """
    table, models_root, output = (
        path.expanduser().resolve() for path in (table, models_root, output)
    )
    if not source_model_root.is_absolute() or not models_root.is_dir():
        raise InputValidationError("Supply the absolute original prefix and current model root.")
    protected = [
        table.parent.parent if table.parent.name == "tables" else table.parent,
        models_root,
    ]
    if config.resource:
        protected.append(Path(config.resource).resolve())
    if output.exists() or any(output.is_relative_to(path) for path in protected):
        raise InputValidationError(
            "Choose a new output directory outside the original results and models."
        )
    before = sha256_file(path=table)
    rows = old_pairs(table=table, clusters=clusters, maximum=maximum)
    executable = shutil.which(config.foldseek) or config.foldseek
    version = lolalign_version(executable=executable)
    output.mkdir(parents=True, mode=0o700)
    results, cache = [], {}
    metric_columns = {
        name: name if name.startswith("lol_") else f"lol_{name}" for name in METRIC_FIELDS
    }
    for index, old in enumerate(rows, start=1):
        result = {f"previous_{key}": value for key, value in old.items()}
        key = (old["reference_model_sha256"], old["mobile_model_sha256"])
        pair = output / f"pair_{index:06d}"
        try:
            query, target = (
                original_model(
                    path=old[f"{name}_model_path"],
                    digest=old[f"{name}_model_sha256"],
                    old_root=source_model_root,
                    root=models_root,
                )
                for name in ("reference", "mobile")
            )
            if key not in cache:
                pair.mkdir(mode=0o700)
                metrics = compare_pair(
                    config=config,
                    query=query,
                    target=target,
                    output=pair,
                    multidomain=multidomain,
                    version=version,
                    retain_models=False,
                )
                cache[key] = pair.name, metrics
            result["pair_artifacts"], metrics = cache[key]
            result.update({metric_columns[name]: value for name, value in metrics.items()})
        except (InputValidationError, ExternalToolError, OSError) as error:
            result["lol_status"] = (
                "NOT_ASSESSED" if isinstance(error, InputValidationError) else "FAILED"
            )
            result["detail"] = str(error)
            result["pair_artifacts"] = pair.name if pair.is_dir() else ""
            LOGGER.warning("Pair %d: %s", index, error)
        results.append(result)
    if sha256_file(path=table) != before:
        raise InputValidationError("Prior alignment table changed during comparison.")
    report = output / "structural_method_comparison.tsv"
    columns = (
        *[f"previous_{name}" for name in OLD_FIELDS],
        *metric_columns.values(),
        "pair_artifacts",
        "detail",
    )
    write_tsv_atomic(
        path=report,
        records=({name: row.get(name) for name in columns} for row in results),
        fieldnames=columns,
    )
    write_json_atomic(
        path=output / "COMPARISON.json",
        value={
            "original_table": str(table),
            "original_table_sha256": before,
            "models_root": str(models_root),
            "source_model_root": str(source_model_root),
            "foldseek_version": version,
            "foldseek_sha256": sha256_file(path=Path(executable)),
            "multidomain": multidomain,
            "selected_clusters": clusters,
            "rows": len(results),
            "comparison_sha256": sha256_file(path=report),
            "status_counts": {
                status: sum(row["lol_status"] == status for row in results)
                for status in sorted({row["lol_status"] for row in results})
            },
            "interpretation": "Exploratory local comparison; original evidence remains separate.",
        },
    )
    return report
