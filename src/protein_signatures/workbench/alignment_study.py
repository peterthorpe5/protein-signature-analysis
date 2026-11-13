"""Exact protein alignment metrics and portable divergence trees."""

from __future__ import annotations

import io
import logging
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio import Phylo
from Bio.Align import substitution_matrices
from Bio.Phylo.NewickIO import NewickError

from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic, write_tsv_atomic

from .catalogue import RESIDUES, identifiers, open_catalogue
from .config import WorkbenchConfig

LOGGER = logging.getLogger(__name__)
CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")
PAIR_ROWS = 120
PAIR_CELLS = 80_000_000
PAIR_FIELDS = (
    "sequence_a",
    "sequence_b",
    "comparable_residue_count",
    "identical_count",
    "identity",
)


def read_alignment(*, path: Path, config: WorkbenchConfig) -> dict[str, str]:
    """Read a bounded rectangular aligned FASTA with distinct exact IDs."""
    records: dict[str, str] = {}
    width, length = 0, 0
    current, parts = None, []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    if not length or (width and length != width):
                        raise InputValidationError(
                            "Alignment rows must be non-empty and rectangular."
                        )
                    records[current] = "".join(parts)
                    width = length
                header = line[1:].split()
                if not header:
                    raise InputValidationError("Alignment contains an empty header.")
                current = identifiers(values=[header[0]])[0]
                if current in records or len(records) >= config.max_alignment_sequences:
                    raise InputValidationError(
                        "Alignment repeats IDs or exceeds the sequence limit."
                    )
                parts, length = [], 0
            else:
                residues = line.upper()
                if current is None or set(residues) - (RESIDUES | {"-"}):
                    raise InputValidationError(
                        "Alignment contains invalid protein residues or syntax."
                    )
                length += len(residues)
                if (len(records) + 1) * max(length, width) > config.max_alignment_cells:
                    raise InputValidationError("Alignment exceeds the configured cell limit.")
                parts.append(residues)
    if current is None or not length or (width and width != length):
        raise InputValidationError("Alignment rows must be non-empty and rectangular.")
    records[current] = "".join(parts)
    if len(records) < 2 or any(not sequence.strip("-") for sequence in records.values()):
        raise InputValidationError("Alignment requires at least two non-empty protein rows.")
    return records


def pairwise_identity(*, records: dict[str, str]) -> list[dict]:
    """Report comparable and identical canonical residues for every unordered pair."""
    count = len(records)
    if not count or count > PAIR_ROWS:
        raise InputValidationError(f"Select 1–{PAIR_ROWS} aligned rows for the identity matrix.")
    labels = tuple(records)
    width = len(records[labels[0]])
    if not width or any(len(row) != width for row in records.values()):
        raise InputValidationError("Pairwise identity needs a rectangular protein alignment.")
    if count * (count + 1) // 2 * width > PAIR_CELLS:
        raise InputValidationError(
            "Identity matrix exceeds the calculation limit. Select fewer rows."
        )
    matrix = np.frombuffer("".join(records.values()).encode("ascii"), dtype=np.uint8).reshape(
        count, width
    )
    known = np.isin(matrix, [ord(value) for value in CANONICAL])
    result = []
    for i, first in enumerate(labels):
        comparable = known[i:] & known[i]
        shared = np.count_nonzero(comparable, axis=1)
        identical = np.count_nonzero((matrix[i:] == matrix[i]) & comparable, axis=1)
        for second, total, equal in zip(labels[i:], shared, identical, strict=True):
            result.append(
                dict(
                    zip(
                        PAIR_FIELDS,
                        (
                            first,
                            second,
                            int(total),
                            int(equal),
                            int(equal) / int(total) if total else None,
                        ),
                        strict=True,
                    )
                )
            )
    return result


def write_alignment_study(*, records: dict[str, str], reference_id: str, output: Path) -> None:
    """Write column and reference metrics, plus all-pair counts when within bounds."""
    if reference_id not in records:
        raise InputValidationError("Choose a reference protein present in the alignment.")
    reference = records[reference_id]
    total, width = len(records), len(reference)
    if not width or any(len(sequence) != width for sequence in records.values()):
        raise InputValidationError("Alignment study requires equal non-empty row lengths.")
    columns, members, position = [], [], 0
    for idx, residues in enumerate(zip(*records.values(), strict=True)):
        counts = Counter(value for value in residues if value in CANONICAL)
        known = sum(counts.values())
        consensus = min(counts, key=lambda aa: (-counts[aa], aa)) if counts else None
        position += reference[idx] != "-"
        columns.append(
            {
                "alignment_column": idx + 1,
                "reference_position": position if reference[idx] != "-" else None,
                "gap_fraction": residues.count("-") / total,
                "canonical_fraction": known / total,
                "consensus_residue": consensus,
                "consensus_fraction_of_known": counts[consensus] / known if known else None,
                "shannon_entropy_bits": (
                    -sum((n / known) * math.log2(n / known) for n in counts.values())
                    if known
                    else None
                ),
            }
        )
    matrix = substitution_matrices.load(name="BLOSUM62")
    for identifier, sequence in records.items():
        overlap, shared, identical, score = 0, 0, 0, 0.0
        for first, second in zip(reference, sequence, strict=True):
            overlap += first != "-" and second != "-"
            if first in CANONICAL and second in CANONICAL:
                shared += 1
                identical += first == second
                score += float(matrix[first, second])
        ref_length, seq_length = len(reference.replace("-", "")), len(sequence.replace("-", ""))
        members.append(
            {
                "sequence_id": identifier,
                "reference_id": reference_id,
                "sequence_length": seq_length,
                "comparable_residue_count": shared,
                "identity_to_reference_fraction": identical / shared if shared else None,
                "reference_overlap_fraction": overlap / ref_length if ref_length else None,
                "member_overlap_fraction": overlap / seq_length if seq_length else None,
                "blosum62_sum_over_comparable_positions": score,
            }
        )
    for name, rows in (("alignment_columns", columns), ("alignment_members", members)):
        write_tsv_atomic(path=output / f"{name}.tsv", fieldnames=tuple(rows[0]), records=rows)
    within_pair_limit = total <= PAIR_ROWS and total * (total + 1) // 2 * width <= PAIR_CELLS
    if within_pair_limit:
        write_tsv_atomic(
            path=output / "pairwise_identity.tsv",
            fieldnames=PAIR_FIELDS,
            records=pairwise_identity(records=records),
        )
    write_json_atomic(
        path=output / "alignment_summary.json",
        value={
            "reference_id": reference_id,
            "aligned_sequence_count": total,
            "alignment_length": width,
            "mean_gap_fraction": sum(row["gap_fraction"] for row in columns) / width,
            "pairwise_table": within_pair_limit,
            "identity_definition": (
                "Identical canonical residues divided by positions canonical in both rows."
            ),
            "overlap_definition": "Positions non-gap in both rows / ungapped length of one row.",
            "blosum_definition": "Raw BLOSUM62 sum over comparable residues, no gap penalty.",
        },
    )


def auspice_dataset(
    *,
    newick: str,
    labels: dict[str, str],
    context: dict[str, dict] | None = None,
) -> dict:
    """Validate inferred tip identities and export cumulative divergence without dates."""
    try:
        tree = Phylo.read(io.StringIO(newick), "newick")
        tips = [tip.name for tip in tree.get_terminals()]
    except (ValueError, NewickError, RecursionError) as error:
        raise InputValidationError("FastTree did not return a valid Newick tree.") from error
    if (
        len(tips) != len(labels)
        or set(tips) != set(labels)
        or len(set(labels.values())) != len(labels)
    ):
        raise InputValidationError("FastTree tips differ from the exact aligned proteins.")
    result: dict = {}
    used = set(labels.values())
    stack = [(tree.root, result, 0.0, 0)]
    node_index = 0
    traits = {
        "released_species": "Released species",
        "reviewed_species": "Reviewed species",
        "reviewed_taxon_id": "Reviewed NCBI taxon",
        "source_group": "Source group",
    }
    while stack:
        clade, node, distance, depth = stack.pop()
        length = clade.branch_length or 0.0
        if depth > 300 or not math.isfinite(length) or length < 0:
            raise InputValidationError("Tree nesting or branch lengths are invalid.")
        distance += length
        node_index += 1
        node_name = labels[clade.name] if clade.is_terminal() else f"INTERNAL_{node_index}"
        while not clade.is_terminal() and node_name in used:
            node_name += "_"
        used.add(node_name)
        attributes = {"div": distance}
        if clade.confidence is not None and math.isfinite(clade.confidence):
            attributes["local_support"] = {"value": clade.confidence}
        node.update(name=node_name, node_attrs=attributes)
        if clade.is_terminal():
            row = (context or {}).get(node_name, {})
            for field in traits:
                if row.get(field):
                    attributes[field] = {"value": str(row[field])}
        else:
            node["children"] = [{} for _ in clade.clades]
            stack.extend(
                (child, target, distance, depth + 1)
                for child, target in zip(clade.clades, node["children"], strict=True)
            )
    return {
        "version": "v2",
        "meta": {
            "title": "Protein family phylogeny",
            "updated": datetime.now(tz=timezone.utc).date().isoformat(),
            "description": "WAG+gamma protein divergence with a display root.",
            "panels": ["tree"],
            "colorings": [
                {"key": field, "title": title, "type": "categorical"}
                for field, title in traits.items()
            ],
            "filters": list(traits),
            "display_defaults": {"distance_measure": "div", "layout": "rect"},
        },
        "tree": result,
    }


def infer_aligned_tree(
    *,
    config: WorkbenchConfig,
    alignment: Path,
    labels: dict[str, str],
    output: Path,
    context: Path | None = None,
) -> None:
    """Run FastTree directly on a saved matrix and export validated Newick and Auspice."""
    from .analysis import run_command
    from .context import open_context

    if len(labels) < 3:
        raise InputValidationError("Tree inference requires at least three sequences.")
    metadata = {}
    if context is not None:
        with open_catalogue(config=config, with_index=False) as connection:
            open_context(connection=connection, path=context)
            metadata = {
                row["sequence_id"]: row
                for row in connection.execute(
                    "SELECT * FROM family_context WHERE sequence_id IN (SELECT unnest(?))",
                    [list(labels.values())],
                )
                .fetchdf()
                .fillna("")
                .to_dict(orient="records")
            }
    run_command(
        config=config,
        output=output,
        stdout_name="tree.nwk",
        arguments=[
            config.fasttree,
            "-wag",
            "-gamma",
            str(alignment),
        ],
    )
    newick = (output / "tree.nwk").read_text(encoding="utf-8").strip()
    if not newick.endswith(";"):
        raise InputValidationError("FastTree did not return a complete Newick tree.")
    write_json_atomic(
        path=output / "auspice.json",
        value=auspice_dataset(
            newick=newick,
            labels=labels,
            context=metadata,
        ),
    )
    LOGGER.info("Verified tree with %d exact protein tips", len(labels))


def tree_from_saved(
    *,
    config: WorkbenchConfig,
    alignment: Path,
    output: Path,
    context: Path | None = None,
) -> None:
    """Infer a tree from a checked matrix, retaining original IDs and exact ungapped rows."""
    from .context import CONTEXT_FILE, carry_context

    records = read_alignment(path=alignment, config=config)
    labels = {f"seq_{idx:06d}": identifier for idx, identifier in enumerate(records)}
    internal = output / "alignment.internal.fasta"
    internal.write_text(
        data="".join(f">{tip}\n{records[identifier]}\n" for tip, identifier in labels.items()),
        encoding="utf-8",
    )
    sequences = output / "sequences.fasta"
    sequences.write_text(
        data="".join(f">{key}\n{row.replace('-', '')}\n" for key, row in records.items()),
        encoding="utf-8",
    )
    write_tsv_atomic(
        path=output / "tip_labels.tsv",
        fieldnames=("tree_tip", "sequence_id"),
        records=(
            {"tree_tip": tip, "sequence_id": identifier} for tip, identifier in labels.items()
        ),
    )
    if context is not None:
        carry_context(config=config, context=context, sequences=sequences, output=output)
    infer_aligned_tree(
        config=config,
        alignment=internal,
        labels=labels,
        output=output,
        context=output / CONTEXT_FILE if context is not None else None,
    )
