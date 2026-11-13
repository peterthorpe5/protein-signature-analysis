"""Reviewed species backbones for MSA display, kept distinct from inferred protein trees."""

from __future__ import annotations

import logging
import math
from io import StringIO
from typing import Mapping, Sequence

from Bio import Phylo
from Bio.Phylo.Newick import Clade
from Bio.Phylo.NewickIO import NewickError

from protein_signatures.errors import InputValidationError

LOGGER = logging.getLogger(__name__)


def species_display_tree(
    *,
    newick: str,
    sequence_ids: Sequence[str],
    metadata: Mapping[str, Mapping[str, object]],
) -> tuple[str, list[dict[str, object]]]:
    """Prune a reviewed species tree and group protein rows at exact species tips for JBrowse.

    Args:
        newick: Reviewed Newick with exact species labels or taxon IDs as unique tip names.
        sequence_ids: All aligned protein IDs, including multiple proteins from one species.
        metadata: Exact sequence_id to reviewed species/taxon metadata.

    Returns:
        A protein-tip display tree and mapping audit. Within-species proteins form a
        zero-length display polytomy; no gene relationships or dates are inferred.

    Raises:
        InputValidationError: If a protein lacks an exact tree mapping or the tree fails.
    """
    if not 1 <= len(sequence_ids) <= 5_000 or len(set(sequence_ids)) != len(sequence_ids):
        raise InputValidationError(
            "Species-tree display requires 1 to 5000 unique aligned protein IDs"
        )
    if not newick.strip().endswith(";") or len(newick.encode()) > 2_000_000:
        raise InputValidationError("Supply a complete reviewed species Newick tree of at most 2 MB")
    try:
        tree = Phylo.read(file=StringIO(newick), format="newick")
    except (ValueError, TypeError, NewickError, RecursionError) as exc:
        raise InputValidationError("Could not parse the reviewed species Newick tree") from exc
    tips = tree.get_terminals()
    labels = [tip.name for tip in tips]
    if any(not name for name in labels) or len(set(labels)) != len(labels):
        raise InputValidationError("Species-tree tip labels must be non-empty and unique")
    if any(
        clade.branch_length is not None
        and (not math.isfinite(clade.branch_length) or clade.branch_length < 0)
        for clade in tree.find_clades()
    ):
        raise InputValidationError("Species-tree branch lengths must be finite and non-negative")
    label_set = set(labels)
    grouped, audit = {}, []
    for identifier in sequence_ids:
        record = metadata.get(identifier, {})
        aliases = {
            str(record.get("released_species", "")),
            str(record.get("reviewed_species", "")),
            str(record.get("reviewed_taxon_id", "")),
        } - {"", "None"}
        matches = label_set & aliases
        if len(matches) != 1:
            raise InputValidationError(
                f"Protein {identifier} needs one exact reviewed species/taxon match in the tree"
            )
        label = next(iter(matches))
        grouped.setdefault(label, []).append(identifier)
        audit.append(
            {
                "sequence_id": identifier,
                "species_tree_tip": label,
                "released_species": record.get("released_species", ""),
                "reviewed_species": record.get("reviewed_species", ""),
                "reviewed_taxon_id": record.get("reviewed_taxon_id", ""),
                "interpretation": "Species backbone; protein grouping is display only",
            }
        )
    for tip in tips:
        if tip.name not in grouped:
            tree.prune(target=tip)
    for tip in tree.get_terminals():
        identifiers = grouped[tip.name]
        tip.name = "species_" + tip.name
        tip.clades = [Clade(name=identifier, branch_length=0.0) for identifier in identifiers]
    final = [tip.name for tip in tree.get_terminals()]
    if set(final) != set(sequence_ids) or len(final) != len(sequence_ids):
        raise InputValidationError(
            "Species display tree did not retain exactly the aligned protein IDs"
        )
    output = StringIO()
    Phylo.write(trees=tree, file=output, format="newick")
    LOGGER.info(
        "Built species-backbone display for %d proteins and %d species",
        len(sequence_ids),
        len(grouped),
    )
    return output.getvalue(), audit
