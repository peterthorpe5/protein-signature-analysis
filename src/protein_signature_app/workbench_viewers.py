"""Bounded alignment and phylogram viewers for completed workbench jobs."""

from __future__ import annotations

import io
import math

import plotly.graph_objects as go
from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError

from protein_signatures.errors import InputValidationError


def alignment_window(*, text: str, start: int, width: int = 80, maximum_rows: int = 50) -> str:
    """Render an explicit cropped alignment view while preserving complete downloads.

    Args:
        text: Complete bounded alignment FASTA.
        start: First displayed alignment position, one-based.
        width: Maximum displayed columns, between 1 and 200.
        maximum_rows: Maximum displayed sequences, between 1 and 100.

    Returns:
        Monospaced alignment window with exact IDs.

    Raises:
        InputValidationError: If the alignment or window controls are invalid.
    """
    if (
        isinstance(start, bool)
        or not isinstance(start, int)
        or start < 1
        or not 1 <= width <= 200
        or not 1 <= maximum_rows <= 100
        or len(text) > 20_000_000
    ):
        raise InputValidationError("Invalid alignment display bounds.")
    records, current = {}, None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:]
            if not current or current in records:
                raise InputValidationError("Invalid alignment header.")
            records[current] = ""
        elif current is None:
            raise InputValidationError("Sequence precedes an alignment header.")
        else:
            records[current] += line
    lengths = {len(value) for value in records.values()}
    if len(lengths) != 1 or not records or start > next(iter(lengths)):
        raise InputValidationError(
            "Alignment is empty, unequal-length or the window starts past its end."
        )
    labels = list(records)[:maximum_rows]
    padding = min(60, max(len(label) for label in labels))
    return "\n".join(
        f"{label[:padding]:<{padding}}  {records[label][start - 1 : start - 1 + width]}"
        for label in labels
    )


def tree_figure(*, newick: str, labels: dict[str, str]) -> go.Figure:
    """Create an interactive rooted-display phylogram without changing tree topology.

    Args:
        newick: Complete FastTree Newick with safe internal tip names.
        labels: Complete safe-tip-to-original-ID map.

    Returns:
        Plotly branch plot with original protein IDs and hover text.

    Raises:
        InputValidationError: If tree tips, labels or branch lengths are invalid.
    """
    if len(newick) > 5_000_000 or not newick.strip().endswith(";"):
        raise InputValidationError("A complete bounded Newick tree is required.")
    try:
        tree = Phylo.read(file=io.StringIO(initial_value=newick), format="newick")
    except (ValueError, RecursionError, NewickError) as error:
        raise InputValidationError("Could not parse the completed Newick tree.") from error
    tips = tree.get_terminals()
    if (
        not 3 <= len(tips) <= 1000
        or len({tip.name for tip in tips}) != len(tips)
        or set(labels) != {tip.name for tip in tips}
    ):
        raise InputValidationError("Tree tips do not match the complete protein label map.")
    depths = tree.depths()
    if any(not math.isfinite(value) or value < 0 for value in depths.values()) or any(
        clade.branch_length is not None
        and (not math.isfinite(clade.branch_length) or clade.branch_length < 0)
        for clade in tree.find_clades()
    ):
        raise InputValidationError("Tree contains invalid branch lengths.")
    # A zero-length tree is displayed cladographically, with this stated on the axis.
    axis = "Inferred divergence (substitutions per site)"
    if max(depths.values()) == 0:
        depths = tree.depths(unit_branch_lengths=True)
        axis = "Branch depth (all inferred branch lengths were zero)"
    positions = {tip: index for index, tip in enumerate(tips)}
    x, y = [], []
    for clade in tree.find_clades(order="postorder"):
        if clade.clades:
            positions[clade] = sum(positions[child] for child in clade.clades) / len(clade.clades)
            child_y = [positions[child] for child in clade.clades]
            x.extend([depths[clade], depths[clade], None])
            y.extend([min(child_y), max(child_y), None])
            for child in clade.clades:
                x.extend([depths[clade], depths[child], None])
                y.extend([positions[child], positions[child], None])
    figure = go.Figure(
        data=[
            go.Scatter(
                x=x, y=y, mode="lines", line={"color": "#1768ac", "width": 1.5}, hoverinfo="skip"
            )
        ]
    )
    figure.add_trace(
        trace=go.Scatter(
            x=[depths[tip] for tip in tips],
            y=[positions[tip] for tip in tips],
            mode="markers",
            marker={"color": "#d47800", "size": 5},
            text=[labels[tip.name] for tip in tips],
            hovertemplate="%{text}<br>Divergence: %{x:.4g}<extra></extra>",
        )
    )
    figure.update_layout(
        showlegend=False,
        height=min(1200, max(400, 17 * len(tips))),
        xaxis_title=axis,
        yaxis={
            "tickvals": list(range(len(tips))),
            "ticktext": [labels[tip.name] for tip in tips],
            "autorange": "reversed",
        },
        margin={"l": 20, "r": 20, "t": 20, "b": 50},
    )
    return figure
