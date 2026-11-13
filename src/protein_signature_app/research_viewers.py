"""Offline JBrowse and Nightingale documents using the E3 viewer assets."""

from __future__ import annotations

import hashlib
import io
import json
import math
import re
from functools import lru_cache
from importlib.resources import files

from Bio import Phylo, SeqIO
from Bio.Phylo.NewickIO import NewickError

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.catalogue import identifiers, normalise_sequence

VIEWER_BYTES = 20 * 1024 * 1024
# Synthetic demonstration from the MIT-licensed E3 v0.23.1 viewer.
EXAMPLE_ALIGNMENT = (
    ">Example_A\nMAAKLVVAGLTQLKEAGLDVAVIDRHGLDVEQLKAALEAGADVVVVTADKETQLAS\n"
    ">Example_B\nMAAKLIVAGLTQLKEAGLDVAVIDRHGLDVEQLKAALEAGADVVIVTADKETQLAS\n"
    ">Example_C\nMSAKLVVAGLTQLKDAGLDVAVIDRHGLDVEQLKSALEAGADIVVVTADKETQLAT\n"
    ">Example_D\nMSAKLIVAGLTQLKDAGLDVAVIDRHGL--EQLKSALEAGADIVIVTADKETQLAT\n"
    ">Example_E\nMNARLIVAGVTQLREAGMDVAIVDRHGLDVEQLKSALESGADVIVVTADRETQLAT\n"
    ">Example_F\nMNARLIVAGVTQLREAGMDVAIVDRHGLDV-QLKSALESGADVIVITADRETQLAT\n"
)
EXAMPLE_TREE = (
    "((Example_A:0.03,Example_B:0.03):0.04,"
    "((Example_C:0.04,Example_D:0.05):0.03,(Example_E:0.03,Example_F:0.04):0.06):0.03);"
)


def fasta_records(*, text: str, aligned: bool = False) -> dict[str, str]:
    """Validate a bounded exact FASTA or a rectangular alignment.

    Args:
        text: UTF-8 decoded FASTA.
        aligned: Permit alignment gaps and require a rectangular matrix.

    Returns:
        Unique sequences keyed by their exact FASTA IDs.
    """
    if len(text.encode()) > VIEWER_BYTES or not text.lstrip().startswith(">"):
        raise InputValidationError("Provide a FASTA file of at most 20 MiB.")
    proteins, residues = {}, 0
    for record in SeqIO.parse(io.StringIO(text.lstrip()), "fasta"):
        name = identifiers(values=[record.id])[0]
        if name in proteins:
            raise InputValidationError(f"Duplicate FASTA identifier: {name}.")
        sequence = str(record.seq).upper()
        if aligned:
            ungapped = sequence.replace("-", "")
            if normalise_sequence(value=ungapped) != ungapped:
                raise InputValidationError("Remove terminal stop symbols from aligned FASTA.")
        else:
            sequence = normalise_sequence(value=sequence)
        proteins[name] = sequence
        residues += len(sequence)
        if len(proteins) > 2000 or residues > 5_000_000:
            raise InputValidationError("Viewer limit: 2000 proteins and five million cells.")
    if not proteins or (aligned and len({len(value) for value in proteins.values()}) != 1):
        raise InputValidationError("The alignment FASTA must be non-empty and rectangular.")
    return proteins


def alignment_payload(
    *,
    fasta: str,
    newick: str = "",
    tip_labels: dict[str, str] | None = None,
) -> dict:
    """Prepare one MSA with an optional exactly matched and validated tree.

    Args:
        fasta: Aligned protein FASTA.
        newick: Optional complete Newick tree.
        tip_labels: Explicit internal tree-tip to exact FASTA identifier mapping.

    Returns:
        Validated JBrowse input. The supplied root and lengths are retained.
    """
    proteins = fasta_records(text=fasta, aligned=True)
    tree_text = ""
    if newick:
        if len(newick.encode()) > 1_000_000 or not newick.strip().endswith(";"):
            raise InputValidationError("Provide one complete Newick tree of at most 1 MB.")
        try:
            tree = Phylo.read(io.StringIO(newick), "newick")
            leaves = tree.get_terminals()
            names = [tip.name for tip in leaves]
            if len(set(names)) != len(leaves) or any(not name for name in names):
                raise InputValidationError("Tree tips must have unique names.")
            if tip_labels is not None and set(names) != set(tip_labels):
                raise InputValidationError("Tree tips do not match their explicit ID mapping.")
            for tip in leaves:
                if tip_labels is not None:
                    tip.name = tip_labels[tip.name]
            mapped = [tip.name for tip in leaves]
            if len(set(mapped)) != len(mapped) or set(mapped) != set(proteins):
                raise InputValidationError("Tree tips must exactly match the alignment IDs.")
            if any(
                clade.branch_length is not None
                and (not math.isfinite(clade.branch_length) or clade.branch_length < 0)
                for clade in tree.find_clades()
            ):
                raise InputValidationError("Tree branch lengths must be finite and non-negative.")
            handle = io.StringIO()
            Phylo.write(tree, handle, "newick", format_branch_length="%1.12g")
            tree_text = handle.getvalue()
        except (ValueError, NewickError, RecursionError) as error:
            raise InputValidationError("Could not read one bounded Newick tree.") from error
    row_height = 26 if len(proteins) <= 16 else 19
    return {
        "msa": "".join(f">{name}\n{sequence}\n" for name, sequence in proteins.items()),
        "tree": tree_text,
        "row_height": row_height,
        "height": min(640, max(250, 118 + row_height * len(proteins))),
        "alignment_length": len(next(iter(proteins.values()))),
    }


def feature_payload(
    *,
    identifier: str,
    sequence: str,
    rows: list[dict],
    evidence_label: str = "predicted annotation",
) -> dict:
    """Map only exact-sequence positional predictions to Nightingale tracks.

    Args:
        identifier: Exact selected protein ID.
        sequence: Its unaligned sequence.
        rows: Normalised evidence rows carrying SHA-256 and inclusive positions.
        evidence_label: Published or predicted scope displayed beside the protein ID.

    Returns:
        Browser-safe data; non-positional results remain in the results table.
    """
    sequence = normalise_sequence(value=sequence)
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    features = []
    for row in rows:
        if row.get("sequence_id") != identifier:
            continue
        if row.get("sequence_sha256") != digest:
            raise InputValidationError("Feature evidence belongs to a different sequence version.")
        if row.get("start") in ("", None) and row.get("end") in ("", None):
            continue
        if not all(re.fullmatch(r"[1-9][0-9]*", str(row.get(key))) for key in ("start", "end")):
            raise InputValidationError("Feature positions must be positive integer coordinates.")
        start, end = int(row["start"]), int(row["end"])
        if not 1 <= start <= end <= len(sequence):
            raise InputValidationError("Feature position is outside the selected protein.")
        features.append(
            {
                "accession": str(row.get("accession", "")),
                "description": str(row.get("description", "")),
                "color": "#346a8a",
                "locations": [{"fragments": [{"start": start, "end": end}]}],
            }
        )
        if len(features) > 5000:
            raise InputValidationError("Viewer limit: 5000 positional predictions per protein.")
    return {
        "sequence": sequence,
        "features": features,
        "coordinate_scope": f"{identifier} · exact sequence · {evidence_label}",
    }


@lru_cache(maxsize=2)
def viewer_script(*, viewer: str) -> str:
    """Load and verify one maintained viewer bundle from the installed package.

    Args:
        viewer: jbrowse or nightingale.

    Returns:
        Trusted script with HTML end tags escaped for inline delivery.
    """
    if viewer not in {"jbrowse", "nightingale"}:
        raise InputValidationError("Unknown research viewer.")
    root = files("protein_signature_app").joinpath("resources", "research_viewers")
    name = f"{viewer}.js"
    try:
        digest = json.loads(root.joinpath("checksums.json").read_text(encoding="utf-8"))[name]
        content = root.joinpath(name).read_bytes()
    except (KeyError, ValueError, OSError) as error:
        raise InputValidationError("Reinstall the package to restore its viewer assets.") from error
    if hashlib.sha256(content).hexdigest() != digest:
        raise InputValidationError("The installed viewer asset checksum differs.")
    return re.sub(r"</script", r"<\\/script", content.decode(), flags=re.IGNORECASE)


def viewer_document(*, viewer: str, payload: dict) -> str:
    """Package a bounded viewer as a standalone document without CDN scripts.

    Args:
        viewer: Bundled viewer name.
        payload: Validated data for the viewer.

    Returns:
        Self-contained HTML with untrusted data encoded as inert JSON.
    """
    data = json.dumps(payload, ensure_ascii=True, allow_nan=False)
    if len(data) > VIEWER_BYTES:
        raise InputValidationError("Viewer payload exceeds its rendering bound.")
    data = data.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    script = viewer_script(viewer=viewer)
    bootstrap = (
        """
const proteinData = JSON.parse(document.getElementById('protein-data').textContent);
const container = document.getElementById('viewer');
const status = document.getElementById('status');
try {
    const width = container.clientWidth || 900;
    const treeAreaWidth = Math.min(235, Math.max(110, Math.floor(width * 0.24)));
    const columns = Math.min(proteinData.alignment_length || 60, 60);
    const colWidth = Math.max(14, Math.min(24, Math.floor((width - treeAreaWidth - 28) / columns)));
    window.proteinMSA = window.ReactMSAView.mount(container, {
        msa: proteinData.msa, tree: proteinData.tree || undefined,
        height: proteinData.height || 640, rowHeight: proteinData.row_height || 22,
        colWidth, treeAreaWidth, colorScheme: 'clustal'
    });
    status.textContent = proteinData.tree ? 'JBrowse MSA · alignment with tree' : 'JBrowse MSA';
} catch (error) {
    status.textContent = 'Viewer error: ' + error.message;
}
"""
        if viewer == "jbrowse"
        else ""
    )
    return (
        """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy"
content="default-src 'none';script-src 'unsafe-inline';style-src 'unsafe-inline';
img-src data: blob:;font-src data:;connect-src 'none';worker-src blob:">
<title>Protein research viewer</title><style>
html,body{margin:0;background:#fff;color:#1b2a3b;font:14px Arial,sans-serif}
#status{padding:9px 12px;background:#eff5fa;font-weight:600}
#viewer{width:100%;position:relative}
.protein-feature-viewer{padding:0 12px 10px}
nightingale-navigation,nightingale-sequence,nightingale-track{display:block;width:100%}
.protein-feature-row{border-top:1px solid #d9e1eb;padding:5px 0}
.protein-feature-label{font-weight:600;overflow-wrap:anywhere}
.protein-feature-label span{font-size:12px;font-weight:400;margin-left:10px;color:#425671}
</style></head><body><div id="status" role="status">Loading viewer…</div><div id="viewer"></div>
<script id="protein-data" type="application/json">"""
        + data
        + "</script><script>"
        + script
        + "</script><script>"
        + bootstrap
        + "</script></body></html>"
    )
