"""Verified coordinate models and sequence overlays for the result viewer.

Structural comparisons in the published tables have aggregate scores, not
residue mappings. Only sequence features with explicit, valid coordinates can
be projected onto an exactly matching full-length model.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

import pandas as pd

from protein_signatures.errors import InputValidationError

LOGGER = logging.getLogger(__name__)
MAX_MODEL_BYTES = 100 * 1024 * 1024
MAX_API_BYTES = 2 * 1024 * 1024
MAX_ALIGNMENT_CELLS = 1_500_000
MAX_ANNOTATION_BYTES = 10 * 1024 * 1024
ANNOTATION_COLUMNS = (
    "comparison_id",
    "protein_id",
    "region_type",
    "region_id",
    "start",
    "end",
    "q_value",
    "prevalence_difference",
    "evidence_source",
)
_ASSET_NAME = re.compile(r"assets/structures/([0-9a-f]{64})\.(pdb|pdb\.gz|cif|cif\.gz)")
_UNIPROT = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})")
_AMINO_ACIDS = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


@dataclass(frozen=True)
class Residue:
    """One C-alpha atom with its model chain, residue number and confidence."""

    chain: str
    position: int
    insertion: str
    amino_acid: str
    x: float
    y: float
    z: float
    confidence: float


@dataclass(frozen=True)
class ModelTrace:
    """The longest protein chain in one bounded PDB model."""

    residues: tuple[Residue, ...]
    chain: str

    def matches_sequence(self, *, sequence: str) -> bool:
        """Require a complete one-to-one mapping to canonical positions."""

        return len(self.residues) == len(sequence) and all(
            residue.position == index and not residue.insertion and residue.amino_acid == amino_acid
            for index, (residue, amino_acid) in enumerate(
                zip(self.residues, sequence, strict=True), start=1
            )
        )


@dataclass(frozen=True)
class AlignmentColumn:
    """One column in an exploratory global *sequence* alignment."""

    reference: str
    comparison: str
    reference_position: int | None
    comparison_position: int | None


def canonical_accession(*, value: object) -> str | None:
    """Accept only a complete canonical UniProt accession."""

    accession = str(value or "").strip().upper()
    return accession if _UNIPROT.fullmatch(accession) else None


def external_links(*, accession: object) -> dict[str, str]:
    """Return exact accession-specific pages for an independently verified ID."""

    value = canonical_accession(value=accession)
    if value is None:
        return {}
    return {
        "UniProt": f"https://www.uniprot.org/uniprotkb/{value}/entry",
        "AlphaFold DB": f"https://alphafold.ebi.ac.uk/entry/{value}",
        "Mol*": f"https://molstar.org/viewer/?{urlencode({'afdb': value})}",
        "InterPro": f"https://www.ebi.ac.uk/interpro/protein/UniProt/{value}/",
    }


def pair_links(*, reference: object, comparison: object) -> dict[str, str]:
    """Build external alignment URLs only for two exact UniProt accessions."""

    left = canonical_accession(value=reference)
    right = canonical_accession(value=comparison)
    if left is None or right is None or left == right:
        return {}
    request = {
        "query": {
            "context": {
                "mode": "pairwise",
                "method": {"name": "fatcat-rigid"},
                "structures": [
                    {"entry_id": f"AF_AF{accession}F1", "selection": {"asym_id": "A"}}
                    for accession in (left, right)
                ],
            }
        }
    }
    return {
        "EMERALD sequence alignment": "https://algbio.github.io/emerald-ui/?"
        + urlencode({"seqA": left, "seqB": right, "alpha": "0.75", "delta": "8"}),
        "RCSB pairwise structure alignment": "https://www.rcsb.org/alignment?"
        + urlencode({"request-body": json.dumps(request, separators=(",", ":"))}),
    }


def read_published_model(*, database: Path, relative_path: str, sha256: str) -> bytes:
    """Read an immutable, checksum-matched model within the verified result."""

    match = _ASSET_NAME.fullmatch(str(relative_path))
    if match is None or match.group(1) != str(sha256):
        raise InputValidationError("Model path or checksum does not match a published asset.")
    root = Path(database).resolve().parent
    path = (root / relative_path).resolve()
    if root not in path.parents or not path.is_file():
        raise InputValidationError("Published model is absent or outside the result bundle.")
    if not 0 < path.stat().st_size <= MAX_MODEL_BYTES:
        raise InputValidationError("Published model exceeds the app model size limit.")
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != sha256:
        raise InputValidationError("Published model checksum does not match the structure table.")
    if relative_path.endswith(".gz"):
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(payload)) as stream:
                payload = stream.read(MAX_MODEL_BYTES + 1)
        except (OSError, EOFError) as error:
            raise InputValidationError("Published model compression is invalid.") from error
    if not 0 < len(payload) <= MAX_MODEL_BYTES:
        raise InputValidationError("Uncompressed model exceeds the app model size limit.")
    return payload


def parse_pdb_trace(*, payload: bytes) -> ModelTrace:
    """Read C-alpha atoms from the first model of a PDB coordinate file."""

    if not 0 < len(payload) <= MAX_MODEL_BYTES:
        raise InputValidationError("PDB model is empty or exceeds the size limit.")
    chains: dict[str, dict[tuple[int, str], Residue]] = {}
    for line in payload.decode("ascii", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  ") or len(line) < 66:
            continue
        if line[12:16].strip() != "CA" or line[16] not in {" ", "A"}:
            continue
        amino_acid = _AMINO_ACIDS.get(line[17:20].strip())
        if amino_acid is None:
            continue
        try:
            position = int(line[22:26])
            x, y, z = (float(line[start : start + 8]) for start in (30, 38, 46))
            confidence = float(line[60:66])
        except ValueError:
            continue
        chain = line[21].strip() or "A"
        insertion = line[26].strip()
        key = position, insertion
        chains.setdefault(chain, {}).setdefault(
            key, Residue(chain, position, insertion, amino_acid, x, y, z, confidence)
        )
    if not chains:
        raise InputValidationError("No protein C-alpha atoms were found in the PDB model.")
    chain = max(sorted(chains), key=lambda value: len(chains[value]))
    residues = tuple(chains[chain].values())
    if len(residues) > 20_000:
        raise InputValidationError("Model has too many residues for the interactive viewer.")
    return ModelTrace(residues=residues, chain=chain)


def feature_intervals(
    *, rows: list[dict[str, object]], sequence_length: int
) -> list[dict[str, object]]:
    """Keep only recorded residue intervals that fit the source sequence."""

    intervals = []
    for row in rows:
        try:
            start, end = int(row["start"]), int(row["end"])
        except (KeyError, TypeError, ValueError):
            continue
        if 1 <= start <= end <= sequence_length:
            intervals.append({**row, "start": start, "end": end})
    return intervals


def parse_annotation_tsv(*, payload: bytes) -> pd.DataFrame:
    """Read an optional bounded, provenance-bearing region annotation table.

    Rows are scoped again to the chosen comparison and protein in the app;
    their coordinates are validated against its canonical sequence there.
    """

    if not 0 < len(payload) <= MAX_ANNOTATION_BYTES:
        raise InputValidationError("Annotation TSV is empty or larger than 10 MiB.")
    try:
        frame = pd.read_csv(io.BytesIO(payload), sep="\t", dtype=str, keep_default_na=False)
    except (UnicodeError, pd.errors.ParserError) as error:
        raise InputValidationError("Annotation TSV could not be parsed.") from error
    if tuple(frame.columns) != ANNOTATION_COLUMNS or len(frame) > 10_000:
        raise InputValidationError(
            "Annotation TSV needs the nine documented columns in order and at most 10,000 rows."
        )
    for index, row in frame.iterrows():
        if any(not str(row[name]).strip() for name in ANNOTATION_COLUMNS):
            raise InputValidationError(f"Annotation row {index + 2} contains a blank field.")
        if any(
            len(str(row[name])) > 255 or re.search(r"[<>\x00-\x1f]", str(row[name]))
            for name in (
                "comparison_id",
                "protein_id",
                "region_type",
                "region_id",
                "evidence_source",
            )
        ):
            raise InputValidationError(f"Annotation row {index + 2} has unsafe text.")
        try:
            start, end = int(row["start"]), int(row["end"])
            q_value = float(row["q_value"])
            difference = float(row["prevalence_difference"])
        except ValueError as error:
            raise InputValidationError(
                f"Annotation row {index + 2} has invalid numbers."
            ) from error
        if start < 1 or end < start or not math.isfinite(q_value) or not 0 <= q_value <= 1:
            raise InputValidationError(f"Annotation row {index + 2} has invalid bounds or q-value.")
        if not math.isfinite(difference):
            raise InputValidationError(f"Annotation row {index + 2} has invalid enrichment.")
    return frame


def enrichment_track(
    *,
    rows: list[dict[str, object]],
    sequence_length: int,
    q_threshold: float = 0.05,
    q_floor: float = 1e-8,
) -> tuple[list[float], list[str]]:
    """Project positive, significant interval enrichment onto sequence positions.

    White (0) means no significant mapped enrichment, including missing data;
    it does not assert that the position was statistically tested. Overlaps use
    the smallest q-value. Blue begins at 0.1, red is 1.0 at the q-value floor.
    """

    if sequence_length < 1 or not 0 < q_floor < q_threshold <= 1:
        raise InputValidationError("Invalid sequence length or colour-scale q-value bounds.")
    q_values: list[float | None] = [None] * sequence_length
    descriptions = ["No significant mapped enrichment"] * sequence_length
    for row in significant_intervals(
        rows=rows, sequence_length=sequence_length, q_threshold=q_threshold
    ):
        q_value = float(row["q_value"])
        difference = float(row["prevalence_difference"])
        label = str(row.get("feature_name") or row.get("region_id") or "Region")
        for index in range(row["start"] - 1, row["end"]):
            if q_values[index] is None or q_value < q_values[index]:
                q_values[index] = q_value
                descriptions[index] = f"{label} · q={q_value:.3g} · Δ={difference:.3g}"
    minimum = -math.log10(q_threshold)
    maximum = -math.log10(q_floor)
    scores = [
        0.0
        if value is None
        else 0.1
        + 0.9
        * min(1.0, max(0.0, (-math.log10(max(value, q_floor)) - minimum) / (maximum - minimum)))
        for value in q_values
    ]
    return scores, descriptions


def significant_intervals(
    *, rows: list[dict[str, object]], sequence_length: int, q_threshold: float = 0.05
) -> list[dict[str, object]]:
    """Select only mapped, positive and statistically significant intervals."""

    selected = []
    for row in feature_intervals(rows=rows, sequence_length=sequence_length):
        try:
            q_value = float(row["q_value"])
            difference = float(row["prevalence_difference"])
        except (KeyError, TypeError, ValueError):
            continue
        if (
            math.isfinite(q_value)
            and math.isfinite(difference)
            and 0 <= q_value <= q_threshold
            and difference > 0
        ):
            selected.append(row)
    return selected


def align_sequences(*, reference: str, comparison: str) -> tuple[AlignmentColumn, ...]:
    """Perform a bounded deterministic global sequence alignment.

    Scores are +2 match, -1 substitution, -2 gap. This exploratory alignment
    is distinct from the published Foldseek structural comparison.
    """

    if (
        not reference
        or not comparison
        or not all(re.fullmatch(r"[A-Z]+", value) for value in (reference, comparison))
    ):
        raise InputValidationError("Sequence alignment requires two amino-acid strings.")
    n, m = len(reference), len(comparison)
    if n * m > MAX_ALIGNMENT_CELLS:
        raise InputValidationError("Pair is too long for the bounded in-app sequence alignment.")
    # Two score rows and a byte per traceback cell keep memory bounded.
    directions = [bytearray(m + 1) for _ in range(n + 1)]
    previous = [-2 * j for j in range(m + 1)]
    for i in range(1, n + 1):
        current = [-2 * i] + [0] * m
        directions[i][0] = 1
        for j in range(1, m + 1):
            diagonal = previous[j - 1] + (2 if reference[i - 1] == comparison[j - 1] else -1)
            above = previous[j] - 2
            left = current[j - 1] - 2
            best = max(diagonal, above, left)
            current[j] = best
            directions[i][j] = 0 if best == diagonal else 1 if best == above else 2
        previous = current
    i, j = n, m
    columns = []
    while i or j:
        direction = directions[i][j] if i and j else 1 if i else 2
        if direction == 0:
            columns.append(AlignmentColumn(reference[i - 1], comparison[j - 1], i, j))
            i -= 1
            j -= 1
        elif direction == 1:
            columns.append(AlignmentColumn(reference[i - 1], "-", i, None))
            i -= 1
        else:
            columns.append(AlignmentColumn("-", comparison[j - 1], None, j))
            j -= 1
    return tuple(reversed(columns))


def alignment_fasta(
    *, reference_id: str, comparison_id: str, columns: tuple[AlignmentColumn, ...]
) -> bytes:
    """Serialise the exact displayed aligned sequences as FASTA."""

    identifiers = (reference_id, comparison_id)
    if any(
        not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:+/-]{0,254}", value) for value in identifiers
    ):
        raise InputValidationError("Alignment FASTA identifiers contain unsafe characters.")
    sequences = (
        "".join(column.reference for column in columns),
        "".join(column.comparison for column in columns),
    )
    return (
        "\n".join(
            "\n".join(
                (f">{identifier}", *(sequence[i : i + 80] for i in range(0, len(sequence), 80)))
            )
            for identifier, sequence in zip(identifiers, sequences, strict=True)
        )
        + "\n"
    ).encode("ascii")


class _RestrictedRedirect(HTTPRedirectHandler):
    """Never send an AlphaFold request to another host or insecure scheme."""

    def redirect_request(self, request, fp, code, message, headers, new_url):
        """Check the destination before urllib follows a redirect."""

        if not _approved_alphafold_url(url=new_url):
            raise InputValidationError("AlphaFold redirected outside its approved host.")
        return super().redirect_request(request, fp, code, message, headers, new_url)


def _approved_alphafold_url(*, url: str) -> bool:
    """Require an HTTPS resource on the exact AlphaFold DB host."""

    parsed = urlsplit(url)
    return parsed.scheme == "https" and parsed.netloc == "alphafold.ebi.ac.uk"


def fetch_alphafold_model(*, accession: str, sequence: str) -> bytes:
    """Fetch a current AFDB PDB on user action, verifying its sequence exactly.

    The response is returned in memory; it does not alter the published result.
    """

    canonical = canonical_accession(value=accession)
    if canonical is None:
        raise InputValidationError("A canonical UniProt accession is required.")
    opener = build_opener(_RestrictedRedirect())
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{canonical}"
    try:
        with opener.open(
            Request(api, headers={"Accept": "application/json"}), timeout=25
        ) as response:
            metadata_bytes = response.read(MAX_API_BYTES + 1)
        if len(metadata_bytes) > MAX_API_BYTES:
            raise InputValidationError("AlphaFold metadata exceeds the size limit.")
        records = json.loads(metadata_bytes)
        if not isinstance(records, list):
            records = [records]
        matching = [
            record
            for record in records
            if isinstance(record, dict)
            and str(record.get("sequence", "")).upper() == sequence.upper()
        ]
        if len(matching) != 1:
            raise InputValidationError("No single AlphaFold model matches this exact sequence.")
        url = matching[0].get("pdbUrl")
        if not isinstance(url, str) or not _approved_alphafold_url(url=url):
            raise InputValidationError("AlphaFold did not provide an approved PDB URL.")
        with opener.open(Request(url), timeout=45) as response:
            payload = response.read(MAX_MODEL_BYTES + 1)
        if not 0 < len(payload) <= MAX_MODEL_BYTES:
            raise InputValidationError("Downloaded AlphaFold model exceeds the size limit.")
        trace = parse_pdb_trace(payload=payload)
        if not trace.matches_sequence(sequence=sequence):
            raise InputValidationError("AlphaFold PDB residues do not match this exact sequence.")
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        LOGGER.warning("AlphaFold retrieval failed for %s: %s", canonical, error)
        raise InputValidationError(
            "AlphaFold model retrieval failed; try the external record."
        ) from error
    return payload
