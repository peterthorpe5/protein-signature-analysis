"""Scientific mapping and safety checks for the model and alignment explorer."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pandas as pd
import pytest

import protein_signature_app.app as app
from protein_signature_app.structure_viewer import (
    align_sequences,
    alignment_fasta,
    canonical_accession,
    enrichment_track,
    external_links,
    fetch_alphafold_model,
    pair_links,
    parse_annotation_tsv,
    parse_mmcif_trace,
    parse_pdb_trace,
    project_kmer_intervals,
    read_published_model,
    significant_intervals,
)
from protein_signatures.errors import InputValidationError


def _pdb() -> bytes:
    """Return a minimal complete PDB with ALA, CYS and GLY."""

    return (
        "HEADER    TEST\n"
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 90.00           C\n"
        "ATOM      2  CA  CYS A   2       2.000   2.000   3.000  1.00 80.00           C\n"
        "ATOM      3  CA  GLY A   3       3.000   2.000   3.000  1.00 70.00           C\n"
        "END\n"
    ).encode("ascii")


def _mmcif() -> bytes:
    """Return a three-residue mmCIF with a skipped alternative conformation."""

    return (
        "data_test\n#\nloop_\n"
        "_atom_site.group_PDB\n_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n"
        "ATOM CA . ALA A 1 ? 1.0 2.0 3.0 90.0 1\n"
        "ATOM CA B CYS A 2 ? 2.0 2.0 3.0 50.0 1\n"
        "ATOM CA . CYS A 2 ? 2.0 2.0 3.0 80.0 1\n"
        "ATOM CA . GLY A 3 ? 3.0 2.0 3.0 70.0 1\n"
        "ATOM CA . ALA A 1 ? 9.0 2.0 3.0 10.0 2\n#\n"
    ).encode("ascii")


def test_mmcif_trace_maps_exact_residues_and_rejects_malformed_loops() -> None:
    """Packaged mmCIF models support exact full-length colouring like PDB models."""

    trace = parse_mmcif_trace(payload=_mmcif())
    assert trace.matches_sequence(sequence="ACG")
    assert trace.residues[1].confidence == 80.0
    assert trace.residues[0].x == 1.0
    assert not trace.matches_sequence(sequence="AGG")
    with pytest.raises(InputValidationError, match="required coordinate columns"):
        parse_mmcif_trace(payload=b"data_bad\nloop_\n_atom_site.group_PDB\nATOM\n#\n")
    with pytest.raises(InputValidationError, match="incomplete columns"):
        parse_mmcif_trace(payload=_mmcif().replace(b"80.0 1\n", b"80.0\n"))


def test_published_kmer_membership_projects_exact_overlapping_positions() -> None:
    """K-mer enrichment can be displayed where the sequence actually contains it."""

    rows = [
        {
            "feature_type": "AMINO_ACID_KMER",
            "feature_id": "k3:AAA",
            "start": None,
            "q_value": 0.001,
            "prevalence_difference": 0.25,
        },
        {
            "feature_type": "AMINO_ACID_KMER",
            "feature_id": "k2:AA",
            "start": pd.NA,
            "q_value": 0.002,
            "prevalence_difference": 0.2,
        },
        {"feature_type": "AMINO_ACID_KMER", "feature_id": "k4:AAA", "start": None},
        {"feature_type": "FOLD", "feature_id": "fold", "start": None},
    ]
    projected = project_kmer_intervals(rows=rows, sequence="AAAA")
    assert [(row["start"], row["end"]) for row in projected if pd.notna(row["start"])] == [
        (1, 3),
        (2, 4),
        (1, 2),
        (2, 3),
        (3, 4),
    ]
    assert projected[-1]["feature_type"] == "FOLD"


def test_model_is_checksum_bound_and_exactly_numbered(tmp_path: Path) -> None:
    """Only bundled models with matching amino acids and numbering get colours."""

    database = tmp_path / "protein_signatures.duckdb"
    database.write_bytes(b"database")
    payload = _pdb()
    digest = hashlib.sha256(payload).hexdigest()
    relative = f"assets/structures/{digest}.pdb.gz"
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    # The published checksum is the checksum of the packaged compressed bytes.
    zipped = gzip.compress(payload, mtime=0)
    zipped_digest = hashlib.sha256(zipped).hexdigest()
    path = tmp_path / f"assets/structures/{zipped_digest}.pdb.gz"
    path.write_bytes(zipped)
    assert (
        read_published_model(
            database=database,
            relative_path=str(path.relative_to(tmp_path)),
            sha256=zipped_digest,
        )
        == payload
    )
    trace = parse_pdb_trace(payload=payload)
    assert trace.matches_sequence(sequence="ACG")
    assert not trace.matches_sequence(sequence="ACA")
    assert trace.residues[1].confidence == 80.0
    with pytest.raises(InputValidationError, match="path or checksum"):
        read_published_model(database=database, relative_path="../../outside.pdb", sha256=digest)
    path.write_bytes(b"tampered")
    with pytest.raises(InputValidationError, match="checksum"):
        read_published_model(
            database=database,
            relative_path=str(path.relative_to(tmp_path)),
            sha256=zipped_digest,
        )


def test_enrichment_track_uses_strongest_q_and_ignores_nonpositive_hits() -> None:
    """White, blue and red encode only mapped positive significant intervals."""

    rows = [
        {
            "start": 2,
            "end": 4,
            "q_value": 0.05,
            "prevalence_difference": 0.3,
            "feature_name": "weak domain",
        },
        {
            "start": 3,
            "end": 3,
            "q_value": 1e-9,
            "prevalence_difference": 0.4,
            "feature_name": "strong motif",
        },
        {
            "start": 1,
            "end": 1,
            "q_value": 0.001,
            "prevalence_difference": -0.2,
            "feature_name": "depleted",
        },
        {
            "start": 5,
            "end": 5,
            "q_value": 0.7,
            "prevalence_difference": 0.4,
            "feature_name": "not significant",
        },
    ]
    values, descriptions = enrichment_track(rows=rows, sequence_length=5)
    assert values == [0.0, 0.1, 1.0, 0.1, 0.0]
    assert "strong motif" in descriptions[2]
    assert "No significant" in descriptions[0]
    assert len(significant_intervals(rows=rows, sequence_length=5)) == 2


def test_held_out_tiers_use_their_own_q_values_and_direction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Validated overlays must be selected from independently validated classes."""

    frame = pd.DataFrame(
        [
            {
                "start": 2,
                "end": 3,
                "discovery_q_value": 0.001,
                "discovery_prevalence_difference": 0.2,
                "validation_q_value": 0.03,
                "validation_study_q_value": 0.08,
                "validation_prevalence_difference": 0.1,
                "evidence_class": "DECISION_CANDIDATE__VALIDATED_WITHIN_COMPARISON",
            },
            {
                "start": 5,
                "end": 6,
                "discovery_q_value": 0.002,
                "discovery_prevalence_difference": 0.3,
                "validation_q_value": 0.001,
                "validation_study_q_value": 0.04,
                "validation_prevalence_difference": 0.15,
                "evidence_class": "DECISION_CANDIDATE__VALIDATED_STUDY_WIDE",
            },
        ]
    )
    monkeypatch.setattr(app, "query_dataframe", lambda **_: frame)
    within = app._enriched_intervals(
        database=Path("unused"),
        protein_id="p1",
        comparison_id="cmp",
        evidence_tier="Validated within comparison",
    )
    assert {row["q_value"] for row in within} == {0.03, 0.001}
    wide = app._enriched_intervals(
        database=Path("unused"),
        protein_id="p1",
        comparison_id="cmp",
        evidence_tier="Validated study-wide",
    )
    assert [row["q_value"] for row in wide] == [0.04]


def test_explorer_projects_published_kmers_and_suggests_target_models(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The app links an enriched motif and packaged models to the selected target class."""

    seen_parameters: list[tuple[str, ...]] = []

    def query(*, sql: str, parameters: tuple[str, ...] = (), **_kwargs: object) -> pd.DataFrame:
        """Return relevant published membership rows for the two queries."""

        if "FROM structures s JOIN label_memberships" in sql:
            assert "IN (?,?)" in sql
            assert "coordinate_path <> ''" in sql
            return pd.DataFrame({"protein_id": ["p1", "p2"]})
        seen_parameters.append(tuple(parameters))
        return pd.DataFrame(
            [
                {
                    "feature_type": "AMINO_ACID_KMER",
                    "feature_id": "k3:AAA",
                    "feature_name": "k3:AAA",
                    "start": None,
                    "end": None,
                    "discovery_q_value": 0.001,
                    "discovery_prevalence_difference": 0.2,
                    "evidence_class": "DISCOVERY_ONLY",
                }
            ]
        )

    monkeypatch.setattr(app, "query_dataframe", query)
    assert app._target_model_proteins(
        database=Path("unused"), target_label_ids="target:one|target:two"
    ) == ("p1", "p2")
    intervals = app._enriched_intervals(
        database=Path("unused"), protein_id="p1", comparison_id="cmp", sequence="AAAA"
    )
    assert [(row["start"], row["end"]) for row in intervals] == [(1, 3), (2, 4)]
    assert seen_parameters == [("p1", "cmp", 50)]


def test_3d_trace_uses_the_same_enrichment_scale_as_the_sequence_track(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The rotatable 3D view must colour actual model residues from mapped scores."""

    rendered = []
    monkeypatch.setattr(app, "_render_plotly_figure", lambda **kwargs: rendered.append(kwargs))
    trace = parse_pdb_trace(payload=_pdb())
    app._render_trace_figure(
        trace=trace,
        scores=[0.0, 0.1, 1.0],
        descriptions=["No significant mapped enrichment", "region A", "region B"],
        download_name="test_model",
    )
    figure = rendered[0]["figure"]
    assert list(figure.data[1].marker.color) == [0.0, 0.1, 1.0]
    assert figure.layout.scene.aspectmode == "data"
    assert "region B" in figure.data[1].text[2]
    rendered.clear()
    app._render_trace_figure(
        trace=trace, scores=None, descriptions=None, download_name="unmapped_model"
    )
    assert rendered[0]["figure"].data[1].marker.showscale is False


def test_annotation_upload_is_validated_and_pair_alignment_preserves_positions() -> None:
    """An optional future pocket can be highlighted without inventing structure pairs."""

    header = (
        "comparison_id\tprotein_id\tregion_type\tregion_id\tstart\tend\t"
        "q_value\tprevalence_difference\tevidence_source\n"
    )
    frame = parse_annotation_tsv(
        payload=(header + "cmp\tp1\tPOCKET\tsite1\t2\t3\t0.001\t0.2\trun42\n").encode()
    )
    assert frame.iloc[0]["region_type"] == "POCKET"
    with pytest.raises(InputValidationError, match="invalid bounds"):
        parse_annotation_tsv(
            payload=(header + "cmp\tp1\tPOCKET\tsite1\t0\t3\t0.001\t0.2\trun42\n").encode()
        )
    columns = align_sequences(reference="ACG", comparison="AG")
    assert "".join(item.reference for item in columns) == "ACG"
    assert "".join(item.comparison for item in columns) == "A-G"
    assert columns[1].comparison_position is None
    assert alignment_fasta(reference_id="p1", comparison_id="p2", columns=columns) == (
        b">p1\nACG\n>p2\nA-G\n"
    )


def test_accession_links_fail_closed_and_remote_requires_exact_sequence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """External sites receive only validated IDs; remote models are sequence checked."""

    assert canonical_accession(value="P12345") == "P12345"
    assert canonical_accession(value="P12345-2") is None
    assert external_links(accession="local-p1") == {}
    assert pair_links(reference="P12345", comparison="Q9XYZ1")[
        "EMERALD sequence alignment"
    ].startswith("https://algbio.github.io/emerald-ui/?seqA=P12345&seqB=Q9XYZ1")

    class FakeResponse:
        """Small bounded HTTP response test double."""

        def __init__(self, payload: bytes) -> None:
            self.payload = payload

        def __enter__(self) -> FakeResponse:
            return self

        def __exit__(self, *_args: object) -> None:
            return None

        def read(self, _limit: int) -> bytes:
            return self.payload

    class FakeOpener:
        """Return an API record followed by a model."""

        def open(self, request: object, timeout: int) -> FakeResponse:
            if "api/prediction" in request.full_url:
                return FakeResponse(
                    b'[{"sequence":"ACG","pdbUrl":"https://alphafold.ebi.ac.uk/files/model.pdb"}]'
                )
            return FakeResponse(_pdb())

    monkeypatch.setattr(
        "protein_signature_app.structure_viewer.build_opener", lambda *_: FakeOpener()
    )
    assert fetch_alphafold_model(accession="P12345", sequence="ACG") == _pdb()
    with pytest.raises(InputValidationError, match="exact sequence"):
        fetch_alphafold_model(accession="P12345", sequence="AAA")
