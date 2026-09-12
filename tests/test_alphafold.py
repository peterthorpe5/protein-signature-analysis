"""Tests for AlphaFold Database mapping, acquisition and confidence parsing."""

from __future__ import annotations

import io
import json
import urllib.error
from pathlib import Path

import pytest

import protein_signatures.alphafold as af
from protein_signatures.checksums import sha256_text
from protein_signatures.errors import InputValidationError
from protein_signatures.models import (
    AlphaFoldRequest,
    AlphaFoldSettings,
    SequenceRecord,
    StructureAnalysisEligibility,
)


def test_read_requests_and_disabled_acquisition(tmp_path: Path) -> None:
    """Mappings should be unique and disabled acquisition should remain explicit."""

    table = tmp_path / "accessions.tsv"
    table.write_text("protein_id\tuniprot_accession\np2\tq22222\np1\tP11111\n", encoding="utf-8")
    requests = af.read_alphafold_requests(path=table, protein_ids=frozenset({"p1", "p2"}))
    assert [(item.protein_id, item.uniprot_accession) for item in requests] == [
        ("p1", "P11111"),
        ("p2", "Q22222"),
    ]
    outcomes, structures = af.acquire_alphafold_models(
        requests=requests,
        sequences=_sequences(),
        settings=_settings(cache=tmp_path / "cache", enabled=False),
    )
    assert structures == ()
    assert {item.acquisition_status for item in outcomes} == {"NOT_SELECTED"}
    table.write_text("protein_id\tuniprot_accession\nunknown\tP1\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="unknown protein"):
        af.read_alphafold_requests(path=table, protein_ids=frozenset({"p1"}))
    table.write_text("protein_id\tuniprot_accession\np1\tP1\np1\tP2\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="multiple"):
        af.read_alphafold_requests(path=table, protein_ids=frozenset({"p1"}))
    with pytest.raises(InputValidationError, match="enabled"):
        af.acquire_alphafold_models(
            requests=(),
            sequences=_sequences(),
            settings=_settings(cache=tmp_path / "cache", enabled=True),
        )


@pytest.mark.parametrize(
    ("remote_sequence", "plddt", "threshold", "expected", "eligibility"),
    [
        ("ACDE", 90.0, 50.0, "ACQUIRED", StructureAnalysisEligibility.ELIGIBLE),
        (
            "ACDE",
            40.0,
            50.0,
            "ACQUIRED_LOW_CONFIDENCE",
            StructureAnalysisEligibility.INELIGIBLE_LOW_CONFIDENCE,
        ),
        (
            "",
            90.0,
            50.0,
            "ACQUIRED_SEQUENCE_UNVERIFIED",
            StructureAnalysisEligibility.INELIGIBLE_SEQUENCE_UNVERIFIED,
        ),
        (
            "ACDE",
            None,
            50.0,
            "ACQUIRED_CONFIDENCE_UNAVAILABLE",
            StructureAnalysisEligibility.INELIGIBLE_CONFIDENCE_UNAVAILABLE,
        ),
    ],
)
def test_acquisition_success_states(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    remote_sequence: str,
    plddt: float | None,
    threshold: float,
    expected: str,
    eligibility: StructureAnalysisEligibility,
) -> None:
    """Acquisition should separate confidence and sequence-verification states."""

    metadata = _metadata(sequence=remote_sequence)
    pdb = _pdb_bytes(plddt=plddt)
    responses = iter((json.dumps([metadata]).encode(), pdb))
    monkeypatch.setattr(af, "_request_bytes", lambda **_kwargs: next(responses))
    outcomes, structures = af.acquire_alphafold_models(
        requests=(AlphaFoldRequest("p1", "P11111"),),
        sequences=_sequences(),
        settings=_settings(cache=tmp_path / "cache", enabled=True, threshold=threshold),
    )
    assert outcomes[0].acquisition_status == expected
    assert outcomes[0].coordinate_path is not None
    assert outcomes[0].coordinate_sha256
    assert len(structures) == 1
    assert structures[0].structure_source == "AlphaFoldDB"
    assert structures[0].fold_evidence_status == "NOT_ASSESSED"
    assert structures[0].analysis_eligibility_status == eligibility
    assert structures[0].is_coordinate_analysis_eligible is (
        eligibility == StructureAnalysisEligibility.ELIGIBLE
    )
    assert structures[0].to_record()["analysis_eligibility_status"] == eligibility.value
    cached = structures[0].coordinate_path
    assert cached is not None and cached.read_bytes() == pdb
    monkeypatch.setattr(
        af,
        "_request_bytes",
        lambda **kwargs: (
            json.dumps([metadata]).encode()
            if kwargs["url"].endswith("P11111")
            else pytest.fail("cached coordinate should not be downloaded again")
        ),
    )
    reused, reused_structures = af.acquire_alphafold_models(
        requests=(AlphaFoldRequest("p1", "P11111"),),
        sequences=_sequences(),
        settings=_settings(cache=tmp_path / "cache", enabled=True, threshold=threshold),
    )
    assert reused[0].coordinate_sha256 == outcomes[0].coordinate_sha256
    assert len(reused_structures) == 1


def test_acquisition_mismatch_http_and_validation_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Sequence mismatch, missing model and malformed metadata should be auditable outcomes."""

    monkeypatch.setattr(
        af,
        "_request_bytes",
        lambda **_kwargs: json.dumps([_metadata(sequence="AAAA")]).encode(),
    )
    mismatch, structures = af.acquire_alphafold_models(
        requests=(AlphaFoldRequest("p1", "P11111"),),
        sequences=_sequences(),
        settings=_settings(cache=tmp_path / "cache", enabled=True),
    )
    assert mismatch[0].acquisition_status == "SEQUENCE_MISMATCH"
    assert structures == ()

    def missing(**_kwargs: object) -> bytes:
        raise urllib.error.HTTPError("url", 404, "missing", {}, io.BytesIO())

    monkeypatch.setattr(af, "_request_bytes", missing)
    unavailable, structures = af.acquire_alphafold_models(
        requests=(AlphaFoldRequest("p1", "P11111"),),
        sequences=_sequences(),
        settings=_settings(cache=tmp_path / "cache2", enabled=True),
    )
    assert unavailable[0].acquisition_status == "MODEL_NOT_AVAILABLE"
    assert structures == ()

    monkeypatch.setattr(af, "_request_bytes", lambda **_kwargs: b"not json")
    failed, structures = af.acquire_alphafold_models(
        requests=(AlphaFoldRequest("p1", "P11111"),),
        sequences=_sequences(),
        settings=_settings(cache=tmp_path / "cache3", enabled=True),
    )
    assert failed[0].acquisition_status == "FAILED"
    assert "metadata JSON" in failed[0].message
    assert structures == ()


def test_plddt_parser_uses_first_model_and_valid_ca_atoms(tmp_path: Path) -> None:
    """Only primary-altloc C-alpha B factors in the first model should be averaged."""

    path = tmp_path / "model.pdb"
    path.write_bytes(
        _pdb_line(atom="CA", alt=" ", plddt=80.0)
        + _pdb_line(atom="CA", alt="A", plddt=60.0)
        + _pdb_line(atom="CB", alt=" ", plddt=99.0)
        + _pdb_line(atom="CA", alt="B", plddt=10.0)
        + b"ATOM      malformed\nENDMDL\n"
        + _pdb_line(atom="CA", alt=" ", plddt=1.0)
    )
    assert af.parse_mean_plddt_from_pdb(path=path) == pytest.approx(70.0)
    empty = tmp_path / "empty.pdb"
    empty.write_text("HEADER\n", encoding="ascii")
    assert af.parse_mean_plddt_from_pdb(path=empty) is None
    with pytest.raises(InputValidationError, match="Could not parse"):
        af.parse_mean_plddt_from_pdb(path=tmp_path / "missing.pdb")


def test_metadata_and_url_validation_helpers() -> None:
    """Official metadata selection and URL allow-listing should fail closed."""

    selected = af._select_metadata(
        payload=json.dumps(
            [
                {"uniprotAccession": "P1", "latestVersion": 1},
                {"uniprotAccession": "P1", "latestVersion": 4},
                {"uniprotAccession": "P2", "latestVersion": 9},
            ]
        ).encode(),
        accession="p1",
    )
    assert selected["latestVersion"] == 4
    for payload in (b"{", b"{}", b"[]"):
        with pytest.raises(InputValidationError):
            af._select_metadata(payload=payload, accession="P1")
    assert af._validated_download_url(
        value="https://alphafold.ebi.ac.uk/files/model.pdb"
    ).startswith("https://")
    assert af._validated_download_url(value="https://ebi.ac.uk/files/model.pdb").startswith(
        "https://"
    )
    for value in ("", "http://alphafold.ebi.ac.uk/x", "https://evil.example/x"):
        with pytest.raises(InputValidationError):
            af._validated_download_url(value=value)
    assert af._metadata_url(accession="P1/unsafe").endswith("P1%2Funsafe")


def test_atomic_byte_writer_and_request_retries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Model caching should be atomic and transient network errors should retry."""

    path = tmp_path / "cache" / "model.pdb"
    af._write_bytes_atomic(path=path, payload=b"MODEL\n")
    assert path.read_bytes() == b"MODEL\n"
    with pytest.raises(InputValidationError):
        af._write_bytes_atomic(path=path, payload=b"")

    calls = {"count": 0}

    class Response:
        """Tiny context-managed URL response."""

        headers: dict[str, str] = {}

        def __enter__(self) -> Response:
            """Return this fake response."""

            return self

        def __exit__(self, *_args: object) -> None:
            """Close the fake response."""

        def read(self, _size: int = -1) -> bytes:
            """Return a non-empty response body."""

            return b"ok"

        def geturl(self) -> str:
            """Return the final response URL."""

            return "https://example.test/result"

    def urlopen(*_args: object, **_kwargs: object) -> Response:
        calls["count"] += 1
        if calls["count"] == 1:
            raise urllib.error.URLError("transient")
        return Response()

    monkeypatch.setattr(af.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(af.time, "sleep", lambda _seconds: None)
    assert af._request_bytes(url="https://example.test", timeout_seconds=1, retries=1) == b"ok"
    assert calls["count"] == 2

    monkeypatch.setattr(
        af.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            urllib.error.HTTPError("url", 400, "bad", {}, None)
        ),
    )
    with pytest.raises(urllib.error.HTTPError):
        af._request_bytes(url="https://example.test", timeout_seconds=1, retries=3)

    class EmptyResponse(Response):
        """Response carrying an invalid empty body."""

        def read(self, _size: int = -1) -> bytes:
            """Return an empty body."""

            return b""

    monkeypatch.setattr(af.urllib.request, "urlopen", lambda *_a, **_k: EmptyResponse())
    with pytest.raises(OSError, match="Empty response"):
        af._request_bytes(url="https://example.test", timeout_seconds=1, retries=0)

    class OversizedResponse(Response):
        """Response that exceeds its caller-provided byte ceiling."""

        def read(self, _size: int = -1) -> bytes:
            """Return more bytes than permitted."""

            return b"too large"

    monkeypatch.setattr(af.urllib.request, "urlopen", lambda *_a, **_k: OversizedResponse())
    with pytest.raises(OSError, match="ceiling"):
        af._request_bytes(
            url="https://example.test",
            timeout_seconds=1,
            retries=0,
            maximum_bytes=2,
        )

    class RedirectedResponse(Response):
        """Response redirected away from the permitted provider."""

        def geturl(self) -> str:
            """Return a disallowed final URL."""

            return "https://example.test/model.pdb"

    monkeypatch.setattr(af.urllib.request, "urlopen", lambda *_a, **_k: RedirectedResponse())
    with pytest.raises(InputValidationError, match="response URL"):
        af._request_bytes(
            url="https://alphafold.ebi.ac.uk/model.pdb",
            timeout_seconds=1,
            retries=0,
            allowed_host_suffix="ebi.ac.uk",
        )
    with pytest.raises(InputValidationError, match="byte ceiling"):
        af._request_bytes(
            url="https://example.test",
            timeout_seconds=1,
            retries=0,
            maximum_bytes=0,
        )


def test_outcome_serialisation() -> None:
    """Outcome helper should retain controlled identifiers and optional values."""

    outcome = af._outcome(
        request=AlphaFoldRequest("p1", "P11111"),
        status="FAILED",
        api_url="https://example.test",
        message="reason",
    )
    assert outcome.to_record()["coordinate_path"] == ""
    assert outcome.acquisition_status == "FAILED"


def _settings(*, cache: Path, enabled: bool, threshold: float = 50.0) -> AlphaFoldSettings:
    """Return deterministic acquisition settings for tests."""

    return AlphaFoldSettings(
        enabled=enabled,
        cache_dir=cache,
        timeout_seconds=1.0,
        retries=0,
        minimum_mean_plddt=threshold,
    )


def _sequences() -> tuple[SequenceRecord, ...]:
    """Return two authoritative test proteins."""

    return (
        SequenceRecord("p1", "", "ACDE", 4, sha256_text(text="ACDE")),
        SequenceRecord("p2", "", "FGHI", 4, sha256_text(text="FGHI")),
    )


def _metadata(*, sequence: str) -> dict[str, object]:
    """Return one official-API-shaped metadata record."""

    return {
        "uniprotAccession": "P11111",
        "entryId": "AF-P11111-F1",
        "latestVersion": 4,
        "uniprotSequence": sequence,
        "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P11111-F1-model_v4.pdb",
    }


def _pdb_bytes(*, plddt: float | None) -> bytes:
    """Return a minimal AlphaFold-style PDB payload."""

    return b"HEADER\n" if plddt is None else _pdb_line(atom="CA", alt=" ", plddt=plddt)


def _pdb_line(*, atom: str, alt: str, plddt: float) -> bytes:
    """Return one fixed-column PDB ATOM record."""

    line = (
        f"ATOM  {1:5d} {atom:>4s}{alt}ALA A{1:4d}    "
        f"{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{plddt:6.2f}          C  \n"
    )
    return line.encode("ascii")
