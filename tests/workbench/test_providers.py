"""Exercise exact-sequence prediction and annotation without transmitting data."""

import csv
import hashlib
import json
from dataclasses import replace
from io import StringIO

import pytest
import requests

from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.workbench import providers
from protein_signatures.workbench.jobs import artifact_path, read_job, run_next, submit_job


def interpro_tsv(*, identifier, sequence, start=1, end=3):
    """Create an InterPro-compatible hit bound to the submitted sequence MD5."""
    digest = hashlib.md5(sequence.encode(), usedforsecurity=False).hexdigest()
    return (
        "\t".join(
            [
                identifier,
                digest,
                str(len(sequence)),
                "Pfam",
                "PF00001",
                "Model domain",
                str(start),
                str(end),
                "1e-8",
                "T",
                "date",
                "IPR1",
                "domain",
                "GO:0000001",
            ]
        )
        + "\n"
    )


def test_selection_is_exact_and_bounded(tmp_path):
    """Reject missing IDs, duplicate IDs and incomplete implicit selection."""
    path = tmp_path / "proteins.fasta"
    path.write_text(">one\nACDEF\n>two\nMKAAA\n")
    kwargs = dict(source=path, maximum=1, scan_limit=10)
    assert providers.select_proteins(selected_ids=["two"], **kwargs) == {"two": "MKAAA"}
    for ids in ([], ["ONE"], ["one", "two"]):
        with pytest.raises(InputValidationError):
            providers.select_proteins(selected_ids=ids, **kwargs)
    path.write_text(">one\nACDEF\n>one\nMKAAA\n")
    with pytest.raises(InputValidationError, match="Duplicate"):
        providers.select_proteins(selected_ids=["one"], **kwargs)
    kwargs["scan_limit"] = 1
    with pytest.raises(InputValidationError, match="scan bound"):
        providers.select_proteins(selected_ids=["one"], **kwargs)


@pytest.mark.parametrize("status,size", [(200, 5), (200, 50), (302, 5), (503, 0)])
def test_provider_transport_is_bounded_and_refuses_redirects(monkeypatch, status, size):
    """Reject non-success responses and excessive bytes with no redirect follow."""

    class Response:
        """Minimal requests stream boundary."""

        status_code = status

        def __enter__(self):
            """Return the fake response."""
            return self

        def __exit__(self, *_args):
            """Close the fake response."""

        def iter_content(self, *, chunk_size):
            """Yield the proposed response bytes."""
            assert chunk_size == 65536
            yield b"x" * size

    def request(**kwargs):
        assert kwargs["allow_redirects"] is False and kwargs["stream"] is True
        return Response()

    monkeypatch.setattr(providers.requests, "request", request)
    if status == 200 and size <= 10:
        assert (
            providers.provider_request(url="https://provider.example/endpoint", maximum=10)
            == b"x" * size
        )
    else:
        with pytest.raises(ExternalToolError):
            providers.provider_request(url="https://provider.example/endpoint", maximum=10)


def test_provider_transport_error_is_reported_without_credentials(monkeypatch):
    """Do not include the body, token or raw HTTP error in the display message."""

    def offline(**_kwargs):
        raise requests.ConnectionError("private service token")

    monkeypatch.setattr(providers.requests, "request", offline)
    with pytest.raises(ExternalToolError, match="Could not reach") as failure:
        providers.provider_request(url="https://provider.example/predict", body=b"private")
    assert "private" not in str(failure.value)


@pytest.mark.parametrize("mode", ["sequence", "json_inputs"])
def test_model_job_validates_sequence_and_retains_model_provenance(
    workbench_config,
    monkeypatch,
    mode,
):
    """The persisted API request is explicitly selected and returns a checked monomer."""
    model = (
        b"ATOM      1  CA  ALA A   1       2.000   1.000   2.000  1.00 85.00           C\n"
        b"ATOM      2  CA  ALA A   2       4.000   1.000   2.000  1.00 86.00           C\nEND\n"
    )
    config = replace(
        workbench_config,
        esm_endpoint="https://provider.example/predict",
        esm_token_env="WORKBENCH_TEST_TOKEN",
        esm_format=mode,
    )
    monkeypatch.setenv("WORKBENCH_TEST_TOKEN", "fixture-token")
    calls = []

    def respond(**kwargs):
        calls.append(kwargs)
        return json.dumps({"pdb": model.decode()}).encode() if mode == "json_inputs" else model

    monkeypatch.setattr(providers, "provider_request", respond)
    identifier = submit_job(
        config=config,
        operation="predict_models",
        uploads={"family": b">seed\nAA\n"},
        parameters={"remote_allowed": True},
    )
    assert run_next(config=config) == identifier
    state = read_job(config=config, job_id=identifier)
    assert state["status"] == "COMPLETE", state
    path = artifact_path(config=config, job_id=identifier, name="models/model_0001.pdb")
    assert path.read_bytes() == model
    assert calls[0]["headers"]["Authorization"] == "Bearer fixture-token"
    body = calls[0]["body"]
    assert (json.loads(body)["inputs"] if mode == "json_inputs" else body.decode()) == "AA"
    assert "fixture-token" not in path.parent.parent.joinpath("esmfold_provenance.json").read_text()


@pytest.mark.parametrize(
    "sequence,parameters,endpoint,token,match",
    [
        ("AA", {}, "https://provider.example/predict", "", "explicit"),
        ("AAA", {"remote_allowed": True}, "https://provider.example/predict", "", "up to 2"),
        ("AA", {"remote_allowed": True}, "", "", "configured"),
        (
            "AA",
            {"remote_allowed": True},
            "https://provider.example/predict",
            "NO_TEST_TOKEN",
            "missing",
        ),
    ],
)
def test_model_preflight_does_not_send_on_invalid_input(
    workbench_config,
    tmp_path,
    monkeypatch,
    sequence,
    parameters,
    endpoint,
    token,
    match,
):
    """Reject missing consent, credentials and excess length before any transfer."""
    config = replace(workbench_config, esm_endpoint=endpoint, esm_max_length=2, esm_token_env=token)
    source = tmp_path / "source.fasta"
    source.write_text(f">protein\n{sequence}\n")
    monkeypatch.delenv("NO_TEST_TOKEN", raising=False)
    monkeypatch.setattr(
        providers, "provider_request", lambda **_: pytest.fail("Unexpected API request")
    )
    with pytest.raises(InputValidationError, match=match):
        providers.predict_models(
            config=config, source=source, parameters=parameters, output=tmp_path
        )


@pytest.mark.parametrize(
    "content",
    [
        "unknown\tbad\trow",
        interpro_tsv(identifier="protein", sequence="GGGG"),
        interpro_tsv(identifier="protein", sequence="AAAA", end=5),
        interpro_tsv(identifier="protein", sequence="AAAA", start="1.5"),
    ],
)
def test_invalid_interpro_coordinates_and_identity(content):
    """Keep incompatible rows away from the residue viewer."""
    with pytest.raises(InputValidationError):
        providers.interpro_rows(content=content, proteins={"protein": "AAAA"})


@pytest.mark.parametrize("engine", ["interpro_local", "interpro_remote", "eggnog"])
def test_annotation_keeps_no_hits_and_nonpositional_evidence(
    workbench_config,
    tmp_path,
    monkeypatch,
    engine,
):
    """Persist selected hits and no-hit rows with exact FASTA identifiers."""
    config = replace(workbench_config, eggnog_data=str(tmp_path), interpro_email="test@example.org")
    source = tmp_path / "proteins.fasta"
    source.write_text(">real|identifier\nAAAA\n>no-hit\nGGGG\n")
    output = tmp_path / "annotation"
    output.mkdir()

    def binary(*, config, arguments, output):
        if engine == "eggnog":
            assert "--tax_scope" in arguments and "--target_taxa" in arguments
            (output / "annotation.emapper.annotations").write_text(
                "#query\tseed_ortholog\tDescription\tscore\tGOs\np000001\tO1\tPredicted\t10\tGO:1\n"
            )
        else:
            (output / "interpro.tsv").write_text(
                interpro_tsv(identifier="p000001", sequence="AAAA")
            )

    counter = []

    def ebi(*, url, **_kwargs):
        if url.endswith("/run"):
            counter.append(1)
            return f"job_{len(counter)}".encode()
        if "/status/" in url:
            return b"FINISHED"
        return (
            interpro_tsv(identifier="p000001", sequence="AAAA").encode() if "job_1" in url else b""
        )

    monkeypatch.setattr(providers, "run_command", binary)
    monkeypatch.setattr(providers, "provider_request", ebi)
    providers.annotate_proteins(
        config=config,
        source=source,
        output=output,
        parameters={"engine": engine, "target_taxon": "2759", "remote_allowed": True},
    )
    rows = list(csv.DictReader(StringIO((output / "annotations.tsv").read_text()), delimiter="\t"))
    assert rows[0]["sequence_id"] == "real|identifier"
    assert rows[0]["sequence_sha256"] == hashlib.sha256(b"AAAA").hexdigest()
    assert rows[0]["start"] == ("" if engine == "eggnog" else "1")
    assert "ASSESSED_NO_HIT" in (output / "annotation_audit.tsv").read_text()
    if engine == "interpro_remote":
        assert "job_2" in (output / "interpro_service_jobs.jsonl").read_text()


def test_annotation_requires_consent_and_configured_taxonomy(
    workbench_config, tmp_path, monkeypatch
):
    """Prevent remote transmission and execution with invalid provider settings."""
    source = tmp_path / "proteins.fasta"
    source.write_text(">protein\nAAAA\n")
    monkeypatch.setattr(
        providers, "provider_request", lambda **_: pytest.fail("Unexpected API request")
    )
    for parameters, match in (
        ({"engine": "interpro_remote"}, "consent"),
        ({"engine": "unknown"}, "supported"),
        ({"engine": "eggnog", "tax_scope": "unverified scope"}, "taxonomy"),
        ({"engine": "interpro_remote", "remote_allowed": True}, "interpro_email"),
    ):
        with pytest.raises(InputValidationError, match=match):
            providers.annotate_proteins(
                config=workbench_config, source=source, output=tmp_path, parameters=parameters
            )
