"""Exercise centroid evidence, remote exchange and exact subcluster handoffs."""

import io
import json
import zipfile
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from protein_signature_app.subcluster_network import (
    load_network_members,
    network_edges,
    network_figure,
    network_graphml,
)
from protein_signatures.errors import InputValidationError
from protein_signatures.workbench import subcluster_remote as remote
from protein_signatures.workbench import subclustering
from protein_signatures.workbench.analysis import iter_fasta
from protein_signatures.workbench.context import CONTEXT_FILE, context_bytes
from protein_signatures.workbench.jobs import artifact_path, read_job, run_next, submit_job
from protein_signatures.workbench.subcluster_review import (
    compare_partition,
    reviewed_taxon_sets,
    taxon_overlap,
)


@pytest.fixture(name="fake_diamond")
def fixture_fake_diamond(monkeypatch):
    """Create native and exact fixture outputs with one rescued singleton."""
    calls = []

    def command(*, config, arguments, output, stdout_name=None):
        """Emulate only the documented DIAMOND command outputs."""
        calls.append(arguments)
        operation = arguments[1]
        if operation == "version":
            (output / stdout_name).write_text("diamond version 2.2.8\n", encoding="utf-8")
        elif operation == "makedb":
            source = Path(arguments[arguments.index("--in") + 1])
            assert len(list(iter_fasta(path=source))) == 3
            assert source.parent.is_relative_to(Path(config.workspace))
        elif operation == "deepclust":
            target = Path(arguments[arguments.index("--out") + 1])
            target.write_text(
                "centroid\tmember\np000000001\tp000000001\n"
                "p000000001\tp000000002\np000000001\tp000000003\n",
                encoding="utf-8",
            )
        elif operation == "realign":
            pairs = Path(arguments[arguments.index("--clusters") + 1])
            assert not pairs.read_text(encoding="utf-8").startswith("centroid")
            target = Path(arguments[arguments.index("--out") + 1])
            target.write_text(
                "p000000001\tp000000002\t100\t9\t9\t1\t9\t1\t9\t9\t1e-20\t65\t8\n"
                "p000000001\tp000000003\t100\t9\t9\t1\t9\t1\t9\t9\t1e-20\t65\t3\n",
                encoding="utf-8",
            )
        else:
            raise AssertionError(operation)

    monkeypatch.setattr(subclustering, "run_command", command)
    return calls


def test_native_and_realignment_validation(tmp_path):
    """Do not accept missing members, rounded identity, invalid lengths or empty evidence."""
    aliases = {key: (key, 100) for key in ("centroid", "rounded", "length_error", "missing")}
    native = tmp_path / "native.tsv"
    native.write_text("centroid\tmember\ncentroid\tcentroid\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="every input"):
        subclustering.native_partition(path=native, aliases=aliases)
    exact = tmp_path / "exact.tsv"
    exact.write_text(
        "centroid\trounded\t70\t100\t100\t1\t99\t1\t99\t99\t1e-30\t60\t69\n"
        "centroid\tlength_error\t100\t90\t100\t1\t90\t1\t90\t90\t1e-30\t60\t90\n",
        encoding="utf-8",
    )
    rows = subclustering.audit_partition(
        native={key: "centroid" for key in aliases},
        aliases=aliases,
        exact=exact,
        controls=subclustering.subcluster_controls(parameters={}),
    )
    statuses = {row["sequence_id"]: row["alignment_status"] for row in rows}
    assert statuses == {
        "centroid": "CENTROID_SELF",
        "rounded": "EXACT_THRESHOLDS_FAIL",
        "length_error": "INVALID_REALIGNMENT",
        "missing": "NO_REALIGNMENT",
    }
    assert len({row["subcluster_id"] for row in rows}) == 4


def test_audited_partition_handoff_and_network(
    workbench_config, config_path, fake_diamond, fake_tools
):
    """Reuse a full audited subgroup with context and keep exactly the measured graph edges."""
    proteins = {"A:1": "ACDEFGHIK", "B|1": "ACDEFGHIL", "C": "MTTTFGHIL"}
    fasta = "".join(f">{key}\n{seq}\n" for key, seq in proteins.items()).encode()
    context = context_bytes(
        rows=[
            {
                "sequence_id": key,
                "sequence_sha256": sha256(seq.encode()).hexdigest(),
                "source_group": "N0.HOG001",
                "group_type": "HOG",
                "released_species": "Fixture species",
                "reviewed_taxon_id": 3702 if key != "C" else 9606,
            }
            for key, seq in proteins.items()
        ]
    )
    first = submit_job(
        config=workbench_config,
        operation="subcluster",
        uploads={"family": fasta, "context": context},
    )
    assert run_next(config=workbench_config) == first
    state = read_job(config=workbench_config, job_id=first)
    assert state["status"] == "COMPLETE", state.get("error")
    assert state["source_groups"] == ["N0.HOG001"]
    membership = artifact_path(config=workbench_config, job_id=first, name="membership.parquet")
    rows = pd.read_parquet(membership)
    chosen = rows.loc[rows.sequence_id.eq("B|1"), "subcluster_id"].item()
    rescued = rows.loc[rows.sequence_id.eq("C"), "subcluster_id"].item()
    assert chosen != rescued
    assert rows.loc[rows.sequence_id.eq("B|1"), "pident"].item() == pytest.approx(800 / 9)
    graph = load_network_members(path=membership, subclusters=[chosen, rescued])
    assert len(graph) == 3
    assert network_edges(members=graph).sequence_id.tolist() == ["B|1"]
    assert network_figure(members=graph).data
    document = ET.fromstring(network_graphml(members=graph))
    namespace = {"g": "http://graphml.graphdrawing.org/xmlns"}
    assert len(document.findall(".//g:node", namespaces=namespace)) == 3
    assert len(document.findall(".//g:edge", namespaces=namespace)) == 1
    context_path = artifact_path(config=workbench_config, job_id=first, name=CONTEXT_FILE)
    comparisons = compare_partition(
        config=workbench_config,
        membership=membership,
        context=context_path,
        subcluster_id=chosen,
    )
    assert comparisons["Size and sequence lengths"].proteins.tolist() == [3, 2]
    taxa = reviewed_taxon_sets(
        config=workbench_config,
        membership=membership,
        context=context_path,
        groups=[chosen, rescued],
    )
    assert taxa == {chosen: {"3702"}, rescued: {"9606"}}
    pairwise, regions = taxon_overlap(groups=taxa)
    assert pairwise[1]["jaccard"] == 0
    assert sorted(row["exclusive_reviewed_taxa"] for row in regions) == [0, 1, 1]

    second = submit_job(
        config=workbench_config,
        operation="align",
        sources={
            "family": {"job_id": first, "name": "sequences.fasta"},
            "context": {"job_id": first, "name": CONTEXT_FILE},
            "subcluster_membership": {
                "job_id": first,
                "name": "membership.parquet",
                "subcluster_id": chosen,
            },
        },
    )
    assert run_next(config=workbench_config) == second
    result = read_job(config=workbench_config, job_id=second)
    assert result["status"] == "COMPLETE", result.get("error")
    alignment_input = artifact_path(config=workbench_config, job_id=second, name="sequences.fasta")
    assert dict(iter_fasta(path=alignment_input)) == {
        "A:1": proteins["A:1"],
        "B|1": proteins["B|1"],
    }
    assert result["source_groups"] == ["N0.HOG001"]
    assert [call[1] for call in fake_diamond] == ["version", "makedb", "deepclust", "realign"]
    app = AppTest.from_string(
        "from pathlib import Path\n"
        "from protein_signature_app.workbench import render_workbench\n"
        f"render_workbench(config_path=Path({str(config_path)!r}))\n",
        default_timeout=15,
    ).run()
    app.radio(key="wb-stage").set_value("5. Research tools").run()
    app.selectbox(key="wb-research-tool").set_value("Audited subclustering").run()
    app.selectbox(key="wb-subcluster-review-id").set_value(chosen).run()
    app.button(key="wb-partition-open-Start MAFFT alignment").click().run()
    assert app.selectbox(key=f"wb-partition-{first}").value == chosen
    assert not app.exception and not app.error


def test_ssh_partition_round_trip_and_bad_archive(
    workbench_config,
    tmp_path,
    fake_diamond,
    monkeypatch,
):
    """Check the finite request, SSH invocation and response binding without a live VM."""
    config = replace(
        workbench_config,
        subcluster_ssh_host="pete@analysis-host",
        subcluster_ssh_executable="/opt/workbench/bin/protein-signature-workbench",
        subcluster_ssh_profile="/data/research/subcluster.toml",
    )
    vm = replace(workbench_config, workspace=str(tmp_path / "vm"))
    fasta = tmp_path / "family.fasta"
    fasta.write_bytes(b">A\nACDEFGHIK\n>B\nACDEFGHIL\n>C\nMTTTFGHIL\n")
    output = tmp_path / "results"
    output.mkdir()
    archive, request = tmp_path / "response.zip", tmp_path / "request.bin"

    def simulate_ssh(*, arguments, stdin_path, stdout_path, stdout_limit, **_kwargs):
        """Run the exact same protocol locally using a separate VM workspace."""
        assert arguments[0] == "ssh"
        assert "pete@analysis-host" in arguments
        assert arguments[-1] == "subcluster-stream"
        assert stdout_limit == remote.RESPONSE_LIMIT
        request.write_bytes(stdin_path.read_bytes())
        with stdin_path.open("rb") as source, stdout_path.open("wb") as destination:
            remote.write_stream_response(config=vm, source=source, destination=destination)
        archive.write_bytes(stdout_path.read_bytes())

    monkeypatch.setattr(remote, "run_command", simulate_ssh)
    controls = subclustering.subcluster_controls(parameters={})
    remote.subcluster_remote(config=config, source=fasta, output=output, controls=controls)
    assert (output / "membership.parquet").is_file()
    response_request = json.loads(request.read_bytes().split(b"\n", 1)[0])
    payload = io.BytesIO()
    with zipfile.ZipFile(archive) as original, zipfile.ZipFile(payload, "w") as changed:
        for item in original.infolist():
            value = original.read(item)
            changed.writestr(
                item.filename, b"x" + value[1:] if item.filename == "membership.tsv" else value
            )
    archive.write_bytes(payload.getvalue())
    tampered = tmp_path / "tampered"
    tampered.mkdir()
    with pytest.raises(InputValidationError, match="checksum"):
        remote._read_response(archive=archive, output=tampered, request=response_request)
    with pytest.raises(InputValidationError, match="truncated"):
        remote._read_request(
            stream=io.BytesIO(request.read_bytes()[:-2]),
            destination=tmp_path / "truncated.fasta",
        )


def test_subclustering_requires_consent_and_valid_thresholds(workbench_config, tmp_path):
    """Fail closed before executing a tool when SSH consent or thresholds are invalid."""
    with pytest.raises(InputValidationError, match="consent"):
        subclustering.subcluster(
            config=workbench_config,
            source=tmp_path / "absent",
            output=tmp_path,
            parameters={"transport": "ssh", "remote_allowed": "yes"},
        )
    with pytest.raises(InputValidationError, match="coverage"):
        subclustering.subcluster_controls(parameters={"mutual_coverage": 101})
