"""Check portable results, explicit server profiles and the packaged viewers."""

import hashlib
import json
import os
import shutil

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from protein_signature_app import research_viewers
from protein_signatures.errors import InputValidationError
from protein_signatures.workbench import exchange
from protein_signatures.workbench.cli import main
from protein_signatures.workbench.config import configure_mount, https_endpoint, load_config
from protein_signatures.workbench.jobs import artifact_path, read_job, run_next, submit_job


@pytest.mark.parametrize(
    "value",
    [
        "http://test.example",
        "https://user:pass@test.example",
        "https://test.example?token=1",
        "https://test.example/#fragment",
        "https://test.example:abc",
        "https:///path",
    ],
)
def test_provider_configuration_rejects_invalid_endpoints(value):
    """Keep credential-bearing and non-HTTPS endpoints out of worker configuration."""
    with pytest.raises(InputValidationError):
        https_endpoint(value=value)


def test_stefan_profile_and_cli_keep_sources_outside_workspace(workbench_config, tmp_path, capsys):
    """Point a separate profile at the mounted release, with no database copying."""
    mount = tmp_path / "ken_deepclust_data"
    mount.mkdir()
    for name, original in (
        ("clust_index_RowGroup.parquet", workbench_config.index),
        ("joined_with_index_RowGroupFinal.parquet", workbench_config.sequences),
        ("SeqIdMapClustId.parquet", workbench_config.sequence_map),
    ):
        shutil.copyfile(original, mount / name)
    workspace, profile = tmp_path / "jobs", tmp_path / "config" / "workbench.toml"
    flags = ["--config", str(profile)]
    assert (
        main(
            argv=flags
            + [
                "--action",
                "configure-mount",
                "--workspace",
                str(workspace),
                "--release",
                "installed fixture",
                "--deepclust-root",
                str(mount),
            ]
        )
        == 0
    )
    config = load_config(path=profile)
    assert config.sequence_map == str(mount / "SeqIdMapClustId.parquet")
    assert main(argv=flags + ["--action", "workspace"]) == 0
    assert capsys.readouterr().out.splitlines()[-1] == str(workspace)
    with pytest.raises(InputValidationError, match="already exists"):
        configure_mount(root=mount, destination=profile, workspace=workspace, release="fixture")
    with pytest.raises(InputValidationError, match="outside"):
        configure_mount(
            root=mount, destination=profile, workspace=mount / "jobs", release="fixture"
        )
    (mount / "SeqIdMapClustId.parquet").unlink()
    with pytest.raises(InputValidationError, match="Missing"):
        configure_mount(
            root=mount, destination=tmp_path / "new.toml", workspace=workspace, release="fixture"
        )


def test_portable_deepclust_export_feeds_complete_fasta(workbench_config, tmp_path):
    """Import an exact cluster selection and keep a portable-source provenance record."""
    source = tmp_path / "export.parquet"
    pq.write_table(
        pa.table(
            {
                "id": ["X", "X", "Y", "Z"],
                "aa": ["AAAA", "AAAA", "GGGG", "VVVV"],
                "cluster": ["c", "c", "c", "d"],
            }
        ),
        source,
    )
    payload = source.read_bytes()
    job = submit_job(
        config=workbench_config,
        operation="import_deepclust",
        uploads={"deepclust_export": payload},
        parameters={
            "id_column": "id",
            "sequence_column": "aa",
            "cluster_column": "cluster",
            "cluster_id": "c",
        },
    )
    assert run_next(config=workbench_config) == job
    assert read_job(config=workbench_config, job_id=job)["status"] == "COMPLETE"
    assert artifact_path(config=workbench_config, job_id=job, name="family.fasta").read_text() == (
        ">X\nAAAA\n>Y\nGGGG\n"
    )
    audit = artifact_path(config=workbench_config, job_id=job, name="import_summary.json")
    assert json.loads(audit.read_text())["duplicate_rows_collapsed"] == 1
    assert source.read_bytes() == payload


def test_portable_import_rejects_invalid_or_ambiguous_records(tmp_path):
    """Do not use corrupt Parquet, conflicting IDs or unmatched cluster selections."""
    with pytest.raises(InputValidationError, match="Parquet"):
        exchange.result_reader(source=b"Not Parquet")
    source = tmp_path / "export.parquet"
    pq.write_table(pa.table({"sequence_id": ["A", "A"], "sequence": ["AA", "GG"]}), source)
    for parameters, match in (
        ({}, "Conflicting"),
        ({"id_column": "missing"}, "absent"),
        ({"id_column": "sequence"}, "distinct"),
    ):
        with pytest.raises(InputValidationError, match=match):
            exchange.import_result(source=source, parameters=parameters, output=tmp_path)


def test_msa_exact_identity_and_explicit_tree_mapping():
    """View legacy FastTree outputs without guessing a relationship between IDs."""
    data = research_viewers.alignment_payload(
        fasta=">UNIPROT_A|v1\nAC-D\n>plant-B\nACED\n",
        newick="(seq_0:0.2,seq_1:0.4);",
        tip_labels={"seq_0": "UNIPROT_A|v1", "seq_1": "plant-B"},
    )
    assert "UNIPROT_A|v1" in data["msa"] and "UNIPROT_A|v1" in data["tree"]
    assert ":0.4" in data["tree"]
    for options in (
        {"fasta": ">a\nAA\n>a\nAG\n"},
        {"fasta": ">a\nAA\n>b\nAGG\n"},
        {"fasta": data["msa"], "newick": "(wrong:0.1,plant-B:0.2);"},
        {"fasta": data["msa"], "newick": "(UNIPROT_A|v1:-1,plant-B:0.2);"},
        {"fasta": data["msa"], "newick": "(seq_0,seq_1);", "tip_labels": {"seq_0": "plant-B"}},
    ):
        with pytest.raises(InputValidationError):
            research_viewers.alignment_payload(**options)


def test_nightingale_exact_hash_and_inert_payload():
    """Keep coordinates tied to their sequence; escape markup and skip orthology rows."""
    digest = hashlib.sha256(b"ACDE").hexdigest()
    rows = [
        {
            "sequence_id": "protein",
            "sequence_sha256": digest,
            "start": "2",
            "end": "3",
            "description": "</script><script>untrusted()</script>",
            "accession": "PF00001",
        },
        {"sequence_id": "protein", "sequence_sha256": digest, "start": "", "end": ""},
    ]
    payload = research_viewers.feature_payload(identifier="protein", sequence="ACDE", rows=rows)
    assert len(payload["features"]) == 1
    document = research_viewers.viewer_document(viewer="nightingale", payload=payload)
    assert "</script><script>untrusted()" not in document and "\\u003c" in document
    assert "connect-src 'none'" in document
    for wrong in ({"sequence_sha256": "wrong"}, {"start": "1.5"}, {"end": 5}):
        with pytest.raises(InputValidationError):
            research_viewers.feature_payload(
                identifier="protein", sequence="ACDE", rows=[{**rows[0], **wrong}]
            )


@pytest.mark.skipif(not os.getenv("PS_BROWSER_EXE"), reason="Requires a local Chromium executable")
def test_real_browser_renders_both_viewers_offline(tmp_path):
    """Exercise the upstream JS bundles in Chromium with network access disabled."""
    from playwright.sync_api import sync_playwright

    data = {
        "jbrowse": research_viewers.alignment_payload(
            fasta=">a\nACDEFGHIK\n>b\nAC-EFGHIK\n", newick="(a:0.1,b:0.2);"
        ),
        "nightingale": research_viewers.feature_payload(
            identifier="a",
            sequence="ACDEFGHIK",
            rows=[
                {
                    "sequence_id": "a",
                    "sequence_sha256": hashlib.sha256(b"ACDEFGHIK").hexdigest(),
                    "start": "2",
                    "end": "5",
                    "accession": "PF00001",
                    "description": "Example domain",
                }
            ],
        ),
    }
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=os.environ["PS_BROWSER_EXE"], headless=True, args=["--no-sandbox"]
        )
        try:
            page = browser.new_page(offline=True)
            for name, payload in data.items():
                path = tmp_path / f"{name}.html"
                path.write_text(research_viewers.viewer_document(viewer=name, payload=payload))
                errors = []
                page.on("pageerror", lambda error: errors.append(str(error)))
                page.goto(path.as_uri())
                if name == "jbrowse":
                    page.wait_for_function("window.proteinMSA !== undefined", timeout=10000)
                    page.locator("#viewer canvas").first.wait_for(timeout=10000)
                else:
                    page.locator("nightingale-track").first.wait_for(timeout=10000)
                    assert "Example domain" in page.locator("#viewer").inner_text()
                assert not errors, errors
        finally:
            browser.close()
