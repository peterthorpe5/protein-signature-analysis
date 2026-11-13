"""Verify released HOG scope, sequence-bound metadata and saved alignment analysis."""

import csv
import io
import json
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import duckdb
import pandas as pd
import pytest
from Bio import Phylo
from streamlit.testing.v1 import AppTest

from protein_signature_app.alignment_ui import context_for_alignment, identity_figure
from protein_signature_app.research_tools import portable_model_request
from protein_signatures.checksums import sha256_file
from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.alignment_study import auspice_dataset, pairwise_identity
from protein_signatures.workbench.context import CONTEXT_FILE, carry_context
from protein_signatures.workbench.hogs import (
    HogScope,
    hog_members,
    hog_snapshot,
    recorded_domains,
    search_hogs,
)
from protein_signatures.workbench.jobs import artifact_path, read_job, run_next, submit_job
from protein_signatures.workbench.species_tree import species_display_tree


@pytest.fixture(name="hog_database")
def fixture_hog_database(tmp_path: Path) -> Path:
    """Build two published runs sharing one HOG ID, with missing sequences and features."""
    database = tmp_path / "published.duckdb"
    with duckdb.connect(str(database)) as connection:
        connection.execute(
            "CREATE TABLE proteins(protein_id VARCHAR, sequence VARCHAR, sequence_sha256 VARCHAR)"
        )
        connection.executemany(
            "INSERT INTO proteins VALUES (?,?,?)",
            [
                (key, seq, sha256(seq.encode()).hexdigest())
                for key, seq in (
                    ("B", "ACDEFGHIK"),
                    ("C", "ACDEFGHIL"),
                    ("UNIPROT_P1", "ACDEFGHIK"),
                )
            ],
        )
        connection.execute("""
            CREATE TABLE orthofinder_memberships(
                run_id VARCHAR, group_type VARCHAR, hierarchy_node VARCHAR,
                group_id VARCHAR, species_label VARCHAR, protein_id VARCHAR
            );
            INSERT INTO orthofinder_memberships VALUES
                ('r1','HOG','N0','N0.HOG001','Published plant','B'),
                ('r1','HOG','N0','N0.HOG001','Published human','C'),
                ('r1','HOG','N0','N0.HOG001','Published plant','UNIPROT_P1'),
                ('r1','HOG','N0','N0.HOG001','Unresolved species','missing'),
                ('r2','HOG','N1','N0.HOG001','Other run','C'),
                ('r1','LEGACY_ORTHOGROUP','','OG001','Published plant','B');
            CREATE TABLE domain_hits(
                protein_id VARCHAR, domain_authority VARCHAR, domain_id VARCHAR,
                domain_name VARCHAR, "start" BIGINT, "end" BIGINT,
                evidence_source VARCHAR, evidence_reference VARCHAR
            );
            INSERT INTO domain_hits VALUES ('B','Pfam','PF00001','Fixture',2,6,'published','record')
        """)
    return database


def test_hog_scope_and_context(workbench_config, hog_database, tmp_path):
    """Use the exact run, and keep evidence only for unchanged IDs and sequences."""
    original = sha256_file(path=hog_database)
    first, second = search_hogs(database=hog_database, text="hog001")
    assert (first.run_id, second.hierarchy_node) == ("r1", "N1")
    frame = hog_members(database=hog_database, scope=first, maximum=10)
    assert list(frame.protein_id) == ["B", "C", "UNIPROT_P1", "missing"]
    assert frame.loc[frame.protein_id.eq("missing"), "sequence"].item() == ""
    assert hog_members(database=hog_database, scope=second, maximum=10).protein_id.tolist() == ["C"]
    with pytest.raises(InputValidationError, match="outside"):
        hog_members(database=hog_database, scope=first, maximum=10, selected_ids=["not_a_member"])
    assert recorded_domains(database=hog_database, identifier="B")[0]["start"] == 2
    assert recorded_domains(database=hog_database, identifier="C") == []
    fasta, data = hog_snapshot(config=workbench_config, scope=first, members=frame)
    assert b">missing" not in fasta
    records = {"B": "ACDEFGHIK", "C": "ACDEFGHIL"}
    assert context_for_alignment(data=data, records=records)["B"]["reviewed_taxon_id"] == "3702"
    assert (
        context_for_alignment(data=data, records=records)["B"]["released_species"]
        == "Published plant"
    )

    context = tmp_path / "context.tsv"
    context.write_bytes(data)
    changed = tmp_path / "changed.fasta"
    changed.write_bytes(b">B\nACDEFGHIK\n>C\nACDEFGHIK\n")
    output = tmp_path / "carry"
    output.mkdir()
    carry_context(config=workbench_config, context=context, sequences=changed, output=output)
    saved = pd.read_csv(output / CONTEXT_FILE, sep="\t")
    assert saved.sequence_id.tolist() == ["B"]
    assert sha256_file(path=hog_database) == original


def test_reviewed_species_mapping_conflicts(workbench_config, hog_database, tmp_path):
    """Accept exact reviewed species labels and leave contradictory taxon evidence unknown."""
    mapping = tmp_path / "reviewed.tsv"
    mapping.write_text(
        "workflow_species_label\taccepted_species_name\tncbi_taxon_id\tmapping_status\n"
        "Published plant\tPlant species\t42\tREVIEWED\n"
        "Published human\tHuman\t9606\tPENDING\n",
        encoding="utf-8",
    )
    config = replace(workbench_config, species_taxonomy=str(mapping))
    scope = search_hogs(database=hog_database, text="001")[0]
    frame = hog_members(database=hog_database, scope=scope, maximum=5)
    _, data = hog_snapshot(config=config, scope=scope, members=frame)
    records = {
        row["sequence_id"]: row
        for row in csv.DictReader(
            io.StringIO(data.decode()),
            delimiter="\t",
        )
    }
    assert records["B"]["reviewed_taxon_id"] == ""
    assert records["C"]["reviewed_taxon_id"] == "9606"
    assert records["B"]["released_species"] == "Published plant"
    no_sequence_taxonomy = replace(config, taxonomy="")
    _, data = hog_snapshot(config=no_sequence_taxonomy, scope=scope, members=frame)
    assert b"Plant species\t42" in data
    assert b"\tHuman\t9606" not in data


def test_saved_alignment_tree_does_not_repeat_mafft(workbench_config, hog_database, fake_tools):
    """Create one MAFFT job and infer a separate tree from its immutable artifact."""
    scope = search_hogs(database=hog_database, text="001")[0]
    frame = hog_members(database=hog_database, scope=scope, maximum=4)
    fasta, context = hog_snapshot(config=workbench_config, scope=scope, members=frame)
    first = submit_job(
        config=workbench_config,
        operation="align",
        uploads={"family": fasta, "context": context},
        parameters={"reference_id": "C"},
    )
    assert run_next(config=workbench_config) == first
    state = read_job(config=workbench_config, job_id=first)
    assert state["status"] == "COMPLETE", state.get("error")
    assert state["source_groups"] == [scope.group_id]
    files = state["artifacts"]
    assert "pairwise_identity.tsv" in files
    alignment = artifact_path(config=workbench_config, job_id=first, name="alignment.fasta")
    metrics = pd.read_csv(alignment.parent / "alignment_members.tsv", sep="\t")
    assert metrics.reference_id.unique().tolist() == ["C"]
    second = submit_job(
        config=workbench_config,
        operation="tree",
        sources={
            "alignment": {"job_id": first, "name": "alignment.fasta"},
            "context": {"job_id": first, "name": CONTEXT_FILE},
        },
    )
    assert run_next(config=workbench_config) == second
    assert read_job(config=workbench_config, job_id=first)["artifacts"] == files
    assert sha256_file(path=alignment) == files["alignment.fasta"]["sha256"]
    assert [call[0] for call in fake_tools] == [workbench_config.mafft, workbench_config.fasttree]
    tree = json.loads(
        artifact_path(
            config=workbench_config,
            job_id=second,
            name="auspice.json",
        ).read_text(encoding="utf-8")
    )
    assert tree["version"] == "v2"
    assert {node["name"] for node in tree["tree"]["children"]} == {"B", "C", "UNIPROT_P1"}
    assert "num_date" not in json.dumps(tree)
    assert "reviewed_taxon_id" in tree["tree"]["children"][0]["node_attrs"]


def test_identity_and_species_display():
    """Keep no-overlap undefined and avoid treating same-species display edges as gene history."""
    records = {"a": "AC-X", "b": "ATX-", "c": "XX-X"}
    pairs = pairwise_identity(records=records)
    by_pair = {(row["sequence_a"], row["sequence_b"]): row for row in pairs}
    assert by_pair["a", "b"]["identity"] == 0.5
    assert by_pair["a", "b"]["comparable_residue_count"] == 2
    assert by_pair["a", "c"]["identity"] is None
    assert identity_figure(selected=list(records), pairs=pairs).data[0].customdata[0][1] == 2
    metadata = {
        "a": {"reviewed_taxon_id": "3702"},
        "b": {"reviewed_taxon_id": "3702"},
        "c": {"released_species": "Human"},
    }
    newick, audit = species_display_tree(
        newick="((3702:1,Human:1):2,unused:3);",
        sequence_ids=list(records),
        metadata=metadata,
    )
    tree = Phylo.read(io.StringIO(newick), "newick")
    assert {tip.name for tip in tree.get_terminals()} == set(records)
    assert all(tip.branch_length == 0 for tip in tree.get_terminals())
    assert audit[0]["species_tree_tip"] == "3702"
    with pytest.raises(InputValidationError, match="one exact"):
        species_display_tree(newick="(3702:1,Human:1);", sequence_ids=list(records), metadata={})
    with pytest.raises(InputValidationError, match="tips"):
        auspice_dataset(newick="(bad:1,b:1);", labels={"a": "a", "b": "b"})


def test_hog_alignment_and_feature_routes(config_path, hog_database, fake_tools):
    """Exercise the published HOG selector and the no-feature path through Streamlit."""
    app = AppTest.from_string(
        "from pathlib import Path\n"
        "from protein_signature_app.workbench import render_workbench\n"
        f"render_workbench(config_path=Path({str(config_path)!r}),"
        f" database=Path({str(hog_database)!r}))\n",
        default_timeout=15,
    ).run()
    app.radio(key="wb-stage").set_value("5. Research tools").run()
    app.radio(key="wb-jb-source").set_value("Align a published HOG").run()
    scope = HogScope(group_id="N0.HOG001", run_id="r1", hierarchy_node="N0")
    app.selectbox(key="wb-jb-hog-hog").set_value(scope).run()
    app.text_input(key="wb-jb-hog-mafft-ref").set_value("B").run()
    app.button(key="wb-jb-hog-mafft-submit").click().run()
    assert not app.exception
    assert not app.error
    assert app.session_state["wb-jb-hog-mafft-last"].startswith("job_")
    app.selectbox(key="wb-research-tool").set_value("Nightingale features").run()
    app.radio(key="wb-feature-source").set_value("Published HOG features").run()
    app.selectbox(key="wb-feature-hog-hog").set_value(scope).run()
    app.selectbox(key="wb-feature-hog-protein").set_value("C").run()
    assert not app.exception
    assert any("No domain intervals" in message.value for message in app.info)


def test_portable_hog_model_includes_only_exact_selected_context(workbench_config, hog_database):
    """Prepare a one-protein ESMFold handoff without contacting a service."""
    import zipfile

    scope = search_hogs(database=hog_database, text="001")[0]
    frame = hog_members(database=hog_database, scope=scope, maximum=4)
    fasta, context = hog_snapshot(config=workbench_config, scope=scope, members=frame)
    selected, bundle = portable_model_request(
        config=workbench_config,
        selected=["C"],
        uploads={"family": fasta, "context": context},
        sources={},
    )
    assert selected == b">C\nACDEFGHIL\n"
    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        assert b"Published human" in archive.read(CONTEXT_FILE)
        assert b"Published plant" not in archive.read(CONTEXT_FILE)


def test_analyse_family_alignment_queues_exact_seed_input(config_path, workbench_config):
    """Keep the original Analyse family route while sharing the new alignment controls."""
    app = AppTest.from_string(
        "from pathlib import Path\n"
        "from protein_signature_app.workbench import render_workbench\n"
        f"render_workbench(config_path=Path({str(config_path)!r}))\n",
        default_timeout=15,
    ).run()
    app.radio(key="wb-stage").set_value("4. Analyse family").run()
    app.selectbox(key="wb-analysis-task").set_value("Alignment and tree").run()
    assert app.button(key="wb-family-mafft-submit").disabled
    app.session_state["wb-seeds"] = b">one\nACD\n>two\nACE\n"
    app.run()
    app.selectbox(key="wb-family-source").set_value("Current seed snapshot").run()
    app.text_input(key="wb-family-mafft-ref").set_value("two").run()
    app.button(key="wb-family-mafft-submit").click().run()
    identifier = app.session_state["wb-family-mafft-last"]
    assert read_job(config=workbench_config, job_id=identifier)["operation"] == "align"
    assert not app.exception
