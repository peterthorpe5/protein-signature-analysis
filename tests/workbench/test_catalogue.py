"""Real Parquet tests for exact mapping, indexed extraction and input corruption."""

from dataclasses import replace

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from protein_signatures.errors import InputValidationError
from protein_signatures.workbench.catalogue import (
    column,
    identifiers,
    normalise_sequence,
    resolve_ids,
    retrieve_clusters,
    row_groups,
)


def test_mapping_preserves_unmatched_ambiguous_unindexed_and_versioned_ids(workbench_config):
    """Exact membership and alias evidence must never be replaced with invented matches."""
    frame = resolve_ids(
        config=workbench_config, values=["sp|P1|NAME", "local_b", "multi", "F", "missing", "NR_D.2"]
    )
    assert frame.loc[frame.input_id == "sp|P1|NAME", "sequence_id"].tolist() == ["UNIPROT_P1"]
    assert frame.loc[frame.input_id == "local_b", "mapping_evidence"].tolist() == [
        "DECLARED_BRIDGE"
    ]
    assert len(frame.loc[frame.input_id == "multi"]) == 2
    assert frame.loc[frame.input_id == "F", "status"].tolist() == ["MAPPED_UNINDEXED"]
    assert set(frame.loc[frame.input_id.isin(["missing", "NR_D.2"]), "status"]) == {"UNMATCHED"}


def test_complete_retrieval_crosses_row_groups_without_unrelated_records(
    workbench_config, tmp_path
):
    """Physical groups can overlap clusters; exact filtering and counts are both necessary."""
    output = tmp_path / "extract"
    output.mkdir()
    retrieve_clusters(config=workbench_config, values=["UNIPROT_P1"], output=output)
    text = (output / "family.fasta").read_text()
    assert text.count(">") == 3 and "*" not in text and "NR_D.1" not in text
    summary = pd.read_csv(output / "clusters.tsv", sep="\t")
    assert summary.loc[0, "retrieved_members"] == summary.loc[0, "published_members"] == 3
    assert summary.loc[0, "physical_row_groups_read"] == 2
    assert (output / "taxonomy_summary.tsv").is_file()


@pytest.mark.parametrize("values", [[], ["unknown"], ["F"]])
def test_unavailable_cluster_is_rejected(workbench_config, tmp_path, values):
    """Unindexed or unknown clusters are explicit failures, never partial extractions."""
    output = tmp_path / "bad"
    output.mkdir()
    with pytest.raises(InputValidationError):
        retrieve_clusters(config=workbench_config, values=values, output=output)


def test_count_bound_and_mismatched_index_fail(workbench_config, tmp_path):
    """Complete counts are enforced before and after reading physical records."""
    output = tmp_path / "bounded"
    output.mkdir()
    with pytest.raises(InputValidationError, match="exceed"):
        retrieve_clusters(
            config=replace(workbench_config, max_sequences=2), values=["UNIPROT_P1"], output=output
        )
    index = tmp_path / "wrong.parquet"
    pq.write_table(
        table=pa.table(data={"CLUSTER": ["UNIPROT_P1"], "NMEMBER": [4], "RowGroup": ["0,1"]}),
        where=index,
    )
    with pytest.raises(InputValidationError, match="disagree"):
        retrieve_clusters(
            config=replace(workbench_config, index=str(index)), values=["UNIPROT_P1"], output=output
        )


@pytest.mark.parametrize("value", [None, "bad", "-1", "8", "[]"])
def test_invalid_row_groups_fail(value):
    """Malformed offsets cannot accidentally trigger a full database scan."""
    with pytest.raises(InputValidationError):
        row_groups(value=value, maximum=3)


def test_validators_and_supported_aliases():
    """Namespace punctuation and accession versions survive strict normalisation."""
    assert row_groups(value="[0, 1;1]", maximum=3) == (0, 1)
    assert identifiers(values=["NR_A.1", "NR_A.1"]) == ("NR_A.1",)
    assert normalise_sequence(value="ac dE*\n") == "ACDE"
    assert column(columns=["f0"], choices=["SEQID", "f0"]) == '"f0"'
    assert column(columns=[], choices=["optional"], optional=True) == "NULL"
    for values in ("abc", [""], ["a b"], ["a\x00"], [None]):
        with pytest.raises(InputValidationError):
            identifiers(values=values)
    for value in (None, "", "A*BC", "A-C", "A1C"):
        with pytest.raises(InputValidationError):
            normalise_sequence(value=value)
    with pytest.raises(InputValidationError):
        column(columns=[], choices=["required"])
