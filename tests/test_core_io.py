"""Unit tests for validation, checksums, text I/O, FASTA and logging."""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import pytest

from protein_signatures.checksums import sha256_file, sha256_json, sha256_text, stable_json
from protein_signatures.errors import ConfigurationError, InputValidationError, PublicationError
from protein_signatures.fasta import (
    _make_record,
    _reject_duplicate,
    read_protein_fasta,
    sequence_kmers,
)
from protein_signatures.io_utils import (
    iter_tsv,
    open_text,
    read_json,
    write_json_atomic,
    write_tsv_atomic,
)
from protein_signatures.logging_config import configure_logging
from protein_signatures.validation import (
    parse_optional_float,
    parse_optional_integer,
    reject_unknown_fields,
    require_mapping,
    require_sequence,
    validate_identifier,
    validate_text,
)


def test_validators_accept_and_normalise_values() -> None:
    """Validators should preserve safe identifiers, text and bounded numbers."""

    assert validate_identifier(value=" A:B/1 ", field_name="id") == "A:B/1"
    assert validate_text(value=" useful text ", field_name="text") == "useful text"
    assert validate_text(value=None, field_name="text", allow_empty=True) == ""
    assert parse_optional_integer(value=" 7 ", field_name="n", minimum=1) == 7
    assert parse_optional_integer(value="", field_name="n") is None
    assert parse_optional_float(value="0.25", field_name="x", minimum=0, maximum=1) == 0.25
    assert parse_optional_float(value=None, field_name="x") is None
    assert require_mapping(value={"x": 1}, field_name="map") == {"x": 1}
    assert require_sequence(value=(1, 2), field_name="list") == [1, 2]
    assert (
        reject_unknown_fields(
            value={"known": 1}, allowed=frozenset({"known"}), field_name="mapping"
        )
        is None
    )
    with pytest.raises(InputValidationError, match="unknown fields"):
        reject_unknown_fields(value={"typo": 1}, allowed=frozenset({"known"}), field_name="mapping")


@pytest.mark.parametrize("value", [None, "", "_bad", "bad space", "x" * 256])
def test_identifier_validator_rejects_unsafe_values(value: object) -> None:
    """Unsafe or overlong identifiers should fail closed."""

    with pytest.raises(InputValidationError):
        validate_identifier(value=value, field_name="identifier")


@pytest.mark.parametrize(
    ("value", "kwargs"),
    [
        ("", {}),
        ("a\tb", {}),
        ("a\nb", {}),
        ("a\x00b", {}),
        ("long", {"maximum_length": 3}),
    ],
)
def test_text_validator_rejects_tsv_breaking_values(value: str, kwargs: dict[str, int]) -> None:
    """Human-readable fields should never corrupt TSV publication."""

    with pytest.raises(InputValidationError):
        validate_text(value=value, field_name="text", **kwargs)


@pytest.mark.parametrize(
    ("function", "value", "kwargs"),
    [
        (parse_optional_integer, "x", {}),
        (parse_optional_integer, "0", {"minimum": 1}),
        (parse_optional_float, "x", {}),
        (parse_optional_float, "nan", {}),
        (parse_optional_float, "-1", {"minimum": 0}),
        (parse_optional_float, "2", {"maximum": 1}),
    ],
)
def test_numeric_validators_reject_invalid_values(
    function: object, value: str, kwargs: dict[str, float]
) -> None:
    """Malformed, non-finite and out-of-range numbers should be rejected."""

    with pytest.raises(InputValidationError):
        function(value=value, field_name="number", **kwargs)  # type: ignore[operator]


def test_container_validators_reject_wrong_types() -> None:
    """Configuration mappings and lists must not accept lookalike scalar values."""

    with pytest.raises(InputValidationError):
        require_mapping(value=[], field_name="map")
    with pytest.raises(InputValidationError):
        require_sequence(value="abc", field_name="list")


def test_checksum_and_stable_json_helpers(tmp_path: Path) -> None:
    """Checksums and canonical JSON should be deterministic and validated."""

    path = tmp_path / "value.txt"
    path.write_text("abc", encoding="utf-8")
    assert sha256_file(path=path, chunk_size=1) == sha256_text(text="abc")
    assert stable_json(value={"b": 2, "a": 1}) == '{"a":1,"b":2}'
    assert sha256_json(value={"a": 1}) == sha256_text(text='{"a":1}')
    with pytest.raises(InputValidationError):
        sha256_file(path=tmp_path / "missing")
    with pytest.raises(InputValidationError):
        sha256_file(path=path, chunk_size=0)


def test_open_text_supports_plain_and_gzip_and_rejects_bad_input(tmp_path: Path) -> None:
    """Text input should support gzip without relaxing validation."""

    plain = tmp_path / "plain.txt"
    plain.write_text("hello\n", encoding="utf-8")
    compressed = tmp_path / "plain.txt.gz"
    with gzip.open(compressed, mode="wt", encoding="utf-8") as handle:
        handle.write("hello\n")
    with open_text(path=plain) as handle:
        assert handle.read() == "hello\n"
    with open_text(path=compressed) as handle:
        assert handle.read() == "hello\n"
    with pytest.raises(InputValidationError):
        with open_text(path=tmp_path / "missing"):
            pass
    bad = tmp_path / "bad.txt"
    bad.write_bytes(b"\xff")
    with pytest.raises(InputValidationError):
        with open_text(path=bad) as handle:
            handle.read()


def test_tsv_reader_and_atomic_writers(tmp_path: Path) -> None:
    """TSV and JSON helpers should publish atomically and retain typed blanks."""

    table = tmp_path / "table.tsv"
    count = write_tsv_atomic(
        path=table,
        fieldnames=("id", "value"),
        records=({"id": "x", "value": 2},),
    )
    assert count == 1
    assert tuple(iter_tsv(path=table, required_fields=("id",))) == ({"id": "x", "value": "2"},)
    document = tmp_path / "value.json"
    write_json_atomic(path=document, value={"b": 2, "a": 1})
    assert read_json(path=document) == {"a": 1, "b": 2}
    with pytest.raises(PublicationError):
        write_tsv_atomic(path=table, fieldnames=(), records=())
    with pytest.raises(PublicationError):
        write_tsv_atomic(path=table, fieldnames=("x",), records=({"y": 1},))
    with pytest.raises(PublicationError):
        write_json_atomic(path=document, value={"bad": object()})
    document.write_text("{", encoding="utf-8")
    with pytest.raises(InputValidationError):
        read_json(path=document)


@pytest.mark.parametrize(
    "content",
    [
        "id\tid\n1\t2\n",
        "id\tvalue\n1\n",
        "id\tvalue\n1\t2\textra\n",
        "other\n1\n",
        "id\n",
    ],
)
def test_tsv_reader_rejects_malformed_tables(tmp_path: Path, content: str) -> None:
    """Duplicate headings, field-count errors, missing columns and emptiness should fail."""

    path = tmp_path / "bad.tsv"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(InputValidationError):
        tuple(iter_tsv(path=path, required_fields=("id",)))
    if content == "id\n":
        assert tuple(iter_tsv(path=path, required_fields=("id",), allow_empty=True)) == ()
    with pytest.raises(InputValidationError):
        tuple(iter_tsv(path=path, required_fields=("id", "id")))


def test_fasta_parsing_and_kmers(tmp_path: Path) -> None:
    """FASTA records, terminal stops and k-mer presence should be deterministic."""

    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">p1 description\nACD\nEF*\n>p2\nAAAA\n", encoding="utf-8")
    records = read_protein_fasta(path=fasta)
    assert [record.protein_id for record in records] == ["p1", "p2"]
    assert records[0].sequence == "ACDEF"
    assert records[0].description == "description"
    assert sequence_kmers(sequence="ABABA", length=2) == frozenset({"AB", "BA"})
    assert sequence_kmers(sequence="AA", length=3) == frozenset()
    with pytest.raises(InputValidationError):
        sequence_kmers(sequence="AA", length=0)
    record = _make_record(header="p3", sequence="ACD")
    assert record.sequence_length == 3
    seen: set[str] = set()
    _reject_duplicate(protein_id="p3", seen=seen)
    with pytest.raises(InputValidationError):
        _reject_duplicate(protein_id="p3", seen=seen)


@pytest.mark.parametrize(
    "content",
    [
        "ACD\n",
        ">\nACD\n",
        ">p1\n",
        ">p1\nAC*D\n",
        ">p1\nAC1\n",
        ">p1\nACD\n>p1\nACE\n",
        "\n",
    ],
)
def test_fasta_rejects_malformed_records(tmp_path: Path, content: str) -> None:
    """Malformed syntax, residues, stops, empty sequences and duplicates should fail."""

    fasta = tmp_path / "bad.faa"
    fasta.write_text(content, encoding="utf-8")
    with pytest.raises(InputValidationError):
        read_protein_fasta(path=fasta)


def test_logging_configuration_validates_level_and_writes_file(tmp_path: Path) -> None:
    """Logging should support explicit levels and persistent UTF-8 output."""

    path = tmp_path / "logs" / "run.log"
    configure_logging(level="debug", log_path=path)
    logging.getLogger("test").warning("visible")
    assert "visible" in path.read_text(encoding="utf-8")
    with pytest.raises(ConfigurationError):
        configure_logging(level="not-a-level")


def test_json_writer_output_is_valid_json(tmp_path: Path) -> None:
    """Atomic JSON should end with a newline and remain machine readable."""

    path = tmp_path / "document.json"
    write_json_atomic(path=path, value=[1, 2])
    assert path.read_text(encoding="utf-8").endswith("\n")
    assert json.loads(path.read_text(encoding="utf-8")) == [1, 2]
