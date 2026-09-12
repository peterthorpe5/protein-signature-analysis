"""Failure-recovery and checksum tests for immutable result publication."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import protein_signatures.publication as publication_module
from protein_signatures.checksums import sha256_file
from protein_signatures.errors import PublicationError
from protein_signatures.publication import (
    _build_input_manifest,
    _copy_assets,
    _create_duckdb,
    _manifest_files,
    _write_table,
    publish_result,
    verify_completed_result,
    verify_input_authorities,
)


def test_completed_result_rejects_marker_and_manifest_corruption(
    completed_result: Path,
) -> None:
    """Every marker and output-manifest invariant should fail closed."""

    marker_path = completed_result / "COMPLETED.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    marker_path.unlink()
    with pytest.raises(PublicationError, match="lacks completion"):
        verify_completed_result(result_dir=completed_result)
    marker_path.write_text(json.dumps({"status": "FAILED"}), encoding="utf-8")
    with pytest.raises(PublicationError, match="Invalid completion"):
        verify_completed_result(result_dir=completed_result)
    marker_path.write_text(json.dumps(marker), encoding="utf-8")
    marker["manifest_sha256"] = "0" * 64
    marker_path.write_text(json.dumps(marker), encoding="utf-8")
    with pytest.raises(PublicationError, match="Manifest checksum"):
        verify_completed_result(result_dir=completed_result)
    _replace_manifest(result_dir=completed_result, value=[])
    with pytest.raises(PublicationError, match="Invalid output manifest"):
        verify_completed_result(result_dir=completed_result)
    _replace_manifest(result_dir=completed_result, value={"outputs": ["bad"]})
    with pytest.raises(PublicationError, match="Invalid output record"):
        verify_completed_result(result_dir=completed_result)
    _replace_manifest(
        result_dir=completed_result,
        value={"outputs": [{"relative_path": "../escape", "size_bytes": 1, "sha256": "0" * 64}]},
    )
    with pytest.raises(PublicationError, match="missing or unsafe"):
        verify_completed_result(result_dir=completed_result)


def test_completed_result_rejects_same_size_checksum_change(completed_result: Path) -> None:
    """A byte change that preserves file size must still invalidate the result."""

    target = completed_result / "run_metadata.json"
    content = target.read_bytes()
    replacement = (b"X" if content[:1] != b"X" else b"Y") + content[1:]
    target.write_bytes(replacement)
    with pytest.raises(PublicationError, match="checksum mismatch"):
        verify_completed_result(result_dir=completed_result)


def test_completed_result_rejects_undeclared_and_duplicate_paths(
    completed_result: Path,
) -> None:
    """An immutable result must contain exactly one copy of every declared file."""

    extra = completed_result / "stale.tsv"
    extra.write_text("undeclared\n", encoding="utf-8")
    with pytest.raises(PublicationError, match="differs from the manifest"):
        verify_completed_result(result_dir=completed_result)
    extra.unlink()

    manifest = json.loads((completed_result / "manifest.json").read_text(encoding="utf-8"))
    manifest["outputs"].append(dict(manifest["outputs"][0]))
    _replace_manifest(result_dir=completed_result, value=manifest)
    with pytest.raises(PublicationError, match="Duplicate output path"):
        verify_completed_result(result_dir=completed_result)


def test_input_authorities_reject_structure_missing_size_and_checksum(
    completed_result: Path,
) -> None:
    """Computational resume should verify the original local input authorities."""

    manifest_path = completed_result / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    original = manifest["inputs"][0]
    manifest["inputs"] = ["bad"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PublicationError, match="Invalid input authority"):
        verify_input_authorities(result_dir=completed_result)
    manifest["inputs"] = [{**original, "path": str(completed_result / "missing")}]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PublicationError, match="no longer available"):
        verify_input_authorities(result_dir=completed_result)
    source = Path(original["path"])
    manifest["inputs"] = [{**original, "size_bytes": source.stat().st_size + 1}]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PublicationError, match="size changed"):
        verify_input_authorities(result_dir=completed_result)
    manifest["inputs"] = [{**original, "sha256": "0" * 64}]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PublicationError, match="checksum changed"):
        verify_input_authorities(result_dir=completed_result)
    manifest_path.write_text("[]", encoding="utf-8")
    with pytest.raises(PublicationError, match="Invalid input manifest"):
        verify_input_authorities(result_dir=completed_result)


def test_input_manifest_fails_closed_for_missing_or_changing_authority(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Publication must never omit an authority or accept one changed during hashing."""

    missing = tmp_path / "missing.tsv"
    with pytest.raises(PublicationError, match="missing or not a file"):
        _build_input_manifest(input_paths=(missing,))
    with pytest.raises(PublicationError, match="missing or not a file"):
        _build_input_manifest(input_paths=(tmp_path,))
    destination = tmp_path / "result"
    with pytest.raises(PublicationError, match="missing or not a file"):
        publish_result(
            output_dir=destination,
            tables={},
            metadata={},
            input_paths=(missing,),
            asset_sources={},
            resume=False,
        )
    assert not destination.exists()
    assert not tuple(tmp_path.glob(".result.staging.*"))

    source = tmp_path / "input.tsv"
    source.write_text("original\n", encoding="utf-8")

    def mutate_while_hashing(*, path: Path) -> str:
        """Change an authority after returning its original-content digest."""

        digest = sha256_file(path=path)
        path.write_text("changed content\n", encoding="utf-8")
        return digest

    monkeypatch.setattr(publication_module, "sha256_file", mutate_while_hashing)
    with pytest.raises(PublicationError, match="changed while checksumming"):
        _build_input_manifest(input_paths=(source,))


def test_table_publication_validates_arrow_schema_and_parquet_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Typed publication should reject incompatible rows and remove failed temporaries."""

    table_dir = tmp_path / "tables"
    table_dir.mkdir()
    with pytest.raises(PublicationError, match="violate schema"):
        _write_table(
            table_name="proteins",
            records=(
                {
                    "protein_id": "p1",
                    "description": "",
                    "sequence": "AAAA",
                    "sequence_length": "not-an-integer",
                    "sequence_sha256": "a" * 64,
                },
            ),
            table_dir=table_dir,
        )

    def fail_parquet(*_args: object, **_kwargs: object) -> None:
        """Simulate a filesystem failure during Parquet publication."""

        raise OSError("simulated")

    monkeypatch.setattr(publication_module.pq, "write_table", fail_parquet)
    with pytest.raises(PublicationError, match="Could not publish Parquet"):
        _write_table(
            table_name="proteins",
            records=(
                {
                    "protein_id": "p1",
                    "description": "",
                    "sequence": "AAAA",
                    "sequence_length": 4,
                    "sequence_sha256": "a" * 64,
                },
            ),
            table_dir=table_dir,
        )
    assert not tuple(table_dir.glob("*.tmp"))


def test_asset_copy_and_database_failure_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Portable asset paths and DuckDB creation should wrap unsafe or failed operations."""

    staging = tmp_path / "staging"
    staging.mkdir()
    source = tmp_path / "model.pdb"
    source.write_text("ATOM\n", encoding="utf-8")
    _copy_assets(staging=staging, asset_sources={"assets/model.pdb": source})
    assert (staging / "assets" / "model.pdb").read_text(encoding="utf-8") == "ATOM\n"
    for relative in ("/absolute", "../escape"):
        with pytest.raises(PublicationError, match="Unsafe"):
            _copy_assets(staging=staging, asset_sources={relative: source})
    with pytest.raises(PublicationError, match="Missing or empty"):
        _copy_assets(staging=staging, asset_sources={"assets/missing": tmp_path / "missing"})

    def fail_copy(*_args: object, **_kwargs: object) -> None:
        """Simulate a failed portable asset copy."""

        raise OSError("simulated")

    monkeypatch.setattr(publication_module.shutil, "copyfile", fail_copy)
    with pytest.raises(PublicationError, match="Could not copy"):
        _copy_assets(staging=staging, asset_sources={"assets/other.pdb": source})
    with pytest.raises(PublicationError, match="Could not create DuckDB"):
        _create_duckdb(path=tmp_path / "database.duckdb", table_dir=tmp_path / "missing_tables")


def test_top_level_publication_wraps_unexpected_failure_and_cleans_staging(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unexpected stage failure should be contextualised and leave no partial result."""

    source = tmp_path / "input.txt"
    source.write_text("input", encoding="utf-8")
    destination = tmp_path / "result"

    def fail_assets(**_kwargs: object) -> None:
        """Simulate an unexpected staging error."""

        raise RuntimeError("simulated")

    monkeypatch.setattr(publication_module, "_copy_assets", fail_assets)
    with pytest.raises(PublicationError, match="Could not publish result"):
        publish_result(
            output_dir=destination,
            tables={},
            metadata={},
            input_paths=(source,),
            asset_sources={},
            resume=False,
        )
    assert not destination.exists()
    assert not tuple(tmp_path.glob(".result.staging.*"))
    destination.mkdir()
    with pytest.raises(PublicationError, match="will not be overwritten"):
        publish_result(
            output_dir=destination,
            tables={},
            metadata={},
            input_paths=(),
            asset_sources={},
            resume=False,
        )


def test_manifest_file_inventory_is_sorted_and_checksummed(tmp_path: Path) -> None:
    """Manifest enumeration should include only files in deterministic path order."""

    (tmp_path / "b").mkdir()
    (tmp_path / "b" / "two.txt").write_text("two", encoding="utf-8")
    (tmp_path / "one.txt").write_text("one", encoding="utf-8")
    rows = _manifest_files(root=tmp_path)
    assert [row["relative_path"] for row in rows] == ["b/two.txt", "one.txt"]
    assert rows[1]["sha256"] == sha256_file(path=tmp_path / "one.txt")


def _replace_manifest(*, result_dir: Path, value: object) -> None:
    """Replace a manifest and refresh its marker checksum for one failure test."""

    manifest_path = result_dir / "manifest.json"
    manifest_path.write_text(json.dumps(value), encoding="utf-8")
    marker_path = result_dir / "COMPLETED.json"
    marker_path.write_text(
        json.dumps(
            {
                "status": "COMPLETE",
                "manifest_sha256": sha256_file(path=manifest_path),
            }
        ),
        encoding="utf-8",
    )
