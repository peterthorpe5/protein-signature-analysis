"""Tests for starter path authority handling."""

from pathlib import Path

import pytest

from protein_signatures.errors import InputValidationError
from protein_signatures.starter import _profile_source


def test_profile_source_preserves_builtin_and_absolutises_custom_yaml(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Custom profile paths should remain valid when the campaign lives elsewhere."""

    profile = tmp_path / "custom.yaml"
    profile.write_text("profile_id: custom\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    assert _profile_source(value="e3") == "e3"
    assert _profile_source(value="custom.yaml") == str(profile.resolve())


def test_profile_source_rejects_missing_or_empty_custom_yaml(tmp_path: Path) -> None:
    """Path-like profile arguments should fail before a campaign YAML is written."""

    with pytest.raises(InputValidationError, match="missing or empty"):
        _profile_source(value=str(tmp_path / "missing.yaml"))
    empty = tmp_path / "empty.yml"
    empty.touch()
    with pytest.raises(InputValidationError, match="missing or empty"):
        _profile_source(value=str(empty))
