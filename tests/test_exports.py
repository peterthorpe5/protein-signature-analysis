"""Unit tests for portable app and pipeline export helpers."""

from __future__ import annotations

import zipfile
from io import BytesIO

import pandas as pd
import pytest

import protein_signatures.exports as exports_module
from protein_signatures.errors import InputValidationError, PublicationError
from protein_signatures.exports import (
    dataframe_to_tsv_bytes,
    dataframe_to_xlsx_bytes,
    normalise_dataframe,
    plotly_figure_to_pdf_bytes,
    safe_download_stem,
)


def test_dataframe_normalisation_and_tsv_are_deterministic() -> None:
    """Frame copies and TSV downloads should preserve typed values and tab delimiters."""

    source = pd.DataFrame({"protein_id": ["p1", "p2"], "score": [1.25, None]})
    normalised = normalise_dataframe(value=source)
    normalised.iloc[0, 0] = "changed"
    assert source.iloc[0, 0] == "p1"
    payload = dataframe_to_tsv_bytes(frame=source)
    assert payload == b"protein_id\tscore\np1\t1.25\np2\t\n"
    with pytest.raises(InputValidationError, match="unique"):
        normalise_dataframe(value=pd.DataFrame([[1, 2]], columns=["same", "same"]))


def test_excel_download_has_formatting_filters_and_safe_text() -> None:
    """Excel downloads should be valid, filterable workbooks with literal formula text."""

    frame = pd.DataFrame(
        {
            "protein_id": ["=unsafe", "p2"],
            "status": ["COMPLETE", "FAILED"],
            "count": [1, 2],
            "p_value": [0.001, 0.5],
        }
    )
    first = dataframe_to_xlsx_bytes(frame=frame, title="Protein evidence")
    second = dataframe_to_xlsx_bytes(frame=frame, title="Protein evidence")
    assert first == second
    with zipfile.ZipFile(BytesIO(first)) as archive:
        names = set(archive.namelist())
        worksheet = archive.read("xl/worksheets/sheet1.xml")
        table = archive.read("xl/tables/table1.xml")
        assert "xl/tables/table1.xml" in names
        assert b"<pane" in worksheet
        assert b"<autoFilter" in table
        assert b"<f>unsafe</f>" not in worksheet
    empty = dataframe_to_xlsx_bytes(
        frame=pd.DataFrame(columns=["protein_id", "status"]),
        title="Empty evidence",
    )
    assert empty.startswith(b"PK")
    with pytest.raises(InputValidationError, match="title"):
        dataframe_to_xlsx_bytes(frame=frame, title="")
    with pytest.raises(InputValidationError, match="at least one column"):
        dataframe_to_xlsx_bytes(frame=pd.DataFrame(), title="Missing columns")


def test_excel_download_splits_rows_across_filterable_sheets(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Large table exports should respect Excel's worksheet row boundary."""

    monkeypatch.setattr(exports_module, "_EXCEL_DATA_ROWS_PER_SHEET", 2)
    payload = dataframe_to_xlsx_bytes(
        frame=pd.DataFrame({"protein_id": ["p1", "p2", "p3"]}),
        title="Split evidence",
    )
    with zipfile.ZipFile(BytesIO(payload)) as archive:
        workbook = archive.read("xl/workbook.xml")
        assert b'Data"' in workbook
        assert b'Data_002"' in workbook
        assert "xl/tables/table1.xml" in archive.namelist()
        assert "xl/tables/table2.xml" in archive.namelist()
    monkeypatch.setattr(exports_module, "_EXCEL_MAX_COLUMNS", 0)
    with pytest.raises(InputValidationError, match="at most"):
        dataframe_to_xlsx_bytes(
            frame=pd.DataFrame({"protein_id": ["p1"]}),
            title="Too wide",
        )


def test_excel_download_preserves_values_longer_than_one_cell() -> None:
    """Long protein text should remain recoverable through continuation sheets."""

    sequence = "ACDE" * 10_000
    payload = dataframe_to_xlsx_bytes(
        frame=pd.DataFrame({"protein_id": ["p1"], "sequence": [sequence]}),
        title="Long sequence",
    )
    with zipfile.ZipFile(BytesIO(payload)) as archive:
        workbook = archive.read("xl/workbook.xml")
        shared = archive.read("xl/sharedStrings.xml")
        assert b"Long_Text" in workbook
        assert shared.count(b"ACDE") == 10_000
        assert b"Full value in Long_Text" in shared
    with pytest.raises(InputValidationError, match="255"):
        dataframe_to_xlsx_bytes(
            frame=pd.DataFrame({"x": [1]}),
            title="x" * 256,
        )
    with pytest.raises(InputValidationError, match="255"):
        normalise_dataframe(value=pd.DataFrame([[1]], columns=["x" * 256]))


def test_plotly_pdf_validation_and_download_names() -> None:
    """PDF validation and portable filename generation should fail closed."""

    class Figure:
        """Minimal Plotly-compatible test double."""

        def __init__(self, payload: bytes | None = None, failure: bool = False) -> None:
            """Store the configured rendering outcome."""

            self.payload = payload
            self.failure = failure

        def to_image(self, *, format: str) -> bytes:
            """Return controlled PDF bytes or fail like an unavailable renderer."""

            assert format == "pdf"
            if self.failure:
                raise RuntimeError("missing browser")
            return self.payload or b""

    document = b"%PDF-1.7\nobject\n%%EOF\n"
    assert plotly_figure_to_pdf_bytes(figure=Figure(document)) == document
    with pytest.raises(PublicationError, match="Kaleido"):
        plotly_figure_to_pdf_bytes(figure=Figure(failure=True))
    with pytest.raises(PublicationError, match="invalid PDF"):
        plotly_figure_to_pdf_bytes(figure=Figure(b"not a pdf"))
    assert safe_download_stem(value="F-box / SHAP plot") == "f_box_shap_plot"
    assert safe_download_stem(value="///") == "export"
    assert len(safe_download_stem(value="x" * 200)) == 100
    with pytest.raises(InputValidationError, match="name"):
        safe_download_stem(value="")
