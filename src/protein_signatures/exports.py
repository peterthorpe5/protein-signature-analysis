"""Portable TSV, formatted Excel and PDF exports for reports and the app."""

from __future__ import annotations

import io
import re
from datetime import UTC, datetime
from typing import Any

import pandas as pd

from protein_signatures.errors import InputValidationError, PublicationError

_EXCEL_MAX_COLUMNS = 16_384
_EXCEL_DATA_ROWS_PER_SHEET = 1_048_573
_EXCEL_MAX_CELL_CHARACTERS = 32_767
_EXCEL_LONG_TEXT_CHUNK = 32_000


def normalise_dataframe(*, value: Any) -> pd.DataFrame:
    """Return a detached frame with unique, non-empty text column names.

    Args:
        value: Pandas frame or other two-dimensional tabular value.

    Returns:
        Defensive frame copy suitable for display and export.

    Raises:
        InputValidationError: If columns are blank or duplicated.
    """

    frame = value.copy(deep=True) if isinstance(value, pd.DataFrame) else pd.DataFrame(value)
    names = tuple(str(column).strip() for column in frame.columns)
    if any(not name for name in names) or len(names) != len(set(names)):
        raise InputValidationError("Downloadable tables require unique, non-empty columns.")
    if any(len(name) > 255 for name in names):
        raise InputValidationError(
            "Downloadable table column names must not exceed 255 characters."
        )
    frame.columns = list(names)
    return frame


def dataframe_to_tsv_bytes(*, frame: pd.DataFrame) -> bytes:
    """Serialise a table as deterministic UTF-8 tab-separated text.

    Args:
        frame: Table to serialise without its pandas index.

    Returns:
        UTF-8 TSV bytes with Unix line endings.
    """

    normalised = normalise_dataframe(value=frame)
    buffer = io.StringIO(newline="")
    normalised.to_csv(
        buffer,
        sep="\t",
        index=False,
        lineterminator="\n",
        na_rep="",
    )
    return buffer.getvalue().encode("utf-8")


def dataframe_to_xlsx_bytes(*, frame: pd.DataFrame, title: str) -> bytes:
    """Build a polished, filterable Excel workbook in memory.

    Args:
        frame: Table to export without its pandas index.
        title: Human-readable workbook and worksheet heading.

    Returns:
        Complete XLSX package bytes.

    Raises:
        InputValidationError: If the title or column contract is invalid.
        PublicationError: If the workbook engine cannot produce valid XLSX bytes.
    """

    heading = str(title).strip()
    if not heading:
        raise InputValidationError("An Excel download title must not be empty.")
    if len(heading) > 255:
        raise InputValidationError("An Excel download title must not exceed 255 characters.")
    normalised = normalise_dataframe(value=frame)
    if not len(normalised.columns):
        raise InputValidationError("An Excel download requires at least one column.")
    if len(normalised.columns) > _EXCEL_MAX_COLUMNS:
        raise InputValidationError(
            f"Excel supports at most {_EXCEL_MAX_COLUMNS:,} columns; "
            f"received {len(normalised.columns):,}."
        )
    excel_frame = normalised.copy(deep=True)
    long_text_rows: list[dict[str, Any]] = []
    for column_index, column in enumerate(excel_frame.columns):
        if not (
            pd.api.types.is_object_dtype(excel_frame[column].dtype)
            or pd.api.types.is_string_dtype(excel_frame[column].dtype)
        ):
            continue
        for row_index, value in enumerate(excel_frame[column].tolist()):
            if not isinstance(value, str) or len(value) <= _EXCEL_MAX_CELL_CHARACTERS:
                continue
            parts = tuple(
                value[offset : offset + _EXCEL_LONG_TEXT_CHUNK]
                for offset in range(0, len(value), _EXCEL_LONG_TEXT_CHUNK)
            )
            excel_frame.iat[row_index, column_index] = (
                f"[Full value in Long_Text: data row {row_index + 1}, "
                f"column {column}, {len(parts)} parts]"
            )
            long_text_rows.extend(
                {
                    "data_row_number": row_index + 1,
                    "column_name": str(column),
                    "part_number": part_number,
                    "part_count": len(parts),
                    "text_part": part,
                }
                for part_number, part in enumerate(parts, start=1)
            )
    output = io.BytesIO()
    try:
        with pd.ExcelWriter(
            output,
            engine="xlsxwriter",
            engine_kwargs={
                "options": {
                    "strings_to_formulas": False,
                    "strings_to_urls": False,
                    "nan_inf_to_errors": False,
                }
            },
        ) as writer:
            workbook = writer.book
            workbook.set_properties(
                {
                    "title": heading,
                    "subject": "Protein Signature Analysis table export",
                    "author": "protein-signature-analysis",
                    "company": "Protein Signature Analysis",
                    "comments": "Generated from a checksum-verified result database.",
                    "created": datetime(1980, 1, 1, tzinfo=UTC),
                }
            )
            title_format = workbook.add_format(
                {
                    "bold": True,
                    "font_color": "#FFFFFF",
                    "bg_color": "#17324D",
                    "font_size": 16,
                    "align": "left",
                    "valign": "vcenter",
                }
            )
            subtitle_format = workbook.add_format(
                {"font_color": "#425466", "italic": True, "align": "left"}
            )
            integer_format = workbook.add_format({"num_format": "#,##0", "valign": "top"})
            float_format = workbook.add_format({"num_format": "0.0000", "valign": "top"})
            scientific_format = workbook.add_format({"num_format": "0.00E+00", "valign": "top"})
            text_format = workbook.add_format({"valign": "top"})
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "font_color": "#FFFFFF",
                    "bg_color": "#1F4E78",
                    "border": 1,
                    "border_color": "#D9E2F3",
                    "text_wrap": True,
                    "align": "center",
                    "valign": "vcenter",
                }
            )
            status_formats = tuple(
                (text_value, workbook.add_format({"bg_color": colour}))
                for text_value, colour in (
                    ("FAILED", "#F4CCCC"),
                    ("COMPLETE", "#D9EAD3"),
                    ("VALIDATED", "#D9EAD3"),
                    ("NOT_", "#FFF2CC"),
                    ("INSUFFICIENT", "#FFF2CC"),
                )
            )
            starts = range(0, max(1, len(excel_frame)), _EXCEL_DATA_ROWS_PER_SHEET)
            for sheet_number, start in enumerate(starts, start=1):
                stop = min(len(excel_frame), start + _EXCEL_DATA_ROWS_PER_SHEET)
                chunk = excel_frame.iloc[start:stop]
                sheet_name = "Data" if sheet_number == 1 else f"Data_{sheet_number:03d}"
                worksheet = workbook.add_worksheet(sheet_name)
                writer.sheets[sheet_name] = worksheet
                final_column = len(chunk.columns) - 1
                if final_column:
                    worksheet.merge_range(0, 0, 0, final_column, heading, title_format)
                else:
                    worksheet.write(0, 0, heading, title_format)
                row_range = (
                    f"rows {start + 1:,}-{stop:,} of {len(excel_frame):,}"
                    if len(excel_frame)
                    else "0 rows"
                )
                worksheet.write(
                    1,
                    0,
                    f"{row_range} | {len(chunk.columns):,} columns",
                    subtitle_format,
                )
                chunk.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    startrow=2,
                    index=False,
                    header=True,
                )
                for index, column in enumerate(chunk.columns):
                    worksheet.write(2, index, str(column), header_format)
                    values = chunk[column].dropna().astype(str).head(250).tolist()
                    width = min(
                        48,
                        max(
                            11,
                            len(str(column)) + 2,
                            *(len(value) + 1 for value in values),
                        ),
                    )
                    lower = str(column).casefold()
                    if pd.api.types.is_integer_dtype(chunk[column].dtype):
                        cell_format = integer_format
                    elif pd.api.types.is_float_dtype(chunk[column].dtype):
                        cell_format = (
                            scientific_format
                            if lower in {"p_value", "q_value", "study_q_value", "e_value"}
                            or lower.endswith(("_p_value", "_q_value", "_e_value"))
                            else float_format
                        )
                    else:
                        cell_format = text_format
                    worksheet.set_column(index, index, width, cell_format)
                    if len(chunk) and (
                        "status" in lower or lower in {"evidence_class", "direction"}
                    ):
                        for text_value, status_format in status_formats:
                            worksheet.conditional_format(
                                3,
                                index,
                                2 + len(chunk),
                                index,
                                {
                                    "type": "text",
                                    "criteria": "containing",
                                    "value": text_value,
                                    "format": status_format,
                                },
                            )
                if len(chunk):
                    worksheet.add_table(
                        2,
                        0,
                        2 + len(chunk),
                        final_column,
                        {
                            "name": f"ProteinSignatureData{sheet_number:03d}",
                            "style": "Table Style Medium 2",
                            "columns": [{"header": str(column)} for column in chunk.columns],
                        },
                    )
                else:
                    worksheet.autofilter(2, 0, 2, final_column)
                worksheet.freeze_panes(3, 0)
                worksheet.hide_gridlines(2)
                worksheet.set_row(0, 28)
                worksheet.set_row(2, 34)
                worksheet.set_zoom(90)
            if long_text_rows:
                long_text = pd.DataFrame.from_records(long_text_rows)
                continuation_starts = range(0, len(long_text), _EXCEL_DATA_ROWS_PER_SHEET)
                for continuation_number, start in enumerate(continuation_starts, start=1):
                    stop = min(len(long_text), start + _EXCEL_DATA_ROWS_PER_SHEET)
                    continuation = long_text.iloc[start:stop]
                    sheet_name = (
                        "Long_Text"
                        if continuation_number == 1
                        else f"Long_Text_{continuation_number:03d}"
                    )
                    worksheet = workbook.add_worksheet(sheet_name)
                    writer.sheets[sheet_name] = worksheet
                    worksheet.merge_range(
                        0,
                        0,
                        0,
                        len(continuation.columns) - 1,
                        "Excel long-text continuation values",
                        title_format,
                    )
                    worksheet.write(
                        1,
                        0,
                        (
                            "Concatenate text_part by data_row_number, column_name and "
                            "part_number to recover each complete value."
                        ),
                        subtitle_format,
                    )
                    continuation.to_excel(
                        writer,
                        sheet_name=sheet_name,
                        startrow=2,
                        index=False,
                        header=True,
                    )
                    for index, column in enumerate(continuation.columns):
                        worksheet.write(2, index, str(column), header_format)
                    worksheet.set_column(0, 0, 18, integer_format)
                    worksheet.set_column(1, 1, 28, text_format)
                    worksheet.set_column(2, 3, 14, integer_format)
                    worksheet.set_column(4, 4, 80, text_format)
                    worksheet.add_table(
                        2,
                        0,
                        2 + len(continuation),
                        len(continuation.columns) - 1,
                        {
                            "name": f"ProteinSignatureLongText{continuation_number:03d}",
                            "style": "Table Style Medium 2",
                            "columns": [{"header": str(column)} for column in continuation.columns],
                        },
                    )
                    worksheet.freeze_panes(3, 0)
                    worksheet.hide_gridlines(2)
                    worksheet.set_row(0, 28)
                    worksheet.set_row(2, 34)
                    worksheet.set_zoom(90)
    except (ImportError, OSError, TypeError, ValueError) as error:
        raise PublicationError(f"Could not create formatted Excel download: {error}") from error
    payload = output.getvalue()
    if len(payload) < 100 or not payload.startswith(b"PK"):
        raise PublicationError("The formatted Excel writer returned an invalid XLSX package.")
    return payload


def plotly_figure_to_pdf_bytes(*, figure: Any) -> bytes:
    """Render one Plotly figure to vector PDF through the declared Kaleido backend.

    Args:
        figure: Plotly-compatible figure exposing ``to_image``.

    Returns:
        PDF bytes.

    Raises:
        PublicationError: If Kaleido/Chrome is unavailable or output is invalid.
    """

    try:
        payload = bytes(figure.to_image(format="pdf"))
    except Exception as error:
        raise PublicationError(
            "Could not render this plot as PDF. Check the declared Kaleido installation "
            f"and its Chrome/Chromium runtime: {error}"
        ) from error
    if not payload.startswith(b"%PDF-") or b"%%EOF" not in payload[-1024:]:
        raise PublicationError("The Plotly renderer returned an invalid PDF document.")
    return payload


def safe_download_stem(*, value: str) -> str:
    """Return a bounded filesystem-safe stem for one app download.

    Args:
        value: Human-readable table or plot identity.

    Returns:
        Lower-case underscore-delimited filename stem.

    Raises:
        InputValidationError: If no usable identifier was supplied.
    """

    text = str(value).strip()
    if not text:
        raise InputValidationError("A download name must not be empty.")
    stem = re.sub(r"[^a-z0-9]+", "_", text.casefold()).strip("_")
    return (stem or "export")[:100]
