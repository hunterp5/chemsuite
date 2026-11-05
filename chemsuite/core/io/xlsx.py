# This file is part of chemsuite.
# Copyright (C) 2026 Hunter Picard
#
# chemsuite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# chemsuite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with chemsuite.  If not, see <https://www.gnu.org/licenses/>.

"""Excel workbooks (.xlsx) via openpyxl, imported only when an Excel file is opened or saved."""

from __future__ import annotations

import datetime as dt
from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.table import NUMBER, Table
from chemsuite.core.table.columns import StructureColumn

from .delimited import AUTO, CHUNK_ROWS
from .formats import Batch, atomic_write, guess_structure_column, unique_headers


def _cell(value: object) -> object:
    if isinstance(value, dt.datetime):
        return value.isoformat(sep=" ", timespec="seconds")
    if isinstance(value, dt.date | dt.time):
        return value.isoformat()
    return value


def read_xlsx(
    path: str | Path,
    *,
    sheet: str | None = None,
    structure_column: str | None = AUTO,
    chunk_rows: int = CHUNK_ROWS,
    cancel: CancelToken | None = None,
) -> Iterator[Batch]:
    """Rows of the first (or named) sheet; the first row is the header."""
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        worksheet = workbook[sheet] if sheet else workbook.worksheets[0]
        rows = worksheet.iter_rows(values_only=True)
        header = next(rows, None)
        if header is None:
            return
        headers = unique_headers(["" if h is None else str(h) for h in header])
        if structure_column == AUTO:
            structure_column = guess_structure_column(headers)
        chunk: list[tuple] = []
        for row in rows:
            if row and any(v is not None and v != "" for v in row):
                chunk.append(row)
            if len(chunk) >= chunk_rows:
                check(cancel)
                yield _batch(headers, chunk, structure_column)
                chunk = []
        if chunk:
            yield _batch(headers, chunk, structure_column)
    finally:
        workbook.close()


def _batch(headers: list[str], rows: list[tuple], structure: str | None) -> Batch:
    batch: dict[str, object] = {}
    for i, name in enumerate(headers):
        values = [_cell(row[i]) if i < len(row) else None for row in rows]
        if name == structure:
            batch[name] = StructureColumn.from_texts(["" if v is None else str(v) for v in values])
        else:
            batch[name] = values
    return batch  # type: ignore[return-value]


def write_xlsx(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | np.ndarray | None = None,
    columns: Sequence[str] | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Numbers are written as numbers (so Excel can compute with them), everything else as text."""
    import openpyxl

    names = list(columns or table.column_names)
    slots = np.arange(len(table)) if oids is None else table.slots_of(oids)
    workbook = openpyxl.Workbook(write_only=True)
    sheet = workbook.create_sheet("chemsuite")
    sheet.append(names)
    for start in range(0, len(slots), 4096):
        check(cancel)
        chunk = slots[start : start + 4096]
        cells = []
        for name in names:
            column = table.column(name)
            if column.kind == NUMBER:
                values = column.values(chunk)
                cells.append([None if np.isnan(v) else float(v) for v in values])
            else:
                cells.append([text or None for text in column.displays(chunk)])
        for row in zip(*cells, strict=True):
            sheet.append(list(row))
    with atomic_write(path, "wb") as handle:
        workbook.save(handle)
