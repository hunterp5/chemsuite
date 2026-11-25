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

"""CSV/TSV and SMILES-list files: text-first reading and streaming writing.

Structures in a tabular file are loaded as text and parsed afterwards (in parallel, see
:func:`chemsuite.io.load.parse_pending_structures`), so a 200k-row CSV opens as fast as the
text can be read.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np

from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.table import Table
from chemsuite.core.table.columns import ReactionColumn, StructureColumn

from .formats import (
    TSV,
    Batch,
    atomic_write,
    file_size,
    guess_structure_column,
    open_text,
    reaction_columns,
    split_suffix,
    unique_headers,
)

AUTO = "auto"
CHUNK_ROWS = 20_000
_DELIMITERS = (",", "\t", ";", "|")
#: The delimiters an export can offer, ``(label, character)``.
EXPORT_DELIMITERS = (
    ("Tab", "\t"),
    ("Comma", ","),
    ("Semicolon", ";"),
    ("Pipe", "|"),
    ("Space", " "),
)
#: Characters that are never part of a SMILES, so in a ``.smi`` line's first word they can
#: only separate fields (a CXSMILES ``|`` comes after a space).
_SMILES_SEPARATORS = (",", ";", "|")


def sniff_delimiter(header_line: str, default: str) -> str:
    """The delimiter that splits the header into the most fields (ties: comma, tab, ; then |).

    ChEMBL exports use ``;``; most other tools use comma or tab.
    """
    best, best_count = default, 1
    for delimiter in _DELIMITERS:
        count = len(next(csv.reader([header_line], delimiter=delimiter), []))
        if count > best_count:
            best, best_count = delimiter, count
    return best


def read_delimited(
    path: str | Path,
    *,
    delimiter: str | None = None,
    structure_column: str | None = AUTO,
    chunk_rows: int = CHUNK_ROWS,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Iterator[Batch]:
    """Yield batches from a delimited text file whose first row is the header.

    *structure_column* names the column holding SMILES/InChI; ``"auto"`` picks it from the
    header (``SMILES``, ``Structure``, ``canonical_smiles``…); None loads every column as data.
    Columns whose first cells are all reaction SMILES become reaction columns (with ``"auto"``,
    a ``rxn_smiles`` column is a reaction column, not the structure column).
    Short rows are padded with blanks; extra fields get ``Column N`` headers.
    """
    fmt, _ = split_suffix(path)
    total = file_size(path)
    with open_text(path) as (text, raw):
        first = text.readline()
        if not first:
            return
        delimiter = delimiter or sniff_delimiter(first.rstrip("\r\n"), "\t" if fmt == TSV else ",")
        header_row = next(csv.reader([first], delimiter=delimiter), [])
        headers = unique_headers(header_row)
        auto = structure_column == AUTO
        if auto:
            structure_column = guess_structure_column(headers)
        elif structure_column is not None and structure_column not in headers:
            raise ValueError(f"no column named {structure_column!r} in {Path(path).name}")
        reactions: set[str] | None = None

        def batch(rows: list[list[str]]) -> Batch:
            nonlocal reactions, structure_column
            if reactions is None:  # decided once, from the first rows
                reactions = reaction_columns(headers, rows)
                if not auto:
                    reactions.discard(structure_column)  # named as structures on purpose
                if auto and structure_column in reactions:
                    structure_column = guess_structure_column(
                        [h for h in headers if h not in reactions]
                    )
            return _batch(headers, rows, structure_column, reactions)

        rows: list[list[str]] = []
        for row in csv.reader(text, delimiter=delimiter):
            if not row or (len(row) == 1 and not row[0].strip()):
                continue  # blank line
            if len(row) > len(headers):
                headers = unique_headers(headers + [""] * (len(row) - len(headers)))[: len(row)]
            rows.append(row)
            if len(rows) >= chunk_rows:
                check(cancel)
                yield batch(rows)
                rows = []
                if progress is not None:
                    progress(raw.tell(), total)
        if rows:
            yield batch(rows)
        if progress is not None:
            progress(total, total)


def _batch(
    headers: list[str], rows: list[list[str]], structure: str | None, reactions: set[str]
) -> Batch:
    width = len(headers)
    columns: list[list[str]] = [[] for _ in range(width)]
    for row in rows:
        for i in range(width):
            columns[i].append(row[i] if i < len(row) else "")
    batch: dict[str, object] = {}
    for name, values in zip(headers, columns, strict=True):
        if name == structure:
            batch[name] = StructureColumn.from_texts(values)
        elif name in reactions:
            batch[name] = ReactionColumn.from_texts(values)
        else:
            batch[name] = values
    return batch  # type: ignore[return-value]


def read_smiles(
    path: str | Path,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int = CHUNK_ROWS,
) -> Iterator[Batch]:
    """A ``.smi`` file: one structure per line, optionally followed by a name, separated by
    whitespace or, when the first line's first word holds one, by a comma, ``;`` or ``|``.
    A first line starting ``SMILES`` is a header."""
    total = file_size(path)
    with open_text(path) as (text, raw):
        structures: list[str] = []
        names: list[str] = []
        first = True
        separator: str | None = None
        for line in text:
            if first:
                word = next(iter(line.split(None, 1)), "")
                found = [(word.index(c), c) for c in _SMILES_SEPARATORS if c in word]
                separator = min(found)[1] if found else None  # the one right after the SMILES
            parts = _smiles_fields(line, separator)
            if not parts:
                continue
            if first and parts[0].lower() in ("smiles", "smi"):
                first = False
                continue  # header line
            first = False
            structures.append(parts[0])
            names.append(parts[1] if len(parts) > 1 else "")
            if len(structures) >= chunk_rows:
                check(cancel)
                yield {"Structure": StructureColumn.from_texts(structures), "Name": names}
                structures, names = [], []
                if progress is not None:
                    progress(raw.tell(), total)
        if structures:
            yield {"Structure": StructureColumn.from_texts(structures), "Name": names}


def _smiles_fields(line: str, separator: str | None) -> list[str]:
    """``[structure, name]``, ``[structure]``, or ``[]`` for a blank line."""
    if separator is None:
        return line.strip().split(None, 1)
    fields = next(csv.reader([line.strip()], delimiter=separator), [])
    if not fields or not fields[0].strip():
        return []
    return [fields[0].strip(), separator.join(fields[1:]).strip()]


def write_delimited(
    table: Table,
    path: str | Path,
    *,
    delimiter: str = ",",
    oids: Sequence[int] | np.ndarray | None = None,
    columns: Sequence[str] | None = None,
    header: bool = True,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Write the cells exactly as displayed; structure columns are written as their SMILES.
    *header* writes the column names as the first line."""
    names = list(columns or table.column_names)
    slots = np.arange(len(table)) if oids is None else table.slots_of(oids)
    with atomic_write(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        if header:
            writer.writerow(names)
        for start in range(0, len(slots), 4096):
            check(cancel)
            chunk = slots[start : start + 4096]
            cells = [table.column(name).displays(chunk) for name in names]
            writer.writerows(zip(*cells, strict=True))
            if progress is not None:
                progress(min(start + 4096, len(slots)), len(slots))


def write_smiles(
    table: Table,
    path: str | Path,
    *,
    structure: str,
    name_column: str | None = None,
    oids: Sequence[int] | np.ndarray | None = None,
    delimiter: str = "\t",
    header: bool = False,
) -> None:
    """One SMILES per line, then *name_column*'s value after *delimiter*. *header* writes
    ``SMILES`` and the name column's name first. After a comma, ``;`` or ``|`` a name is
    quoted where it needs to be; after a tab or space the rest of the line is the name."""
    column = table.column(structure)
    names = table.column(name_column) if name_column else None
    slots = np.arange(len(table)) if oids is None else table.slots_of(oids)
    quoted = delimiter not in ("\t", " ")
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")

        def write(smiles: str, name: str) -> None:
            if not name:
                handle.write(f"{smiles}\n")
            elif quoted:
                writer.writerow([smiles, name])
            else:
                handle.write(f"{smiles}{delimiter}{name}\n")

        if header:
            write("SMILES", name_column or "")
        for slot in slots.tolist():
            smiles = column.display(slot)
            if "\n" in smiles:  # an unparsed molfile cannot be written to a SMILES list
                smiles = ""
            name = names.display(slot).replace("\n", " ") if names is not None else ""
            write(smiles, name)
