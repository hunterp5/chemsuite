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
"""Structural alerts: RDKit's FilterCatalog sets matched against many structures, in chunks,
optionally in worker processes.

Each set gives a count column and a names column. A structure is matched as drawn and without
its explicit hydrogens: some sets' patterns count hydrogens and would miss on a structure that
carries them as atoms.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass
from functools import cache

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import as_drawn, mol_from_blob, mol_from_text, without_hydrogens
from chemsuite.core.table.columns import NumberColumn, TextColumn

#: A structure cell for a worker: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
CHUNK_ROWS = 500  # rows per task in a worker process
INLINE_CHUNK_ROWS = 50  # in the job's own thread: a cancel is noticed quickly
SEPARATOR = "; "


@dataclass(frozen=True)
class AlertSet:
    key: str  # RDKit's FilterCatalogs name, and the key saved in the dialog's settings
    title: str  # the column name stem
    help: str


ALERT_SETS: dict[str, AlertSet] = {
    s.key: s
    for s in (
        AlertSet("PAINS", "PAINS", "Pan-assay interference substructures, families A, B and C "
                 "(Baell & Holloway 2010): 480 alerts."),
        AlertSet("BRENK", "Brenk", "Reactive, toxic or metabolically unstable groups unwanted in "
                 "a screening library (Brenk et al. 2008): 105 alerts."),
        AlertSet("NIH", "NIH", "Reactive and assay-interfering groups from the NIH Chemical "
                 "Genomics Center: 180 alerts."),
        AlertSet("ZINC", "ZINC", "Drug-likeness and reactivity filters used to build the ZINC "
                 "database: 50 alerts."),
        AlertSet("CHEMBL_Glaxo", "Glaxo", "Glaxo Wellcome hard filters for unsuitable leads "
                 "(Hann et al. 1999): 55 alerts."),
        AlertSet("CHEMBL_Dundee", "Dundee", "University of Dundee NTD screening-library filters: "
                 "105 alerts."),
        AlertSet("CHEMBL_BMS", "BMS", "Bristol-Myers Squibb HTS deck filters (Pearce et al. "
                 "2006): 180 alerts."),
        AlertSet("CHEMBL_MLSMR", "MLSMR", "NIH Molecular Libraries Small Molecule Repository "
                 "exclusion filters: 116 alerts."),
        AlertSet("CHEMBL_Inpharmatica", "Inpharmatica", "Inpharmatica's unwanted-fragment "
                 "filters: 91 alerts."),
        AlertSet("CHEMBL_SureChEMBL", "SureChEMBL", "SureChEMBL's non-MedChem-friendly filters: "
                 "166 alerts."),
        AlertSet("CHEMBL_LINT", "LINT", "Pfizer's LINT filters: 57 alerts."),
    )
}  # fmt: skip

#: The sets ticked when the dialog is first opened.
COMMON: tuple[str, ...] = ("PAINS", "BRENK", "NIH")

_TRAILING_COUNT = re.compile(r"\(\d+\)$")  # PAINS names carry their hit count: "catechol_A(92)"


def count_column(title: str) -> str:
    return f"{title} Alerts"


def names_column(title: str) -> str:
    return f"{title} Alert Names"


@cache
def catalog(key: str) -> FilterCatalog:
    """RDKit's catalog for *key*, built once per process. Unknown keys raise ``KeyError``."""
    if key not in ALERT_SETS:
        raise KeyError(f"unknown alert set: {key}")
    params = FilterCatalogParams()
    params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, key))
    return FilterCatalog(params)


def alert_name(description: str) -> str:
    """An alert's name as shown: PAINS' trailing hit count dropped, underscores spaced."""
    return _TRAILING_COUNT.sub("", description).strip().replace("_", " ")


def alerts_of(mol: Chem.Mol, key: str) -> list[str]:
    """The names of *key*'s alerts *mol* hits, in the catalog's order, without duplicates."""
    matched = without_hydrogens(as_drawn(mol))
    names: list[str] = []
    for entry in catalog(key).GetMatches(matched):
        name = alert_name(entry.GetDescription())
        if name not in names:
            names.append(name)
    return names


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def match_chunk(args: tuple[tuple[str, ...], Sequence[Cell]]) -> dict[str, list]:
    """``{key: [names per cell]}`` for one chunk (a picklable worker entry point).

    A cell with no structure, or one that cannot be parsed, gives None rather than no alerts.
    """
    keys, cells = args
    out: dict[str, list] = {key: [] for key in keys}
    with rdBase.BlockLogs():
        for cell in cells:
            mol = _mol(cell)
            for key in keys:
                out[key].append(None if mol is None else alerts_of(mol, key))
    return out


def _chunks(keys: tuple[str, ...], cells: Sequence[Cell], size: int) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        yield keys, list(cells[start : start + size])


def compute_alerts(
    cells: Sequence[Cell],
    keys: Sequence[str],
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int | None = None,
) -> dict[str, object]:
    """``{column name: Column}``: a count and a names column per set, cell for cell.

    Unknown keys raise ``KeyError``. If *cancel* is set, raises :class:`Cancelled` whose
    ``partial`` holds the columns for the rows finished so far, or None if none were.
    """
    keys = tuple(keys)
    unknown = [key for key in keys if key not in ALERT_SETS]
    if unknown:  # fail before starting any work
        raise KeyError(f"unknown alert set(s): {', '.join(unknown)}")
    total = len(cells)
    size = chunk_rows or (CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS)
    collected: dict[str, list] = {key: [] for key in keys}
    done = 0
    try:
        for part in ordered_map(match_chunk, _chunks(keys, cells, size), executor, cancel=cancel):
            for key in keys:
                collected[key].extend(part[key])
            done = min(total, done + size)
            if progress is not None:
                progress(done, total)
    except Cancelled:
        raise Cancelled(_columns(keys, collected) if done else None) from None
    return _columns(keys, collected)


def _columns(keys: tuple[str, ...], collected: dict[str, list]) -> dict[str, object]:
    columns: dict[str, object] = {}
    for key in keys:
        title = ALERT_SETS[key].title
        hits = collected[key]
        counts = np.array([np.nan if names is None else len(names) for names in hits], float)
        columns[count_column(title)] = NumberColumn.from_floats(counts, 0)
        columns[names_column(title)] = TextColumn.from_texts(
            ["" if names is None else SEPARATOR.join(names) for names in hits]
        )
    return columns
