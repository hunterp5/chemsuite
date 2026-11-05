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
"""Computing descriptors for many structures, in chunks, optionally in worker processes."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor

import numpy as np
from rdkit import Chem, rdBase

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text
from chemsuite.core.table.columns import Column, NumberColumn, TextColumn

from .catalog import CATALOG, SHAPE_GROUP

#: A structure cell for a worker: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
#: rows per task in a worker process; in the job's own thread chunks are smaller, so a cancel
#: is noticed within a fraction of a second
CHUNK_ROWS = 500
INLINE_CHUNK_ROWS = 50

# Expected per-molecule failures (odd valences, InChI errors…); the cell is left blank.
_FAILURES = (ValueError, RuntimeError, ArithmeticError, KeyError, IndexError)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def _as_drawn(mol: Chem.Mol) -> Chem.Mol:
    """A 3D structure carries explicit hydrogens; 2D descriptors (H counts, SMILES…) see it
    without them, as drawn, so they match the same compound's 2D cell."""
    if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
        return mol
    try:
        return Chem.RemoveHs(mol)
    except _FAILURES:
        return mol


def compute_chunk(args: tuple[tuple[str, ...], Sequence[Cell]]) -> dict[str, list]:
    """``{key: values}`` for one chunk (a picklable worker entry point).

    Numbers are floats (NaN for a failure or a missing structure); text is ``""`` instead.
    """
    keys, cells = args
    specs = [CATALOG[key] for key in keys]
    out: dict[str, list] = {key: [] for key in keys}
    with rdBase.BlockLogs():
        for cell in cells:
            mol = _mol(cell)
            drawn = _as_drawn(mol) if mol is not None else None
            for spec in specs:
                value: object
                if mol is None:
                    value = "" if spec.is_text else np.nan
                else:
                    try:
                        value = spec.fn(mol if spec.group == SHAPE_GROUP else drawn)
                    except _FAILURES:
                        value = "" if spec.is_text else np.nan
                    if not spec.is_text:
                        value = float(value)
                    elif value is None:
                        value = ""
                out[spec.key].append(value)
    return out


def _chunks(keys: tuple[str, ...], cells: Sequence[Cell], size: int) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        yield keys, list(cells[start : start + size])


def compute_descriptors(
    cells: Sequence[Cell],
    keys: Sequence[str],
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
    chunk_rows: int | None = None,
) -> dict[str, Column]:
    """One result column per key, cell for cell. Unknown keys raise ``KeyError``.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` holds the columns for the
    rows finished so far (the first rows, in order), or None if none were.
    """
    keys = tuple(keys)
    unknown = [key for key in keys if key not in CATALOG]
    if unknown:  # fail before starting any work
        raise KeyError(f"unknown descriptor(s): {', '.join(unknown)}")
    total = len(cells)
    size = chunk_rows or (CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS)
    collected: dict[str, list] = {key: [] for key in keys}
    done = 0
    try:
        chunks = _chunks(keys, cells, size)
        for part in ordered_map(compute_chunk, chunks, executor, cancel=cancel):
            for key in keys:
                collected[key].extend(part[key])
            done = min(total, done + size)
            if progress is not None:
                progress(done, total)
    except Cancelled:
        raise Cancelled(_columns(keys, collected) if done else None) from None
    return _columns(keys, collected)


def _columns(keys: tuple[str, ...], collected: dict[str, list]) -> dict[str, Column]:
    columns: dict[str, Column] = {}
    for key in keys:
        spec = CATALOG[key]
        values = collected[key]
        if spec.is_text:
            columns[key] = TextColumn.from_texts([str(v) for v in values])
        else:
            numbers = np.array(values, dtype=np.float64)
            # A tiny negative value would display as "-0.00": show (and store) it as 0.
            with np.errstate(invalid="ignore"):
                numbers[(numbers < 0) & (np.round(numbers, spec.decimals) == 0)] = 0.0
            columns[key] = NumberColumn.from_floats(numbers, spec.decimals)
    return columns
