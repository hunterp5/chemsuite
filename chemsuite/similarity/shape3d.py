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
"""3D similarity on the superposition engine: each molecule is fitted onto every reference
(its conformers generated once and tried against them all) and scored by the fit it gets.

The score is the fit's own: a shape Tanimoto, a shape + features combo, or a field similarity.
Only the best reference's score is kept, so this is already fused: a row's "nearest reference"
is the one it fits best.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor

import numpy as np
from rdkit import rdBase

from chemsuite.analysis.fingerprints import Cell, cell_mol
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob
from chemsuite.superpose import model as sp
from chemsuite.superpose.run import best_over_references, prepare_reference

from .model import COLOR, FIELD, SHAPE, Method, Reference

NAN = float("nan")
#: the superpose method each 3D similarity method scores with
METHODS = {SHAPE: sp.SHAPE, COLOR: sp.COLOR, FIELD: sp.FIELD}
#: the fit score each one reports (``COLOR`` adds the two Tanimotos, as Superpose does)
SCORES = {SHAPE: sp.SHAPE_TANIMOTO, COLOR: sp.COMBO, FIELD: sp.FIELD_SIMILARITY}


def settings_for(kind: str, conformers: int, color_weight: float, field_weight: float):
    """The superpose settings a 3D similarity method runs with."""
    return sp.Settings(
        method=METHODS[kind],
        conformers=conformers,
        color_weight=color_weight,
        field_weight=field_weight,
    )


def score_chunk(args: tuple[str, object, Sequence[bytes], Sequence[int], Sequence[Cell]]):
    """Per cell, ``(score, the reference it fits best)`` (a picklable worker entry point).

    *args* is ``(method, superpose settings, the prepared references as blobs, the reference
    each cell is (-1: none), cells)``. ``(NaN, -1)`` for a blank or invalid cell, for a cell
    that will not embed in 3D, and where every reference fails.
    """
    kind, settings, reference_blobs, own, cells = args
    name = SCORES[kind]
    out: list[tuple[float, int]] = []
    with rdBase.BlockLogs():
        references = [mol_from_blob(blob) for blob in reference_blobs]
        for cell, skip in zip(cells, own, strict=True):
            mol = cell_mol(cell)
            if mol is None:
                out.append((NAN, -1))
                continue
            usable = [r for i, r in enumerate(references) if i != skip and r is not None]
            keep = [i for i in range(len(references)) if i != skip and references[i] is not None]
            if not usable:
                out.append((NAN, -1))
                continue
            found = best_over_references(mol, usable, settings)
            if found is None:
                out.append((NAN, -1))
                continue
            index, fit = found
            out.append((float(fit.scores[name]), keep[index]))
    return out


def prepared_references(mols: Sequence[object]) -> list[bytes]:
    """The references as the aligners take them, as blobs for the workers (a reference without
    3D coordinates raises; :func:`~chemsuite.similarity.model.check_references` says so first)."""
    from chemsuite.core.mols import mol_to_blob

    return [mol_to_blob(prepare_reference(mol)) for mol in mols]


#: cells per worker task: a row costs conformers plus one alignment per reference
CHUNK_ROWS = 2


def _chunks(
    kind: str, settings, blobs: list[bytes], own: np.ndarray, cells: Sequence[Cell], size: int
) -> Iterator[tuple]:
    for start in range(0, len(cells), size):
        stop = start + size
        yield kind, settings, blobs, own[start:stop].tolist(), list(cells[start:stop])


def score(
    cells: Sequence[Cell],
    references: Sequence[Reference],
    method: Method,
    own: np.ndarray,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(score, nearest reference)`` per cell by *method* (a 3D method).

    *own* gives, per cell, the reference that is that cell (-1: none), so it is not compared
    with itself. On cancel, raises ``Cancelled`` whose partial is that pair for the cells
    finished so far, or None if none were.
    """
    settings = settings_for(method.kind, method.conformers, method.color_weight,
                            method.field_weight)  # fmt: skip
    blobs = prepared_references([r.mol for r in references])
    rows: list[tuple[float, int]] = []
    done = [0]

    def finished(part: list) -> None:
        done[0] += len(part)
        if progress is not None:
            progress(done[0], len(cells))

    def split(found: list[tuple[float, int]]) -> tuple[np.ndarray, np.ndarray]:
        if not found:
            return np.empty(0), np.empty(0, dtype=np.int64)
        values = np.array([f[0] for f in found], dtype=np.float64)
        which = np.array([f[1] for f in found], dtype=np.int64)
        return values, which

    chunks = _chunks(method.kind, settings, blobs, own, cells, CHUNK_ROWS)
    try:
        for part in ordered_map(score_chunk, chunks, executor, cancel=cancel, finished=finished):
            rows.extend(part)
    except Cancelled:
        raise Cancelled(split(rows) if rows else None) from None
    return split(rows)
