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

"""Cutting molecules into fragments: BRICS and RECAP.

Both rules cut bonds a chemist could make again:

* **BRICS** (Degen et al., *ChemMedChem* 2008, 3, 1503) cuts 16 kinds of bond and labels each
  open end with the environment it came from (``[4*]``, ``[16*]``), so the pieces can only be
  put back together where the rules allow. :mod:`chemsuite.analysis.recombine` needs those
  labels.
* **RECAP** (Lewell et al., *J. Chem. Inf. Comput. Sci.* 1998, 38, 511) cuts 11 retrosynthetic
  bonds and writes every open end as a plain ``*``. Its result is a tree; the leaves are the
  pieces that cannot be cut further.

Each fragment is canonical SMILES. With *plain* the dummy labels are dropped from both
methods, so the same piece counts once however it was cut; the labels are kept by default
because recombination needs them.

Rows are fragmented in chunks, in worker processes when given an executor, and a cancel
raises :class:`Cancelled` carrying the rows already done.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import BRICS, Recap

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text

#: A structure cell as sent to a worker: an RDKit blob, text not parsed yet, or None.
Cell = bytes | str | None
Progress = Callable[[int, int], None]

CHUNK_ROWS = 200
INLINE_CHUNK_ROWS = 20
POOL_MIN_ROWS = 500
#: the methods :func:`fragment` takes
METHODS = ("brics", "recap")
#: fragments smaller than this are not worth keeping (heavy atoms, the dummy not counted)
MIN_SIZE = 3
_LABELLED_DUMMY = re.compile(r"\[\d+\*\]")
_FAILURES = (ValueError, RuntimeError)


def plain_dummies(smiles: str) -> str:
    """*smiles* with every labelled attachment point written as a plain ``*``."""
    return _LABELLED_DUMMY.sub("[*]", smiles)


def _heavy(smiles: str) -> int:
    """Heavy atoms in a fragment, not counting its attachment points."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return 0
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def _brics(mol: Chem.Mol) -> list[str]:
    return sorted(BRICS.BRICSDecompose(mol, returnMols=False, keepNonLeafNodes=False))


def _recap(mol: Chem.Mol) -> list[str]:
    tree = Recap.RecapDecompose(mol)
    leaves = tree.GetLeaves()
    if leaves:
        return sorted(leaves)
    return [Chem.MolToSmiles(mol)]  # nothing could be cut: the molecule is its own fragment


def fragment(
    mol: Chem.Mol, method: str = "brics", *, min_size: int = MIN_SIZE, plain: bool = False
) -> list[str]:
    """*mol* cut into fragments by *method*, largest first.

    Fragments with fewer than *min_size* heavy atoms are left out. With *plain*, attachment
    points lose their BRICS labels.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, not {method!r}")
    with rdBase.BlockLogs():
        try:
            found = _brics(mol) if method == "brics" else _recap(mol)
        except _FAILURES:
            return []
    if plain:
        seen: dict[str, None] = {}
        for smiles in found:
            seen.setdefault(plain_dummies(smiles), None)
        found = list(seen)
    kept = [s for s in found if _heavy(s) >= min_size]
    return sorted(kept, key=lambda s: (-_heavy(s), s))


@dataclass
class FragmentRows:
    """The fragments of each row, in row order. ``[]`` where a row has none."""

    rows: list[list[str]] = field(default_factory=list)
    #: rows whose cell held something that is not a structure
    failed: int = 0

    def __len__(self) -> int:
        return len(self.rows)

    def extend(self, other: FragmentRows) -> None:
        self.rows.extend(other.rows)
        self.failed += other.failed

    def widest(self) -> int:
        """The most fragments any row gave."""
        return max((len(r) for r in self.rows), default=0)

    def fragmented(self) -> int:
        """How many rows gave at least one fragment."""
        return sum(1 for row in self.rows if row)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def fragment_chunk(args: tuple[str, int, bool, Sequence[Cell]]) -> FragmentRows:
    """One chunk of cells (a picklable worker entry point)."""
    method, min_size, plain, cells = args
    out = FragmentRows()
    with rdBase.BlockLogs():
        for cell in cells:
            mol = _mol(cell)
            if mol is None:
                out.rows.append([])
                out.failed += cell is not None
                continue
            out.rows.append(fragment(mol, method, min_size=min_size, plain=plain))
    return out


def _chunked(total: int, executor: Executor | None) -> Iterator[slice]:
    size = CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS
    for start in range(0, total, size):
        yield slice(start, min(start + size, total))


def fragment_rows(
    cells: Sequence[Cell],
    method: str = "brics",
    *,
    min_size: int = MIN_SIZE,
    plain: bool = False,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> FragmentRows:
    """Fragment every cell. Cancelled: ``partial`` holds the rows already finished."""
    cells = list(cells)
    total = len(cells)
    chunks = ((method, min_size, plain, cells[s]) for s in _chunked(total, executor))
    out = FragmentRows()
    try:
        for part in ordered_map(fragment_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(len(out), total)
    except Cancelled:
        raise Cancelled(out if len(out) else None) from None
    return out


@dataclass(frozen=True)
class FragmentStat:
    """One distinct fragment: where it is, and the activity of the rows that hold it."""

    smiles: str
    rows: tuple[int, ...]  # indices into the oids the run covered
    mean: float
    sd: float

    @property
    def count(self) -> int:
        return len(self.rows)


def fragment_stats(found: FragmentRows, activity: np.ndarray | None = None) -> list[FragmentStat]:
    """Each distinct fragment with its compound count and activity, commonest first.

    *activity* is one value per row (NaN where it is missing); without it the mean and SD are
    NaN. Rows are counted once per fragment however often it appears in them.
    """
    where: dict[str, list[int]] = {}
    for index, fragments in enumerate(found.rows):
        for smiles in dict.fromkeys(fragments):  # once per row
            where.setdefault(smiles, []).append(index)
    values = np.asarray(activity, dtype=np.float64) if activity is not None else None
    out: list[FragmentStat] = []
    for smiles, rows in where.items():
        mean = sd = float("nan")
        if values is not None and len(values):
            taken = values[np.asarray(rows, dtype=np.int64)]
            taken = taken[~np.isnan(taken)]
            if len(taken):
                mean = float(taken.mean())
                sd = float(taken.std(ddof=1)) if len(taken) > 1 else 0.0
        out.append(FragmentStat(smiles, tuple(rows), mean, sd))
    out.sort(key=lambda s: (-s.count, s.smiles))
    return out


__all__ = [
    "METHODS",
    "MIN_SIZE",
    "POOL_MIN_ROWS",
    "Cell",
    "FragmentRows",
    "FragmentStat",
    "fragment",
    "fragment_chunk",
    "fragment_rows",
    "fragment_stats",
    "plain_dummies",
]
