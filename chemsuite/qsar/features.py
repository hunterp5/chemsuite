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
"""The feature pool: a recipe and some structures in, one matrix out.

Nothing here is new chemistry. The descriptors come from
:mod:`chemsuite.descriptors.compute` and the fingerprints from
:mod:`chemsuite.analysis.fingerprints`, exactly as the Calculate Descriptors and similarity
tools get them, so a feature in a model means the same number as the column in the table.

The pool is deliberately wide — Morgan counts at 2048 bits beside every numeric descriptor is
of order 2,500 columns, on perhaps a few hundred rows. That is safe only because selection
(:mod:`chemsuite.qsar.select`) runs inside every fold; it is the recipe, not restraint, that
keeps a fit honest.

**Column order is part of the recipe**: descriptors in the recipe's order, then each
fingerprint variant's bits in order, then the extra table columns. A model archive stores the
recipe and the names, and rebuilding the pool from that recipe gives the same columns in the
same places — which is what lets a saved model be applied to rows it has never seen.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from chemsuite.analysis import fingerprints as fp
from chemsuite.core.jobs import CancelToken, check
from chemsuite.descriptors.compute import Cell, compute_descriptors

from .recipe import Recipe

if TYPE_CHECKING:  # rdkit is imported where it is used, so the pool costs nothing to import
    from rdkit import Chem

#: Called with (rows done, rows expected) over the whole pool.
Progress = Callable[[int, int], None]

#: What a column of the pool came from, which is also the prefix of its name.
DESCRIPTOR, FINGERPRINT, COLUMN = "d", "f", "c"


@dataclass(frozen=True)
class Pool:
    """The feature matrix, its column names, and which rows could be featurized at all."""

    x: np.ndarray  # rows × columns, float64, NaN where a feature could not be computed
    names: tuple[str, ...]
    kinds: tuple[str, ...]  # DESCRIPTOR, FINGERPRINT or COLUMN, per column
    usable: np.ndarray  # bool per row

    def __post_init__(self) -> None:
        if not (self.x.shape[1] == len(self.names) == len(self.kinds)):
            raise ValueError("the pool's columns, names and kinds disagree")
        if self.x.shape[0] != len(self.usable):
            raise ValueError("the pool's rows and usable flags disagree")

    @property
    def rows(self) -> int:
        return self.x.shape[0]

    @property
    def columns(self) -> int:
        return self.x.shape[1]

    def title_of(self, column: int) -> str:
        """A column's name as a person reads it."""
        return title_of(self.names[column])


def title_of(name: str) -> str:
    """A pool name (``d:MolWt``, ``f:morgan2@c:17``, ``c:Measured LogD``) as a title."""
    kind, _, rest = name.partition(":")
    if kind == DESCRIPTOR:
        from chemsuite.descriptors.catalog import CATALOG

        return CATALOG[rest].title if rest in CATALOG else rest
    if kind == FINGERPRINT:
        key, _, bit = rest.rpartition(":")
        return f"{fp.kind_of(key).column} bit {bit}"
    return rest


def pool_names(recipe: Recipe) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """The names and kinds of every column the recipe asks for, in order — without computing
    anything, so a dialog can say how wide the pool will be."""
    names: list[str] = []
    kinds: list[str] = []
    for key in recipe.descriptors:
        names.append(f"{DESCRIPTOR}:{key}")
        kinds.append(DESCRIPTOR)
    for key in recipe.fingerprints:
        bits = fp.kind_of(key).bits
        names.extend(f"{FINGERPRINT}:{key}:{bit}" for bit in range(bits))
        kinds.extend([FINGERPRINT] * bits)
    for name in recipe.columns:
        names.append(f"{COLUMN}:{name}")
        kinds.append(COLUMN)
    return tuple(names), tuple(kinds)


def width(recipe: Recipe) -> int:
    """How many columns the pool will have."""
    return (
        len(recipe.descriptors)
        + sum(fp.kind_of(key).bits for key in recipe.fingerprints)
        + len(recipe.columns)
    )


class _Share:
    """Spreads one progress count over several blocks of work."""

    def __init__(self, blocks: int, rows: int, progress: Progress | None) -> None:
        self.blocks = max(blocks, 1)
        self.rows = rows
        self.progress = progress
        self.done = 0

    def block(self) -> Progress | None:
        """A progress callback for one block, mapped onto the whole."""
        if self.progress is None:
            return None
        start = self.done

        def report(done: int, total: int) -> None:
            if self.progress is not None:
                fraction = done / total if total else 1.0
                self.progress(int(start + fraction * self.rows), self.rows * self.blocks)

        return report

    def finished(self) -> None:
        self.done += self.rows
        if self.progress is not None:
            self.progress(self.done, self.rows * self.blocks)


def _descriptor_block(
    cells: Sequence[Cell],
    recipe: Recipe,
    executor: Executor | None,
    progress: Progress | None,
    cancel: CancelToken | None,
) -> np.ndarray:
    computed = compute_descriptors(
        cells, recipe.descriptors, executor=executor, progress=progress, cancel=cancel
    )
    return np.column_stack([computed[key].values() for key in recipe.descriptors])


def _fingerprint_block(
    cells: Sequence[Cell],
    key: str,
    executor: Executor | None,
    progress: Progress | None,
    cancel: CancelToken | None,
) -> tuple[np.ndarray, np.ndarray]:
    """One variant's columns, and which rows it could fingerprint."""
    if fp.parse(key).counts:
        counts, valid = fp.compute_counts(
            cells, key, executor=executor, progress=progress, cancel=cancel
        )
        return counts.astype(np.float64), valid
    bits = fp.compute_fingerprints(cells, key, executor=executor, progress=progress, cancel=cancel)
    return bits.unpacked().astype(np.float64), bits.valid


def build(
    cells: Sequence[Cell],
    recipe: Recipe,
    *,
    columns: Mapping[str, np.ndarray] | None = None,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Pool:
    """The pool a recipe asks for, row for row with *cells*.

    *columns* holds the table columns the recipe names, as float arrays with NaN for a blank
    cell — the caller reads them, because this package never sees a table.
    """
    rows = len(cells)
    columns = columns or {}
    missing = [name for name in recipe.columns if name not in columns]
    if missing:
        raise KeyError(f"the recipe needs column(s) the caller did not pass: {missing}")
    blocks = bool(recipe.descriptors) + len(recipe.fingerprints) + bool(recipe.columns)
    share = _Share(blocks, rows, progress)
    parts: list[np.ndarray] = []
    usable = np.ones(rows, dtype=np.bool_)
    parsed: np.ndarray | None = None

    if recipe.descriptors:
        check(cancel)
        block = _descriptor_block(cells, recipe, executor, share.block(), cancel)
        parts.append(block)
        share.finished()
    for key in recipe.fingerprints:
        check(cancel)
        block, valid = _fingerprint_block(cells, key, executor, share.block(), cancel)
        parts.append(block)
        parsed = valid if parsed is None else (parsed & valid)
        share.finished()
    if recipe.columns:
        check(cancel)
        parts.append(
            np.column_stack(
                [
                    np.asarray(columns[name], dtype=np.float64).reshape(rows)
                    for name in recipe.columns
                ]
            )
        )
        share.finished()

    if parsed is not None:
        usable &= parsed
    elif recipe.descriptors:  # no fingerprint said so: a row with no descriptor at all failed
        computed = len(recipe.descriptors)
        usable &= np.isfinite(parts[0][:, :computed]).any(axis=1)
    names, kinds = pool_names(recipe)
    x = np.column_stack(parts) if parts else np.zeros((rows, 0))
    return Pool(np.ascontiguousarray(x, dtype=np.float64), names, kinds, usable)


def mols_for(cells: Sequence[Cell]) -> list[Chem.Mol | None]:
    """A molecule per cell, or None where RDKit cannot read the structure.

    A structure column holds a binary blob or text, never a molecule, so anything that works
    on molecules (the scaffolds a split groups by, the keys a refit matches on) has to read
    them first. RDKit's parse warnings are blocked, so an unreadable row becomes None rather
    than a line of console noise per row.
    """
    from rdkit import rdBase

    from chemsuite.core.mols import mol_from_blob, mol_from_text

    out: list[Chem.Mol | None] = []
    with rdBase.BlockLogs():
        for cell in cells:
            if isinstance(cell, bytes):
                out.append(mol_from_blob(cell) if cell else None)
            elif isinstance(cell, str):
                out.append(mol_from_text(cell))
            else:
                out.append(None)
    return out


def keys_for(cells: Sequence[Cell]) -> list[str]:
    """An InChIKey per cell, or "" where RDKit cannot read the structure.

    This is how a refit tells a measurement that changed from a compound it has never seen:
    the key, not the row number, which moves when rows are sorted or deleted. A tautomer or a
    different salt gives a different key, which is the conservative answer — such a row joins
    the training set as a new compound rather than silently replacing one.
    """
    from rdkit import Chem, rdBase

    with rdBase.BlockLogs():
        return ["" if mol is None else (Chem.MolToInchiKey(mol) or "") for mol in mols_for(cells)]
