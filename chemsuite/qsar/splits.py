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
"""How the rows are divided, for the outer folds that judge a model and the inner ones that
choose it.

A split decides what a score means. A random split on a congeneric series puts a test
compound's near neighbors in the training set, and every model looks good; a scaffold split
asks the question a project actually faces, which is whether the model carries to the next
series. So **scaffold is the default**, and a grouped split never puts the same scaffold or
cluster on both sides of a fold.

:func:`make_folds` takes the rows to divide and returns folds of indices *into the full set*,
so the inner folds of an outer fold are the same call again over that fold's training rows.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

#: A group of -1 is "no group of its own": such rows are spread over the folds like ungrouped
#: rows, never pooled into one group that would then land whole in one fold.
UNGROUPED = -1


@dataclass(frozen=True)
class Fold:
    """One division: the rows a model is fitted on, and the rows that judge it."""

    train: np.ndarray  # int64 indices into the full set
    test: np.ndarray

    def __post_init__(self) -> None:
        if np.intersect1d(self.train, self.test).size:
            raise ValueError("a fold's training and test rows overlap")


def scaffold_groups(mols: Sequence[Chem.Mol | None]) -> np.ndarray:
    """One group per Bemis-Murcko scaffold, as a text array; a row without a molecule, or one
    with no ring system, gets no group of its own (it is acyclic, so it has no scaffold to
    share) and is marked with the empty string."""
    from chemsuite.descriptors.scaffold import murcko

    out = np.full(len(mols), "", dtype=object)
    for i, mol in enumerate(mols):
        if mol is None:
            continue
        core = murcko(mol)
        if core is not None and core.GetNumAtoms():
            out[i] = Chem.MolToSmiles(core)
    return out.astype(str)


def cluster_groups(groups: np.ndarray) -> np.ndarray:
    """Butina labels (``analysis.clustering.butina``) as split groups: a row that clustered
    keeps its label, and an unclustered one (-1) is ungrouped rather than all one group."""
    labels = np.asarray(groups, dtype=np.int64)
    return np.where(labels < 0, UNGROUPED, labels)


def as_groups(values: Sequence[object]) -> np.ndarray:
    """Any hashable column as integer group ids; "" and None are ungrouped."""
    out = np.full(len(values), UNGROUPED, dtype=np.int64)
    seen: dict[object, int] = {}
    for i, value in enumerate(values):
        if value is None or value == "":
            continue
        if value not in seen:
            seen[value] = len(seen)
        out[i] = seen[value]
    return out


def _deal(sizes: Sequence[int], k: int) -> list[list[int]]:
    """Largest-first: give each item to the fold holding fewest rows so far. The classic
    balanced grouped split — with groups of wildly different size, exact equality is not
    available and this is the closest that keeps groups whole."""
    folds: list[list[int]] = [[] for _ in range(k)]
    load = np.zeros(k, dtype=np.int64)
    for item in np.argsort(-np.asarray(sizes, dtype=np.int64), kind="stable").tolist():
        fold = int(np.argmin(load))
        folds[fold].append(item)
        load[fold] += sizes[item]
    return folds


def _ungrouped_mask(labels: np.ndarray) -> np.ndarray:
    """Which rows have no group of their own: -1 for integer labels, "" for text ones."""
    if np.issubdtype(labels.dtype, np.number):
        return labels == UNGROUPED
    return labels == ""


def _grouped_folds(rows: np.ndarray, groups: np.ndarray, k: int, rng: np.random.Generator):
    """Folds that keep a group whole. Ungrouped rows are dealt out one by one afterwards, so
    they balance the folds instead of forming a group of their own."""
    labels = groups[rows]
    spare = rows[_ungrouped_mask(labels)]
    grouped = rows[~_ungrouped_mask(labels)]
    _, inverse, counts = np.unique(groups[grouped], return_inverse=True, return_counts=True)
    order = rng.permutation(len(counts))  # the seed only breaks ties between equal sizes
    parts: list[list[np.ndarray]] = [[] for _ in range(k)]
    for fold, items in enumerate(_deal(counts[order].tolist(), k)):
        parts[fold].extend(grouped[inverse == order[item]] for item in items)
    for i, row in enumerate(rng.permutation(spare).tolist()):
        parts[i % k].append(np.array([row], dtype=np.int64))
    return [
        np.sort(np.concatenate(part).astype(np.int64)) if part else np.zeros(0, dtype=np.int64)
        for part in parts
    ]


def _stratified_parts(rows: np.ndarray, y: np.ndarray, k: int, rng: np.random.Generator):
    """Equal shares of every class in every fold: deal each class round-robin."""
    parts: list[list[int]] = [[] for _ in range(k)]
    start = 0
    for value in np.unique(y[rows]):
        mine = rng.permutation(rows[y[rows] == value])
        for i, row in enumerate(mine.tolist()):
            parts[(start + i) % k].append(row)
        start = (start + len(mine)) % k  # carry on round-robin, so small classes still spread
    return [np.sort(np.asarray(part, dtype=np.int64)) for part in parts]


def _random_parts(rows: np.ndarray, k: int, rng: np.random.Generator):
    return [np.sort(part) for part in np.array_split(rng.permutation(rows), k)]


def _forward_chaining(rows: np.ndarray, order: np.ndarray, k: int) -> list[Fold]:
    """Train on the past, test on what came next: *k* contiguous blocks in the order of a
    date or series column, each tested against everything before it. No fold ever trains on a
    row that came after the rows it is judged on."""
    ranked = rows[np.argsort(order[rows], kind="stable")]
    blocks = [block for block in np.array_split(ranked, k + 1) if len(block)]
    folds = []
    for i in range(1, len(blocks)):
        train = np.sort(np.concatenate(blocks[:i]))
        folds.append(Fold(train.astype(np.int64), np.sort(blocks[i]).astype(np.int64)))
    return folds


def make_folds(
    kind: str,
    k: int,
    rows: np.ndarray,
    *,
    y: np.ndarray | None = None,
    groups: np.ndarray | None = None,
    order: np.ndarray | None = None,
    seed: int = 0,
    stratify: bool = False,
) -> list[Fold]:
    """*k* folds over *rows*, as indices into the full set.

    *kind* is a key of :data:`~.recipe.SPLITS`. A grouped kind ("scaffold", "cluster") needs
    *groups*, one entry per row of the full set; "column" needs *order*. *stratify* keeps each
    class's share equal in every fold, for a classification on a random split.
    """
    rows = np.sort(np.asarray(rows, dtype=np.int64))
    if k < 2:
        raise ValueError(f"a split needs two folds or more, not {k}")
    if len(rows) < k:
        raise ValueError(f"{len(rows)} rows cannot make {k} folds")
    rng = np.random.default_rng(seed)
    if kind == "column":
        if order is None:
            raise ValueError("the column split needs the column to order by")
        return _forward_chaining(rows, np.asarray(order), k)
    if kind in ("scaffold", "cluster"):
        if groups is None:
            raise ValueError(f"the {kind} split needs a group per row")
        parts = _grouped_folds(rows, np.asarray(groups), k, rng)
    elif kind == "random":
        if stratify and y is not None:
            parts = _stratified_parts(rows, np.asarray(y), k, rng)
        else:
            parts = _random_parts(rows, k, rng)
    else:
        raise ValueError(f"unknown split {kind!r}")
    return [Fold(np.setdiff1d(rows, test), test) for test in parts if len(test)]
