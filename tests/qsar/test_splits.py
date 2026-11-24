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
"""The folds. The test that matters is that a grouped split keeps a scaffold on one side: a
fold that leaks a scaffold into its own training set is what makes a QSAR score a lie.
"""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.qsar.splits import (
    UNGROUPED,
    Fold,
    as_groups,
    cluster_groups,
    make_folds,
    scaffold_groups,
)

SERIES = [
    "c1ccccc1C(=O)NC",
    "c1ccccc1C(=O)NCC",
    "c1ccccc1C(=O)NCCC",
    "c1ccncc1C(=O)NC",
    "c1ccncc1C(=O)NCC",
    "C1CCCCC1C(=O)NC",
    "C1CCCCC1C(=O)NCC",
    "CCCCCCC",
    "CCCCCCCC",
]


def mols():
    return [Chem.MolFromSmiles(smiles) for smiles in SERIES]


def rows(n: int) -> np.ndarray:
    return np.arange(n, dtype=np.int64)


# -- the point of the module ----------------------------------------------------------------
def test_a_grouped_fold_never_sees_its_own_scaffold():
    groups = scaffold_groups(mols())
    for fold in make_folds("scaffold", 3, rows(len(SERIES)), groups=groups, seed=0):
        shared = set(groups[fold.train]) & set(groups[fold.test])
        assert not (shared - {""}), f"a scaffold is on both sides: {shared}"


def test_a_random_split_does_leak_the_scaffold():
    """Not a wish — the contrast is why scaffold is the default."""
    groups = scaffold_groups(mols())
    leaks = [
        bool((set(groups[f.train]) & set(groups[f.test])) - {""})
        for f in make_folds("random", 3, rows(len(SERIES)), seed=1)
    ]
    assert any(leaks)


def test_clusters_stay_whole():
    labels = np.array([0, 0, 1, 1, 1, 2, -1, -1, 3, 3])
    groups = cluster_groups(labels)
    for fold in make_folds("cluster", 3, rows(len(labels)), groups=groups, seed=0):
        shared = set(groups[fold.train]) & set(groups[fold.test])
        assert not (shared - {UNGROUPED})


# -- the mechanics --------------------------------------------------------------------------
@pytest.mark.parametrize("kind", ["random", "scaffold", "cluster"])
def test_every_row_is_tested_exactly_once(kind):
    n = len(SERIES)
    groups = scaffold_groups(mols()) if kind == "scaffold" else as_groups(list(range(n)))
    folds = make_folds(kind, 3, rows(n), groups=groups, seed=0)
    tested = np.concatenate([fold.test for fold in folds])
    assert sorted(tested.tolist()) == list(range(n))


@pytest.mark.parametrize("kind", ["random", "scaffold", "cluster"])
def test_train_and_test_together_are_the_whole_set(kind):
    n = len(SERIES)
    groups = scaffold_groups(mols()) if kind == "scaffold" else as_groups(list(range(n)))
    for fold in make_folds(kind, 3, rows(n), groups=groups, seed=0):
        assert sorted(fold.train.tolist() + fold.test.tolist()) == list(range(n))


def test_inner_folds_stay_inside_their_outer_fold():
    groups = scaffold_groups(mols())
    outer = make_folds("scaffold", 3, rows(len(SERIES)), groups=groups, seed=0)[0]
    for inner in make_folds("scaffold", 2, outer.train, groups=groups, seed=0):
        assert set(inner.train) | set(inner.test) <= set(outer.train.tolist())
        assert not set(inner.test) & set(outer.test.tolist())


def test_the_same_seed_gives_the_same_folds():
    groups = scaffold_groups(mols())
    first = make_folds("scaffold", 3, rows(len(SERIES)), groups=groups, seed=4)
    again = make_folds("scaffold", 3, rows(len(SERIES)), groups=groups, seed=4)
    assert [f.test.tolist() for f in first] == [f.test.tolist() for f in again]


def test_stratified_folds_keep_each_class():
    y = np.array([0] * 20 + [1] * 6)
    folds = make_folds("random", 3, rows(len(y)), y=y, seed=0, stratify=True)
    for fold in folds:
        assert set(y[fold.test].tolist()) == {0, 1}
        assert set(y[fold.train].tolist()) == {0, 1}


def test_a_column_split_never_trains_on_the_future():
    order = np.arange(20)[::-1]  # the column need not be sorted
    folds = make_folds("column", 4, rows(20), order=order)
    assert folds
    for fold in folds:
        assert order[fold.train].max() < order[fold.test].min()
    assert len(folds[0].train) < len(folds[-1].train), "training grows as time passes"


def test_an_acyclic_row_has_no_scaffold_to_share():
    groups = scaffold_groups([Chem.MolFromSmiles("CCCCC"), None])
    assert groups.tolist() == ["", ""]


def test_ungrouped_rows_are_spread_not_pooled():
    groups = np.full(12, UNGROUPED)
    folds = make_folds("cluster", 3, rows(12), groups=groups, seed=0)
    assert len(folds) == 3
    assert all(len(fold.test) == 4 for fold in folds)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"kind": "scaffold", "k": 3}, "needs a group"),
        ({"kind": "column", "k": 3}, "column to order by"),
        ({"kind": "tea-leaves", "k": 3}, "unknown split"),
        ({"kind": "random", "k": 1}, "two folds"),
    ],
)
def test_nonsense_is_refused(kwargs, message):
    with pytest.raises(ValueError, match=message):
        make_folds(rows=rows(9), **kwargs)


def test_too_few_rows_for_the_folds_is_refused():
    with pytest.raises(ValueError, match="cannot make"):
        make_folds("random", 5, rows(4))


def test_a_fold_cannot_overlap_itself():
    with pytest.raises(ValueError, match="overlap"):
        Fold(np.array([0, 1, 2]), np.array([2, 3]))
