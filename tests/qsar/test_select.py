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
"""The selection cascade: it finds the columns that carry the signal out of a pool of noise,
and it refuses to run on rows whose holes have not been filled from this same fold.
"""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.qsar.prepare import Prepare
from chemsuite.qsar.select import (
    DEFAULT_FILTER,
    FILTERS,
    Selection,
    decorrelated,
    embedded,
    flat_columns,
    select,
    univariate_scores,
)

#: Columns 0, 1 and 2 carry the signal; the rest is noise.
SIGNAL = (0, 1, 2)


def a_pool(rows: int = 200, noise: int = 60, seed: int = 0):
    """A pool where three columns make y and the others are noise."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(rows, 3 + noise))
    y = 3.0 * x[:, 0] - 2.0 * x[:, 1] + x[:, 2] + rng.normal(scale=0.1, size=rows)
    return x, y


# -- the point of the module ----------------------------------------------------------------
def test_the_cascade_finds_the_columns_that_made_y():
    x, y = a_pool()
    kept = select(x, y, task="regression", max_features=5).kept
    assert set(SIGNAL) <= set(kept.tolist())
    assert len(kept) <= 5


def test_a_shuffled_y_keeps_no_column_in_particular():
    """Selection must not invent signal: with y shuffled, the three real columns are no more
    likely to be kept than any others."""
    x, y = a_pool()
    shuffled = np.random.default_rng(1).permutation(y)
    kept = set(select(x, shuffled, task="regression", max_features=5).kept.tolist())
    assert len(kept & set(SIGNAL)) <= 1


def test_the_holes_have_to_be_filled_first():
    x, y = a_pool(rows=40, noise=3)
    x[0, 0] = np.nan
    with pytest.raises(ValueError, match="holes filled first"):
        select(x, y, task="regression")
    filled = Prepare.fit(x).filled(x)
    assert len(select(filled, y, task="regression").kept)


# -- the stages -----------------------------------------------------------------------------
def test_flat_columns_are_the_ones_with_nothing_in_them():
    x = np.array([
        [1.0, 5.0, np.nan, 2.0],
        [2.0, 5.0, np.nan, 3.0],
        [3.0, 5.0, np.nan, 4.0],
    ])  # fmt: skip
    assert flat_columns(x).tolist() == [1, 2]


def test_a_flat_column_never_survives_the_cascade():
    x, y = a_pool(rows=60, noise=5)
    x[:, 4] = 7.0
    assert 4 not in select(x, y, task="regression", max_features=8).kept.tolist()


@pytest.mark.parametrize("kind", list(FILTERS))
def test_every_filter_ranks_the_signal_above_the_noise(kind):
    x, y = a_pool(rows=120, noise=20)
    scores = univariate_scores(x, y, "regression", kind)
    assert scores[list(SIGNAL)].min() > np.median(scores[3:])


def test_an_unknown_filter_is_refused():
    x, y = a_pool(rows=30, noise=2)
    with pytest.raises(ValueError, match="unknown filter"):
        univariate_scores(x, y, "regression", "astrology")


def test_a_copy_of_a_column_is_dropped_and_the_better_one_kept():
    rng = np.random.default_rng(2)
    first = rng.normal(size=80)
    x = np.column_stack([first, first + 1e-9, rng.normal(size=80)])
    kept = decorrelated(x, np.array([0, 1, 2]), threshold=0.95)
    assert kept.tolist() == [0, 2], "the column offered first wins its pair"
    assert decorrelated(x, np.array([1, 0, 2]), threshold=0.95).tolist() == [1, 2]


def test_decorrelation_keeps_columns_that_merely_resemble_each_other():
    rng = np.random.default_rng(3)
    first = rng.normal(size=200)
    close = 0.5 * first + rng.normal(size=200)  # correlated, but not the same information
    kept = decorrelated(np.column_stack([first, close]), np.array([0, 1]), threshold=0.95)
    assert kept.tolist() == [0, 1]


def test_the_embedded_pass_cuts_to_the_columns_a_forest_uses():
    x, y = a_pool(rows=150, noise=20)
    columns = np.arange(x.shape[1])
    kept = embedded(x, y, columns, "regression", keep=4, seed=0)
    assert len(kept) == 4
    assert set(SIGNAL) <= set(kept.tolist())


def test_the_embedded_pass_does_nothing_when_there_is_nothing_to_cut():
    x, y = a_pool(rows=40, noise=1)
    columns = np.arange(x.shape[1])
    assert embedded(x, y, columns, "regression", keep=99, seed=0).tolist() == columns.tolist()


# -- what it reports ------------------------------------------------------------------------
def test_every_stage_reports_what_it_dropped():
    x, y = a_pool(rows=100, noise=40)
    x[:, 5] = 0.0  # flat
    x[:, 6] = x[:, 0]  # a copy of a signal column
    chosen = select(x, y, task="regression", max_features=4)
    assert set(chosen.dropped) == {"flat", "filter", "correlated", "unused"}
    assert chosen.dropped["flat"] == 1
    assert sum(chosen.dropped.values()) + len(chosen) == x.shape[1]


def test_the_scores_cover_the_whole_pool():
    x, y = a_pool(rows=60, noise=10)
    chosen = select(x, y, task="regression", max_features=5)
    assert chosen.scores is not None
    assert len(chosen.scores) == x.shape[1]
    assert np.isfinite(chosen.scores).all()


def test_kept_columns_come_back_in_order():
    x, y = a_pool(rows=80, noise=30)
    kept = select(x, y, task="regression", max_features=6).kept
    assert kept.tolist() == sorted(kept.tolist())


def test_the_same_seed_selects_the_same_columns():
    x, y = a_pool(rows=100, noise=30)
    first = select(x, y, task="regression", max_features=5, seed=3).kept
    again = select(x, y, task="regression", max_features=5, seed=3).kept
    assert first.tolist() == again.tolist()


def test_a_classification_selects_too():
    x, _ = a_pool(rows=160, noise=30)
    y = (x[:, 0] + x[:, 1] > 0).astype(np.int64)
    kept = select(x, y, task="classification", max_features=4, kind=DEFAULT_FILTER).kept
    assert {0, 1} <= set(kept.tolist())


def test_an_empty_pool_selects_nothing_rather_than_failing():
    chosen = select(np.zeros((10, 3)), np.arange(10.0), task="regression")
    assert len(chosen) == 0
    assert chosen.dropped["flat"] == 3


def test_mismatched_rows_are_refused():
    with pytest.raises(ValueError, match="rows of features"):
        select(np.zeros((5, 2)), np.zeros(4), task="regression")


def test_a_selection_knows_its_size():
    assert len(Selection(np.array([1, 4, 9]))) == 3
