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
"""Refining a model on the data measured since.

What has to hold: a repeated measurement replaces the old value rather than being added
twice; a new compound joins; a compound the table no longer has is reported, not silently
dropped; and the refitted model's score is its own, measured the same way the first was,
with one line of history per round.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from chemsuite.qsar.fit import fit
from chemsuite.qsar.recipe import Recipe
from chemsuite.qsar.refine import merge, refine
from chemsuite.qsar.search import NoModel

ROWS = 80


def a_recipe(**changes) -> Recipe:
    values = {
        "task": "regression",
        "target": "pIC50",
        "columns": tuple(f"c{i}" for i in range(4)),
        "split": "random",
        "budget": "quick",
        "max_features": 4,
        "seed": 0,
        "without": ("rf", "mlp", "knn"),
    }
    values.update(changes)
    return Recipe(**values)


def a_dataset(rows: int = ROWS, seed: int = 0):
    rng = np.random.default_rng(seed)
    pool = rng.normal(size=(rows, 8))
    y = 2 * pool[:, 0] - pool[:, 1] + rng.normal(scale=0.2, size=rows)
    names = [f"c:col{i}" for i in range(8)]
    keys = [f"KEY{i:04d}" for i in range(rows)]
    return pool, y, names, keys


def a_fit(rows: int = ROWS, **changes):
    pool, y, names, keys = a_dataset(rows)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fit(pool, names, y, np.arange(rows), a_recipe(**changes), keys=keys), pool, y, keys


def a_refit(fitted, pool, y, keys, rows=None, **changes):
    names = [f"c:col{i}" for i in range(pool.shape[1])]
    rows = np.arange(len(y)) if rows is None else rows
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return refine(fitted.model, pool, names, keys, y, rows, recipe=a_recipe(**changes))


# -- merging ---------------------------------------------------------------------------------
def test_the_same_compound_measured_again_replaces_its_old_value() -> None:
    fitted, pool, y, keys = a_fit()
    changed = y.copy()
    changed[:5] += 1.0  # the first five were remeasured
    merged = merge(fitted.model, keys, changed, np.arange(len(y)))
    assert len(merged) == ROWS, "nothing was added, so the set is the same size"
    assert merged.added == 0
    assert merged.changed == 5
    assert merged.carried == 0
    assert np.allclose(merged.y, changed)


def test_a_new_compound_joins_the_training_set() -> None:
    fitted, pool, y, keys = a_fit()
    extra_y = np.concatenate([y, y[:4]])
    extra_keys = [*keys, "NEW00001", "NEW00002", "NEW00003", "NEW00004"]
    merged = merge(fitted.model, extra_keys, extra_y, np.arange(len(extra_y)))
    assert merged.added == 4
    assert len(merged) == ROWS + 4
    assert "4 new" in merged.summary


def test_a_compound_the_table_no_longer_has_is_carried_and_counted() -> None:
    """Its features are not in this pool, so it cannot be refitted — and that is reported."""
    fitted, pool, y, keys = a_fit()
    kept = np.arange(ROWS - 10)  # ten compounds dropped out of the table
    merged = merge(fitted.model, keys, y, kept)
    assert merged.carried == 10
    assert merged.added == 0
    assert "10 carried over" in merged.summary
    assert (merged.rows[-10:] < 0).all(), "a carried row has no place in the pool"


def test_a_structure_that_cannot_be_read_counts_as_new() -> None:
    """With no key there is nothing to match on, so it is never mistaken for a known one."""
    fitted, pool, y, keys = a_fit()
    blanked = ["", *keys[1:]]
    merged = merge(fitted.model, blanked, y, np.arange(len(y)))
    assert merged.added == 1
    assert merged.carried == 1, "the compound it should have matched is carried instead"


def test_an_unchanged_refit_reports_no_change() -> None:
    fitted, pool, y, keys = a_fit()
    merged = merge(fitted.model, keys, y, np.arange(len(y)))
    assert merged.summary == "no change"
    assert (merged.added, merged.changed, merged.carried) == (0, 0, 0)


# -- refitting -------------------------------------------------------------------------------
def test_a_refit_adds_exactly_one_line_of_history() -> None:
    fitted, pool, y, keys = a_fit()
    again = a_refit(fitted, pool, y, keys)
    assert len(fitted.model.lineage) == 1
    assert len(again.model.lineage) == 2, "one line per round, not two"
    last = again.model.lineage[-1]
    assert last["note"].startswith("refit")
    assert last["rows"] == len(again.rows)
    third = a_refit(again, pool, y, keys)
    assert len(third.model.lineage) == 3


def test_a_refit_reports_its_own_score_not_the_old_one() -> None:
    fitted, pool, y, keys = a_fit()
    more_pool, more_y, _names, more_keys = a_dataset(ROWS + 40, seed=1)
    again = a_refit(fitted, more_pool, more_y, more_keys)
    assert again.model.scores == again.outcome.scores
    assert again.model.lineage[-1]["scores"] == again.scores
    assert again.model.scores["r2"] > 0.8, "the planted signal is still there"
    assert len(again.model.training) == len(again.rows)


def test_a_refit_keeps_the_recipe_and_the_name() -> None:
    fitted, pool, y, keys = a_fit()
    named = fitted.model.__class__(**{**fitted.model.__dict__, "label": "March model"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        again = refine(
            named,
            pool,
            [f"c:col{i}" for i in range(8)],
            keys,
            y,
            np.arange(len(y)),
        )
    assert again.model.label == "March model"
    assert again.model.recipe.target == "pIC50"
    assert again.model.created == named.created


def test_more_effort_can_be_asked_for_on_a_refit() -> None:
    """A model first fitted Quick can be refitted harder once the data is worth it.

    Only ridge is left in the zoo here: a full Standard search over five candidates is
    thousands of estimator fits, which is the engine's own tests' business, not this one's.
    """
    only_ridge = ("hgb", "svr", "pls", "elasticnet", "rf", "mlp", "knn")
    fitted, pool, y, keys = a_fit(budget="quick", without=only_ridge)
    again = a_refit(fitted, pool, y, keys, budget="standard", without=only_ridge)
    assert fitted.model.recipe.budget == "quick"
    assert again.model.recipe.budget == "standard"
    assert again.model.candidate == "ridge"


def test_refining_against_a_table_of_strangers_is_refused() -> None:
    """Nothing to fit: a model cannot be refined on compounds it has no features for."""
    fitted, pool, y, keys = a_fit()
    with pytest.raises(NoModel, match="nothing to fit"):
        a_refit(fitted, pool, y, keys, rows=np.array([], dtype=np.int64))
