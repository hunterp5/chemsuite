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
"""The search.

Two tests carry the module. **Recovery**: when y really is a function of three columns out of
sixty, the search finds them and scores well. **Leakage**: when y is shuffled, it scores
nothing — R² about zero, ROC-AUC about a half. A pipeline that selects features or tunes
settings on all the rows passes the first and fails the second, which is exactly why the
second one is here.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.qsar import search, zoo
from chemsuite.qsar.recipe import Recipe
from chemsuite.qsar.splits import as_groups

SIGNAL = (0, 1, 2)


def a_recipe(**changes) -> Recipe:
    """A quick recipe without the two candidates that dominate the clock.

    A random forest of 500 trees is most of the cost of a run, so the tests that are about the
    loops leave it out and one test below puts the whole Quick zoo back.
    """
    values = {
        "task": "regression",
        "target": "y",
        "columns": tuple(f"c{i}" for i in range(5)),
        "split": "random",
        "budget": "quick",
        "max_features": 6,
        "seed": 0,
        "without": ("rf", "mlp"),
    }
    values.update(changes)
    return Recipe(**values)


def a_pool(rows: int = 160, noise: int = 40, seed: int = 0, strength: float = 1.0):
    """Three columns make y; the rest is noise."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(rows, 3 + noise))
    y = strength * (3.0 * x[:, 0] - 2.0 * x[:, 1] + x[:, 2]) + rng.normal(scale=0.3, size=rows)
    return x, y


def classes_pool(rows: int = 200, noise: int = 30, seed: int = 0):
    x, raw = a_pool(rows, noise, seed)
    y = (raw > np.median(raw)).astype(np.int64)
    return x, y, np.array(["low", "high"])


# -- the two that matter --------------------------------------------------------------------
def test_it_recovers_a_model_that_is_really_there():
    x, y = a_pool()
    outcome = search.run(x, y, np.arange(len(y)), a_recipe())
    assert outcome.scores["r2"] > 0.8, outcome.scores
    assert set(SIGNAL) <= set(outcome.model.selection.kept.tolist())


def test_a_shuffled_target_scores_nothing():
    """The leakage test. Selection and tuning sit inside the folds, so there is nothing to
    find here and the search must report that it found nothing."""
    x, y = a_pool()
    shuffled = np.random.default_rng(99).permutation(y)
    outcome = search.run(x, shuffled, np.arange(len(y)), a_recipe())
    assert outcome.scores["r2"] < 0.15, outcome.scores
    assert abs(outcome.scores["spearman"]) < 0.3, outcome.scores


def test_a_shuffled_class_label_scores_a_coin_toss():
    x, y, classes = classes_pool()
    shuffled = np.random.default_rng(7).permutation(y)
    recipe = a_recipe(task="classification")
    outcome = search.run(x, shuffled, np.arange(len(y)), recipe, classes=classes)
    assert 0.3 < outcome.scores["roc_auc"] < 0.7, outcome.scores


# -- what comes back ------------------------------------------------------------------------
def test_every_row_is_predicted_out_of_fold_exactly_once():
    x, y = a_pool(rows=120, noise=10)
    outcome = search.run(x, y, np.arange(len(y)), a_recipe())
    assert (outcome.fold >= 0).all()
    assert sorted(set(outcome.fold.tolist())) == list(range(a_recipe().budget_spec.outer))
    assert np.isfinite(outcome.predicted).all()


def test_the_headline_score_is_the_out_of_fold_score_not_the_fit():
    x, y = a_pool(rows=120, noise=10)
    outcome = search.run(x, y, np.arange(len(y)), a_recipe())
    from chemsuite.qsar.metrics import regression_scores

    own_fit = regression_scores(y, search.predict(outcome.model, x, np.arange(len(y))))
    assert outcome.scores["r2"] == pytest.approx(
        regression_scores(y, outcome.predicted)["r2"], abs=1e-9
    )
    assert own_fit["rmse"] <= outcome.scores["rmse"] + 1e-9, "a fit flatters itself"


def test_the_leaderboard_covers_the_candidates_and_is_sorted():
    """The one test on the whole Quick zoo, random forest and neural network included."""
    x, y = a_pool(rows=140, noise=10)
    recipe = a_recipe(without=())
    outcome = search.run(x, y, np.arange(len(y)), recipe)
    tried = {c.key for c in zoo.for_task("regression", most=recipe.budget_spec.candidates)}
    assert "rf" in tried
    assert {s.candidate for s in outcome.leaderboard} <= tried
    objective = outcome.objective
    values = [s.scores[objective] for s in outcome.leaderboard]
    assert values == sorted(values), "the objective is an error, so lower comes first"
    assert all(s.title for s in outcome.leaderboard)


def test_the_kept_model_is_one_of_the_candidates_and_can_predict_again():
    x, y = a_pool(rows=120, noise=10)
    outcome = search.run(x, y, np.arange(len(y)), a_recipe())
    assert outcome.model.candidate in {c.key for c in zoo.for_task("regression")}
    again = search.predict(outcome.model, x, np.arange(10))
    assert again.shape == (10,)
    assert np.isfinite(again).all()


def test_a_classification_comes_back_with_class_scores():
    x, y, classes = classes_pool()
    outcome = search.run(x, y, np.arange(len(y)), a_recipe(task="classification"), classes=classes)
    assert outcome.class_scores is not None
    assert outcome.class_scores.shape == (len(y), 2)
    assert outcome.scores["roc_auc"] > 0.7
    assert set(outcome.predicted.tolist()) <= {0, 1}


def test_the_same_seed_gives_the_same_model():
    x, y = a_pool(rows=120, noise=10)
    first = search.run(x, y, np.arange(len(y)), a_recipe(seed=5))
    again = search.run(x, y, np.arange(len(y)), a_recipe(seed=5))
    assert first.model.candidate == again.model.candidate
    assert first.model.params == again.model.params
    assert first.scores == again.scores


def test_a_hole_in_the_pool_is_filled_not_fatal():
    x, y = a_pool(rows=140, noise=10)
    x[::7, 0] = np.nan
    outcome = search.run(x, y, np.arange(len(y)), a_recipe())
    assert math.isfinite(outcome.scores["rmse"])


def test_only_the_rows_given_are_used():
    x, y = a_pool(rows=160, noise=10)
    rows = np.arange(0, 120)
    y[120:] = np.nan  # the rows left out have no value at all
    outcome = search.run(x, y, rows, a_recipe())
    assert (outcome.fold[120:] == -1).all()
    assert math.isfinite(outcome.scores["r2"])


# -- the inner loop -------------------------------------------------------------------------
def test_choosing_returns_one_best_trial_per_candidate_best_first():
    x, y = a_pool(rows=120, noise=10)
    recipe = a_recipe()
    candidates = zoo.for_task("regression", most=3)
    trials, spent = search.choose(x, y, np.arange(len(y)), recipe, candidates=candidates)
    assert [t.candidate for t in trials] == sorted(
        {t.candidate for t in trials}, key=[t.candidate for t in trials].index
    )
    assert len({t.candidate for t in trials}) == len(trials) <= 3
    assert trials[0].score <= trials[-1].score  # RMSE: lower is better
    assert spent > 0


def test_halving_rounds_double_up_to_the_fold_count():
    assert search.halving_rounds(5) == [1, 2, 4, 5]
    assert search.halving_rounds(4) == [1, 2, 4]
    assert search.halving_rounds(1) == [1]


def test_the_basis_is_fitted_on_the_training_rows_alone():
    x, y = a_pool(rows=120, noise=10)
    rows = np.arange(60)
    basis = search.basis_for(x, y, rows, a_recipe())
    only = search.basis_for(x[:60], y[:60], rows, a_recipe())
    assert basis.prepare.fill.tolist() == only.prepare.fill.tolist()
    assert basis.selection.kept.tolist() == only.selection.kept.tolist()


def test_progress_counts_up_and_cancel_stops_the_search():
    x, y = a_pool(rows=120, noise=10)
    seen: list[tuple[int, int]] = []
    search.run(x, y, np.arange(len(y)), a_recipe(), progress=lambda d, t: seen.append((d, t)))
    assert seen and [d for d, _ in seen] == sorted(d for d, _ in seen)
    assert seen[-1][0] <= seen[-1][1]

    token = CancelToken()

    def stop(done: int, _total: int) -> None:
        if done >= 3:
            token.cancel()

    with pytest.raises(Cancelled):
        search.run(x, y, np.arange(len(y)), a_recipe(), progress=stop, cancel=token)


def test_a_scaffold_like_split_is_honored():
    x, y = a_pool(rows=150, noise=10)
    groups = as_groups([f"series{i % 5}" for i in range(len(y))])
    outcome = search.run(x, y, np.arange(len(y)), a_recipe(split="scaffold"), groups=groups)
    assert math.isfinite(outcome.scores["rmse"])


# -- refusals -------------------------------------------------------------------------------
def test_too_few_rows_is_refused_with_a_reason():
    x, y = a_pool(rows=12, noise=3)
    with pytest.raises(search.NoModel, match="fewest a model can be judged on"):
        search.run(x, y, np.arange(len(y)), a_recipe())


def test_a_pool_with_no_information_is_refused():
    x, y = a_pool(rows=80, noise=5)
    flat = np.ones_like(x)
    with pytest.raises(search.NoModel, match="no feature"):
        search.run(flat, y, np.arange(len(y)), a_recipe())


def test_a_classification_without_its_labels_is_refused():
    x, y, _ = classes_pool(rows=80, noise=5)
    with pytest.raises(ValueError, match="class labels"):
        search.run(x, y, np.arange(len(y)), a_recipe(task="classification"))
