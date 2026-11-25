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
"""Fitted with scikit-learn, applied in numpy — and the two must agree exactly.

The test that carries the module runs **every candidate in the zoo**, for regression and for
classification, and demands that the exported arrays reproduce scikit-learn to within 1e-9.
That is what lets a saved model be arrays instead of a pickle: the numpy path is not an
approximation of the estimator, it is the estimator.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from chemsuite.qsar import apply, export, search, zoo
from chemsuite.qsar.apply import ModelError
from chemsuite.qsar.recipe import Recipe


def a_recipe(task: str = "regression", **changes) -> Recipe:
    values = {
        "task": task,
        "target": "y",
        "columns": tuple(f"c{i}" for i in range(4)),
        "split": "random",
        "budget": "quick",
        "max_features": 8,
        "seed": 0,
    }
    values.update(changes)
    return Recipe(**values)


def a_dataset(rows: int = 120, features: int = 8, classes: int = 0, seed: int = 0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(rows, features))
    raw = 2.0 * x[:, 0] - x[:, 1] + 0.5 * x[:, 2] + rng.normal(scale=0.2, size=rows)
    if not classes:
        return x, raw, None
    edges = np.quantile(raw, np.linspace(0, 1, classes + 1)[1:-1])
    y = np.digitize(raw, edges).astype(np.int64)
    return x, y, np.array([f"class{i}" for i in range(classes)])


def fit(
    candidate: zoo.Candidate, x, y, classes, recipe: Recipe, params: dict | None = None
) -> search.Pipeline:
    rows = np.arange(len(y))
    basis = search.basis_for(x, y, rows, recipe)
    rng = np.random.default_rng(1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return search.fit_pipeline(
            x, y, rows, candidate, params or candidate.draw(rng), recipe, basis, classes
        )


def design_of(model: search.Pipeline, x) -> np.ndarray:
    return model.prepare.design(x, model.selection.kept)


# -- the point of the module ----------------------------------------------------------------
@pytest.mark.parametrize("key", [c.key for c in zoo.for_task("regression")])
def test_every_regression_candidate_survives_the_round_trip(key):
    x, y, _ = a_dataset()
    recipe = a_recipe()
    model = fit(zoo.candidate("regression", key), x, y, None, recipe)
    weights, difference = export.exported(model, design_of(model, x))
    assert difference <= export.TOLERANCE, f"{key} is off by {difference:.3g}"
    assert weights.family in apply.FAMILIES


@pytest.mark.parametrize("key", [c.key for c in zoo.for_task("classification")])
def test_every_classification_candidate_survives_the_round_trip(key):
    x, y, classes = a_dataset(classes=2)
    recipe = a_recipe("classification")
    model = fit(zoo.candidate("classification", key), x, y, classes, recipe)
    weights, difference = export.exported(model, design_of(model, x))
    assert difference <= export.TOLERANCE, f"{key} is off by {difference:.3g}"
    scores = apply.values(weights, design_of(model, x))
    assert scores.shape == (len(y), 2)


@pytest.mark.parametrize("key", ["hgb", "rf", "logistic", "knn", "mlp"])
def test_three_classes_survive_the_round_trip(key):
    """Boosting grows one tree per class, a forest's leaf holds all three shares, and the
    linear family turns from a logistic into a softmax."""
    x, y, classes = a_dataset(rows=180, classes=3)
    recipe = a_recipe("classification")
    model = fit(zoo.candidate("classification", key), x, y, classes, recipe)
    weights, difference = export.exported(model, design_of(model, x))
    assert difference <= export.TOLERANCE, f"{key} is off by {difference:.3g}"
    assert apply.values(weights, design_of(model, x)).shape == (len(y), 3)


def test_support_vectors_are_refused_above_two_classes():
    """Not a silent approximation: scikit-learn's multiclass margin is a vote between pairs."""
    x, y, classes = a_dataset(rows=150, classes=3)
    model = fit(zoo.candidate("classification", "svc"), x, y, classes, a_recipe("classification"))
    with pytest.raises(ModelError, match="two classes only"):
        export.export(model)


def test_a_wrong_array_is_caught_rather_than_written():
    x, y, _ = a_dataset()
    recipe = a_recipe()
    model = fit(zoo.candidate("regression", "ridge"), x, y, None, recipe)
    weights = export.export(model)
    broken = apply.Weights(
        weights.family,
        weights.link,
        {**weights.arrays, "intercept": weights.arrays["intercept"] + 1.0},
        weights.meta,
    )
    with pytest.raises(ModelError, match="disagrees with scikit-learn"):
        export.verify(model, broken, design_of(model, x))


def test_the_predicted_class_is_the_one_scikit_learn_picks():
    x, y, classes = a_dataset(rows=150, classes=3)
    recipe = a_recipe("classification")
    for key in ("hgb", "rf", "logistic", "knn"):
        model = fit(zoo.candidate("classification", key), x, y, classes, recipe)
        weights = export.export(model)
        design = design_of(model, x)
        mine = apply.predicted(weights, design)
        theirs = model.estimator.predict(design)
        assert mine.tolist() == theirs.tolist(), key


# -- the arrays themselves ------------------------------------------------------------------
def test_partial_least_squares_folds_its_centering_into_the_intercept():
    x, y, _ = a_dataset()
    recipe = a_recipe()
    model = fit(zoo.candidate("regression", "pls"), x, y, None, recipe)
    weights = export.export(model)
    assert weights.family == apply.LINEAR, "once centered, PLS is a straight line"
    assert weights.link == apply.IDENTITY


def test_a_forest_leaf_holds_one_share_per_class():
    x, y, classes = a_dataset(rows=150, classes=3)
    model = fit(zoo.candidate("classification", "rf"), x, y, classes, a_recipe("classification"))
    weights = export.export(model)
    assert weights.arrays["node_value"].shape[1] == 3
    assert weights.meta["classes"] == 3


def test_boosting_records_which_class_each_tree_is_for():
    x, y, classes = a_dataset(rows=150, classes=3)
    model = fit(zoo.candidate("classification", "hgb"), x, y, classes, a_recipe("classification"))
    weights = export.export(model)
    assert set(weights.arrays["tree_output"].tolist()) == {0, 1, 2}
    assert len(weights.arrays["baseline"]) == 3


def test_a_neighborhood_model_carries_its_training_rows():
    x, y, _ = a_dataset()
    recipe = a_recipe()
    model = fit(zoo.candidate("regression", "knn"), x, y, None, recipe)
    weights = export.export(model)
    assert weights.arrays["design"].shape == (len(y), len(model.selection.kept))
    assert weights.meta["weights"] in ("uniform", "distance")


def test_an_unknown_family_is_refused():
    x, y, _ = a_dataset()
    recipe = a_recipe()
    model = fit(zoo.candidate("regression", "ridge"), x, y, None, recipe)
    pretend = search.Pipeline(
        model.task, model.candidate, model.params, model.prepare, model.selection,
        model.estimator, model.classes,
    )  # fmt: skip
    assert export.export(pretend).family == apply.LINEAR  # the real one still works
    with pytest.raises(ModelError, match="unknown model family"):
        apply.Weights("astrology", apply.IDENTITY)


def test_features_must_be_a_matrix():
    with pytest.raises(ModelError, match="two-dimensional"):
        apply.outputs(
            apply.Weights(apply.LINEAR, apply.IDENTITY,
                          {"coef": np.ones((1, 2)), "intercept": np.zeros(1)}),
            np.ones(2),
        )  # fmt: skip


def test_a_missing_array_says_which_one():
    with pytest.raises(ModelError, match="'coef'"):
        apply.outputs(apply.Weights(apply.LINEAR, apply.IDENTITY, {}), np.ones((2, 2)))


def test_the_search_never_offers_a_model_it_cannot_write_out():
    """The guard that keeps a three-class fit from searching for minutes and then failing:
    support vectors are simply not among the candidates."""
    two = {c.key for c in zoo.for_task("classification", classes=2)}
    three = {c.key for c in zoo.for_task("classification", classes=3)}
    assert "svc" in two
    assert "svc" not in three
    assert three == two - {"svc"}


def test_a_forest_keeps_scikit_learns_missing_value_rule() -> None:
    """The export carries the estimator's own missing-value routing, not zeros.

    ``Prepare.design`` fills every hole before the estimator sees it, so no prediction made
    through the pipeline reaches this branch today, and writing zeros was harmless. It is
    still the estimator's rule and not ours to invent: scikit-learn sets ``missing_go_to_left``
    on its nodes whatever the training data held, and an archive that says something different
    from the forest it was exported from is a trap for whoever next feeds a design through
    ``apply`` directly.
    """
    rng = np.random.default_rng(0)
    x = rng.normal(size=(120, 12))
    y = 2 * x[:, 0] - x[:, 1] + rng.normal(scale=0.2, size=120)
    x[rng.random(x.shape) < 0.08] = np.nan

    model = fit(zoo.candidate("regression", "rf"), x, y, None, a_recipe("regression"))
    weights, difference = export.exported(model, design_of(model, x))

    assert difference == 0.0, "the exported forest is not scikit-learn's forest"
    kept = weights.arrays["node_missing_left"]
    theirs = np.concatenate([t.tree_.missing_go_to_left for t in model.estimator.estimators_])
    assert np.array_equal(kept, theirs.astype(np.int8)), "the routing is not the estimator's"


def test_a_forest_walks_the_way_scikit_learn_walks_it() -> None:
    """A scikit-learn tree casts its input to float32 before walking it, so the numpy walk
    has to compare the same way.

    The probe puts a row exactly on each threshold the forest uses, at the float64 value just
    below it: wherever that threshold is not float32-representable, the two precisions take
    different branches. Comparing in float64 sent 309 of these rows the other way, by up to
    0.045 — rare and small, which is why ordinary random data never caught it.
    """
    rng = np.random.default_rng(3)
    x = rng.normal(size=(200, 6))
    y = 2 * x[:, 0] - x[:, 1] + rng.normal(scale=0.1, size=200)

    params = {"max_depth": 8, "max_features": "sqrt", "min_samples_leaf": 2}
    model = fit(zoo.candidate("regression", "rf"), x, y, None, a_recipe("regression"), params)
    design = design_of(model, x)
    weights = export.export(model)

    inner = np.flatnonzero(weights.arrays["node_leaf"] == 0)
    feature = weights.arrays["node_feature"]
    threshold = weights.arrays["node_threshold"]
    probe = np.repeat(design[:1], len(inner), axis=0)
    for i, node in enumerate(inner):
        probe[i, feature[node]] = np.nextafter(threshold[node], -np.inf)

    difference = np.abs(apply.values(weights, probe) - model.estimator.predict(probe)).max()

    assert difference == 0.0, "the exported forest walks a row the estimator would not"
