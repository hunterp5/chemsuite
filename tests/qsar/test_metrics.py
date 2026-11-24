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
"""The metrics, and which direction counts as better."""

from __future__ import annotations

import math

import numpy as np
import pytest

from chemsuite.qsar import metrics


def test_regression_metrics_agree_with_scikit_learn():
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    rng = np.random.default_rng(0)
    y = rng.normal(size=60)
    predicted = y + rng.normal(scale=0.3, size=60)
    got = metrics.regression_scores(y, predicted)
    assert got["r2"] == pytest.approx(r2_score(y, predicted))
    assert got["mae"] == pytest.approx(mean_absolute_error(y, predicted))
    assert got["rmse"] == pytest.approx(math.sqrt(mean_squared_error(y, predicted)))
    assert got["q2"] == pytest.approx(got["r2"])  # out of fold already


def test_a_perfect_regression_and_a_useless_one():
    y = np.linspace(0, 1, 30)
    perfect = metrics.regression_scores(y, y.copy())
    assert perfect["rmse"] == pytest.approx(0.0)
    assert perfect["r2"] == pytest.approx(1.0)
    assert perfect["spearman"] == pytest.approx(1.0)
    mean_only = metrics.regression_scores(y, np.full_like(y, y.mean()))
    assert mean_only["r2"] == pytest.approx(0.0)
    assert math.isnan(mean_only["spearman"])  # a flat prediction has no rank order


def test_rows_without_a_value_are_left_out():
    y = np.array([1.0, 2.0, 3.0, np.nan])
    predicted = np.array([1.0, 2.0, 3.0, 99.0])
    assert metrics.regression_scores(y, predicted)["rmse"] == pytest.approx(0.0)
    missed = np.array([1.0, 2.0, 3.0, np.nan])
    assert metrics.regression_scores(y, missed)["rmse"] == pytest.approx(0.0)


def test_too_few_rows_to_score_is_nan_not_an_error():
    got = metrics.regression_scores(np.array([1.0]), np.array([1.0]))
    assert set(got) == set(metrics.REGRESSION_METRICS)
    assert all(math.isnan(value) for value in got.values())


def test_classification_metrics_agree_with_scikit_learn():
    from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

    rng = np.random.default_rng(1)
    y = (rng.random(80) < 0.3).astype(np.int64)
    probability = np.clip(0.2 + 0.6 * y + rng.normal(scale=0.2, size=80), 0.01, 0.99)
    scores = np.column_stack([1 - probability, probability])
    predicted = (probability > 0.5).astype(np.int64)
    got = metrics.classification_scores(y, predicted, scores, np.array(["no", "yes"]))
    assert got["roc_auc"] == pytest.approx(roc_auc_score(y, probability))
    assert got["pr_auc"] == pytest.approx(average_precision_score(y, probability))
    assert got["mcc"] == pytest.approx(matthews_corrcoef(y, predicted))
    assert 0.0 <= got["calibration"] <= 1.0


def test_a_calibrated_model_scores_better_than_a_confident_wrong_one():
    y = np.array([0] * 50 + [1] * 50)
    honest = np.where(y == 1, 0.8, 0.2)
    overconfident = np.where(y == 1, 0.999, 0.4)
    classes = np.array(["no", "yes"])
    good = metrics.classification_scores(
        y, (honest > 0.5).astype(np.int64), np.column_stack([1 - honest, honest]), classes
    )
    bad = metrics.classification_scores(
        y,
        (overconfident > 0.5).astype(np.int64),
        np.column_stack([1 - overconfident, overconfident]),
        classes,
    )
    assert good["calibration"] < bad["calibration"]
    assert good["roc_auc"] == pytest.approx(bad["roc_auc"])  # the ranking is the same


def test_a_decision_value_has_no_calibration():
    y = np.array([0, 0, 1, 1])
    decision = np.array([-2.0, -1.0, 1.0, 2.5])
    got = metrics.classification_scores(
        y, (decision > 0).astype(np.int64), np.column_stack([-decision, decision]),
        np.array(["no", "yes"]),
    )  # fmt: skip
    assert math.isnan(got["calibration"])
    assert got["roc_auc"] == pytest.approx(1.0)


def test_three_classes_are_scored_one_against_the_rest():
    rng = np.random.default_rng(2)
    y = np.repeat([0, 1, 2], 20)
    scores = rng.random((60, 3))
    scores[np.arange(60), y] += 1.5  # mostly right
    scores /= scores.sum(axis=1, keepdims=True)
    predicted = scores.argmax(axis=1)
    got = metrics.classification_scores(y, predicted, scores, np.array(["a", "b", "c"]))
    assert got["roc_auc"] > 0.8
    assert got["pr_auc"] > 0.6
    assert got["mcc"] > 0.5


def test_one_class_in_the_fold_cannot_be_scored():
    y = np.zeros(10, dtype=np.int64)
    scores = np.column_stack([np.ones(10), np.zeros(10)])
    got = metrics.classification_scores(y, y, scores, np.array(["no", "yes"]))
    assert all(math.isnan(value) for value in got.values())


def test_the_objective_follows_the_task():
    assert metrics.objective_of("regression") == "rmse"
    assert metrics.objective_of("classification") == "roc_auc"
    assert metrics.objective_of("classification", classes=3) == "mcc"


@pytest.mark.parametrize(
    ("metric", "value", "than", "better"),
    [
        ("rmse", 0.4, 0.5, True),
        ("rmse", 0.6, 0.5, False),
        ("r2", 0.9, 0.5, True),
        ("roc_auc", 0.5, 0.9, False),
        ("rmse", math.nan, 0.5, False),
        ("rmse", 0.5, math.nan, True),
        ("rmse", math.nan, math.nan, False),
    ],
)
def test_which_way_is_better(metric, value, than, better):
    assert metrics.is_better(metric, value, than) is better


def test_the_worst_value_is_beaten_by_anything():
    for metric in ("rmse", "r2"):
        assert metrics.is_better(metric, 0.42, metrics.worst(metric))


def test_every_metric_has_a_title():
    for metric in metrics.REGRESSION_METRICS + metrics.CLASSIFICATION_METRICS:
        assert metrics.TITLES[metric]


def test_scores_for_needs_the_classes_of_a_classification():
    with pytest.raises(ValueError, match="class scores"):
        metrics.scores_for("classification", np.array([0, 1]), np.array([0, 1]))
