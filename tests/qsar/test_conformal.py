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
"""Intervals, class sets and the applicability domain.

The test that carries the module is coverage: a 90% interval must cover about 90% of rows the
calibration never saw. That is the whole promise of conformal prediction, and it is checked
here on data where the right answer is known — including the hard case, where the noise is
bigger in one corner than another, and a single interval width would be a lie in both.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from chemsuite.qsar.conformal import (
    DOMAIN_QUANTILE,
    MIN_FOR_LOCAL,
    Reliability,
    conformal_quantile,
)


def a_model(rows: int = 400, features: int = 3, noise: float = 1.0, seed: int = 0):
    """Rows in feature space, the truth, and a model that is wrong by *noise*."""
    rng = np.random.default_rng(seed)
    design = rng.normal(size=(rows, features))
    truth = design[:, 0] * 2.0
    predicted = truth + rng.normal(scale=noise, size=rows)
    return design, predicted, truth


def uneven_noise(rows: int = 600, seed: int = 0):
    """Half the rows are easy, half are hard: the point of a local interval."""
    rng = np.random.default_rng(seed)
    design = np.column_stack([rng.normal(size=rows), rng.normal(size=rows)])
    hard = design[:, 0] > 0
    scale = np.where(hard, 3.0, 0.2)
    truth = design[:, 1]
    predicted = truth + rng.normal(scale=scale, size=rows)
    return design, predicted, truth, hard


# -- the promise ----------------------------------------------------------------------------
@pytest.mark.parametrize("confidence", [0.8, 0.9, 0.95])
def test_an_interval_covers_what_it_says_it_covers(confidence):
    design, predicted, truth = a_model(rows=500)
    fit = Reliability.for_regression(design[:250], predicted[:250], truth[:250],
                                     confidence=confidence)  # fmt: skip
    low, high = fit.interval(predicted[250:], design[250:])
    covered = float(np.mean((truth[250:] >= low) & (truth[250:] <= high)))
    assert covered == pytest.approx(confidence, abs=0.07), covered


def test_the_interval_is_wider_where_the_model_is_worse():
    design, predicted, truth, hard = uneven_noise()
    half = len(truth) // 2
    fit = Reliability.for_regression(design[:half], predicted[:half], truth[:half], confidence=0.9)
    low, high = fit.interval(predicted[half:], design[half:])
    width = high - low
    assert width[hard[half:]].mean() > 2 * width[~hard[half:]].mean()


def test_coverage_holds_in_both_corners_not_just_on_average():
    design, predicted, truth, hard = uneven_noise()
    half = len(truth) // 2
    fit = Reliability.for_regression(design[:half], predicted[:half], truth[:half], confidence=0.9)
    low, high = fit.interval(predicted[half:], design[half:])
    inside = (truth[half:] >= low) & (truth[half:] <= high)
    assert float(inside[hard[half:]].mean()) > 0.75
    assert float(inside[~hard[half:]].mean()) > 0.75


def test_a_class_set_covers_the_true_class():
    rng = np.random.default_rng(1)
    rows = 500
    design = rng.normal(size=(rows, 2))
    y = (design[:, 0] > 0).astype(np.int64)
    confident = np.clip(0.5 + 0.35 * np.where(y == 1, 1, -1) + rng.normal(scale=0.2, size=rows),
                        0.01, 0.99)  # fmt: skip
    scores = np.column_stack([1 - confident, confident])
    half = rows // 2
    fit = Reliability.for_classification(
        design[:half], scores[:half], y[:half], classes=2, confidence=0.9
    )
    sets = fit.class_sets(scores[half:])
    held = y[half:]
    assert float(np.mean(sets[np.arange(len(held)), held])) > 0.85
    assert sets.sum(axis=1).mean() < 2.0, "a useful set is not just both classes"


def test_a_rare_class_gets_its_own_quantile():
    rng = np.random.default_rng(2)
    y = np.array([0] * 270 + [1] * 30)
    design = rng.normal(size=(len(y), 2))
    confident = np.where(y == 1, 0.45, 0.2) + rng.normal(scale=0.05, size=len(y))
    scores = np.column_stack([1 - confident, confident])
    fit = Reliability.for_classification(design, scores, y, classes=2, confidence=0.9)
    assert fit.class_quantiles is not None
    assert fit.class_quantiles[1] > fit.class_quantiles[0], "the rare class is given more room"


# -- the quantile ---------------------------------------------------------------------------
def test_the_quantile_is_corrected_for_a_small_calibration_set():
    scores = np.arange(1.0, 11.0)  # ten scores, 1 to 10
    # ceil(11 * 0.9) / 10 = 1.0 -> the largest score
    assert conformal_quantile(scores, 0.9) == 10.0
    # ceil(11 * 0.8) / 10 = 0.9 -> the ninth
    assert conformal_quantile(scores, 0.8) == 9.0


def test_too_few_rows_to_promise_anything_gives_an_endless_interval():
    assert conformal_quantile(np.array([]), 0.9) == math.inf


def test_a_not_a_number_score_is_left_out():
    assert conformal_quantile(np.array([1.0, 2.0, np.nan, 3.0]), 0.5) == 2.0


# -- the applicability domain ---------------------------------------------------------------
def test_a_far_away_row_is_outside_the_domain():
    design, predicted, truth = a_model(rows=200)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    near = design[:20]
    far = design[:20] + 50.0
    assert fit.inside(near).mean() > 0.9
    assert not fit.inside(far).any()
    assert (fit.distance(far) > fit.distance(near)).all()


def test_the_domain_edge_is_where_the_training_set_itself_stops():
    design, predicted, truth = a_model(rows=300)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    own = fit.distance(design)
    assert float(np.mean(own <= fit.domain_edge)) == pytest.approx(DOMAIN_QUANTILE, abs=0.06)


def test_a_row_is_not_its_own_neighbor():
    """Otherwise every training row would sit at distance zero and the domain would be a
    point, and every row's difficulty would be its own error."""
    design, predicted, truth = a_model(rows=120)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    assert fit.distance(design).min() > 0.0
    assert math.isfinite(fit.domain_edge)


# -- the small-data path --------------------------------------------------------------------
def test_with_few_rows_the_width_is_constant_rather_than_pretend_local():
    design, predicted, truth = a_model(rows=MIN_FOR_LOCAL - 10)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    assert not fit.local
    low, high = fit.interval(predicted, design)
    width = high - low
    assert np.allclose(width, width[0])


def test_with_enough_rows_the_width_varies():
    design, predicted, truth = a_model(rows=MIN_FOR_LOCAL + 150)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    assert fit.local
    low, high = fit.interval(predicted, design)
    assert np.ptp(high - low) > 0.0


# -- the archive ----------------------------------------------------------------------------
def test_a_regression_calibration_round_trips_as_arrays():
    design, predicted, truth = a_model(rows=200)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    back = Reliability.from_saved(fit.as_arrays(), fit.as_json())
    assert np.allclose(back.interval(predicted, design), fit.interval(predicted, design))
    assert back.inside(design).tolist() == fit.inside(design).tolist()


def test_a_classification_calibration_round_trips_as_arrays():
    rng = np.random.default_rng(3)
    design = rng.normal(size=(120, 2))
    y = (design[:, 0] > 0).astype(np.int64)
    scores = np.column_stack([1 - y * 0.8 - 0.1, y * 0.8 + 0.1])
    fit = Reliability.for_classification(design, scores, y, classes=2, confidence=0.9)
    back = Reliability.from_saved(fit.as_arrays(), fit.as_json())
    assert back.class_sets(scores).tolist() == fit.class_sets(scores).tolist()


def test_an_archive_without_calibration_rows_is_refused():
    with pytest.raises(ValueError, match="no calibration rows"):
        Reliability.from_saved({}, {})


def test_class_sets_need_a_classification():
    design, predicted, truth = a_model(rows=120)
    fit = Reliability.for_regression(design, predicted, truth, confidence=0.9)
    with pytest.raises(ValueError, match="not calibrated for classes"):
        fit.class_sets(np.zeros((4, 2)))


def test_a_decision_value_still_gets_a_class_set():
    rng = np.random.default_rng(4)
    design = rng.normal(size=(200, 2))
    y = (design[:, 0] > 0).astype(np.int64)
    margin = np.where(y == 1, 1.5, -1.5) + rng.normal(scale=0.4, size=len(y))
    scores = np.column_stack([-margin, margin])  # a decision value, not a probability
    fit = Reliability.for_classification(design, scores, y, classes=2, confidence=0.9)
    sets = fit.class_sets(scores)
    assert float(np.mean(sets[np.arange(len(y)), y])) > 0.85
