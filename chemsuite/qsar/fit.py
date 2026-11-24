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
"""The one call a tool makes: a feature pool and some values in, a saved model out.

It runs the search, writes the winner out as arrays and checks them against scikit-learn,
calibrates the intervals on the predictions the search already made out of fold, and packs the
lot into a :class:`~.model.Model`. Everything the tool needs for its columns comes back
alongside, in the :class:`~.search.Outcome`.

This is also the only module in the package that knows about all the others, which keeps the
dependency arrows pointing one way: the apply path (``model``, ``apply``, ``trees``,
``prepare``, ``neighbors``, ``conformal``) never imports the fitting side.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import date

import numpy as np

from chemsuite.core.jobs import CancelToken

from . import export, search
from .conformal import Reliability
from .model import Model, Prediction, Training
from .recipe import Recipe
from .search import NoModel, Outcome, Progress
from .splits import as_groups, cluster_groups, scaffold_groups

__all__ = [
    "Fitted",
    "NoModel",
    "as_groups",
    "cluster_groups",
    "fit",
    "scaffold_groups",
]


@dataclass(frozen=True)
class Fitted:
    """A finished fit: the model to keep, the run it came out of, and the honest predictions.

    *held* is what the model would have said about each training row had it never seen it —
    the out-of-fold prediction, calibrated. It is what a table should show, and it is not the
    same as applying the finished model to those rows, which would flatter it.
    """

    model: Model
    outcome: Outcome
    #: the rows that were fitted on, as indices into the feature pool
    rows: np.ndarray
    #: the out-of-fold prediction for :attr:`rows`, in that order
    held: Prediction
    #: the measured value for :attr:`rows` (a class index for a classification)
    measured: np.ndarray

    @property
    def scores(self) -> dict[str, float]:
        return self.outcome.scores

    @property
    def fold(self) -> np.ndarray:
        """Which outer fold held each of :attr:`rows` out (-1 for none)."""
        return self.outcome.fold[self.rows]


def fit(
    pool: np.ndarray,
    names: Sequence[str],
    y: np.ndarray,
    rows: np.ndarray,
    recipe: Recipe,
    *,
    keys: Sequence[str],
    classes: np.ndarray | None = None,
    groups: np.ndarray | None = None,
    order: np.ndarray | None = None,
    label: str = "",
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Fitted:
    """Fit a model on *rows* of *pool*, and come back with something that can be saved.

    *pool* is the whole feature matrix (:func:`chemsuite.qsar.features.build`), *names* its
    column names, *y* the target — class indices with *classes* for a classification — and
    *keys* one identifier per row, so the model can be refined later with more data.
    """
    rows = np.sort(np.asarray(rows, dtype=np.int64))
    values = np.asarray(y, dtype=np.float64)
    outcome = search.run(
        pool,
        y,
        rows,
        recipe,
        classes=classes,
        groups=groups,
        order=order,
        progress=progress,
        cancel=cancel,
    )
    design = outcome.model.prepare.design(pool, outcome.model.selection.kept)
    weights, _difference = export.exported(outcome.model, design[rows])
    reliability = _calibrate(outcome, design, y, rows, recipe, classes)
    keys = np.asarray(keys, dtype=np.str_)
    if len(keys) != len(y):
        raise ValueError(f"{len(keys)} keys for {len(y)} rows")
    model = Model(
        recipe=recipe,
        weights=weights,
        prepare=outcome.model.prepare,
        kept=outcome.model.selection.kept,
        names=tuple(names),
        candidate=outcome.model.candidate,
        params=outcome.model.params,
        scores=outcome.scores,
        reliability=reliability,
        training=Training(keys[rows], np.asarray(y, dtype=np.float64)[rows]),
        classes=classes,
        leaderboard=tuple(
            {"candidate": s.candidate, "title": s.title, "params": s.params, "scores": s.scores}
            for s in outcome.leaderboard
        ),
        lineage=(
            {
                "date": date.today().isoformat(),
                "note": "fitted",
                "rows": len(rows),
                "scores": outcome.scores,
            },
        ),
        created=date.today().isoformat(),
        objective=outcome.objective,
        label=label,
    )
    return Fitted(model, outcome, rows, _held(model, outcome, design, rows), values[rows])


def _held(model: Model, outcome: Outcome, design: np.ndarray, rows: np.ndarray) -> Prediction:
    """The out-of-fold predictions the search made, wrapped with the model's calibration.

    Deliberately not a re-prediction: these numbers come from folds that did not contain
    their own row, which is the whole basis of the reported score.
    """
    held = design[rows]
    reliability = model.reliability
    distance = reliability.distance(held)
    inside = distance <= reliability.domain_edge
    if model.recipe.task == "regression":
        predicted = outcome.predicted[rows]
        low, high = reliability.interval(predicted, held)
        return Prediction(predicted, low=low, high=high, distance=distance, inside=inside)
    if outcome.class_scores is None:
        raise NoModel("a classification needs its class scores")
    scores = outcome.class_scores[rows]
    return Prediction(
        np.argmax(scores, axis=1),
        scores=scores,
        sets=reliability.class_sets(scores),
        distance=distance,
        inside=inside,
    )


def _calibrate(
    outcome: Outcome,
    design: np.ndarray,
    y: np.ndarray,
    rows: np.ndarray,
    recipe: Recipe,
    classes: np.ndarray | None,
) -> Reliability:
    """Calibrate on the out-of-fold predictions the search already made.

    Rows no fold held out — only a column split leaves any — cannot calibrate anything, so
    they are left out of the quantile but stay in the domain's own reckoning.
    """
    held = rows[outcome.fold[rows] >= 0]
    if not len(held):
        raise NoModel("no row was held out, so nothing can be calibrated")
    if recipe.task == "regression":
        return Reliability.for_regression(
            design[held],
            outcome.predicted[held],
            np.asarray(y, dtype=np.float64)[held],
            confidence=recipe.confidence,
        )
    if outcome.class_scores is None or classes is None:
        raise NoModel("a classification needs its class scores to be calibrated")
    return Reliability.for_classification(
        design[held],
        outcome.class_scores[held],
        np.asarray(y, dtype=np.int64)[held],
        classes=len(classes),
        confidence=recipe.confidence,
    )
