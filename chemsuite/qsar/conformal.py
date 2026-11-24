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
"""How much to trust one prediction: an interval, a class set, and the distance to what the
model has seen.

A number with no error bar is the thing chemists rightly distrust about QSAR. Conformal
prediction gives the bar a meaning that holds whatever the model is: calibrate on predictions
the model did not fit, take a quantile of how wrong it was, and an interval built from that
quantile covers the true value at least as often as asked — distribution-free, with no
assumption that the errors are normal.

Three things live here.

- **Regression intervals**, normalized: the nonconformity of a calibration row is its error
  divided by a local difficulty — the mean error the model made on that row's nearest
  training neighbors. So the interval is wide in a sparse, badly modelled corner and narrow
  where the model has company, instead of one width everywhere.
- **Class sets**, Mondrian: a quantile per class rather than one overall, so a rare class is
  not swallowed by a common one. A row's set can hold one class, several (the model cannot
  separate them at this confidence) or none (it looks like nothing the model was trained on).
- **The applicability domain**: the distance from a row to its nearest training neighbors, and
  whether that is further than the training set's own rows usually sit. This is not conformal
  and makes no guarantee; it answers a different question — "is this compound the kind of
  thing the model knows about" — and the honest reading is that a row outside the domain has
  an interval the calibration never covered.

The calibration rows are the **out-of-fold** predictions the outer loop already produced
(:mod:`chemsuite.qsar.search`), so nothing is held back from fitting and nothing a prediction
is calibrated on was fitted on. The one compromise is that those predictions came from the
outer folds' models while the interval is put around the final model's: standard practice for
cross-conformal prediction, and the reason the coverage test allows a little slack.

Everything is arrays and arithmetic, so :mod:`chemsuite.qsar.apply` replays it without
scikit-learn.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .neighbors import mean_distance, mean_of_neighbors

#: Neighbors used for the local difficulty and the domain distance.
NEIGHBORS = 10
#: The fewest calibration rows a normalized interval needs; below this the width is constant,
#: because a difficulty estimated from a handful of neighbors is noise.
MIN_FOR_LOCAL = 50
#: The domain's edge, as a quantile of the training set's own neighbor distances.
DOMAIN_QUANTILE = 0.95
#: Keeps a difficulty of zero from making an interval of zero width.
BETA_SHARE = 0.1


def conformal_quantile(scores: np.ndarray, confidence: float) -> float:
    """The nonconformity quantile that covers *confidence* of future rows.

    It is the ``ceil((n + 1) · confidence)``-th smallest score — an order statistic, not an
    interpolated quantile. The ``n + 1`` is the finite-sample correction that makes the
    guarantee hold for the *next* row rather than on average over an infinite calibration set,
    and interpolating between two scores would quietly give it away again.

    With too few rows for that rank to exist, the largest score is the quantile and coverage is
    as good as a small set allows.
    """
    scores = np.asarray(scores, dtype=np.float64)
    scores = scores[np.isfinite(scores)]
    n = len(scores)
    if n == 0:
        return math.inf
    rank = math.ceil((n + 1) * confidence)
    if rank > n:
        return float(scores.max())
    return float(np.partition(scores, rank - 1)[rank - 1])


@dataclass(frozen=True)
class Reliability:
    """What a fit learned about its own errors, as arrays a model archive can carry.

    *design* is the training rows as the final model sees them — filled, selected and scaled —
    which is the space both the difficulty and the domain distance are measured in.
    """

    confidence: float
    k: int
    design: np.ndarray
    #: regression: the out-of-fold absolute error of each training row
    residual: np.ndarray | None = None
    #: regression: the nonconformity quantile, and the floor under the local difficulty
    quantile: float = math.nan
    beta: float = 0.0
    local: bool = True
    #: classification: one nonconformity quantile per class
    class_quantiles: np.ndarray | None = None
    #: the distance beyond which a row is outside the domain
    domain_edge: float = math.inf

    # -- building -------------------------------------------------------------------------
    @classmethod
    def for_regression(
        cls,
        design: np.ndarray,
        predicted: np.ndarray,
        y: np.ndarray,
        *,
        confidence: float,
        k: int = NEIGHBORS,
    ) -> Reliability:
        """Calibrate on out-of-fold predictions: *design*, *predicted* and *y* are the
        training rows, in the same order."""
        design = np.ascontiguousarray(design, dtype=np.float64)
        error = np.abs(np.asarray(y, dtype=np.float64) - np.asarray(predicted, dtype=np.float64))
        local = len(design) >= MIN_FOR_LOCAL
        beta = 0.0
        if local:
            difficulty = mean_of_neighbors(design, design, error, k, skip=1)  # not itself
            beta = BETA_SHARE * float(np.nanstd(error))
            scores = error / np.maximum(difficulty + beta, 1e-12)
        else:
            scores = error
        return cls(
            confidence=confidence,
            k=k,
            design=design,
            residual=error,
            quantile=conformal_quantile(scores, confidence),
            beta=beta,
            local=local,
            domain_edge=_domain_edge(design, k),
        )

    @classmethod
    def for_classification(
        cls,
        design: np.ndarray,
        scores: np.ndarray,
        y: np.ndarray,
        *,
        classes: int,
        confidence: float,
        k: int = NEIGHBORS,
    ) -> Reliability:
        """Mondrian calibration: one quantile per class, from the out-of-fold score each row's
        true class was given."""
        design = np.ascontiguousarray(design, dtype=np.float64)
        scores = np.asarray(scores, dtype=np.float64).reshape(len(y), -1)
        y = np.asarray(y, dtype=np.int64)
        nonconformity = 1.0 - _as_weights(scores)
        quantiles = np.full(classes, math.inf)
        for label in range(classes):
            mine = y == label
            if mine.any():
                quantiles[label] = conformal_quantile(nonconformity[mine, label], confidence)
        return cls(
            confidence=confidence,
            k=k,
            design=design,
            class_quantiles=quantiles,
            domain_edge=_domain_edge(design, k),
        )

    # -- using ----------------------------------------------------------------------------
    def difficulty(self, design: np.ndarray) -> np.ndarray:
        """The local difficulty of each row: how wrong the model was near it."""
        if self.residual is None or not self.local:
            return np.zeros(len(design))
        return mean_of_neighbors(
            np.asarray(design, dtype=np.float64), self.design, self.residual, self.k
        )

    def interval(self, predicted: np.ndarray, design: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """The low and high end of each prediction's interval."""
        predicted = np.asarray(predicted, dtype=np.float64)
        if self.local:
            half = self.quantile * (self.difficulty(design) + self.beta)
        else:
            half = np.full(len(predicted), self.quantile)
        return predicted - half, predicted + half

    def class_sets(self, scores: np.ndarray) -> np.ndarray:
        """Rows × classes: whether each class stays in the row's prediction set."""
        if self.class_quantiles is None:
            raise ValueError("this model was not calibrated for classes")
        nonconformity = 1.0 - _as_weights(np.asarray(scores, dtype=np.float64))
        return nonconformity <= self.class_quantiles[None, :]

    def distance(self, design: np.ndarray) -> np.ndarray:
        """How far each row sits from the training set."""
        return mean_distance(np.asarray(design, dtype=np.float64), self.design, self.k)

    def inside(self, design: np.ndarray) -> np.ndarray:
        """Whether each row is close enough to the training set to be in the domain."""
        return self.distance(design) <= self.domain_edge

    # -- the archive ----------------------------------------------------------------------
    def as_arrays(self) -> dict[str, np.ndarray]:
        out = {"reliability_design": self.design}
        if self.residual is not None:
            out["reliability_residual"] = self.residual
        if self.class_quantiles is not None:
            out["reliability_class_quantiles"] = self.class_quantiles
        return out

    def as_json(self) -> dict[str, float | int | bool]:
        return {
            "confidence": self.confidence,
            "k": self.k,
            "quantile": self.quantile,
            "beta": self.beta,
            "local": self.local,
            "domain_edge": self.domain_edge,
        }

    @classmethod
    def from_saved(
        cls, arrays: dict[str, np.ndarray], meta: dict[str, float | int | bool]
    ) -> Reliability:
        if "reliability_design" not in arrays:
            raise ValueError("the model archive has no calibration rows")
        return cls(
            confidence=float(meta["confidence"]),
            k=int(meta["k"]),
            design=np.asarray(arrays["reliability_design"], dtype=np.float64),
            residual=arrays.get("reliability_residual"),
            quantile=float(meta["quantile"]),
            beta=float(meta["beta"]),
            local=bool(meta["local"]),
            class_quantiles=arrays.get("reliability_class_quantiles"),
            domain_edge=float(meta["domain_edge"]),
        )


def _as_weights(scores: np.ndarray) -> np.ndarray:
    """Score columns as rows that sum to one. A decision value becomes a weight by a softmax,
    so support vectors get prediction sets too."""
    if scores.size and np.all(scores >= 0) and np.allclose(scores.sum(axis=1), 1.0):
        return scores
    shifted = scores - scores.max(axis=1, keepdims=True)
    weights = np.exp(shifted)
    return weights / weights.sum(axis=1, keepdims=True)


def _domain_edge(design: np.ndarray, k: int) -> float:
    """Where the training set itself stops: the 95th percentile of its own rows' distances to
    their neighbors, each row's own place in the set left out."""
    if len(design) < 2:
        return math.inf
    gaps = mean_distance(design, design, k, skip=1)  # a row is not its own neighbor
    return float(np.quantile(gaps, DOMAIN_QUANTILE))
