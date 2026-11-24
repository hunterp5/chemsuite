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
"""How a model is scored, and which way is better.

Every number here is computed from **held-out** predictions: the search scores a candidate on
inner folds it was not fitted on, and the leaderboard reports the outer folds. R² on a model's
own training rows is not a metric and is never reported.

Q² is the same arithmetic as R², kept apart by name because chemists read it as the
cross-validated one: R² on out-of-fold predictions over the whole set.
"""

from __future__ import annotations

import math

import numpy as np

#: Reported for a regression, in the order the leaderboard shows them.
REGRESSION_METRICS = ("rmse", "mae", "r2", "q2", "spearman")
#: Reported for a classification.
CLASSIFICATION_METRICS = ("roc_auc", "pr_auc", "mcc", "balanced_accuracy", "calibration")
#: What the search minimizes or maximizes, per task.
OBJECTIVE = {"regression": "rmse", "classification": "roc_auc"}
#: The objective above two classes: ROC-AUC needs a one-against-rest average that MCC avoids.
MULTICLASS_OBJECTIVE = "mcc"
#: Metrics a higher value is better for; the rest are errors.
HIGHER_IS_BETTER = frozenset({"r2", "q2", "spearman", "roc_auc", "pr_auc", "mcc",
                              "balanced_accuracy"})  # fmt: skip

TITLES = {
    "rmse": "RMSE",
    "mae": "MAE",
    "r2": "R²",
    "q2": "Q²",
    "spearman": "Spearman ρ",
    "roc_auc": "ROC-AUC",
    "pr_auc": "PR-AUC",
    "mcc": "MCC",
    "balanced_accuracy": "Balanced accuracy",
    "calibration": "Calibration error",
}


def objective_of(task: str, classes: int = 2) -> str:
    """The metric a search for *task* optimizes."""
    if task == "classification" and classes > 2:
        return MULTICLASS_OBJECTIVE
    return OBJECTIVE[task]


def is_better(metric: str, value: float, than: float) -> bool:
    """Whether *value* beats *than* for *metric*. A NaN never beats a number."""
    if math.isnan(value):
        return False
    if math.isnan(than):
        return True
    return value > than if metric in HIGHER_IS_BETTER else value < than


def worst(metric: str) -> float:
    """The starting value a search can only improve on."""
    return -math.inf if metric in HIGHER_IS_BETTER else math.inf


def regression_scores(y: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    """Every regression metric from held-out predictions. Rows whose y or prediction is NaN
    are left out; fewer than two rows leaves every metric NaN."""
    from scipy import stats
    from sklearn import metrics as skm

    y = np.asarray(y, dtype=np.float64)
    predicted = np.asarray(predicted, dtype=np.float64)
    keep = np.isfinite(y) & np.isfinite(predicted)
    y, predicted = y[keep], predicted[keep]
    if len(y) < 2:
        return dict.fromkeys(REGRESSION_METRICS, math.nan)
    rmse = float(np.sqrt(np.mean((y - predicted) ** 2)))
    r2 = float(skm.r2_score(y, predicted)) if np.ptp(y) > 0 else math.nan
    rho = stats.spearmanr(y, predicted).statistic if np.ptp(predicted) > 0 else math.nan
    return {
        "rmse": rmse,
        "mae": float(np.mean(np.abs(y - predicted))),
        "r2": r2,
        "q2": r2,  # the same number: these predictions are already out of fold
        "spearman": float(rho),
    }


def classification_scores(
    y: np.ndarray, predicted: np.ndarray, scores: np.ndarray, classes: np.ndarray
) -> dict[str, float]:
    """Every classification metric from held-out predictions.

    *y* and *predicted* are class indices; *scores* is one column per class (a probability or a
    decision value). ROC-AUC and PR-AUC need two classes; above two they are the one-against-
    rest average, weighted by how many rows each class has.
    """
    from sklearn import metrics as skm

    y = np.asarray(y, dtype=np.int64)
    predicted = np.asarray(predicted, dtype=np.int64)
    scores = np.asarray(scores, dtype=np.float64).reshape(len(y), -1)
    out = dict.fromkeys(CLASSIFICATION_METRICS, math.nan)
    if len(y) < 2 or len(np.unique(y)) < 2:
        return out
    out["mcc"] = float(skm.matthews_corrcoef(y, predicted))
    out["balanced_accuracy"] = float(skm.balanced_accuracy_score(y, predicted))
    labels = np.arange(len(classes))
    if len(classes) == 2:
        positive = scores[:, -1]
        out["roc_auc"] = float(skm.roc_auc_score(y, positive))
        out["pr_auc"] = float(skm.average_precision_score(y, positive))
        out["calibration"] = _calibration(y, positive)
    else:
        present = np.isin(labels, np.unique(y))
        if present.all():
            out["roc_auc"] = float(
                skm.roc_auc_score(y, _normalized(scores), multi_class="ovr", average="weighted")
            )
            out["pr_auc"] = float(
                skm.average_precision_score(
                    np.eye(len(classes), dtype=np.int64)[y], scores, average="weighted"
                )
            )
    return out


def _normalized(scores: np.ndarray) -> np.ndarray:
    """Score columns as rows that sum to one, which is what ROC-AUC's one-against-rest
    average wants. Decision values are turned into weights by a softmax."""
    if np.all(scores >= 0) and np.allclose(scores.sum(axis=1), 1.0):
        return scores
    shifted = scores - scores.max(axis=1, keepdims=True)
    weights = np.exp(shifted)
    return weights / weights.sum(axis=1, keepdims=True)


def _calibration(y: np.ndarray, positive: np.ndarray, bins: int = 10) -> float:
    """Expected calibration error: how far predicted probability sits from the share that
    actually belong to the class, averaged over equal-width bins weighted by their size.

    NaN when the scores are not probabilities (a decision value has no calibration to check).
    """
    if positive.min() < 0.0 or positive.max() > 1.0:
        return math.nan
    edges = np.linspace(0.0, 1.0, bins + 1)
    which = np.clip(np.digitize(positive, edges[1:-1]), 0, bins - 1)
    error = 0.0
    for b in range(bins):
        rows = which == b
        n = int(rows.sum())
        if n:
            error += n * abs(float(positive[rows].mean()) - float(y[rows].mean()))
    return error / len(y)


def scores_for(
    task: str,
    y: np.ndarray,
    predicted: np.ndarray,
    *,
    scores: np.ndarray | None = None,
    classes: np.ndarray | None = None,
) -> dict[str, float]:
    """The metrics for whichever task, from held-out predictions."""
    if task == "regression":
        return regression_scores(y, predicted)
    if scores is None or classes is None:
        raise ValueError("a classification needs its class scores and labels to be scored")
    return classification_scores(y, predicted, scores, classes)
