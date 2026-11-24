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
"""A fitted scikit-learn estimator, written out as flat arrays.

Every candidate in the zoo has to come through here, because an estimator that cannot be
written out is an estimator the app could only keep as a pickle.

**Nothing is written that cannot be reproduced.** :func:`verify` applies the arrays with
:mod:`chemsuite.qsar.apply` and compares them against the estimator itself; past
:data:`TOLERANCE` the export is refused. The two sum the same float64 numbers in a different
order, so a disagreement of that size is round-off and anything larger is a misreading of the
model. ``scripts/som/export_gbm.py`` has held the site-of-metabolism trees to the same rule
since they were trained.

Reading another library's private attributes is the price of not pickling. They are named
explicitly here, never reached for with ``getattr``, so a scikit-learn release that moves one
breaks a test rather than quietly exporting nonsense.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from . import apply
from .apply import ModelError, Weights
from .search import Pipeline

#: The largest disagreement accepted between scikit-learn and the exported arrays.
TOLERANCE = 1e-9


def export(model: Pipeline) -> Weights:
    """*model*'s estimator as arrays. Raises :class:`~.apply.ModelError` for a family that
    cannot be written out."""
    family = model.family
    if family in ("linear", "pls"):
        return _linear(model)
    if family == "svm":
        return _svm(model)
    if family == "forest":
        return _forest(model)
    if family == "boosting":
        return _boosting(model)
    if family == "knn":
        return _knn(model)
    if family == "mlp":
        return _mlp(model)
    raise ModelError(f"no way to write out a {family} model")


def _classes(model: Pipeline) -> int:
    return 0 if model.classes is None else len(model.classes)


def _linear(model: Pipeline) -> Weights:
    """Ridge, the elastic net, logistic regression — and partial least squares, which only
    centers its columns first, so the centering folds into the intercept."""
    estimator = model.estimator
    coef = np.atleast_2d(np.asarray(estimator.coef_, dtype=np.float64))
    intercept = np.atleast_1d(np.asarray(estimator.intercept_, dtype=np.float64)).astype(np.float64)
    if model.family == "pls":
        intercept = intercept - np.asarray(estimator._x_mean, dtype=np.float64) @ coef.T
    classes = _classes(model)
    link = apply.IDENTITY
    if classes:
        link = apply.LOGISTIC if len(coef) == 1 else apply.SOFTMAX
    return Weights(
        apply.LINEAR,
        link,
        {"coef": coef, "intercept": intercept.reshape(len(coef))},
        {"classes": classes},
    )


def _svm(model: Pipeline) -> Weights:
    """Support vectors with an RBF kernel. Only two classes: scikit-learn's multiclass
    decision function is a vote between pairs, which is not worth reproducing by hand."""
    estimator = model.estimator
    classes = _classes(model)
    if classes > 2:
        raise ModelError("a support vector machine is written out for two classes only")
    return Weights(
        apply.SVM,
        apply.DECISION if classes else apply.IDENTITY,
        {
            "support": np.asarray(estimator.support_vectors_, dtype=np.float64),
            "dual": np.asarray(estimator.dual_coef_, dtype=np.float64).reshape(-1),
            "intercept": np.asarray(estimator.intercept_, dtype=np.float64).reshape(-1),
        },
        {"gamma": float(estimator._gamma), "classes": classes},
    )


def _tree_nodes(tree: Any, offset: int) -> dict[str, np.ndarray]:
    """One scikit-learn decision tree's nodes, its child indices shifted into flat arrays.

    The holes in the pool are filled before any estimator sees them, so no tree in the zoo
    carries a missing-value rule and ``missing_left`` is written as zeros.
    """
    leaf = (tree.children_left == -1).astype(np.int8)
    value = np.asarray(tree.value, dtype=np.float64)
    value = value.reshape(len(leaf), -1)
    if value.shape[1] > 1:  # a classifier's leaf: class shares, normalized across versions
        totals = value.sum(axis=1, keepdims=True)
        value = np.divide(value, totals, out=np.zeros_like(value), where=totals > 0)
    children = np.where(tree.children_left == -1, 0, tree.children_left + offset)
    right = np.where(tree.children_right == -1, 0, tree.children_right + offset)
    return {
        "node_feature": np.maximum(np.asarray(tree.feature, dtype=np.int32), 0),
        "node_threshold": np.asarray(tree.threshold, dtype=np.float64),
        "node_left": children.astype(np.int32),
        "node_right": right.astype(np.int32),
        "node_leaf": leaf,
        "node_value": value,
        "node_missing_left": np.zeros(len(leaf), dtype=np.int8),
    }


def _concatenate(parts: list[dict[str, np.ndarray]], starts: list[int], outputs: list[int]):
    joined = {name: np.concatenate([part[name] for part in parts]) for name in parts[0]}
    joined["tree_start"] = np.asarray(starts, dtype=np.int64)
    joined["tree_output"] = np.asarray(outputs, dtype=np.int32)
    return joined


def _forest(model: Pipeline) -> Weights:
    """A random forest: every tree votes, and the votes are averaged."""
    classes = _classes(model)
    width = max(classes, 1)
    parts: list[dict[str, np.ndarray]] = []
    starts: list[int] = []
    offset = 0
    for fitted in model.estimator.estimators_:
        nodes = _tree_nodes(fitted.tree_, offset)
        if nodes["node_value"].shape[1] != width:
            raise ModelError("a forest's leaf does not hold one value per output")
        starts.append(offset)
        parts.append(nodes)
        offset += len(nodes["node_leaf"])
    arrays = _concatenate(parts, starts, [0] * len(starts))
    arrays["baseline"] = np.zeros(width)
    return Weights(apply.FOREST, apply.IDENTITY, arrays, {"classes": classes})


def _boosting(model: Pipeline) -> Weights:
    """Histogram gradient boosting: a baseline, and every tree adds its leaf onto it.

    One iteration is one tree for a regression or two classes, and one tree per class above
    that — which is what ``tree_output`` records.
    """
    estimator = model.estimator
    classes = _classes(model)
    per_iteration = int(estimator.n_trees_per_iteration_)
    parts: list[dict[str, np.ndarray]] = []
    starts: list[int] = []
    outputs: list[int] = []
    offset = 0
    for iteration in estimator._predictors:
        for which, predictor in enumerate(iteration):
            nodes = predictor.nodes
            if nodes["is_categorical"].any():
                raise ModelError("categorical splits are not written out; features are numeric")
            starts.append(offset)
            outputs.append(which)
            parts.append(
                {
                    "node_feature": nodes["feature_idx"].astype(np.int32),
                    "node_threshold": nodes["num_threshold"].astype(np.float64),
                    "node_left": (nodes["left"].astype(np.int64) + offset).astype(np.int32),
                    "node_right": (nodes["right"].astype(np.int64) + offset).astype(np.int32),
                    "node_leaf": nodes["is_leaf"].astype(np.int8),
                    "node_value": nodes["value"].astype(np.float64).reshape(-1, 1),
                    "node_missing_left": nodes["missing_go_to_left"].astype(np.int8),
                }
            )
            offset += len(nodes)
    arrays = _concatenate(parts, starts, outputs)
    arrays["baseline"] = np.ravel(
        np.asarray(estimator._baseline_prediction, dtype=np.float64)
    ).copy()
    if len(arrays["baseline"]) != per_iteration:
        raise ModelError("the baseline does not match the trees per iteration")
    link = apply.IDENTITY
    if classes:
        link = apply.LOGISTIC if per_iteration == 1 else apply.SOFTMAX
    return Weights(apply.BOOSTING, link, arrays, {"classes": classes})


def _knn(model: Pipeline) -> Weights:
    """The training rows themselves: a neighborhood model is its data."""
    estimator = model.estimator
    classes = _classes(model)
    return Weights(
        apply.KNN,
        apply.IDENTITY,
        {
            "design": np.asarray(estimator._fit_X, dtype=np.float64),
            "target": np.asarray(estimator._y, dtype=np.float64),
        },
        {
            "n_neighbors": int(estimator.n_neighbors),
            "weights": str(estimator.weights),
            "classes": classes,
        },
    )


def _mlp(model: Pipeline) -> Weights:
    """The weight matrices and biases, layer by layer."""
    estimator = model.estimator
    classes = _classes(model)
    arrays: dict[str, np.ndarray] = {}
    for layer, (coef, bias) in enumerate(zip(estimator.coefs_, estimator.intercepts_, strict=True)):
        arrays[f"coef_{layer}"] = np.asarray(coef, dtype=np.float64)
        arrays[f"intercept_{layer}"] = np.asarray(bias, dtype=np.float64)
    out = str(estimator.out_activation_)
    link = {"identity": apply.IDENTITY, "logistic": apply.LOGISTIC, "softmax": apply.SOFTMAX}
    if out not in link:
        raise ModelError(f"a neural network with a {out!r} output is not written out")
    if str(estimator.activation) not in (apply.RELU, apply.TANH):
        raise ModelError(f"a {estimator.activation!r} activation is not written out")
    return Weights(
        apply.MLP,
        link[out],
        arrays,
        {
            "activation": str(estimator.activation),
            "layers": len(estimator.coefs_),
            "classes": classes,
        },
    )


def verify(model: Pipeline, weights: Weights, design: np.ndarray) -> float:
    """How far the arrays are from the estimator on *design*, and a refusal past the tolerance.

    *design* is the matrix the estimator was fitted on — the rows most likely to show up a
    misread tree or a dropped intercept.
    """
    theirs = _from_estimator(model, design)
    ours = apply.values(weights, design)
    if ours.shape != theirs.shape:
        raise ModelError(
            f"the exported arrays give {ours.shape} where scikit-learn gives {theirs.shape}"
        )
    difference = float(np.abs(ours - theirs).max()) if ours.size else 0.0
    if difference > TOLERANCE:
        raise ModelError(
            f"the exported {model.candidate} disagrees with scikit-learn by {difference:.3g} "
            f"(tolerance {TOLERANCE:.0e}); the export is wrong, not the model"
        )
    return difference


def _from_estimator(model: Pipeline, design: np.ndarray) -> np.ndarray:
    """What the estimator itself says, in the same shape :func:`~.apply.values` returns."""
    if model.task == "regression":
        return np.asarray(model.estimator.predict(design), dtype=np.float64).reshape(len(design))
    if model.family == "svm":
        margin = np.asarray(model.estimator.decision_function(design), dtype=np.float64).reshape(
            len(design)
        )
        return np.column_stack([-margin, margin])
    return np.asarray(model.estimator.predict_proba(design), dtype=np.float64)


def exported(model: Pipeline, design: np.ndarray) -> tuple[Weights, float]:
    """Write *model* out and check it in one step: the arrays, and how far off they were."""
    weights = export(model)
    return weights, verify(model, weights, design)
