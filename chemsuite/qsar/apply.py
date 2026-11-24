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
"""Applying a saved model, in numpy and nothing else.

**This module must not import scikit-learn.** That is the whole point of it: a model the app
applies is arrays and arithmetic, so a saved model needs no pickle, survives a scikit-learn
upgrade, and cannot run code when it is opened. :mod:`chemsuite.qsar.export` writes the arrays
and refuses to write any it cannot reproduce here to within 1e-9.

A model is a *family* plus arrays plus a little JSON. The families, and the arithmetic each
one is:

| family | arrays | what it computes |
|---|---|---|
| ``linear`` | ``coef`` (outputs × features), ``intercept`` | ``x @ coefᵀ + intercept`` |
| ``svm`` | ``support``, ``dual``, ``intercept`` | an RBF kernel against the support vectors |
| ``forest`` | node arrays | the average of the trees' leaves |
| ``boosting`` | node arrays | a baseline plus every tree's leaf |
| ``knn`` | ``design``, ``target`` | the mean of the nearest neighbors, perhaps weighted |
| ``mlp`` | ``coef_0…``, ``intercept_0…`` | the layers, in order |

Ridge, the elastic net, logistic regression and partial least squares are all the ``linear``
family: PLS only centers its columns before applying its coefficients, so the centering folds
into the intercept and what is left is a straight line.

The *link* turns that output into what a person reads: ``identity`` leaves it alone,
``logistic`` and ``softmax`` turn log-odds into probability, and ``decision`` marks a support
vector machine's margin — a score that ranks rows without being a probability.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .neighbors import nearest, tile_rows
from .trees import Forest, logistic, softmax

LINEAR, SVM, FOREST, BOOSTING, KNN, MLP = "linear", "svm", "forest", "boosting", "knn", "mlp"
FAMILIES = (LINEAR, SVM, FOREST, BOOSTING, KNN, MLP)

IDENTITY, LOGISTIC, SOFTMAX, DECISION = "identity", "logistic", "softmax", "decision"
LINKS = (IDENTITY, LOGISTIC, SOFTMAX, DECISION)

#: Hidden-layer activations a neural network may use.
RELU, TANH = "relu", "tanh"


class ModelError(ValueError):
    """The arrays are not a model this version can apply."""


@dataclass(frozen=True)
class Weights:
    """A fitted model as arrays: what :mod:`chemsuite.qsar.export` writes and this applies."""

    family: str
    link: str
    arrays: dict[str, np.ndarray] = field(repr=False, default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.family not in FAMILIES:
            raise ModelError(f"unknown model family {self.family!r}")
        if self.link not in LINKS:
            raise ModelError(f"unknown link {self.link!r}")

    def array(self, name: str) -> np.ndarray:
        try:
            return np.asarray(self.arrays[name], dtype=np.float64)
        except KeyError:
            raise ModelError(f"a {self.family} model needs the {name!r} array") from None

    @property
    def classification(self) -> bool:
        return self.link in (LOGISTIC, SOFTMAX, DECISION) or bool(self.meta.get("classes"))


# -- the families ---------------------------------------------------------------------------
def _linear(w: Weights, x: np.ndarray) -> np.ndarray:
    return x @ w.array("coef").T + w.array("intercept")


def _svm(w: Weights, x: np.ndarray) -> np.ndarray:
    """One RBF kernel row per support vector: ``exp(-gamma · ‖x − v‖²)``, dual-weighted."""
    support = w.array("support")
    gamma = float(w.meta["gamma"])
    out = np.empty((len(x), 1))
    for start in range(0, len(x), 512):  # in tiles: rows × support vectors can be large
        stop = min(start + 512, len(x))
        square = (
            np.square(x[start:stop]).sum(axis=1)[:, None]
            - 2.0 * x[start:stop] @ support.T
            + np.square(support).sum(axis=1)[None, :]
        )
        kernel = np.exp(-gamma * np.maximum(square, 0.0))
        out[start:stop, 0] = kernel @ w.array("dual").reshape(-1)
    return out + w.array("intercept").reshape(-1)[0]


def _trees(w: Weights, x: np.ndarray) -> np.ndarray:
    from .trees import MEAN, SUM, raw

    return raw(x, Forest.from_arrays(w.arrays, MEAN if w.family == FOREST else SUM))


def _knn_weights(gaps: np.ndarray, kind: str) -> np.ndarray:
    """Neighbor weights. With inverse-distance weighting, a row that sits exactly on a
    training point takes all of the weight — scikit-learn's own rule, and the reason applying
    a model to its own training rows reproduces it instead of smoothing it."""
    if kind == "uniform":
        return np.ones_like(gaps)
    out = np.zeros_like(gaps)
    exact = gaps == 0.0
    coincide = exact.any(axis=1)
    with np.errstate(divide="ignore"):
        out[~coincide] = 1.0 / gaps[~coincide]
    out[coincide] = exact[coincide].astype(np.float64)
    return out


def _knn(w: Weights, x: np.ndarray) -> np.ndarray:
    design = w.array("design")
    target = w.array("target")
    k = int(w.meta["n_neighbors"])
    kind = str(w.meta["weights"])
    classes = int(w.meta.get("classes", 0))
    outputs = max(classes, 1)
    out = np.zeros((len(x), outputs))
    step = tile_rows(len(design), design.shape[1])
    for start in range(0, len(x), step):
        stop = min(start + step, len(x))
        near, chosen = nearest(x[start:stop], design, k)
        share = _knn_weights(chosen, kind)
        share /= share.sum(axis=1, keepdims=True)
        if classes:
            labels = target[near].astype(np.int64)
            for label in range(classes):
                out[start:stop, label] = (share * (labels == label)).sum(axis=1)
        else:
            out[start:stop, 0] = (share * target[near]).sum(axis=1)
    return out


def _mlp(w: Weights, x: np.ndarray) -> np.ndarray:
    activation = str(w.meta.get("activation", RELU))
    layers = int(w.meta["layers"])
    out = x
    for layer in range(layers):
        out = out @ w.array(f"coef_{layer}") + w.array(f"intercept_{layer}")
        if layer < layers - 1:
            out = np.maximum(out, 0.0) if activation == RELU else np.tanh(out)
    return out


_BY_FAMILY = {
    LINEAR: _linear,
    SVM: _svm,
    FOREST: _trees,
    BOOSTING: _trees,
    KNN: _knn,
    MLP: _mlp,
}


# -- what a caller asks for -----------------------------------------------------------------
def outputs(weights: Weights, x: np.ndarray) -> np.ndarray:
    """The model's raw output, as (rows, outputs), before the link."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 2:
        raise ModelError("features must be a two-dimensional array")
    out = np.asarray(_BY_FAMILY[weights.family](weights, x), dtype=np.float64)
    return out.reshape(len(x), -1)


def values(weights: Weights, x: np.ndarray) -> np.ndarray:
    """What the model says: numbers for a regression, one score column per class otherwise.

    The score columns are the same thing the search measured — a probability where the model
    has one, a margin for support vectors — so an interval or a class set calibrated during
    the fit means the same here.
    """
    raw_out = outputs(weights, x)
    if weights.link == IDENTITY and not weights.meta.get("classes"):
        return raw_out.reshape(len(x))
    if weights.link == LOGISTIC:
        probability = logistic(raw_out.reshape(len(x)))
        return np.column_stack([1.0 - probability, probability])
    if weights.link == SOFTMAX:
        return softmax(raw_out)
    if weights.link == DECISION:
        margin = raw_out.reshape(len(x))
        return np.column_stack([-margin, margin])
    return raw_out  # a forest's or a neighborhood's class shares are already probabilities


def predicted(weights: Weights, x: np.ndarray) -> np.ndarray:
    """A number per row, or the index of the class the model picks."""
    said = values(weights, x)
    if said.ndim == 1:
        return said
    return np.argmax(said, axis=1).astype(np.int64)
