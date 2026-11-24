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
"""An ensemble of trees as flat arrays, walked in numpy.

This is how a QSAR model travels without a pickle. A fitted forest or boosted ensemble is
written out as node arrays — every tree concatenated, ``tree_start`` marking each one's first
node — and :func:`raw` walks them. :mod:`chemsuite.qsar.export` refuses to write an archive
whose predictions here differ from scikit-learn's by more than 1e-9, so this file is not an
approximation of the estimator; it *is* the estimator.

Two shapes of ensemble, told apart by :attr:`Forest.combine`:

- **boosting** (``"sum"``): every tree adds its leaf onto a baseline. One tree holds one
  value per leaf, and ``tree_output`` says which output column it contributes to — a
  multiclass ``HistGradientBoosting`` grows one tree per class per iteration.
- **a forest** (``"mean"``): every tree votes and the votes are averaged. A leaf holds one
  value per output at once, so a classifier's leaf carries its whole class distribution.

The split rule is scikit-learn's own: a row goes **left** when its value is less than or equal
to the node's threshold, and a missing value follows the node's ``missing_left`` flag.

A near relation lives in :mod:`chemsuite.metabolism.trees`, which walks the sites-of-metabolism
model. That one stays specialized — a binary classifier with a scalar baseline and a logistic —
because its arrays are pinned by a shipped model and its golden test; the node descent is
written twice on purpose.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: The arrays an ensemble cannot be rebuilt without.
REQUIRED = (
    "tree_start",
    "tree_output",
    "node_feature",
    "node_threshold",
    "node_left",
    "node_right",
    "node_leaf",
    "node_value",
    "node_missing_left",
    "baseline",
)
#: How a tree's leaves are put together.
SUM, MEAN = "sum", "mean"


@dataclass(frozen=True)
class Forest:
    """Every tree of one ensemble, concatenated into flat node arrays."""

    tree_start: np.ndarray  # int64, the first node of each tree
    tree_output: np.ndarray  # int32, the output column each tree adds to
    feature: np.ndarray  # int32, the column a node splits on
    threshold: np.ndarray  # float64
    left: np.ndarray  # int32, an index into these same arrays
    right: np.ndarray  # int32
    leaf: np.ndarray  # int8, 1 at a leaf
    value: np.ndarray  # float64, (nodes, width): what a leaf is worth
    missing_left: np.ndarray  # int8, where a NaN goes
    baseline: np.ndarray  # float64, (outputs,): what every tree adds to
    combine: str = SUM

    def __post_init__(self) -> None:
        if self.combine not in (SUM, MEAN):
            raise ValueError(f"an ensemble is combined by {SUM} or {MEAN}, not {self.combine!r}")
        if self.value.ndim != 2:
            raise ValueError("node values must be (nodes, width)")
        if self.combine == MEAN and self.width != self.outputs:
            raise ValueError("a forest's leaf holds one value per output")
        if self.combine == SUM and self.width != 1:
            raise ValueError("a boosted tree's leaf holds one value")
        if len(self.tree_output) != len(self.tree_start):
            raise ValueError("every tree needs an output column")

    @property
    def trees(self) -> int:
        return int(self.tree_start.size)

    @property
    def outputs(self) -> int:
        return int(self.baseline.size)

    @property
    def width(self) -> int:
        return int(self.value.shape[1])

    def leaves(self, rows: np.ndarray, start: int) -> np.ndarray:
        """The leaf each row of *rows* reaches in the tree beginning at *start*.

        Every row walks down together, so the loop runs the tree's depth, not its size.
        """
        node = np.full(len(rows), start, dtype=np.int64)
        walking = self.leaf[node] == 0
        while walking.any():
            here = node[walking]
            values = rows[walking, self.feature[here]]
            go_left = np.where(
                np.isnan(values),
                self.missing_left[here] == 1,
                values <= self.threshold[here],
            )
            node[walking] = np.where(go_left, self.left[here], self.right[here])
            walking = self.leaf[node] == 0
        return node

    def as_arrays(self) -> dict[str, np.ndarray]:
        return {
            "tree_start": self.tree_start,
            "tree_output": self.tree_output,
            "node_feature": self.feature,
            "node_threshold": self.threshold,
            "node_left": self.left,
            "node_right": self.right,
            "node_leaf": self.leaf,
            "node_value": self.value,
            "node_missing_left": self.missing_left,
            "baseline": self.baseline,
        }

    @classmethod
    def from_arrays(cls, arrays: dict[str, np.ndarray], combine: str = SUM) -> Forest:
        missing = [name for name in REQUIRED if name not in arrays]
        if missing:
            raise ValueError(f"the ensemble is missing {', '.join(missing)}")
        value = np.asarray(arrays["node_value"], dtype=np.float64)
        return cls(
            tree_start=np.asarray(arrays["tree_start"], dtype=np.int64),
            tree_output=np.asarray(arrays["tree_output"], dtype=np.int32),
            feature=np.asarray(arrays["node_feature"], dtype=np.int32),
            threshold=np.asarray(arrays["node_threshold"], dtype=np.float64),
            left=np.asarray(arrays["node_left"], dtype=np.int32),
            right=np.asarray(arrays["node_right"], dtype=np.int32),
            leaf=np.asarray(arrays["node_leaf"], dtype=np.int8),
            value=value.reshape(len(value), -1),
            missing_left=np.asarray(arrays["node_missing_left"], dtype=np.int8),
            baseline=np.atleast_1d(np.asarray(arrays["baseline"], dtype=np.float64)),
            combine=combine,
        )


def raw(x: np.ndarray, model: Forest) -> np.ndarray:
    """The ensemble's output for every row of *x*, as (rows, outputs).

    Before any logistic or softmax: the sum of what each tree contributes for boosting, or the
    average of what each tree holds for a forest.
    """
    rows = np.asarray(x, dtype=np.float64)
    if rows.ndim != 2:
        raise ValueError("features must be a two-dimensional array")
    total = np.tile(model.baseline, (len(rows), 1))
    if not len(rows) or not model.trees:
        return total
    for tree, start in enumerate(model.tree_start.tolist()):
        reached = model.leaves(rows, start)
        if model.width == 1:
            total[:, int(model.tree_output[tree])] += model.value[reached, 0]
        else:
            total += model.value[reached]
    if model.combine == MEAN:
        total /= model.trees
    return total


def logistic(scores: np.ndarray) -> np.ndarray:
    """A logistic that cannot overflow, for a boosted binary classifier's log-odds."""
    positive = scores >= 0
    small = np.exp(-np.abs(scores))
    return np.where(positive, 1.0 / (1.0 + small), small / (1.0 + small))


def softmax(scores: np.ndarray) -> np.ndarray:
    """Rows of log-odds as rows of probability, for a boosted multiclass classifier."""
    shifted = scores - scores.max(axis=1, keepdims=True)
    weights = np.exp(shifted)
    return weights / weights.sum(axis=1, keepdims=True)
