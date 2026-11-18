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

"""A gradient boosted ensemble evaluated in numpy, from the arrays ``scripts/som`` exports.

The app does not depend on scikit-learn's pickle format: a fitted
``HistGradientBoostingClassifier`` is exported as flat node arrays (every tree
concatenated, ``tree_start`` marking each one's first node) and walked here. One tree is
one boosting iteration; its leaf values are log-odds contributions that sum onto the
baseline, and a logistic turns the total into a probability.

The split rule is scikit-learn's own: a row goes **left** when its feature value is less
than or equal to the node's threshold, and a missing value follows the node's
``missing_go_to_left`` flag. ``scripts/som/export_gbm.py`` refuses to write an archive
whose predictions here differ from scikit-learn's by more than 1e-9, so this file and
that one are checked against each other every time a model is trained.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: The array names an ensemble needs; an archive may carry more (the isoform list, the
#: feature names, the applicability-domain index).
REQUIRED = (
    "tree_start",
    "node_feature",
    "node_threshold",
    "node_left",
    "node_right",
    "node_leaf",
    "node_value",
    "node_missing_left",
    "baseline",
)


@dataclass(frozen=True)
class Ensemble:
    """Every tree of one boosted model, concatenated into flat node arrays."""

    tree_start: np.ndarray  # int64, first node of each tree
    feature: np.ndarray  # int32, the column a node splits on
    threshold: np.ndarray  # float64
    left: np.ndarray  # int32, index into these same arrays
    right: np.ndarray  # int32
    leaf: np.ndarray  # int8, 1 at a leaf
    value: np.ndarray  # float64, a leaf's log-odds contribution
    missing_left: np.ndarray  # int8, where a NaN goes
    baseline: float  # the log-odds every tree adds to

    @classmethod
    def from_arrays(cls, arrays: dict[str, np.ndarray]) -> Ensemble:
        """An ensemble from a mapping of the exported array names."""
        missing = [name for name in REQUIRED if name not in arrays]
        if missing:
            raise ValueError(f"the model archive is missing {', '.join(missing)}")
        return cls(
            tree_start=np.asarray(arrays["tree_start"], dtype=np.int64),
            feature=np.asarray(arrays["node_feature"], dtype=np.int32),
            threshold=np.asarray(arrays["node_threshold"], dtype=np.float64),
            left=np.asarray(arrays["node_left"], dtype=np.int32),
            right=np.asarray(arrays["node_right"], dtype=np.int32),
            leaf=np.asarray(arrays["node_leaf"], dtype=np.int8),
            value=np.asarray(arrays["node_value"], dtype=np.float64),
            missing_left=np.asarray(arrays["node_missing_left"], dtype=np.int8),
            baseline=float(np.ravel(arrays["baseline"])[0]),
        )

    @classmethod
    def load(cls, path: Path) -> Ensemble:
        """The ensemble in the ``.npz`` archive at *path*."""
        with np.load(path) as archive:
            return cls.from_arrays({name: archive[name] for name in REQUIRED})

    @property
    def trees(self) -> int:
        return int(self.tree_start.size)


def raw_scores(features: np.ndarray, model: Ensemble) -> np.ndarray:
    """The log-odds of each row of *features* (one row per atom and isoform)."""
    rows = np.asarray(features, dtype=np.float64)
    if rows.ndim != 2:
        raise ValueError("features must be a two-dimensional array")
    total = np.full(rows.shape[0], model.baseline, dtype=np.float64)
    if rows.shape[0] == 0:
        return total

    for start in model.tree_start:
        # Every row starts at the tree's root and walks down together; the loop ends
        # when they have all reached a leaf, so it runs the tree's depth, not its size.
        node = np.full(rows.shape[0], start, dtype=np.int64)
        walking = model.leaf[node] == 0
        while walking.any():
            here = node[walking]
            values = rows[walking, model.feature[here]]
            go_left = np.where(
                np.isnan(values),
                model.missing_left[here] == 1,
                values <= model.threshold[here],
            )
            node[walking] = np.where(go_left, model.left[here], model.right[here])
            walking = model.leaf[node] == 0
        total += model.value[node]
    return total


def probabilities(features: np.ndarray, model: Ensemble) -> np.ndarray:
    """:func:`raw_scores` through a logistic, written so it cannot overflow."""
    scores = raw_scores(features, model)
    return np.where(
        scores >= 0,
        1.0 / (1.0 + np.exp(-np.abs(scores))),
        np.exp(-np.abs(scores)) / (1.0 + np.exp(-np.abs(scores))),
    )
