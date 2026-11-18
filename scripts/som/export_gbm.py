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
#
"""A fitted ``HistGradientBoostingClassifier`` to the flat arrays the app evaluates.

docs/SOM_PLAN.md: the app must not depend on scikit-learn's pickle format, so the
ensemble ships as one ``.npz`` of tree arrays that ``chemsuite/metabolism/trees.py``
walks in numpy. Deliberately not pickle: a pickled estimator is opaque, tied to one
scikit-learn version, and a load-time hazard.

Used by ``train_gbm.py``; nothing here ships. The export is verified against
scikit-learn's own ``predict_proba`` before it is written, and the archive is written
only if they agree to :data:`TOLERANCE`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

#: Largest probability difference accepted between scikit-learn and the exported arrays.
#: The two evaluate the same float64 sums in a different order, so this is round-off, and
#: anything larger means the export misread the trees.
TOLERANCE = 1e-9


def tree_arrays(model: HistGradientBoostingClassifier) -> dict[str, np.ndarray]:
    """The ensemble's nodes, every tree concatenated, with ``tree_start`` to split them.

    One iteration of a binary ``HistGradientBoostingClassifier`` is one tree, and its
    leaves hold log-odds contributions that sum onto ``baseline``.
    """
    if model.n_trees_per_iteration_ != 1:
        raise ValueError("only the binary classifier is supported")

    trees = [iteration[0] for iteration in model._predictors]
    feature, threshold, left, right, is_leaf, value, missing_left = [], [], [], [], [], [], []
    starts = []
    offset = 0
    for tree in trees:
        nodes = tree.nodes
        if nodes["is_categorical"].any():
            raise ValueError("categorical splits are not supported; features are numeric")
        starts.append(offset)
        feature.append(nodes["feature_idx"].astype(np.int32))
        threshold.append(nodes["num_threshold"].astype(np.float64))
        # The child indices are per-tree, so they are shifted into the flat arrays here
        # rather than in the app.
        left.append(nodes["left"].astype(np.int64) + offset)
        right.append(nodes["right"].astype(np.int64) + offset)
        is_leaf.append(nodes["is_leaf"].astype(np.int8))
        value.append(nodes["value"].astype(np.float64))
        missing_left.append(nodes["missing_go_to_left"].astype(np.int8))
        offset += len(nodes)

    return {
        "tree_start": np.asarray(starts, dtype=np.int64),
        "node_feature": np.concatenate(feature),
        "node_threshold": np.concatenate(threshold),
        "node_left": np.concatenate(left).astype(np.int32),
        "node_right": np.concatenate(right).astype(np.int32),
        "node_leaf": np.concatenate(is_leaf),
        "node_value": np.concatenate(value),
        "node_missing_left": np.concatenate(missing_left),
        "baseline": np.asarray([float(np.ravel(model._baseline_prediction)[0])]),
    }


def write(
    path: Path,
    model: HistGradientBoostingClassifier,
    *,
    check_on: np.ndarray,
    extra: dict[str, np.ndarray],
) -> dict[str, float]:
    """Export *model* to *path*, with *extra* alongside, after checking it on *check_on*.

    Returns the check's largest differences. Raises if they exceed :data:`TOLERANCE`.
    """
    arrays = tree_arrays(model)

    from chemsuite.metabolism import trees

    ours = trees.probabilities(check_on, trees.Ensemble.from_arrays(arrays))
    theirs = model.predict_proba(check_on)[:, 1]
    difference = float(np.abs(ours - theirs).max())
    if difference > TOLERANCE:
        raise ValueError(
            f"the exported trees disagree with scikit-learn by {difference:.3g} "
            f"(tolerance {TOLERANCE:.0e}); the export is wrong, not the model"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays, **extra)
    return {
        "probability_difference": difference,
        "nodes": float(arrays["node_feature"].size),
        "trees": float(arrays["tree_start"].size),
        "megabytes": path.stat().st_size / 1e6,
    }
