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

"""The exported ensemble, walked in numpy.

The test that matters is the last one: a model fitted here with scikit-learn, exported by
``scripts/som/export_gbm.py`` and evaluated by :mod:`chemsuite.metabolism.trees` must give
scikit-learn's own probabilities back. That is what lets the app ship trees instead of a
pickle.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from chemsuite.metabolism import trees

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "som"


def one_split(threshold: float = 0.5, left: float = -1.0, right: float = 2.0) -> dict:
    """One tree: a root splitting on column 0, and two leaves."""
    return {
        "tree_start": np.array([0]),
        "node_feature": np.array([0, 0, 0], dtype=np.int32),
        "node_threshold": np.array([threshold, 0.0, 0.0]),
        "node_left": np.array([1, 0, 0], dtype=np.int32),
        "node_right": np.array([2, 0, 0], dtype=np.int32),
        "node_leaf": np.array([0, 1, 1], dtype=np.int8),
        "node_value": np.array([0.0, left, right]),
        "node_missing_left": np.array([1, 0, 0], dtype=np.int8),
        "baseline": np.array([0.25]),
    }


def test_a_row_goes_left_when_it_is_at_or_below_the_threshold():
    model = trees.Ensemble.from_arrays(one_split())
    scores = trees.raw_scores(np.array([[0.0], [0.5], [0.5001]]), model)
    assert scores == pytest.approx([0.25 - 1.0, 0.25 - 1.0, 0.25 + 2.0])


def test_a_missing_value_follows_the_nodes_own_flag():
    arrays = one_split()
    model = trees.Ensemble.from_arrays(arrays)
    assert trees.raw_scores(np.array([[np.nan]]), model) == pytest.approx([0.25 - 1.0])
    arrays["node_missing_left"] = np.array([0, 0, 0], dtype=np.int8)
    model = trees.Ensemble.from_arrays(arrays)
    assert trees.raw_scores(np.array([[np.nan]]), model) == pytest.approx([0.25 + 2.0])


def test_every_tree_adds_to_the_baseline():
    arrays = one_split()
    # the same tree twice, concatenated as the exporter concatenates them
    doubled = dict(arrays)
    doubled["tree_start"] = np.array([0, 3])
    for name in ("node_feature", "node_threshold", "node_leaf", "node_value",
                 "node_missing_left"):  # fmt: skip
        doubled[name] = np.concatenate([arrays[name], arrays[name]])
    doubled["node_left"] = np.array([1, 0, 0, 4, 3, 3], dtype=np.int32)
    doubled["node_right"] = np.array([2, 0, 0, 5, 3, 3], dtype=np.int32)
    model = trees.Ensemble.from_arrays(doubled)
    assert model.trees == 2
    assert trees.raw_scores(np.array([[1.0]]), model) == pytest.approx([0.25 + 4.0])


def test_probabilities_are_a_logistic_and_survive_a_large_score():
    model = trees.Ensemble.from_arrays(one_split(left=-800.0, right=800.0))
    found = trees.probabilities(np.array([[0.0], [1.0]]), model)
    assert found[0] == pytest.approx(0.0, abs=1e-300)
    assert found[1] == pytest.approx(1.0)
    assert np.isfinite(found).all()


def test_an_archive_missing_an_array_says_which():
    arrays = one_split()
    del arrays["node_value"]
    with pytest.raises(ValueError, match="node_value"):
        trees.Ensemble.from_arrays(arrays)


def test_no_rows_is_an_empty_answer_not_a_failure():
    model = trees.Ensemble.from_arrays(one_split())
    assert trees.probabilities(np.zeros((0, 1)), model).shape == (0,)


def test_the_export_reproduces_scikit_learns_own_probabilities(tmp_path):
    sklearn_ensemble = pytest.importorskip("sklearn.ensemble")
    sys.path.insert(0, str(SCRIPTS))
    try:
        import export_gbm
    finally:
        sys.path.remove(str(SCRIPTS))

    rng = np.random.default_rng(0)
    x = rng.normal(size=(400, 6))
    y = (x[:, 0] + 0.5 * x[:, 2] - x[:, 4] ** 2 > 0).astype(int)
    model = sklearn_ensemble.HistGradientBoostingClassifier(
        max_iter=40, early_stopping=False, max_leaf_nodes=15, random_state=0
    ).fit(x, y)

    path = tmp_path / "check.npz"
    written = export_gbm.write(path, model, check_on=x, extra={})
    assert written["trees"] == 40
    assert written["probability_difference"] < export_gbm.TOLERANCE

    # and again from the written archive, which is how the app reads it
    loaded = trees.Ensemble.load(path)
    ours = trees.probabilities(x, loaded)
    assert ours == pytest.approx(model.predict_proba(x)[:, 1], abs=1e-12)
