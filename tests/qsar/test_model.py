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
"""The saved model.

Three tests carry the module. A model **round trips** — saved, reopened, and predicting the
same numbers. Reopening it **never unpickles**, checked by handing the loader a file with a
pickled object in it. And a fresh interpreter can apply it **without scikit-learn importable
at all**, which is the property that makes the archive worth its complexity.
"""

from __future__ import annotations

import io
import json
import math
import pickle
import subprocess
import sys
import warnings

import numpy as np
import pytest

from chemsuite.qsar import model as model_module
from chemsuite.qsar.apply import ModelError
from chemsuite.qsar.fit import fit
from chemsuite.qsar.model import SUFFIX, Model, Training, load, unpack
from chemsuite.qsar.recipe import Recipe


def a_recipe(task: str = "regression", **changes) -> Recipe:
    values = {
        "task": task,
        "target": "y",
        "columns": tuple(f"c{i}" for i in range(6)),
        "split": "random",
        "budget": "quick",
        "max_features": 6,
        "seed": 0,
        "without": ("rf", "mlp", "knn"),
    }
    values.update(changes)
    return Recipe(**values)


def a_fit(task: str = "regression", rows: int = 140, **changes):
    rng = np.random.default_rng(0)
    pool = rng.normal(size=(rows, 12))
    raw = 2.0 * pool[:, 0] - pool[:, 1] + rng.normal(scale=0.2, size=rows)
    names = [f"c:col{i}" for i in range(pool.shape[1])]
    keys = [f"KEY{i:04d}" for i in range(rows)]
    classes = None
    y = raw
    if task == "classification":
        y = (raw > np.median(raw)).astype(np.int64)
        classes = np.array(["low", "high"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted = fit(
            pool, names, y, np.arange(rows), a_recipe(task, **changes),
            keys=keys, classes=classes,
        )  # fmt: skip
    return fitted, pool


# -- the three that matter ------------------------------------------------------------------
def test_a_saved_model_predicts_the_same_numbers_when_reopened(tmp_path):
    fitted, pool = a_fit()
    path = fitted.model.save(tmp_path / "model")
    assert path.suffix == SUFFIX
    before = fitted.model.predict(pool)
    after = load(path).predict(pool)
    assert np.allclose(before.predicted, after.predicted)
    assert np.allclose(before.low, after.low)
    assert np.allclose(before.high, after.high)
    assert np.allclose(before.distance, after.distance)


def test_reopening_a_model_never_unpickles(tmp_path):
    """A model file can be mailed around, so the loader must not be a way to run code."""
    buffer = io.BytesIO()
    np.savez(buffer, _format=np.array(1, dtype=np.int32), evil=np.array(pickle.dumps(len)))
    with pytest.raises(ModelError):
        unpack(buffer.getvalue())

    fitted, _pool = a_fit()
    path = fitted.model.save(tmp_path / "model.csqm")
    assert b"pickle" not in path.read_bytes().lower()


def test_a_fresh_interpreter_applies_it_without_scikit_learn(tmp_path):
    """The point of the arrays: no scikit-learn on the apply path, and nothing imports it."""
    fitted, pool = a_fit()
    path = fitted.model.save(tmp_path / "model.csqm")
    np.save(tmp_path / "pool.npy", pool)
    expected = fitted.model.predict(pool).predicted
    script = tmp_path / "apply_it.py"
    script.write_text(
        "import sys\n"
        "class Blocked:\n"
        "    def find_module(self, name, path=None):\n"
        "        if name == 'sklearn' or name.startswith('sklearn.'):\n"
        "            raise ImportError('scikit-learn must not be needed to apply a model')\n"
        "sys.meta_path.insert(0, Blocked())\n"
        "import numpy as np\n"
        "from chemsuite.qsar.model import load\n"
        f"model = load(r'{path}')\n"
        f"pool = np.load(r'{tmp_path / 'pool.npy'}')\n"
        "np.save(r'" + str(tmp_path / "said.npy") + "', model.predict(pool).predicted)\n"
        "assert 'sklearn' not in sys.modules, sorted(m for m in sys.modules if 'sk' in m)\n",
        encoding="utf-8",
    )
    done = subprocess.run(
        [sys.executable, "-I", str(script)], capture_output=True, text=True, timeout=300
    )
    assert done.returncode == 0, done.stderr
    assert np.allclose(np.load(tmp_path / "said.npy"), expected)


# -- what the archive carries ---------------------------------------------------------------
def test_the_model_knows_what_it_was_trained_on():
    fitted, _pool = a_fit(rows=140)
    model = fitted.model
    assert len(model.training) == 140
    assert model.training.keys[0] == "KEY0000"
    assert model.recipe == a_recipe()
    assert model.features == len(model.kept) <= 6
    assert model.kept_names()[0].startswith("c:col")


def test_the_scores_are_the_out_of_fold_ones():
    fitted, _pool = a_fit()
    assert fitted.model.scores == fitted.outcome.scores
    assert fitted.model.scores["r2"] > 0.8
    assert fitted.model.objective == "rmse"


def test_the_leaderboard_and_the_lineage_are_saved(tmp_path):
    fitted, _pool = a_fit()
    back = load(fitted.model.save(tmp_path / "m.csqm"))
    assert len(back.leaderboard) == len(fitted.outcome.leaderboard)
    assert {entry["candidate"] for entry in back.leaderboard}
    assert len(back.lineage) == 1
    assert back.lineage[0]["note"] == "fitted"
    assert back.lineage[0]["rows"] == len(back.training)


def test_a_refit_adds_a_line_to_the_lineage():
    fitted, _pool = a_fit()
    again = fitted.model.with_lineage("refit", {"rmse": 0.3}, rows=200)
    assert len(again.lineage) == 2
    assert again.lineage[-1]["note"] == "refit"
    assert len(fitted.model.lineage) == 1, "the model itself did not change"


def test_a_not_a_number_score_survives_the_manifest(tmp_path):
    fitted, _pool = a_fit()
    broken = model_module.Model(
        **{**fitted.model.__dict__, "scores": {"rmse": 0.5, "spearman": math.nan}}
    )
    back = load(broken.save(tmp_path / "m.csqm"))
    assert back.scores["rmse"] == 0.5
    assert math.isnan(back.scores["spearman"]), "JSON has no NaN, so it travels as text"


# -- classification -------------------------------------------------------------------------
def test_a_classification_model_round_trips_with_its_class_sets(tmp_path):
    fitted, pool = a_fit("classification")
    back = load(fitted.model.save(tmp_path / "m.csqm"))
    before = fitted.model.predict(pool)
    after = back.predict(pool)
    assert np.allclose(before.scores, after.scores)
    assert before.sets.tolist() == after.sets.tolist()
    assert back.classes is not None and back.classes.tolist() == ["low", "high"]
    names = after.set_names(back.classes)
    assert len(names) == len(pool)
    assert any(name in ("low", "high") for name in names)


def test_a_prediction_set_can_be_empty_or_both():
    fitted, pool = a_fit("classification")
    said = fitted.model.predict(pool)
    sizes = said.sets.sum(axis=1)
    assert sizes.min() >= 0 and sizes.max() <= 2


def test_asking_a_regression_for_class_names_is_refused():
    fitted, pool = a_fit()
    with pytest.raises(ValueError, match="no class sets"):
        fitted.model.predict(pool).set_names(np.array(["a"]))


# -- refusals -------------------------------------------------------------------------------
def test_a_pool_of_the_wrong_width_is_refused_with_advice():
    fitted, pool = a_fit()
    with pytest.raises(ModelError, match="Rebuild it from the model's own recipe"):
        fitted.model.predict(pool[:, :4])


def test_rubbish_is_not_a_model():
    with pytest.raises(ModelError, match="not a chemsuite model"):
        unpack(b"this is not an npz file at all")


def test_a_model_without_a_format_number_is_refused():
    buffer = io.BytesIO()
    np.savez(buffer, kept=np.arange(3))
    with pytest.raises(ModelError, match="no format number"):
        unpack(buffer.getvalue())


def test_a_model_from_a_newer_chemsuite_says_so():
    buffer = io.BytesIO()
    np.savez(buffer, _format=np.array(99, dtype=np.int32))
    with pytest.raises(ModelError, match="newer chemsuite"):
        unpack(buffer.getvalue())


def test_a_model_without_a_manifest_is_refused():
    buffer = io.BytesIO()
    np.savez(buffer, _format=np.array(1, dtype=np.int32), kept=np.arange(3))
    with pytest.raises(ModelError, match="no manifest"):
        unpack(buffer.getvalue())


def test_a_manifest_without_its_arrays_is_refused():
    buffer = io.BytesIO()
    meta = json.dumps({"format": 1, "family": "linear", "link": "identity"}).encode()
    np.savez(
        buffer,
        _format=np.array(1, dtype=np.int32),
        _meta=np.frombuffer(meta, dtype=np.uint8),
    )
    with pytest.raises(ModelError, match="missing"):
        unpack(buffer.getvalue())


def test_a_training_row_needs_both_a_key_and_a_value():
    with pytest.raises(ValueError, match="key and a value"):
        Training(np.array(["a", "b"]), np.array([1.0]))


def test_a_model_is_named_for_what_it_predicts():
    fitted, _pool = a_fit()
    assert fitted.model.title.startswith("y (")
    named = Model(**{**fitted.model.__dict__, "label": "Solubility, March"})
    assert named.title == "Solubility, March"


def test_a_model_can_be_kept_as_a_document_asset(tmp_path):
    """The document stores a model as an asset, so it must satisfy core.assets.Packable:
    able to turn itself back into bytes, and able to say how much it is holding."""
    from chemsuite.core.assets import Packable

    fitted, pool = a_fit()
    model = fitted.model
    assert isinstance(model, Packable)
    assert model.nbytes > 0
    back = unpack(model.to_bytes())
    assert np.allclose(back.predict(pool).predicted, model.predict(pool).predicted)
