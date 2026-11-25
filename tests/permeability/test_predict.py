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

"""The batch driver: what it does with a blank cell, a bad structure, and a cancel.

These tests fit their own small model rather than reading the shipped archives, so they
measure the driver and not the models — the shipped archives are checked by
``test_models.py``, and what they predict is pinned by ``test_golden.py``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.permeability import Endpoint
from chemsuite.permeability import predict as driver
from chemsuite.qsar.fit import fit
from chemsuite.qsar.recipe import Recipe

#: enough molecules to fit on, varied enough that the descriptors are not all the same number
SMILES = [
    "CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccccc1O", "c1ccccc1N", "CC(=O)O", "CC(=O)N",
    "CCN(CC)CC", "COc1ccccc1", "Clc1ccccc1", "Fc1ccccc1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "CC(=O)Nc1ccc(O)cc1", "OC(=O)c1ccccc1O",
    "CCCCCCCC", "CCCCCCCCO", "c1ccc2ccccc2c1", "C1CCCCC1", "C1CCNCC1", "CCOC(=O)C",
    "NCCc1ccc(O)cc1", "CSCC", "CC#N", "CCBr", "CCI", "OCC(O)CO", "NCCO", "CC(N)C(=O)O",
]  # fmt: skip

DESCRIPTORS = ("MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "RingCount")


def a_model(task: str = "regression"):
    """A small model whose recipe reads structures, so the driver can build its own pool."""
    from chemsuite.qsar import features

    recipe = Recipe(
        task=task,
        target="y",
        structures="Structure",
        descriptors=DESCRIPTORS,
        split="random",
        budget="quick",
        max_features=len(DESCRIPTORS),
        seed=0,
        # the driver does not care which estimator won, so the search is cut to the two
        # cheapest candidates per task: a boosted ensemble here would cost minutes
        without=("hgb", "rf", "mlp", "knn", "pls", "svr", "svc"),
    )
    pool = features.build(SMILES, recipe)
    rng = np.random.default_rng(0)
    raw = pool.x[:, 1] - 0.01 * pool.x[:, 2] + rng.normal(scale=0.1, size=pool.rows)
    classes = None
    y = raw
    if task == "classification":
        y = (raw > np.median(raw)).astype(np.float64)
        classes = np.array(["no", "yes"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted = fit(
            pool.x,
            pool.names,
            y,
            np.arange(pool.rows),
            recipe,
            keys=[f"KEY{i:04d}" for i in range(pool.rows)],
            classes=classes,
        )
    return fitted.model


@pytest.fixture(scope="module")
def regression_model():
    return a_model("regression")


@pytest.fixture(scope="module")
def classification_model():
    return a_model("classification")


def an_endpoint(task: str = "regression") -> Endpoint:
    return Endpoint("test", "Test endpoint", "Test", task, "units", "test.csqm", "a tooltip")


@pytest.fixture
def one_model(monkeypatch):
    """Hand the driver one model whatever endpoint it asks for."""

    def use(model):
        monkeypatch.setattr(driver, "shipped_model", lambda _endpoint: model)

    return use


def test_a_prediction_per_row_with_an_interval(one_model, regression_model):
    one_model(regression_model)
    endpoint = an_endpoint()
    results = driver.predict_permeability(SMILES[:8], [endpoint])

    assert len(results) == 8
    assert results.failed == 0
    said = results.per_endpoint["test"].as_prediction()
    assert np.isfinite(said.predicted).all()
    # an interval that contains its own prediction, and a domain distance for every row
    assert (said.low <= said.predicted).all()
    assert (said.predicted <= said.high).all()
    assert np.isfinite(said.distance).all()


def test_a_blank_cell_is_not_a_failure(one_model, regression_model):
    one_model(regression_model)
    results = driver.predict_permeability(["CCO", "", None, "   "], [an_endpoint()])

    said = results.per_endpoint["test"].as_prediction()
    assert results.failed == 0, "a row with no structure is not a structure that failed"
    assert np.isfinite(said.predicted[0])
    assert np.isnan(said.predicted[1:]).all()


def test_a_structure_that_cannot_be_read_is_counted_and_left_blank(one_model, regression_model):
    one_model(regression_model)
    results = driver.predict_permeability(["CCO", "not a structure", "CCCO"], [an_endpoint()])

    said = results.per_endpoint["test"].as_prediction()
    assert results.failed == 1
    assert np.isfinite(said.predicted[[0, 2]]).all()
    assert np.isnan(said.predicted[1])
    assert not said.inside[1], "a row with no prediction is not inside the domain"


def test_a_cancel_keeps_the_rows_already_finished(one_model, regression_model, monkeypatch):
    one_model(regression_model)
    monkeypatch.setattr(driver, "CHUNK_ROWS", 4)
    cells = SMILES[:16]

    class StopsAfterTwoChunks(CancelToken):
        def __init__(self) -> None:
            super().__init__()
            self.asked = 0

        @property
        def cancelled(self) -> bool:
            self.asked += 1
            return self.asked > 2

    with pytest.raises(Cancelled) as raised:
        driver.predict_permeability(cells, [an_endpoint()], cancel=StopsAfterTwoChunks())

    partial = raised.value.partial
    assert partial is not None
    assert len(partial) == 8, "the finished chunks are kept, in order"
    kept = partial.per_endpoint["test"].as_prediction()
    whole = driver.predict_permeability(cells[:8], [an_endpoint()])
    assert kept.predicted == pytest.approx(whole.per_endpoint["test"].as_prediction().predicted)


def test_a_classification_predicts_a_label_and_a_class_set(one_model, classification_model):
    one_model(classification_model)
    endpoint = an_endpoint("classification")
    results = driver.predict_permeability(["CCO", "", "CCCCO"], [endpoint])

    said = results.per_endpoint["test"].as_prediction()
    assert said.scores is not None and said.scores.shape == (3, 2)
    assert said.sets is not None and said.sets.shape == (3, 2)
    # a blank row gets an index no class name can be read from, which writes as ""
    assert said.predicted[1] == -1
    assert said.predicted[0] in (0, 1)
    assert not said.sets[1].any()


def test_every_endpoint_asked_for_gets_its_own_columns(one_model, regression_model):
    one_model(regression_model)
    other = Endpoint("other", "Other", "Other", "regression", "u", "o.csqm", "t")
    endpoints = [an_endpoint(), other]
    results = driver.predict_permeability(SMILES[:5], endpoints)

    assert sorted(results.per_endpoint) == ["other", "test"]
    assert all(len(r) == 5 for r in results.per_endpoint.values())
