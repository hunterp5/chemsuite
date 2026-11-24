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
"""The recipe: validated on the way in, and the same recipe after a JSON round trip."""

from __future__ import annotations

import pytest

from chemsuite.descriptors.catalog import CATALOG, PKA_KEYS, SHAPE_GROUP
from chemsuite.qsar.recipe import (
    BUDGETS,
    MAX_CLASSES,
    Recipe,
    default_descriptors,
    default_fingerprints,
    detect_task,
)


def a_recipe(**changes) -> Recipe:
    values = {
        "task": "regression",
        "target": "pIC50",
        "structures": "Structure",
        "descriptors": ("MolWt", "MolLogP"),
    }
    values.update(changes)
    return Recipe(**values)


def test_round_trip_keeps_every_setting():
    recipe = a_recipe(
        fingerprints=default_fingerprints(),
        columns=("Measured LogD",),
        split="column",
        split_column="Date",
        budget="thorough",
        confidence=0.95,
        seed=7,
        max_features=50,
        without=("mlp",),
    )
    assert Recipe.from_json(recipe.as_json()) == recipe


def test_json_is_plain_types():
    data = a_recipe().as_json()
    for value in data.values():
        assert isinstance(value, (str, int, float, list)), value


def test_a_setting_from_a_newer_chemsuite_is_refused_not_guessed():
    data = a_recipe().as_json()
    data["quantum_bogosity"] = 3
    with pytest.raises(ValueError, match="does not know"):
        Recipe.from_json(data)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"task": "divination"}, "unknown task"),
        ({"target": ""}, "needs a target"),
        ({"split": "vibes"}, "unknown split"),
        ({"split": "column"}, "column to order by"),
        ({"budget": "lavish"}, "unknown budget"),
        ({"confidence": 1.0}, "confidence"),
        ({"confidence": 0.1}, "confidence"),
        ({"correlation": 0.0}, "correlation"),
        ({"max_features": 0}, "at least one feature"),
        ({"descriptors": ("NoSuchDescriptor",)}, "unknown descriptor"),
        ({"fingerprints": ("morgan9",)}, "unknown fingerprint"),
        ({"descriptors": (), "columns": ()}, "at least one kind of feature"),
        ({"structures": ""}, "need a structure column"),
    ],
)
def test_nonsense_is_refused(changes, message):
    with pytest.raises(ValueError, match=message):
        a_recipe(**changes)


def test_columns_alone_need_no_structures():
    recipe = Recipe(task="regression", target="pIC50", descriptors=(), columns=("LogD",))
    assert not recipe.needs_structures


def test_default_pool_leaves_out_conformers_and_pka():
    keys = default_descriptors()
    assert keys, "the default pool is not empty"
    assert not PKA_KEYS.intersection(keys)
    assert not [k for k in keys if CATALOG[k].group == SHAPE_GROUP]
    assert not [k for k in keys if CATALOG[k].is_text]
    assert "MolWt" in keys and "TPSA" in keys


def test_default_fingerprints_count_features():
    from chemsuite.analysis.fingerprints import parse

    variants = [parse(key) for key in default_fingerprints()]
    assert variants and all(v.counts for v in variants)


@pytest.mark.parametrize(
    ("numeric", "distinct", "task"),
    [
        (False, 2, "classification"),
        (False, 40, "classification"),
        (True, 2, "classification"),
        (True, MAX_CLASSES, "classification"),
        (True, MAX_CLASSES + 1, "regression"),
        (True, 300, "regression"),
        (True, 1, "regression"),
    ],
)
def test_task_detection(numeric, distinct, task):
    assert detect_task(numeric, distinct) == task


def test_budgets_grow():
    order = ["quick", "standard", "thorough"]
    fits = [BUDGETS[key].fits for key in order]
    assert fits == sorted(fits)
    assert BUDGETS["quick"].candidates is not None
    assert BUDGETS["thorough"].candidates is None
