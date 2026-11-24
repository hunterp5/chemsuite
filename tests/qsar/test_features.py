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
"""The feature pool.

What matters here is that the columns are the same numbers the table's own tools compute, and
that they land in the order the recipe implies — a saved model rebuilds its pool from the
recipe alone, so a column that moves is a model that silently predicts nonsense.
"""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.analysis import fingerprints as fp
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.descriptors.compute import compute_descriptors
from chemsuite.qsar import features
from chemsuite.qsar.recipe import Recipe

SMILES = ["CCO", "c1ccccc1", "CC(=O)O", "CCN(CC)CC", "C1CCCCC1O"]
COUNTS = fp.variant("morgan2", counts=True)


def a_recipe(**changes) -> Recipe:
    values = {
        "task": "regression",
        "target": "y",
        "structures": "Structure",
        "descriptors": ("MolWt", "MolLogP", "TPSA"),
    }
    values.update(changes)
    return Recipe(**values)


def test_the_pool_is_as_wide_as_the_recipe_says_before_anything_is_computed():
    recipe = a_recipe(fingerprints=(COUNTS,), columns=("LogD",))
    expected = 3 + fp.kind_of(COUNTS).bits + 1
    assert features.width(recipe) == expected
    names, kinds = features.pool_names(recipe)
    assert len(names) == len(kinds) == expected
    pool = features.build(SMILES, recipe, columns={"LogD": np.arange(5.0)})
    assert pool.columns == expected
    assert pool.names == names and pool.kinds == kinds


def test_the_columns_come_in_the_recipe_order():
    recipe = a_recipe(fingerprints=(COUNTS,), columns=("LogD",))
    pool = features.build(SMILES, recipe, columns={"LogD": np.arange(5.0)})
    assert pool.kinds[0] == features.DESCRIPTOR
    assert pool.kinds[3] == features.FINGERPRINT
    assert pool.kinds[-1] == features.COLUMN
    assert pool.names[0] == "d:MolWt"
    assert pool.names[-1] == "c:LogD"


def test_the_descriptors_are_the_table_tools_own_numbers():
    recipe = a_recipe()
    pool = features.build(SMILES, recipe)
    computed = compute_descriptors(SMILES, recipe.descriptors)
    for i, key in enumerate(recipe.descriptors):
        assert np.allclose(pool.x[:, i], computed[key].values(), equal_nan=True)


def test_the_fingerprint_block_is_the_similarity_tools_own_counts():
    recipe = a_recipe(descriptors=(), fingerprints=(COUNTS,))
    pool = features.build(SMILES, recipe)
    counts, valid = fp.compute_counts(SMILES, COUNTS)
    assert np.allclose(pool.x, counts.astype(float))
    assert pool.usable.tolist() == valid.tolist()


def test_fingerprint_bits_are_zero_or_one():
    bits = fp.variant("morgan2", bits=512)
    pool = features.build(SMILES, a_recipe(descriptors=(), fingerprints=(bits,)))
    assert set(np.unique(pool.x).tolist()) <= {0.0, 1.0}
    assert pool.columns == 512


def test_a_table_column_is_carried_through_blanks_and_all():
    values = np.array([1.0, np.nan, 3.0, 4.0, 5.0])
    pool = features.build(SMILES, a_recipe(descriptors=(), columns=("LogD",),
                                           structures=""), columns={"LogD": values})  # fmt: skip
    assert np.allclose(pool.x[:, 0], values, equal_nan=True)
    assert pool.usable.all(), "a blank number is a hole to fill, not an unusable row"


def test_a_structure_that_cannot_be_read_is_marked_unusable():
    cells = [*SMILES, "not a molecule at all", None]
    pool = features.build(cells, a_recipe(fingerprints=(COUNTS,)))
    assert pool.usable[:5].all()
    assert not pool.usable[5] and not pool.usable[6]


def test_with_descriptors_alone_a_row_with_no_value_at_all_is_unusable():
    cells = [*SMILES, "not a molecule at all"]
    pool = features.build(cells, a_recipe())
    assert pool.usable[:5].all()
    assert not pool.usable[5]


def test_a_column_the_caller_forgot_is_an_error_not_a_hole():
    with pytest.raises(KeyError, match="did not pass"):
        features.build(SMILES, a_recipe(columns=("LogD",)))


def test_titles_read_as_a_person_would_say_them():
    assert features.title_of("d:MolWt") == "MW"  # the catalog's own title
    assert features.title_of(f"f:{COUNTS}:17").endswith("bit 17")
    assert features.title_of("c:Measured LogD") == "Measured LogD"
    recipe = a_recipe(fingerprints=(COUNTS,))
    pool = features.build(SMILES, recipe)
    assert pool.title_of(0) == "MW"


def test_progress_counts_up_to_the_end():
    seen: list[tuple[int, int]] = []
    recipe = a_recipe(fingerprints=(COUNTS,), columns=("LogD",))
    features.build(
        SMILES, recipe, columns={"LogD": np.arange(5.0)},
        progress=lambda done, total: seen.append((done, total)),
    )  # fmt: skip
    assert seen
    assert [done for done, _ in seen] == sorted(done for done, _ in seen)
    assert seen[-1][0] == seen[-1][1]


def test_cancel_stops_the_pool():
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        features.build(SMILES, a_recipe(), cancel=token)


def test_a_pool_whose_parts_disagree_is_refused():
    with pytest.raises(ValueError, match="columns, names and kinds"):
        features.Pool(np.zeros((3, 2)), ("one",), ("d",), np.ones(3, dtype=bool))
    with pytest.raises(ValueError, match="rows and usable"):
        features.Pool(np.zeros((3, 1)), ("one",), ("d",), np.ones(2, dtype=bool))
