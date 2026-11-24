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
"""Filling and scaling, as arrays a model archive can carry."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.qsar.prepare import Prepare


def test_a_hole_is_filled_with_the_training_median():
    x = np.array([[1.0, 10.0], [3.0, np.inf], [5.0, 30.0], [np.nan, 50.0]])
    prepare = Prepare.fit(x)
    assert prepare.fill.tolist() == [3.0, 30.0]
    filled = prepare.filled(x)
    assert np.isfinite(filled).all()
    assert filled[3, 0] == 3.0
    assert filled[1, 1] == 30.0


def test_filling_matches_scikit_learns_imputer():
    from sklearn.impute import SimpleImputer

    rng = np.random.default_rng(0)
    x = rng.normal(size=(40, 6))
    x[rng.random(x.shape) < 0.2] = np.nan
    mine = Prepare.fit(x).filled(x)
    theirs = SimpleImputer(strategy="median").fit_transform(x)
    assert np.allclose(mine, theirs)


def test_a_column_with_no_value_at_all_is_filled_with_zero():
    x = np.array([[np.nan, 1.0], [np.nan, 2.0]])
    assert Prepare.fit(x).fill.tolist() == [0.0, 1.5]  # wholly empty, then the median of 1 and 2


def test_new_rows_are_filled_from_the_training_fold_not_themselves():
    train = np.array([[1.0], [1.0], [1.0]])
    prepare = Prepare.fit(train)
    later = np.array([[np.nan], [99.0]])
    assert prepare.filled(later)[0, 0] == 1.0


def test_the_design_matrix_fills_takes_and_scales_in_that_order():
    x = np.array([[1.0, 100.0, 7.0], [3.0, 300.0, 7.0], [5.0, np.nan, 7.0]])
    kept = np.array([0, 1])
    prepare = Prepare.fit(x)
    scaled = prepare.with_scaling(prepare.filled(x), kept, scaled=True)
    design = scaled.design(x, kept)
    assert design.shape == (3, 2)
    assert np.allclose(design.mean(axis=0), 0.0)
    assert np.allclose(design.std(axis=0), 1.0)


def test_unscaled_keeps_the_numbers_a_tree_wants():
    x = np.array([[1.0, 100.0], [3.0, 300.0]])
    prepare = Prepare.fit(x)
    plain = prepare.with_scaling(prepare.filled(x), np.array([0, 1]), scaled=False)
    assert plain.center is None and plain.scale is None
    assert np.allclose(plain.design(x, np.array([0, 1])), x)


def test_a_flat_column_is_not_divided_by_its_noise():
    x = np.array([[2.0, 1.0], [2.0, 5.0], [2.0, 9.0]])
    prepare = Prepare.fit(x)
    scaled = prepare.with_scaling(prepare.filled(x), np.array([0, 1]), scaled=True)
    assert scaled.scale[0] == 1.0
    assert np.isfinite(scaled.design(x, np.array([0, 1]))).all()


def test_arrays_round_trip():
    x = np.array([[1.0, 10.0], [3.0, 30.0]])
    prepare = Prepare.fit(x)
    for scaled in (False, True):
        fitted = prepare.with_scaling(prepare.filled(x), np.array([0, 1]), scaled=scaled)
        back = Prepare.from_arrays(fitted.as_arrays())
        assert np.allclose(back.design(x, np.array([0, 1])), fitted.design(x, np.array([0, 1])))


def test_the_wrong_number_of_columns_is_refused():
    prepare = Prepare.fit(np.zeros((3, 4)))
    with pytest.raises(ValueError, match="the fit saw 4"):
        prepare.filled(np.zeros((3, 5)))


def test_an_archive_without_fill_values_is_refused():
    with pytest.raises(ValueError, match="no fill values"):
        Prepare.from_arrays({})
