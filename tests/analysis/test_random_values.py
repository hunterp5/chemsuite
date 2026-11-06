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

"""Random numbers and samples: ranges, seeds and settings checks."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.analysis.random_values import (
    INTEGER,
    LOGNORMAL,
    NORMAL,
    ORDER,
    UNIFORM,
    RandomSpec,
    random_numbers,
    random_sample,
    sample_size,
)


def test_uniform_stays_in_range_after_rounding():
    values = random_numbers(RandomSpec(UNIFORM, low=0, high=1, decimals=0, seed=1), 2000)
    assert set(np.unique(values)) <= {0.0, 1.0}


def test_integers_include_both_ends():
    values = random_numbers(RandomSpec(INTEGER, low=1, high=3, seed=1), 3000)
    assert set(np.unique(values)) == {1.0, 2.0, 3.0}


def test_normal_and_lognormal_have_the_asked_centre():
    normal = random_numbers(RandomSpec(NORMAL, mean=5, sd=2, decimals=None, seed=2), 20_000)
    assert normal.mean() == pytest.approx(5, abs=0.05)
    assert normal.std() == pytest.approx(2, abs=0.05)
    lognormal = random_numbers(RandomSpec(LOGNORMAL, mean=10, sd=2, decimals=None, seed=2), 20_000)
    assert np.median(lognormal) == pytest.approx(10, rel=0.03)


def test_clip_and_order():
    clipped = random_numbers(RandomSpec(NORMAL, mean=0, sd=10, clip=True, low=-1, high=1), 500)
    assert clipped.min() >= -1 and clipped.max() <= 1
    order = random_numbers(RandomSpec(ORDER, seed=3), 10)
    assert sorted(order) == list(range(1, 11))


def test_same_seed_same_numbers():
    spec = RandomSpec(UNIFORM, seed=42)
    assert np.array_equal(random_numbers(spec, 50), random_numbers(spec, 50))
    assert not np.array_equal(random_numbers(RandomSpec(), 50), random_numbers(RandomSpec(), 50))


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        (RandomSpec(UNIFORM, low=2, high=1), "Low must not"),
        (RandomSpec(INTEGER, low=1.2, high=1.8), "no whole number"),
        (RandomSpec(LOGNORMAL, mean=0), "median"),
        (RandomSpec(LOGNORMAL, mean=1, sd=0.5), "spread factor"),
    ],
)
def test_problems(spec, message):
    assert message in spec.problem()
    with pytest.raises(ValueError):
        random_numbers(spec, 3)


def test_sample():
    oids = np.arange(100, 200)
    picked = random_sample(oids, 10, seed=5)
    assert len(set(picked)) == 10 and set(picked) <= set(oids)
    assert np.array_equal(picked, random_sample(oids, 10, seed=5))
    assert len(random_sample(oids, 500)) == 100
    assert sample_size(37, 10, True) == 4 and sample_size(37, 50, False) == 37
