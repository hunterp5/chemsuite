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

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from chemsuite.analysis.statistics import (
    correlation_matrix,
    find_outliers,
    fit_curve,
    run_test,
    summarize,
)

RNG = np.random.default_rng(7)


def test_summary_matches_scipy() -> None:
    x = RNG.normal(5, 2, 500)
    values = np.r_[x, np.nan, np.nan]
    s = summarize("x", values)
    assert (s.n, s.missing) == (500, 2)
    assert s.mean == pytest.approx(x.mean())
    assert s.sd == pytest.approx(x.std(ddof=1))
    assert s.median == pytest.approx(np.median(x))
    assert s.skew == pytest.approx(stats.skew(x, bias=False))
    assert s.kurtosis == pytest.approx(stats.kurtosis(x, bias=False))


def test_summary_of_tiny_and_empty_columns() -> None:
    assert summarize("e", np.array([np.nan])).n == 0
    one = summarize("one", np.array([3.0]))
    assert one.mean == 3 and math.isnan(one.sd) and math.isnan(one.skew)
    flat = summarize("flat", np.full(10, 2.0))
    assert flat.sd == 0 and math.isnan(flat.skew)


def test_correlations_are_pairwise_complete() -> None:
    a = RNG.normal(size=200)
    b = 2 * a + RNG.normal(scale=0.5, size=200)
    c = np.exp(a)
    b_missing = b.copy()
    b_missing[:50] = np.nan
    pearson = correlation_matrix(["a", "b", "c"], [a, b_missing, c])
    assert pearson.n[0, 1] == 150 and pearson.n[0, 2] == 200
    assert pearson.r[0, 1] == pytest.approx(stats.pearsonr(a[50:], b[50:])[0])
    spearman = correlation_matrix(["a", "c"], [a, c], "spearman")
    assert spearman.r[0, 1] == pytest.approx(1.0)
    assert np.allclose(np.diag(spearman.r), 1.0)
    constant = correlation_matrix(["a", "k"], [a, np.ones(200)])
    assert math.isnan(constant.r[0, 1]) and math.isnan(constant.r[1, 1])
    with pytest.raises(ValueError):
        correlation_matrix(["a"], [a], "kendall")


def test_outlier_methods() -> None:
    values = np.r_[RNG.normal(0, 1, 100), 12.0, np.nan]
    for method in ("iqr", "zscore", "modified_z"):
        result = find_outliers(values, method)
        assert result.mask[100] and not result.mask[101], method
        assert result.mask.sum() < 5
    assert not find_outliers(np.array([1.0, 2.0]), "iqr").mask.any()
    assert find_outliers(values, "zscore", threshold=100).mask.sum() == 0


def test_curve_fits_recover_parameters() -> None:
    x = np.linspace(1, 10, 40)
    linear = fit_curve(x, 3 * x - 2, x_name="MW", y_name="pIC50")
    assert linear.coefficients == pytest.approx((3, -2))
    assert linear.r2 == pytest.approx(1) and linear.equation == "pIC50 = 3·MW − 2"
    assert fit_curve(x, 2 * np.exp(0.3 * x), "exponential").coefficients == pytest.approx((2, 0.3))
    assert fit_curve(x, 1 + 4 * np.log(x), "logx").coefficients == pytest.approx((1, 4))
    assert fit_curve(x, 5 * x**1.5, "power").coefficients == pytest.approx((5, 1.5))
    quadratic = fit_curve(x, x**2 - x, degree=2)
    assert quadratic.predict(np.array([2.0]))[0] == pytest.approx(2)
    assert not fit_curve(x - 5, x, "logx").ok
    assert not fit_curve(x[:2], x[:2]).ok
    assert "Needs" in fit_curve(x[:2], x[:2]).report()


def test_hypothesis_tests_match_scipy() -> None:
    a = RNG.normal(0, 1, 40)
    b = RNG.normal(0.8, 1, 50)
    welch = run_test("welch_t", a, b)
    assert welch.p_value == pytest.approx(stats.ttest_ind(a, b, equal_var=False).pvalue)
    paired = run_test("paired_t", a, a + 0.5 + RNG.normal(0, 0.1, 40))
    assert paired.ok and paired.p_value < 1e-6 and "< 0.0001" in paired.report()
    assert run_test("mann_whitney", a, b).ok
    one = run_test("one_sample_t", a, mean=0.0)
    assert one.p_value == pytest.approx(stats.ttest_1samp(a, 0).pvalue)
    assert run_test("shapiro", a).p_value > 0.01
    assert not run_test("welch_t", a).ok
    assert not run_test("shapiro", np.array([1.0, np.nan])).ok
    big = run_test("shapiro", RNG.normal(size=6000))
    assert big.ok and "5,000" in big.report()
