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
"""Column statistics: summaries, correlations, outliers, curve fits and hypothesis tests.

Every function takes float arrays with NaN for missing values. SciPy is imported only by the
functions that need it.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from chemsuite.core.jobs import CancelToken, check

# -- summaries ---------------------------------------------------------------------------------


@dataclass(frozen=True)
class Summary:
    name: str
    n: int
    missing: int
    mean: float
    sd: float
    minimum: float
    p5: float
    q1: float
    median: float
    q3: float
    p95: float
    maximum: float
    skew: float
    kurtosis: float  # excess: 0 for a normal distribution


#: (heading, Summary attribute, help) in display order.
SUMMARY_FIELDS: tuple[tuple[str, str, str], ...] = (
    ("N", "n", "Rows with a number"),
    ("Missing", "missing", "Rows that are blank or not a number"),
    ("Mean", "mean", ""),
    ("SD", "sd", "Sample standard deviation"),
    ("Min", "minimum", ""),
    ("5%", "p5", "5th percentile"),
    ("Q1", "q1", "25th percentile"),
    ("Median", "median", ""),
    ("Q3", "q3", "75th percentile"),
    ("95%", "p95", "95th percentile"),
    ("Max", "maximum", ""),
    ("Skew", "skew", "Sample skewness (0 for a symmetric distribution)"),
    ("Kurtosis", "kurtosis", "Excess kurtosis (0 for a normal distribution)"),
)


def summarize(name: str, values: np.ndarray) -> Summary:
    values = np.asarray(values, dtype=np.float64)
    x = values[np.isfinite(values)]
    n = len(x)
    nan = math.nan
    if n == 0:
        return Summary(name, 0, len(values), *([nan] * 11))
    mean = float(x.mean())
    sd = float(x.std(ddof=1)) if n > 1 else nan
    p5, q1, median, q3, p95 = (float(v) for v in np.percentile(x, [5, 25, 50, 75, 95]))
    skew = kurtosis = nan
    if n > 2:
        m2 = float(np.mean((x - mean) ** 2))
        if m2 > 0:
            m3 = float(np.mean((x - mean) ** 3))
            # Adjusted Fisher-Pearson skewness, as Excel and pandas report it.
            skew = math.sqrt(n * (n - 1)) / (n - 2) * m3 / m2**1.5
            if n > 3:
                m4 = float(np.mean((x - mean) ** 4))
                g2 = m4 / m2**2 - 3
                kurtosis = (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6)
    return Summary(
        name,
        n,
        len(values) - n,
        mean,
        sd,
        float(x.min()),
        p5,
        q1,
        median,
        q3,
        p95,
        float(x.max()),
        skew,
        kurtosis,
    )


# -- correlation -------------------------------------------------------------------------------

CORRELATION_METHODS = {"pearson": "Pearson r", "spearman": "Spearman ρ"}


@dataclass(frozen=True)
class Correlation:
    names: tuple[str, ...]
    r: np.ndarray  # k × k, NaN where fewer than 3 pairs or a constant column
    n: np.ndarray  # k × k pair counts


def _ranks(x: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return rankdata(x)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    dx, dy = x - x.mean(), y - y.mean()
    denominator = math.sqrt(float(dx @ dx) * float(dy @ dy))
    return float(dx @ dy) / denominator if denominator > 0 else math.nan


def _pairwise_pearson(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pearson r and pair counts for every column pair at once (rows × k, NaN = missing).

    Sums over the rows both columns share come from matrix products, so 50 columns of 200k
    rows take well under a second.
    """
    present = np.isfinite(data)
    # Centre each column on its own mean first: keeps the sums small and accurate.
    with np.errstate(invalid="ignore"):
        centred = np.where(present, data - np.nanmean(np.where(present, data, np.nan), axis=0), 0.0)
    m = present.astype(np.float64)
    n = m.T @ m
    sums = centred.T @ m  # [i, j]: sum of column i over the rows j also has
    squares = (centred**2).T @ m
    products = centred.T @ centred
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = products - sums * sums.T / n
        var = squares - sums**2 / n
        r = cov / np.sqrt(var * var.T)
    r[(n < 3) | ~np.isfinite(r)] = math.nan
    return np.clip(r, -1.0, 1.0), n.astype(np.int64)


def correlation_matrix(
    names: Sequence[str],
    columns: Sequence[np.ndarray],
    method: str = "pearson",
    *,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> Correlation:
    """Pairwise-complete correlations: each pair uses the rows where both have a number.

    Spearman re-ranks each pair whose columns have blanks in different rows, which is slow on
    big tables: *progress* and *cancel* cover that part.
    """
    if method not in CORRELATION_METHODS:
        raise ValueError(f"unknown correlation method: {method!r}")
    k = len(columns)
    if not k:
        return Correlation((), np.empty((0, 0)), np.empty((0, 0), dtype=np.int64))
    data = np.column_stack([np.asarray(c, dtype=np.float64) for c in columns])
    if method == "pearson":
        r, n = _pairwise_pearson(data)
        return Correlation(tuple(names), r, n)
    # Spearman: rank each column over its own values; a pair whose rows differ is re-ranked.
    present = np.isfinite(data)
    ranked = np.full_like(data, math.nan)
    for i in range(k):
        ranked[present[:, i], i] = _ranks(data[present[:, i], i])
    r, n = _pairwise_pearson(ranked)
    counts = present.sum(axis=0)
    ragged = [
        (i, j)
        for i in range(k)
        for j in range(i + 1, k)
        if n[i, j] >= 3 and (n[i, j] != counts[i] or n[i, j] != counts[j])
    ]
    for done, (i, j) in enumerate(ragged):
        check(cancel)
        both = present[:, i] & present[:, j]
        r[i, j] = r[j, i] = _pearson(_ranks(data[both, i]), _ranks(data[both, j]))
        if progress is not None:
            progress(done + 1, len(ragged))
    return Correlation(tuple(names), r, n)


# -- outliers ----------------------------------------------------------------------------------

OUTLIER_METHODS = {
    "iqr": ("IQR (Tukey fences)", "k", 1.5),
    "zscore": ("Z-score", "|Z| above", 3.0),
    "modified_z": ("Modified Z (median and MAD)", "score above", 3.5),
}


@dataclass(frozen=True)
class Outliers:
    mask: np.ndarray  # bool per input value
    description: str  # how the threshold was applied, with the numbers used


def find_outliers(
    values: np.ndarray, method: str = "iqr", threshold: float | None = None
) -> Outliers:
    """Flag outliers; *threshold* defaults to the method's usual value (see OUTLIER_METHODS)."""
    if method not in OUTLIER_METHODS:
        raise ValueError(f"unknown outlier method: {method!r}")
    title, _label, default = OUTLIER_METHODS[method]
    t = default if threshold is None else float(threshold)
    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    x = v[ok]
    mask = np.zeros(len(v), dtype=np.bool_)
    if method == "iqr":
        if len(x) < 4:
            return Outliers(mask, "Needs at least 4 numbers.")
        q1, q3 = np.percentile(x, [25, 75])
        iqr = q3 - q1
        lo, hi = q1 - t * iqr, q3 + t * iqr
        if iqr > 0:
            mask[ok] = (x < lo) | (x > hi)
        return Outliers(
            mask, f"{title}, k = {t:g}: Q1 {q1:.4g}, Q3 {q3:.4g}, fences [{lo:.4g}, {hi:.4g}]"
        )
    if method == "zscore":
        if len(x) < 2:
            return Outliers(mask, "Needs at least 2 numbers.")
        mean, sd = x.mean(), x.std(ddof=1)
        if sd > 0:
            mask[ok] = np.abs(x - mean) / sd > t
        return Outliers(mask, f"{title}, |Z| > {t:g}: mean {mean:.4g}, SD {sd:.4g}")
    if len(x) < 3:
        return Outliers(mask, "Needs at least 3 numbers.")
    median = np.median(x)
    mad = np.median(np.abs(x - median))
    if mad > 0:
        mask[ok] = 0.6745 * np.abs(x - median) / mad > t
    return Outliers(mask, f"{title}, score > {t:g}: median {median:.4g}, MAD {mad:.4g}")


# -- curve fits --------------------------------------------------------------------------------

CURVE_MODELS = {
    "polynomial": "Polynomial",
    "exponential": "Exponential  y = a·e^(bx)",
    "logx": "Logarithmic  y = a + b·ln x",
    "power": "Power  y = a·x^b",
}


@dataclass(frozen=True)
class CurveFit:
    ok: bool
    equation: str = ""
    n: int = 0
    r2: float = math.nan
    rmse: float = math.nan
    coefficients: tuple[float, ...] = ()
    message: str = ""
    model: str = ""
    degree: int = 1

    def predict(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        c = self.coefficients
        with np.errstate(all="ignore"):
            if self.model == "polynomial":
                return np.polyval(c, x)
            if self.model == "exponential":
                return c[0] * np.exp(c[1] * x)
            if self.model == "logx":
                return c[0] + c[1] * np.log(x)
            if self.model == "power":
                return c[0] * x ** c[1]
        raise ValueError(f"unknown model: {self.model!r}")

    def report(self) -> str:
        if not self.ok:
            return self.message
        return "\n".join(
            [
                self.equation,
                "",
                f"N (pairs):  {self.n:,}",
                f"R²:         {self.r2:.4f}",
                f"RMSE:       {self.rmse:.4g}",
            ]
            + ([self.message] if self.message else [])
        )


def _r2_rmse(y: np.ndarray, fitted: np.ndarray) -> tuple[float, float]:
    residual = float(np.sum((y - fitted) ** 2))
    total = float(np.sum((y - y.mean()) ** 2))
    return (1 - residual / total if total > 0 else math.nan), math.sqrt(residual / len(y))


def _term(coefficient: float, power: int, x: str) -> str:
    if power == 0:
        return f"{coefficient:.4g}"
    return f"{coefficient:.4g}·{x}" + (f"^{power}" if power > 1 else "")


def fit_curve(
    x: np.ndarray,
    y: np.ndarray,
    model: str = "polynomial",
    *,
    degree: int = 1,
    x_name: str = "x",
    y_name: str = "y",
) -> CurveFit:
    """Least-squares fit of *model* to the rows where both values are numbers."""
    if model not in CURVE_MODELS:
        raise ValueError(f"unknown model: {model!r}")
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    needed = degree + 1 if model == "polynomial" else 2
    if n < max(needed, 2) + 1:
        return CurveFit(
            False, message=f"Needs at least {max(needed, 2) + 1} rows with both values."
        )
    try:
        if model == "polynomial":
            coefficients = tuple(float(c) for c in np.polyfit(x, y, degree))
            terms = [_term(c, degree - i, x_name) for i, c in enumerate(coefficients)]
            equation = f"{y_name} = " + " + ".join(terms).replace("+ -", "− ")
        elif model == "logx":
            if (x <= 0).any():
                return CurveFit(False, message="A logarithmic fit needs every x above 0.")
            b, a = np.polyfit(np.log(x), y, 1)
            coefficients = (float(a), float(b))
            equation = f"{y_name} = {a:.4g} + {b:.4g}·ln({x_name})"
        elif model == "power":
            if (x <= 0).any() or (y <= 0).any():
                return CurveFit(False, message="A power fit needs every x and y above 0.")
            b, log_a = np.polyfit(np.log(x), np.log(y), 1)
            coefficients = (float(np.exp(log_a)), float(b))
            equation = f"{y_name} = {coefficients[0]:.4g}·{x_name}^{b:.4g}"
        else:
            if (y <= 0).any():
                return CurveFit(False, message="An exponential fit needs every y above 0.")
            from scipy.optimize import curve_fit

            b0, log_a0 = np.polyfit(x, np.log(y), 1)  # log-linear start, then refine
            popt, _ = curve_fit(
                lambda t, a, b: a * np.exp(b * t), x, y, p0=(np.exp(log_a0), b0), maxfev=20_000
            )
            coefficients = (float(popt[0]), float(popt[1]))
            equation = f"{y_name} = {coefficients[0]:.4g}·e^({coefficients[1]:.4g}·{x_name})"
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as err:
        return CurveFit(False, message=f"The fit failed: {err}")
    fit = CurveFit(True, equation, n, coefficients=coefficients, model=model, degree=degree)
    r2, rmse = _r2_rmse(y, fit.predict(x))
    return CurveFit(True, equation, n, r2, rmse, coefficients, model=model, degree=degree)


# -- hypothesis tests --------------------------------------------------------------------------

TESTS = {
    "welch_t": ("Welch t-test (two columns, independent)", 2),
    "paired_t": ("Paired t-test (two columns, same rows)", 2),
    "mann_whitney": ("Mann-Whitney U (two columns, independent)", 2),
    "one_sample_t": ("One-sample t-test (against a mean)", 1),
    "shapiro": ("Shapiro-Wilk normality test", 1),
}
SHAPIRO_MAX = 5000


@dataclass(frozen=True)
class TestResult:
    ok: bool
    lines: list[str] = field(default_factory=list)
    statistic: float = math.nan
    p_value: float = math.nan

    def report(self) -> str:
        return "\n".join(self.lines)


def _p(p: float) -> str:
    return "< 0.0001" if p < 1e-4 else f"{p:.4f}"


def run_test(
    test: str,
    a: np.ndarray,
    b: np.ndarray | None = None,
    *,
    a_name: str = "A",
    b_name: str = "B",
    mean: float = 0.0,
) -> TestResult:
    """Run one test on finite values; two-sided p-values."""
    from scipy import stats

    if test not in TESTS:
        raise ValueError(f"unknown test: {test!r}")
    a = np.asarray(a, dtype=np.float64)
    fa = a[np.isfinite(a)]
    if TESTS[test][1] == 2:
        if b is None:
            return TestResult(False, ["Choose a second column."])
        b = np.asarray(b, dtype=np.float64)
        fb = b[np.isfinite(b)]
    lines: list[str] = []
    if test == "shapiro":
        if len(fa) < 3:
            return TestResult(False, ["Needs at least 3 numbers."])
        if len(fa) > SHAPIRO_MAX:
            lines.append(
                f"Using a random {SHAPIRO_MAX:,} of {len(fa):,} values (the test's limit)."
            )
            fa = np.random.default_rng(0).choice(fa, SHAPIRO_MAX, replace=False)
        stat, p = stats.shapiro(fa)
        lines += [
            f"{a_name}: W = {stat:.4f}, p = {_p(p)} (n = {len(fa):,})",
            "A small p (say below 0.05) suggests the values are not normally distributed.",
        ]
    elif test == "one_sample_t":
        if len(fa) < 2:
            return TestResult(False, ["Needs at least 2 numbers."])
        stat, p = stats.ttest_1samp(fa, mean)
        lines += [
            f"{a_name}: mean {fa.mean():.4g} vs {mean:g}",
            f"t = {stat:.4f}, p = {_p(p)} (n = {len(fa):,})",
        ]
    elif test == "paired_t":
        both = np.isfinite(a) & np.isfinite(b)
        if both.sum() < 2:
            return TestResult(False, ["Needs at least 2 rows with both values."])
        stat, p = stats.ttest_rel(a[both], b[both])
        diff = float(np.mean(a[both] - b[both]))
        lines += [
            f"{a_name} − {b_name}: mean difference {diff:.4g}",
            f"t = {stat:.4f}, p = {_p(p)} (pairs = {int(both.sum()):,})",
        ]
    else:
        minimum = 2 if test == "welch_t" else 1
        if len(fa) < minimum or len(fb) < minimum:
            return TestResult(False, [f"Each column needs at least {minimum} number(s)."])
        if test == "welch_t":
            stat, p = stats.ttest_ind(fa, fb, equal_var=False)
            lines.append(f"Means: {a_name} {fa.mean():.4g}, {b_name} {fb.mean():.4g}")
            lines.append(f"t = {stat:.4f}, p = {_p(p)}")
        else:
            stat, p = stats.mannwhitneyu(fa, fb, alternative="two-sided")
            lines.append(f"Medians: {a_name} {np.median(fa):.4g}, {b_name} {np.median(fb):.4g}")
            lines.append(f"U = {stat:.4g}, p = {_p(p)}")
        lines.append(f"n = {len(fa):,} and {len(fb):,}")
    return TestResult(True, lines, float(stat), float(p))
