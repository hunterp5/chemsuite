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

from chemsuite.analysis.binning import BinSpec, assign, bin_values, parse_names, parse_numbers
from chemsuite.analysis.mpo import Desirability, MpoProfile, combine
from chemsuite.analysis.selectivity import looks_like_assay, selectivity

NAN = math.nan


# -- MPO ---------------------------------------------------------------------------------------
def test_linear_shapes() -> None:
    x = np.array([0, 2, 5, 8, 10, NAN])
    up = Desirability("a", "maximize", 2, 8).score(x)
    assert np.allclose(up[:5], [0, 0, 0.5, 1, 1]) and math.isnan(up[5])
    down = Desirability("a", "minimize", 2, 8).score(x)
    assert np.allclose(down[:5], [1, 1, 0.5, 0, 0])
    target = Desirability("a", "target", 2, 10, target=4).score(x)
    assert np.allclose(target[:5], [0, 0, 5 / 6, 2 / 6, 0])
    hard = Desirability("a", "maximize", 2, 8, hard=True).score(x)
    assert np.allclose(hard[:5], [0, 0, 0, 1, 1])
    window = Desirability("a", "target", 2, 8, hard=True).score(x)
    assert np.allclose(window[:5], [0, 1, 1, 1, 0])


def test_gaussian_and_validation() -> None:
    d = Desirability("a", "gaussian", target=3, width=1).score(np.array([3.0, 4.0]))
    assert d == pytest.approx([1, math.exp(-0.5)])
    with pytest.raises(ValueError):
        Desirability("a", "sideways")
    with pytest.raises(ValueError):
        Desirability("a", weight=-1)


def test_combine_methods_and_missing_values() -> None:
    a = np.array([1.0, 0.5, 0.0, NAN])
    b = np.array([0.25, 0.5, 1.0, 0.8])
    geometric = combine([a, b], [1, 1])
    assert np.allclose(geometric[:3], [0.5, 0.5, 0]) and math.isnan(geometric[3])
    arithmetic = combine([a, b], [3, 1], "arithmetic")
    assert np.allclose(arithmetic[:3], [0.8125, 0.5, 0.25])
    skipping = combine([a, b], [1, 1], skip_missing=True)
    assert skipping[3] == pytest.approx(0.8)
    weighted = combine([np.array([0.25]), np.array([1.0])], [2, 0])
    assert weighted[0] == pytest.approx(0.25)  # weight 0: ignored


def test_profile_round_trips_through_json() -> None:
    profile = MpoProfile(
        (Desirability("cLogP", "target", 1, 4), Desirability("MW", "minimize", 350, 500, weight=2)),
        "arithmetic",
        True,
    )
    assert MpoProfile.from_json(profile.to_json()) == profile
    overall, parts = profile.score({"cLogP": np.array([2.5]), "MW": np.array([425.0])})
    assert len(parts) == 2 and overall[0] == pytest.approx((1 + 2 * 0.5) / 3)
    with pytest.raises(ValueError):
        MpoProfile.from_json('{"format": "other"}')
    with pytest.raises(ValueError):
        MpoProfile.from_json("not json")


# -- binning -----------------------------------------------------------------------------------
X = np.array([1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, NAN])


def test_equal_width_and_count() -> None:
    width = bin_values(X, BinSpec("equal_width", 3))
    assert width.range_names == ("[1, 4)", "[4, 7)", "[7, 10]")
    assert width.index.tolist() == [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, -1]
    assert width.labels(width.class_names)[0] == "Low" and width.labels(width.range_names)[-1] == ""
    count = bin_values(X, BinSpec("equal_count", 3))
    assert count.counts().tolist() == [4, 4, 3]
    ties = bin_values(np.array([1, 1, 1, 1, 2]), BinSpec("equal_count", 4))
    assert len(ties.intervals) == 2  # equal values share a bin


def test_width_and_limits() -> None:
    aligned = bin_values(X, BinSpec("equal_width", width=2.5))
    assert aligned.range_names == ("[0, 2.5)", "[2.5, 5)", "[5, 7.5)", "[7.5, 10]")
    limited = bin_values(X, BinSpec("equal_width", 2, minimum=2, maximum=6))
    assert limited.index.tolist()[:8] == [-1, 0, 0, 0, 1, 1, 1, -1]


def test_custom_edges_and_names() -> None:
    open_ends = bin_values(X, BinSpec("custom", edges=(3, 6), names=("A", "B", "C")))
    assert open_ends.range_names == ("< 3", "[3, 6)", "≥ 6")
    assert open_ends.labels(open_ends.class_names)[:5] == ["A", "A", "A", "B", "B"]
    closed = bin_values(X, BinSpec("custom", edges=(3, 6), open_ends=False))
    assert closed.index.tolist()[:8] == [-1, -1, -1, 0, 0, 0, 0, -1]
    with pytest.raises(ValueError, match="3 bins"):
        bin_values(X, BinSpec("custom", edges=(3, 6), names=("A", "B")))
    assert parse_numbers("5, 6;7 8") == (5, 6, 7, 8)
    for bad in ("", "5, x", "6, 5"):
        with pytest.raises(ValueError):
            parse_numbers(bad)
    assert parse_names(" Low, High ") == ("Low", "High") and parse_names("") == ()
    with pytest.raises(ValueError):
        parse_names("a,,b")


def test_assign_edges_and_errors() -> None:
    spec_bins = bin_values(np.array([0.0, 10.0]), BinSpec("equal_width", 2)).intervals
    assert assign(np.array([-1, 0, 4.99, 5, 10, 10.01]), spec_bins).tolist() == [-1, 0, 0, 1, 1, -1]
    assert len(bin_values(np.array([4.0, 4.0]), BinSpec("equal_width", 5)).intervals) == 1
    with pytest.raises(ValueError):
        bin_values(np.array([NAN]), BinSpec("equal_width"))
    with pytest.raises(ValueError):
        bin_values(X, BinSpec("equal_width", width=0.001))


# -- selectivity -------------------------------------------------------------------------------
def test_selectivity_on_p_scale() -> None:
    primary = np.array([8.0, 7.0, NAN, 6.0])
    result = selectivity(
        primary,
        {"T2": np.array([6, 6.5, 5, NAN]), "T3": np.array([5, NAN, 1, NAN])},
        potency_cutoff=7,
        window=1,
    )
    assert np.allclose(result.by_target["T2"][:2], [2, 0.5])
    assert result.minimum[:2].tolist() == [2, 0.5] and math.isnan(result.minimum[3])
    assert result.narrow == ["no", "yes", "", ""]


def test_selectivity_from_concentrations() -> None:
    result = selectivity(
        np.array([10.0, 100.0, 0.0]),
        {"hERG": np.array([1000.0, 200.0, 5.0])},
        concentrations=True,
        potency_cutoff=150,
        window=1,
    )
    assert np.allclose(result.by_target["hERG"][:2], [2, math.log10(2)])
    assert math.isnan(result.by_target["hERG"][2])  # zero is not a concentration
    assert result.narrow == ["no", "yes", ""]


def test_assay_names() -> None:
    assert all(looks_like_assay(n) for n in ("pIC50", "hERG IC50 (uM)", "Ki", "target_pKi"))
    assert not any(looks_like_assay(n) for n in ("MW", "Kinase family", "Notes"))
