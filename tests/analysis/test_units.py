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

import numpy as np
import pytest

from chemsuite.analysis.units import (
    convert,
    molecular_weights,
    needs_weight,
    parse_unit,
    parse_value,
)


@pytest.mark.parametrize(
    ("typed", "unit"),
    [
        ("nM", "nM"), ("nm", "nM"), (" NM ", "nM"), ("uM", "µM"), ("µM", "µM"), ("μM", "µM"),
        ("mcM", "µM"), ("mM", "mM"), ("M", "M"), ("mol/L", "M"), ("pM", "pM"),
        ("micromolar", "µM"), ("mg/mL", "mg/mL"), ("ug/ml", "µg/mL"), ("ng / mL", "ng/mL"),
        ("mg/L", "mg/L"), ("g/L", "g/L"),
        ("", None), (None, None), ("%", None), ("mg/kg", None), ("mg/dL", None),
    ],
)  # fmt: skip
def test_units_are_read_whatever_their_case_or_spelling(typed, unit) -> None:
    assert parse_unit(typed) == unit


def test_values_keep_their_qualifiers() -> None:
    assert parse_value(12.5) == ("", 12.5)
    assert parse_value(">10") == (">", 10.0)
    assert parse_value("≤ 0.5") == ("<=", 0.5)
    assert parse_value("1e-3") == ("", 0.001)
    q, v = parse_value("n.d.")
    assert q == "" and np.isnan(v)


def test_molar_units_convert_by_powers_of_ten() -> None:
    result = convert([1.0, 2.5, 500, ">10"], ["µM", "nM", "pM", "uM"], "nM")
    np.testing.assert_allclose(result.values, [1000.0, 2.5, 0.5, 10000.0])
    assert result.qualifiers == ["", "", "", ">"]
    assert result.unknown_units == 0 and result.no_weight == 0


def test_mass_units_convert_through_the_molecular_weight() -> None:
    """Aspirin, 180.16 g/mol: 1 mg/mL is 5.55 mM, and 1 µM is 0.18 µg/mL."""
    weights = np.array([180.16, 180.16, np.nan])
    result = convert([1.0, 1.0, 1.0], ["mg/mL", "µg/mL", "mg/mL"], "mM", weights)
    np.testing.assert_allclose(result.values[:2], [1 / 180.16 * 1e3, 1 / 180.16], rtol=1e-9)
    assert np.isnan(result.values[2]) and result.no_weight == 1
    back = convert([1.0], ["µM"], "µg/mL", weights[:1])
    np.testing.assert_allclose(back.values, [0.18016])
    assert needs_weight(["nM", "mg/mL"], "nM") and not needs_weight(["nM", "µM"], "pM")


def test_rows_without_a_unit_or_a_number_stay_blank() -> None:
    result = convert([5.0, "abc", None, 3.0], ["parsecs", "nM", "nM", None], "nM")
    assert np.isnan(result.values).all()
    assert result.unknown_units == 2  # a number with no unit read; blanks are not counted


def test_molecular_weights_of_structure_cells() -> None:
    weights = molecular_weights(["CC(=O)Oc1ccccc1C(=O)O", "", None, "not a smiles"])
    assert weights[0] == pytest.approx(180.159, abs=0.01)
    assert np.isnan(weights[1:]).all()


def test_cancel_keeps_the_rows_already_converted(monkeypatch) -> None:
    from chemsuite.analysis import units
    from chemsuite.core.jobs import Cancelled, CancelToken

    monkeypatch.setattr(units, "CHUNK_ROWS", 2)
    token = CancelToken()

    def progress(done: int, _total: int) -> None:
        if done == 2:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        units.convert_rows([1, 2, 3, 4, 5], ["µM"] * 5, "nM", progress=progress, cancel=token)
    np.testing.assert_allclose(stopped.value.partial.values, [1000.0, 2000.0])
    done = units.convert_rows([1, ">2"], ["nM", "mg/mL"], "nM", ["", "CCO"])
    assert done.values[0] == 1 and done.qualifiers == ["", ">"] and done.values[1] > 0
