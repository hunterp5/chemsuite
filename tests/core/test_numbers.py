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
from hypothesis import given
from hypothesis import strategies as st

from chemsuite.core.table.numbers import (
    GENERAL,
    INTEGER,
    classify_number,
    decimals_for_value,
    format_general,
    format_number,
)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("0", (0.0, INTEGER)),
        ("12", (12.0, INTEGER)),
        ("-12", (-12.0, INTEGER)),
        ("1.50", (1.5, 2)),
        ("0.001", (0.001, 3)),
        ("-0.5", (-0.5, 1)),
        ("123456789012345", (123456789012345.0, INTEGER)),
    ],
)
def test_plain_numbers_are_classified(text: str, expected: tuple[float, int]) -> None:
    assert classify_number(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        "",
        "007",
        "+1",
        "1e5",
        "1.",
        ".5",
        "nan",
        "inf",
        ">10",
        "1,000",
        " 1",
        "1 ",
        "-0",
        "-0.00",
        "1234567890123456",
        "0x10",
        "١٢",
    ],
)
def test_text_that_would_not_round_trip_stays_text(text: str) -> None:
    assert classify_number(text) is None


@given(
    st.integers(0, 6).flatmap(lambda places: st.decimals(-1e9, 1e9, allow_nan=False, places=places))
)
def test_classified_numbers_round_trip_exactly(value) -> None:
    text = str(value)
    parsed = classify_number(text)
    if parsed is not None:
        assert format_number(*parsed) == text


def test_format_general_is_compact_and_never_reads_as_zero() -> None:
    assert format_general(1.0) == "1"
    assert format_general(3.14159265) == "3.1416"
    assert format_general(1.2e-7) == "1.200e-07"
    assert format_general(-0.00001) == "-1.000e-05"
    assert format_general(0.0) == "0"


def test_values_from_code_are_numbers() -> None:
    assert decimals_for_value(3) == (3.0, INTEGER)
    assert decimals_for_value(np.int64(3)) == (3.0, INTEGER)
    assert decimals_for_value(np.float32(0.5)) == (0.5, GENERAL)
    assert decimals_for_value(float("nan")) is None
    assert decimals_for_value(True) is None
    assert decimals_for_value(np.bool_(True)) is None
    assert decimals_for_value("1") is None
