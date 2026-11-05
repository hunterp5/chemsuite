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

"""Lossless parsing and display of plain numbers in table cells.

A cell typed as ``1.50`` must still read ``1.50`` after a save, a sort or an undo, so a number
column keeps the value as ``float64`` and, per cell, how it was written:

* ``decimals >= 0``  fixed-point with that many decimals (``"1.50"`` -> 2)
* ``INTEGER`` (-1)   written as an integer (``"12"``)
* ``GENERAL`` (-2)   produced by a computation, not typed; shown compactly (:func:`format_general`)

Text that would not survive that round trip stays text: leading zeros (``"00123"``), a leading
``+``, scientific notation, more than 15 significant digits, ``nan``/``inf``, qualifiers such as
``">10"``. Converting such a column to numbers is an explicit user action.
"""

from __future__ import annotations

import math
import numbers

import numpy as np

INTEGER = -1
GENERAL = -2

_MAX_SIGNIFICANT = 15  # float64 round-trips any decimal string with <= 15 significant digits
_DIGITS = frozenset("0123456789")


def classify_number(text: str) -> tuple[float, int] | None:
    """Return ``(value, decimals_code)`` if *text* is a plain number that round-trips exactly.

    >>> classify_number("1.50")
    (1.5, 2)
    >>> classify_number("-12")
    (-12.0, -1)
    >>> classify_number("007") is None
    True
    """
    if not text:
        return None
    body = text[1:] if text[0] == "-" else text
    whole, dot, frac = body.partition(".")
    if not whole or not _DIGITS.issuperset(whole):
        return None
    if dot and (not frac or not _DIGITS.issuperset(frac)):
        return None
    if len(whole) > 1 and whole[0] == "0":
        return None
    significant = (whole + frac).lstrip("0")
    if len(significant) > _MAX_SIGNIFICANT:
        return None
    if text == "-0" or (text.startswith("-") and not significant):
        return None  # "-0" / "-0.00" would display without the sign
    return float(text), (len(frac) if dot else INTEGER)


def format_number(value: float, decimals: int) -> str:
    """Display *value* with its stored *decimals* code (see module docstring)."""
    if decimals == INTEGER:
        return str(int(value))
    if decimals == GENERAL:
        return format_general(value)
    return f"{value:.{decimals}f}"


def format_general(value: float) -> str:
    """Compact display for computed values: 4 decimals, trailing zeros trimmed.

    Very small or very large magnitudes switch to scientific notation so they do not read as 0.
    """
    if not math.isfinite(value):
        return ""
    magnitude = abs(value)
    if magnitude != 0.0 and (magnitude < 1e-3 or magnitude >= 1e12):
        return f"{value:.3e}"
    text = f"{value:.4f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def decimals_for_value(value: object) -> tuple[float, int] | None:
    """``(value, code)`` for a Python number produced by code rather than typed by a user."""
    # numbers.Integral / numbers.Real also cover numpy scalars (np.int64, np.float32, ...).
    if isinstance(value, (bool, np.bool_)):
        return None
    if isinstance(value, numbers.Integral):
        return (float(value), INTEGER) if abs(int(value)) <= 2**53 else None
    if isinstance(value, numbers.Real):
        as_float = float(value)
        return (as_float, GENERAL) if math.isfinite(as_float) else None
    return None
