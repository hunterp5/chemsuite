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

"""A column at a glance, for its header tooltip: counts, and the spread of its numbers.

Computed when the tooltip is asked for, over every row (not just the shown ones), with whole-
column numpy operations, so it stays quick at 200k rows.
"""

from __future__ import annotations

import numpy as np

from .columns import Column, NumberColumn, TextColumn
from .style import FIXED, NUMBER_FORMATS, PERCENT, SCIENTIFIC, SIGNIFICANT, THOUSANDS, ColumnStyle

FORMAT_TITLES = {
    FIXED: "Fixed decimals",
    THOUSANDS: "Thousands separator",
    SCIENTIFIC: "Scientific",
    SIGNIFICANT: "Significant figures",
    PERCENT: "Percent",
}
assert set(FORMAT_TITLES) == set(NUMBER_FORMATS)


def format_title(style: ColumnStyle) -> str:
    """How the column shows its numbers, e.g. ``Significant figures (3)``."""
    unit = "significant figures" if style.number_format == SIGNIFICANT else "decimals"
    return f"{FORMAT_TITLES[style.number_format]} ({style.digits} {unit})"


def _number(value: float) -> str:
    if not np.isfinite(value):
        return "—"
    return f"{value:.4g}"


def column_summary(column: Column) -> list[str]:
    """Lines describing *column*: rows filled and blank, then by kind:

    - numbers: lowest, highest, mean ± standard deviation and median;
    - text: how many distinct values, and the commonest one when it repeats.
    """
    total = len(column)
    blanks = int(np.count_nonzero(column.blank_mask())) if total else 0
    lines = [f"{total - blanks:,} filled · {blanks:,} blank"]
    if isinstance(column, NumberColumn):
        values = column.values()
        values = values[np.isfinite(values)]
        if len(values):
            lines.append(f"Range {_number(values.min())} to {_number(values.max())}")
            spread = f" ± {_number(values.std(ddof=1))}" if len(values) > 1 else ""
            lines.append(f"Mean {_number(values.mean())}{spread}")
            lines.append(f"Median {_number(float(np.median(values)))}")
    elif isinstance(column, TextColumn):
        counts = [(value, count) for value, count in column.value_counts() if value != ""]
        if counts:
            lines.append(f"{len(counts):,} distinct value{'s' if len(counts) != 1 else ''}")
            value, count = max(counts, key=lambda pair: pair[1])
            if count > 1:
                shown = value if len(value) <= 30 else value[:29] + "…"
                lines.append(f"Commonest: {shown} ({count:,})")
    return lines
