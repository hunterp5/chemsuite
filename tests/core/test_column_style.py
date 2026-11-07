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

"""Column styles: precision, Logarithmic, colours, duplicates; all undoable and saved."""

from __future__ import annotations

import math

import numpy as np
import pytest

from chemsuite.core.document import load_document, save_document
from chemsuite.core.table import NUMBER, TEXT, ColourRule, ColumnStyle, Table
from chemsuite.core.table.style import CATEGORIES, GRADIENT, GRADIENT3, format_fixed


def _table() -> Table:
    table = Table()
    table.add_rows({"IC50": ["1.5", "10", "", "250"], "Series": ["A", "B", "A", "C"]})
    return table


def _shown(table: Table, name: str) -> list[str]:
    column = table.column(name)
    return column.displays(np.arange(len(column)))


def test_precision_changes_the_display_not_the_values():
    table = _table()
    before = table.column("IC50").values().copy()
    table.set_style("IC50", ColumnStyle(precision=2))
    assert _shown(table, "IC50") == ["1.50", "10.00", "", "250.00"]
    np.testing.assert_array_equal(table.column("IC50").values(), before)
    table.undo()
    assert _shown(table, "IC50") == ["1.5", "10", "", "250"]


def test_precision_applies_to_search_text():
    table = _table()
    table.set_style("IC50", ColumnStyle(precision=1))
    hits = table.column("IC50").match_text(lambda text: text == "10.0")
    assert hits.tolist() == [False, True, False, False]


@pytest.mark.parametrize(
    ("value", "decimals", "text"),
    [(-0.001, 2, "0.00"), (-0.4, 0, "0"), (-1.25, 1, "-1.2"), (2.675, 2, "2.67"), (7, 0, "7")],
)
def test_format_fixed(value, decimals, text):
    # Python rounds the binary value: 2.675 is 2.67499999… and -1.25 ties to even.
    assert format_fixed(value, decimals) == text


def test_logarithmic_converts_and_back():
    table = _table()
    assert table.set_log10("IC50", True) == 3
    column = table.column("IC50")
    assert column.style.log10
    np.testing.assert_allclose(
        column.values(), [math.log10(1.5), 1.0, np.nan, math.log10(250)], equal_nan=True
    )
    table.set_log10("IC50", False)
    assert not table.column("IC50").style.log10
    assert _shown(table, "IC50") == ["1.5", "10", "", "250"]
    np.testing.assert_array_equal(table.column("IC50").values()[[0, 1, 3]], [1.5, 10.0, 250.0])


def test_logarithmic_is_one_undo_step():
    table = _table()
    table.set_log10("IC50", True)
    table.undo()
    assert not table.column("IC50").style.log10
    assert _shown(table, "IC50") == ["1.5", "10", "", "250"]


def test_logarithmic_refuses_zero_and_negative():
    table = Table()
    table.add_rows({"x": ["1", "0"]})
    with pytest.raises(ValueError, match="zero or negative"):
        table.set_log10("x", True)
    table.add_rows({"t": ["a", "b", "c"]})
    with pytest.raises(ValueError, match="not a number"):
        table.set_log10("t", True)


def test_duplicate_column_copies_cells_and_style_to_the_right():
    table = _table()
    table.set_style("IC50", ColumnStyle(precision=1, colour=ColourRule()))
    name = table.duplicate_column("IC50")
    assert name == "IC50 copy"
    assert table.column_names == ["IC50", "IC50 copy", "Series"]
    assert _shown(table, name) == _shown(table, "IC50")
    assert table.column(name).style == table.column("IC50").style
    assert table.duplicate_column("IC50") == "IC50 copy (2)"
    table.undo()
    table.undo()
    assert table.column_names == ["IC50", "Series"]


def test_style_follows_rename_and_delete_undo():
    table = _table()
    style = ColumnStyle(precision=3)
    table.set_style("IC50", style)
    table.rename_column("IC50", "IC50 (nM)")
    assert table.column("IC50 (nM)").style == style
    table.remove_column("IC50 (nM)")
    table.undo()
    assert table.column("IC50 (nM)").style == style


def test_convert_keeps_only_what_still_applies():
    table = _table()
    table.set_style("IC50", ColumnStyle(precision=2, colour=ColourRule(mode=GRADIENT)))
    table.convert_column("IC50", TEXT)
    assert table.column("IC50").style == ColumnStyle()
    table.set_style("Series", ColumnStyle(colour=ColourRule(mode=CATEGORIES)))
    table.convert_column("Series", NUMBER, force=True)
    assert table.column("Series").style.colour.mode == CATEGORIES


def test_snapshot_keeps_styles():
    table = _table()
    table.set_style("IC50", ColumnStyle(precision=2))
    assert table.snapshot().column("IC50").style.precision == 2


def test_document_round_trips_styles(tmp_path):
    table = _table()
    rule = ColourRule(mode=GRADIENT3, minimum=1.0, maximum=None, low=0x112233, opacity=80)
    table.set_style("IC50", ColumnStyle(precision=1, colour=rule))
    table.set_log10("IC50", True)
    table.set_style("Series", ColumnStyle(colour=ColourRule(mode=CATEGORIES)))
    path = tmp_path / "styled.csd"
    save_document(table, path)
    back, _extras = load_document(path)
    assert back.column("IC50").style == table.column("IC50").style
    assert back.column("Series").style == table.column("Series").style
    assert _shown(back, "IC50") == _shown(table, "IC50")


# -- colours ------------------------------------------------------------------------------


def test_gradient_ends_and_middle():
    rule = ColourRule(mode=GRADIENT, low=0x000000, high=0xFFFFFF, opacity=255)
    bounds = rule.bounds((0.0, 10.0))
    assert bounds == (0.0, 5.0, 10.0)
    assert rule.number_colour(0.0, bounds) == 0xFF000000
    assert rule.number_colour(10.0, bounds) == 0xFFFFFFFF
    assert rule.number_colour(5.0, bounds) == 0xFF808080  # round(127.5) == 128
    assert rule.number_colour(-3.0, bounds) == 0xFF000000  # clamped
    assert rule.number_colour(99.0, bounds) == 0xFFFFFFFF


def test_three_colour_gradient_passes_through_the_middle_colour():
    rule = ColourRule(mode=GRADIENT3, low=0xFF0000, mid=0x00FF00, high=0x0000FF, opacity=10)
    bounds = rule.bounds((0.0, 100.0))
    assert rule.number_colour(50.0, bounds) == (10 << 24) | 0x00FF00
    assert rule.number_colour(0.0, bounds) == (10 << 24) | 0xFF0000
    assert rule.number_colour(100.0, bounds) == (10 << 24) | 0x0000FF


def test_fixed_bounds_override_the_data():
    rule = ColourRule(minimum=2.0, maximum=4.0, middle=3.5)
    assert rule.bounds((0.0, 100.0)) == (2.0, 3.5, 4.0)
    assert ColourRule().bounds(None) is None


def test_categories_are_stable_and_distinct():
    rule = ColourRule(mode=CATEGORIES)
    assert rule.text_colour("A") == rule.text_colour("A")
    assert len({rule.text_colour(t) for t in ("A", "B", "C", "D")}) == 4
    assert rule.text_colour("A") >> 24 == rule.opacity


def test_bad_values_are_rejected():
    with pytest.raises(ValueError):
        ColourRule(mode="rainbow")
    with pytest.raises(ValueError):
        ColumnStyle(precision=13)


def test_number_formats_show_every_value_the_columns_way():
    from chemsuite.core.table.style import format_with

    assert format_with(1234567.891, "thousands", 2) == "1,234,567.89"
    assert format_with(0.000123456, "scientific", 2) == "1.23e-04"
    assert format_with(1234.5678, "significant", 3) == "1230"
    assert format_with(0.0012345, "significant", 2) == "0.0012"
    assert format_with(0.4567, "percent", 1) == "45.7%"
    assert format_with(-0.0001, "thousands", 2) == "0.00"  # never -0.00
    table = _table()
    table.set_style("IC50", ColumnStyle(number_format="percent"))  # default 1 decimal
    shown = [table.cell("IC50", int(o)) for o in table.oids]
    assert all(text.endswith("%") for text in shown if text)
    with pytest.raises(ValueError):
        ColumnStyle(number_format="roman")
    with pytest.raises(ValueError):
        ColumnStyle(precision=0, number_format="significant")


def test_number_format_is_saved_and_searched():
    from chemsuite.core.table.search import Query, find

    table = _table()
    style = ColumnStyle(precision=1, number_format="percent")
    table.set_style("IC50", style)
    assert ColumnStyle.from_json(style.to_json()) == style
    assert ColumnStyle.from_json({"number_format": "nonsense"}).number_format == "fixed"
    assert len(find(table, Query("%"), ["IC50"])) == sum(
        1 for o in table.oids if table.cell("IC50", int(o))
    )
