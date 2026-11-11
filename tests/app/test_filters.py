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

import pytest

from .conftest import gui

pytestmark = gui


@pytest.fixture
def panel(qapp):
    from chemsuite.app.filters.panel import FilterPanel
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": [f"CPD-{i}" for i in range(1, 9)],
            "pIC50": ["5", "6", "7", "8", "", "6.5", "9", "5.5"],
            "Series": ["A", "B", "A", "C", "", "B", "A", "A"],
        },
        undoable=False,
    )
    p = FilterPanel(table)
    yield p
    p.close()


def shown_ids(panel) -> list[str]:
    table = panel.table
    return [table.cell("ID", int(o)) for o in table.visible_oids()]


def test_card_types_follow_the_column(panel) -> None:
    from chemsuite.app.filters.cards import CategoryCard, RangeCard, TextCard

    assert isinstance(panel.add_filter("pIC50"), RangeCard)
    assert isinstance(panel.add_filter("Series"), CategoryCard)
    assert isinstance(panel.add_filter("ID", "text"), TextCard)


def test_range_card_filters_and_negates(panel) -> None:
    card = panel.add_filter("pIC50")
    assert "8 numbers" not in card.details and "7 numbers from 5 to 9" in card.details
    card.set_bounds(6, 8)
    panel.apply_now()
    assert shown_ids(panel) == ["CPD-2", "CPD-3", "CPD-4", "CPD-5", "CPD-6"]  # blank kept
    card.blanks.setChecked(False)
    panel.apply_now()
    assert "CPD-5" not in shown_ids(panel)
    card.not_button.setChecked(True)
    panel.apply_now()
    assert shown_ids(panel) == ["CPD-1", "CPD-7", "CPD-8"]
    assert panel.table.history.is_clean  # filtering is not an edit
    assert panel.summary().startswith("Showing 3 of 8 rows")


def test_range_slider_and_edits_stay_in_step(panel) -> None:
    card = panel.add_filter("pIC50")
    card.slider.set_values(0, card.slider.steps // 2)  # 5 .. 7
    assert card.bounds() == (None, 7.0) and card.high_edit.text() == "7"
    card.low_edit.setText("6")
    card.low_edit.editingFinished.emit()
    assert card.bounds() == (6.0, 7.0)
    assert card.slider.values()[0] == card.slider.steps // 4


def test_category_card_and_unlisted_values(panel, monkeypatch) -> None:
    import chemsuite.app.filters.cards as cards

    card = panel.add_filter("Series")
    assert card.listed_values()[:2] == ["A", "B"]  # most frequent first
    card.set_checked({"A", "C"})
    panel.apply_now()
    assert shown_ids(panel) == ["CPD-1", "CPD-3", "CPD-4", "CPD-7", "CPD-8"]
    card.set_checked(None)
    panel.apply_now()
    assert len(panel.table.view()) == 8 and not panel.table.filters

    monkeypatch.setattr(cards, "CATEGORY_LIST_LIMIT", 2)  # list only A and B
    small = panel.add_filter("Series")
    assert small.listed_values() == ["A", "B"] and "most common" in small.details
    small.set_checked({"A"})  # untick B; C and blank are not listed and stay visible
    panel.apply_now()
    assert "CPD-4" in shown_ids(panel) and "CPD-2" not in shown_ids(panel)


def test_text_card_modes_and_bad_regex(panel) -> None:
    card = panel.add_filter("ID", "text")
    card.query.setText("cpd-1")
    panel.apply_now()
    assert shown_ids(panel) == ["CPD-1"]
    card.mode.setCurrentIndex(card.mode.findData("regex"))
    card.query.setText("CPD-[")  # invalid: filter is off, not an error
    panel.apply_now()
    assert len(panel.table.view()) == 8 and "Invalid" in card.query.toolTip()
    card.query.setText(r"CPD-[78]$")
    panel.apply_now()
    assert shown_ids(panel) == ["CPD-7", "CPD-8"]


def test_disable_remove_and_follow_column_changes(panel) -> None:
    card = panel.add_filter("pIC50")
    card.set_bounds(8, None)
    panel.apply_now()
    assert len(panel.table.view()) == 3  # 8, 9 and the blank
    card.set_enabled(False)
    panel.apply_now()
    assert len(panel.table.view()) == 8
    card.set_enabled(True)
    panel.apply_now()
    panel.table.rename_column("pIC50", "Potency")
    assert card.column == "Potency" and len(panel.table.view()) == 3
    panel.table.remove_column("Potency")
    assert panel.cards == [] and len(panel.table.view()) == 8


def test_cards_refresh_when_data_changes(qtbot, panel) -> None:
    card = panel.add_filter("Series")
    panel.table.set_cells("Series", [5], ["D"])
    qtbot.waitUntil(lambda: "D" in card.listed_values(), timeout=2000)


def test_header_menu_adds_a_filter(window) -> None:
    from .test_column_actions import trigger

    trigger(window.columns.build_menu("pIC50", [], 1), "Filter")
    assert [c.column for c in window.filters.cards] == ["pIC50"]
    assert window.filter_dock.isVisibleTo(window)


def test_panel_tools_sit_in_the_title_bar_and_cards_carry_no_text(window) -> None:
    from PySide6.QtWidgets import QLabel

    panel = window._filters
    header = window.filter_dock.titleBarWidget()
    assert header.toolbar is panel.toolbar  # Add, Sets and Remove All are header icons
    assert not panel.clear_button.isEnabled()  # nothing to remove yet
    card = panel.add_filter("pIC50")
    assert panel.clear_button.isEnabled() and "Showing 3 of 3 rows" in panel.clear_button.toolTip()
    # no helper or summary text anywhere in the pane: the card's details are in its tooltip
    labels = [w for w in panel.findChildren(QLabel) if w.isVisibleTo(panel) and w is not card.title]
    assert [w.text() for w in labels if w.text() not in ("to",)] == []
    assert "3 numbers from 5 to 7" in card.title.toolTip()
    old_toolbar = panel.toolbar
    window._install(type(window.workspace)(window.workspace.table))  # a fresh document
    assert header.toolbar is window._filters.toolbar is not old_toolbar
    assert old_toolbar.parent() is None


def test_cards_are_striped_compact_and_use_glyphs(panel) -> None:
    cards = [panel.add_filter(c) for c in ("pIC50", "Series", "ID")]
    assert [c._alternate for c in cards] == [False, True, False]
    background = [c.styleSheet().split("background: ")[1][:7] for c in cards]
    assert background[0] == background[2] != background[1]
    panel.remove_filter(cards[0])  # the stripes follow the new order
    assert [c._alternate for c in panel.cards] == [False, True]
    series = cards[1]
    for button in (series.not_button, series.remove_button, series.all_button):
        assert not button.icon().isNull() and not button.text()
    # The list is as tall as its four values (A, B, C, blank), not a fixed 180 px.
    assert series.values.height() < 4 * series.values.sizeHintForRow(0) + 20
    series.none_button.click()
    assert series.checked_values() == set()
    series.all_button.click()
    assert len(series.checked_values()) == 4


def test_add_by_type_then_choose_the_column(panel) -> None:
    from chemsuite.app.filters.cards import CategoryCard, RangeCard, TextCard, sources_for

    panel._fill_add_menu()
    titles = [a.text() for a in panel.add_menu.actions()]
    assert titles == ["Number Range", "Category", "Text Match"]  # no searcher: no structures
    assert sources_for(panel.table, "range") == ["pIC50"]
    card = panel.add_filter_of_type("category")
    assert isinstance(card, CategoryCard) and card.column == "Series"  # not the ID column
    card = panel.change_source(card, "ID")
    assert [card.title.itemText(i) for i in range(card.title.count())] == ["ID", "pIC50", "Series"]
    card.not_button.setChecked(True)
    other = panel.add_filter_of_type("text")
    new = panel.change_source(card, "Series")
    assert isinstance(new, CategoryCard) and new.column == "Series" and new.negated
    assert panel.cards == [new, other] and new.title.currentText() == "Series"
    assert [c._alternate for c in panel.cards] == [False, True]
    new.set_checked({"A"})
    panel.apply_now()
    assert "CPD-1" not in shown_ids(panel) and "CPD-2" in shown_ids(panel)  # Not A
    new.title.activated.emit(new.title.findText("pIC50"))  # picked in the card itself
    assert isinstance(panel.cards[0], CategoryCard) and panel.cards[0].column == "pIC50"
    assert isinstance(other, TextCard) and not isinstance(panel.cards[0], RangeCard)
    panel.table.add_column("Notes", ["x"] * len(panel.table))
    assert "Notes" in [other.title.itemText(i) for i in range(other.title.count())]
