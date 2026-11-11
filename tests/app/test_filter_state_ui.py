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
"""Filters saved with the document, named filter sets, and copying shown rows."""

from __future__ import annotations

import pytest

from .conftest import gui

pytestmark = gui


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["C0", "C1", "C2", "C3", "C4"],
            "SMILES": ["c1ccccc1O", "CCO", "Cc1ccccc1", "c1ccncc1", "CCCC"],
            "pIC50": ["7.5", "6.0", "8.2", "5.1", ""],
            "Series": ["A", "B", "A", "B", "C"],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.convert_column("pIC50", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def shown(window) -> list[str]:
    table = window.workspace.table
    return [table.cell("ID", int(o)) for o in table.oids[table.view()]]


def add_cards(window) -> None:
    panel = window.filters
    panel.add_filter("pIC50").set_bounds(6.0, None)  # range; blanks are kept by default
    panel.add_filter("Series").set_checked({"A", "C"})  # category
    text = panel.add_filter("ID", "text")
    text.query.setText("C[0-3]")
    text.mode.setCurrentIndex(text.mode.findData("regex"))
    text.not_button.setChecked(True)
    panel.apply_now()


def test_filters_are_saved_with_the_document(chem, tmp_path, qtbot) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    add_cards(chem)
    assert shown(chem) == ["C4"]
    structure = chem.filters.add_filter("SMILES")
    structure.query.set_text("c1ccccc1")
    structure.enabled_box.setChecked(False)
    chem.filter_dock.show()
    target = tmp_path / "filters.csd"
    assert chem._save_to(target)

    chem._install(Workspace(Table()))
    assert chem.filters.cards == []
    chem.open_path(target)
    qtbot.waitUntil(lambda: chem.workspace.path == target, timeout=10_000)
    cards = chem.filters.cards
    assert [c.kind for c in cards] == ["range", "category", "text", "structure"]
    assert cards[0].bounds() == (6.0, None) and cards[0].blanks.isChecked()
    assert cards[1].checked_values() == {"A", "C"}
    assert cards[2].negated and cards[2].mode.currentData() == "regex"
    assert not cards[3].is_enabled and cards[3].query.text() == "c1ccccc1"
    assert shown(chem) == ["C4"]
    assert chem.filter_dock.isVisibleTo(chem)  # open, as when saved
    assert not chem.isWindowModified()


def test_filter_sets(chem) -> None:
    panel = chem.filters
    add_cards(chem)
    panel.save_set("Actives")
    panel.remove_all()
    panel.add_filter("Series").set_checked({"B"})
    panel.apply_now()
    panel.save_set("Series B")
    assert shown(chem) == ["C1", "C3"]
    assert panel.apply_set("Actives") == 0
    assert shown(chem) == ["C4"] and len(panel.cards) == 3
    assert panel.apply_set("Series B") == 0
    assert shown(chem) == ["C1", "C3"]
    # A set whose column is gone restores what it can.
    chem.workspace.table.remove_columns(["Series"])
    assert panel.apply_set("Actives") == 1 and len(panel.cards) == 2
    panel.delete_set("Series B")
    assert list(panel.sets) == ["Actives"]
    assert panel.state()["sets"]["Actives"][0]["column"] == "pIC50"


def test_copy_shown_rows_to_a_new_document(chem) -> None:
    add_cards(chem)
    chem.workspace.table.set_sort([("ID", False)])
    chem.filters.cards[0].set_bounds(None, None)
    chem.filters.cards[2].enabled_box.setChecked(False)
    chem.filters.apply_now()
    assert shown(chem) == ["C4", "C2", "C0"]
    new = chem.copy_shown_rows_to_new_document()
    try:
        table = new.workspace.table
        assert [table.cell("ID", int(o)) for o in table.oids] == ["C4", "C2", "C0"]
        assert table.column_names == ["ID", "SMILES", "pIC50", "Series"]
        assert table.kind("SMILES") == "structure" and table.kind("pIC50") == "number"
        assert not table.filters and new.workspace.is_modified  # never saved: ask on close
        assert new.document_name == "Untitled"
    finally:
        new.workspace.mark_saved(table.history.checkpoint(), new.workspace.lists.revision)
        new.close()
