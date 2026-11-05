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

SMILES = ["c1ccccc1O", "CCO", "Cc1ccccc1", "c1ccncc1", "CCCC", ""]
AROMATIC = [0, 2]  # rows with a benzene ring


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": [f"C{i}" for i in range(len(SMILES))], "SMILES": SMILES},
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _oids(window, rows):
    return [int(window.workspace.table.oids[r]) for r in rows]


def _shown_ids(window):
    table = window.workspace.table
    return [table.cell("ID", int(o)) for o in table.oids[table.view()]]


def _wait(qtbot, job):
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)


def test_dialog_selects_matches_and_highlights(chem, qtbot) -> None:
    dialog = chem.structure_search()
    dialog.query.set_text("c1ccccc1")
    dialog.actions["select"].setChecked(True)
    assert dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()
    assert dialog.start()
    _wait(qtbot, dialog.job)
    assert sorted(chem.workspace.selection.oids) == _oids(chem, AROMATIC)
    assert "SMILES" in chem.view.highlights()
    chem.clear_structure_highlights()
    assert not chem.view.highlights()
    # The next dialog starts from the last query.
    again = chem.structure_search()
    assert again.query.text() == "c1ccccc1"
    again.reject()
    dialog.reject()


def test_dialog_rejects_bad_queries(chem) -> None:
    dialog = chem.structure_search()
    ok = dialog.buttons.button(dialog.buttons.StandardButton.Ok)
    dialog.query.set_text("C((((")
    assert not ok.isEnabled()
    assert "SMILES or SMARTS" in dialog.query_info.text()
    dialog.query.set_text("[C,N]")  # SMARTS: fine for substructures only
    assert ok.isEnabled()
    dialog.mode.setCurrentIndex(dialog.mode.findData("exact"))
    assert not ok.isEnabled()
    dialog.reject()


def test_dialog_adds_matches_to_a_list_within_the_selection(chem, qtbot) -> None:
    chem.workspace.selection.set(_oids(chem, [0, 1, 3]), "test")
    dialog = chem.structure_search()
    dialog.query.set_text("c1ccccc1")
    dialog.actions["list"].setChecked(True)
    dialog.list_name.setCurrentText("Aromatics")
    dialog.rows.setCurrentIndex(dialog.rows.findData("selected"))
    assert dialog.start()
    _wait(qtbot, dialog.job)
    assert chem.workspace.lists.oids("Aromatics").tolist() == _oids(chem, [0])
    dialog.reject()


def test_dialog_show_only_adds_a_structure_filter(chem, qtbot) -> None:
    dialog = chem.structure_search()
    dialog.query.set_text("c1ccccc1")
    assert dialog.action() == "show"
    assert dialog.start()
    card = chem.filters.cards[0]
    assert card.kind == "structure"
    _wait(qtbot, card.job or dialog.job)
    qtbot.waitUntil(lambda: _shown_ids(chem) == ["C0", "C2"], timeout=5000)
    dialog.reject()


def test_structure_card(chem, qtbot) -> None:
    card = chem.filter_column("SMILES")
    assert card.kind == "structure"
    card.query.set_text("c1ccccc1")
    card.search_now()
    _wait(qtbot, card.job)
    qtbot.waitUntil(lambda: _shown_ids(chem) == ["C0", "C2"], timeout=5000)
    assert card.status.text() == "2 match(es)"
    assert "SMILES" in chem.view.highlights()

    card.not_button.setChecked(True)
    qtbot.waitUntil(lambda: _shown_ids(chem) == ["C1", "C3", "C4", "C5"], timeout=5000)
    card.not_button.setChecked(False)

    # Editing the column searches again.
    table = chem.workspace.table
    table.set_cells("SMILES", _oids(chem, [4]), ["c1ccccc1CC"])
    qtbot.waitUntil(lambda: _shown_ids(chem) == ["C0", "C2", "C4"], timeout=20_000)

    # Similarity mode.
    card.set_query(
        card.current_query().__class__("c1ccccc1O", "similarity", threshold=0.99),
    )
    _wait(qtbot, card.job)
    qtbot.waitUntil(lambda: _shown_ids(chem) == ["C0"], timeout=5000)
    assert "SMILES" not in chem.view.highlights()  # nothing to highlight for similarity

    # An invalid query stops filtering and says why.
    card.mode.setCurrentIndex(card.mode.findData("substructure"))
    card.query.set_text("C((((")
    card.search_now()
    assert "SMILES or SMARTS" in card.status.text()
    qtbot.waitUntil(lambda: len(_shown_ids(chem)) == len(SMILES), timeout=5000)

    card.query.set_text("c1ccccc1")
    card.search_now()
    _wait(qtbot, card.job)
    qtbot.waitUntil(lambda: "SMILES" in chem.view.highlights(), timeout=5000)
    chem.filters.remove_filter(card)
    assert not chem.view.highlights()
    assert len(_shown_ids(chem)) == len(SMILES)


def test_highlighted_drawing(chem) -> None:
    from PySide6.QtCore import QSize

    from chemsuite.analysis.structure_search import StructureQuery
    from chemsuite.app.structure_search import highlight_for

    table = chem.workspace.table
    column = table.column("SMILES")
    images = chem._images
    plain = images.image(column, 0, QSize(120, 90), 1.0, False)
    marked = images.image(
        column, 0, QSize(120, 90), 1.0, False, highlight_for(StructureQuery("c1ccccc1"))
    )
    assert plain is not None and marked is not None
    assert plain != marked
