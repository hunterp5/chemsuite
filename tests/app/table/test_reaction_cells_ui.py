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

"""Reaction cells in the window (REACTIONS_PLAN R3): drawing, the sketcher, copy and paste."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QMimeData, QSize, Qt
from PySide6.QtGui import QGuiApplication

from chemsuite.app.table.table_model import KIND_ROLE, STATUS_ROLE
from chemsuite.core.reactions import reaction_from_text
from chemsuite.core.table.columns import REACTION
from chemsuite.io.rdf import rxn_block
from chemsuite.sketch import convert, reaction
from tests.gui import gui

pytestmark = gui

ESTER = "CC(=O)O.OCC>[H+]>CC(=O)OCC"
AMIDE = "[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]"


@pytest.fixture
def reactions(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": ["A", "B", "C"], "Rxn": [ESTER, AMIDE, "not>>a reaction"]})
    table.convert_column("Rxn", REACTION)
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _index(window, row: int):
    model = window.workspace.model
    return model.index(row, model.column_index("Rxn"))


def _sketcher(window):
    from chemsuite.app.sketcher import SketcherDialog

    return window.findChild(SketcherDialog)


def test_reaction_cells_report_kind_and_status(reactions) -> None:
    assert _index(reactions, 0).data(KIND_ROLE) == REACTION
    assert _index(reactions, 0).data(STATUS_ROLE) in ("ok", "pending")
    reactions.workspace.table.column("Rxn").reaction(2)  # parse the bad one
    assert _index(reactions, 2).data(STATUS_ROLE) == "invalid"
    assert _index(reactions, 2).data(0x03).startswith("Invalid reaction")  # ToolTipRole


def test_reaction_cells_draw_as_pictures(reactions) -> None:
    from chemsuite.app.table.table_view import StructureImages

    images = StructureImages()
    column = reactions.workspace.table.column("Rxn")
    image = images.image(column, 0, QSize(400, 120), 1.0, False)
    assert image is not None and not image.isNull()
    assert images.image(column, 2, QSize(400, 120), 1.0, False) is None


def test_double_click_opens_the_sketcher_and_ok_writes_one_undo_step(reactions) -> None:
    window = reactions
    table = window.workspace.table
    oid = window.workspace.model.oid_at(0)
    view = window.view
    assert view.edit(_index(window, 0), view.EditTrigger.DoubleClicked, None) is False
    sketcher = _sketcher(window)
    assert sketcher is not None and sketcher.reaction_only
    assert reaction.to_smiles(sketcher.canvas.sketch).count(">") == 2
    sketcher.canvas.set_sketch(reaction.from_text("C=C>>CC"), "Paste")
    sketcher.accept()
    assert table.cell("Rxn", oid) == "C=C>>CC"
    assert table.column("Rxn").reaction(0) is not None  # stored parsed
    assert table.history.undo_label == "Draw Rxn"
    table.undo()
    assert table.cell("Rxn", oid) == ESTER


def test_f2_opens_the_sketcher_but_typing_edits_in_place(reactions, qtbot) -> None:
    window = reactions
    view = window.view
    assert view.edit(_index(window, 1), view.EditTrigger.EditKeyPressed, None) is False
    assert _sketcher(window) is not None
    _sketcher(window).reject()
    view.setCurrentIndex(_index(window, 1))
    qtbot.keyClick(view, Qt.Key.Key_C)
    assert view.state() == view.State.EditingState  # an inline editor, not a dialog


def test_a_reaction_cell_does_not_take_a_lone_structure(reactions) -> None:
    window = reactions
    dialog = window.draw_reaction(window.workspace.model.oid_at(0), "Rxn")
    dialog.canvas.set_sketch(convert.from_text("CCO"), "Paste")
    dialog.accept()
    assert dialog.isVisible() and "arrow" in dialog.info.text()
    dialog.reject()
    assert window.workspace.table.cell("Rxn", window.workspace.model.oid_at(0)) == ESTER


def test_invalid_reaction_text_opens_an_inline_editor(reactions) -> None:
    window = reactions
    assert window.draw_reaction(window.workspace.model.oid_at(2), "Rxn") is None
    assert _sketcher(window) is None
    assert window.view.state() == window.view.State.EditingState


def test_cell_menu_offers_sketcher_and_reaction_copy_formats(reactions) -> None:
    menu = reactions.cell_menu(_index(reactions, 0))
    titles = [a.text() for a in menu.actions()]
    assert "Edit in Sketcher…" in titles and "Edit as Text…" not in titles
    copy_as = next(a.menu() for a in menu.actions() if a.text() == "Copy As")
    assert [a.text() for a in copy_as.actions()] == ["RXN File", "Reaction SMILES", "Image"]


def test_copy_as_rxn_and_reaction_smiles(reactions) -> None:
    window = reactions
    assert "Copied reaction as RXN File" in window.copy_reactions_as(_index(window, 1), "rxn")
    mime = QGuiApplication.clipboard().mimeData()
    assert bytes(mime.data("chemical/x-mdl-rxnfile").data()).startswith(b"$RXN")
    window.copy_reactions_as(_index(window, 0), "reaction-smiles")
    assert QGuiApplication.clipboard().text() == "CC(=O)O.CCO>[H+]>CCOC(C)=O"  # canonical
    assert "image" in window.copy_reactions_as(_index(window, 0), "image")


def test_paste_an_rxn_file_into_a_reaction_cell(reactions) -> None:
    window = reactions
    table = window.workspace.table
    mime = QMimeData()
    mime.setData("chemical/x-mdl-rxnfile", rxn_block(reaction_from_text(AMIDE)).encode())
    QGuiApplication.clipboard().setMimeData(mime)
    window.view.setCurrentIndex(_index(window, 2))
    assert "1 cell" in window.view.paste()
    assert table.column("Rxn").reaction(2) is not None
    model = window.workspace.model
    window.view.setCurrentIndex(model.index(0, model.column_index("ID")))
    window.view.clearSelection()
    assert "Select a reaction cell" in window.view.paste()


def test_column_type_menu_converts_to_reaction(window) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"Text": [ESTER]}, undoable=False)
    window._install(Workspace(table))
    window._columns.convert("Text", REACTION)
    assert table.kind("Text") == REACTION
