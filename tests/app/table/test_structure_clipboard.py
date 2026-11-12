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
from PySide6.QtCore import QMimeData
from PySide6.QtGui import QGuiApplication
from rdkit import Chem

from chemsuite.core.structure_input import molblock_of
from chemsuite.sketch import convert
from tests.gui import gui

pytestmark = gui

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"


@pytest.fixture
def structures(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": ["A", "B"], "SMILES": ["CCO", "c1ccccc1"]}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def set_clipboard_molfile(smiles: str) -> str:
    block = molblock_of(Chem.MolFromSmiles(smiles))
    mime = QMimeData()
    mime.setData("chemical/x-mdl-molfile", block.encode())
    mime.setText("ignored plain text")
    QGuiApplication.clipboard().setMimeData(mime)
    return block


def test_clipboard_reads_a_molfile_before_plain_text(qapp) -> None:
    from chemsuite.app.table.structure_clipboard import clipboard_structure_text

    block = set_clipboard_molfile(ASPIRIN_SMILES)
    assert clipboard_structure_text() == block


def test_mdlct_clipboard_format() -> None:
    from chemsuite.app.table.structure_clipboard import _mdlct_to_text

    lines = ["", "  ChemDraw", "", "  0  0  0  0  0  0  0  0  0  0999 V2000", "M  END"]
    data = b"".join(bytes([len(line)]) + line.encode() for line in lines)
    assert _mdlct_to_text(data) == "\n".join(lines)


def _sketcher(window):
    from chemsuite.app.sketcher import SketcherDialog

    return window.findChild(SketcherDialog)


def test_f2_opens_sketcher_and_ok_writes_undoable_edit(structures) -> None:
    window = structures
    table = window.workspace.table
    model = window.workspace.model
    index = model.index(0, model.column_index("SMILES"))
    opened = []
    window.view.structureEditRequested.connect(lambda i, t: opened.append((i.row(), t)))
    assert window.view.edit(index, window.view.EditTrigger.EditKeyPressed, None) is False
    assert opened == [(0, "")]
    sketcher = _sketcher(window)
    assert sketcher is not None and sketcher.isVisible()
    assert convert.to_smiles(sketcher.canvas.history.sketch) == "CCO"
    sketcher.canvas.set_sketch(convert.from_text(ASPIRIN_SMILES), "Paste")
    sketcher.accept()
    assert table.cell("SMILES", model.oid_at(0)) == Chem.CanonSmiles(ASPIRIN_SMILES)
    assert table.history.undo_label == "Draw SMILES"
    assert table.column("SMILES").mol(0) is not None  # stored parsed, no background job
    table.undo()
    assert table.cell("SMILES", model.oid_at(0)) == "CCO"


def test_typing_in_a_structure_cell_opens_the_sketcher(structures, qtbot) -> None:
    window = structures
    model = window.workspace.model
    window.show()
    qtbot.waitExposed(window)
    window.view.setFocus()
    window.view.setCurrentIndex(model.index(1, model.column_index("SMILES")))
    qtbot.keyClicks(window.view, "N")
    sketcher = _sketcher(window)
    assert sketcher is not None and convert.to_smiles(sketcher.canvas.sketch) == "c1ccccc1"
    sketcher.reject()


def test_cancel_changes_nothing(structures) -> None:
    window = structures
    model = window.workspace.model
    window.draw_structure(model.oid_at(0), "SMILES")
    sketcher = _sketcher(window)
    assert convert.to_smiles(sketcher.canvas.history.sketch) == "CCO"
    sketcher.canvas.set_sketch(convert.from_text("CCN"), "Paste")
    sketcher.reject()
    assert window.workspace.table.cell("SMILES", model.oid_at(0)) == "CCO"
    assert not window.isWindowModified()


def test_add_compound_appends_a_row(structures) -> None:
    window = structures
    table = window.workspace.table
    sketcher = window.add_compound()
    sketcher.canvas.set_sketch(convert.from_text("CCN"), "Paste")
    sketcher.accept()
    assert len(table) == 3
    new = table.oids[-1]
    assert table.cell("SMILES", int(new)) == "CCN"
    assert window.workspace.selection.oids == {int(new)}


def test_add_compound_to_a_table_without_structures(window) -> None:
    table = window.workspace.table
    sketcher = window.add_compound()
    sketcher.canvas.set_sketch(convert.from_text("CCN"), "Paste")
    sketcher.accept()
    assert table.kind("Structure") == "structure"
    assert len(table) == 4


def test_paste_molfile_into_structure_cell(structures) -> None:
    window = structures
    table = window.workspace.table
    model = window.workspace.model
    set_clipboard_molfile(ASPIRIN_SMILES)
    window.view.setCurrentIndex(model.index(1, model.column_index("SMILES")))
    message = window.view.paste()
    assert "1 cell" in message
    assert table.cell("SMILES", model.oid_at(1)) == Chem.MolToSmiles(
        Chem.MolFromSmiles(ASPIRIN_SMILES)
    )
    window.view.setCurrentIndex(model.index(0, model.column_index("ID")))
    window.view.clearSelection()
    assert "Select a structure cell" in window.view.paste()


def test_copy_structure_as_molfile(structures) -> None:
    window = structures
    model = window.workspace.model
    window.view.setCurrentIndex(model.index(1, model.column_index("SMILES")))
    window.copy_structure()
    mime = QGuiApplication.clipboard().mimeData()
    assert b"M  END" in bytes(mime.data("chemical/x-mdl-molfile").data())
    assert mime.text() == "c1ccccc1"
