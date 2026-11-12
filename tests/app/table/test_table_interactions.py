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

"""Table interactions: header selection, column delete, structure double-click and menus,
and the column menu's Duplicate, Logarithmic, Precision and Colour."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QMessageBox
from rdkit import Chem

from chemsuite.core.table.style import CATEGORIES, GRADIENT
from tests.gui import gui

from .test_column_actions import menu_titles, trigger

pytestmark = gui


@pytest.fixture
def structures(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["A", "B", "C"], "SMILES": ["CCO", "c1ccccc1", "C1CC"], "IC50": ["1", "10", "100"]},
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def click_header(qtbot, view, column: int, modifier=Qt.KeyboardModifier.NoModifier) -> None:
    header = view.horizontalHeader()
    x = header.sectionViewportPosition(column) + header.sectionSize(column) // 2
    qtbot.mouseClick(
        header.viewport(), Qt.MouseButton.LeftButton, modifier, QPoint(x, header.height() // 2)
    )


# -- header selection ---------------------------------------------------------------------------


def test_header_click_selects_the_column_without_sorting(qtbot, window) -> None:
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    click_header(qtbot, view, 1)
    assert view.whole_columns() == ["pIC50"]
    assert window.workspace.table.sort_keys == ()
    # Whole columns are not a compound selection.
    assert not window.workspace.selection.oids
    assert len(view.selected_oids()) == 0


def test_shift_and_ctrl_click_select_several_columns(qtbot, structures) -> None:
    window = structures
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    click_header(qtbot, view, 0)
    click_header(qtbot, view, 2, Qt.KeyboardModifier.ShiftModifier)
    assert view.whole_columns() == ["ID", "SMILES", "IC50"]
    click_header(qtbot, view, 1, Qt.KeyboardModifier.ControlModifier)
    assert view.whole_columns() == ["ID", "IC50"]


def test_delete_key_deletes_selected_columns_in_one_step(qtbot, structures) -> None:
    window = structures
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    table = window.workspace.table
    click_header(qtbot, view, 0)
    click_header(qtbot, view, 2, Qt.KeyboardModifier.ControlModifier)
    view.setFocus()
    qtbot.keyClick(view, Qt.Key.Key_Delete)
    assert table.column_names == ["SMILES"]
    table.undo()
    assert table.column_names == ["ID", "SMILES", "IC50"]


def test_delete_key_on_cells_still_clears_them(qtbot, window) -> None:
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    model = window.workspace.model
    view.setFocus()
    view.setCurrentIndex(model.index(0, 1))
    view.selectionModel().select(model.index(0, 1), view.selectionModel().SelectionFlag.Select)
    qtbot.keyClick(view, Qt.Key.Key_Delete)
    assert window.workspace.table.column_names == ["ID", "pIC50"]
    assert window.workspace.table.cell("pIC50", model.oid_at(0)) == ""


def test_select_all_after_a_column_click_is_not_a_column_selection(qtbot, window) -> None:
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    click_header(qtbot, view, 0)
    view.setFocus()
    qtbot.keyClick(view, Qt.Key.Key_A, Qt.KeyboardModifier.ControlModifier)
    assert view.whole_columns() == []
    qtbot.keyClick(view, Qt.Key.Key_Delete)
    assert window.workspace.table.column_names == ["ID", "pIC50"]  # cells cleared, not columns


def test_sorting_from_the_menu_keeps_a_column_selection(qtbot, window) -> None:
    window.show()
    qtbot.waitExposed(window)
    click_header(qtbot, window.view, 1)
    trigger(window.columns.build_menu("pIC50", ["pIC50"], 1), "Descending")
    assert window.workspace.table.sort_keys == (("pIC50", False),)
    assert window.view.whole_columns() == ["pIC50"]


# -- structures -----------------------------------------------------------------------------------


def test_double_click_on_a_structure_opens_the_sketcher(structures) -> None:
    from chemsuite.app.sketcher import SketcherDialog

    window = structures
    table = window.workspace.table
    model = window.workspace.model
    index = model.index(0, model.column_index("SMILES"))
    assert window.view.edit(index, window.view.EditTrigger.DoubleClicked, None) is False
    sketcher = window.findChild(SketcherDialog)
    assert sketcher is not None and sketcher.isVisible()
    assert sketcher.canvas.sketch.atoms  # ethanol is on the canvas
    sketcher.accept()
    assert table.history.undo_label != "Draw SMILES"  # OK without a change writes nothing

    from chemsuite.sketch.convert import from_text

    window.view.edit(index, window.view.EditTrigger.DoubleClicked, None)
    sketcher = window.findChild(SketcherDialog)
    sketcher.canvas.set_sketch(from_text("CCN"), "Replace")
    sketcher.accept()
    assert table.cell("SMILES", model.oid_at(0)) == "CCN"
    assert table.column("SMILES").mol(0).GetNumConformers() == 1  # the drawing's layout
    assert table.history.undo_label == "Draw SMILES"
    table.undo()
    assert table.cell("SMILES", model.oid_at(0)) == "CCO"


def test_double_click_on_an_invalid_structure_opens_an_empty_sketcher(structures) -> None:
    from chemsuite.app.sketcher import SketcherDialog

    window = structures
    model = window.workspace.model
    table = window.workspace.table
    index = model.index(2, model.column_index("SMILES"))  # "C1CC" does not parse
    window.view.edit(index, window.view.EditTrigger.DoubleClicked, None)
    sketcher = window.findChild(SketcherDialog)
    assert sketcher is not None and sketcher.canvas.sketch.is_empty()
    sketcher.reject()
    assert table.cell("SMILES", model.oid_at(2)) == "C1CC"  # kept until a drawing replaces it


def test_structure_cell_menu(structures) -> None:
    window = structures
    model = window.workspace.model
    menu = window.cell_menu(model.index(0, model.column_index("SMILES")))
    assert menu_titles(menu) == [  # no plain Copy, no text editor
        "Browse",
        "Edit in Sketcher…",
        "Copy As",
        "Paste",
        "Clear Cells",
    ]
    copy_as = next(a.menu() for a in menu.actions() if a.text() == "Copy As")
    assert menu_titles(copy_as) == ["Molfile", "SMILES", "InChI", "InChIKey", "Image"]
    assert menu_titles(window.cell_menu(model.index(0, 0))) == [
        "Browse",
        "Copy",
        "Paste",
        "Clear Cells",
    ]


@pytest.mark.parametrize(
    ("fmt", "expected"),
    [
        ("smiles", "CCO"),
        ("inchi", "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
        ("inchikey", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"),
    ],
)
def test_copy_as_text_formats(structures, fmt, expected) -> None:
    window = structures
    model = window.workspace.model
    window.copy_structures_as(model.index(0, model.column_index("SMILES")), fmt)
    assert QGuiApplication.clipboard().text() == expected


def test_copy_as_molfile_offers_the_drawing_format(structures) -> None:
    window = structures
    model = window.workspace.model
    window.copy_structures_as(model.index(1, model.column_index("SMILES")), "molfile")
    mime = QGuiApplication.clipboard().mimeData()
    block = bytes(mime.data("chemical/x-mdl-molfile").data()).decode()
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(block)) == "c1ccccc1"


def test_copy_as_takes_the_selected_structures_in_shown_order(structures) -> None:
    window = structures
    view = window.view
    model = window.workspace.model
    col = model.column_index("SMILES")
    flags = view.selectionModel().SelectionFlag.Select
    view.selectionModel().select(model.index(1, col), flags)
    view.selectionModel().select(model.index(0, col), flags)
    message = window.copy_structures_as(model.index(0, col), "smiles")
    assert QGuiApplication.clipboard().text() == "CCO\nc1ccccc1"
    assert message == "Copied 2 structures as SMILES"
    window.copy_structures_as(model.index(0, col), "molfile")
    text = QGuiApplication.clipboard().text()
    assert text.count("$$$$") == 2  # an SD file
    assert not QGuiApplication.clipboard().mimeData().hasFormat("chemical/x-mdl-molfile")


def test_copy_as_image(structures) -> None:
    window = structures
    model = window.workspace.model
    window.copy_structures_as(model.index(0, model.column_index("SMILES")), "image")
    assert not QGuiApplication.clipboard().image().isNull()


def test_copy_as_reports_an_invalid_structure(structures) -> None:
    window = structures
    model = window.workspace.model
    message = window.copy_structures_as(model.index(2, model.column_index("SMILES")), "smiles")
    assert message == "No valid structure to copy"


# -- column menu ----------------------------------------------------------------------------------


def test_menu_has_the_new_column_items(window) -> None:
    titles = menu_titles(window.columns.build_menu("pIC50", [], 1))
    for title in ("Duplicate Column", "Logarithmic (log₁₀)", "Number Format…", "Colour Column…"):
        assert title in titles


def test_duplicate_column(window) -> None:
    table = window.workspace.table
    trigger(window.columns.build_menu("pIC50", [], 1), "Duplicate Column")
    assert table.column_names == ["ID", "pIC50", "pIC50 copy"]


def test_logarithmic_toggles_and_marks_the_header(structures) -> None:
    window = structures
    table = window.workspace.table
    model = window.workspace.model
    col = model.column_index("IC50")
    menu = window.columns.build_menu("IC50", [], col)
    trigger(menu, "Logarithmic (log₁₀)")
    assert [table.cell("IC50", oid) for oid in model.visible_oids] == ["0", "1", "2"]
    assert model.headerData(col, Qt.Orientation.Horizontal) == "IC50 (log₁₀)"
    menu = window.columns.build_menu("IC50", [], col)
    action = next(a for a in menu.actions() if a.text() == "Logarithmic (log₁₀)")
    assert action.isChecked()
    action.trigger()
    assert [table.cell("IC50", oid) for oid in model.visible_oids] == ["1", "10", "100"]
    assert model.headerData(col, Qt.Orientation.Horizontal) == "IC50"


def test_logarithmic_is_disabled_with_a_reason(window, monkeypatch) -> None:
    table = window.workspace.table
    table.set_cells("pIC50", [table.oids[0]], ["-1"])
    menu = window.columns.build_menu("pIC50", [], 1)
    action = next(a for a in menu.actions() if a.text() == "Logarithmic (log₁₀)")
    assert not action.isEnabled() and "above zero" in action.toolTip()
    menu = window.columns.build_menu("ID", [], 0)
    action = next(a for a in menu.actions() if a.text() == "Logarithmic (log₁₀)")
    assert not action.isEnabled() and "Number" in action.toolTip()
    warned = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: warned.append(a[2]))
    assert not window.columns.toggle_log("pIC50")
    assert warned


def test_number_format_is_display_only_and_undoable(window) -> None:
    table = window.workspace.table
    model = window.workspace.model
    table.set_cells("pIC50", [int(model.visible_oids[0])], ["1234.5678"])
    shown = lambda: [model.data(model.index(r, 1)) for r in range(3)]  # noqa: E731
    dialog = window.columns.number_format("pIC50")
    dialog.digits.setValue(2)  # Fixed decimals
    dialog.format.setCurrentIndex(dialog.format.findData("fixed"))
    dialog.digits.setValue(2)
    dialog.accept()
    assert shown() == ["1234.57", "6.00", "7.00"]
    assert table.column("pIC50").values().tolist() == [1234.5678, 6.0, 7.0]

    dialog = window.columns.number_format("pIC50")
    assert dialog.chosen() == "fixed" and dialog.digits.value() == 2
    dialog.format.setCurrentIndex(dialog.format.findData("significant"))
    assert dialog.digits_label.text() == "Significant figures:" and dialog.digits.value() == 3
    dialog.accept()
    assert shown() == ["1230", "6.00", "7.00"]
    # editing starts from the number as typed, not the formatted text
    assert model.data(model.index(0, 1), Qt.ItemDataRole.EditRole) == "1234.5678"

    dialog = window.columns.number_format("pIC50")
    dialog.format.setCurrentIndex(dialog.format.findData("thousands"))
    dialog.digits.setValue(1)
    dialog.accept()
    assert shown() == ["1,234.6", "6.0", "7.0"]
    header = model.headerData(1, Qt.Orientation.Horizontal, Qt.ItemDataRole.ToolTipRole)
    assert "Shown as Thousands separator (1 decimals)" in header
    assert "3 filled · 0 blank" in header and "Median 7" in header

    dialog = window.columns.number_format("pIC50")
    dialog.format.setCurrentIndex(dialog.format.findData("as_typed"))
    assert not dialog.digits.isEnabled()
    dialog.accept()
    assert shown() == ["1234.5678", "6", "7"]
    table.undo()
    assert shown() == ["1,234.6", "6.0", "7.0"]


def test_colour_dialog_gradient(window) -> None:
    from chemsuite.app.table.column_colour import ColumnColourDialog

    table = window.workspace.table
    model = window.workspace.model
    dialog = window.columns.colour("pIC50")
    assert isinstance(dialog, ColumnColourDialog)
    assert dialog.chosen_mode() == GRADIENT
    assert not dialog.status.text().strip()
    dialog.minimum.setText("9")
    assert "below high" in dialog.status.text() and not dialog.primary.isEnabled()
    dialog.minimum.setText("")
    dialog.accept()
    rule = table.column("pIC50").style.colour
    assert rule.mode == GRADIENT and rule.minimum is None
    low = model.data(model.index(0, 1), Qt.ItemDataRole.BackgroundRole)
    high = model.data(model.index(2, 1), Qt.ItemDataRole.BackgroundRole)
    assert low.rgb() & 0xFFFFFF == rule.low and high.rgb() & 0xFFFFFF == rule.high
    assert model.data(model.index(0, 0), Qt.ItemDataRole.BackgroundRole) is None
    table.undo()
    assert model.data(model.index(0, 1), Qt.ItemDataRole.BackgroundRole) is None


def test_colour_dialog_categories_for_text(window) -> None:
    table = window.workspace.table
    model = window.workspace.model
    dialog = window.columns.colour("ID")
    modes = [dialog.mode.itemData(i) for i in range(dialog.mode.count())]
    assert modes == ["off", CATEGORIES]  # no gradients for text
    dialog.mode.setCurrentIndex(1)
    dialog.accept()
    assert table.column("ID").style.colour.mode == CATEGORIES
    colours = {
        model.data(model.index(r, 0), Qt.ItemDataRole.BackgroundRole).rgba() for r in range(3)
    }
    assert len(colours) == 3


def test_colour_is_not_offered_for_structures(structures) -> None:
    model = structures.workspace.model
    col = model.column_index("SMILES")
    menu = structures.columns.build_menu("SMILES", [], col)
    action = next(a for a in menu.actions() if a.text() == "Colour Column…")
    assert not action.isEnabled()


# -- Tools ▸ Sketcher -------------------------------------------------------------------------


def test_chemistry_sketcher_adds_the_drawing_as_a_row(window) -> None:
    from PySide6.QtWidgets import QPushButton

    from chemsuite.app.sketcher import SketcherDialog
    from chemsuite.sketch.convert import from_text

    table = window.workspace.table
    dialog = window.open_tool("sketcher")
    assert isinstance(dialog, SketcherDialog) and dialog.canvas.sketch.is_empty()
    assert dialog.buttons.button(dialog.buttons.StandardButton.Ok).text() == "Add to Table"
    # Only OK and Cancel: undo, copy and the rest are on the keyboard and right-click menu.
    assert len(dialog.findChildren(QPushButton)) == 2
    assert not dialog.info.isVisible()
    # Elements run down the left of the canvas.
    assert dialog.toolbar.side.geometry().right() < dialog.canvas.geometry().left()
    dialog.canvas.set_sketch(from_text("c1ccccc1O"), "Draw")
    dialog.accept()
    assert len(table) == 4
    assert "Structure" in table.column_names
    assert Chem.MolToSmiles(table.column("Structure").mol(3)) == "Oc1ccccc1"
    assert window.workspace.selection.oids == {int(table.oids[3])}


def test_frozen_columns_stay_put_and_follow_the_table(window, qtbot) -> None:
    from chemsuite.app.document import session

    table = window.workspace.table
    table.add_rows({"ID": [f"X{i}" for i in range(60)], "pIC50": ["1"] * 60}, undoable=False)
    for i in range(6):
        table.add_column(f"Wide {i}", ["text"] * len(table))
    window.resize(700, 400)
    window.show()
    qtbot.waitExposed(window)
    view = window.view
    for col in range(view.table_model.columnCount()):
        view.setColumnWidth(col, 300)
    qtbot.waitUntil(lambda: view.horizontalScrollBar().maximum() > 600)
    titles = menu_titles(window.columns.build_menu("pIC50", [], 1))
    assert "Freeze Columns up to Here" in titles and "Unfreeze Columns" not in titles
    trigger(window.columns.build_menu("pIC50", [], 1), "Freeze Columns up to Here")
    frozen = view._frozen
    assert view.frozen_count == 2 and frozen.isVisible()
    assert frozen.width() == 2 * 300 + 1
    assert [frozen.isColumnHidden(c) for c in range(4)] == [False, False, True, True]

    view.verticalScrollBar().setValue(20)  # up and down: both move
    assert frozen.verticalScrollBar().value() == 20
    view.setColumnWidth(0, 120)  # a frozen column resized on the table
    assert frozen.columnWidth(0) == 120 and frozen.width() == 120 + 300 + 1
    last = view.table_model.columnCount() - 1
    view.scrollTo(view.table_model.index(25, last))
    view.scrollTo(view.table_model.index(25, 2))  # back left: not hidden under the frozen part
    assert view.columnViewportPosition(2) >= frozen.frozen_width()
    view.horizontalScrollBar().setValue(view.horizontalScrollBar().maximum())
    sideways = view.horizontalScrollBar().value()
    assert sideways > 0
    view.scrollTo(view.table_model.index(30, 0))  # a frozen cell: no sideways jump
    assert view.horizontalScrollBar().value() == sideways

    menus = []
    view.cellMenuRequested.disconnect()  # listen instead of showing the window's menu
    view.cellMenuRequested.connect(lambda index, _pos: menus.append(index.column()))
    from PySide6.QtGui import QContextMenuEvent

    at = frozen.visualRect(frozen.model().index(21, 0)).center()
    event = QContextMenuEvent(QContextMenuEvent.Reason.Mouse, at, frozen.mapToGlobal(at))
    frozen.contextMenuEvent(event)
    assert menus == [0]

    extras = session.capture(window)
    view.set_frozen(0)
    assert not frozen.isVisible()
    session.restore(window, extras)
    assert window.view.frozen_count == 2
    assert "Unfreeze Columns" in menu_titles(window.columns.build_menu("ID", [], 0))


# -- aligned drawings and previews --------------------------------------------------------------


def test_align_drawings_from_the_column_menu_and_saved_with_the_layout(structures) -> None:
    from chemsuite.app.document import session

    window = structures
    view = window.view
    assert "Align Drawings…" in menu_titles(window.columns.build_menu("SMILES", [], 1))
    assert "Align Drawings…" not in menu_titles(window.columns.build_menu("IC50", [], 2))
    assert "Stop Aligning Drawings" not in menu_titles(window.columns.build_menu("SMILES", [], 1))

    dialog = window.align_drawings("SMILES")
    ok = dialog.primary
    assert not ok.isEnabled()  # nothing typed yet
    dialog.core.set_text("not [[ a structure")
    assert not ok.isEnabled() and "not SMILES" in dialog.status.text()
    dialog.core.set_text("c1ccccc1")
    assert ok.isEnabled() and not dialog.footer.isVisibleTo(dialog)
    ok.click()
    assert view.alignments() == {"SMILES": "c1ccccc1"}
    assert "Stop Aligning Drawings" in menu_titles(window.columns.build_menu("SMILES", [], 1))

    images = view.itemDelegate()._images
    column = window.workspace.table.column("SMILES")
    alignment = view.itemDelegate().alignments["SMILES"]
    from PySide6.QtCore import QSize

    plain = images.image(column, 1, QSize(120, 90), 1.0, False)
    aligned = images.image(column, 1, QSize(120, 90), 1.0, False, None, alignment)
    assert plain is not aligned  # separate cache entries

    extras = session.capture(window)
    trigger(window.columns.build_menu("SMILES", [], 1), "Stop Aligning Drawings")
    assert view.alignments() == {}
    session.restore(window, extras)
    assert window.view.alignments() == {"SMILES": "c1ccccc1"}
    assert not window.view.set_alignment("SMILES", "not [[ a structure")


def test_structure_cells_preview_a_large_drawing_on_hover(structures) -> None:
    model = structures.view.table_model
    column = structures.workspace.table.column("SMILES")
    assert model.data(model.index(1, 1), Qt.ItemDataRole.ToolTipRole) == "c1ccccc1"  # unparsed
    for row in range(3):
        column.obj(int(model._slots[row]))  # what drawing the cells does
    tip = model.data(model.index(1, 1), Qt.ItemDataRole.ToolTipRole)
    assert tip.startswith('<img src="data:image/png;base64,')
    assert 'width="360" height="270"' in tip and "c1ccccc1" in tip
    assert model.data(model.index(2, 1), Qt.ItemDataRole.ToolTipRole).startswith("Invalid")


# -- the row header -----------------------------------------------------------------------------


def test_short_row_numbers() -> None:
    from chemsuite.app.table.table_model import short_row_number

    assert [short_row_number(n) for n in (7, 999, 1000, 1234, 12_345, 123_456, 2_500_000)] == [
        "7",
        "999",
        "1k",
        "1.2k",
        "12k",
        "123k",
        "2.5M",
    ]


def test_cell_selection_survives_sorting_and_row_header_selects_rows(window, qtbot) -> None:
    table = window.workspace.table
    view = window.view
    window.show()
    qtbot.waitExposed(window)
    model = view.table_model
    qtbot.mouseClick(
        view.viewport(),
        Qt.MouseButton.LeftButton,
        pos=view.visualRect(model.index(0, 1)).center(),  # A's pIC50
    )
    table.set_sort([("pIC50", False)])  # A moves to the bottom
    ranges = [(r.top(), r.left(), r.bottom(), r.right()) for r in view.selectionModel().selection()]
    assert ranges == [(2, 1, 2, 1)]  # still that one cell, not the whole row
    header = view.verticalHeader()
    qtbot.mouseClick(
        header.viewport(),
        Qt.MouseButton.LeftButton,
        pos=QPoint(header.width() // 2, header.sectionViewportPosition(0) + 4),
    )
    assert view.selected_whole_rows().tolist() == [0]
    assert window.workspace.selection.oids == {int(table.oids[2])}  # C, now first


def test_row_menu_deletes_duplicates_and_freezes(window, qtbot) -> None:
    from chemsuite.app.document import session

    table = window.workspace.table
    view = window.view
    table.add_rows({"ID": [f"X{i}" for i in range(40)], "pIC50": ["1"] * 40}, undoable=False)
    view.selectRow(1)
    titles = menu_titles(window.row_menu(1))
    assert titles[:2] == ["Delete Row", "Duplicate Row"] and "Unfreeze Rows" not in titles
    trigger(window.row_menu(1), "Duplicate Row")
    assert len(table) == 44 and table.cell("ID", int(table.oids[-1])) == "B"
    assert window.workspace.selection.oids == {int(table.oids[-1])}  # the copy is selected
    table.undo()
    assert len(table) == 43

    window.resize(700, 400)
    window.show()
    qtbot.waitExposed(window)
    trigger(window.row_menu(1), "Freeze Rows up to Here")
    rows = view._frozen_rows
    assert view.frozen_row_count == 2 and rows.isVisible()
    assert rows.height() == 2 * view.rowHeight(0) + 1
    view.scrollTo(model_index := view.table_model.index(30, 0))
    assert view.rowViewportPosition(model_index.row()) >= rows.frozen_height()
    view.set_frozen(1)  # columns too: the corner appears
    assert view._frozen_corner is not None and view._frozen_corner.isVisible()

    extras = session.capture(window)
    trigger(window.row_menu(0), "Unfreeze Rows")
    assert view.frozen_row_count == 0 and not rows.isVisible()
    assert not view._frozen_corner.isVisible()
    session.restore(window, extras)
    assert window.view.frozen_row_count == 2
