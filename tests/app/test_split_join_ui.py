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
def assays(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["A", "B", "C"],
            "Assay": ["IC50; 5.2", "Ki; 7", "EC50; 6.1"],
            "Series": ["S1", "", "S3"],
        },
        undoable=False,
    )
    window._install(Workspace(table))
    return window


def preview_headers(dialog) -> list[str]:
    table = dialog.preview
    return [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]


def test_split_dialog_previews_and_splits(assays) -> None:
    window = assays
    dialog = window.split_column("Assay")
    assert preview_headers(dialog) == ["Assay", "Assay_1", "Assay_2"]
    assert dialog.preview.item(0, 2).text() == "5.2"
    assert not dialog.status.text().strip()
    dialog.prefix.setText("Data")
    dialog.update_preview()
    assert preview_headers(dialog)[1:] == ["Data_1", "Data_2"]
    dialog.accept()
    table = window.workspace.table
    assert table.column_names == ["ID", "Assay", "Data_1", "Data_2", "Series"]
    assert table.kind("Data_2") == "number"


def test_split_dialog_custom_delimiter_and_errors(assays) -> None:
    dialog = assays.split_column("Assay")
    dialog.delimiter.setCurrentIndex(dialog.delimiter.findData("custom"))
    dialog.update_preview()
    assert "custom delimiter" in dialog.status.text()
    ok = dialog.primary
    assert not ok.isEnabled()
    dialog.custom.setText("C")
    dialog.update_preview()
    assert ok.isEnabled() and not dialog.status.text().strip()
    dialog.keep.setCurrentIndex(dialog.keep.findData("last"))
    dialog.update_preview()
    assert preview_headers(dialog) == ["Assay", "Assay_last"]
    dialog.reject()


def test_split_selected_rows_only(assays) -> None:
    window = assays
    table = window.workspace.table
    window.workspace.selection.set([int(table.oids[1])], "test")
    dialog = window.split_column("Assay")
    dialog.rows.set_value("selected")
    dialog.update_preview()
    assert dialog.preview.rowCount() == 1
    dialog.accept()
    assert [table.cell("Assay_1", o) for o in table.oids] == ["", "Ki", ""]


def test_join_dialog_uses_selected_columns_and_order(assays, qtbot) -> None:
    window = assays
    table = window.workspace.table
    dialog = window.join_columns(["Series", "ID"])
    assert dialog.chosen_columns() == ["Series", "ID"]  # given order first
    dialog.delimiter.setCurrentIndex(dialog.delimiter.findData("none"))
    dialog.update_preview()
    assert preview_headers(dialog) == ["Series", "ID", "Series_ID"]
    assert dialog.preview.item(0, 2).text() == "S1A"
    assert dialog.preview.item(1, 2).text() == "B"
    dialog.name.setText("Key")
    dialog.accept()
    assert table.column_names == ["ID", "Assay", "Series", "Key"]
    assert table.history.undo_label == "Join 2 columns"


def test_join_needs_two_ticked_columns(assays) -> None:
    dialog = assays.join_columns()
    dialog.set_checked(["ID"])
    assert "at least two" in dialog.status.text()
    dialog.accept()  # refused
    assert dialog.isVisible()
    dialog.reject()
