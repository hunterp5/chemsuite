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

"""Tools ▸ Table ▸ Split Column and Join Columns, with a live preview of the result."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.selection import Selection
from chemsuite.core.table import Table
from chemsuite.core.table.split_join import (
    JOIN_DELIMITERS,
    MAX_SPLIT_COLUMNS,
    SPLIT_DELIMITERS,
    join_columns,
    plan_join,
    plan_split,
    split_column,
    split_names,
    unescape,
)

PREVIEW_ROWS = 8
SPLIT_TITLES = {
    "auto": "Auto (detect)",
    "comma": "Comma (,)",
    "semicolon": "Semicolon (;)",
    "tab": "Tab",
    "pipe": "Pipe (|)",
    "space": "Spaces",
    "custom": "Custom",
}
JOIN_TITLES = {
    "comma_space": "Comma and space (, )",
    "comma": "Comma (,)",
    "semicolon": "Semicolon (; )",
    "tab": "Tab",
    "pipe": "Pipe (|)",
    "space": "Space",
    "dot": "Dot (.), joins SMILES",
    "none": "Nothing (run together)",
    "custom": "Custom",
}
KEEP_TITLES = {
    "all": "All fields, one column each",
    "first": "First field only",
    "last": "Last field only",
    "largest": "Largest value only",
    "smallest": "Smallest value only",
}


def _describe(delimiter: str) -> str:
    return {"\t": "tab", " ": "spaces"}.get(delimiter, f"“{delimiter}”")


class _Base(QDialog):
    """Shared pieces: delimiter choice with a custom field, scope, preview, debounce."""

    def __init__(self, table: Table, selection: Selection, title: str, parent: QWidget | None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.table = table
        self.selection = selection
        self.delimiter = QComboBox()
        self.custom = QLineEdit()
        self.custom.setPlaceholderText(r"e.g. :  or  \t")
        self.custom.setMaximumWidth(120)
        self.selected_only = QCheckBox(f"Selected rows only ({len(selection):,})")
        self.selected_only.setEnabled(len(selection) > 0)
        self.preview = QTableWidget()
        self.preview.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.preview.verticalHeader().hide()
        self.preview.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.preview.setMinimumHeight(180)
        self.note = QLabel(" ")
        self.note.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(120)
        self._timer.timeout.connect(self.update_preview)
        self.delimiter.currentIndexChanged.connect(self._delimiter_changed)
        self.custom.textChanged.connect(self.schedule)
        self.selected_only.toggled.connect(self.schedule)

    def _delimiter_row(self) -> QHBoxLayout:
        row = QHBoxLayout()
        row.addWidget(self.delimiter, 1)
        row.addWidget(self.custom)
        return row

    def _delimiter_changed(self) -> None:
        self.custom.setVisible(self.delimiter.currentData() == "custom")
        self.schedule()

    def schedule(self, *_args) -> None:
        self._timer.start()

    def oids(self):
        """Rows to work on: the selection (in storage order) or every row."""
        return self.selection.as_array() if self.selected_only.isChecked() else None

    def _sample(self) -> np.ndarray:
        """Oids of the first rows shown (and in scope), for the preview."""
        oids = self.table.oids[self.table.view()]
        if self.selected_only.isChecked():
            oids = oids[np.isin(oids, self.selection.as_array())]
        return oids[:PREVIEW_ROWS]

    def _say(self, text: str, error: bool = False) -> None:
        self.note.setText(text or " ")
        self.note.setStyleSheet("color: #c62828" if error else "")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not error)

    def _fill(self, headers: list[str], rows: list[list[str]]) -> None:
        self.preview.clear()
        self.preview.setColumnCount(len(headers))
        self.preview.setRowCount(len(rows))
        self.preview.setHorizontalHeaderLabels(headers)
        for r, row in enumerate(rows):
            for c, text in enumerate(row):
                self.preview.setItem(r, c, QTableWidgetItem(text))
        self.preview.resizeColumnsToContents()
        for c in range(len(headers)):  # one long cell must not push the results out of view
            self.preview.setColumnWidth(c, min(self.preview.columnWidth(c), 200))

    def update_preview(self) -> None:
        raise NotImplementedError

    def accept(self) -> None:
        self.update_preview()
        if not self.buttons.button(QDialogButtonBox.StandardButton.Ok).isEnabled():
            return
        try:
            self.apply()
        except ValueError as err:
            QMessageBox.warning(self, self.windowTitle(), str(err))
            return
        super().accept()

    def apply(self) -> None:
        raise NotImplementedError


class SplitColumnDialog(_Base):
    def __init__(
        self,
        table: Table,
        selection: Selection,
        column: str | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(table, selection, "Split Column", parent)
        self.result_names: list[str] = []
        self.column = QComboBox()
        self.column.addItems(table.column_names)
        if column in table:
            self.column.setCurrentText(column)
        for key, title in SPLIT_TITLES.items():
            self.delimiter.addItem(title, key)
        self.keep = QComboBox()
        for key, title in KEEP_TITLES.items():
            self.keep.addItem(title, key)
        self.prefix = QLineEdit()

        form = QFormLayout()
        form.addRow("Column:", self.column)
        form.addRow("Split on:", self._delimiter_row())
        form.addRow("Keep:", self.keep)
        form.addRow("New column names:", self.prefix)
        form.addRow("", self.selected_only)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(QLabel("Preview (first shown rows):"))
        layout.addWidget(self.preview, 1)
        layout.addWidget(self.note)
        layout.addWidget(self.buttons)
        self.resize(620, 480)

        self.column.currentIndexChanged.connect(self._column_changed)
        self.keep.currentIndexChanged.connect(self.schedule)
        self.prefix.textChanged.connect(self.schedule)
        self._delimiter_changed()
        self._column_changed()
        self.update_preview()

    def _column_changed(self, *_args) -> None:
        self.prefix.setPlaceholderText(self.column.currentText())
        self.schedule()

    def chosen_delimiter(self) -> str | None:
        """None for auto; raises ValueError for an empty custom delimiter."""
        key = self.delimiter.currentData()
        if key == "auto":
            return None
        if key == "custom":
            text = unescape(self.custom.text())
            if not text:
                raise ValueError("Type the custom delimiter.")
            return text
        return SPLIT_DELIMITERS[key]

    def update_preview(self) -> None:
        self._timer.stop()
        name = self.column.currentText()
        if not name:
            self._fill([], [])
            self._say("The table has no columns.", error=True)
            return
        try:
            delimiter = self.chosen_delimiter()
            # Detect on every row in scope, so the preview matches what OK will do.
            full = plan_split(
                self.table, name, delimiter, keep=self.keep.currentData(), oids=self.oids()
            )
            sample = plan_split(
                self.table, name, full.delimiter, keep=self.keep.currentData(), oids=self._sample()
            )
        except ValueError as err:
            self._fill([], [])
            self._say(str(err), error=True)
            return
        prefix = self.prefix.text().strip() or name
        keep = self.keep.currentData()
        width = full.width
        headers = split_names(prefix, width) if keep == "all" else [f"{prefix}_{keep}"]
        rows = [
            [self.table.cell(name, int(oid))]
            + [f[i] if i < len(f) else "" for i in range(len(headers))]
            for oid, f in zip(sample.oids.tolist(), sample.fields, strict=True)
        ]
        self._fill([name, *headers], rows)
        if not width:
            self._say("Nothing to split: the cells are blank.", error=True)
        elif width > MAX_SPLIT_COLUMNS:
            self._say(
                f"Splitting on {_describe(full.delimiter)} would make {width:,} columns; "
                "check the delimiter.",
                error=True,
            )
        else:
            self._say("")

    def apply(self) -> None:
        self.result_names = split_column(
            self.table,
            self.column.currentText(),
            self.chosen_delimiter(),
            prefix=self.prefix.text(),
            keep=self.keep.currentData(),
            oids=self.oids(),
        )


class JoinColumnsDialog(_Base):
    def __init__(
        self,
        table: Table,
        selection: Selection,
        columns: list[str] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(table, selection, "Join Columns", parent)
        self.result_name = ""
        chosen = list(columns or [])
        self.columns = QListWidget()
        self.columns.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.columns.setToolTip("Tick the columns to join; drag to change their order.")
        order = [n for n in chosen if n in table] + [
            n for n in table.column_names if n not in chosen
        ]
        for name in order:
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if name in chosen else Qt.CheckState.Unchecked)
            self.columns.addItem(item)
        self.columns.setMaximumHeight(150)
        for key, title in JOIN_TITLES.items():
            self.delimiter.addItem(title, key)
        self.name = QLineEdit()
        self.skip_blank = QCheckBox("Skip blank cells (no doubled delimiters)")
        self.skip_blank.setChecked(True)

        form = QFormLayout()
        form.addRow("Columns:", self.columns)
        form.addRow("Join with:", self._delimiter_row())
        form.addRow("New column name:", self.name)
        form.addRow("", self.skip_blank)
        form.addRow("", self.selected_only)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(QLabel("Preview (first shown rows):"))
        layout.addWidget(self.preview, 1)
        layout.addWidget(self.note)
        layout.addWidget(self.buttons)
        self.resize(560, 560)

        self.columns.itemChanged.connect(self.schedule)
        self.columns.model().rowsMoved.connect(self.schedule)
        self.name.textChanged.connect(self.schedule)
        self.skip_blank.toggled.connect(self.schedule)
        self._delimiter_changed()
        self.update_preview()

    def chosen_columns(self) -> list[str]:
        items = (self.columns.item(i) for i in range(self.columns.count()))
        return [i.text() for i in items if i.checkState() == Qt.CheckState.Checked]

    def set_checked(self, names: list[str]) -> None:
        for i in range(self.columns.count()):
            item = self.columns.item(i)
            item.setCheckState(
                Qt.CheckState.Checked if item.text() in names else Qt.CheckState.Unchecked
            )
        self.update_preview()

    def chosen_delimiter(self) -> str:
        key = self.delimiter.currentData()
        return unescape(self.custom.text()) if key == "custom" else JOIN_DELIMITERS[key]

    def update_preview(self) -> None:
        self._timer.stop()
        names = self.chosen_columns()
        self.name.setPlaceholderText("_".join(names) if names else "New column")
        if len(names) < 2:
            self._fill(names, [])
            self._say("Tick at least two columns to join.", error=True)
            return
        oids, joined = plan_join(
            self.table,
            names,
            self.chosen_delimiter(),
            skip_blank=self.skip_blank.isChecked(),
            oids=self._sample(),
        )
        result = self.name.text().strip() or "_".join(names)
        rows = [
            [self.table.cell(n, int(oid)) for n in names] + [text]
            for oid, text in zip(oids.tolist(), joined, strict=True)
        ]
        self._fill([*names, result], rows)
        self._say("")

    def apply(self) -> None:
        self.result_name = join_columns(
            self.table,
            self.chosen_columns(),
            self.chosen_delimiter(),
            name=self.name.text(),
            skip_blank=self.skip_blank.isChecked(),
            oids=self.oids(),
        )
