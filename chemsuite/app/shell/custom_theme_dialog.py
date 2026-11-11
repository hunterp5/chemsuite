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
"""Create, edit and delete named custom colour themes, with a live preview."""

from __future__ import annotations

from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from . import theme

NEW_THEME = "(New theme)"
STARTS = {
    "Light": theme.light_palette,
    "Dark": theme.dark_palette,
    "Groovy (random)": theme.groovy_palette,
}


class CustomThemeDialog(QDialog):
    """Edits a copy of *themes* (name → colours).

    After Accept, ``themes`` holds the edited set and ``chosen`` the theme to switch to (None:
    keep the current theme, unless it was deleted).
    """

    def __init__(
        self, themes: dict[str, dict[str, str]], current: str | None, parent: QWidget | None = None
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Custom Themes")
        self.themes = {name: dict(colours) for name, colours in themes.items()}
        self.chosen: str | None = None
        self.colours = theme.default_custom_colours()
        self.swatches: dict[str, QPushButton] = {}

        self.combo = QComboBox()
        self.combo.currentIndexChanged.connect(self._load_selected)
        self.delete_button = QPushButton("Delete")
        self.delete_button.clicked.connect(self.delete_current)
        start = QComboBox()
        start.addItem("Start from…")
        start.addItems(list(STARTS))
        start.activated.connect(lambda i: self._start_from(start, i))
        top = QHBoxLayout()
        top.addWidget(QLabel("Theme:"))
        top.addWidget(self.combo, 1)
        top.addWidget(self.delete_button)
        top.addWidget(start)

        grid = QGridLayout()
        for i, (key, label, _role) in enumerate(theme.CUSTOM_ROLES):
            swatch = QPushButton()
            swatch.setFixedSize(84, 24)
            swatch.setToolTip(f"Choose the colour for {label.lower()}")
            swatch.clicked.connect(lambda _c=False, k=key: self.pick(k))
            self.swatches[key] = swatch
            row, col = i % 8, (i // 8) * 2
            grid.addWidget(QLabel(label + ":"), row, col)
            grid.addWidget(swatch, row, col + 1)
        colours_box = QGroupBox("Colours")
        colours_box.setLayout(grid)

        self.preview = _Preview()
        middle = QHBoxLayout()
        middle.addWidget(colours_box)
        middle.addWidget(self.preview, 1)

        self.name = QLineEdit()
        self.name.setPlaceholderText("Theme name")
        self.name.textChanged.connect(self._update_buttons)
        buttons = QDialogButtonBox()
        self.save_button = buttons.addButton("Save and Use", QDialogButtonBox.ButtonRole.AcceptRole)
        self.save_button.setDefault(True)
        close = buttons.addButton(QDialogButtonBox.StandardButton.Close)
        self.save_button.clicked.connect(self.save)
        close.clicked.connect(self.close_editor)
        bottom = QHBoxLayout()
        bottom.addWidget(QLabel("Name:"))
        bottom.addWidget(self.name, 1)
        bottom.addWidget(buttons)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addLayout(middle, 1)
        layout.addLayout(bottom)
        self._original = {name: dict(c) for name, c in self.themes.items()}
        self._fill_combo(current if current in self.themes else None)

    # -- theme list --------------------------------------------------------------------------------
    def _fill_combo(self, select: str | None) -> None:
        self.combo.blockSignals(True)
        self.combo.clear()
        self.combo.addItem(NEW_THEME, "")
        for name in sorted(self.themes, key=str.lower):
            self.combo.addItem(name, name)
        self.combo.setCurrentIndex(max(0, self.combo.findData(select or "")))
        self.combo.blockSignals(False)
        self._load_selected()

    def _load_selected(self, *_args) -> None:
        name = self.combo.currentData() or ""
        self.colours = theme.default_custom_colours()  # fills any colour a saved theme lacks
        if name:
            self.colours.update(theme.palette_colours(theme.custom_palette(self.themes[name])))
        self.name.setText(name)
        self.delete_button.setEnabled(bool(name))
        self._refresh()

    def _start_from(self, combo: QComboBox, index: int) -> None:
        make = STARTS.get(combo.itemText(index))
        combo.setCurrentIndex(0)
        if make is not None:
            self.colours = theme.palette_colours(make())
            self._refresh()

    # -- colours -----------------------------------------------------------------------------------
    def pick(self, key: str) -> None:
        chosen = QColorDialog.getColor(QColor(self.colours[key]), self, "Choose Colour")
        if chosen.isValid():
            self.set_colour(key, chosen.name())

    def set_colour(self, key: str, colour: str) -> None:
        self.colours[key] = QColor(colour).name()
        self._refresh()

    def _refresh(self) -> None:
        for key, swatch in self.swatches.items():
            colour = QColor(self.colours[key])
            text = theme.contrasting_text(colour).name()
            swatch.setText(colour.name())
            swatch.setStyleSheet(
                f"background-color: {colour.name()}; color: {text};"
                " border: 1px solid #777; border-radius: 3px;"
            )
        self.preview.show_palette(theme.custom_palette(self.colours))
        self._update_buttons()

    def _update_buttons(self, *_args) -> None:
        self.save_button.setEnabled(bool(_clean(self.name.text())))

    # -- saving ------------------------------------------------------------------------------------
    def save(self) -> None:
        name = _clean(self.name.text())
        if not name:
            return
        old = self.combo.currentData() or ""
        if old and old != name:  # renamed
            self.themes.pop(old, None)
        self.themes[name] = dict(self.colours)
        self.chosen = name
        self.accept()

    def confirm_delete(self, name: str) -> bool:
        answer = QMessageBox.question(self, "Delete Theme", f"Delete the custom theme “{name}”?")
        return answer == QMessageBox.StandardButton.Yes

    def delete_current(self) -> None:
        name = self.combo.currentData() or ""
        if name and self.confirm_delete(name):
            self.themes.pop(name, None)
            self._fill_combo(None)

    def close_editor(self) -> None:
        """Close; deletions still count, so accept when the saved set changed."""
        if self.themes != self._original:
            self.accept()
        else:
            self.reject()


def _clean(name: str) -> str:
    return " ".join(name.split())


class _Preview(QGroupBox):
    """A few ordinary widgets drawn in the edited palette."""

    def __init__(self) -> None:
        super().__init__("Preview")
        self.setAutoFillBackground(True)
        table = QTableWidget(4, 2)
        table.setHorizontalHeaderLabels(["ID", "pIC50"])
        table.setAlternatingRowColors(True)
        table.horizontalHeader().setStretchLastSection(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        for r, (cid, value) in enumerate((("CPD-1", "6.2"), ("CPD-2", "7.9"),
                                          ("CPD-3", "5.4"), ("CPD-4", "8.1"))):  # fmt: skip
            table.setItem(r, 0, QTableWidgetItem(cid))
            table.setItem(r, 1, QTableWidgetItem(value))
        table.selectRow(1)
        field = QLineEdit("Text field")
        link = QLabel('<a href="#">A link</a>')
        disabled = QPushButton("Disabled")
        disabled.setEnabled(False)
        buttons = QHBoxLayout()
        buttons.addWidget(QPushButton("Button"))
        buttons.addWidget(disabled)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Window text"))
        layout.addWidget(field)
        layout.addWidget(QCheckBox("Check box"))
        layout.addLayout(buttons)
        layout.addWidget(link)
        layout.addWidget(table, 1)
        self.setToolTip("Tooltips use the tooltip colours.")

    def show_palette(self, palette) -> None:
        self.setPalette(palette)
        for child in self.findChildren(QWidget):
            child.setPalette(palette)
