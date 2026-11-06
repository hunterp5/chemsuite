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

"""The sketcher window: tool bar, canvas and a status line; OK gives back a molfile."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from rdkit.Chem import Descriptors, rdMolDescriptors

from chemsuite.app.sketcher.canvas import SketchCanvas
from chemsuite.app.sketcher.toolbar import SketchToolBar
from chemsuite.sketch import convert
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import Sketch


class SketcherDialog(QDialog):
    """Draw a structure. :attr:`molblock` holds the accepted drawing ("" when empty)."""

    def __init__(
        self, text: str = "", title: str = "Draw Structure", parent: QWidget | None = None
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        start = convert.from_text(text) if text.strip() else None
        self.canvas = SketchCanvas(History(start) if start is not None else None)
        self.toolbar = SketchToolBar(self)  # holds the strips laid out below; not shown itself
        self.toolbar.hide()
        self.info = QLabel(" ")
        self.undo_button = QPushButton("Undo")
        self.redo_button = QPushButton("Redo")
        self.paste_button = QPushButton("Paste")
        self.copy_button = QPushButton("Copy")
        self.clear_button = QPushButton("Clear")
        self.paste_button.setToolTip("Add the structure on the clipboard (molfile or SMILES)")
        self.copy_button.setToolTip("Copy as molfile and SMILES (the selection, if any)")
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        commands = QHBoxLayout()
        for button in (
            self.undo_button,
            self.redo_button,
            self.paste_button,
            self.copy_button,
            self.clear_button,
        ):
            button.setAutoDefault(False)
            commands.addWidget(button)
        commands.addStretch(1)

        layout = QVBoxLayout(self)
        layout.addWidget(self.toolbar.top)
        middle = QHBoxLayout()
        middle.addWidget(self.canvas, 1)
        middle.addWidget(self.toolbar.side)
        layout.addLayout(middle, 1)
        layout.addWidget(self.toolbar.bottom)
        layout.addLayout(commands)
        layout.addWidget(self.info)
        layout.addWidget(self.buttons)
        self.resize(760, 560)

        self.toolbar.picked.connect(self.canvas.set_tool)
        self.canvas.tool_changed.connect(lambda _n: self.toolbar.show_tool(self.canvas.tool))
        self.canvas.changed.connect(self._update)
        self.canvas.message.connect(self._show_message)
        self.undo_button.clicked.connect(self.canvas.history.undo)
        self.redo_button.clicked.connect(self.canvas.history.redo)
        self.paste_button.clicked.connect(self.canvas.paste)
        self.copy_button.clicked.connect(self.canvas.copy)
        self.clear_button.clicked.connect(self._clear)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.molblock = ""
        self._update()
        self.canvas.setFocus()

    def _clear(self) -> None:
        if not self.canvas.sketch.is_empty():
            self.canvas.set_sketch(Sketch(), "Clear")

    def _update(self) -> None:
        history = self.canvas.history
        self.undo_button.setEnabled(history.can_undo)
        self.redo_button.setEnabled(history.can_redo)
        self.undo_button.setToolTip(f"Undo {history.undo_name}" if history.can_undo else "")
        self.redo_button.setToolTip(f"Redo {history.redo_name}" if history.can_redo else "")
        report = self.canvas.report
        ok = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        ok.setEnabled(report.valid)
        self.copy_button.setEnabled(not self.canvas.sketch.is_empty())
        if self.canvas.sketch.is_empty():
            self.info.setText("Empty")
            self.info.setStyleSheet("color: palette(placeholder-text)")
        elif not report.valid:
            n = len(report.problems)
            self.info.setText(
                f"{n} atom{'s' if n != 1 else ''} with a valence problem (marked; hover for why)"
            )
            self.info.setStyleSheet("color: #c62828")
        else:
            mol = convert.to_mol(self.canvas.sketch)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            self.info.setText(f"{formula} · MW {Descriptors.MolWt(mol):.2f}")
            self.info.setStyleSheet("color: palette(placeholder-text)")

    def _show_message(self, text: str) -> None:
        self.info.setText(text)
        self.info.setStyleSheet("color: #c62828")

    def accept(self) -> None:
        self.canvas.cancel_preview()
        if not self.canvas.report.valid:
            return
        sketch = self.canvas.sketch
        self.molblock = "" if sketch.is_empty() else convert.to_molblock(sketch)
        super().accept()
