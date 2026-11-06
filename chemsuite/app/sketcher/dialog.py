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

"""The sketcher window: tool strips around the canvas; OK gives back a molfile.

Tools run along the top, rings along the bottom and elements down the left, as in Ketcher.
Undo, redo, cut, copy, paste and clear are on the keyboard and the canvas's right-click menu,
so the window holds only the drawing. The line under the canvas appears only to explain a
problem (an atom with a bad valence, a label that could not be read).

A text line above the drawing takes SMILES, InChI or a molfile: Enter draws it. It is the way
in for anyone who would rather type, and what typing into a structure cell opens with. Text
that was typed and then left undrawn-over is given back exactly as typed (``result_text``).
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.sketcher.canvas import SketchCanvas
from chemsuite.app.sketcher.toolbar import SketchToolBar
from chemsuite.sketch import convert
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import Sketch


class SketcherDialog(QDialog):
    """Draw a structure. :attr:`molblock` holds the accepted drawing ("" when empty)."""

    def __init__(
        self,
        text: str = "",
        title: str = "Draw Structure",
        parent: QWidget | None = None,
        *,
        ok_text: str = "",
        typed: str = "",
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        start = convert.from_text(text) if text.strip() else None
        self.canvas = SketchCanvas(History(start) if start is not None else None)
        self.entry = QLineEdit()
        self.entry.setObjectName("sketch-text-entry")
        self.entry.setPlaceholderText("Type or paste SMILES, InChI or a molfile, then Enter")
        self.entry.setClearButtonEnabled(True)
        self.entry.returnPressed.connect(self.load_text)
        #: (text, the edit that drew it) for the text last drawn from the entry line
        self._typed: tuple[str, object] | None = None
        self.toolbar = SketchToolBar(self)  # holds the strips laid out below; not shown itself
        self.toolbar.hide()
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setStyleSheet("color: #c62828")
        self.info.hide()
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        if ok_text:
            self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText(ok_text)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        grid = QGridLayout()  # rings line up with the canvas edge, elements to its left
        grid.setSpacing(6)
        grid.addWidget(self.toolbar.top, 0, 0, 1, 2)
        grid.addWidget(self.entry, 1, 1)
        grid.addWidget(self.toolbar.side, 2, 0)
        grid.addWidget(self.canvas, 2, 1)
        grid.addWidget(self.toolbar.bottom, 3, 1)
        grid.setColumnStretch(1, 1)
        grid.setRowStretch(2, 1)
        layout.addLayout(grid, 1)
        bottom = QHBoxLayout()
        bottom.addWidget(self.info, 1)
        bottom.addWidget(self.buttons)
        layout.addLayout(bottom)
        self.resize(760, 560)

        self.toolbar.picked.connect(self.canvas.set_tool)
        self.canvas.tool_changed.connect(lambda _n: self.toolbar.show_tool(self.canvas.tool))
        self.canvas.changed.connect(self._update)
        self.canvas.message.connect(self._show_message)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.molblock = ""
        self._update()
        if typed or (text.strip() and start is None):
            # typed into a cell, or text that cannot be drawn (SMARTS, a typo): keep it to edit
            self.entry.setText(typed or text.strip())
            if start is None and not typed:
                self._show_message("This text is not a structure the sketcher can draw.")
            self.entry.setFocus()
            self.entry.end(False)
        else:
            self.canvas.setFocus()

    TYPED_EDIT = "Type structure"

    def load_text(self) -> bool:
        """Draw the entry line's text, replacing the drawing (one undo step)."""
        text = self.entry.text().strip()
        if not text:
            return False
        sketch = convert.from_text(text)
        if sketch is None:
            self._show_message(f"Not SMILES, InChI or a molfile the sketcher can read: {text}")
            return False
        self.canvas.set_sketch(sketch, self.TYPED_EDIT)
        self._typed = (text, self.canvas.history.last_edit)
        self.canvas.setFocus()
        return True

    def result_text(self) -> str:
        """What the dialog gives back: the typed text exactly as typed when nothing was drawn
        after it, otherwise the drawing as a molfile ("" when empty)."""
        if self._typed is not None and self.canvas.history.last_edit is self._typed[1]:
            return self._typed[0]
        return self.molblock

    def clear(self) -> None:
        if not self.canvas.sketch.is_empty():
            self.canvas.set_sketch(Sketch(), "Clear")

    def _update(self) -> None:
        report = self.canvas.report
        ok = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        ok.setEnabled(report.valid)
        wedges = len(self.canvas.wedge_problems)
        if report.valid and wedges:
            self._show_message(
                f"{wedges} wedge{'s' if wedges != 1 else ''} that may not mean what you intend "
                "(marked; hover for why)"
            )
        elif self.canvas.sketch.is_empty() or report.valid:
            self.info.hide()
            self.info.clear()
        else:
            n = len(report.problems)
            self._show_message(
                f"{n} atom{'s' if n != 1 else ''} with a valence problem (marked; hover for why)"
            )

    def _show_message(self, text: str) -> None:
        self.info.setText(text)
        self.info.setVisible(bool(text))

    def keyPressEvent(self, event) -> None:
        # Enter in the text line draws the text (returnPressed); it must not also press OK.
        if self.entry.hasFocus() and event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            event.accept()
            return
        super().keyPressEvent(event)

    def accept(self) -> None:
        self.canvas.cancel_preview()
        pending = self.entry.text().strip()
        undrawn = pending and (self._typed is None or pending != self._typed[0])
        if undrawn and not self.load_text():  # OK with undrawn text: draw it first
            return
        if not self.canvas.report.valid:
            return
        sketch = self.canvas.sketch
        self.molblock = "" if sketch.is_empty() else convert.to_molblock(sketch)
        super().accept()
