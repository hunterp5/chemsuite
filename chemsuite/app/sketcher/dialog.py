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

SMILES, InChI or a molfile pastes straight onto the canvas (Ctrl+V). The sketcher is the only
way to edit a structure in the app; there is no text editor.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.sketcher.canvas import SketchCanvas
from chemsuite.app.sketcher.toolbar import SketchToolBar
from chemsuite.sketch import convert, reaction
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import Sketch


class SketcherDialog(QDialog):
    """Draw a structure. :attr:`molblock` holds the accepted drawing ("" when empty).

    With ``reactions=True`` a drawing with an arrow is accepted as a reaction: :attr:`rxn`
    holds it as an RXN file and :attr:`molblock` stays "". Otherwise OK on a reaction says
    that only a structure can go there.
    """

    def __init__(
        self,
        text: str = "",
        title: str = "Draw Structure",
        parent: QWidget | None = None,
        *,
        ok_text: str = "",
        reactions: bool = False,
        reaction_only: bool = False,
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.reactions = reactions or reaction_only
        #: a reaction cell: OK needs a reaction (or an empty drawing), not a lone structure
        self.reaction_only = reaction_only
        self.rxn = ""
        start = None
        if text.strip() and reaction.looks_like_reaction(text):
            start = reaction.from_text(text)
        start = start or (convert.from_text(text) if text.strip() else None)
        self.canvas = SketchCanvas(History(start) if start is not None else None)
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
        grid.addWidget(self.toolbar.side, 1, 0)
        grid.addWidget(self.canvas, 1, 1)
        grid.addWidget(self.toolbar.bottom, 2, 1)
        grid.setColumnStretch(1, 1)
        grid.setRowStretch(1, 1)
        layout.addLayout(grid, 1)
        bottom = QHBoxLayout()
        bottom.addWidget(self.info, 1)
        bottom.addWidget(self.buttons)
        layout.addLayout(bottom)
        self.resize(760, 560)

        self.toolbar.picked.connect(self.canvas.set_tool)
        self.toolbar.library_requested.connect(self.open_library)
        self.canvas.library_requested.connect(self.open_library)
        self.canvas.tool_changed.connect(lambda _n: self.toolbar.show_tool(self.canvas.tool))
        self.canvas.changed.connect(self._update)
        self.canvas.message.connect(self._show_message)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.molblock = ""
        self._update()
        if text.strip() and start is None:  # SMARTS or a typo: start blank and say why
            self._show_message("This text is not a structure the sketcher can draw.")
        self.canvas.setFocus()

    def clear(self) -> None:
        if not self.canvas.sketch.is_empty():
            self.canvas.set_sketch(Sketch(), "Clear")

    def open_library(self):
        """The template library; the template picked becomes the drawing tool."""
        from chemsuite.app.sketcher.library import TemplateLibrary
        from chemsuite.app.sketcher.tools import TemplateTool
        from chemsuite.sketch import library

        window = TemplateLibrary(self)
        window.chosen.connect(
            lambda name: self.canvas.set_tool(TemplateTool(library.template(name)))
        )
        window.rejected.connect(lambda: self.toolbar.show_tool(self.canvas.tool))
        window.finished.connect(lambda _r: window.deleteLater())
        window.open()
        return window

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

    def accept(self) -> None:
        self.canvas.cancel_preview()
        if not self.canvas.report.valid:
            return
        sketch = self.canvas.sketch
        if sketch.is_reaction or sketch.pluses:
            if not self.reactions:
                self._show_message(
                    "This holds one structure, not a reaction: delete the arrow and plus signs."
                )
                return
            try:
                self.rxn = reaction.to_rxn_block(sketch)
            except (reaction.ReactionError, convert.SketchError) as err:
                self._show_message(str(err))
                return
            self.molblock = ""
            super().accept()
            return
        if self.reaction_only and not sketch.is_empty():
            self._show_message(
                "This holds a reaction: draw an arrow from the reactants to the products."
            )
            return
        self.molblock = "" if sketch.is_empty() else convert.to_molblock(sketch)
        super().accept()
