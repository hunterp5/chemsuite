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
"""Tools ▸ Minimise…: a prepared protein minimised, kept as a new protein
(:mod:`.minimise`)."""

from __future__ import annotations

from PySide6.QtWidgets import QComboBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import ToolDialog

from . import minimise
from .minimise_options import MinimiseFields
from .prepare_dialog import included_of, proteins, shown_protein
from .tool import dialog_parent


class MinimiseDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(dialog_parent(ctx), "Minimise Protein", guide="minimise")
        self.ctx = ctx
        inputs = self.section("Input")
        self.protein = QComboBox()
        for id, name in proteins(ctx):
            self.protein.addItem(name, id)
        chosen = shown_protein(ctx)
        if chosen is not None and self.protein.findData(chosen) >= 0:
            self.protein.setCurrentIndex(self.protein.findData(chosen))
        inputs.row("Protein", self.protein)
        self.fields = MinimiseFields(self.section("Method"), minimise.last_settings())
        self.protein.currentIndexChanged.connect(lambda *_: self.refresh())
        self.fields.forcefield.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def refresh(self) -> None:
        protein = self.protein.currentData()
        s = included_of(self.ctx, protein) if isinstance(protein, str) else None
        if s is not None:
            self.fields.allow_waters(minimise.has_waters(s))
        super().refresh()

    def problem(self) -> str:
        protein = self.protein.currentData()
        if not isinstance(protein, str):
            return "There are no proteins in this document"
        return minimise.problem(self.ctx, protein, self.fields.settings())

    def verb(self) -> str:
        name = self.protein.currentText()
        return f"Minimise {name}" if name else "Minimise"

    def run(self):
        protein = self.protein.currentData()
        if not isinstance(protein, str):
            return None
        return minimise.minimise_protein(self.ctx, protein, self.fields.settings())


def open_minimise(ctx: ToolContext):
    return ctx.show_dialog(MinimiseDialog(ctx))
