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

"""Align Drawings…: lay out a structure column's drawings to match a common core.

Display only: the structures and their saved coordinates are unchanged, the setting is saved
with the document's layout, and Stop Aligning Drawings brings the stored layouts back.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtWidgets import QWidget

from chemsuite.core.depict import alignment_core

from .dialogs import FIELD_WIDTH, ToolDialog
from .structure_search import QueryInput


class AlignDrawingsDialog(ToolDialog):
    """Ask for the core (typed, drawn or from the current row) and apply it on OK."""

    def __init__(
        self,
        parent: QWidget | None,
        column: str,
        apply: Callable[[str, str], bool],
        current_structure: Callable[[], str] | None = None,
        start: str = "",
    ) -> None:
        super().__init__(parent, f"Align Drawings — {column}", guide="columns")
        self.column = column
        self._apply = apply
        section = self.section("Input")
        self.core = QueryInput(current_structure, compact=True)
        self.core.set_text(start)
        self.core.setFixedWidth(FIELD_WIDTH)
        section.row(
            "Core",
            self.core,
            help="Rows containing this core are drawn with it in the same place and "
            "orientation; other rows keep their usual drawing.",
        )
        self.core.changed.connect(self.refresh)
        self.refresh()

    def problem(self) -> str:
        text = self.core.text().strip()
        if not text:
            return "Enter, draw or take a core from the current row."
        if alignment_core(text) is None:
            return "The core is not SMILES or SMARTS."
        return ""

    def verb(self) -> str:
        return "Align"

    def accept(self) -> None:
        if self._problem():
            self.refresh()
            return
        if self._apply(self.column, self.core.text().strip()):
            super(ToolDialog, self).accept()
