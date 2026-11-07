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

from PySide6.QtWidgets import QDialog, QDialogButtonBox, QFormLayout, QLabel, QVBoxLayout, QWidget

from chemsuite.core.depict import alignment_core

from .structure_search import QueryInput


class AlignDrawingsDialog(QDialog):
    """Ask for the core (typed, drawn or from the current row) and apply it on OK."""

    def __init__(
        self,
        parent: QWidget | None,
        column: str,
        apply: Callable[[str, str], bool],
        current_structure: Callable[[], str] | None = None,
        start: str = "",
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"Align Drawings: {column}")
        self.column = column
        self._apply = apply
        self.core = QueryInput(current_structure)
        self.core.set_text(start)
        self.core.setToolTip(
            "Rows containing this core are drawn with it in the same place and orientation; "
            "other rows keep their usual drawing."
        )
        self.problem = QLabel()
        self.problem.hide()
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        form = QFormLayout()
        form.addRow("Core:", self.core)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.problem)
        layout.addWidget(self.buttons)
        self.core.changed.connect(self._check)
        self._check()
        self.resize(460, self.sizeHint().height())

    def _check(self) -> None:
        text = self.core.text().strip()
        problem = "" if not text or alignment_core(text) is not None else "⚠ Not SMILES or SMARTS"
        self.problem.setText(problem)
        self.problem.setVisible(bool(problem))
        ok = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        ok.setEnabled(bool(text) and not problem)

    def accept(self) -> None:
        if self._apply(self.column, self.core.text().strip()):
            super().accept()
