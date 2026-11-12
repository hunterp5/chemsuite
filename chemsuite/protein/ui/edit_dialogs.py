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

"""Rename a residue or a chain, renumber residues: small dialogs on the selection.

Each asks for one value and makes its edit on OK; an edit that cannot be made (a name another
chain has, a number already taken) is shown in the dialog, which stays open to fix it.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

from PySide6.QtWidgets import QLineEdit, QSpinBox, QWidget

from chemsuite.app.widgets.dialogs.kit import ToolDialog
from chemsuite.protein.structure.edit import EditError

if TYPE_CHECKING:
    from .panel import ProteinPanel


class EditDialog(ToolDialog):
    """One field and the edit it makes (*apply*, which raises EditError)."""

    def __init__(
        self,
        parent: QWidget,
        title: str,
        label: str,
        field: QLineEdit | QSpinBox,
        apply: Callable[[str], None],
        help: str = "",
    ) -> None:
        super().__init__(parent, title, guide="proteins")
        self.field = field
        self._apply = apply
        self.section("Edit").row(label, field, help)
        if isinstance(field, QLineEdit):
            field.textChanged.connect(self.refresh)
        self.refresh()

    def value(self) -> str:
        return self.field.text() if isinstance(self.field, QLineEdit) else str(self.field.value())

    def problem(self) -> str:
        return "" if self.value().strip() else "Type a value"

    def verb(self) -> str:
        return "Apply"

    def accept(self) -> None:
        if self.problem():
            self.refresh()
            return
        try:
            self._apply(self.value())
        except EditError as error:
            self.report(str(error))
            return
        self.done(ToolDialog.DialogCode.Accepted)


def rename_residue(panel: ProteinPanel) -> EditDialog | None:
    found = panel.session.selected_residues()
    if found is None:
        return None
    shown, residues = found
    names = {str(n) for n in shown.structure.res_name[residues]}
    field = QLineEdit(names.pop() if len(names) == 1 else "")
    field.setMaxLength(5)
    title = "Rename Residue" if len(residues) == 1 else f"Rename {len(residues)} Residues"

    def apply(name: str) -> None:
        panel.session.rename_residues(shown.id, residues, name)

    return _show(EditDialog(panel, title, "Name", field, apply,
                            "1 to 5 letters or digits, as in the PDB (HOH, LIG)"))  # fmt: skip


def rename_chain(panel: ProteinPanel) -> EditDialog | None:
    found = panel.session.selected_chain()
    if found is None:
        return None
    shown, chain = found
    field = QLineEdit(str(shown.structure.chain_name[chain]))
    field.setMaxLength(4)

    def apply(name: str) -> None:
        panel.session.rename_chain(shown.id, chain, name)

    return _show(EditDialog(panel, "Rename Chain", "Name", field, apply,
                            "1 to 4 letters or digits; no other chain may have it"))  # fmt: skip


def renumber(panel: ProteinPanel) -> EditDialog | None:
    found = panel.session.selected_residues()
    if found is None:
        return None
    shown, residues = found
    field = QSpinBox()
    field.setRange(-999_999, 9_999_999)
    field.setValue(int(shown.structure.res_seq[residues[0]]))

    def apply(first: str) -> None:
        panel.session.renumber(shown.id, residues, int(first))

    title = "Renumber Residue" if len(residues) == 1 else f"Renumber {len(residues)} Residues"
    help = "The selected residues are numbered from here up, in order"
    return _show(EditDialog(panel, title, "First number", field, apply, help))


def _show(dialog: EditDialog) -> EditDialog:
    dialog.show()
    return dialog
