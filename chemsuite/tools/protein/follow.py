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

"""The workbench follows the table's current row: its own pose, or the ligand it names.

A row's pose is its ``Ligand`` when that is 3D and is not just the protein's residue (a ligand
moved in the sketcher, or docked); otherwise the residue its ``Ligand ID`` names is focused.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from chemsuite.app.document.tools import ToolContext
from chemsuite.protein.chem import ligands, poses

from .rows import LIGAND_COLUMN, LIGAND_ID_COLUMN, PROTEIN_COLUMN, ROW_COLUMNS, mol_of

if TYPE_CHECKING:
    from chemsuite.protein.ui.session import ProteinSession


class TableRow:
    """The document's current row, for the workbench (a :class:`~chemsuite.protein.ui.host.
    LinkedRow`)."""

    columns = ROW_COLUMNS

    def __init__(self, ctx: ToolContext) -> None:
        self._ctx = ctx

    def follow(self, session: ProteinSession, oid: int | None) -> bool:
        table = self._ctx.table
        if oid is None or PROTEIN_COLUMN not in table or oid not in table.oids:
            return False
        protein = str(table.cell(PROTEIN_COLUMN, oid) or "")
        shown = session.display(protein)
        if shown is None:
            _leave_pose(session)
            return False
        label = str(table.cell(LIGAND_ID_COLUMN, oid) or "") if LIGAND_ID_COLUMN in table else ""
        residue = ligands.residue_named(shown.structure, label) if label else -1
        pose = self._pose(oid)
        if pose is not None and (residue < 0 or not poses.same_pose(
                shown.structure, residue, pose, shown.active, shown.model)):  # fmt: skip
            session.set_visible(protein, True)
            return session.show_pose(protein, pose, oid)
        if residue < 0:
            _leave_pose(session)
            return False
        session.set_visible(protein, True)
        if session.pose is not None or session.focused_residue() != (protein, residue):
            session.focus(protein, residue)
        return True

    def missing(self, session: ProteinSession, oid: int | None) -> str:
        table = self._ctx.table
        if oid is None or PROTEIN_COLUMN not in table or oid not in table.oids:
            return ""
        protein = str(table.cell(PROTEIN_COLUMN, oid) or "")
        return protein if protein and session.display(protein) is None else ""

    def _pose(self, oid: int):
        table = self._ctx.table
        if LIGAND_COLUMN not in table or table.kind(LIGAND_COLUMN) != "structure":
            return None
        (given,) = self._ctx.structure_inputs(LIGAND_COLUMN, [oid])
        return poses.pose_of(mol_of(given))


def _leave_pose(session: ProteinSession) -> None:
    """A row that shows nothing here takes away another row's pose."""
    if session.pose is not None:
        session.clear_focus()
