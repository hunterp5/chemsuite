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
"""Protein ▸ Dock… and Rescore Poses…: AutoDock Vina on the document's ligands.

The box comes from what the workbench is looking at — the focused ligand, the selection, or a
ligand row — or from numbers typed in. The receptor must be prepared first (Vina types atoms
by their hydrogens), which :meth:`DockDialog.problem` says when it is not.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QComboBox, QDoubleSpinBox, QSpinBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog, inline
from chemsuite.protein.dock import box as box_tools
from chemsuite.protein.dock import vina
from chemsuite.protein.dock.box import PADDING, Box
from chemsuite.protein.structure.model import Structure

from . import dock
from .prepare_dialog import proteins, shown_protein, structure_of
from .tool import dialog_parent, protein_workspace

#: where the search box comes from
FOCUSED = "focused"
SELECTED = "selected"
WHOLE = "whole"
TYPED = "typed"
SOURCES = (
    (FOCUSED, "Focused Ligand"),
    (SELECTED, "Selection"),
    (WHOLE, "Whole Protein"),
    (TYPED, "Typed"),
)
SCORING = ((vina.VINA, "Vina"), (vina.VINARDO, "Vinardo"))
_last = vina.Options()  # remembered for the session


class DockDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(dialog_parent(ctx), "Dock Ligands", guide="dock")
        self.ctx = ctx

        inputs = self.section("Input")
        self.protein = QComboBox()
        for id, name in proteins(ctx):
            self.protein.addItem(name, id)
        chosen = shown_protein(ctx)
        if chosen is not None and self.protein.findData(chosen) >= 0:
            self.protein.setCurrentIndex(self.protein.findData(chosen))
        inputs.row("Protein", self.protein)
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row(
            "Ligands",
            self.structure,
            help="A row whose structure has no 3D coordinates gets one conformer first.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx, "selected" if len(ctx.selection) else "all"))

        where = self.section("Box")
        self.source = Segmented(SOURCES, FOCUSED)
        where.row(
            "Around",
            self.source,
            help="The ligand the viewer is focused on, whatever is selected there, the whole "
            "protein, or a centre and size of your own.",
        )
        self.padding = QDoubleSpinBox()
        self.padding.setRange(0.0, 20.0)
        self.padding.setSingleStep(0.5)
        self.padding.setSuffix(" Å")
        self.padding.setValue(PADDING)
        where.row("Padding", self.padding, help="Room left around those atoms on every side.")
        self.centre = [_number(-999.0, 999.0) for _ in range(3)]
        where.row("Centre", inline(*self.centre))
        self.size = [_number(1.0, 100.0, 20.0) for _ in range(3)]
        where.row("Size", inline(*self.size))

        method = self.section("Method")
        self.scoring = Segmented(SCORING, _last.scoring)
        method.row(
            "Scoring",
            self.scoring,
            help="Vina's own function, or Vinardo (Quiroga & Villarreal, 2016), which scores "
            "hydrophobic contacts and hydrogen bonds differently.",
        )
        self.exhaustiveness = QSpinBox()
        self.exhaustiveness.setRange(1, 64)
        self.exhaustiveness.setValue(_last.exhaustiveness)
        method.row(
            "Exhaustiveness",
            self.exhaustiveness,
            help="How hard the search tries. Doubling it roughly doubles the time.",
        )
        self.poses = QSpinBox()
        self.poses.setRange(1, 20)
        self.poses.setValue(_last.poses)
        method.row("Poses", self.poses, help="At most this many per ligand, best first.")

        more = self.section("More Options", collapsed=True)
        self.energy_range = QDoubleSpinBox()
        self.energy_range.setRange(0.5, 20.0)
        self.energy_range.setSingleStep(0.5)
        self.energy_range.setSuffix(" kcal/mol")
        self.energy_range.setValue(_last.energy_range)
        more.row(
            "Energy range",
            self.energy_range,
            help="Keep only poses this close to the best one found.",
        )
        self.seed = QSpinBox()
        self.seed.setRange(0, 2**31 - 1)
        self.seed.setValue(_last.seed)
        more.row("Seed", self.seed, help="The same seed docks the same way every time.")

        self.protein.currentIndexChanged.connect(lambda *_: self._update())
        self.source.changed.connect(lambda *_: self._update())
        self.padding.valueChanged.connect(lambda *_: self._update())
        self.rows.changed.connect(lambda *_: self.refresh())
        self._update()

    # -- the box ------------------------------------------------------------------------------
    def _atoms(self) -> tuple[Structure | None, np.ndarray | None]:
        """The atoms the box is built around: the structure and a mask, or (structure, None)
        when the source has nothing to offer."""
        protein = self.protein.currentData()
        s = structure_of(self.ctx, protein) if isinstance(protein, str) else None
        if s is None:
            return None, None
        source = self.source.value()
        if source == WHOLE:
            return s, np.ones(s.n_atoms, bool)
        workspace = protein_workspace(self.ctx)
        if not workspace.started:
            return s, None
        session = workspace.session()
        if source == FOCUSED:
            focused = session.focused_residue()
            if focused is None or focused[0] != protein:
                return s, None
            return s, s.residue_of == focused[1]
        shown = session.display(protein)
        if shown is None or not shown.selected.any():
            return s, None
        return s, np.asarray(shown.selected, bool)

    def box(self) -> Box | None:
        """The search box as the dialog stands, or None when it cannot be worked out yet."""
        if self.source.value() == TYPED:
            cx, cy, cz = (float(field.value()) for field in self.centre)
            sx, sy, sz = (float(field.value()) for field in self.size)
            return Box((cx, cy, cz), (sx, sy, sz))
        s, atoms = self._atoms()
        if s is None or atoms is None or not atoms.any():
            return None
        return box_tools.from_atoms(s, atoms, padding=float(self.padding.value()))

    def _update(self) -> None:
        typed = self.source.value() == TYPED
        box = self.box()
        if box is not None and not typed:  # show what the chosen atoms give
            for field, value in zip(self.centre, box.centre, strict=True):
                field.setValue(value)
            for field, value in zip(self.size, box.size, strict=True):
                field.setValue(value)
        for field in (*self.centre, *self.size):
            field.setReadOnly(not typed)
            field.setEnabled(typed)
        self.padding.setEnabled(not typed)
        self.refresh()

    # -- running ------------------------------------------------------------------------------
    def options(self) -> vina.Options:
        return vina.Options(
            scoring=self.scoring.value(),
            exhaustiveness=int(self.exhaustiveness.value()),
            poses=int(self.poses.value()),
            energy_range=float(self.energy_range.value()),
            seed=int(self.seed.value()),
        )

    def problem(self) -> str:
        if self.protein.count() == 0:
            return "There are no proteins in this document"
        if not self.structure.currentText():
            return "There is no structure column"
        if not len(self.rows.oids()):
            return "There are no rows to dock"
        trouble = dock.docking_problem(structure_of(self.ctx, self.protein.currentData()))
        if trouble:
            return trouble
        if self.box() is None:
            source = dict(SOURCES)[self.source.value()]
            return f"There is nothing to put a box around ({source.lower()})"
        return ""

    def verb(self) -> str:
        return f"Dock {len(self.rows.oids()):,} Rows"

    def run(self):
        global _last
        _last = self.options()
        protein = self.protein.currentData()
        s = structure_of(self.ctx, protein)
        box = self.box()
        if not isinstance(protein, str) or s is None or box is None:
            return None
        return dock.dock_rows(
            self.ctx, protein, s, self.structure.currentText(), self.rows.oids(), box, _last
        )


class RescoreDialog(ToolDialog):
    """Score poses that already exist, where they stand."""

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(dialog_parent(ctx), "Rescore Poses", guide="rescore-poses")
        self.ctx = ctx

        inputs = self.section("Input")
        self.protein = QComboBox()
        for id, name in proteins(ctx):
            self.protein.addItem(name, id)
        chosen = shown_protein(ctx)
        if chosen is not None and self.protein.findData(chosen) >= 0:
            self.protein.setCurrentIndex(self.protein.findData(chosen))
        inputs.row("Protein", self.protein)
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Poses", self.structure, help="Each one is scored exactly where it sits.")
        self.rows = inputs.row("Rows", RowsChoice(ctx, "selected" if len(ctx.selection) else "all"))

        method = self.section("Method")
        self.scoring = Segmented(SCORING, vina.VINA)
        method.row("Scoring", self.scoring)
        self.protein.currentIndexChanged.connect(lambda *_: self.refresh())
        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if self.protein.count() == 0:
            return "There are no proteins in this document"
        if not self.structure.currentText():
            return "There is no structure column"
        if not len(self.rows.oids()):
            return "There are no rows to score"
        return dock.docking_problem(structure_of(self.ctx, self.protein.currentData()))

    def verb(self) -> str:
        return f"Rescore {len(self.rows.oids()):,} Rows"

    def run(self):
        protein = self.protein.currentData()
        s = structure_of(self.ctx, protein)
        if not isinstance(protein, str) or s is None:
            return None
        options = vina.Options(scoring=self.scoring.value())
        return dock.dock_rows(
            self.ctx,
            protein,
            s,
            self.structure.currentText(),
            self.rows.oids(),
            None,
            options,
            rescore=True,
        )


def _number(low: float, high: float, value: float = 0.0) -> QDoubleSpinBox:
    field = QDoubleSpinBox()
    field.setRange(low, high)
    field.setDecimals(2)
    field.setSuffix(" Å")
    field.setValue(value)
    return field


def open_dock(ctx: ToolContext):
    return ctx.show_dialog(DockDialog(ctx))


def open_rescore(ctx: ToolContext):
    return ctx.show_dialog(RescoreDialog(ctx))
