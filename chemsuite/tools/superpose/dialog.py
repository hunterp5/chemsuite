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
"""The Superpose dialog: rows onto a reference (a row or a protein's ligand), or each row's
conformers onto its first. The reference is kept with the document as an asset that each
superposed structure names, so the Browser can show it."""

from __future__ import annotations

import uuid

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QLineEdit,
    QPushButton,
    QSpinBox,
)
from rdkit import Chem

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog, inline
from chemsuite.core.assets import Asset, PutAsset
from chemsuite.core.table.columns import NumberColumn, StructureColumn
from chemsuite.core.wording import plural
from chemsuite.superpose import model as sp
from chemsuite.superpose.run import align_conformers_rows, mol_of_cell, superpose_rows
from chemsuite.superpose.substructure import parse_core

#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 20
#: the asset kind of a superposition's reference (a MOL block)
REFERENCE_KIND = "superpose-reference"
SUFFIX = " superposed"
BEST_CONFORMER = "Best Conformer"
ONTO_REFERENCE = "reference"
ONTO_FIRST = "conformers"
CURRENT_ROW = "row"
METHOD_TITLES = {
    sp.MCS: "Substructure (MCS)",
    sp.CORE: "Substructure (core)",
    sp.O3A: "Atoms (O3A)",
    sp.SHAPE: "Shape",
    sp.COLOR: "Shape + Features",
    sp.FIELD: "Field (Electrostatics + Shape)",
}
METHOD_HELP = (
    "Substructure: fit the atoms of a common core, the largest one (MCS) or one you give. "
    "Atoms: Open3DALIGN pairs atoms by type without a core. Shape: overlap of the "
    "molecules' volumes. Shape + Features: volume plus donors, acceptors, rings and "
    "charges (Combo is the two Tanimotos added, 0–2). Field: electrostatic potential "
    "overlap plus shape."
)

_last: tuple[str, sp.Settings] | None = None  # (mode, settings), remembered for the session


def output_name(structure: str) -> str:
    return structure if structure.endswith(SUFFIX) else structure + SUFFIX


class SuperposeDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Superpose", guide="superpose")
        self.ctx = ctx
        self._refs = self._ligand_refs()
        self._reference_cache: tuple[object, Chem.Mol | None] | None = None
        names = ctx.structure_columns()
        three_d = [n for n in names if n.endswith(" 3D")]
        default = three_d[-1] if three_d else ctx.default_structure_column()

        inputs = self.inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(names)
        if default:
            self.structure.setCurrentText(default)
        inputs.row(
            "Structures",
            self.structure,
            help="Each row's stored 3D conformers are all tried and the best fit kept; a row "
            "with only 2D gets conformers first.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx))
        self.mode = QComboBox()
        self.mode.addItem("Onto a reference", ONTO_REFERENCE)
        self.mode.addItem("Each row's conformers onto its first", ONTO_FIRST)
        inputs.row("Superpose", self.mode)
        self.reference = QComboBox()
        self.reference.addItem("Current row", CURRENT_ROW)
        for index, ref in enumerate(self._refs):
            self.reference.addItem(ref.title, index)
        inputs.row(
            "Reference",
            self.reference,
            help="The table's current row (click another row to change it), or a ligand of a "
            "protein in the document, in its bound pose.",
        )
        self.reference_column = QComboBox()
        self.reference_column.addItems(names)
        if default:
            self.reference_column.setCurrentText(default)
        inputs.row(
            "Reference from",
            self.reference_column,
            help="The column the reference row's 3D structure is read from (its first conformer).",
        )

        self.method_section = self.section("Method")
        self.method = QComboBox()
        for key, title in METHOD_TITLES.items():
            self.method.addItem(title, key)
        self.method_section.row("Method", self.method, help=METHOD_HELP)
        self.core = QLineEdit()
        self.core.setPlaceholderText("SMILES or SMARTS")
        self.core.setClearButtonEnabled(True)
        draw = QPushButton("Draw…")
        draw.clicked.connect(self.draw_core)
        self.core_row = inline(self.core, draw, stretch=0)
        self.method_section.row("Core", self.core_row)

        more = self.section("More Options", collapsed=True)
        self.ring_only = QCheckBox()
        self.complete_rings = QCheckBox()
        self.elements = QCheckBox()
        self.timeout = QSpinBox()
        self.timeout.setRange(1, 600)
        self.timeout.setSuffix(" s")
        more.row("Rings match rings only", self.ring_only)
        more.row("Complete rings only", self.complete_rings)
        more.row(
            "Compare elements",
            self.elements,
            help="Off: any atom matches any other, so only the bond pattern counts.",
        )
        more.row("MCS time limit", self.timeout, help="Per row: the best core found so far.")
        self.mcs_options = (self.ring_only, self.complete_rings, self.elements, self.timeout)
        self.typing = QComboBox()
        self.typing.addItem("MMFF94 types and charges", "mmff")
        self.typing.addItem("Crippen logP contributions", "crippen")
        more.row("Atom typing", self.typing)
        self.color_weight = QDoubleSpinBox()
        self.color_weight.setRange(0.0, 1.0)
        self.color_weight.setSingleStep(0.1)
        more.row(
            "Feature weight",
            self.color_weight,
            help="The share of the features in the optimization: 0 is shape only, 1 features only.",
        )
        self.charges = QComboBox()
        self.charges.addItem("Gasteiger", "gasteiger")
        self.charges.addItem("MMFF94", "mmff")
        more.row("Charges", self.charges)
        self.field_weight = QDoubleSpinBox()
        self.field_weight.setRange(0.0, 1.0)
        self.field_weight.setSingleStep(0.1)
        more.row(
            "Field weight",
            self.field_weight,
            help="The share of the field in the optimization: 0 is shape only, 1 field only.",
        )
        self.conformers = QSpinBox()
        self.conformers.setRange(1, 300)
        more.row(
            "Conformers for 2D rows",
            self.conformers,
            help="Generated as Generate Conformers does (CONFORGE where installed, else ETKDG).",
        )
        self.heavy_only = QCheckBox()
        more.row("Heavy atoms only", self.heavy_only, help="Fit conformers on heavy atoms only.")
        self.more = more

        mode, settings = _last or (ONTO_REFERENCE, sp.Settings())
        self.mode.setCurrentIndex(max(0, self.mode.findData(mode)))
        self.set_settings(settings)
        for signal in (
            self.structure.currentIndexChanged,
            self.mode.currentIndexChanged,
            self.reference.currentIndexChanged,
            self.reference_column.currentIndexChanged,
            self.method.currentIndexChanged,
            self.rows.changed,
            self.core.textChanged,
        ):
            signal.connect(lambda *_: self._update())
        stop = ctx.watch_current(lambda _oid: self._update())
        self.destroyed.connect(lambda *_: stop())
        self._update()

    def _ligand_refs(self) -> list:
        from chemsuite.tools.macromolecule.ligand_refs import ligand_refs

        return ligand_refs(self.ctx)

    # -- settings -----------------------------------------------------------------------------
    def set_settings(self, settings: sp.Settings) -> None:
        self.method.setCurrentIndex(max(0, self.method.findData(settings.method)))
        self.core.setText(settings.core)
        self.ring_only.setChecked(settings.mcs_ring_matches_ring)
        self.complete_rings.setChecked(settings.mcs_complete_rings)
        self.elements.setChecked(settings.mcs_compare_elements)
        self.timeout.setValue(settings.mcs_timeout)
        self.typing.setCurrentIndex(max(0, self.typing.findData(settings.o3a_typing)))
        self.color_weight.setValue(settings.color_weight)
        self.charges.setCurrentIndex(max(0, self.charges.findData(settings.charges)))
        self.field_weight.setValue(settings.field_weight)
        self.conformers.setValue(settings.conformers)
        self.heavy_only.setChecked(settings.heavy_only)

    def settings(self) -> sp.Settings:
        return sp.Settings(
            method=self.method.currentData(),
            core=self.core.text().strip(),
            mcs_ring_matches_ring=self.ring_only.isChecked(),
            mcs_complete_rings=self.complete_rings.isChecked(),
            mcs_compare_elements=self.elements.isChecked(),
            mcs_timeout=self.timeout.value(),
            o3a_typing=self.typing.currentData(),
            color_weight=self.color_weight.value(),
            charges=self.charges.currentData(),
            field_weight=self.field_weight.value(),
            conformers=self.conformers.value(),
            heavy_only=self.heavy_only.isChecked(),
        )

    def _onto_reference(self) -> bool:
        return self.mode.currentData() == ONTO_REFERENCE

    def _update(self) -> None:
        onto = self._onto_reference()
        method = self.method.currentData()
        row = self.reference.currentData() == CURRENT_ROW
        inputs = self.inputs
        inputs.set_row_visible(self.reference, onto)
        inputs.set_row_visible(self.reference_column, onto and row)
        self.method_section.setVisible(onto)
        self.method_section.set_row_visible(self.core_row, method == sp.CORE)
        more = self.more
        for option in self.mcs_options:
            more.set_row_visible(option, onto and method == sp.MCS)
        more.set_row_visible(self.typing, onto and method == sp.O3A)
        more.set_row_visible(self.color_weight, onto and method == sp.COLOR)
        more.set_row_visible(self.charges, onto and method == sp.FIELD)
        more.set_row_visible(self.field_weight, onto and method == sp.FIELD)
        more.set_row_visible(self.conformers, onto)
        more.set_row_visible(self.heavy_only, not onto)
        self.refresh()

    # -- the reference ------------------------------------------------------------------------
    def _reference_key(self) -> tuple:
        choice = self.reference.currentData()
        if choice == CURRENT_ROW:
            return (CURRENT_ROW, self.reference_column.currentText(), self.ctx.current_oid())
        return ("ligand", choice)

    def reference_mol(self) -> Chem.Mol | None:
        """The chosen reference as it is now (cached until the choice changes)."""
        key = self._reference_key()
        if self._reference_cache is not None and self._reference_cache[0] == key:
            return self._reference_cache[1]
        if key[0] == CURRENT_ROW:
            _kind, column, oid = key
            mol = None
            if column and oid is not None:
                mol = mol_of_cell(self.ctx.structure_inputs(column, [oid])[0])
        else:
            from chemsuite.tools.macromolecule.ligand_refs import ligand_mol

            mol = ligand_mol(self.ctx, self._refs[key[1]])
        self._reference_cache = (key, mol)
        return mol

    def _reference_asset(self, mol: Chem.Mol, settings: sp.Settings) -> Asset:
        key = self._reference_key()
        if key[0] == CURRENT_ROW:
            _kind, column, oid = key
            position = int(np.flatnonzero(self.ctx.table.oids == oid)[0]) + 1
            name = f"Row {position} ({column})"
            meta = {"source": "row", "oid": str(oid), "column": column}
        else:
            ref = self._refs[key[1]]
            name = ref.title
            meta = {"source": "protein", "protein": ref.protein, "ligand": ref.label}
        meta["method"] = settings.method
        asset_id = self.ctx.table.unique_asset_id(f"superpose-{uuid.uuid4().hex[:8]}")
        block = Chem.MolToMolBlock(sp.first_conformer(mol))
        return Asset(asset_id, REFERENCE_KIND, name, block.encode(), meta)

    def draw_core(self) -> None:
        from chemsuite.app.sketcher import SketcherDialog

        dialog = SketcherDialog(self.core.text(), "Core", self)
        if not dialog.exec() or not dialog.molblock:
            return
        mol = Chem.MolFromMolBlock(dialog.molblock)
        if mol is not None:
            self.core.setText(Chem.MolToSmiles(mol))

    # -- checks -------------------------------------------------------------------------------
    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column."
        if self.rows.problem():
            return self.rows.problem()
        if not len(self.rows.oids()):
            return "There are no rows to work on."
        if not self._onto_reference():
            return ""
        if self.reference.currentData() == CURRENT_ROW and self.ctx.current_oid() is None:
            return "There is no current row: click the reference row in the table."
        mol = self.reference_mol()
        if mol is None:
            return "The reference could not be read."
        if not sp.has_3d(mol):
            column = self.reference_column.currentText()
            return f"The reference row has no 3D structure in “{column}”."
        if self.method.currentData() == sp.CORE:
            core = parse_core(self.core.text())
            if core is None:
                return "Enter or draw the core (SMILES or SMARTS)."
            if not mol.HasSubstructMatch(core):
                return "The core is not in the reference."
        return ""

    def verb(self) -> str:
        return f"Superpose {len(self.rows.oids()):,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        global _last
        settings = self.settings()
        _last = (self.mode.currentData(), settings)
        if self._onto_reference():
            return self.superpose(settings)
        return self.align_conformers(settings.heavy_only)

    def superpose(self, settings: sp.Settings):
        ctx = self.ctx
        name = self.structure.currentText()
        oids = self.rows.oids()
        cells = ctx.structure_inputs(name, oids)
        reference = self.reference_mol()
        asset = self._reference_asset(reference, settings)
        target = output_name(name)

        def work(job):
            return superpose_rows(
                cells,
                reference,
                asset.id,
                settings,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(batch) -> None:
            done_oids = oids[: len(batch.blobs)]  # after a cancel: only the first rows
            results = {
                target: (done_oids, StructureColumn.from_parsed(batch.blobs, batch.displays)),
                BEST_CONFORMER: (done_oids, NumberColumn.from_floats(batch.conformers, 0)),
            }
            for score, values in batch.scores.items():
                decimals = sp.DECIMALS.get(score, 3)
                results[score] = (done_oids, NumberColumn.from_floats(values, decimals))
            names = ctx.apply_results(results, label="Superpose",
                                      extra_ops=[PutAsset(asset.id, None, asset)])  # fmt: skip
            message = (
                f"Superposed {plural(len(done_oids) - batch.failed, 'row')} onto "
                f"{asset.name} in “{names[target]}”"
            )
            if batch.failed:
                message += f"; {plural(batch.failed, 'row')} could not be superposed"
            ctx.show_message(message)

        return ctx.run("Superposing…", work, done, keep_partial=True)

    def align_conformers(self, heavy_only: bool):
        ctx = self.ctx
        name = self.structure.currentText()
        oids = self.rows.oids()
        cells = ctx.structure_inputs(name, oids)
        target = output_name(name)

        def work(job):
            return align_conformers_rows(
                cells,
                heavy_only,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS * 5),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(batch) -> None:
            done_oids = oids[: len(batch.blobs)]
            results = {
                target: (done_oids, StructureColumn.from_parsed(batch.blobs, batch.displays))
            }
            for score, values in batch.scores.items():
                results[score] = (done_oids, NumberColumn.from_floats(values, 3))
            names = ctx.apply_results(results, label="Superpose conformers")
            message = f"Aligned the conformers of {plural(len(done_oids), 'row')} in "
            ctx.show_message(message + f"“{names[target]}”")

        return ctx.run("Aligning conformers…", work, done, keep_partial=True)
