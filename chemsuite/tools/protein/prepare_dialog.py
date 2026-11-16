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
"""Protein ▸ Prepare…: hydrogens, flips, His states and protonation, as one undoable edit.

The work itself is :func:`chemsuite.protein.prepare.pipeline.prepare` (repair, waters,
hydrogens, then flips, His states and rotors at their best). This module is the bridge: the
dialog, the background job, replacing the protein in one step
(:meth:`~chemsuite.protein.workspace.session.WorkbenchSession.edit`) and opening the report
panel. What the preparation did is kept with the document, so the report comes back with it
(:mod:`chemsuite.protein.workspace.reports`).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import Segmented, ToolDialog
from chemsuite.core.wording import plural
from chemsuite.protein.prepare import pipeline, sites
from chemsuite.protein.structure import store
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.workspace.reports import PrepareReport, Row

from .tool import KIND, TOOL_ID, ccd, dialog_parent, protein_workspace

WATERS = ((pipeline.ALL, "All"), (pipeline.NEAR_LIGANDS, "Near Ligands"), (pipeline.NONE, "None"))
TITRATE = ((sites.OFF, "Off"), (sites.POCKET, "Pocket"), (sites.ALL, "All"))
HYDROGEN = 1


class PrepareDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(dialog_parent(ctx), "Prepare Protein", guide="prepare")
        self.ctx = ctx

        inputs = self.section("Input")
        self.protein = QComboBox()
        for id, name in proteins(ctx):
            self.protein.addItem(name, id)
        chosen = shown_protein(ctx)
        if chosen is not None and self.protein.findData(chosen) >= 0:
            self.protein.setCurrentIndex(self.protein.findData(chosen))
        inputs.row("Protein", self.protein)

        method = self.section("Method")
        self.ph = QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setSingleStep(0.5)
        self.ph.setDecimals(1)
        self.ph.setValue(7.0)
        method.row("pH", self.ph, help="Which form each acid, base and ligand takes.")
        self.waters = Segmented(WATERS, pipeline.ALL)
        method.row(
            "Waters",
            self.waters,
            help="Keep every water, only those within 5 Å of a ligand, or none.",
        )
        self.repair = QCheckBox()
        self.repair.setChecked(True)
        method.row(
            "Repair",
            self.repair,
            help="Build missing side-chain atoms and terminal oxygens from the component "
            "templates. No loops are built.",
        )
        self.optimise = QCheckBox()
        self.optimise.setChecked(True)
        method.row(
            "Optimise",
            self.optimise,
            help="Choose each Asn/Gln/His flip, each His tautomer and each hydroxyl and "
            "amine rotor by the hydrogen bonds they make.",
        )

        more = self.section("More Options", collapsed=True)
        self.titrate = Segmented(TITRATE, sites.OFF)
        more.row(
            "Titrate",
            self.titrate,
            help="Work out each group's pKa and take the form it mostly holds at this pH, "
            "instead of the textbook form. Less accurate on average than the textbook forms.",
        )
        self.protein.currentIndexChanged.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if self.protein.count() == 0:
            return "There are no proteins in this document"
        return ""

    def verb(self) -> str:
        name = self.protein.currentText()
        return f"Prepare {name}" if name else "Prepare"

    def options(self) -> pipeline.Options:
        return pipeline.Options(
            ph=float(self.ph.value()),
            waters=self.waters.value(),
            repair=self.repair.isChecked(),
            optimise=self.optimise.isChecked(),
            titrate=self.titrate.value(),
        )

    def run(self):
        protein = self.protein.currentData()
        if not isinstance(protein, str):
            return None
        return prepare_protein(self.ctx, protein, self.options())


# -- the proteins -------------------------------------------------------------------------------
def proteins(ctx: ToolContext) -> list[tuple[str, str]]:
    """Every protein in the document, as (id, name)."""
    return [(asset.id, asset.name) for asset in ctx.table.assets.values() if asset.kind == KIND]


def shown_protein(ctx: ToolContext) -> str | None:
    """The protein the viewer is looking at (the focused one, else the first shown), or None
    when the viewer has not been opened in this document."""
    workspace = protein_workspace(ctx)
    if not workspace.started:
        return None
    session = workspace.session()
    focused = session.focused_residue()
    if focused is not None:
        return focused[0]
    for shown in session.displays:
        if shown.visible:
            return shown.id
    return None


def structure_of(ctx: ToolContext, protein: str) -> Structure | None:
    asset = ctx.table.assets.get(protein)
    if asset is None or asset.kind != KIND:
        return None
    return store.as_structure(asset.payload)


# -- running ------------------------------------------------------------------------------------
def prepare_protein(ctx: ToolContext, protein: str, options: pipeline.Options):
    """Prepare *protein* in the background, apply it as one undoable edit and show the
    report. Returns the job, or None when there is no such protein."""
    s = structure_of(ctx, protein)
    if s is None:
        return None
    cache = ccd(ctx)
    before = int((s.element == HYDROGEN).sum())

    def work(job):
        return pipeline.prepare(s, options, cache, job.cancel_token)

    def done(prepared: pipeline.Prepared) -> None:
        if structure_of(ctx, protein) is None:  # removed while it ran
            ctx.show_message(f"{protein} is no longer in this document")
            return
        workspace = protein_workspace(ctx)
        workspace.session().edit(f"Prepare {protein}", {protein: prepared.structure})
        report = summary(protein, prepared, options.ph)
        workspace.set_report(report)
        show_report(ctx, report)
        added = int((prepared.structure.element == HYDROGEN).sum()) - before
        warnings = f", {report.warnings} to check" if report.warnings else ""
        ctx.show_message(
            f"Prepared {protein}: {plural(added, 'hydrogen')} added, "
            f"net charge {prepared.net_charge:+d}{warnings}"
        )

    return ctx.run(f"Preparing {protein}…", work, done)


def summary(protein: str, prepared: pipeline.Prepared, ph: float) -> PrepareReport:
    """What the preparation did, in the compact form the document keeps."""
    rows = tuple(
        Row(finding.label, finding.text, finding.warn, finding.residue)
        for finding in prepared.report
    )
    return PrepareReport(protein, rows, ph=ph, net_charge=prepared.net_charge)


# -- the report panel ---------------------------------------------------------------------------
def show_report(ctx: ToolContext, report: PrepareReport):
    """Open the report panel on *report*, or fill the open one.

    It docks in the viewer's window, beside the structure it describes; with the viewer
    closed it goes beside the table, so a report is never out of reach."""
    from chemsuite.protein.ui.report import PrepareReportPanel

    panel = ctx.find_panel(PrepareReportPanel.tool_id)
    if isinstance(panel, PrepareReportPanel):
        panel.fill(report)
        return panel
    panel = PrepareReportPanel(report, select=lambda id, residue: select(ctx, id, residue))
    ctx.add_panel(panel, panel.title(), beside=ctx.find_panel(TOOL_ID))
    return panel


def open_report(ctx: ToolContext, state: dict):
    """Reopen the report panel with the document (the report itself was saved with the
    protein workspace)."""
    protein = state.get("protein")
    if not isinstance(protein, str):
        return None
    report = protein_workspace(ctx).report(protein)
    return show_report(ctx, report) if report is not None else None


def select(ctx: ToolContext, protein: str, residue: int) -> None:
    """Select *residue* in the workbench (a report line was clicked)."""
    workspace = protein_workspace(ctx)
    session = workspace.session()
    shown = session.display(protein)
    if shown is None:
        return
    mask = shown.structure.residue_of == residue
    if mask.any():
        session.select({protein: np.asarray(mask, dtype=bool)})


def select_named(ctx: ToolContext, protein: str, label: str) -> None:
    """Select the residue the workbench calls *label* (``A/HIS 57``) in *protein*: for a line
    found on another structure than the one shown (Inspect's, which is prepared, with caps and
    built residues put in, while the viewer still shows the file's)."""
    shown = protein_workspace(ctx).session().display(protein)
    if shown is None:
        return
    s = shown.structure
    found = [r for r in range(s.n_residues) if s.residue_label(r) == label]
    if found:
        select(ctx, protein, found[0])


def open_prepare(ctx: ToolContext):
    return ctx.show_dialog(PrepareDialog(ctx))
