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
from PySide6.QtWidgets import QComboBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import ToolDialog
from chemsuite.core.wording import plural
from chemsuite.protein.prepare import pipeline
from chemsuite.protein.structure import store
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.workspace.prepared import PreparedRecord, carried_built
from chemsuite.protein.workspace.reports import PrepareReport, Row

from .prepare_options import OptionsFields, describe, settings_of
from .tool import KIND, TOOL_ID, ccd, dialog_parent, protein_workspace

HYDROGEN = 1
#: between the two residues a report line names (a gap, a break, a disulfide)
SPAN = " – "
#: the report line naming the settings changed from the defaults
SETTINGS = "Settings"


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

        self.fields = OptionsFields(
            self.section("Method"), self.section("More Options", collapsed=True)
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
        return self.fields.options()

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
        report = summary(protein, prepared, options)
        workspace.set_report(report)
        workspace.set_prepared(record_of(protein, prepared, options, workspace.prepared(protein)))
        show_report(ctx, report)
        added = int((prepared.structure.element == HYDROGEN).sum()) - before
        warnings = f", {report.warnings} to check" if report.warnings else ""
        ctx.show_message(
            f"Prepared {protein}: {plural(added, 'hydrogen')} added, "
            f"net charge {prepared.net_charge:+d}{warnings}"
        )

    return ctx.run(f"Preparing {protein}…", work, done)


def summary(protein: str, prepared: pipeline.Prepared, options: pipeline.Options) -> PrepareReport:
    """What the preparation did, in the compact form the document keeps; led by the settings
    changed from the defaults, so the report says how it was made."""
    rows = tuple(
        Row(finding.label, finding.text, finding.warn, finding.residue)
        for finding in prepared.report
    )
    if changed := describe(options):
        rows = (Row(SETTINGS, changed), *rows)
    return PrepareReport(protein, rows, ph=options.ph, net_charge=prepared.net_charge)


def record_of(
    protein: str,
    prepared: pipeline.Prepared,
    options: pipeline.Options,
    before: PreparedRecord | None = None,
) -> PreparedRecord:
    """What the preparation left besides the structure, kept with the document; the atoms an
    earlier preparation of the same protein built (*before*) are still counted built."""
    s = prepared.structure
    return PreparedRecord.of(
        protein, s, ph=options.ph, settings=settings_of(options), states=prepared.states,
        pka=prepared.pka, population=prepared.population,
        built=carried_built(prepared.built, before, s),
    )  # fmt: skip


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
    panel = PrepareReportPanel(
        report, select=lambda id, residue, label: select(ctx, id, residue, label)
    )
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


def select(ctx: ToolContext, protein: str, residue: int, line: str = "") -> None:
    """Select *residue* in the workbench and bring it into view (a report line was clicked).

    *line* is the line's label: one that names two residues (``A/LYS 10 – A/GLY 15``, a gap,
    a break or a disulfide) frames both, so the gap between them is seen; the first is
    selected."""
    session = protein_workspace(ctx).session()
    shown = session.display(protein)
    if shown is None:
        return
    s = shown.structure
    mask = np.asarray(s.residue_of == residue, dtype=bool)
    if not mask.any():
        return
    session.select({protein: mask})
    around = mask.copy()
    if SPAN in line:
        other = pipeline.residues_named(s, [line.split(SPAN, 1)[1].strip()])
        if other:
            around |= s.residue_of == other[0]
    frame(ctx, protein, around)


def frame(ctx: ToolContext, protein: str, atoms: np.ndarray) -> None:
    """Zoom the viewer onto these atoms of *protein* (not when it is hidden or the viewer is
    closed)."""
    from chemsuite.protein.ui.panel import ProteinPanel

    panel = ctx.find_panel(TOOL_ID)
    shown = protein_workspace(ctx).session().display(protein)
    if isinstance(panel, ProteinPanel) and shown is not None and shown.visible:
        panel.frame_atoms(shown.xyz[atoms])


def select_named(ctx: ToolContext, protein: str, label: str, line: str = "") -> None:
    """Select the residue the workbench calls *label* (``A/HIS 57``) in *protein*: for a line
    found on another structure than the one shown (Inspect's, which is prepared, with caps and
    built residues put in, while the viewer still shows the file's)."""
    shown = protein_workspace(ctx).session().display(protein)
    if shown is None:
        return
    found = pipeline.residues_named(shown.structure, [label])
    if found:
        select(ctx, protein, found[0], line)


def open_prepare(ctx: ToolContext):
    return ctx.show_dialog(PrepareDialog(ctx))
