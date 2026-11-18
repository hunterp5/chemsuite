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
"""Tools ▸ Prepare…: hydrogens, flips, His states and protonation, as one undoable edit.

The work itself is :func:`chemsuite.macromolecule.prepare.pipeline.prepare` (repair, waters,
hydrogens, then flips, His states and rotors at their best). This module is the bridge: the
dialog, the background job, replacing the protein in one step
(:meth:`~chemsuite.macromolecule.workspace.session.WorkbenchSession.edit`) and opening the report
panel. What the preparation did is kept with the document, so the report comes back with it
(:mod:`chemsuite.macromolecule.workspace.reports`).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import ToolDialog
from chemsuite.core.wording import plural
from chemsuite.macromolecule.prepare import pipeline
from chemsuite.macromolecule.structure import store
from chemsuite.macromolecule.structure.model import Structure
from chemsuite.macromolecule.workspace.prepared import PreparedRecord, carried_built
from chemsuite.macromolecule.workspace.reports import PrepareReport, Row
from chemsuite.simulate.spec import MinimiseSettings

from ..pka import model as pka_model
from . import minimise
from .minimise_options import MinimiseFields
from .prepare_options import OptionsFields, describe, settings_of
from .tool import KIND, TOOL_ID, ccd, dialog_parent, macromolecule_workspace

HYDROGEN = 1
#: between the two residues a report line names (a gap, a break, a disulfide)
SPAN = " – "
#: the report line naming the settings changed from the defaults
SETTINGS = "Settings"
_minimise_on = False  # Minimise, remembered for the session (off until M5's benchmark passes)


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

        self.method = method = self.section("Method")
        self.fields = OptionsFields(method, self.section("More Options", collapsed=True))
        self.minimise = method.row("Minimise", QCheckBox())
        self.minimise.setChecked(_minimise_on)
        self.minimising = MinimiseFields(
            self.section("Minimisation", collapsed=True), minimise.last_settings()
        )
        self._minimise_help = (
            "Relax the prepared structure in OpenMM, held near the deposited coordinates. The "
            "protein is left as it was; the prepared and the minimised structures are added "
            "as new proteins (-prep and -min)."
        )
        for box in (self.fields.repair, self.fields.side_chains, self.fields.caps):
            box.toggled.connect(lambda *_: self.refresh())
        self.fields.waters.changed.connect(lambda *_: self.refresh())
        self.minimise.toggled.connect(lambda *_: self.refresh())
        self.minimising.forcefield.changed.connect(lambda *_: self.refresh())
        self.protein.currentIndexChanged.connect(lambda *_: self.refresh())
        self.refresh()

    def refresh(self) -> None:
        why = self.minimise_blocked()
        self.minimise.setEnabled(not why)
        if why:
            self.minimise.setChecked(False)
        self.minimise.setToolTip(why or self._minimise_help)
        self.minimising.set_enabled(self.minimise.isChecked())
        protein = self.protein.currentData()
        s = structure_of(self.ctx, protein) if isinstance(protein, str) else None
        kept = self.fields.waters.value() != pipeline.NONE
        self.minimising.allow_waters(kept and s is not None and minimise.has_waters(s))
        super().refresh()

    def minimise_blocked(self) -> str:
        """Why Minimise cannot be ticked, or ``""``: OpenMM set up, and a preparation that
        leaves nothing missing (a force field has no template for half a residue)."""
        if why := minimise.ready():
            return why
        options = self.fields.options()
        if not (options.repair and options.side_chains and options.caps):
            return "Minimising needs Repair, with Side Chains and Cap Ends on"
        return ""

    def problem(self) -> str:
        if self.protein.count() == 0:
            return "There are no proteins in this document"
        protein = self.protein.currentData()
        if isinstance(protein, str) and (why := minimise.minimised_already(self.ctx, protein)):
            return why.replace(": minimise", ": prepare", 1)
        return ""

    def minimise_settings(self) -> MinimiseSettings | None:
        """The minimisation to run after the preparation, or None."""
        return self.minimising.settings() if self.minimise.isChecked() else None

    def verb(self) -> str:
        name = self.protein.currentText()
        return f"Prepare {name}" if name else "Prepare"

    def options(self) -> pipeline.Options:
        return self.fields.options()

    def run(self):
        global _minimise_on
        protein = self.protein.currentData()
        if not isinstance(protein, str):
            return None
        if self.minimise.isEnabled():
            _minimise_on = self.minimise.isChecked()
        settings = self.minimise_settings()
        if settings is not None:
            minimise.remember(settings)
        return prepare_protein(self.ctx, protein, self.options(), settings)


# -- the proteins -------------------------------------------------------------------------------
def proteins(ctx: ToolContext) -> list[tuple[str, str]]:
    """Every protein in the document, as (id, name)."""
    return [(asset.id, asset.name) for asset in ctx.table.assets.values() if asset.kind == KIND]


def shown_protein(ctx: ToolContext) -> str | None:
    """The protein the viewer is looking at (the focused one, else the first shown), or None
    when the viewer has not been opened in this document."""
    workspace = macromolecule_workspace(ctx)
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
def prepare_protein(
    ctx: ToolContext,
    protein: str,
    options: pipeline.Options,
    minimising: MinimiseSettings | None = None,
):
    """Prepare *protein* in the background and show the report. Returns the job, or None when
    there is no such protein.

    Without *minimising*, the prepared structure replaces the protein in one undoable edit.
    With it, the protein is left as it was: the prepared structure and the minimised one are
    added as new proteins (``X-prep``, ``X-min``) in one undoable step. A minimisation that
    cannot run still adds ``X-prep``, its report saying why it was not minimised."""
    s = structure_of(ctx, protein)
    if s is None:
        return None
    cache = ccd(ctx)
    before = int((s.element == HYDROGEN).sum())
    earlier = macromolecule_workspace(ctx).prepared(protein)
    worker = minimise.simulator() if minimising is not None else None

    def work(job) -> Outcome:
        prepared = pipeline.prepare(s, options, cache, job.cancel_token,
                                    ligand_model=pka_model.ligand_pka(options))  # fmt: skip
        if worker is None or minimising is None:
            return Outcome(prepared)
        record = record_of(protein, prepared, options, earlier)
        try:
            relaxed = minimise.relax(worker, prepared.structure, record, minimising, cache,
                                     report=job.report, cancel=job.cancel_token)  # fmt: skip
        except minimise.FAILURES as error:
            return Outcome(prepared, failure=error)
        return Outcome(prepared, relaxed)

    def done(outcome: Outcome) -> None:
        if structure_of(ctx, protein) is None:  # removed while it ran
            ctx.show_message(f"{protein} is no longer in this document")
            return
        prepared = outcome.prepared
        if minimising is None:
            workspace = macromolecule_workspace(ctx)
            workspace.session().edit(f"Prepare {protein}", {protein: prepared.structure})
            report = summary(protein, prepared, options)
            workspace.set_report(report)
            workspace.set_prepared(record_of(protein, prepared, options, earlier))
            show_report(ctx, report)
        else:
            report = _add_copies(ctx, protein, options, outcome, earlier)
        added = int((prepared.structure.element == HYDROGEN).sum()) - before
        warnings = f", {report.warnings} to check" if report.warnings else ""
        ctx.show_message(
            f"Prepared {report.protein}: {plural(added, 'hydrogen')} added, "
            f"net charge {prepared.net_charge:+d}{warnings}"
        )

    return ctx.run(f"Preparing {protein}…", work, done)


@dataclass(frozen=True)
class Outcome:
    """A preparation, and the minimisation after it (or why there was none)."""

    prepared: pipeline.Prepared
    relaxed: minimise.Relaxed | None = None
    failure: BaseException | None = None


def _add_copies(
    ctx: ToolContext,
    protein: str,
    options: pipeline.Options,
    outcome: Outcome,
    earlier: PreparedRecord | None,
) -> PrepareReport:
    """Add the prepared protein, and the minimised one when there is one, in one step; keep
    their records and reports, and show the last. Returns the report shown."""
    workspace = macromolecule_workspace(ctx)
    meta = ctx.table.assets[protein].meta
    prepared, relaxed = outcome.prepared, outcome.relaxed
    prep = minimise.copy_id(ctx, protein, minimise.PREPARED)
    ops = [minimise.copy_op(prep, protein, prepared.structure, minimise.PREPARED, meta)]
    made = None
    if relaxed is not None:
        made = minimise.copy_id(ctx, protein, minimise.MINIMISED, {prep})
        ops.append(minimise.copy_op(made, prep, relaxed.structure, minimise.MINIMISED, meta))
    label = f"Prepare and minimise {protein}" if made else f"Prepare {protein} as {prep}"
    ctx.edit(label, ops)
    record = record_of(prep, prepared, options, earlier)
    workspace.set_prepared(record)
    report = summary(prep, prepared, options)
    if outcome.failure is not None:
        report = replace(report, rows=(minimise.failed_row(outcome.failure), *report.rows))
    workspace.set_report(report)
    if made is not None and relaxed is not None:
        workspace.set_prepared(minimise.record_of(made, prep, record, relaxed))
        report = minimise.report_of(made, relaxed, options.ph)
        workspace.set_report(report)
    show_report(ctx, report)
    return report


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
    from chemsuite.macromolecule.ui.report import PrepareReportPanel

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
    report = macromolecule_workspace(ctx).report(protein)
    return show_report(ctx, report) if report is not None else None


def select(ctx: ToolContext, protein: str, residue: int, line: str = "") -> None:
    """Select *residue* in the workbench and bring it into view (a report line was clicked).

    *line* is the line's label: one that names two residues (``A/LYS 10 – A/GLY 15``, a gap,
    a break or a disulfide) frames both, so the gap between them is seen; the first is
    selected."""
    session = macromolecule_workspace(ctx).session()
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
    from chemsuite.macromolecule.ui.panel import MacromoleculePanel

    panel = ctx.find_panel(TOOL_ID)
    shown = macromolecule_workspace(ctx).session().display(protein)
    if isinstance(panel, MacromoleculePanel) and shown is not None and shown.visible:
        panel.frame_atoms(shown.xyz[atoms])


def select_named(ctx: ToolContext, protein: str, label: str, line: str = "") -> None:
    """Select the residue the workbench calls *label* (``A/HIS 57``) in *protein*: for a line
    found on another structure than the one shown (Inspect's, which is prepared, with caps and
    built residues put in, while the viewer still shows the file's)."""
    shown = macromolecule_workspace(ctx).session().display(protein)
    if shown is None:
        return
    found = pipeline.residues_named(shown.structure, [label])
    if found:
        select(ctx, protein, found[0], line)


def open_prepare(ctx: ToolContext):
    return ctx.show_dialog(PrepareDialog(ctx))
