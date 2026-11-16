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
"""Protein ▸ Inspect…: a preparation surveyed, reviewed line by line, then applied.

Prepare… decides everything and applies it at once. Inspect… runs the same preparation
(:mod:`chemsuite.protein.prepare.review`) without committing it and opens the Inspect panel
in the viewer, where each issue can be left alone and each choice overruled. Re-survey solves
again with the hand-made choices held; Apply prepares the structure as reviewed, as one
undoable edit, and turns the review into the protein's preparation report — which says what
was chosen by hand, so the document records how its structure came to be.

What has been decided so far is kept with the document (:mod:`~chemsuite.protein.workspace.
reviews`) until it is applied; reopening the document surveys the structure again and puts
those decisions back.
"""

from __future__ import annotations

from collections.abc import Mapping

from chemsuite.app.document.tools import ToolContext
from chemsuite.core.wording import plural
from chemsuite.protein.prepare import pipeline, review
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.workspace.reports import PrepareReport, Row
from chemsuite.protein.workspace.reviews import Review

from .prepare_dialog import (
    PrepareDialog,
    record_of,
    select_named,
    show_report,
    structure_of,
    summary,
)
from .prepare_options import options_from, settings_of
from .tool import TOOL_ID, ccd, protein_workspace

HYDROGEN = 1


class InspectDialog(PrepareDialog):
    """The same settings as Prepare…; the survey opens for review instead of being applied."""

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx)
        self.setWindowTitle("Inspect Protein")

    def verb(self) -> str:
        name = self.protein.currentText()
        return f"Inspect {name}" if name else "Inspect"

    def run(self):
        protein = self.protein.currentData()
        if not isinstance(protein, str):
            return None
        return inspect_protein(self.ctx, protein, self.options())


# -- surveying ------------------------------------------------------------------------------------
def _survey(
    protein: str,
    s: Structure,
    options: pipeline.Options,
    cache,
    cancel,
    skipped: tuple[str, ...],
    overrides: Mapping[str, str],
    base: review.Survey | None = None,
) -> tuple[review.Survey, review.Survey]:
    """The survey of *s* with nothing held (the proposals), and with what was decided held.
    Decisions are named by the lines of the survey with nothing held, so without a *base* it
    is surveyed once plainly first."""
    if base is None:
        base = review.survey(protein, s, options, cache, cancel)
    if not skipped and not overrides:
        return base, base
    decisions = review.decisions_from(
        base, skip_issues=dict.fromkeys(skipped, True), overrides=overrides
    )
    return base, review.survey(protein, s, options, cache, cancel, decisions=decisions)


def inspect_protein(
    ctx: ToolContext,
    protein: str,
    options: pipeline.Options,
    skipped: tuple[str, ...] = (),
    overrides: Mapping[str, str] | None = None,
    base: review.Survey | None = None,
):
    """Survey *protein* in the background and open the Inspect panel on it (or fill the open
    one). Returns the job, or None when there is no such protein."""
    s = structure_of(ctx, protein)
    if s is None:
        return None
    cache = ccd(ctx)
    held = dict(overrides or {})

    def work(job):
        return _survey(protein, s, options, cache, job.cancel_token, skipped, held, base)

    def done(surveyed: tuple[review.Survey, review.Survey]) -> None:
        if structure_of(ctx, protein) is None:
            ctx.show_message(f"{protein} is no longer in this document")
            return
        proposals, found = surveyed
        panel = show_inspector(ctx, found, skipped, held, base=proposals)
        panel.set_busy(False)
        remember(ctx, panel)
        warned = f", {found.warnings} to check" if found.warnings else ""
        ctx.show_message(f"Surveyed {protein}: {plural(len(found.issues), 'issue')}, "
                         f"{plural(len(found.choices), 'choice')}{warned}")  # fmt: skip

    return ctx.run(f"Surveying {protein}…", work, done)


def show_inspector(
    ctx: ToolContext,
    found: review.Survey,
    skipped: tuple[str, ...] = (),
    overrides: Mapping[str, str] | None = None,
    base: review.Survey | None = None,
):
    """Open the Inspect panel on *found*, or fill the open one (*base*: the survey with
    nothing held). It docks in the viewer's window, beside the structure it describes."""
    from chemsuite.protein.ui.inspector import InspectorPanel

    panel = ctx.find_panel(InspectorPanel.tool_id)
    if isinstance(panel, InspectorPanel):
        panel.fill(found, skipped, overrides, base)
        return panel
    panel = InspectorPanel(
        found,
        select=lambda id, residue, line: select_prepared(ctx, panel, id, residue, line),
        resurvey=lambda skip, held: _resurvey(ctx, panel, skip, held),
        apply=lambda skip, held: apply_review(ctx, panel, skip, held),
    )
    panel.fill(found, skipped, overrides, base)
    panel.decided.connect(lambda: remember(ctx, panel))
    ctx.add_panel(panel, panel.title(), beside=ctx.find_panel(TOOL_ID))
    return panel


def select_prepared(ctx: ToolContext, panel, protein: str, residue: int, line: str = "") -> None:
    """Select the residue an Inspect line names. The line's *residue* indexes the survey's
    prepared structure, which has caps and built residues the viewer's structure does not yet
    have, so the residue is found again by its label."""
    label = panel.survey.prepared.structure.residue_label(residue)
    select_named(ctx, protein, label, line)


def remember(ctx: ToolContext, panel) -> None:
    """Keep what has been decided so far with the document."""
    found = panel.survey
    skipped, overrides = panel.decisions()
    protein_workspace(ctx).set_review(
        Review(found.protein, settings_of(found.options), skipped, overrides)
    )


def _resurvey(ctx: ToolContext, panel, skipped, overrides):
    panel.set_busy(True)
    proposals = panel.base  # the survey with nothing held: the decisions are read against it
    job = inspect_protein(ctx, proposals.protein, proposals.options, skipped, overrides,
                          base=proposals)  # fmt: skip
    if job is None:
        panel.set_busy(False)
    return job


def open_saved(ctx: ToolContext, state: dict):
    """Reopen the Inspect panel with the document: the structure is surveyed again and the
    decisions saved with it put back."""
    protein = state.get("protein")
    if not isinstance(protein, str):
        return None
    saved = protein_workspace(ctx).review(protein)
    if saved is None:
        return None
    return inspect_protein(ctx, protein, options_from(saved.settings), saved.skipped,
                           saved.overrides)  # fmt: skip


# -- applying -------------------------------------------------------------------------------------
def apply_review(ctx: ToolContext, panel, skipped, overrides):
    """Prepare the protein as reviewed, as one undoable edit; the review becomes its report."""
    found: review.Survey = panel.base  # the proposals: read the answers against them
    protein = found.protein
    s = structure_of(ctx, protein)
    if s is None:
        return None
    decisions = review.decisions_from(
        found, skip_issues=dict.fromkeys(skipped, True), overrides=overrides
    )
    options = found.options
    cache = ccd(ctx)
    before = int((s.element == HYDROGEN).sum())
    panel.set_busy(True)

    def work(job):
        return pipeline.prepare(s, options, cache, job.cancel_token, decisions)

    def done(prepared: pipeline.Prepared) -> None:
        panel.set_busy(False)
        if structure_of(ctx, protein) is None:
            ctx.show_message(f"{protein} is no longer in this document")
            return
        workspace = protein_workspace(ctx)
        workspace.session().edit(f"Prepare {protein} (reviewed)", {protein: prepared.structure})
        report = _with_provenance(summary(protein, prepared, options), found, skipped,
                                  overrides, prepared.structure)  # fmt: skip
        workspace.set_report(report)
        workspace.set_prepared(record_of(protein, prepared, options, workspace.prepared(protein)))
        workspace.drop_review(protein)
        show_report(ctx, report)
        panel.close_panel()
        dock = panel.parent()
        if dock is not None:
            dock.close()
        added = int((prepared.structure.element == HYDROGEN).sum()) - before
        by_hand = len(skipped) + len(overrides)
        ctx.show_message(
            f"Prepared {protein}: {plural(added, 'hydrogen')} added, "
            f"net charge {prepared.net_charge:+d}, {plural(by_hand, 'line')} decided by hand"
        )

    return ctx.run(f"Preparing {protein}…", work, done)


def _with_provenance(
    report: PrepareReport,
    found: review.Survey,
    skipped,
    overrides: Mapping[str, str],
    made: Structure,
) -> PrepareReport:
    """*report*, led by what was decided by hand: the record of how this structure was made.
    Each row points at its residue in *made*, the structure applied (a dropped water is not
    in it, and moves the residues after it)."""
    where = {made.residue_label(r): r for r in range(made.n_residues)}
    rows: list[Row] = []
    for issue in found.issues:
        if issue.id in skipped:
            rows.append(Row(issue.label, f"Left as in the file, by hand: {issue.text}",
                            False, where.get(issue.label, -1)))  # fmt: skip
    for id, label in overrides.items():
        choice = found.choice(id)
        if choice is not None:
            proposed = choice.state.label if choice.state else "?"
            rows.append(Row(choice.label, f"Chosen by hand: {choice.what} {label} "
                            f"(proposed: {proposed})", False,
                            where.get(choice.label, -1)))  # fmt: skip
    return PrepareReport(report.protein, (*rows, *report.rows), report.ph, report.net_charge)


def open_inspect(ctx: ToolContext):
    return ctx.show_dialog(InspectDialog(ctx))
