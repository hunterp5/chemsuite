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

"""Chemistry ▸ Enumerate ▸ Prepare Reagents: reagents as substituents for a core.

The engine is :mod:`chemsuite.analysis.reagents`. A new window holds one row per fragment, so
it can be saved and handed back to R-Groups on a Core, which takes a file per attachment
point; columns hold one fragment per row, beside the reagent they came from.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QWidget

from chemsuite.analysis.reagents import (
    FIT_SAMPLE,
    GROUPS,
    POOL_MIN_ROWS,
    PreparedRows,
    Reagent,
    any_fits,
    by_id,
    in_group,
    prepare_cells,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.core.mols import mol_from_text, parse_structure
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.wording import plural
from chemsuite.tools.reactions.enumerate import new_table

#: the columns a prepared set is written to; Fragment comes first so a file of them is read
#: back as the substituents, not as the reagents
FRAGMENT = "Fragment"
REAGENT = "Reagent"
SITES = "Sites Found"


class ReagentChoice(QComboBox):
    """Pick a kind of reagent, grouped by the bond it makes, as a menu separates its
    sections."""

    picked = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        for group in GROUPS:
            reagents = in_group(group)
            if not reagents:
                continue
            self.addItem(group, None)
            self.model().item(self.count() - 1).setEnabled(False)
            for reagent in reagents:
                self.addItem(f"    {reagent.name}", reagent.id)
                self.setItemData(self.count() - 1, reagent.describe(), Qt.ItemDataRole.ToolTipRole)
        if self.count() > 1:
            self.setCurrentIndex(1)  # the first kind, not its group heading
        self.activated.connect(lambda _index: self.picked.emit())

    def reagent(self) -> Reagent | None:
        """The chosen kind, or None while a group heading is current."""
        return by_id(self.currentData() or "")

    def set_reagent(self, reagent_id: str) -> None:
        index = self.findData(reagent_id)
        if index >= 0:
            self.setCurrentIndex(index)


class PrepareReagentsDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Prepare Reagents", guide="prepare-reagents")
        self.ctx = ctx
        self._built = False

        inputs = self.section("Input")
        self.source = inputs.row(
            "Reagents",
            QComboBox(),
            help="The structures to prepare: a catalogue of acyl chlorides, amines, boronic "
            "acids and the like. Open a supplier's file with File ▸ Open first.",
        )
        for name in ctx.structure_columns():
            self.source.addItem(name, name)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.choice = method.row(
            "Reagent",
            ReagentChoice(),
            help="What kind of reagent these are. Its leaving group goes, and the atom that "
            "makes the new bond keeps an attachment point.",
        )

        output = self.section("Output")
        self.output = output.row(
            "Result",
            Segmented((("new", "New Window"), ("columns", "Columns")), "new"),
            help="New Window: one row per fragment, to save as a reagent file for R-Groups on "
            "a Core. Columns: the first fragment of each row, beside the rows it came from.",
        )

        more = self.section("More Options", collapsed=True)
        self.clean = more.row(
            "Clean up",
            QCheckBox("Strip salts and neutralise"),
            help="Prepare a hydrochloride or a carboxylate the way its free base or acid "
            "would be prepared.",
        )
        self.clean.setChecked(True)

        self.source.currentIndexChanged.connect(lambda _index: self.refresh())
        self.choice.currentIndexChanged.connect(lambda _index: self.refresh())
        self.rows.changed.connect(lambda *_: self.refresh())
        self.clean.toggled.connect(lambda *_: self.refresh())
        self._built = True
        self.refresh()

    # -- what was chosen --------------------------------------------------------------------
    def chosen_column(self) -> str:
        chosen = self.source.currentData()
        return str(chosen) if chosen else ""

    def _sample(self, column: str, oids: np.ndarray) -> list:
        texts = self.ctx.distinct_values(column, oids)[:FIT_SAMPLE]
        return [mol_from_text(text) for text in texts]

    # -- ToolDialog -------------------------------------------------------------------------
    def refresh(self) -> None:
        if self._built:
            super().refresh()

    def problem(self) -> str:
        column = self.chosen_column()
        if not column:
            return "There are no structures to prepare."
        reagent = self.choice.reagent()
        if reagent is None:
            return "Choose a kind of reagent."
        oids = self.rows.oids()
        if not len(oids):
            return "There are no rows to take reagents from."
        mols = self._sample(column, oids)
        if mols and not any_fits(reagent.id, mols, clean=self.clean.isChecked()):
            return f"Nothing in {column} is a {reagent.name.lower()}."
        return ""

    def verb(self) -> str:
        count = len(self.rows.oids())
        return f"Prepare {count:,} Reagents" if count else "Prepare"

    def run(self):
        ctx = self.ctx
        column = self.chosen_column()
        reagent = self.choice.reagent()
        if reagent is None:
            return None
        oids = self.rows.oids()
        cells = ctx.structure_inputs(column, oids)
        clean = self.clean.isChecked()
        to_columns = self.output.value() == "columns"

        def work(job):
            return prepare_cells(
                reagent.id,
                cells,
                clean=clean,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(rows: PreparedRows) -> None:
            if to_columns:
                finish_columns(ctx, oids[: len(rows)], rows)
            else:
                finish_window(ctx, rows)

        return ctx.run("Preparing reagents…", work, done, keep_partial=True)


def _message(rows: PreparedRows) -> str:
    message = f"Prepared {plural(rows.fragments, 'fragment')} from {plural(len(rows), 'reagent')}"
    if rows.missing:
        message += f"; {plural(rows.missing, 'reagent')} did not match"
    several = sum(1 for row in rows if row.sites > 1)
    if several:
        message += f"; {plural(several, 'reagent')} had more than one site"
    return message


def finish_columns(ctx: ToolContext, oids: np.ndarray, rows: PreparedRows) -> dict[str, str]:
    """The first fragment of each row, and how many sites it had."""
    texts = [row.fragment for row in rows]
    blobs = [parse_structure(text)[0] if text else b"" for text in texts]
    written = ctx.apply_results(
        {
            FRAGMENT: (oids, StructureColumn.from_parsed(blobs, texts)),
            SITES: (oids, [float(row.sites) for row in rows]),
        },
        label="Prepare Reagents",
    )
    ctx.show_message(_message(rows))
    return written


def finish_window(ctx: ToolContext, rows: PreparedRows) -> None:
    """One row per prepared fragment, in a new window to save as a reagent file."""
    fragments: list[str] = []
    reagents: list[str] = []
    counts: list[float] = []
    for row in rows:
        for fragment in row.fragments:
            fragments.append(fragment)
            reagents.append(row.reagent)
            counts.append(float(row.sites))
    if not fragments:
        ctx.show_message(
            f"No fragment came from {plural(len(rows), 'reagent')}: check the kind of reagent "
            "chosen"
        )
        return
    table = new_table(
        {
            FRAGMENT: StructureColumn.from_texts(fragments),
            REAGENT: StructureColumn.from_texts(reagents),
            SITES: counts,
        }
    )
    ctx.open_table(table, f"{len(fragments):,} fragments from {len(rows):,} reagents")
    ctx.show_message(_message(rows))


def open_prepare(ctx: ToolContext):
    return ctx.show_dialog(PrepareReagentsDialog(ctx))
