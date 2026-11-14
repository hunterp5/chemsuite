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

"""Chemistry ▸ Decompose ▸ Central Ring: R-groups around each compound's own core."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QComboBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.wording import plural
from chemsuite.sar.central_ring import POOL_MIN_ROWS, CentralRings, central_ring_rows
from chemsuite.tools.sar.common import structure_cells

CORE_COLUMN = "Core"


class CentralRingDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Central Ring", guide="central-ring")
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row(
            "Structures",
            self.structure,
            help="Each compound is split around its own most central ring system, so no core "
            "is needed. The Core column says which one each row got.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        self.structure.currentTextChanged.connect(lambda *_: self.refresh())
        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to decompose."
        if not self._rows():
            return "There are no rows to decompose."
        return ""

    def verb(self) -> str:
        return f"Decompose {self._rows():,} Rows"

    def run(self):
        ctx = self.ctx
        oids = self.rows.oids()
        cells = structure_cells(ctx, self.structure.currentText(), oids)

        def work(job):
            return central_ring_rows(
                cells,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(found: CentralRings) -> None:
            finish(ctx, oids[: len(found)], found)

        return ctx.run("Finding central rings…", work, done, keep_partial=True)


def finish(ctx: ToolContext, oids: np.ndarray, found: CentralRings) -> dict[str, str]:
    """Write the core and R-group columns for the rows *found* covers."""
    results: dict[str, tuple[np.ndarray, object]] = {
        CORE_COLUMN: (oids, StructureColumn.from_texts(found.cores))
    }
    for position in found.positions():
        results[position] = (oids, found.labels[position])
    written = ctx.apply_results(results, label="Central Ring")
    cores = len({core for core in found.cores if core})
    message = (
        f"Decomposed {plural(found.decomposed(), 'compound')} around "
        f"{plural(cores, 'distinct core')}"
    )
    if found.failed:
        message += f"; {found.failed:,} could not be read"
    ctx.show_message(message)
    return written
