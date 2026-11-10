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
"""R-Group Decomposition…: write a series' substituents as Core, R1, R2… columns."""

from __future__ import annotations

import numpy as np

from chemsuite.app.dialogs import RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.sar.rgroups import Decomposition

from .common import CoreChooser

MATCH_COLUMN = "Core Match"


class RGroupDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "R-Group Decomposition", guide="r-group-decomposition")
        self.ctx = ctx
        inputs = self.section("Input")
        self.core = CoreChooser(ctx, inputs)
        self.rows = inputs.row("Rows", RowsChoice(ctx, "shown"))
        self.core.core.textChanged.connect(lambda *_: self.refresh())
        self.core.structure.currentIndexChanged.connect(lambda *_: self.refresh())
        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.core.core.text().strip():
            return "Enter or suggest a core."
        return self.core.problem()

    def verb(self) -> str:
        return f"Decompose {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.decompose()

    def decompose(self):
        oids = self.rows.oids()
        work = self.core.task(oids)
        # After a cancel, d covers the first rows only (a numbered core decomposes in chunks).
        return self.ctx.run(
            "Decomposing R-groups…",
            work,
            lambda d: self.finish(oids[: len(d)], d),
            keep_partial=True,
        )

    def finish(self, oids: np.ndarray, d: Decomposition) -> dict[str, str]:
        ctx = self.ctx
        matched = d.matched.tolist()
        results = {MATCH_COLUMN: (oids, ["yes" if m else "" for m in matched])}
        for position in d.positions:
            results[position] = (oids, d.labels[position])
        names = ctx.apply_results(results, label="R-Group Decomposition")
        count = int(d.matched.sum())
        varying = ", ".join(names[p] for p in d.varying()) or "none"
        ctx.show_message(
            f"{count:,} of {len(oids):,} row(s) match the core; varying positions: {varying}",
            10_000,
        )
        return names
