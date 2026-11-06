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
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QFormLayout, QLabel, QVBoxLayout

from chemsuite.app.tools import ToolContext
from chemsuite.sar.rgroups import Decomposition
from chemsuite.tools.similarity.common import scope_combo, scope_oids

from .common import CoreChooser

MATCH_COLUMN = "Core Match"


class RGroupDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("R-Group Decomposition")
        self.ctx = ctx
        self.job = None
        self.core = CoreChooser(ctx)
        self.rows = scope_combo(ctx, "shown")
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Decompose")
        form = QFormLayout()
        form.addRow("Rows:", self.rows)
        layout = QVBoxLayout(self)
        layout.addWidget(self.core)
        layout.addLayout(form)
        layout.addWidget(
            QLabel("Writes a column per position (R1, R2…). A missing substituent is H.")
        )
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.resize(620, 0)
        self.core.core.textChanged.connect(lambda *_: self._update())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self._update()

    def _update(self) -> None:
        problem = (
            self.core.problem() if self.core.core.text().strip() else "Enter or suggest a core."
        )
        self.status.setText(problem)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not problem)

    def accept(self) -> None:
        if self.core.problem():
            return
        self.job = self.decompose()
        super().accept()

    def decompose(self):
        oids = scope_oids(self.ctx, self.rows.currentData())
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
        names = ctx.apply_results(results, label="R-Group Decomposition", overwrite=True)
        count = int(d.matched.sum())
        varying = ", ".join(names[p] for p in d.varying()) or "none"
        ctx.show_message(
            f"{count:,} of {len(oids):,} row(s) match the core; varying positions: {varying}",
            10_000,
        )
        return names
