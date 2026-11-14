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

"""The fragment panel: each distinct fragment, how many compounds hold it, and their activity."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QHBoxLayout, QPushButton, QVBoxLayout

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.sar.fragments import FragmentStat
from chemsuite.tools.sar.common import ResultTable

NO_ACTIVITY = "(none)"


class FragmentPanel(ToolPanel):
    """One row per distinct fragment, with the compounds it came from."""

    tool_id = "fragments"

    def __init__(
        self,
        ctx: ToolContext,
        oids,
        stats: list[FragmentStat],
        method: str,
        activity: str = NO_ACTIVITY,
    ) -> None:
        super().__init__()
        self.ctx = ctx
        self.oids = np.asarray(oids)
        self.stats = stats
        self.method = method
        self.activity = activity
        plain = activity == NO_ACTIVITY
        headers = ["Fragment", "Compounds"] if plain else [
            "Fragment",
            "Compounds",
            f"Mean {activity}",
            "SD",
        ]  # fmt: skip
        self.table = ResultTable(headers)
        # fragment_stats gives the commonest first, and ResultTable keeps that until a header
        # is clicked
        self.table.fill(
            [[s.smiles, s.count] if plain else [s.smiles, s.count, s.mean, s.sd] for s in stats],
            stats,
        )
        select = QPushButton("Select Compounds")
        select.setToolTip("Select the compounds holding the chosen fragments.")
        select.clicked.connect(lambda: self.select_rows())
        copy = QPushButton("Copy Fragments")
        copy.clicked.connect(lambda: self.table.copy_text())
        buttons = QHBoxLayout()
        buttons.addWidget(select)
        buttons.addStretch(1)
        buttons.addWidget(copy)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)

    def title(self) -> str:
        return f"{self.method.upper()} Fragments"

    def select_rows(self) -> list[int]:
        """Select the compounds holding the chosen fragments."""
        chosen = [s for s in self.table.selected_data() if s is not None]
        rows = sorted({r for s in chosen for r in s.rows})
        oids = self.oids[rows].tolist() if rows else []
        if oids:
            self.ctx.selection.set(oids, "fragments")
        return oids

    def state(self) -> dict:
        return {"method": self.method, "activity": self.activity, "oids": self.oids.tolist()}

    def close_panel(self) -> None:
        pass
