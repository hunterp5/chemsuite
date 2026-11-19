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
"""PCA loadings: what each number column contributes to each principal component, as a
sortable table beside the map. Saved with the document."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from PySide6.QtWidgets import QVBoxLayout

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.widgets.result_table import ResultTable

TOOL_ID = "pca_loadings"


class LoadingsPanel(ToolPanel):
    tool_id = TOOL_ID

    def __init__(self, columns: Sequence[str], loadings: np.ndarray) -> None:
        super().__init__()
        self.columns = list(columns)
        self.loadings = np.asarray(loadings, dtype=np.float64)
        k = self.loadings.shape[1]
        self.table = ResultTable(["Column", *(f"PC{i + 1}" for i in range(k))])
        self.table.fill(
            [
                [name, *(round(float(v), 4) for v in row)]
                for name, row in zip(self.columns, self.loadings, strict=True)
            ]  # fmt: skip
        )
        self.table.setToolTip(
            "Each column's weight on each component: the larger its size, the more the column "
            "moves compounds along that axis; the sign gives the direction."
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.table)

    def title(self) -> str:
        return "PCA Loadings"

    def state(self) -> dict:
        return {"columns": self.columns, "loadings": self.loadings.round(6).tolist()}


def open_loadings(ctx: ToolContext, columns: Sequence[str], loadings: np.ndarray):
    panel = LoadingsPanel(columns, loadings)
    ctx.add_panel(panel, panel.title())
    return panel


def restore(ctx: ToolContext, state: dict):
    columns = state.get("columns")
    loadings = state.get("loadings")
    if not isinstance(columns, list) or not isinstance(loadings, list) or not loadings:
        return None
    values = np.asarray(loadings, dtype=np.float64)
    if values.ndim != 2 or len(values) != len(columns):
        return None
    return open_loadings(ctx, [str(c) for c in columns], values)
