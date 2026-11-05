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
"""Data ▸ Statistics…: a panel of summaries, correlations, outliers, curve fits and tests."""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext


def _open(ctx: ToolContext, state: dict | None = None):
    from .panel import StatisticsPanel

    panel = StatisticsPanel(ctx)
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


TOOL = Tool(
    id="statistics",
    title="Statistics…",
    menu=("Data",),
    open=_open,
    restore=_open,
    requires=("scipy",),
    tip="Summaries, correlations, outliers, curve fits and hypothesis tests for number columns.",
)
