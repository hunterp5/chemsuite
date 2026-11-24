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
"""Data ▸ QSAR ▸ Auto QSAR…, Apply Model…, Refine Model…"""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

QSAR_MENU = ("Data", "QSAR")
#: every candidate here needs these, and the engine imports them lazily
NEEDS = ("rdkit", "sklearn", "scipy")


def _fit(ctx: ToolContext):
    from .dialog import AutoQsarDialog

    dialog = AutoQsarDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _apply(ctx: ToolContext):
    from .apply_dialog import ApplyModelDialog

    dialog = ApplyModelDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _refine(ctx: ToolContext):
    from .refine_dialog import RefineModelDialog

    dialog = RefineModelDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _panel(ctx: ToolContext, state: dict | None = None):
    from .panel import ModelsPanel

    panel = ModelsPanel(ctx)
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


TOOLS = (
    Tool(
        id="auto_qsar",
        title="Auto QSAR…",
        menu=QSAR_MENU,
        open=_fit,
        requires=NEEDS,
        tip="Fit the best model it can for a property or an activity: it picks the features "
        "out of descriptors and fingerprints, picks and tunes the model, and scores itself "
        "on compounds it did not learn from.",
    ),
    Tool(
        id="qsar_apply",
        title="Apply Model…",
        menu=QSAR_MENU,
        open=_apply,
        requires=("rdkit",),
        tip="Predict with a model this document already holds, or one saved to a file.",
    ),
    Tool(
        id="qsar_refine",
        title="Refine Model…",
        menu=QSAR_MENU,
        open=_refine,
        requires=NEEDS,
        tip="Fit a model again on its own rows plus the data measured since.",
    ),
    Tool(
        id="qsar_models",
        title="QSAR Models",
        menu=(),  # a fit opens it; this entry reopens it with the document
        open=lambda ctx: None,
        restore=_panel,
        tip="The models this document holds, how well each one scored, and what it used.",
    ),
)
