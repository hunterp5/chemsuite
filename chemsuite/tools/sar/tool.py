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
"""Data ▸ SAR: R-Group Decomposition, SAR Matrix, Free-Wilson, Hansch, Matched Pairs, Cliffs."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

MENU = ("Data", "SAR")


def _dialog(module: str, name: str):
    def open_dialog(ctx: ToolContext):
        from importlib import import_module

        cls = getattr(import_module(f"{__package__}.{module}"), name)
        dialog = cls(ctx)
        ctx.show_dialog(dialog)
        return dialog

    return open_dialog


def _restore(module: str):
    def restore(ctx: ToolContext, state: dict):
        from importlib import import_module

        return import_module(f"{__package__}.{module}").restore(ctx, state)

    return restore


def _matrix(ctx: ToolContext, state: dict | None = None):
    from .matrix_panel import SarMatrixPanel

    panel = SarMatrixPanel(ctx)
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


TOOLS = (
    Tool(
        id="rgroups",
        title="R-Group Decomposition…",
        menu=("Chemistry",),
        open=_dialog("rgroup_dialog", "RGroupDialog"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Split a series onto a common core, writing the substituents as R1, R2… columns.",
    ),
    Tool(
        id="sar_matrix",
        title="SAR Matrix…",
        menu=MENU,
        open=_matrix,
        restore=_matrix,
        tip="Activity by substituent pairs (R1 × R2), with Free-Wilson estimates for gaps.",
    ),
    Tool(
        id="free_wilson",
        title="Free-Wilson…",
        menu=MENU,
        open=_dialog("free_wilson_ui", "FreeWilsonDialog"),
        restore=_restore("free_wilson_ui"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Additive substituent contributions, and the best combinations not yet made.",
    ),
    Tool(
        id="hansch",
        title="Hansch…",
        menu=MENU,
        open=_dialog("hansch_ui", "HanschDialog"),
        restore=_restore("hansch_ui"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Regress activity on substituent constants (π, σ, steric) from your own table.",
    ),
    Tool(
        id="mmp",
        title="Matched Pairs…",
        menu=MENU,
        open=_dialog("mmp_ui", "MmpDialog"),
        restore=_restore("mmp_ui"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Pairs that differ by one change, grouped by transformation, with activity shifts.",
    ),
    Tool(
        id="sali",
        title="Activity Cliffs…",
        menu=MENU,
        open=_dialog("sali_ui", "SaliDialog"),
        restore=_restore("sali_ui"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Similar compounds with very different activity (structure-activity landscape index).",
    ),
)
