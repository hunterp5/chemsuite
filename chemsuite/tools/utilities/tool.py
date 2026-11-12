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

"""Tools (utilities): everyday helpers in the spirit of MolManager's Tools ▸ Utilities."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

MENU = ("Tools",)


def _calculator(ctx: ToolContext):
    from .calculator import open_calculator

    return open_calculator(ctx)


def _dialog(name: str):
    def open_(ctx: ToolContext):
        from . import random_dialogs

        dialog = getattr(random_dialogs, name)(ctx)
        ctx.show_dialog(dialog)
        return dialog

    return open_


TOOLS = (
    Tool(
        id="calculator",
        title="Calculator…",
        menu=MENU,
        open=_calculator,
        shortcut="Ctrl+Shift+=",
        tip="A new column from a formula: math, text, statistics and chemistry functions "
        "over other columns, e.g. 9 - log10([IC50 (nM)]) or MW([Structure]).",
        separator_after=True,
    ),
    Tool(
        id="random_numbers",
        title="Random Numbers…",
        menu=MENU,
        open=_dialog("RandomNumbersDialog"),
        tip="Fill a column with random numbers: uniform, whole numbers, normal, log-normal or "
        "a random order, with an optional seed.",
    ),
    Tool(
        id="random_rows",
        title="Random Rows…",
        menu=MENU,
        open=_dialog("RandomRowsDialog"),
        tip="Pick a random sample of rows (a number or a percentage) and select them or save "
        "them as a list.",
    ),
)
