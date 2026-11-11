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

"""What the reaction tools share: the current row's reaction, and reaction cells."""

from __future__ import annotations

from collections.abc import Callable

from chemsuite.app.widgets.reaction_text import ReactionText as ReactionText


def current_reaction(ctx, column: Callable[[], str]) -> Callable[[], str]:
    """A function giving the reaction text of *column()* in the table's current row."""

    def text() -> str:
        name, oid = column(), ctx.current_oid()
        table = ctx.table
        if not name or name not in table or oid is None or oid not in table.oids:
            return ""
        return table.cell(name, oid)

    return text


def cells(table, column: str, oids):
    """Cells to send to workers: the parsed blob, else the text (workers parse it)."""
    data = table.column(column)
    slots = table.slots_of(oids)
    return [
        blob if blob is not None else (text or None)
        for blob, text in zip(data.blobs(slots), data.displays(slots), strict=True)
    ]
