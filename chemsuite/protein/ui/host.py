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

"""What the workbench asks of whatever hosts it, so it never knows what that is.

A table document hosts it through its bridge (:mod:`chemsuite.tools.protein`); later the
workbench runs on its own, with no table and so no :class:`LinkedRow`.
"""

from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from .session import ProteinSession


class LinkedRow(Protocol):
    """The host's current row, as the workbench follows it: a row may name a protein, one of
    its ligands, or a pose of its own."""

    #: the columns whose cells say what a row shows: an edit to one follows the row again
    columns: Collection[str]

    def follow(self, session: ProteinSession, oid: int | None) -> bool:
        """Focus what row *oid* names in a protein *session* shows (its pose, or its ligand).
        False when it names nothing shown there."""
        ...

    def missing(self, session: ProteinSession, oid: int | None) -> str:
        """The protein row *oid* names when the document has no such protein, else ""."""
        ...
