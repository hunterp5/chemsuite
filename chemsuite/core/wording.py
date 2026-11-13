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

"""Wording shared by the tools' status bar messages.

Every tool words its messages the same way:

- **While it runs**, the job's title says what it is doing: "Generating conformers…",
  "Calculating UMAP…". A verb in -ing, sentence case, no counts (the status bar adds the
  progress, "(120/500)"), ending in "…".
- **When it is done**, one sentence that starts with what was done, in the past tense:
  "Generated 1,234 conformers for 100 rows in “Structure 3D”". Counts carry thousands
  separators and real plurals (:func:`plural`, never "row(s)"); the column written is named in
  quotes; anything that went wrong follows a semicolon. A message that only reports a problem
  starts "Could not…".
"""

from __future__ import annotations


def plural(n: int, noun: str, many: str | None = None) -> str:
    """``1 row``, ``1,234 rows``; *many* for a noun that does not just take an "s"."""
    return f"{n:,} {noun if n == 1 else (many or noun + 's')}"
