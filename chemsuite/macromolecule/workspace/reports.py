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
"""What a preparation did, kept with the document so its report comes back with it.

:class:`PrepareReport` is what the Prepare report panel shows: one
:class:`Row` per thing found or decided, each naming the residue it concerns so the panel can
select it in the workbench. It holds no structure and no :mod:`~chemsuite.macromolecule.prepare`
types, so the workspace can save and reload it without preparing anything.

The saved form is plain JSON (:meth:`PrepareReport.to_json`), under ``reports`` of the protein
workspace's state. A report whose protein has gone is dropped on the way out.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Row:
    """One line of a report: what was found, and where."""

    label: str  # the residue (or residues) as the viewer names them
    text: str
    warn: bool = False  # needs the user's eye (⚠)
    residue: int = -1  # the residue it concerns (-1 for none)

    def to_json(self) -> list[Any]:
        return [self.label, self.text, self.warn, self.residue]

    @classmethod
    def from_json(cls, saved: Any) -> Row | None:
        if not isinstance(saved, list) or len(saved) != 4:
            return None
        label, text, warn, residue = saved
        if not isinstance(label, str) or not isinstance(text, str):
            return None
        if not isinstance(residue, int) or isinstance(residue, bool):
            return None
        return cls(label, text, bool(warn), residue)


@dataclass(frozen=True)
class PrepareReport:
    protein: str
    rows: tuple[Row, ...]
    ph: float = 7.0
    net_charge: int = 0

    @property
    def warnings(self) -> int:
        return sum(1 for row in self.rows if row.warn)

    def to_json(self) -> dict[str, Any]:
        return {
            "protein": self.protein,
            "ph": self.ph,
            "net_charge": self.net_charge,
            "rows": [row.to_json() for row in self.rows],
        }

    @classmethod
    def from_json(cls, saved: Any) -> PrepareReport | None:
        """The report *saved* held, or None when it is not one."""
        if not isinstance(saved, dict) or not isinstance(saved.get("protein"), str):
            return None
        rows = saved.get("rows")
        if not isinstance(rows, list):
            return None
        made = [Row.from_json(row) for row in rows]
        if any(row is None for row in made):
            return None
        ph = saved.get("ph")
        charge = saved.get("net_charge")
        return cls(
            saved["protein"],
            tuple(row for row in made if row is not None),
            float(ph) if isinstance(ph, (int, float)) and not isinstance(ph, bool) else 7.0,
            int(charge) if isinstance(charge, int) and not isinstance(charge, bool) else 0,
        )
