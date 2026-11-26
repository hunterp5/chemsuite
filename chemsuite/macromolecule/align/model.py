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
"""An alignment of chains as the document keeps it.

Each :class:`AlignedRow` is one chain of one protein; its ``positions`` say, for every column,
which position of the chain's sequence (:func:`~..structure.sequence.chain_sequences`, the full
SEQRES sequence where the file has one) sits there, or -1 for a gap. Positions rather than
residue indices are kept so the alignment survives edits that renumber residues. The saved form
is plain JSON (as :mod:`..workspace.poses` is); anything unreadable is dropped, not raised.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

#: how an alignment was made
PAIRWISE, MSA = "pairwise", "msa"


@dataclass(frozen=True, eq=False)
class AlignedRow:
    protein: str  # the protein's asset id
    chain: str  # the chain's name in that protein
    positions: np.ndarray  # int32: sequence position at each column, -1 for a gap

    def to_json(self) -> dict[str, Any]:
        return {"protein": self.protein, "chain": self.chain,
                "positions": self.positions.tolist()}  # fmt: skip

    @classmethod
    def from_json(cls, saved: Any) -> AlignedRow | None:
        if not isinstance(saved, dict):
            return None
        protein, chain, positions = saved.get("protein"), saved.get("chain"), saved.get("positions")
        if not isinstance(protein, str) or not isinstance(chain, str):
            return None
        if not isinstance(positions, list) or not all(isinstance(p, int) for p in positions):
            return None
        return cls(protein, chain, np.asarray(positions, dtype=np.int32))


@dataclass(frozen=True, eq=False)
class Alignment:
    rows: tuple[AlignedRow, ...]
    kind: str = PAIRWISE
    #: how it was made: mode, matrix, gap costs, ... (plain JSON values)
    parameters: dict[str, Any] = field(default_factory=dict)
    title: str = ""

    @property
    def width(self) -> int:
        return len(self.rows[0].positions) if self.rows else 0

    def proteins(self) -> set[str]:
        return {row.protein for row in self.rows}

    def to_json(self) -> dict[str, Any]:
        return {"title": self.title, "kind": self.kind, "parameters": dict(self.parameters),
                "rows": [row.to_json() for row in self.rows]}  # fmt: skip

    @classmethod
    def from_json(cls, saved: Any) -> Alignment | None:
        if not isinstance(saved, dict) or not isinstance(saved.get("rows"), list):
            return None
        rows = tuple(r for r in map(AlignedRow.from_json, saved["rows"]) if r is not None)
        if not rows or len({len(r.positions) for r in rows}) != 1:
            return None
        kind = saved.get("kind")
        parameters = saved.get("parameters")
        title = saved.get("title")
        return cls(
            rows,
            kind if kind in (PAIRWISE, MSA) else PAIRWISE,
            parameters if isinstance(parameters, dict) else {},
            title if isinstance(title, str) else "",
        )
