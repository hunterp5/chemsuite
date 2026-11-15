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
"""The poses a docking run produced, kept with the document.

A docking run is work done in the Protein Viewer, so its poses stay there: the Poses panel
lists them beside the pocket they were docked into, and they reach the table only when the
chemist asks for them (Poses to Table). Before this they went straight into the table, which
put a pocket's results in another window and filled the table with rows nobody had asked for.

Each :class:`Pose` holds its own 3D molblock, so a pose can be drawn, scored again or copied
to the table without the docking run. :class:`DockedPoses` is one run against one protein. The
saved form is plain JSON (as :mod:`.reports` is), under ``poses`` of the protein workspace's
state; a run whose protein has gone is dropped on the way out.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

#: what a run was: a dock, or a rescore of poses that already existed
DOCKED, RESCORED = "docked", "rescored"


def _number(value: Any, fallback: float = 0.0) -> float:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    return fallback


@dataclass(frozen=True)
class Pose:
    """One pose of one ligand: its scores and the pose itself."""

    ligand: str  # the row it was docked from, as the table names it
    rank: int  # 1 is Vina's best
    affinity: float  # kcal/mol
    rmsd_lb: float
    rmsd_ub: float
    molblock: str  # the pose in 3D

    def to_json(self) -> list[Any]:
        return [self.ligand, self.rank, self.affinity, self.rmsd_lb, self.rmsd_ub, self.molblock]

    @classmethod
    def from_json(cls, saved: Any) -> Pose | None:
        if not isinstance(saved, list) or len(saved) != 6:
            return None
        ligand, rank, affinity, lb, ub, molblock = saved
        if not isinstance(ligand, str) or not isinstance(molblock, str):
            return None
        if not isinstance(rank, int) or isinstance(rank, bool):
            return None
        return cls(ligand, rank, _number(affinity), _number(lb), _number(ub), molblock)


@dataclass(frozen=True)
class DockedPoses:
    """One docking run against one protein."""

    protein: str
    poses: tuple[Pose, ...]
    kind: str = DOCKED  # DOCKED or RESCORED
    #: the ligands that could not be docked, and the rows that held no structure
    failed: int = 0
    skipped: int = 0

    @property
    def ligands(self) -> tuple[str, ...]:
        """The ligands with a pose, in the order they were docked."""
        return tuple(dict.fromkeys(pose.ligand for pose in self.poses))

    @property
    def best(self) -> float | None:
        """The best (lowest) affinity of the run, if it produced anything."""
        return min((pose.affinity for pose in self.poses), default=None)

    def to_json(self) -> dict[str, Any]:
        return {
            "protein": self.protein,
            "kind": self.kind,
            "failed": self.failed,
            "skipped": self.skipped,
            "poses": [pose.to_json() for pose in self.poses],
        }

    @classmethod
    def from_json(cls, saved: Any) -> DockedPoses | None:
        """The run *saved* held, or None when it is not one."""
        if not isinstance(saved, dict) or not isinstance(saved.get("protein"), str):
            return None
        poses = saved.get("poses")
        if not isinstance(poses, list):
            return None
        made = [Pose.from_json(pose) for pose in poses]
        if any(pose is None for pose in made):
            return None
        kind = saved.get("kind")
        counts = [saved.get("failed"), saved.get("skipped")]
        failed, skipped = (int(n) if isinstance(n, int) and not isinstance(n, bool) else 0
                           for n in counts)  # fmt: skip
        return cls(
            saved["protein"],
            tuple(pose for pose in made if pose is not None),
            kind if kind in (DOCKED, RESCORED) else DOCKED,
            failed,
            skipped,
        )
