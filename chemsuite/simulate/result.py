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

"""What a minimisation gives back, as the worker reports it.

The new coordinates come by uid, so they land on the right atoms of the structure whatever
order it holds them in (:meth:`MinimiseResult.coordinates`). Alongside: the energies, how far
the heavy atoms moved (overall, in the pocket, the ligand), the worst contact before and after,
the largest force left, the residues that moved most, the charges and force-field versions
used, the platform, and how long each stage took. :meth:`MinimiseResult.warnings` says what
deserves a ⚠: heavy atoms further than the limit from where they started, or a clash left.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from chemsuite.macromolecule.structure.model import Structure

from .spec import SpecError

#: closer than this (Å) between heavy atoms of different residues, not bonded, is a clash
CLASH = 2.2


@dataclass(frozen=True)
class Contact:
    """The closest approach of two atoms (their labels) of different residues."""

    a: str
    b: str
    distance: float  # Å


@dataclass(frozen=True)
class MinimiseResult:
    digest: str  # the SystemSpec's, which these coordinates belong to
    uid: tuple[int, ...]
    xyz: tuple[tuple[float, float, float], ...]  # Å, by uid
    energy_before: float  # kJ/mol
    energy_after: float
    rmsd: float  # Å, every heavy atom
    rmsd_pocket: float | None = None  # heavy atoms within 5 Å of a ligand
    rmsd_ligand: float | None = None
    contact_before: Contact | None = None
    contact_after: Contact | None = None
    max_force: float = 0.0  # kJ/mol/nm, the largest force left
    max_force_atom: str = ""
    moved_most: tuple[tuple[str, float], ...] = ()  # residues by label, largest first, Å
    charges: tuple[tuple[str, str, tuple[float, ...]], ...] = ()  # (residue, model, charges)
    versions: tuple[tuple[str, str], ...] = ()  # (package, version)
    platform: str = ""
    seconds: tuple[tuple[str, float], ...] = ()  # (stage, seconds)

    def __post_init__(self) -> None:
        if len(self.uid) != len(self.xyz):
            raise SpecError("One position per atom")

    def coordinates(self, s: Structure, model: int = 0) -> np.ndarray:
        """*s*'s *model* (n, 3) with the minimised positions put in, by uid; atoms the
        minimisation did not hold keep theirs. Raises :class:`~.spec.SpecError` when an atom
        it moved is no longer in *s*."""
        found = s.indices(np.array(self.uid, dtype=np.int64))
        if (found < 0).any():
            raise SpecError("Atoms of the minimisation are no longer in the structure")
        out = s.xyz(model).astype(np.float32).copy()
        out[found] = np.array(self.xyz, dtype=np.float32).reshape(-1, 3)
        return out

    def warnings(self, rmsd_limit: float) -> list[str]:
        """What deserves a ⚠, worded for the report (empty when nothing does)."""
        out = []
        if self.rmsd > rmsd_limit:
            out.append(f"Heavy atoms moved {self.rmsd:.2f} Å, more than {rmsd_limit:.2f} Å")
        if self.contact_after is not None and self.contact_after.distance < CLASH:
            c = self.contact_after
            out.append(f"A clash is left: {c.a} – {c.b}, {c.distance:.2f} Å")
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            "digest": self.digest,
            "uid": list(self.uid),
            "xyz": [round(v, 4) for p in self.xyz for v in p],
            "energy_before": self.energy_before,
            "energy_after": self.energy_after,
            "rmsd": self.rmsd,
            "rmsd_pocket": self.rmsd_pocket,
            "rmsd_ligand": self.rmsd_ligand,
            "contact_before": _contact_json(self.contact_before),
            "contact_after": _contact_json(self.contact_after),
            "max_force": self.max_force,
            "max_force_atom": self.max_force_atom,
            "moved_most": [list(m) for m in self.moved_most],
            "charges": [[label, model, list(q)] for label, model, q in self.charges],
            "versions": dict(self.versions),
            "platform": self.platform,
            "seconds": [list(s) for s in self.seconds],
        }

    @classmethod
    def from_json(cls, saved: Any) -> MinimiseResult:
        """The result *saved* holds (from the worker or the document); :class:`SpecError`
        when it is not one."""
        if not isinstance(saved, dict):
            raise SpecError("A minimisation result is a JSON object")
        uid = saved.get("uid")
        flat = saved.get("xyz")
        if (not isinstance(uid, list) or not all(_int(u) for u in uid)
                or not isinstance(flat, list) or not all(_number(v) for v in flat)
                or len(flat) != 3 * len(uid)):  # fmt: skip
            raise SpecError("A result has one uid and three coordinates per atom")
        digest = saved.get("digest")
        if not isinstance(digest, str):
            raise SpecError("A result names the system it belongs to")
        numbers = {k: saved.get(k) for k in ("energy_before", "energy_after", "rmsd")}
        if not all(_number(v) for v in numbers.values()):
            raise SpecError("A result has its energies and RMSD")
        xyz = tuple((float(flat[k]), float(flat[k + 1]), float(flat[k + 2]))
                    for k in range(0, len(flat), 3))  # fmt: skip
        return cls(
            digest=digest,
            uid=tuple(int(u) for u in uid),
            xyz=xyz,
            energy_before=float(numbers["energy_before"]),
            energy_after=float(numbers["energy_after"]),
            rmsd=float(numbers["rmsd"]),
            rmsd_pocket=_optional(saved.get("rmsd_pocket")),
            rmsd_ligand=_optional(saved.get("rmsd_ligand")),
            contact_before=_contact(saved.get("contact_before")),
            contact_after=_contact(saved.get("contact_after")),
            max_force=_optional(saved.get("max_force")) or 0.0,
            max_force_atom=_text(saved.get("max_force_atom")),
            moved_most=tuple((str(m[0]), float(m[1])) for m in _pairs(saved.get("moved_most"))),
            charges=_charges(saved.get("charges")),
            versions=tuple(sorted((str(k), str(v)) for k, v in _mapping(saved.get("versions")))),
            platform=_text(saved.get("platform")),
            seconds=tuple((str(m[0]), float(m[1])) for m in _pairs(saved.get("seconds"))),
        )


def _int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _optional(value: Any) -> float | None:
    return float(value) if _number(value) else None


def _text(value: Any) -> str:
    return value if isinstance(value, str) else ""


def _pairs(value: Any) -> list[list]:
    if not isinstance(value, list):
        return []
    return [
        v
        for v in value
        if isinstance(v, list) and len(v) == 2 and isinstance(v[0], str) and _number(v[1])
    ]


def _mapping(value: Any) -> list[tuple[Any, Any]]:
    return list(value.items()) if isinstance(value, dict) else []


def _contact_json(contact: Contact | None) -> dict[str, Any] | None:
    if contact is None:
        return None
    return {"a": contact.a, "b": contact.b, "distance": contact.distance}


def _contact(value: Any) -> Contact | None:
    if (
        not isinstance(value, dict)
        or not isinstance(value.get("a"), str)
        or not isinstance(value.get("b"), str)
        or not _number(value.get("distance"))
    ):
        return None
    return Contact(value["a"], value["b"], float(value["distance"]))


def _charges(value: Any) -> tuple[tuple[str, str, tuple[float, ...]], ...]:
    if not isinstance(value, list):
        return ()
    out = []
    for v in value:
        if (isinstance(v, list) and len(v) == 3 and isinstance(v[0], str)
                and isinstance(v[1], str) and isinstance(v[2], list)
                and all(_number(q) for q in v[2])):  # fmt: skip
            out.append((v[0], v[1], tuple(float(q) for q in v[2])))
    return tuple(out)
