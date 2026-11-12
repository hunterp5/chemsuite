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

"""Residue templates with ideal geometry, from the wwPDB Chemical Component Dictionary.

A :class:`Template` is one CCD entry: every atom (hydrogens too) with its ideal coordinates
(``pdbx_model_Cartn_*_ideal``), formal charge and bonds. The CCD is public domain (CC0 1.0).

Atoms are built from a template by internal coordinates (:class:`Placement`): an atom sits at
a bond length from its *parent*, at an angle to a second atom and a dihedral to a third, all
measured on the ideal geometry. The second and third atoms are heavy atoms bonded to the
parent where it has two, so a CH₂ hydrogen or an aromatic hydrogen follows the real heavy
atoms rather than the ideal side-chain angles. A hydrogen on a parent with one heavy neighbour
is a *rotor* when it can turn (OH, SH, NH₃⁺, CH₃) and fixed when its parent is planar (an amide
or guanidinium NH₂); the dihedral of a rotor is then a free choice.

The standard amino acids, the nucleotides and the ACE/NME caps ship with the app
(``data/templates.json``, rebuilt by ``scripts/build_prepare_templates.py``); any other residue
comes from the CCD through :class:`~chemsuite.protein.structure.components.CcdCache`.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from functools import cache, cached_property
from pathlib import Path

import gemmi
import numpy as np

from chemsuite.core.jobs import CancelToken
from chemsuite.protein.structure.components import CcdCache

DATA = Path(__file__).with_name("data") / "templates.json"
HYDROGEN = 1
_ORDERS = {"SING": 1, "DOUB": 2, "TRIP": 3}
#: what ships with the app: the amino acids (standard and common modified ones), the caps and
#: the nucleotides
BUNDLED = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "MSE", "SEP", "TPO", "PTR",
    "HYP", "MLY", "M3L", "CSO", "CSD", "OCS", "CME", "KCX", "LLP", "ALY", "CGU", "FME",
    "ACE", "NH2", "NME", "A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI",
)  # fmt: skip


@dataclass(frozen=True, eq=False)
class Template:
    name: str
    atoms: tuple[str, ...]
    elements: tuple[int, ...]  # atomic numbers
    charges: tuple[int, ...]
    xyz: np.ndarray  # (n, 3) ideal coordinates
    bonds: tuple[tuple[int, int, int], ...]  # (i, j, order)
    leaving: frozenset[str] = frozenset()

    @cached_property
    def index(self) -> dict[str, int]:
        return {name: i for i, name in enumerate(self.atoms)}

    @cached_property
    def neighbours(self) -> tuple[tuple[int, ...], ...]:
        found: list[list[int]] = [[] for _ in self.atoms]
        for i, j, _ in self.bonds:
            found[i].append(j)
            found[j].append(i)
        return tuple(tuple(sorted(n)) for n in found)

    def is_hydrogen(self, i: int) -> bool:
        return self.elements[i] == HYDROGEN

    def heavy_neighbours(self, i: int) -> tuple[int, ...]:
        return tuple(j for j in self.neighbours[i] if not self.is_hydrogen(j))

    def hydrogens_on(self, name: str) -> tuple[str, ...]:
        i = self.index[name]
        return tuple(self.atoms[j] for j in self.neighbours[i] if self.is_hydrogen(j))

    @property
    def heavy(self) -> tuple[str, ...]:
        return tuple(a for i, a in enumerate(self.atoms) if not self.is_hydrogen(i))

    # variants -----------------------------------------------------------------------------
    def with_hydrogen(self, parent: str, name: str) -> Template:
        """One more hydrogen on *parent*, where a tetrahedral (sp³) atom has room for it: the
        CCD draws Lys and the N-terminus neutral, NH₂ not NH₃⁺."""
        p = self.index[parent]
        around = self.xyz[list(self.neighbours[p])] - self.xyz[p]
        direction = -(around / np.linalg.norm(around, axis=1)[:, None]).sum(axis=0)
        length = x_h_length(self.elements[p])
        if np.linalg.norm(direction) < 1e-3:
            raise ValueError(f"{self.name} {parent} has no room for another hydrogen")
        at = self.xyz[p] + length * direction / np.linalg.norm(direction)
        i = len(self.atoms)
        return replace(
            self,
            atoms=(*self.atoms, name),
            elements=(*self.elements, HYDROGEN),
            charges=(*self.charges, 0),
            xyz=np.vstack([self.xyz, at]),
            bonds=(*self.bonds, (p, i, 1)),
        )

    def without(self, names: Iterable[str]) -> Template:
        gone = set(names)
        keep = [i for i, a in enumerate(self.atoms) if a not in gone]
        new = {old: k for k, old in enumerate(keep)}
        return replace(
            self,
            atoms=tuple(self.atoms[i] for i in keep),
            elements=tuple(self.elements[i] for i in keep),
            charges=tuple(self.charges[i] for i in keep),
            xyz=self.xyz[keep],
            bonds=tuple((new[i], new[j], o) for i, j, o in self.bonds if i in new and j in new),
            leaving=self.leaving - gone,
        )

    # internal coordinates -----------------------------------------------------------------
    @cached_property
    def hydrogens(self) -> tuple[Placement, ...]:
        """How to place each hydrogen from the heavy atoms (those that cannot be placed, on a
        lone atom such as a water oxygen, are left out)."""
        found = []
        for h in range(len(self.atoms)):
            if not self.is_hydrogen(h) or not self.neighbours[h]:
                continue
            placement = self._placement(h, self.neighbours[h][0])
            if placement is not None:
                found.append(placement)
        return tuple(found)

    def placement(self, atom: str, parent: str, angle: str, dihedral: str) -> Placement:
        """Atom *atom* placed from three named atoms, measured on the ideal geometry."""
        i, p, a, d = (self.index[n] for n in (atom, parent, angle, dihedral))
        x = self.xyz
        return Placement(atom, parent, angle, dihedral, float(np.linalg.norm(x[i] - x[p])),
                         angle_of(x[a], x[p], x[i]), dihedral_of(x[d], x[a], x[p], x[i]),
                         rotor=False)  # fmt: skip

    def _placement(self, h: int, p: int) -> Placement | None:
        heavy = self.heavy_neighbours(p)
        if len(heavy) >= 2:
            a, d = heavy[0], heavy[1]
            rotor = False
        elif len(heavy) == 1:
            a = heavy[0]
            beyond = [n for n in self.heavy_neighbours(a) if n != p]
            if not beyond:
                return None
            d = beyond[0]
            symbol = self.elements[p]
            rotor = symbol in (8, 16, 34) or len(self.neighbours[p]) == 4
        else:
            return None
        placed = self.placement(self.atoms[h], self.atoms[p], self.atoms[a], self.atoms[d])
        return replace(placed, rotor=rotor)

    # storage ------------------------------------------------------------------------------
    def as_json(self) -> dict:
        symbols = [gemmi.Element(z).name for z in self.elements]
        return {
            "atoms": [[a, s, c, *(round(float(v), 3) for v in xyz)]
                      for a, s, c, xyz in zip(self.atoms, symbols, self.charges, self.xyz,
                                              strict=True)],
            "bonds": [[self.atoms[i], self.atoms[j], o] for i, j, o in self.bonds],
            "leaving": sorted(self.leaving),
        }  # fmt: skip

    @classmethod
    def from_json(cls, name: str, entry: Mapping) -> Template:
        rows = entry["atoms"]
        atoms = tuple(r[0] for r in rows)
        index = {a: i for i, a in enumerate(atoms)}
        return cls(
            name=name,
            atoms=atoms,
            elements=tuple(gemmi.Element(r[1]).atomic_number for r in rows),
            charges=tuple(int(r[2]) for r in rows),
            xyz=np.array([r[3:6] for r in rows], dtype=np.float64).reshape(-1, 3),
            bonds=tuple((index[a], index[b], int(o)) for a, b, o in entry["bonds"]),
            leaving=frozenset(entry.get("leaving", ())),
        )


@dataclass(frozen=True)
class Placement:
    """Where *atom* goes: *length* Å from *parent*, at *angle* (radians) to *parent*–*angle*
    and at *dihedral* about the *angle*–*parent* bond, measured from *dihedral_atom*."""

    atom: str
    parent: str
    angle_atom: str
    dihedral_atom: str
    length: float
    angle: float
    dihedral: float
    rotor: bool


def from_block(block: gemmi.cif.Block) -> Template | None:
    """The template a CCD entry defines (None when it has no coordinates)."""
    atoms = block.find("_chem_comp_atom.", [
        "comp_id", "atom_id", "type_symbol", "?charge", "?pdbx_leaving_atom_flag",
        "?pdbx_model_Cartn_x_ideal", "?pdbx_model_Cartn_y_ideal", "?pdbx_model_Cartn_z_ideal",
        "?model_Cartn_x", "?model_Cartn_y", "?model_Cartn_z",
    ])  # fmt: skip
    rows = list(atoms)
    if not rows:
        return None
    ideal = all(r.has(5) and _number(r.str(5)) is not None for r in rows)
    first = 5 if ideal else 8
    xyz = []
    for r in rows:
        point = [_number(r.str(first + k)) if r.has(first + k) else None for k in range(3)]
        if any(v is None for v in point):
            return None
        xyz.append(point)
    names = tuple(r.str(1) for r in rows)
    index = {a: i for i, a in enumerate(names)}
    bonds = []
    for r in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        a, b = r.str(0), r.str(1)
        if a in index and b in index:
            bonds.append((index[a], index[b], _ORDERS.get(r.str(2).upper(), 1)))
    return Template(
        name=rows[0].str(0).upper(),
        atoms=names,
        elements=tuple(gemmi.Element(r.str(2).capitalize()).atomic_number for r in rows),
        charges=tuple(int(_number(r.str(3)) or 0) if r.has(3) else 0 for r in rows),
        xyz=np.array(xyz, dtype=np.float64),
        bonds=tuple(bonds),
        leaving=frozenset(r.str(1) for r in rows if r.has(4) and r.str(4).upper() == "Y"),
    )


def _number(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


@cache
def bundled() -> dict[str, Template]:
    entries = json.loads(DATA.read_text(encoding="utf-8"))
    return {name: Template.from_json(name, entry) for name, entry in entries.items()}


class Library:
    """Templates by residue name: the bundled ones, then the CCD (when a cache is given)."""

    def __init__(self, ccd: CcdCache | None = None) -> None:
        self.ccd = ccd
        self._found: dict[str, Template | None] = {}

    def get(self, name: str, cancel: CancelToken | None = None) -> Template | None:
        name = name.upper()
        if name in bundled():
            return bundled()[name]
        if name not in self._found:
            block = self.ccd.block(name, cancel) if self.ccd is not None else None
            self._found[name] = from_block(block) if block is not None else None
        return self._found[name]


# geometry ---------------------------------------------------------------------------------
def x_h_length(element: int) -> float:
    return {6: 1.09, 7: 1.01, 8: 0.96, 16: 1.34}.get(element, 1.0)


def angle_of(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """The angle a–b–c in radians."""
    u, v = a - b, c - b
    cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.arccos(np.clip(cos, -1.0, 1.0)))


def dihedral_of(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    """The dihedral a–b–c–d in radians."""
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def place(d: np.ndarray, a: np.ndarray, p: np.ndarray, length, angle, dihedral) -> np.ndarray:
    """The atom at *length* from *p*, *angle* to *a*–*p* and *dihedral* d–a–p–atom (NeRF).
    Every argument may carry leading batch dimensions."""
    bc = p - a
    bc = bc / np.linalg.norm(bc, axis=-1, keepdims=True)
    n = np.cross(a - d, bc)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    m = np.cross(n, bc)
    length, angle, dihedral = (np.asarray(v, dtype=np.float64)[..., None]
                               for v in (length, angle, dihedral))  # fmt: skip
    local = (
        -length * np.cos(angle) * bc
        + length * np.sin(angle) * np.cos(dihedral) * m
        + length * np.sin(angle) * np.sin(dihedral) * n
    )
    return p + local
