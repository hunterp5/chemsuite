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

"""Hydrogens on a structure's heavy atoms, every group in its textbook state.

Polymer residues are built from their templates (:mod:`.templates`):

- Asp and Glu charged, Lys and Arg charged, His neutral with H on NE2 (HIE), Cys, Ser, Thr
  and Tyr neutral; nucleotides with a charged phosphate;
- the N-terminus NH₃⁺ and the C-terminus a carboxylate;
- an atom bonded to another residue (a peptide or phosphodiester link, a disulfide, a glycan,
  a covalent ligand) carries one hydrogen fewer for each such bond, and the hydrogens it keeps
  are placed from the real neighbours on both sides (a backbone NH lies in the peptide plane);
- an N, O or S within :data:`METAL_CONTACT` of a transition metal gives up one hydrogen (a Zn
  thiolate); a His on a metal keeps its hydrogen on the other ring N.

Ligands, additives and glycans get hydrogens from RDKit (``AddHs(addCoords=True)``) on the
molecule :func:`~chemsuite.protein.chem.ligands.to_mol` builds, so their state is the one their
component defines. Waters and ions get none. A residue that already has hydrogens is left as it
is. Flips, rotor angles and other states are chosen later (``sites``/``solve``); this module
records the rotors they turn (:class:`Rotor`).

Hydrogens are placed on the atoms :meth:`~.model.Structure.active` picks (one location each).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

from ..chem import ligands
from ..structure import edit
from ..structure.components import CcdCache
from ..structure.model import AMINO, ION, NUCLEIC, WATER, Structure
from . import templates
from .templates import HYDROGEN, Library, Template

#: metals that take a proton's place on a coordinating N, O or S
SOFT_METALS = frozenset([23, 24, 25, 26, 27, 28, 29, 30, 42, 47, 48, 74, 78, 79, 80])
METAL_CONTACT = 2.9  # Å
#: hydrogens a residue loses in its textbook state at neutral pH
TEXTBOOK = {"ASP": ("HD2",), "GLU": ("HE2",)}
PHOSPHATE = ("HOP2", "HOP3")
POLAR = frozenset([7, 8, 16])
KEPT = "Kept the hydrogens it has"  # the note on a residue that came with hydrogens
_TETRAHEDRAL = np.radians(109.47)
_TRIGONAL = np.radians(120.0)


@dataclass(frozen=True)
class Rotor:
    """Hydrogens that turn together about the bond *axis*–*parent* (all uids)."""

    parent: int
    axis: int
    hydrogens: tuple[int, ...]


@dataclass(frozen=True)
class Note:
    residue: int
    text: str


@dataclass(frozen=True)
class Protonated:
    structure: Structure
    #: the state picked for a residue where there was a choice: HID, HIE, CYM (thiolate)…
    states: Mapping[int, str]
    rotors: tuple[Rotor, ...]
    notes: tuple[Note, ...]


def add_hydrogens(
    s: Structure,
    library: Library | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
) -> Protonated:
    """*s* with hydrogens on every residue that has none (uids of the old atoms kept)."""
    library = library or Library(ccd)
    plan = _Plan(s, library, ccd, cancel)
    for residue in range(s.n_residues):
        check(cancel)
        plan.residue(residue)
    return plan.result()


# -- the plan --------------------------------------------------------------------------------------
@dataclass
class _New:
    residue: int
    name: str
    parent: int
    xyz: np.ndarray  # (k, 3)
    rotor_axis: int = -1  # the heavy atom it turns about, -1 when fixed


@dataclass
class _Plan:
    s: Structure
    library: Library
    ccd: CcdCache | None
    cancel: CancelToken | None
    new: list[_New] = field(default_factory=list)
    charges: dict[int, int] = field(default_factory=dict)
    states: dict[int, str] = field(default_factory=dict)
    notes: list[Note] = field(default_factory=list)
    active: np.ndarray = field(init=False)
    coords: np.ndarray = field(init=False)
    neighbours: list[list[int]] = field(init=False)
    metals: np.ndarray = field(init=False)
    metal_tree: cKDTree | None = field(init=False)

    def __post_init__(self) -> None:
        s = self.s
        self.active = s.active()
        self.coords = s.coords.astype(np.float64)
        heavy = self.active & (s.element != HYDROGEN)
        self.neighbours: list[list[int]] = [[] for _ in range(s.n_atoms)]
        for a, b in s.bonds.tolist():
            if heavy[a] and heavy[b]:
                self.neighbours[a].append(b)
                self.neighbours[b].append(a)
        metals = np.flatnonzero(self.active & np.isin(s.element, list(SOFT_METALS))
                                & (s.atom_kind == ION))  # fmt: skip
        self.metals = metals
        self.metal_tree = cKDTree(self.coords[0][metals]) if len(metals) else None

    # -- residues ---------------------------------------------------------------------------
    def residue(self, r: int) -> None:
        s = self.s
        kind = int(s.res_kind[r])
        if kind in (WATER, ION):
            return
        span = s.atoms_of(r)
        atoms = [i for i in range(span.start, span.stop) if self.active[i]]
        if not atoms:
            return
        if (s.element[atoms] == HYDROGEN).any():
            self.notes.append(Note(r, KEPT))
            return
        if kind in (AMINO, NUCLEIC):
            template = self.library.get(str(s.res_name[r]), self.cancel)
            if template is not None:
                self.from_template(r, atoms, template)
                return
            self.notes.append(Note(r, "No template: hydrogens from RDKit"))
        self.from_rdkit(r)

    def external(self, i: int) -> list[int]:
        r = self.s.residue_of[i]
        return [j for j in self.neighbours[i] if self.s.residue_of[j] != r]

    def metals_near(self, i: int) -> list[int]:
        if self.metal_tree is None or int(self.s.element[i]) not in POLAR:
            return []
        found = self.metal_tree.query_ball_point(self.coords[0][i], METAL_CONTACT)
        return [int(self.metals[k]) for k in sorted(found)]

    # -- polymer residues from templates ----------------------------------------------------
    def from_template(self, r: int, atoms: list[int], t: Template) -> None:
        s = self.s
        here = {str(s.name[i]): i for i in atoms}
        unknown = sorted(set(here) - set(t.atoms))
        if unknown:
            self.notes.append(Note(r, f"Atoms not in the {t.name} template: {', '.join(unknown)}"))
        delta: dict[str, int] = {}
        removed: set[str] = set()
        kind = int(s.res_kind[r])
        if kind == AMINO:
            removed |= {"HXT"} & set(t.atoms)  # a carboxylate C-terminus
            n = here.get("N")
            if n is not None and not self.external(n) and "N" in t.index and t.hydrogens_on("N"):
                free = next(h for h in ("H", "H2", "H3", "HN3") if h not in t.index)
                try:
                    t = t.with_hydrogen("N", free)  # NH3+
                    delta["N"] = 1
                except ValueError:
                    pass  # a planar N (N-formylmethionine's amide) keeps its hydrogen
        elif kind == NUCLEIC:
            removed |= set(PHOSPHATE) & set(t.atoms)
        removed |= set(TEXTBOOK.get(t.name, ())) & set(t.atoms)
        ring: dict[str, int] = {}
        if t.name == "HIS":
            gone, ring = self.histidine(r, here)
            removed |= gone & set(t.atoms)
        for h in removed:
            parent = t.atoms[t.neighbours[t.index[h]][0]]
            delta[parent] = delta.get(parent, 0) - 1
        t = t.without(removed)
        for name, i in here.items():
            if name in t.index:
                self.charges[i] = ring.get(name, t.charges[t.index[name]] + delta.get(name, 0))
        by_name = {p.atom: p for p in t.hydrogens}
        for name, i in here.items():
            if name in t.index and t.hydrogens_on(name):
                self.place_on(r, i, t, name, here, by_name)

    def histidine(self, r: int, here: dict[str, int]) -> tuple[set[str], dict[str, int]]:
        """The tautomer, and the ring N charges: H on NE2 (HIE) unless NE2 binds a metal or
        another residue (then HID). A ring N bonded to another residue next to an NH is
        charged (an alkylated imidazolium)."""
        ring = {"ND1": 0, "NE2": 0}
        if "ND1" not in here or "NE2" not in here:
            return {"HD1"}, ring
        nd1, ne2 = here["ND1"], here["NE2"]
        taken_d = bool(self.metals_near(nd1) or self.external(nd1))
        taken_e = bool(self.metals_near(ne2) or self.external(ne2))
        if taken_d and taken_e:
            state, gone = "HIS-", {"HD1", "HE2"}
            ring["ND1"] = -1
        elif taken_e:
            state, gone = "HID", {"HE2"}
            ring["NE2"] = 1 if self.external(ne2) else 0
        else:
            state, gone = "HIE", {"HD1"}
            ring["ND1"] = 1 if self.external(nd1) else 0
        self.states[r] = state
        return gone, ring

    def place_on(
        self, r: int, i: int, t: Template, name: str, here: dict[str, int], by_name: dict
    ) -> None:
        s = self.s
        hs = t.hydrogens_on(name)
        partners = self.external(i)
        metals = [] if t.name == "HIS" else self.metals_near(i)
        lost = min(len(hs), len(metals))
        if lost:
            self.charges[i] = self.charges.get(i, 0) - lost
            if t.name == "CYS":
                self.states[r] = "CYM"
        keep = len(hs) - len(partners) - lost
        if keep <= 0:
            return
        x = self.coords
        placed = [by_name.get(h) for h in hs]
        fixed = all(
            p is not None and all(a in here for a in (p.parent, p.angle_atom, p.dihedral_atom))
            for p in placed
        )
        if fixed and not partners:
            spots = [
                (h, templates.place(x[:, here[p.dihedral_atom]], x[:, here[p.angle_atom]],
                                    x[:, i], p.length, p.angle, p.dihedral),
                 here[p.angle_atom] if p.rotor else -1)
                for h, p in zip(hs, placed, strict=True)
            ]  # fmt: skip
            avoid = [x[0, j] for j in metals]
            for h, at, axis in _nearest_dropped(spots, avoid, len(hs) - keep):
                self.new.append(_New(r, h, i, at, axis))
            return
        # placed from the real neighbours, those in other residues too
        names = sorted(hs, key=lambda h: h in t.leaving)[:keep]
        heavy = [j for j in self.neighbours[i] if j != i]
        planar = _planar(t, name) or (int(s.element[i]) == 7
                                      and any(str(s.name[j]) == "C" for j in partners))  # fmt: skip
        positions = _geometric(self, i, heavy, keep, planar)
        if positions is None:
            self.notes.append(Note(r, f"No room to place hydrogens on {name}"))
            return
        rotor = (
            len(heavy) == 1
            and not planar
            and (int(s.element[i]) in (8, 16, 34) or len(heavy) + keep == 4)
        )
        for h, at in zip(names, positions, strict=True):
            self.new.append(_New(r, h, i, at, heavy[0] if rotor else -1))

    # -- ligands from RDKit -----------------------------------------------------------------
    def from_rdkit(self, r: int) -> None:
        s = self.s
        try:
            first = ligands.to_mol(s, r, active=self.active, ccd=self.ccd, cancel=self.cancel)
        except (ValueError, RuntimeError) as error:
            self.notes.append(Note(r, f"No hydrogens: {error}"))
            return
        if first.guessed:
            self.notes.append(Note(r, "Bond orders guessed from the geometry"))
        mols = [first.mol] + [
            ligands.to_mol(s, r, active=self.active, model=k, ccd=self.ccd).mol
            for k in range(1, s.n_models)
        ]
        full = [Chem.AddHs(m, addCoords=True) for m in mols]
        template = self.library.get(str(s.res_name[r]), self.cancel)
        taken = {str(n) for n in s.name[s.atoms_of(r)]}
        lookup = s.index_of_uid
        mol = full[0]
        positions = [m.GetConformer().GetPositions() for m in full]
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == HYDROGEN:
                continue
            i = lookup[atom.GetIntProp("uid")]
            self.charges[i] = atom.GetFormalCharge()
            hs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == HYDROGEN]
            if not hs:
                continue
            partners = self.external(i)
            spots = [(h, np.array([at[h] for at in positions]), -1) for h in hs]
            kept = _nearest_dropped(spots, [self.coords[0][j] for j in partners], len(partners))
            heavy = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != HYDROGEN]
            rotor = (len(heavy) + len(partners) == 1
                     and atom.GetHybridization() == Chem.HybridizationType.SP3)  # fmt: skip
            axis = -1
            if rotor:
                axis = lookup[heavy[0].GetIntProp("uid")] if heavy else partners[0]
            names = ligand_names(template, str(s.name[i]), len(kept), taken)
            for name, (_h, at, _a) in zip(names, kept, strict=True):
                self.new.append(_New(r, name, i, at, axis))

    # -- the result -------------------------------------------------------------------------
    def result(self) -> Protonated:
        s = self.s
        charge = s.charge.copy()
        for i, c in self.charges.items():
            charge[i] = c
        s = replace(s, charge=charge)
        n, m = s.n_atoms, len(self.new)
        if not m:
            return Protonated(s, dict(self.states), (), tuple(self.notes))
        grown = edit.add_atoms(
            s,
            residues=np.array([h.residue for h in self.new]),
            names=np.array([h.name for h in self.new]),
            elements=np.full(m, HYDROGEN),
            xyz=np.stack([h.xyz for h in self.new], axis=1),
            bonds=np.array([[h.parent, n + j] for j, h in enumerate(self.new)]),
        )
        groups: dict[tuple[int, int], list[int]] = {}
        for j, h in enumerate(self.new):
            if h.rotor_axis >= 0:
                groups.setdefault((h.parent, h.rotor_axis), []).append(s.next_uid + j)
        rotors = tuple(Rotor(int(s.uid[p]), int(s.uid[a]), tuple(u))
                       for (p, a), u in groups.items())  # fmt: skip
        return Protonated(grown, dict(self.states), rotors, tuple(self.notes))


# -- geometry --------------------------------------------------------------------------------------
def _nearest_dropped(spots: list, avoid: list[np.ndarray], drop: int) -> list:
    """*spots* less the *drop* nearest to any point of *avoid* (judged in the first model)."""
    if drop <= 0 or not avoid:
        return spots[: len(spots) - max(drop, 0)]
    points = np.array(avoid)
    distance = [float(np.min(np.linalg.norm(points - at[0], axis=1))) for _h, at, _a in spots]
    gone = set(np.argsort(distance)[:drop].tolist())
    return [spot for k, spot in enumerate(spots) if k not in gone]


def _planar(t: Template, name: str) -> bool:
    """Whether *name* is sp²: it has a multiple bond, or it is an N conjugated with one."""
    i = t.index[name]

    def multiple(k: int) -> bool:
        return any(o >= 2 and k in (a, b) for a, b, o in t.bonds)

    if multiple(i):
        return True
    return t.elements[i] == 7 and any(multiple(k) for k in t.heavy_neighbours(i))


def _geometric(plan: _Plan, i: int, heavy: list[int], count: int, planar: bool) -> list | None:
    """*count* hydrogens on atom *i* from its heavy neighbours, in each model (VSEPR)."""
    x = plan.coords
    length = templates.x_h_length(int(plan.s.element[i]))
    p = x[:, i]
    units = [_unit(x[:, j] - p) for j in heavy]
    if len(heavy) >= 2:
        bisector = _unit(-sum(units))
        if len(heavy) >= 3 or planar:
            return [p + length * bisector][:count]
        normal = _unit(np.cross(units[0], units[1]))
        half = _TETRAHEDRAL / 2
        return [p + length * (np.cos(half) * bisector + sign * np.sin(half) * normal)
                for sign in (1, -1)][:count]  # fmt: skip
    if len(heavy) == 1:
        a = heavy[0]
        beyond = [j for j in plan.neighbours[a] if j != i]
        ref = x[:, beyond[0]] if beyond else x[:, a] + _perpendicular(units[0])
        if planar:
            spots = [templates.place(ref, x[:, a], p, length, _TRIGONAL, d) for d in (np.pi, 0.0)]
        else:
            spots = [templates.place(ref, x[:, a], p, length, _TETRAHEDRAL, d)
                     for d in (np.pi, np.pi / 3, -np.pi / 3)]  # fmt: skip
        return spots[:count]
    return None


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _perpendicular(u: np.ndarray) -> np.ndarray:
    axis = np.array([1.0, 0.0, 0.0]) if abs(u[0, 0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    return _unit(np.cross(u, axis))


def ligand_names(template: Template | None, parent: str, count: int, taken: set[str]) -> list:
    """The component's names for the hydrogens on *parent*, else H1, H2… not yet taken."""
    own = () if template is None or parent not in template.index else template.hydrogens_on(parent)
    names = [h for h in own if h not in taken][:count]
    k = 1
    while len(names) < count:
        if f"H{k}" not in taken and f"H{k}" not in names:
            names.append(f"H{k}")
        k += 1
    taken.update(names)
    return names
