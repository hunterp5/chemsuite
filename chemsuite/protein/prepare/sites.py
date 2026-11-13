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

"""Sites: the groups whose hydrogens (or whose orientation) are a choice, each with its states.

- **Asn and Gln**: the amide as built, or flipped (O and N swapped). X-ray density cannot tell
  O from N, so either may be right.
- **His**: HID, HIE or HIP, each with the ring as built or flipped (ND1↔CD2, CE1↔NE2). A His
  whose ring N is bonded to another residue keeps what :mod:`.hydrogens` gave it; one on a
  metal is a site like the others, which the energy's metal term orients.
- **Rotors**: the polar hydrogens that turn (OH, SH, NH₃⁺, on the protein and the ligands),
  in 15° steps; a phenol or acid OH stays in the plane of its neighbour.
- **Titration** (``titrate=True``): Asp and Glu also neutral, with H on either O (syn or
  anti); Lys
  and the N-terminus also neutral NH₂ (its lone pair in 15° steps); free Cys also S⁻; Tyr
  also O⁻.

Every state lists the atoms that move with it (heavy atoms and hydrogens, in every model), so
the energy can score it and the solver can write the choice back. Each carries its protons
(relative to the site's least protonated form) and its intrinsic free energy in water
(:data:`MICRO_PKA`, kcal/mol, the pH part left to the solver). The flips follow Reduce
(Word et al., *J. Mol. Biol.* 285, 1735–1747, 1999); the model pKa values are PROPKA 3's
(Olsson et al., *J. Chem. Theory Comput.* 7, 525–537, 2011).
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

import numpy as np

from ..structure.model import AMINO, AROMATIC, Structure
from . import hydrogens, templates
from .hydrogens import HYDROGEN, Protonated

AMIDES = {"ASN": ("CG", "OD1", "ND2", "CB"), "GLN": ("CD", "OE1", "NE2", "CG")}
ACIDS = {"ASP": ("CG", "OD1", "OD2", "HD1", "HD2"), "GLU": ("CD", "OE1", "OE2", "HE1", "HE2")}
RING = ("ND1", "CD2", "CE1", "NE2")
FLIPPED_RING = ("CD2", "ND1", "NE2", "CE1")  # where each ring atom goes when the ring flips
HIS_HYDROGENS = {"HID": ("HD1",), "HIE": ("HE2",), "HIP": ("HD1", "HE2")}
STEP = np.radians(15.0)
POLAR = frozenset([7, 8, 16])
EXCLUDED_BONDS = 3  # atoms this many bonds or fewer from a site's anchors are not scored
RT_LN10 = 1.364  # kcal/mol per pH unit at 298 K
KT = RT_LN10 / np.log(10)  # kcal/mol at 298 K
#: model pKa values (PROPKA 3)
MODEL_PKA = {"ASP": 3.8, "GLU": 4.5, "HIS": 6.5, "CYS": 9.0, "TYR": 10.0, "LYS": 10.5, "NTERM": 8.0}
#: His in water is about 4 parts HIE to 1 part HID
HIE_PER_HID = 4.0
#: (residue, parent atom) of the rotors that titrate, by residue name; NTERM is any N-terminus
TITRATING = {("LYS", "NZ"), ("CYS", "SG"), ("TYR", "OH")}
OH_LENGTH, SYN_ANGLE = 0.97, np.radians(107.0)
NH_LENGTH, TETRAHEDRAL = 1.01, np.radians(109.47)


@dataclass(frozen=True)
class State:
    label: str
    names: tuple[str, ...]  # the atoms of the state
    parents: tuple[str, ...]  # each hydrogen's heavy atom ("" for a heavy atom)
    elements: np.ndarray  # (m,)
    xyz: np.ndarray  # (k, m, 3): every model
    donor: np.ndarray  # (m, 3): the heavy atom a polar hydrogen sits on, NaN otherwise
    flipped: bool = False
    anti: bool = False  # an acid's OH pointing away from its other O
    tautomer: str = ""  # HID, HIE or HIP for a His
    form: str = ""  # the protonation form: ASP/ASH, GLU/GLH, LYS/LYN, CYS/CYM, TYR/TYM, HI?
    protons: int = 0  # relative to the site's least protonated form
    intrinsic: float = 0.0  # kcal/mol in water, less the pH term (RT ln10 · protons · pH)
    formal: tuple[int, ...] = ()  # formal charge of each atom (written back)
    partial: tuple[float, ...] = ()  # charge of each atom for the Coulomb term

    @property
    def polar_h(self) -> np.ndarray:
        return ~np.isnan(self.donor[:, 0])

    @property
    def acceptor(self) -> np.ndarray:
        """O, a thiolate S, and N with no hydrogen in this state."""
        has_h = {p for p in self.parents if p}
        charge = self.formal or (0,) * len(self.names)
        atoms = zip(self.names, self.elements, charge, strict=True)
        return np.array([e == 8 or (e == 7 and n not in has_h) or (e == 16 and q < 0)
                         for n, e, q in atoms])  # fmt: skip

    @property
    def lone(self) -> np.ndarray:
        """Acceptors with no hydrogen of their own."""
        has_h = {p for p in self.parents if p}
        return self.acceptor & np.array([n not in has_h for n in self.names])

    @property
    def charges(self) -> np.ndarray:
        return np.array(self.partial or (0.0,) * len(self.names))


@dataclass(frozen=True)
class Site:
    residue: int
    kind: str  # AMIDE, HIS, ACID or ROTOR
    moving: np.ndarray  # structure atoms the states replace
    anchors: np.ndarray  # structure atoms the moving ones hang from
    states: tuple[State, ...]
    excluded: frozenset[int] = field(default=frozenset())  # atoms not scored against it

    @property
    def titrates(self) -> bool:
        return len({st.protons for st in self.states}) > 1


AMIDE, HIS, ACID, ROTOR = "amide", "his", "acid", "rotor"


def find_sites(p: Protonated, titrate: bool = False) -> list[Site]:
    """Every site of the protonated structure *p*; with *titrate*, every titratable group
    also takes its other protonation form."""
    s = p.structure
    graph = _graph(s)
    found: list[Site] = []
    for r in np.flatnonzero(s.res_kind == AMINO).tolist():
        name = str(s.res_name[r])
        here = _named(s, r)
        site = None
        if name in AMIDES:
            site = _amide(s, r, here, graph)
        elif name == "HIS":
            site = _histidine(s, r, here, graph)
        elif name in ACIDS and titrate:
            site = _acid(s, r, here, graph)
        if site is not None:
            found.append(site)
    for rotor in p.rotors:
        site = _rotor(s, rotor, graph, titrate)
        if site is not None:
            found.append(site)
    return found


# -- amides ----------------------------------------------------------------------------------------
def _amide(s: Structure, r: int, here: dict[str, int], graph: list[list[int]]) -> Site | None:
    carbon, oxygen, nitrogen, before = AMIDES[str(s.res_name[r])]
    if not all(a in here for a in (carbon, oxygen, nitrogen, before)):
        return None
    if any(s.residue_of[j] != r for j in graph[here[nitrogen]] if s.element[j] != HYDROGEN):
        return None  # a glycan on the N
    t = templates.bundled()[str(s.res_name[r])]
    placements = [p for p in t.hydrogens if p.parent == nitrogen]
    x = s.coords.astype(np.float64)
    states = []
    for flipped in (False, True):
        at = {a: x[:, here[a]] for a in (carbon, before)}
        at[oxygen], at[nitrogen] = (
            (x[:, here[nitrogen]], x[:, here[oxygen]])
            if flipped
            else (x[:, here[oxygen]], x[:, here[nitrogen]])
        )
        names = [oxygen, nitrogen]
        parents = ["", ""]
        xyz = [at[oxygen], at[nitrogen]]
        for pl in placements:
            names.append(pl.atom)
            parents.append(nitrogen)
            xyz.append(_placed(pl, at))
        states.append(_state("flipped" if flipped else "as built", names, parents, xyz,
                             flipped=flipped))  # fmt: skip
    moving = _with_hydrogens(s, graph, [here[oxygen], here[nitrogen]])
    return _site(s, r, AMIDE, moving, [here[carbon]], states, graph)


# -- histidine -------------------------------------------------------------------------------------
def _histidine(s: Structure, r: int, here: dict[str, int], graph: list[list[int]]) -> Site | None:
    if not all(a in here for a in ("CB", "CG", *RING)):
        return None
    ring_n = (here["ND1"], here["NE2"])
    if any(s.residue_of[j] != r for n in ring_n for j in graph[n] if s.element[j] != HYDROGEN):
        return None  # bonded to another residue: hydrogens.py chose
    t = templates.bundled()["HIS"]
    by_name = {p.atom: p for p in t.hydrogens}
    x = s.coords.astype(np.float64)
    # free energies in water relative to HIE: HID by the tautomer ratio, HIP by HIE's pKa
    pka_hie = MODEL_PKA["HIS"] + np.log10(1 + 1 / HIE_PER_HID)
    intrinsic = {"HIE": 0.0, "HID": RT_LN10 * np.log10(HIE_PER_HID), "HIP": -RT_LN10 * pka_hie}
    states = []
    for flipped in (False, True):
        at = {"CB": x[:, here["CB"]], "CG": x[:, here["CG"]]}
        for name, source in zip(RING, FLIPPED_RING if flipped else RING, strict=True):
            at[name] = x[:, here[source]]
        for tautomer, polar in HIS_HYDROGENS.items():
            names, parents, xyz = list(RING), ["", "", "", ""], [at[a] for a in RING]
            for h in ("HD2", "HE1", *polar):
                pl = by_name[h]
                names.append(h)
                parents.append(pl.parent)
                xyz.append(_placed(pl, at))
            hip = tautomer == "HIP"
            formal = [1 if hip and n == "ND1" else 0 for n in names]
            partial = [0.5 if hip and n in ("ND1", "NE2") else 0.0 for n in names]
            label = f"{tautomer} flipped" if flipped else tautomer
            states.append(_state(label, names, parents, xyz, flipped=flipped, tautomer=tautomer,
                                 form=tautomer, protons=int(hip),
                                 intrinsic=intrinsic[tautomer] + KT * np.log(2),
                                 formal=formal, partial=partial))  # fmt: skip
    moving = _with_hydrogens(s, graph, [here[a] for a in RING])
    return _site(s, r, HIS, moving, [here["CG"]], states, graph)


# -- Asp and Glu -----------------------------------------------------------------------------------
def _acid(s: Structure, r: int, here: dict[str, int], graph: list[list[int]]) -> Site | None:
    name = str(s.res_name[r])
    carbon, first, second, h_first, h_second = ACIDS[name]
    if not all(a in here for a in (carbon, first, second)):
        return None
    if any(s.residue_of[j] != r for o in (first, second) for j in graph[here[o]]
           if s.element[j] != HYDROGEN):  # fmt: skip
        return None  # an ester or another link
    x = s.coords.astype(np.float64)
    at = {a: x[:, here[a]] for a in (carbon, first, second)}
    neutral = MODEL_PKA[name] - np.log10(2)  # two oxygens to carry the H
    states = [_state(name, [first, second], ["", ""], [at[first], at[second]], form=name,
                     formal=[0, -1], partial=[-0.5, -0.5])]  # fmt: skip
    form = {"ASP": "ASH", "GLU": "GLH"}[name]
    for oxygen, other, h in ((first, second, h_first), (second, first, h_second)):
        for anti in (False, True):
            xyz = templates.place(at[other], at[carbon], at[oxygen], OH_LENGTH, SYN_ANGLE,
                                  np.pi if anti else 0.0)  # fmt: skip
            label = f"{form} ({oxygen}{', anti' if anti else ''})"
            states.append(_state(label, [first, second, h], ["", "", oxygen],
                                 [at[first], at[second], xyz], anti=anti, form=form, protons=1,
                                 intrinsic=-RT_LN10 * neutral,
                                 formal=[0, 0, 0], partial=[0.0, 0.0, 0.0]))  # fmt: skip
    moving = _with_hydrogens(s, graph, [here[first], here[second]])
    return _site(s, r, ACID, moving, [here[carbon]], states, graph)


# -- rotors ----------------------------------------------------------------------------------------
def _rotor(
    s: Structure, rotor: hydrogens.Rotor, graph: list[list[int]], titrate: bool
) -> Site | None:
    parent, axis = (int(i) for i in s.indices(np.array([rotor.parent, rotor.axis])))
    hs = [int(i) for i in s.indices(np.array(rotor.hydrogens))]
    if min(parent, axis, *hs) < 0 or int(s.element[parent]) not in POLAR:
        return None
    others = [j for j in graph[axis] if j != parent and s.element[j] != HYDROGEN]
    if not others:
        return None
    ref = others[0]
    x = s.coords.astype(np.float64)
    first = x[0]
    length = [float(np.linalg.norm(first[h] - first[parent])) for h in hs]
    angle = [templates.angle_of(first[axis], first[parent], first[h]) for h in hs]
    dihedral = [templates.dihedral_of(first[ref], first[axis], first[parent], first[h]) for h in hs]
    offsets = [d - dihedral[0] for d in dihedral]
    if len(hs) == 1 and _planar(s, axis):
        starts = [0.0, np.pi]  # a phenol or acid OH lies in the plane of its neighbours
    else:
        period = 2 * np.pi / 3 if len(hs) == 3 else 2 * np.pi
        starts = list(np.arange(0.0, period - 1e-6, STEP))
    names = [str(s.name[h]) for h in hs]
    r = int(s.residue_of[parent])
    kind = _titrating(s, r, parent, len(hs)) if titrate else ""
    lead = [str(s.name[parent])] if kind else []  # a titrating rotor carries its heavy atom
    at_parent = [x[:, parent]] if kind else []
    on = [str(s.name[parent])] * len(hs)
    states = []
    for start in starts:
        xyz = [templates.place(x[:, ref], x[:, axis], x[:, parent], length[k], angle[k],
                               start + offsets[k]) for k in range(len(hs))]  # fmt: skip
        states.append(_state(f"{np.degrees(start):.0f}°", lead + names, [""] * len(lead) + on,
                             at_parent + xyz, donor_xyz=x[0, parent]))  # fmt: skip
    moving = [parent, *hs] if kind else hs
    if kind:
        states = _titrated(states, kind, s, parent, axis, ref, names, x)
    return _site(s, r, ROTOR, np.array(moving), [parent, axis], states, graph)


def _titrating(s: Structure, r: int, parent: int, n_h: int) -> str:
    """The residue kind (LYS, CYS, TYR or NTERM) when this rotor titrates, else ""."""
    name, atom = str(s.res_name[r]), str(s.name[parent])
    if (name, atom) in TITRATING:
        return name
    if atom == "N" and n_h == 3 and s.res_kind[r] == AMINO:
        return "NTERM"
    return ""


def _titrated(
    states: list[State],
    kind: str,
    s: Structure,
    parent: int,
    axis: int,
    ref: int,
    names: list[str],
    x: np.ndarray,
) -> list[State]:
    """The rotor's protonated states with their charges and free energies, and the
    deprotonated form: S⁻ or O⁻ (one state), or NH₂ with its lone pair in 15° steps."""
    heavy = str(s.name[parent])
    pka = MODEL_PKA[kind]
    base = kind in ("LYS", "NTERM")
    up, down = {"LYS": ("LYS", "LYN"), "NTERM": ("NH3+", "NH2"), "CYS": ("CYS", "CYM"),
                "TYR": ("TYR", "TYM")}[kind]  # fmt: skip
    q = 1 if base else 0
    variants = KT * np.log(len(states))
    protonated = [
        _state(st.label, list(st.names), list(st.parents), list(np.moveaxis(st.xyz, 1, 0)),
               donor_xyz=x[0, parent], form=up, protons=1,
               intrinsic=-RT_LN10 * pka + variants,
               formal=[q] + [0] * (len(st.names) - 1),
               partial=[float(q)] + [0.0] * (len(st.names) - 1))
        for st in states
    ]  # fmt: skip
    if not base:
        lone = [_state(down, [heavy], [""], [x[:, parent]], form=down, formal=[-1],
                       partial=[-1.0])]  # fmt: skip
        return protonated + lone
    neutral = []
    steps = np.arange(0.0, 2 * np.pi - 1e-6, STEP)
    for phi in steps:  # the lone pair at phi, the two hydrogens 120° either side
        xyz = [templates.place(x[:, ref], x[:, axis], x[:, parent], NH_LENGTH, TETRAHEDRAL,
                               phi + turn) for turn in (2 * np.pi / 3, 4 * np.pi / 3)]  # fmt: skip
        neutral.append(_state(f"{down} {np.degrees(phi):.0f}°", [heavy, *names[:2]],
                              ["", heavy, heavy], [x[:, parent], *xyz], donor_xyz=x[0, parent],
                              form=down, intrinsic=KT * np.log(len(steps)), formal=[0, 0, 0],
                              partial=[0.0, 0.0, 0.0]))  # fmt: skip
    return protonated + neutral


def _planar(s: Structure, atom: int) -> bool:
    on = (s.bonds == atom).any(axis=1)
    return bool(((s.order[on] >= 2) | (s.order[on] == AROMATIC)).any())


# -- helpers ---------------------------------------------------------------------------------------
def _state(
    label: str,
    names: list[str],
    parents: list[str],
    xyz: list[np.ndarray],
    flipped: bool = False,
    tautomer: str = "",
    anti: bool = False,
    donor_xyz: np.ndarray | None = None,
    form: str = "",
    protons: int = 0,
    intrinsic: float = 0.0,
    formal: list[int] | None = None,
    partial: list[float] | None = None,
) -> State:
    stacked = np.stack(xyz, axis=1)  # (k, m, 3)
    elements = np.array(
        [HYDROGEN if p else _element(n) for n, p in zip(names, parents, strict=True)]
    )
    donor = np.full((len(names), 3), np.nan)
    for k, parent in enumerate(parents):
        if not parent:
            continue
        if donor_xyz is not None:
            donor[k] = donor_xyz
        elif parent in names and _element(parent) in POLAR:
            donor[k] = stacked[0, names.index(parent)]
    return State(label, tuple(names), tuple(parents), elements, stacked, donor, flipped, anti,
                 tautomer, form, protons, float(intrinsic), tuple(formal or ()),
                 tuple(partial or ()))  # fmt: skip


def _placed(pl: templates.Placement, at: dict[str, np.ndarray]) -> np.ndarray:
    """Hydrogen *pl* from the positions *at* (in every model)."""
    return templates.place(
        at[pl.dihedral_atom], at[pl.angle_atom], at[pl.parent], pl.length, pl.angle, pl.dihedral
    )


def _element(name: str) -> int:
    return {"N": 7, "O": 8, "S": 16, "C": 6}[name[0]]


def _named(s: Structure, r: int) -> dict[str, int]:
    active = s.active()
    span = s.atoms_of(r)
    return {str(s.name[i]): i for i in range(span.start, span.stop) if active[i]}


def _graph(s: Structure) -> list[list[int]]:
    graph: list[list[int]] = [[] for _ in range(s.n_atoms)]
    for a, b in s.bonds.tolist():
        graph[a].append(b)
        graph[b].append(a)
    return graph


def _with_hydrogens(s: Structure, graph: list[list[int]], heavy: list[int]) -> np.ndarray:
    hs = [j for i in heavy for j in graph[i] if s.element[j] == HYDROGEN]
    return np.array(heavy + hs, dtype=np.int64)


def _site(
    s: Structure,
    r: int,
    kind: str,
    moving: np.ndarray,
    anchors: list[int],
    states: list[State],
    graph: list[list[int]],
) -> Site:
    """The site, with the atoms near its anchors in the bond graph left out of its score."""
    seen = {a: 0 for a in anchors}
    queue = deque(anchors)
    while queue:
        a = queue.popleft()
        if seen[a] >= EXCLUDED_BONDS:
            continue
        for b in graph[a]:
            if b not in seen:
                seen[b] = seen[a] + 1
                queue.append(b)
    return Site(r, kind, np.asarray(moving, dtype=np.int64), np.array(anchors), tuple(states),
                frozenset(seen))  # fmt: skip
