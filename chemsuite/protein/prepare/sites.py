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

Every state lists the atoms that move with it (heavy atoms and hydrogens, in every model), so
the energy can score it and the solver can write the choice back. The flips follow Reduce
(Word et al., *J. Mol. Biol.* 285, 1735–1747, 1999).
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

import numpy as np

from ..structure.model import AMINO, AROMATIC, Structure
from . import hydrogens, templates
from .hydrogens import HYDROGEN, Protonated

AMIDES = {"ASN": ("CG", "OD1", "ND2", "CB"), "GLN": ("CD", "OE1", "NE2", "CG")}
RING = ("ND1", "CD2", "CE1", "NE2")
FLIPPED_RING = ("CD2", "ND1", "NE2", "CE1")  # where each ring atom goes when the ring flips
HIS_HYDROGENS = {"HID": ("HD1",), "HIE": ("HE2",), "HIP": ("HD1", "HE2")}
STEP = np.radians(15.0)
POLAR = frozenset([7, 8, 16])
EXCLUDED_BONDS = 3  # atoms this many bonds or fewer from a site's anchors are not scored


@dataclass(frozen=True)
class State:
    label: str
    names: tuple[str, ...]  # the atoms of the state
    parents: tuple[str, ...]  # each hydrogen's heavy atom ("" for a heavy atom)
    elements: np.ndarray  # (m,)
    xyz: np.ndarray  # (k, m, 3): every model
    donor: np.ndarray  # (m, 3): the heavy atom a polar hydrogen sits on, NaN otherwise
    flipped: bool = False
    tautomer: str = ""  # HID, HIE or HIP for a His

    @property
    def polar_h(self) -> np.ndarray:
        return ~np.isnan(self.donor[:, 0])

    @property
    def acceptor(self) -> np.ndarray:
        """O, and N with no hydrogen in this state."""
        has_h = {p for p in self.parents if p}
        return np.array([(e == 8 or e == 7) and n not in has_h
                         for n, e in zip(self.names, self.elements, strict=True)])  # fmt: skip


@dataclass(frozen=True)
class Site:
    residue: int
    kind: str  # AMIDE, HIS or ROTOR
    moving: np.ndarray  # structure atoms the states replace
    anchors: np.ndarray  # structure atoms the moving ones hang from
    states: tuple[State, ...]
    excluded: frozenset[int] = field(default=frozenset())  # atoms not scored against it


AMIDE, HIS, ROTOR = "amide", "his", "rotor"


def find_sites(p: Protonated) -> list[Site]:
    """Every site of the protonated structure *p*."""
    s = p.structure
    graph = _graph(s)
    found: list[Site] = []
    for r in np.flatnonzero(s.res_kind == AMINO).tolist():
        name = str(s.res_name[r])
        here = _named(s, r)
        if name in AMIDES:
            site = _amide(s, r, here, graph)
        elif name == "HIS":
            site = _histidine(s, r, here, graph)
        else:
            continue
        if site is not None:
            found.append(site)
    for rotor in p.rotors:
        site = _rotor(s, rotor, graph)
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
            label = f"{tautomer} flipped" if flipped else tautomer
            states.append(_state(label, names, parents, xyz, flipped=flipped, tautomer=tautomer))
    moving = _with_hydrogens(s, graph, [here[a] for a in RING])
    return _site(s, r, HIS, moving, [here["CG"]], states, graph)


# -- rotors ----------------------------------------------------------------------------------------
def _rotor(s: Structure, rotor: hydrogens.Rotor, graph: list[list[int]]) -> Site | None:
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
    states = []
    for start in starts:
        xyz = [templates.place(x[:, ref], x[:, axis], x[:, parent], length[k], angle[k],
                               start + offsets[k]) for k in range(len(hs))]  # fmt: skip
        states.append(_state(f"{np.degrees(start):.0f}°", names, [str(s.name[parent])] * len(hs),
                             xyz, donor_xyz=x[0, parent]))  # fmt: skip
    r = int(s.residue_of[parent])
    return _site(s, r, ROTOR, np.array(hs), [parent, axis], states, graph)


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
    donor_xyz: np.ndarray | None = None,
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
    return State(label, tuple(names), tuple(parents), elements, stacked, donor, flipped, tautomer)


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
