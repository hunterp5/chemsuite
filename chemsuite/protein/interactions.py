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
"""Non-covalent interactions between a query (a ligand, a peptide, any selection) and the atoms
around it: hydrogen bonds, salt bridges, π-stacking, cation–π, hydrophobic contacts, halogen
bonds, metal contacts and water bridges.

Structures from the PDB rarely carry hydrogens, so everything is judged from heavy atoms.
Amino acids get their features from atom names; ligands and other small molecules from their
RDKit molecule with bond orders (:func:`.ligands.to_mol`); waters are donors and acceptors.
A hydrogen bond needs the donor–acceptor distance and, at both ends, the partner where a
hydrogen or lone pair can point (:func:`_points_at`); one to a water is shown as part of a
water bridge, one between salt-bridge atoms as the salt bridge. Cut-offs follow PLIP (Salentin
et al. 2015). Hydrophobic contacts and charged pairs keep only the closest pair per residue
pair, so a pocket does not fill with lines.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from . import ligands
from .model import ADDITIVE, AMINO, ION, WATER, Structure

# Kinds
HBOND = 0
SALT_BRIDGE = 1
PI_STACK = 2
CATION_PI = 3
HYDROPHOBIC = 4
HALOGEN = 5
METAL = 6
WATER_BRIDGE = 7
NAMES = {
    HBOND: "Hydrogen Bond",
    SALT_BRIDGE: "Salt Bridge",
    PI_STACK: "π-Stacking",
    CATION_PI: "Cation–π",
    HYDROPHOBIC: "Hydrophobic",
    HALOGEN: "Halogen Bond",
    METAL: "Metal",
    WATER_BRIDGE: "Water Bridge",
}
COLOURS = {
    HBOND: (0.20, 0.50, 1.00),
    SALT_BRIDGE: (0.90, 0.30, 0.80),
    PI_STACK: (0.20, 0.75, 0.35),
    CATION_PI: (1.00, 0.55, 0.10),
    HYDROPHOBIC: (0.60, 0.60, 0.60),
    HALOGEN: (0.10, 0.75, 0.80),
    METAL: (0.55, 0.35, 0.85),
    WATER_BRIDGE: (0.45, 0.75, 1.00),
}

# Cut-offs, Å and degrees
HBOND_MIN, HBOND_MAX = 2.5, 3.5
HBOND_NEIGHBOUR_ANGLE = 90.0
HBOND_BISECTOR_ANGLE = 50.0
SALT_BRIDGE_MAX = 4.0
PI_STACK_MAX = 5.5
PI_OFFSET_MAX = 2.0
PI_PARALLEL, PI_T_SHAPED = 30.0, 60.0
CATION_PI_MAX = 6.0
HYDROPHOBIC_MAX = 4.0
HALOGEN_MAX = 4.0
HALOGEN_ANGLE = 140.0  # C–X···A
METAL_MAX = 2.8
WATER_MAX = 3.5
NEAR = 7.0  # partners are looked for this far from the query

# Amino-acid features by atom name ("*": every amino acid)
_DONORS = {"*": {"N"}, "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
           "HIS": {"ND1", "NE2"}, "LYS": {"NZ"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
           "TRP": {"NE1"}, "CYS": {"SG"}}  # fmt: skip
_ACCEPTORS = {"*": {"O", "OXT"}, "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"},
              "GLN": {"OE1"}, "HIS": {"ND1", "NE2"}, "SER": {"OG"}, "THR": {"OG1"},
              "TYR": {"OH"}, "MET": {"SD"}, "CYS": {"SG"}}  # fmt: skip
_POSITIVE = {"ARG": {"NE", "NH1", "NH2"}, "LYS": {"NZ"}}
_NEGATIVE = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "*": {"OXT"}}
_BENZENE = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
_RINGS = {
    "PHE": (_BENZENE,),
    "TYR": (_BENZENE,),
    "HIS": (("CG", "ND1", "CD2", "CE1", "NE2"),),
    "TRP": (("CG", "CD1", "NE1", "CE2", "CD2"), ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")),
}
#: Li Na Mg K Ca Mn Fe Co Ni Cu Zn Cd Hg
_METALS = {3, 11, 12, 19, 20, 25, 26, 27, 28, 29, 30, 48, 80}
_HALOGENS = {17, 35, 53}  # Cl Br I (fluorine makes no halogen bond)
_POSITIVE_SMARTS = [
    Chem.MolFromSmarts(s)
    for s in (
        "[NX3;H2,H1;!$(N-[C,S,P]=[O,S,N]);!$(N-a);!$(N-[#7,#8])]",  # primary, secondary amines
        "[NX3;H0;!$(N-[C,S,P]=[O,S,N]);!$(N-a);!$(N-[#7,#8]);!$(N-C=C)]",  # tertiary amines
        "[NX3][CX3](=[NX2;!a])[NX3]",  # guanidine
        "[CX3;!$(C-[NX3]-[#6]=O)](=[NX2;!a])[NX3;!a]",  # amidine
        "[#7+;!$([#7+]-[O-]);!$([#7+]=O)]",
    )
]
_NEGATIVE_SMARTS = [
    Chem.MolFromSmarts(s)
    for s in (
        "[CX3](=O)[OX2H1,OX1-]",  # carboxylic acid
        "[SX4,PX4](=O)[OX2H1,OX1-]",  # sulfonic, phosphonic, phosphate
        "c1nn[nH]n1",  # tetrazole
        "[#8-,#16-;!$([#8-]-[#7+])]",
    )
]


@dataclass
class Features:
    """What a set of atoms can do (atom indices into the structure)."""

    donors: set[int] = field(default_factory=set)
    acceptors: set[int] = field(default_factory=set)
    positive: dict[int, int] = field(default_factory=dict)  # atom → charged group
    negative: dict[int, int] = field(default_factory=dict)
    rings: list[tuple[int, ...]] = field(default_factory=list)
    hydrophobic: set[int] = field(default_factory=set)
    halogens: set[int] = field(default_factory=set)
    metals: set[int] = field(default_factory=set)
    waters: set[int] = field(default_factory=set)


@dataclass(frozen=True)
class Interactions:
    """``m`` interactions; ``a`` is the query's atom (a ring's first atom), ``b`` the partner's."""

    kind: np.ndarray  # int8 (m,)
    a: np.ndarray  # int64 (m,)
    b: np.ndarray  # int64 (m,)
    start: np.ndarray  # float32 (m, 3): the query end (an atom or a ring centre)
    end: np.ndarray  # float32 (m, 3)
    distance: np.ndarray  # float32 (m,)

    @classmethod
    def empty(cls) -> Interactions:
        return cls(np.zeros(0, np.int8), np.zeros(0, np.int64), np.zeros(0, np.int64),
                   np.zeros((0, 3), np.float32), np.zeros((0, 3), np.float32),
                   np.zeros(0, np.float32))  # fmt: skip

    def __len__(self) -> int:
        return len(self.kind)

    def counts(self) -> dict[int, int]:
        kinds, counts = np.unique(self.kind, return_counts=True)
        return {int(k): int(c) for k, c in zip(kinds, counts, strict=True)}


def describe(structure: Structure, found: Interactions, k: int) -> str:
    """One interaction in words: ``Hydrogen Bond  ATP E355 N6 – GLU E121 O  2.9 Å``."""
    return (f"{NAMES[int(found.kind[k])]}  {structure.atom_label(int(found.a[k]))} – "
            f"{structure.atom_label(int(found.b[k]))}  {found.distance[k]:.1f} Å")  # fmt: skip


# -- features -----------------------------------------------------------------------------------


def features(
    structure: Structure,
    atoms: np.ndarray,
    active: np.ndarray | None = None,
    model: int = 0,
    smiles: dict[int, str] | None = None,
) -> Features:
    """The features of ``atoms`` (a bool mask): amino acids by name, waters, ions, and other
    residues through RDKit. ``smiles`` gives a residue's bond orders when known."""
    s = structure
    if active is None:
        active = s.active()
    atoms = atoms & active & (s.element != 1)
    out = Features()
    neighbours = _heavy_neighbours(s)
    for residue in np.unique(s.residue_of[atoms]).tolist():
        span = s.atoms_of(residue)
        mine = span.start + np.flatnonzero(atoms[span])
        kind = int(s.res_kind[residue])
        if kind == AMINO:
            _amino(s, residue, mine, out)
        elif kind == WATER:
            out.waters.update(i for i in mine.tolist() if s.element[i] == 8)
        elif kind == ION:
            out.metals.update(i for i in mine.tolist() if int(s.element[i]) in _METALS)
        else:
            _small_molecule(s, residue, mine, active, model, (smiles or {}).get(residue), out)
    for i in np.flatnonzero(atoms).tolist():
        element = int(s.element[i])
        near = neighbours.get(i, ())
        if element == 6 and all(s.element[j] in (1, 6) for j in near):
            out.hydrophobic.add(i)
        elif element in _HALOGENS:
            out.hydrophobic.add(i)
            if any(s.element[j] == 6 for j in near):
                out.halogens.add(i)
    return out


def _amino(s: Structure, residue: int, mine: np.ndarray, out: Features) -> None:
    name = str(s.res_name[residue])
    by_name = {str(s.name[i]): i for i in mine.tolist()}
    for table, target in ((_DONORS, out.donors), (_ACCEPTORS, out.acceptors)):
        for atom_name in table.get("*", set()) | table.get(name, set()):
            if atom_name in by_name and not (name == "PRO" and atom_name == "N"):
                target.add(by_name[atom_name])
    for table, charged in ((_POSITIVE, out.positive), (_NEGATIVE, out.negative)):
        for atom_name in table.get("*", set()) | table.get(name, set()):
            if atom_name in by_name:
                charged[by_name[atom_name]] = residue
    for ring in _RINGS.get(name, ()):
        if all(n in by_name for n in ring):
            out.rings.append(tuple(by_name[n] for n in ring))


def _small_molecule(s, residue, mine, active, model, smiles, out: Features) -> None:
    try:
        mol = ligands.to_mol(s, residue, smiles=smiles, active=active, model=model).mol
    except (ValueError, RuntimeError):
        return
    index = [s.index_of_uid[a.GetIntProp("uid")] for a in mol.GetAtoms()]
    keep = set(mine.tolist())
    _polar(mol, index, keep, out)
    _charged_groups(mol, index, keep, residue, out)
    for ring in mol.GetRingInfo().AtomRings():
        members = tuple(index[m] for m in ring)
        aromatic = all(mol.GetAtomWithIdx(m).GetIsAromatic() for m in ring)
        if len(ring) in (5, 6) and aromatic and all(m in keep for m in members):
            out.rings.append(members)


def _polar(mol: Chem.Mol, index: list[int], keep: set[int], out: Features) -> None:
    """Donors: N and O with a hydrogen. Acceptors: O, and N with a free lone pair (no
    hydrogen, fewer than three neighbours), neither positively charged."""
    for atom in mol.GetAtoms():
        i = index[atom.GetIdx()]
        element, hydrogens = atom.GetAtomicNum(), atom.GetTotalNumHs()
        if i not in keep or element not in (7, 8):
            continue
        if hydrogens:
            out.donors.add(i)
        lone_pair = element == 8 or (not hydrogens and atom.GetDegree() < 3)
        if lone_pair and atom.GetFormalCharge() <= 0:
            out.acceptors.add(i)


def _charged_groups(mol, index, keep, residue: int, out: Features) -> None:
    """Groups charged at pH 7 (amines, amidines, acids…), each its own group number
    (negative, so never equal to an amino acid's residue index)."""
    group = 0
    for patterns, charged in ((_POSITIVE_SMARTS, out.positive), (_NEGATIVE_SMARTS, out.negative)):
        for pattern in patterns:
            for match in mol.GetSubstructMatches(pattern):
                group += 1
                for m in match:
                    element = mol.GetAtomWithIdx(m).GetAtomicNum()
                    if element in (7, 8, 16) and index[m] in keep:
                        charged.setdefault(index[m], -(residue * 64 + group) - 1)


def _heavy_neighbours(s: Structure) -> dict[int, list[int]]:
    out: dict[int, list[int]] = {}
    heavy = s.element != 1
    for a, b in s.bonds.tolist():
        if heavy[a] and heavy[b]:
            out.setdefault(a, []).append(b)
            out.setdefault(b, []).append(a)
    return out


# -- interactions -------------------------------------------------------------------------------


def find(
    structure: Structure,
    query: np.ndarray,
    active: np.ndarray | None = None,
    model: int = 0,
    smiles: dict[int, str] | None = None,
) -> Interactions:
    """The interactions between the ``query`` atoms (a bool mask) and every other atom near."""
    s = structure
    if active is None:
        active = s.active()
    query = query & active & (s.element != 1)
    if not query.any():
        return Interactions.empty()
    xyz = s.xyz(model).astype(np.float64)
    tree = cKDTree(xyz[~query & active])
    others = np.flatnonzero(~query & active)
    near = np.zeros(s.n_atoms, bool)
    for hits in tree.query_ball_point(xyz[query], NEAR):
        near[others[hits]] = True
    q = features(s, query, active, model, smiles)
    p = features(s, near, active, model, smiles)
    bonded = {tuple(sorted(pair)) for pair in s.bonds.tolist()}
    neighbours = _heavy_neighbours(s)
    found = _Found(s, xyz)
    _charged(found, q.positive, p.negative, SALT_BRIDGE)
    _charged(found, q.negative, p.positive, SALT_BRIDGE)
    _hbonds(found, q, p, neighbours, bonded)
    _stacking(found, q.rings, p.rings)
    _cation_pi(found, q.positive, p.rings, query_is_cation=True)
    _cation_pi(found, p.positive, q.rings, query_is_cation=False)
    _hydrophobic(found, q.hydrophobic, p.hydrophobic, bonded)
    _halogens(found, q.halogens, p.acceptors | p.waters, neighbours, query_is_donor=True)
    _halogens(found, p.halogens, q.acceptors | q.waters, neighbours, query_is_donor=False)
    _metals(found, q, p)
    _water_bridges(found, q, p, query)
    return found.result()


class _Found:
    def __init__(self, structure: Structure, xyz: np.ndarray) -> None:
        self.s = structure
        self.xyz = xyz
        self.rows: list[tuple[int, int, int, np.ndarray, np.ndarray]] = []

    def add(self, kind: int, a: int, b: int, start=None, end=None) -> None:
        start = self.xyz[a] if start is None else start
        end = self.xyz[b] if end is None else end
        self.rows.append((kind, a, b, start, end))

    def result(self) -> Interactions:
        if not self.rows:
            return Interactions.empty()
        kind, a, b, start, end = zip(*self.rows, strict=True)
        start, end = np.array(start, np.float32), np.array(end, np.float32)
        return Interactions(np.array(kind, np.int8), np.array(a, np.int64),
                            np.array(b, np.int64), start, end,
                            np.linalg.norm(end - start, axis=1).astype(np.float32))  # fmt: skip


def _pairs(xyz, first, second, cutoff: float):
    """(i, j, distance) for atoms of ``first`` and ``second`` closer than ``cutoff``."""
    first, second = sorted(first), sorted(second)
    if not first or not second:
        return
    a, b = np.array(first), np.array(second)
    tree = cKDTree(xyz[b])
    for k, hits in enumerate(tree.query_ball_point(xyz[a], cutoff)):
        for h in hits:
            yield int(a[k]), int(b[h]), float(np.linalg.norm(xyz[a[k]] - xyz[b[h]]))


def _angle(a, b, c) -> float:
    """The angle at ``b`` in degrees."""
    u, v = a - b, c - b
    cos = u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-9)
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _points_at(xyz, neighbours, atom: int, partner: int) -> bool:
    """Whether ``atom``'s hydrogen or lone pair can point at ``partner``: with one heavy
    neighbour, the partner is more than 90° from it; with more, the partner is within 50° of
    the direction away from them all (the bisector, where an sp2 N–H or lone pair points)."""
    near = neighbours.get(atom, ())
    if len(near) == 1:
        return _angle(xyz[near[0]], xyz[atom], xyz[partner]) > HBOND_NEIGHBOUR_ANGLE
    if not near:
        return True
    away = xyz[atom] - xyz[list(near)]
    away /= np.maximum(np.linalg.norm(away, axis=1, keepdims=True), 1e-9)
    direction = away.sum(axis=0)
    if np.linalg.norm(direction) < 0.2:  # planar with three neighbours: nothing to give
        return False
    return _angle(xyz[atom] + direction, xyz[atom], xyz[partner]) <= HBOND_BISECTOR_ANGLE


def _hbonds(found: _Found, q: Features, p: Features, neighbours, bonded) -> None:
    xyz = found.xyz
    seen = {(a, b) for kind, a, b, _s, _e in found.rows if kind == SALT_BRIDGE}
    for donors, acceptors in ((q.donors, p.acceptors), (q.acceptors, p.donors)):
        for i, j, d in _pairs(xyz, donors, acceptors, HBOND_MAX):
            if (d < HBOND_MIN or (min(i, j), max(i, j)) in bonded or (i, j) in seen
                    or found.s.res_kind[found.s.residue_of[j]] == WATER):  # fmt: skip
                continue
            if _points_at(xyz, neighbours, i, j) and _points_at(xyz, neighbours, j, i):
                seen.add((i, j))
                found.add(HBOND, i, j)


def _charged(found: _Found, mine: dict[int, int], theirs: dict[int, int], kind: int) -> None:
    best: dict[tuple[int, int], tuple[float, int, int]] = {}
    for i, j, d in _pairs(found.xyz, mine, theirs, SALT_BRIDGE_MAX):
        key = (mine[i], theirs[j])
        if key not in best or d < best[key][0]:
            best[key] = (d, i, j)
    for _d, i, j in best.values():
        found.add(kind, i, j)


def _ring_geometry(xyz, ring: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    points = xyz[list(ring)]
    centre = points.mean(axis=0)
    normal = np.linalg.svd(points - centre)[2][2]
    return centre, normal


def _offset(point, centre, normal) -> float:
    """How far ``point`` is from the ring's axis."""
    v = point - centre
    return float(np.linalg.norm(v - (v @ normal) * normal))


def _stacking(found: _Found, mine, theirs) -> None:
    xyz = found.xyz
    rings = [_ring_geometry(xyz, r) for r in theirs]
    for ring in mine:
        c1, n1 = _ring_geometry(xyz, ring)
        for other, (c2, n2) in zip(theirs, rings, strict=True):
            d = float(np.linalg.norm(c2 - c1))
            if d > PI_STACK_MAX:
                continue
            angle = float(np.degrees(np.arccos(min(abs(float(n1 @ n2)), 1.0))))
            offset = min(_offset(c2, c1, n1), _offset(c1, c2, n2))
            if offset <= PI_OFFSET_MAX and (angle < PI_PARALLEL or angle > PI_T_SHAPED):
                found.add(PI_STACK, ring[0], other[0], c1, c2)


def _cation_pi(found: _Found, cations: dict[int, int], rings, query_is_cation: bool) -> None:
    xyz = found.xyz
    best: dict[tuple[int, tuple], tuple[float, int]] = {}
    for ring in rings:
        centre, normal = _ring_geometry(xyz, ring)
        for atom, group in cations.items():
            d = float(np.linalg.norm(xyz[atom] - centre))
            if d <= CATION_PI_MAX and _offset(xyz[atom], centre, normal) <= PI_OFFSET_MAX:
                key = (group, ring)
                if key not in best or d < best[key][0]:
                    best[key] = (d, atom)
    for (_group, ring), (_d, atom) in best.items():
        centre = _ring_geometry(xyz, ring)[0]
        if query_is_cation:
            found.add(CATION_PI, atom, ring[0], None, centre)
        else:
            found.add(CATION_PI, ring[0], atom, centre, None)


def _hydrophobic(found: _Found, mine, theirs, bonded) -> None:
    s = found.s
    best: dict[tuple[int, int], tuple[float, int, int]] = {}
    for i, j, d in _pairs(found.xyz, mine, theirs, HYDROPHOBIC_MAX):
        if (min(i, j), max(i, j)) in bonded:
            continue
        key = (int(s.residue_of[i]), int(s.residue_of[j]))
        if key not in best or d < best[key][0]:
            best[key] = (d, i, j)
    for _d, i, j in best.values():
        found.add(HYDROPHOBIC, i, j)


def _halogens(found: _Found, halogens, acceptors, neighbours, query_is_donor: bool) -> None:
    xyz = found.xyz
    for x, a, _d in _pairs(xyz, halogens, acceptors, HALOGEN_MAX):
        carbons = [n for n in neighbours.get(x, ()) if found.s.element[n] == 6]
        if carbons and _angle(xyz[carbons[0]], xyz[x], xyz[a]) >= HALOGEN_ANGLE:
            found.add(HALOGEN, *((x, a) if query_is_donor else (a, x)))


def _metals(found: _Found, q: Features, p: Features) -> None:
    s = found.s
    ligating = (s.element == 7) | (s.element == 8) | (s.element == 16)
    mine = {i for i in q.donors | q.acceptors | q.waters | q.negative.keys() if ligating[i]}
    theirs = {i for i in p.donors | p.acceptors | p.waters | p.negative.keys() if ligating[i]}
    for i, j, _d in _pairs(found.xyz, q.metals, theirs, METAL_MAX):
        found.add(METAL, i, j)
    for i, j, _d in _pairs(found.xyz, mine, p.metals, METAL_MAX):
        found.add(METAL, i, j)


def _water_bridges(found: _Found, q: Features, p: Features, query: np.ndarray) -> None:
    """A water within hydrogen-bond reach of both the query and a non-water partner."""
    s = found.s
    xyz = found.xyz
    polar = q.donors | q.acceptors
    partners = {i for i in p.donors | p.acceptors
                if s.res_kind[s.residue_of[i]] not in (WATER, ADDITIVE)}  # fmt: skip
    reach: dict[int, list[int]] = {}
    for i, w, d in _pairs(xyz, polar, p.waters, WATER_MAX):
        if d >= HBOND_MIN:
            reach.setdefault(w, []).append(i)
    if not reach:
        return
    for w, j, d in _pairs(xyz, set(reach), partners, WATER_MAX):
        if d < HBOND_MIN or query[j]:
            continue
        for i in reach.pop(w, []):
            found.add(WATER_BRIDGE, i, w)
        found.add(WATER_BRIDGE, w, j)


def ligand_query(structure: Structure, residue: int) -> np.ndarray:
    """The mask of one residue's atoms, the usual query."""
    mask = np.zeros(structure.n_atoms, bool)
    mask[structure.atoms_of(residue)] = True
    return mask
