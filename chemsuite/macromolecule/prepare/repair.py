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

"""Repairs before hydrogens: one location per atom, and missing heavy atoms built.

- **Alternate locations**: the location shown (the highest occupancy, or the user's pick) is
  kept and the others dropped.
- **Missing atoms**: a polymer residue's missing side-chain atoms, and the OXT of a chain's
  last residue, are built from the residue's template, one at a time, each from three atoms
  already there by the template's internal coordinates (bond length, angle and dihedral of
  the CCD's ideal geometry). The side chain takes the ideal conformation, which is not
  checked for clashes. No loops are built.
- **Stubs**: an amino acid missing CA or C (or, at the C-side end of a stretch, its O) has
  nothing to build on or cap from, so it is removed (:func:`remove_stubs`). Its neighbours
  then end at a gap, which is filled or capped like any other (4AGC's last residue, GLN 1169,
  is deposited as its N alone).
- **Caps**: a chain end that the deposit cut short (at a gap, or before or after a stretch of
  the sequence that was not modelled) takes a neutral ACE or NME cap (:func:`cap_ends`), as
  a force field needs; true termini stay charged.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.spatial import cKDTree

from chemsuite.core.jobs import CancelToken, check

from ..structure import edit
from ..structure.model import AMINO, ION, Structure
from .hydrogens import HYDROGEN
from .inspect import BACKBONE, CAPS, POLYMER, ends, last_in_chain, missing_atoms
from .templates import Library, Template, place

ACE, NME = "ACE", "NME"
PEPTIDE = 1.335  # Å, the C–N peptide bond
N_CA, CA_C = 1.458, 1.525  # Å
SPAN = 2.46  # Å, N···C across a residue's CA (N–CA–C at 111°)
SPAN_SLACK = 0.4  # Å: a gap that cannot be spanned this closely is capped instead
ROOM = 2.8  # Å: a built CA or CB closer than this to another atom counts against its place
CONTACT = 3.2  # Å: a built side-chain atom closer than this to another residue's (not an
# ion's) is crowded
TURN = 30  # degrees, the steps a built side chain is turned through
BACKBONE_ATOMS = frozenset(["N", "CA", "C", "O", "OXT"])


def one_location(s: Structure, altlocs: dict[int, str] | None = None) -> Structure:
    """*s* with only the atoms shown, none of them marked as an alternate location."""
    if (s.altloc == "").all():
        return s
    kept = s.subset(s.active(altlocs))
    return replace(kept, altloc=np.full(kept.n_atoms, "", dtype=kept.altloc.dtype),
                   occupancy=np.ones(kept.n_atoms, dtype=kept.occupancy.dtype))  # fmt: skip


def remove_stubs(s: Structure) -> tuple[Structure, list[str]]:
    """*s* without the amino acids too incomplete to build on or cap: those missing CA or C,
    and, at the C-side end of a stretch (which takes an OXT or an NME, placed from it), O.
    Inside a chain only a residue of the deposited sequence is removed, since its gap is then
    marked by the sequence and filled or capped; without one, only chain ends are. Repeats
    until the new ends are whole. Returns the labels of the residues removed."""
    removed: list[str] = []
    while True:
        found = ends(s)
        tails = found.c_termini | found.open_c
        stretch_ends = tails | found.n_termini | found.open_n
        gone = np.zeros(s.n_atoms, dtype=bool)
        now: list[str] = []
        for r in range(s.n_residues):
            if s.res_kind[r] != AMINO or str(s.res_name[r]) in CAPS:
                continue
            if r not in stretch_ends and s.seqres[r] < 0:
                continue
            span = s.atoms_of(r)
            names = {str(s.name[i]) for i in range(span.start, span.stop)
                     if s.element[i] != HYDROGEN}  # fmt: skip
            needed = ("CA", "C", "O") if r in tails else ("CA", "C")
            if any(name not in names for name in needed):
                gone[span] = True
                now.append(s.residue_label(r))
        if not now or gone.all():
            return s, removed
        removed += now
        s = edit.delete(s, gone)


def build_missing(
    s: Structure,
    library: Library,
    cancel: CancelToken | None = None,
    only: frozenset[int] | None = None,
) -> Structure:
    """*s* with the missing heavy atoms of its polymer residues built (uids of the old atoms
    kept). With *only*, just those residues are built: the rest are left as the file has
    them, for a chemist who would rather see a gap than a modelled side chain."""
    ends = last_in_chain(s)
    x = s.coords.astype(np.float64)
    residues, names, elements, coords, bonds, orders = [], [], [], [], [], []
    n = s.n_atoms
    heavy = np.flatnonzero(s.element != HYDROGEN)
    tree = cKDTree(x[0, heavy]) if len(heavy) else None
    for r in range(s.n_residues):
        if s.res_kind[r] not in POLYMER or (only is not None and r not in only):
            continue
        check(cancel)
        lacking = missing_atoms(s, r, library, ends, cancel)
        if not lacking or set(BACKBONE[int(s.res_kind[r])]) & set(lacking):
            continue
        t = library.get(str(s.res_name[r]), cancel)
        if t is None:
            continue
        span = s.atoms_of(r)
        index = {str(s.name[i]): i for i in range(span.start, span.stop)
                 if s.element[i] != HYDROGEN}  # fmt: skip
        at = {name: x[:, i] for name, i in index.items()}
        built = _build(t, at, lacking)
        if tree is not None:
            # not ions: a Cys or His that binds a metal is meant to sit close to it
            near = s.residue_of[heavy]
            others = heavy[(near != r) & (s.res_kind[near] != ION)]
            built = _turned_clear(t, at, built,
                                  lambda p, o=others: _crowding(tree, heavy, o, p))  # fmt: skip
        start = n + len(names)
        new_index = {name: start + k for k, name in enumerate(built)}
        for name, xyz in built.items():
            residues.append(r)
            names.append(name)
            elements.append(t.elements[t.index[name]])
            coords.append(xyz)
        everywhere = {**index, **new_index}
        for i, j, order in t.bonds:
            a, b = t.atoms[i], t.atoms[j]
            if (a in new_index or b in new_index) and a in everywhere and b in everywhere:
                bonds.append((everywhere[a], everywhere[b]))
                orders.append(order)
    if not names:
        return s
    return edit.add_atoms(
        s,
        residues=np.array(residues),
        names=np.array(names),
        elements=np.array(elements),
        xyz=np.stack(coords, axis=1),
        bonds=np.array(bonds, dtype=np.int64).reshape(-1, 2),
        order=np.array(orders),
    )


def _crowding(tree: cKDTree, heavy: np.ndarray, others: np.ndarray, p: np.ndarray) -> float:
    """How far the points *p* (m, 3) come inside :data:`CONTACT` of the atoms *others* (heavy
    atoms of other residues; *tree* holds every heavy atom, in the order of *heavy*):
    the sum of the squared overlaps."""
    allowed = np.zeros(len(heavy), dtype=bool)
    allowed[np.searchsorted(heavy, others)] = True
    total = 0.0
    for point, near in zip(p, tree.query_ball_point(p, CONTACT), strict=True):
        hits = [k for k in near if allowed[k]]
        if hits:
            d = np.linalg.norm(tree.data[hits] - point, axis=1)
            total += float(np.sum((CONTACT - d) ** 2))
    return total


def _turned_clear(
    t: Template, at: dict[str, np.ndarray], built: dict[str, np.ndarray], crowding
) -> dict[str, np.ndarray]:
    """*built* with its side chain turned about CA–CB, then about the bond beyond CB, in
    :data:`TURN` steps, to where it comes least close to other residues' atoms (first model;
    every model turned alike). Only atoms that were built move; a bond is turned only when all
    the atoms beyond it were built."""
    out = dict(built)
    for axis in _axes(t):
        a, b = axis
        beyond = _beyond(t, a, b)
        if not beyond or any(name not in out for name in beyond):
            continue
        placed = {**at, **out}
        if a not in placed or b not in placed:
            continue
        origin, direction = placed[b], placed[b] - placed[a]
        direction = direction / np.linalg.norm(direction, axis=-1, keepdims=True)
        best, chosen = np.inf, None
        for angle in np.radians(np.arange(0, 360, TURN)):
            moved = {name: _rotate(out[name], origin, direction, angle) for name in beyond}
            score = crowding(np.array([moved[name][0] for name in beyond]))
            if score < best - 1e-9:
                best, chosen = score, moved
        if chosen is not None:
            out.update(chosen)
    return out


def _axes(t: Template) -> list[tuple[str, str]]:
    """CA–CB, and CB to its one heavy neighbour other than CA, where the template has them."""
    if "CA" not in t.index or "CB" not in t.index:
        return []
    cb = t.index["CB"]
    out = [("CA", "CB")]
    beyond = [t.atoms[j] for j in t.heavy_neighbours(cb) if t.atoms[j] != "CA"]
    if len(beyond) == 1:
        out.append(("CB", beyond[0]))
    return out


def _beyond(t: Template, a: str, b: str) -> list[str]:
    """The heavy atoms reached from *b* without crossing back to *a* (a ring that closes on
    the backbone, as Pro's, gives none: it cannot turn)."""
    seen, todo = {t.index[b]}, [t.index[b]]
    while todo:
        i = todo.pop()
        for j in t.heavy_neighbours(i):
            if t.atoms[j] == a and i != t.index[b]:
                return []
            if t.atoms[j] != a and j not in seen:
                seen.add(j)
                todo.append(j)
    if any(t.atoms[j] in BACKBONE_ATOMS for j in seen if t.atoms[j] != b):
        return []
    return [t.atoms[j] for j in seen if t.atoms[j] != b]


def _rotate(p: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Points *p* (k, 3) turned by *angle* about the line through *origin* along *axis*
    (Rodrigues)."""
    v = p - origin
    cos, sin = np.cos(angle), np.sin(angle)
    dot = np.sum(v * axis, axis=-1, keepdims=True)
    return origin + v * cos + np.cross(axis, v) * sin + axis * dot * (1 - cos)


def fill_single_gaps(s: Structure) -> tuple[Structure, list[str]]:
    """*s* with the backbone (N, CA, C, O) of each residue that is alone in its gap in the
    deposited sequence, so its two neighbours need no caps (which would meet where it is);
    :func:`build_missing` then builds its side chain. Returns the labels of the residues built.

    The new N lies in the peptide plane of the residue before (trans), the new C in that of the
    residue after, turned about the latter's N–CA to the φ (every 5°) that puts the two
    :data:`SPAN` apart; CA closes the triangle at the angle nearest the ideal ones. A gap whose
    ends are too far apart or too close for one residue is left to the caps."""
    found = ends(s)
    pairs = [(r, r + 1) for r in sorted(found.open_c)
             if r + 1 in found.open_n and s.chain_of[r] == s.chain_of[r + 1]
             and _one_missing(s, r, r + 1)]  # fmt: skip
    labels = [(s.residue_label(a), s.residue_label(b)) for a, b in pairs]
    built: list[str] = []
    for before, after in labels:
        where = {s.residue_label(r): r for r in range(s.n_residues)}
        filled = _fill(s, where[before], where[after])
        if filled is not None:
            s = filled
            built.append(s.residue_label(where[after]))  # the new residue sits where *after* was
    return s, built


def _one_missing(s: Structure, a: int, b: int) -> bool:
    e = int(s.entity[a])
    return (s.seqres[a] >= 0 and s.seqres[b] == s.seqres[a] + 2 and e == s.entity[b]
            and 0 <= e < len(s.entity_id))  # fmt: skip


def _fill(s: Structure, a: int, b: int) -> Structure | None:
    heavy = [{str(s.name[i]): i for i in range(*s.atoms_of(r).indices(s.n_atoms))
              if s.element[i] != HYDROGEN} for r in (a, b)]  # fmt: skip
    if any(n not in heavy[0] for n in ("CA", "C", "O")) or any(
        n not in heavy[1] for n in ("N", "CA", "C")
    ):
        return None
    if "OXT" in heavy[0]:
        return None
    x = s.coords.astype(np.float64)
    ca0, c0, o0 = (x[:, heavy[0][n]] for n in ("CA", "C", "O"))
    n1, ca1, c1 = (x[:, heavy[1][n]] for n in ("N", "CA", "C"))
    nitrogen = place(o0, ca0, c0, PEPTIDE, np.radians(116.2), np.pi)
    options = []
    for phi in np.radians(np.arange(-180, 180, 5)):
        carbon = place(c1, ca1, n1, PEPTIDE, np.radians(121.7), phi)
        options.append((abs(float(np.linalg.norm(carbon[0] - nitrogen[0])) - SPAN), carbon))
    miss, carbon = min(options, key=lambda option: option[0])
    if miss > SPAN_SLACK:
        return None
    alpha = _apex(s, nitrogen, carbon, heavy[0]["C"], heavy[1]["N"])
    oxygen = place(n1, alpha, carbon, 1.229, np.radians(120.5), np.pi)
    seqres = int(s.seqres[a]) + 1
    entity = int(s.entity[a])
    name = str(s.sequence_of(entity)[seqres])
    seq = int(s.res_seq[a]) + 1
    chain = int(s.chain_of[a])
    n = s.n_atoms
    bonds = [[heavy[0]["C"], n], [n, n + 1], [n + 1, n + 2], [n + 2, n + 3], [n + 2, heavy[1]["N"]]]
    return edit.insert_residue(s, b, chain, name, seq, np.array(["N", "CA", "C", "O"]),
                               np.array([7, 6, 6, 8]),
                               np.stack([nitrogen, alpha, carbon, oxygen], axis=1),
                               np.array(bonds), np.array([1, 1, 1, 2, 1]),
                               icode=_free_icode(s, chain, seq), seqres=seqres,
                               entity=entity)  # fmt: skip


def _apex(s: Structure, n: np.ndarray, c: np.ndarray, before: int, after: int) -> np.ndarray:
    """CA (k, 3) at N–CA and CA–C bond lengths from *n* and *c*, turned about N···C to the
    angles C(*before*)–N–CA and CA–C–N(*after*) (atom indices) nearest their ideal values, each
    degree off costing one and each tenth of an ångström that CA or an L-amino acid's CB comes
    closer than :data:`ROOM` to another atom costing ten (first model)."""
    others = np.ones(s.n_atoms, dtype=bool)
    others[[before, after]] = False
    tree = cKDTree(s.coords[0, others].astype(np.float64))
    before_xyz, after_xyz = s.coords[0, before], s.coords[0, after]
    axis = c - n
    d = np.linalg.norm(axis, axis=-1, keepdims=True)
    unit = axis / d
    # the circle of points at N_CA from N and CA_C from C
    along = (d**2 + N_CA**2 - CA_C**2) / (2 * d)
    radius = np.sqrt(np.clip(N_CA**2 - along**2, 0.0, None))
    helper = np.where(np.abs(unit[..., :1]) < 0.9, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    u = np.cross(unit, helper)
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    v = np.cross(unit, u)
    best, chosen = np.inf, n
    for turn in np.radians(np.arange(0, 360, 5)):
        alpha = n + along * unit + radius * (np.cos(turn) * u + np.sin(turn) * v)
        beta = place(c[0], n[0], alpha[0], 1.53, np.radians(109.5), np.radians(-120.0))
        cost = (abs(_angle(before_xyz, n[0], alpha[0]) - 121.7)
                + abs(_angle(alpha[0], c[0], after_xyz) - 116.2))  # fmt: skip
        for atom in (alpha[0], beta):
            near = tree.query(atom)[0]
            cost += 100.0 * max(0.0, ROOM - float(near))
        if cost < best:
            best, chosen = cost, alpha
    return chosen


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u, w = a - b, c - b
    cos = float(np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


@dataclass(frozen=True)
class Cap:
    residue: str  # the label of the residue capped
    name: str  # ACE or NME


def cap_ends(s: Structure) -> tuple[Structure, list[Cap]]:
    """*s* with an ACE before every stretch of amino acids that starts after a gap or a
    missing N-terminal stretch, and an NME after every one that stops before one
    (:func:`.inspect.ends`), so no chain end the deposit cut short is left a charged terminus.

    A cap sits where the missing neighbour's backbone would: the NME's N and C in the peptide
    plane of the residue's C=O (trans), the ACE's C, O and CH3 in a trans peptide with the
    residue's N, turned about N–CA to the φ that comes nearest no other atom. A residue that
    lacks the backbone atoms for it, has an OXT, or is bonded onward already is left alone."""
    found = ends(s)
    plans: list[tuple[int, str]] = [(r, NME) for r in found.open_c]
    plans += [(r, ACE) for r in found.open_n]
    # in chain order, an NME before the ACE across its gap, so the ACE can turn away from it
    plans.sort(key=lambda plan: (plan[0], plan[1] == ACE))
    labels = [(s.residue_label(r), kind) for r, kind in plans]
    caps: list[Cap] = []
    for label, kind in labels:
        r = next(r for r in range(s.n_residues) if s.residue_label(r) == label)
        built = _cap(s, r, kind)
        if built is not None:
            caps.append(Cap(label, kind))
            s = built
    return s, caps


def _cap(s: Structure, r: int, kind: str) -> Structure | None:
    span = s.atoms_of(r)
    index = {str(s.name[i]): i for i in range(span.start, span.stop) if s.element[i] != HYDROGEN}
    needed = ("N", "CA", "C") if kind == ACE else ("N", "CA", "C", "O")
    if any(a not in index for a in needed) or (kind == NME and "OXT" in index):
        return None
    end = index["N"] if kind == ACE else index["C"]
    partners = s.bonds[(s.bonds == end).any(axis=1)].ravel()
    if (s.residue_of[partners] != r).any():
        return None  # bonded onward already (a cyclic peptide, an isopeptide)
    if kind == NME:
        x = s.coords.astype(np.float64)
        ca, c = x[:, index["CA"]], x[:, index["C"]]
        nn = place(x[:, index["O"]], ca, c, PEPTIDE, np.radians(116.2), np.pi)
        methyl = place(ca, c, nn, 1.458, np.radians(121.7), np.pi)
        names, elements, xyz = ["N", "C"], [7, 6], np.stack([nn, methyl], axis=1)
        bonds, order = [[end, s.n_atoms], [s.n_atoms, s.n_atoms + 1]], [1, 1]
        at, seq = r + 1, int(s.res_seq[r]) + 1
    else:
        xyz = _ace(s, index)
        names, elements = ["C", "O", "CH3"], [6, 8, 6]
        bonds = [[s.n_atoms, end], [s.n_atoms, s.n_atoms + 1], [s.n_atoms, s.n_atoms + 2]]
        order = [1, 2, 1]
        at, seq = r, int(s.res_seq[r]) - 1
    chain = int(s.chain_of[r])
    icode = _free_icode(s, chain, seq)
    return edit.insert_residue(s, at, chain, kind, seq, np.array(names), np.array(elements), xyz,
                               np.array(bonds), np.array(order), icode=icode)  # fmt: skip


def _ace(s: Structure, index: dict[str, int]) -> np.ndarray:
    """An ACE's C, O and CH3 (k, 3, 3) on the N of a residue (its heavy atoms by name in
    *index*), at the φ (every 30°) whose atoms come nearest no atom but that N and CA, in the
    first model."""
    x = s.coords.astype(np.float64)
    n, ca, c = (x[:, index[a]] for a in ("N", "CA", "C"))
    others = np.ones(s.n_atoms, dtype=bool)
    others[[index["N"], index["CA"]]] = False
    tree = cKDTree(x[0, others])
    options = []
    for phi in np.radians(np.arange(-180, 180, 30)):
        carbon = place(c, ca, n, PEPTIDE, np.radians(121.7), phi)
        oxygen = place(ca, n, carbon, 1.229, np.radians(122.7), 0.0)
        methyl = place(ca, n, carbon, 1.52, np.radians(116.5), np.pi)
        xyz = np.stack([carbon, oxygen, methyl], axis=1)
        options.append((float(tree.query(xyz[0])[0].min()), xyz))
    return max(options, key=lambda option: option[0])[1]


def _free_icode(s: Structure, chain: int, seq: int) -> str:
    """No insertion code when *seq* is free in *chain*, else the first letter that is."""
    span = s.residues_of(chain)
    taken = {str(s.icode[r]) for r in range(span.start, span.stop) if s.res_seq[r] == seq}
    return next(code for code in ["", *"ABCDEFGHIJKLMNOPQRSTUVWXYZ"] if code not in taken)


def _build(t: Template, at: dict[str, np.ndarray], lacking: list[str]) -> dict[str, np.ndarray]:
    """Positions (every model) for the atoms *lacking*, each from three placed atoms."""
    placed = dict(at)
    built: dict[str, np.ndarray] = {}
    todo = [a for a in lacking if a in t.index]
    while todo:
        progress = False
        for name in list(todo):
            refs = _references(t, name, placed)
            if refs is None:
                continue
            parent, angle, dihedral = refs
            p = t.placement(name, parent, angle, dihedral)
            position = place(placed[dihedral], placed[angle], placed[parent], p.length, p.angle,
                             p.dihedral)  # fmt: skip
            placed[name] = built[name] = position
            todo.remove(name)
            progress = True
        if not progress:
            break
    return built


def _references(t: Template, name: str, placed: dict) -> tuple[str, str, str] | None:
    """A placed parent of *name*, a placed neighbour of the parent, and a third placed atom
    bonded to either (heavy atoms only)."""
    i = t.index[name]
    for p in t.heavy_neighbours(i):
        parent = t.atoms[p]
        if parent not in placed:
            continue
        for a in t.heavy_neighbours(p):
            angle = t.atoms[a]
            if a == i or angle not in placed:
                continue
            thirds = [t.atoms[d] for d in (*t.heavy_neighbours(a), *t.heavy_neighbours(p))
                      if d not in (i, p, a) and t.atoms[d] in placed]  # fmt: skip
            if thirds:
                return parent, angle, thirds[0]
    return None
