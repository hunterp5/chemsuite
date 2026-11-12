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
  checked for clashes. Residues missing N, CA or C are left as they are; no loops are built.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from ..structure import edit
from ..structure.model import Structure
from .hydrogens import HYDROGEN
from .inspect import BACKBONE, POLYMER, last_in_chain, missing_atoms
from .templates import Library, Template, place


def one_location(s: Structure, altlocs: dict[int, str] | None = None) -> Structure:
    """*s* with only the atoms shown, none of them marked as an alternate location."""
    if (s.altloc == "").all():
        return s
    kept = s.subset(s.active(altlocs))
    return replace(kept, altloc=np.full(kept.n_atoms, "", dtype=kept.altloc.dtype),
                   occupancy=np.ones(kept.n_atoms, dtype=kept.occupancy.dtype))  # fmt: skip


def build_missing(s: Structure, library: Library, cancel: CancelToken | None = None) -> Structure:
    """*s* with the missing heavy atoms of its polymer residues built (uids of the old atoms
    kept)."""
    ends = last_in_chain(s)
    x = s.coords.astype(np.float64)
    residues, names, elements, coords, bonds, orders = [], [], [], [], [], []
    n = s.n_atoms
    for r in range(s.n_residues):
        if s.res_kind[r] not in POLYMER:
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
