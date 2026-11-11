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
"""Bonds of a structure.

- **Within a residue**, from its template (``data/residues.json``: the bonds of the standard
  residues and common cofactors, from the wwPDB Chemical Component Dictionary), with their
  orders. Residues without a template (most ligands) and atoms a template does not name (a
  hydrogen named by other software) are bonded by distance: closer than the sum of the
  covalent radii plus 0.45 Å.
- **Between residues**: the peptide bond C–N and the nucleic-acid link O3'–P of neighbours in
  a chain, disulfides SG–SG, and the covalent links the file declares (``struct_conn``), which
  is how a covalent inhibitor stays attached.

Two atoms in different alternate locations of the same residue are never bonded.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from functools import cache
from importlib import resources

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from .model import AMINO, AROMATIC, NUCLEIC

#: added to the two covalent radii for a distance bond
TOLERANCE = 0.45
#: closer than this is not a bond but a clash or a duplicate atom
TOO_CLOSE = 0.4
#: longest peptide, nucleic-acid or disulfide bond accepted between residues
LINK = 2.0
DISULFIDE = 2.5
#: and of a link the file declares (it may be a little stretched)
DECLARED = 3.5


@cache
def templates() -> dict[str, list[tuple[str, str, int]]]:
    """Residue name -> bonds ``(atom, atom, order)`` with AROMATIC as the order of aromatic
    bonds."""
    text = resources.files("chemsuite.protein.structure").joinpath("data/residues.json").read_text()
    out: dict[str, list[tuple[str, str, int]]] = {}
    for name, entry in json.loads(text).items():
        out[name] = [(a, b, AROMATIC if aromatic else order) for a, b, order, aromatic
                     in entry["bonds"]]  # fmt: skip
    return out


@cache
def template_smiles(name: str) -> str:
    text = resources.files("chemsuite.protein.structure").joinpath("data/residues.json").read_text()
    entry = json.loads(text).get(name)
    return entry["smiles"] if entry else ""


@cache
def _covalent_radii() -> np.ndarray:
    table = Chem.GetPeriodicTable()
    radii = np.full(119, 0.77, dtype=np.float32)
    for z in range(1, 119):
        radii[z] = table.GetRcovalent(z)
    radii[0] = 0.77
    return radii


def compatible(altloc: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Pairs that can be bonded: same alternate location, or at least one has none."""
    a, b = altloc[i], altloc[j]
    return (a == "") | (b == "") | (a == b)


def perceive(
    names: np.ndarray,
    elements: np.ndarray,
    altloc: np.ndarray,
    xyz: np.ndarray,
    res_start: np.ndarray,
    res_name: np.ndarray,
    res_kind: np.ndarray,
    chain_start: np.ndarray,
    links: Iterable[tuple[int, int]] = (),
    library: Mapping[str, list[tuple[str, str, int]]] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Bonds ``(m, 2)`` int32 with ``i < j`` and their orders ``(m,)`` int8. ``library`` adds
    residue templates (the file's own components) to the bundled ones."""
    pairs: list[np.ndarray] = []
    orders: list[np.ndarray] = []
    templated = np.zeros(len(names), dtype=bool)
    library = {**templates(), **(library or {})}
    for residue in range(len(res_name)):
        bonds = library.get(str(res_name[residue]))
        if not bonds:
            continue
        start, end = int(res_start[residue]), int(res_start[residue + 1])
        found = _template_bonds(bonds, names[start:end], altloc[start:end])
        if found is not None:
            pairs.append(found[0] + start)
            orders.append(found[1])
            templated[start:end][np.unique(found[0])] = True
    residue_of = np.repeat(np.arange(len(res_name)), np.diff(res_start))
    distance = _distance_bonds(elements, altloc, xyz, residue_of, templated)
    pairs.append(distance)
    orders.append(np.ones(len(distance), dtype=np.int8))
    between = _links(names, altloc, xyz, res_start, res_kind, chain_start, residue_of, links)
    pairs.append(between)
    orders.append(np.ones(len(between), dtype=np.int8))
    return _unique(np.concatenate(pairs).reshape(-1, 2), np.concatenate(orders))


def _template_bonds(
    bonds: list[tuple[str, str, int]], names: np.ndarray, altloc: np.ndarray
) -> tuple[np.ndarray, np.ndarray] | None:
    where: dict[str, list[int]] = {}
    for i, name in enumerate(names):
        where.setdefault(str(name), []).append(i)
    pairs: list[tuple[int, int]] = []
    orders: list[int] = []
    for a, b, order in bonds:
        for i in where.get(a, ()):
            for j in where.get(b, ()):
                if altloc[i] == "" or altloc[j] == "" or altloc[i] == altloc[j]:
                    pairs.append((i, j))
                    orders.append(order)
    if not pairs:
        return None
    return np.array(pairs, dtype=np.int64), np.array(orders, dtype=np.int8)


def _distance_bonds(
    elements: np.ndarray,
    altloc: np.ndarray,
    xyz: np.ndarray,
    residue_of: np.ndarray,
    templated: np.ndarray,
) -> np.ndarray:
    """Bonds by distance inside each residue, for atoms its template did not bond."""
    loose = np.flatnonzero(~templated)
    if not len(loose):
        return np.zeros((0, 2), dtype=np.int64)
    # pairs with at least one loose atom, both in the same residue
    candidates = np.flatnonzero(np.isin(residue_of, np.unique(residue_of[loose])))
    radii = _covalent_radii()[elements[candidates]]
    tree = cKDTree(xyz[candidates])
    reach = 2 * float(radii.max()) + TOLERANCE
    pairs = tree.query_pairs(reach, output_type="ndarray")
    if not len(pairs):
        return np.zeros((0, 2), dtype=np.int64)
    i, j = candidates[pairs[:, 0]], candidates[pairs[:, 1]]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    limit = radii[pairs[:, 0]] + radii[pairs[:, 1]] + TOLERANCE
    keep = (
        (residue_of[i] == residue_of[j])
        & ~(templated[i] & templated[j])
        & (d > TOO_CLOSE)
        & (d < limit)
        & compatible(altloc, i, j)
        & ~((elements[i] == 1) & (elements[j] == 1))
    )
    return np.stack([i[keep], j[keep]], axis=1)


def _links(
    names: np.ndarray,
    altloc: np.ndarray,
    xyz: np.ndarray,
    res_start: np.ndarray,
    res_kind: np.ndarray,
    chain_start: np.ndarray,
    residue_of: np.ndarray,
    declared: Iterable[tuple[int, int]],
) -> np.ndarray:
    found = [np.array(list(declared), dtype=np.int64).reshape(-1, 2)]
    limits = [DECLARED]
    chain_of_residue = np.repeat(np.arange(len(chain_start) - 1), np.diff(chain_start))
    for kind, head, tail in ((AMINO, "C", "N"), (NUCLEIC, "O3'", "P")):
        first = _atoms_named(names, head, residue_of, res_kind == kind)
        second = _atoms_named(names, tail, residue_of, res_kind == kind)
        if not len(first) or not len(second):
            continue
        # a head of residue r with a tail of residue r + 1 in the same chain
        by_residue: dict[int, list[int]] = {}
        for atom in second:
            by_residue.setdefault(int(residue_of[atom]), []).append(int(atom))
        pairs = [
            (int(a), b)
            for a in first
            for b in by_residue.get(int(residue_of[a]) + 1, ())
            if chain_of_residue[residue_of[a]] == chain_of_residue[residue_of[b]]
        ]
        found.append(np.array(pairs, dtype=np.int64).reshape(-1, 2))
        limits.append(LINK)
    sulfur = np.flatnonzero(names == "SG")
    if len(sulfur) > 1:
        pairs = cKDTree(xyz[sulfur]).query_pairs(DISULFIDE, output_type="ndarray")
        found.append(sulfur[pairs].reshape(-1, 2))
        limits.append(DISULFIDE)
    pairs = np.concatenate(found)
    if not len(pairs):
        return pairs
    limit = np.repeat(limits, [len(f) for f in found])
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    return pairs[(d < limit) & (d > TOO_CLOSE) & compatible(altloc, i, j)]


def _atoms_named(
    names: np.ndarray, name: str, residue_of: np.ndarray, residue_ok: np.ndarray
) -> np.ndarray:
    atoms = np.flatnonzero(names == name)
    return atoms[residue_ok[residue_of[atoms]]]


def _unique(pairs: np.ndarray, orders: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if not len(pairs):
        return np.zeros((0, 2), np.int32), np.zeros(0, np.int8)
    pairs = np.sort(pairs, axis=1)
    # sorted by atom, so the same atoms give the same bonds whatever found them first
    unique, first = np.unique(pairs, axis=0, return_index=True)
    return unique.astype(np.int32), orders[first].astype(np.int8)
