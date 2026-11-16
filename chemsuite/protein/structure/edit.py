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

"""Edits to a structure: add or delete atoms, change an element or a charge, rename and
renumber.

Each edit returns a new :class:`~.model.Structure` (structures are never changed in place) and
keeps every atom's ``uid``, so selections, labels and measurements follow the atoms. An edit
that cannot be made raises :class:`EditError`, whose message says why.
"""

from __future__ import annotations

import re
from dataclasses import replace

import gemmi
import numpy as np

from .classify import residue_kind
from .model import AMINO, NUCLEIC, Structure

_NAME = re.compile(r"[A-Z0-9]{1,5}")
_CHAIN = re.compile(r"[A-Za-z0-9]{1,4}")
MAX_CHARGE = 8


class EditError(ValueError):
    """The edit cannot be made; the message says why."""


def delete(s: Structure, atoms: np.ndarray) -> Structure:
    """Without the atoms in the mask *atoms*. An atom goes with all its alternate locations
    (the same name in the same residue), so no hidden location is left behind."""
    gone = with_alternates(s, atoms)
    if not gone.any():
        return s
    if gone.all():
        raise EditError("Deleting every atom would leave nothing: use Remove Protein")
    return s.subset(~gone)


def with_alternates(s: Structure, atoms: np.ndarray) -> np.ndarray:
    """The mask *atoms* widened to the other locations of each atom it holds."""
    atoms = np.asarray(atoms, dtype=bool).copy()
    alternate = s.altloc != ""
    chosen = atoms & alternate
    if chosen.any():
        keys = np.char.add(s.residue_of.astype("<U10"), np.char.add(":", s.name))
        atoms |= alternate & np.isin(keys, keys[chosen])
    return atoms


def set_element(s: Structure, atoms: np.ndarray, symbol: str) -> Structure:
    """The atoms in the mask *atoms* made *symbol* (``Cl``, ``N``…); bonds are kept."""
    element = gemmi.Element(symbol.strip().capitalize())
    if element.atomic_number == 0:
        raise EditError(f"'{symbol}' is not an element")
    values = s.element.copy()
    values[np.asarray(atoms, dtype=bool)] = element.atomic_number
    return replace(s, element=values)


def set_charge(s: Structure, atoms: np.ndarray, charge: int) -> Structure:
    """The atoms in the mask *atoms* given the formal charge *charge*."""
    if not -MAX_CHARGE <= charge <= MAX_CHARGE:
        raise EditError(f"A formal charge is between -{MAX_CHARGE} and +{MAX_CHARGE}")
    values = s.charge.copy()
    values[np.asarray(atoms, dtype=bool)] = charge
    return replace(s, charge=values)


def rename_residues(s: Structure, residues: np.ndarray, name: str) -> Structure:
    """The residues *residues* (indices) renamed *name*. A residue outside a polymer takes
    the kind its new name gives it: a water renamed ``NA`` becomes an ion."""
    name = name.strip().upper()
    if not _NAME.fullmatch(name):
        raise EditError("A residue name is 1 to 5 letters or digits")
    residues = np.asarray(residues, dtype=np.int64)
    names = s.res_name.copy()
    names[residues] = name
    kinds = s.res_kind.copy()
    for residue in residues:
        if kinds[residue] in (AMINO, NUCLEIC):
            continue
        span = s.atoms_of(int(residue))
        heavy = tuple(
            gemmi.Element(int(z)).name.upper() for z in s.element[span] if z not in (0, 1)
        )
        kinds[residue] = residue_kind(name, gemmi.EntityType.NonPolymer, heavy)
    return replace(s, res_name=names, res_kind=kinds)


def add_atoms(
    s: Structure,
    residues: np.ndarray,
    names: np.ndarray,
    elements: np.ndarray,
    xyz: np.ndarray,
    bonds: np.ndarray | None = None,
    order: np.ndarray | None = None,
    charge: np.ndarray | None = None,
    b: np.ndarray | None = None,
) -> Structure:
    """New atoms added to existing residues, each at the end of its residue's atoms.

    Atom *j* goes into residue ``residues[j]`` as ``names[j]`` (element ``elements[j]``, an
    atomic number) at ``xyz[j]``, or ``xyz[k, j]`` in model *k*. *bonds* are pairs of atom
    indices, where ``n_atoms + j`` is new atom *j*; *order* is each bond's order (single when
    left out). The new atoms take uids from ``next_uid`` up, in the order given; every old atom
    keeps its uid."""
    residues = np.asarray(residues, dtype=np.int64)
    m = len(residues)
    if m == 0:
        return s
    if residues.min() < 0 or residues.max() >= s.n_residues:
        raise EditError("Atoms can only be added to a residue of the structure")
    xyz = np.asarray(xyz, dtype=np.float32)
    if xyz.ndim == 2:
        xyz = np.broadcast_to(xyz, (s.n_models, m, 3))
    if xyz.shape != (s.n_models, m, 3):
        raise EditError("The new atoms need a position in every model")
    n = s.n_atoms
    added = np.bincount(residues, minlength=s.n_residues)
    before = np.concatenate([[0], np.cumsum(added)])  # atoms added ahead of each residue
    res_start = (s.res_start + before).astype(np.int32)
    # where every old and every new atom lands
    old_at = np.arange(n) + before[s.residue_of]
    sorted_new = np.argsort(residues, kind="stable")
    rank = np.empty(m, dtype=np.int64)
    rank[sorted_new] = np.arange(m) - before[residues[sorted_new]]
    new_at = res_start[residues + 1] - added[residues] + rank
    where = np.concatenate([old_at, new_at])

    def grown(old: np.ndarray, new: np.ndarray, axis: int = 0) -> np.ndarray:
        both = np.concatenate([old, np.asarray(new, dtype=old.dtype)], axis=axis)
        out = np.empty_like(both)
        if axis == 0:
            out[where] = both
        else:
            out[:, where] = both
        return out

    zeros = np.zeros(m)
    links = np.asarray(bonds if bonds is not None else np.zeros((0, 2)), dtype=np.int64)
    links = where[links.reshape(-1, 2)]
    kinds = np.asarray(order if order is not None else np.ones(len(links)), dtype=np.int8)
    all_bonds = np.concatenate([where[s.bonds.astype(np.int64)], links]) if len(s.bonds) else links
    all_bonds = np.sort(all_bonds, axis=1).astype(np.int32)
    return replace(
        s,
        uid=grown(s.uid, s.next_uid + np.arange(m)),
        coords=grown(s.coords, xyz, axis=1),
        element=grown(s.element, elements),
        name=grown(s.name, names),
        altloc=grown(s.altloc, np.full(m, "")),
        occupancy=grown(s.occupancy, np.ones(m)),
        b=grown(s.b, b if b is not None else zeros),
        charge=grown(s.charge, charge if charge is not None else zeros),
        res_start=res_start,
        bonds=all_bonds.reshape(-1, 2),
        order=np.concatenate([s.order, kinds]),
        next_uid=s.next_uid + m,
    )


def insert_residue(
    s: Structure,
    at: int,
    chain: int,
    name: str,
    seq: int,
    names: np.ndarray,
    elements: np.ndarray,
    xyz: np.ndarray,
    bonds: np.ndarray | None = None,
    order: np.ndarray | None = None,
    icode: str = "",
    seqres: int = -1,
    entity: int = -1,
) -> Structure:
    """A new residue *name* numbered *seq* (*icode*) put before residue *at* (an index;
    ``n_residues`` for the very end), in chain *chain*, which must hold residue *at* or end
    there. Its atoms are *names*, *elements* and *xyz* as for :func:`add_atoms`, as are *bonds*
    (``n_atoms + j`` is new atom *j*) and *order*. Its kind follows from its name as in a
    polymer; it is at *seqres* in *entity*'s sequence (-1: in none). The new atoms take uids
    from ``next_uid`` up."""
    first, last = int(s.chain_start[chain]), int(s.chain_start[chain + 1])
    if not first <= at <= last:
        raise EditError("A residue goes inside its chain")
    if not _NAME.fullmatch(name):
        raise EditError("A residue name is 1 to 5 letters or digits")
    m = len(names)
    xyz = np.asarray(xyz, dtype=np.float32)
    if xyz.ndim == 2:
        xyz = np.broadcast_to(xyz, (s.n_models, m, 3))
    if m == 0 or xyz.shape != (s.n_models, m, 3):
        raise EditError("The new atoms need a position in every model")
    p = int(s.res_start[at])  # where the new atoms go
    n = s.n_atoms
    old_to_new = np.arange(n) + (np.arange(n) >= p) * m
    where = np.concatenate([old_to_new, p + np.arange(m)])
    links = np.asarray(bonds if bonds is not None else np.zeros((0, 2)), dtype=np.int64)
    links = where[links.reshape(-1, 2)]
    kinds = np.asarray(order if order is not None else np.ones(len(links)), dtype=np.int8)
    old_bonds = old_to_new[s.bonds.astype(np.int64)] if len(s.bonds) else np.zeros((0, 2))
    all_bonds = np.sort(np.concatenate([old_bonds, links]), axis=1).astype(np.int32)
    heavy = tuple(gemmi.Element(int(z)).name.upper() for z in elements if z not in (0, 1))
    kind = residue_kind(name, gemmi.EntityType.Polymer, heavy)
    neighbour = at if at < last else at - 1  # a residue of the chain, for its subchain
    res_start = np.concatenate([s.res_start[: at + 1], s.res_start[at:] + m]).astype(np.int32)
    chain_start = s.chain_start.copy()
    chain_start[chain + 1 :] += 1

    def put(values: np.ndarray, new: object) -> np.ndarray:
        return np.insert(values, at, np.asarray(new, dtype=values.dtype))

    def atoms(values: np.ndarray, new: object, axis: int = 0) -> np.ndarray:
        return np.insert(values, [p] * m, np.asarray(new, dtype=values.dtype), axis=axis)

    return replace(
        s,
        uid=atoms(s.uid, s.next_uid + np.arange(m)),
        coords=np.concatenate([s.coords[:, :p], xyz, s.coords[:, p:]], axis=1),
        element=atoms(s.element, elements),
        name=atoms(s.name, names),
        altloc=atoms(s.altloc, np.full(m, "")),
        occupancy=atoms(s.occupancy, np.ones(m)),
        b=atoms(s.b, np.zeros(m)),
        charge=atoms(s.charge, np.zeros(m)),
        res_start=res_start,
        res_name=put(s.res_name, name),
        res_seq=put(s.res_seq, seq),
        icode=put(s.icode, icode),
        res_kind=put(s.res_kind, kind),
        ss=put(s.ss, 0),
        seqres=put(s.seqres, seqres),
        label_asym=put(s.label_asym, s.label_asym[neighbour] if s.n_residues else ""),
        entity=put(s.entity, entity),
        chain_start=chain_start,
        bonds=all_bonds.reshape(-1, 2),
        order=np.concatenate([s.order, kinds]),
        next_uid=s.next_uid + m,
    )


def rename_chain(s: Structure, chain: int, name: str) -> Structure:
    """Chain *chain* (an index) renamed *name*, which no other chain may have."""
    name = name.strip()
    if not _CHAIN.fullmatch(name):
        raise EditError("A chain name is 1 to 4 letters or digits")
    others = np.delete(s.chain_name, chain)
    if name in others:
        raise EditError(f"There is already a chain {name}")
    names = s.chain_name.copy()
    names[chain] = name
    return replace(s, chain_name=names)


def renumber(s: Structure, residues: np.ndarray, first: int) -> Structure:
    """The residues *residues* (indices, in one chain) numbered from *first* up, in their
    order, without insertion codes. Refused when a new number is already another residue's."""
    residues = np.sort(np.asarray(residues, dtype=np.int64))
    if not len(residues):
        return s
    chains = np.unique(s.chain_of[residues])
    if len(chains) != 1:
        raise EditError("Renumber residues of one chain at a time")
    numbers = np.arange(first, first + len(residues), dtype=np.int64)
    if numbers[-1] > np.iinfo(np.int32).max or numbers[0] < np.iinfo(np.int32).min:
        raise EditError("Those numbers are out of range")
    span = s.residues_of(int(chains[0]))
    rest = np.setdiff1d(np.arange(span.start, span.stop), residues)
    taken = rest[(s.icode[rest] == "") & np.isin(s.res_seq[rest], numbers)]
    if len(taken):
        raise EditError(f"{s.residue_label(int(taken[0]))} already has that number")
    seq = s.res_seq.copy()
    seq[residues] = numbers
    icode = s.icode.copy()
    icode[residues] = ""
    return replace(s, res_seq=seq, icode=icode)
