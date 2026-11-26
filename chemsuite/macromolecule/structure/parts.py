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

"""Parts of a structure: which atoms are left out of calculations, and whole residues copied
from one structure into another.

A structure's excluded atoms are kept beside it (in its asset's meta, under :data:`EXCLUDED`)
as uids, so they survive edits and leave the structure, and so its preparation, as it was. A
calculation runs on :func:`included` and :func:`restore` puts the parts it left out back
unchanged. :func:`copy_residues` duplicates a ligand, the waters or a chain.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace

import numpy as np

from .edit import EditError
from .model import AMINO, KIND_NAMES, LIGAND, NUCLEIC, Structure

#: the asset meta key holding the excluded atoms' uids (:func:`write_excluded`)
EXCLUDED = "excluded"


# -- excluded atoms -------------------------------------------------------------------------------
def write_excluded(uids: np.ndarray) -> str:
    """*uids* as runs: ``"812-1040 1502"`` (``""`` for none)."""
    values = np.unique(np.asarray(uids, dtype=np.int64))
    if not len(values):
        return ""
    breaks = np.flatnonzero(np.diff(values) != 1) + 1
    runs = []
    for run in np.split(values, breaks):
        first, last = int(run[0]), int(run[-1])
        runs.append(str(first) if first == last else f"{first}-{last}")
    return " ".join(runs)


def read_uids(text: str) -> np.ndarray:
    """The uids :func:`write_excluded` wrote (anything unreadable is passed over)."""
    found: list[np.ndarray] = []
    for word in text.split():
        first, _, last = word.partition("-")
        if not first.isdigit() or (last and not last.isdigit()):
            continue
        a, b = int(first), int(last or first)
        if a <= b:
            found.append(np.arange(a, b + 1, dtype=np.int64))
    return np.concatenate(found) if found else np.zeros(0, dtype=np.int64)


def excluded(meta: Mapping[str, str], s: Structure) -> np.ndarray:
    """The atoms of *s* left out of calculations (bool per atom); uids *s* no longer has are
    passed over."""
    text = meta.get(EXCLUDED, "")
    if not text:
        return np.zeros(s.n_atoms, dtype=bool)
    return np.isin(s.uid, read_uids(text))


def with_excluded(meta: Mapping[str, str], s: Structure, mask: np.ndarray) -> dict[str, str]:
    """*meta* with *s*'s atoms *mask* the excluded ones."""
    out = {k: v for k, v in meta.items() if k != EXCLUDED}
    if mask.any():
        out[EXCLUDED] = write_excluded(s.uid[mask])
    return out


def included(s: Structure, left_out: np.ndarray) -> Structure:
    """*s* without the atoms *left_out*."""
    return s.subset(~np.asarray(left_out, dtype=bool))


def whole_residues(s: Structure, atoms: np.ndarray) -> np.ndarray:
    """*atoms* grown to the whole residues they touch."""
    out = np.zeros(s.n_atoms, dtype=bool)
    if atoms.any():
        out = np.isin(s.residue_of, np.unique(s.residue_of[atoms]))
    return out


def restore(done: Structure, original: Structure, left_out: np.ndarray) -> Structure:
    """*done* (a calculation's result on :func:`included`) with *original*'s atoms *left_out*
    put back as they were, uids kept: each residue at the end of its chain, a chain left out
    whole as a chain of its own, and *original*'s bonds between them and *done*'s atoms."""
    residues = np.unique(original.residue_of[np.asarray(left_out, dtype=bool)])
    if not len(residues):
        return done
    taken = set(done.uid.tolist())
    if taken & set(original.uid[np.isin(original.residue_of, residues)].tolist()):
        raise EditError("The parts left out are already in the structure")
    return _graft(done, original, residues, keep_uids=True)


def copy_residues(s: Structure, residues: np.ndarray, chain_name: str = "") -> Structure:
    """*s* with copies of *residues* (new uids, the same positions): each at the end of its
    chain, numbered on from the chain's last, or all in a new chain *chain_name*."""
    residues = np.unique(np.asarray(residues, dtype=np.int64))
    if not len(residues):
        return s
    if chain_name:
        segments = s.chain_segment[s.chain_of[residues]]
        if (s.chain_name[s.chain_segment == segments[0]] == chain_name).any():
            raise EditError(f"There is already a chain {chain_name}")
    return _graft(s, s, residues, keep_uids=False, chain_name=chain_name)


def free_chain_name(s: Structure) -> str:
    """The first chain name *s* does not use: A–Z, then a–z, then 0–9, then two letters."""
    used = set(s.chain_name.tolist())
    letters = [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    for name in [*letters, *(c.lower() for c in letters), *map(str, range(10))]:
        if name not in used:
            return name
    return next(a + b for a in letters for b in letters if a + b not in used)


# -- the grafting ---------------------------------------------------------------------------------
def _graft(
    s: Structure, source: Structure, residues: np.ndarray, *, keep_uids: bool, chain_name: str = ""
) -> Structure:
    """*s* with *source*'s *residues* added (see :func:`restore`, :func:`copy_residues`)."""
    piece = source.subset(np.isin(source.residue_of, residues))
    n, m = s.n_atoms, piece.n_atoms
    # the chain of s each piece residue joins (-1: a new chain, numbered after s's)
    where = {(str(s.chain_name[c]), str(s.chain_segment[c])): c for c in range(s.n_chains)}
    piece_chain = piece.chain_of
    target = np.empty(piece.n_residues, dtype=np.int64)
    new_chains: list[tuple[str, str]] = []
    for c in range(piece.n_chains):
        key = (str(piece.chain_name[c]), str(piece.chain_segment[c]))
        if chain_name:
            key = (chain_name, key[1])
        if key in where and not chain_name:
            target[piece_chain == c] = where[key]
        else:
            if key not in new_chains:
                new_chains.append(key)
            target[piece_chain == c] = s.n_chains + new_chains.index(key)
    # numbers: a copy into its own chain goes on from the chain's last
    seq = piece.res_seq.astype(np.int64)
    icode = piece.icode
    if not keep_uids and not chain_name:
        seq = seq.copy()
        for c in np.unique(target).tolist():
            span = s.residues_of(c)
            last = int(s.res_seq[span].max()) if span.stop > span.start else 0
            mine = np.flatnonzero(target == c)
            seq[mine] = last + 1 + np.arange(len(mine))
        icode = np.full(piece.n_residues, "", dtype=piece.icode.dtype)
    # residues in their new order: by chain (s's, then the new ones); within a chain a copy goes
    # last, a residue put back right after the residue it followed in *source*
    old_chain = s.chain_of.astype(np.int64)
    order_chain = np.concatenate([old_chain, target])
    place = np.concatenate([np.arange(s.n_residues, dtype=np.float64),
                            _places(s, source, residues, target, keep_uids)])  # fmt: skip
    order = np.lexsort((place, order_chain))
    counts = np.concatenate([np.diff(s.res_start), np.diff(piece.res_start)])
    atom_start = np.concatenate([s.res_start[:-1], n + piece.res_start[:-1]])
    atoms = np.concatenate([np.arange(atom_start[r], atom_start[r] + counts[r]) for r in order])
    new_index = np.empty(n + m, dtype=np.int64)
    new_index[atoms] = np.arange(n + m)
    res_start = np.concatenate([[0], np.cumsum(counts[order])]).astype(np.int32)
    chain_count = np.bincount(order_chain, minlength=s.n_chains + len(new_chains))
    chain_start = np.concatenate([[0], np.cumsum(chain_count)]).astype(np.int32)

    uid = piece.uid if keep_uids else s.next_uid + np.arange(m, dtype=np.int64)
    next_uid = max(s.next_uid, source.next_uid) if keep_uids else s.next_uid + m
    entity, seqres = _entities(s, piece)

    def atom(a: np.ndarray, b: np.ndarray, axis: int = 0) -> np.ndarray:
        both = np.concatenate([a, b.astype(a.dtype) if a.dtype.kind != "U" else b], axis=axis)
        return np.take(both, atoms, axis=axis)

    def residue(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        both = np.concatenate([a, b.astype(a.dtype) if a.dtype.kind != "U" else b])
        return both[order]

    names = np.array([k[0] for k in new_chains], dtype=str)
    segments = np.array([k[1] for k in new_chains], dtype=str)
    bonds, bond_order = _bonds(s, source, piece, keep_uids)
    bonds = np.sort(new_index[bonds], axis=1).astype(np.int32) if len(bonds) else bonds
    return replace(
        s,
        uid=atom(s.uid, uid),
        coords=atom(s.coords, _models(piece.coords, s.n_models), axis=1),
        element=atom(s.element, piece.element),
        name=atom(s.name, piece.name),
        altloc=atom(s.altloc, piece.altloc),
        occupancy=atom(s.occupancy, piece.occupancy),
        b=atom(s.b, piece.b),
        charge=atom(s.charge, piece.charge),
        res_start=res_start,
        res_name=residue(s.res_name, piece.res_name),
        res_seq=residue(s.res_seq, seq),
        icode=residue(s.icode, icode),
        res_kind=residue(s.res_kind, piece.res_kind),
        ss=residue(s.ss, piece.ss),
        seqres=residue(s.seqres, seqres),
        label_asym=residue(s.label_asym, piece.label_asym),
        entity=residue(s.entity, entity),
        chain_start=chain_start,
        chain_name=np.concatenate([s.chain_name, names]) if len(names) else s.chain_name,
        chain_segment=np.concatenate([s.chain_segment, segments])
        if len(names)
        else s.chain_segment,
        bonds=bonds.reshape(-1, 2),
        order=bond_order,
        next_uid=int(next_uid),
    )


def _places(
    s: Structure, source: Structure, residues: np.ndarray, target: np.ndarray, keep_uids: bool
) -> np.ndarray:
    """Where each added residue goes among *s*'s (residue positions, a fraction between):
    after the residue of its chain it followed in *source* that *s* still has, else first in
    the chain; last in the chain for a copy."""
    k = len(residues)
    tie = (np.arange(k) + 1) / (k + 1)  # keeps added residues in their order, between two of s
    if not keep_uids:
        return s.n_residues + tie
    # each source residue's place in s, by its first atom (-1: not in s)
    found = s.indices(source.uid[source.res_start[:-1]])
    present = np.full(len(found), -1, dtype=np.int64)
    present[found >= 0] = s.residue_of[found[found >= 0]]
    out = np.empty(k)
    for j, r in enumerate(residues.tolist()):
        chain = int(target[j])
        anchor = -1.0
        if chain < s.n_chains:
            span = s.residues_of(chain)
            anchor = span.start - 1.0
            p = r - 1
            while p >= 0 and source.chain_of[p] == source.chain_of[r]:
                if span.start <= present[p] < span.stop:
                    anchor = float(present[p])
                    break
                p -= 1
        out[j] = anchor + tie[j]
    return out


def _models(coords: np.ndarray, k: int) -> np.ndarray:
    """*coords* in *k* models: the first one repeated when they have fewer."""
    if len(coords) == k:
        return coords
    return np.repeat(coords[:1], k, axis=0)


def _entities(s: Structure, piece: Structure) -> tuple[np.ndarray, np.ndarray]:
    """The piece's entities as *s* numbers them (by entity id; -1 and no sequence position
    where *s* has no such entity)."""
    ids = {str(e): i for i, e in enumerate(s.entity_id.tolist())}
    entity = np.full(piece.n_residues, -1, dtype=s.entity.dtype)
    seqres = np.full(piece.n_residues, -1, dtype=s.seqres.dtype)
    for r in range(piece.n_residues):
        e = int(piece.entity[r])
        if 0 <= e < len(piece.entity_id) and str(piece.entity_id[e]) in ids:
            entity[r] = ids[str(piece.entity_id[e])]
            seqres[r] = piece.seqres[r]
    return entity, seqres


def _bonds(
    s: Structure, source: Structure, piece: Structure, keep_uids: bool
) -> tuple[np.ndarray, np.ndarray]:
    """The bonds of *s* and of the piece (atoms numbered *s*'s then the piece's), and, with
    uids kept, *source*'s bonds between a piece atom and an atom of *s*."""
    n = s.n_atoms
    bonds = [s.bonds.astype(np.int64).reshape(-1, 2), piece.bonds.astype(np.int64) + n]
    orders = [s.order, piece.order]
    if keep_uids and len(source.bonds):
        in_piece = np.isin(source.uid, piece.uid)
        crossing = in_piece[source.bonds].sum(axis=1) == 1
        here = {int(u): i for i, u in enumerate(s.uid.tolist())}
        there = {int(u): n + i for i, u in enumerate(piece.uid.tolist())}
        extra, kinds = [], []
        for (a, b), kind in zip(source.bonds[crossing].tolist(), source.order[crossing].tolist(),
                                strict=True):  # fmt: skip
            ua, ub = int(source.uid[a]), int(source.uid[b])
            i = there.get(ua, here.get(ua))
            j = there.get(ub, here.get(ub))
            if i is not None and j is not None:
                extra.append((i, j))
                kinds.append(kind)
        if extra:
            bonds.append(np.array(extra, dtype=np.int64))
            orders.append(np.array(kinds, dtype=np.int8))
    return np.concatenate(bonds), np.concatenate(orders).astype(np.int8)


def describe(s: Structure, left_out: np.ndarray) -> str:
    """The parts *left_out* in words: ``chain I, E/ATP 355, 103 waters`` (``""`` for none)."""
    if not left_out.any():
        return ""
    residues = np.unique(s.residue_of[left_out])
    polymer = np.isin(s.res_kind, (AMINO, NUCLEIC))
    said: list[str] = []
    loose = 0
    for c in np.unique(s.chain_of[residues[polymer[residues]]]).tolist():
        span = np.arange(*s.residues_of(c).indices(s.n_residues))
        mine = span[polymer[span]]
        if np.isin(mine, residues).all():
            said.append(f"chain {s.chain_label(c)}")
        else:
            loose += int(np.isin(mine, residues).sum())
    if loose:
        said.append(f"{loose} polymer residue{'s' if loose != 1 else ''}")
    said += [s.residue_label(int(r)) for r in residues if s.res_kind[r] == LIGAND]
    for kind, name in KIND_NAMES.items():
        if kind in (AMINO, NUCLEIC, LIGAND):
            continue
        count = int((s.res_kind[residues] == kind).sum())
        if count:
            said.append(f"{count} {name.lower()}{'s' if count != 1 else ''}")
    return ", ".join(said)
