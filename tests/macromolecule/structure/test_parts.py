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

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.macromolecule.structure import edit, parts
from chemsuite.macromolecule.structure.model import LIGAND, WATER, Structure
from tests.macromolecule.conftest import only_ligand


def _bonds_by_uid(s: Structure) -> set[tuple[int, int, int]]:
    return {(*sorted((int(s.uid[a]), int(s.uid[b]))), int(o))
            for (a, b), o in zip(s.bonds.tolist(), s.order.tolist(), strict=True)}  # fmt: skip


def _same(a: Structure, b: Structure) -> None:
    """*a* and *b* hold the same atoms, residues, chains and bonds, in the same order."""
    assert a.uid.tolist() == b.uid.tolist()
    assert np.array_equal(a.coords, b.coords)
    assert a.name.tolist() == b.name.tolist() and a.charge.tolist() == b.charge.tolist()
    assert a.res_start.tolist() == b.res_start.tolist()
    assert [a.residue_label(r) for r in range(a.n_residues)] == [
        b.residue_label(r) for r in range(b.n_residues)
    ]
    assert a.res_kind.tolist() == b.res_kind.tolist()
    assert a.entity.tolist() == b.entity.tolist() and a.seqres.tolist() == b.seqres.tolist()
    assert a.chain_start.tolist() == b.chain_start.tolist()
    assert _bonds_by_uid(a) == _bonds_by_uid(b)


def test_excluded_atoms_are_kept_as_uid_runs(lysozyme) -> None:
    s = lysozyme
    mask = s.atom_kind == WATER
    meta = parts.with_excluded({"title": "x"}, s, mask)
    assert meta["title"] == "x" and "-" in meta[parts.EXCLUDED]
    assert np.array_equal(parts.excluded(meta, s), mask)
    assert parts.EXCLUDED not in parts.with_excluded(meta, s, np.zeros(s.n_atoms, bool))
    assert parts.write_excluded(np.array([5, 1, 2, 3, 9])) == "1-3 5 9"
    assert parts.read_uids("1-3 x 5 9-8").tolist() == [1, 2, 3, 5]
    # an atom deleted since is passed over
    gone = np.zeros(s.n_atoms, bool)
    gone[np.flatnonzero(mask)[0]] = True
    edited = edit.delete(s, gone)
    assert parts.excluded(meta, edited).sum() == mask.sum() - 1


@pytest.mark.parametrize("kind", ["waters", "ligand", "chain"])
def test_restore_puts_back_exactly_what_was_left_out(lysozyme, kind) -> None:
    s = lysozyme
    if kind == "waters":
        left_out = s.atom_kind == WATER
    elif kind == "ligand":
        left_out = s.residue_of == only_ligand(s)
    else:
        left_out = s.atom_chain == 0
    sub = parts.included(s, left_out)
    assert sub.n_atoms == s.n_atoms - left_out.sum()
    back = parts.restore(sub, s, left_out)
    assert sorted(back.uid.tolist()) == sorted(s.uid.tolist())
    assert _bonds_by_uid(back) == _bonds_by_uid(s)
    found = back.indices(s.uid)
    assert np.array_equal(back.coords[:, found], s.coords)
    labels = {back.residue_label(r) for r in range(back.n_residues)}
    assert labels == {s.residue_label(r) for r in range(s.n_residues)}
    if kind != "chain":  # each residue back in its place
        _same(back, s)


def test_restore_keeps_what_the_calculation_changed(lysozyme) -> None:
    s = lysozyme
    waters = s.atom_kind == WATER
    sub = parts.included(s, waters)
    # a calculation adds an atom and moves another
    lig = only_ligand(sub)
    atom = sub.atoms_of(lig).start
    xyz = sub.xyz()[atom][None] + 1.0
    done = edit.add_atoms(sub, np.array([lig]), np.array(["HX"]), np.array([1]), xyz,
                          bonds=np.array([[atom, sub.n_atoms]]))  # fmt: skip
    done = done.with_coords(done.coords + 0.5)
    back = parts.restore(done, s, waters)
    assert back.n_atoms == s.n_atoms + 1 and back.next_uid == done.next_uid
    found = back.indices(s.uid[waters])
    assert np.array_equal(back.coords[:, found], s.coords[:, waters])  # untouched
    assert (back.res_kind[back.residue_of[found]] == WATER).all()
    with pytest.raises(edit.EditError):
        parts.restore(back, s, waters)


def test_a_ligand_copy_joins_its_chain_numbered_on(lysozyme) -> None:
    s = lysozyme
    lig = only_ligand(s)
    copy = parts.copy_residues(s, np.array([lig]))
    assert copy.n_atoms == s.n_atoms + len(range(*s.atoms_of(lig).indices(s.n_atoms)))
    assert copy.n_residues == s.n_residues + 1 and copy.n_chains == s.n_chains
    new = np.flatnonzero(~np.isin(copy.uid, s.uid))
    assert copy.uid[new].min() >= s.next_uid and copy.next_uid == s.next_uid + len(new)
    (r,) = np.unique(copy.residue_of[new])
    assert copy.res_kind[r] == LIGAND and copy.res_name[r] == s.res_name[lig]
    chain = copy.chain_of[r]
    assert copy.chain_name[chain] == s.chain_name[s.chain_of[lig]]
    assert copy.res_seq[r] == s.res_seq[s.residues_of(int(s.chain_of[lig]))].max() + 1
    assert np.array_equal(copy.xyz()[new], s.xyz()[s.atoms_of(lig)])
    inside = sum(1 for a, b in copy.bonds.tolist() if a in new and b in new)
    assert inside == sum(
        1 for a, b in s.bonds.tolist() if s.residue_of[a] == s.residue_of[b] == lig
    )
    # the original is unchanged within the copy
    old = copy.indices(s.uid)
    assert np.array_equal(copy.xyz()[old], s.xyz())


def test_a_chain_copy_is_a_new_chain(crambin) -> None:
    s = crambin
    span = s.residues_of(0)
    name = parts.free_chain_name(s)
    assert name not in s.chain_name.tolist()
    copy = parts.copy_residues(s, np.arange(span.start, span.stop), chain_name=name)
    assert copy.n_chains == s.n_chains + 1 and copy.chain_name[-1] == name
    last = copy.residues_of(copy.n_chains - 1)
    assert copy.res_seq[last].tolist() == s.res_seq[span].tolist()
    assert copy.entity[last].tolist() == s.entity[span].tolist()
    assert len(copy.bonds) == len(s.bonds) * 2 - sum(
        1 for a, b in s.bonds.tolist() if s.atom_chain[a] != 0 or s.atom_chain[b] != 0
    )
    with pytest.raises(edit.EditError):
        parts.copy_residues(s, np.array([0]), chain_name=str(s.chain_name[0]))


def test_a_chain_left_out_whole_comes_back_as_its_own_chain(kinase) -> None:
    s = kinase
    (i,) = [c for c in range(s.n_chains) if s.chain_name[c] == "I"]
    assert i == s.n_chains - 1  # the last chain, so it comes back where it was
    left_out = s.atom_chain == i
    sub = parts.included(s, left_out)
    assert "I" not in sub.chain_name.tolist()
    _same(parts.restore(sub, s, left_out), s)
