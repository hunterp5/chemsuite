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
from rdkit import Chem

from chemsuite.macromolecule.prepare import sites
from chemsuite.macromolecule.structure import edit
from chemsuite.macromolecule.structure.model import Structure
from chemsuite.simulate import build, spec


def _texts(found: list[build.Problem]) -> list[str]:
    return [str(p) for p in found]


def _atom(s: Structure, label: str) -> int:
    (i,) = [i for i in range(s.n_atoms) if s.atom_label(i) == label]
    return i


def _without_hydrogen_on(s: Structure, label: str) -> Structure:
    """*s* less one hydrogen of the atom *label*."""
    parent = _atom(s, label)
    partners = [b if a == parent else a for a, b in s.bonds.tolist() if parent in (a, b)]
    gone = np.zeros(s.n_atoms, dtype=bool)
    gone[[j for j in partners if s.element[j] == 1][0]] = True
    return edit.delete(s, gone)


def test_a_protein_becomes_its_amber_forms(crambin) -> None:
    s = crambin.structure
    system = build.system_spec(s, crambin.states)
    assert system.n_atoms == s.n_atoms and set(system.uid) == set(s.uid.tolist())
    names = [r.name for r in system.residues]
    assert names[0] == "NTHR" and names[-1] == "CASN"
    assert names.count("CYX") == 6 and len(system.disulfides) == 3
    assert all(r.kind == spec.PROTEIN for r in system.residues)
    for i, j in system.disulfides:
        assert system.name[i] == system.name[j] == "SG"
        assert (i, j, 1) in system.bonds
    assert system.net_charge == int(s.charge.sum())
    # every atom of a residue together, in residue order
    of = system.residue_of()
    assert all(system.residues[of[k]].label == s.residue_label(int(s.residue_of[s.index_of_uid[u]]))
               for k, u in enumerate(system.uid))  # fmt: skip
    xyz = build.positions(s, system)
    assert np.allclose(xyz, s.xyz()[s.indices(np.array(system.uid))])


def _without(prepared, names: tuple[str, ...]) -> tuple[Structure, dict[int, str]]:
    """The prepared structure less the residues named *names*, its states renumbered."""
    s = prepared.structure
    keep = ~np.isin(s.res_name[s.residue_of], names)
    out = s.subset(keep)
    where = {s.residue_label(r): r for r in range(s.n_residues)}
    labels = {out.residue_label(r): r for r in range(out.n_residues)}
    states = {labels[label]: prepared.states[where[label]] for label in labels
              if where[label] in prepared.states}  # fmt: skip
    return out, states


def test_bond_orders_guessed_from_the_geometry_block_the_build(lysozyme) -> None:
    found = _texts(build.problems(lysozyme.structure, lysozyme.states))
    assert {t.split(":")[0] for t in found} == {"A/PO4 165", "A/PO4 166", "A/BME 168"}
    assert all("guess" in t for t in found)


def test_a_ligand_is_a_molecule_with_the_hydrogens_prepare_gave_it(lysozyme) -> None:
    s, states = _without(lysozyme, ("PO4", "BME"))
    system = build.system_spec(s, states)
    (molecule,) = system.molecules
    residue = system.residues[molecule.residue]
    assert residue.kind == spec.MOLECULE
    mol = Chem.MolFromMolBlock(molecule.molblock, removeHs=False)
    assert mol.GetNumAtoms() == residue.count
    symbols = [a.GetAtomicNum() for a in mol.GetAtoms()]
    assert symbols == list(system.element[residue.first : residue.first + residue.count])
    assert any(e == 1 for e in symbols)
    assert not system.charged  # until the worker computes them
    # the residue capped where the deposit stops short of the sequence's end
    assert [r.name for r in system.residues if r.kind == spec.CAP] == ["NME"]
    assert system.residues[[r.name for r in system.residues].index("NME") - 1].name == "ASN"


def test_the_states_prepare_recorded_must_match_the_hydrogens(lysozyme) -> None:
    s, states = _without(lysozyme, ("PO4", "BME"))
    his = next(r for r in range(s.n_residues) if s.res_name[r] == "HIS")
    system = build.system_spec(s, states)
    form = system.residues[[r.label for r in system.residues].index(s.residue_label(his))].name
    assert form in ("HID", "HIE", "HIP") and states.get(his) == form
    wrong = {**states, his: "HIP" if form != "HIP" else "HID"}
    with pytest.raises(build.StateMismatch, match="HIS 31"):
        build.system_spec(s, wrong)


def test_what_no_template_covers_blocks_the_build(kinase, egfr, crambin) -> None:
    found = _texts(build.problems(kinase.structure, kinase.states))
    assert any("TPO 197" in t and "20 amino acids" in t for t in found)
    assert any("SEP 338" in t for t in found)
    covalent = _texts(build.problems(egfr.structure, egfr.states))
    assert any("80U" in t and "CYS 797" in t for t in covalent)
    with pytest.raises(build.BuildError) as raised:
        build.system_spec(kinase.structure, kinase.states)
    assert len(raised.value.problems) == len(found)
    # tyrosinate, a neutral N-terminus, an ion with no charge, a water without hydrogens
    s = crambin.structure
    tym = _without_hydrogen_on(s, "A/TYR 29/OH")
    assert any("Tyrosinate" in t for t in _texts(build.problems(tym)))
    neutral = _without_hydrogen_on(s, "A/THR 1/N")
    assert any("neutral N-terminus" in t for t in _texts(build.problems(neutral)))
    mn = next(i for i in range(kinase.structure.n_atoms) if kinase.structure.res_name[
        kinase.structure.residue_of[i]] == "MN")  # fmt: skip
    mask = np.zeros(kinase.structure.n_atoms, dtype=bool)
    mask[mn] = True
    uncharged = edit.set_charge(kinase.structure, mask, 0)
    assert any("unknown charge" in t for t in _texts(build.problems(uncharged)))
    w = egfr.structure
    water = next(r for r in range(w.n_residues) if w.res_name[r] == "HOH")
    one_h = [i for i in range(*w.atoms_of(water).indices(w.n_atoms)) if w.element[i] == 1][0]
    mask = np.zeros(w.n_atoms, dtype=bool)
    mask[one_h] = True
    assert any("two hydrogens" in t for t in _texts(build.problems(edit.delete(w, mask))))


def test_a_mid_chain_form_at_a_chain_end_blocks(crambin) -> None:
    s = crambin.structure
    last = s.n_residues - 1
    assert s.res_name[last] == "ASN"
    # crambin's C-terminal Asn with a proton on OXT is a neutral C-terminus: no template
    oxt = _atom(s, f"{s.residue_label(last)}/OXT")
    mask = np.zeros(s.n_atoms, dtype=bool)
    mask[oxt] = True
    charge = edit.set_charge(s, mask, 0)
    proton = edit.add_atoms(charge, np.array([last]), np.array(["HXT"]), np.array([1]),
                            s.xyz()[oxt][None] + [0.0, 0.0, 0.97],
                            bonds=np.array([[oxt, s.n_atoms]]))  # fmt: skip
    assert any("neutral C-terminus" in t for t in _texts(build.problems(proton)))


def test_the_neutral_n_terminus_state_is_named_alike_on_both_sides() -> None:
    assert build.NEUTRAL_N == sites.NEUTRAL_N
    from chemsuite.macromolecule.prepare import solve

    assert build.PROTEIN_STATES == solve.RECORDED
