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
from scipy.spatial import cKDTree

from chemsuite.protein.prepare import hydrogens, templates
from chemsuite.protein.structure import io
from chemsuite.protein.structure.model import Structure
from tests.protein.conftest import load


def _heavy(s: Structure) -> Structure:
    return s.subset(s.element != 1)


def _on(s: Structure, residue: int, parent: str) -> list[str]:
    """The hydrogens bonded to *parent* in *residue*."""
    span = s.atoms_of(residue)
    (p,) = [i for i in range(span.start, span.stop) if s.name[i] == parent]
    partners = [b if a == p else a for a, b in s.bonds.tolist() if p in (a, b)]
    return sorted(str(s.name[j]) for j in partners if s.element[j] == 1)


def _positions(s: Structure, residue: int) -> dict[str, np.ndarray]:
    span = s.atoms_of(residue)
    return {str(s.name[i]): s.xyz()[i] for i in range(span.start, span.stop)}


def _residue(s: Structure, label: str) -> int:
    (r,) = [r for r in range(s.n_residues) if s.residue_label(r) == label]
    return r


def test_hydrogens_land_where_the_nmr_structure_has_them() -> None:
    nmr = load("1d3z_2models.cif.gz")  # ubiquitin, two models, refined with hydrogens
    out = hydrogens.add_hydrogens(_heavy(nmr))
    s = out.structure
    assert (s.element == 1).sum() == (nmr.element == 1).sum()
    assert s.n_models == 2 and not out.notes
    turning = set(s.indices(np.array([u for r in out.rotors for u in r.hydrogens])).tolist())
    placed = [i for i in np.flatnonzero(s.element == 1) if i not in turning
              and s.res_name[s.residue_of[i]] != "HIS"]  # fmt: skip
    deposited = np.flatnonzero(nmr.element == 1)
    for model in range(2):  # judged by position, since names of CH2 pairs vary by convention
        d, _ = cKDTree(nmr.coords[model][deposited]).query(s.coords[model][placed])
        assert d.mean() < 0.05 and d.max() < 0.25, (model, d.mean(), d.max())


def test_textbook_states_and_termini(crambin) -> None:
    out = hydrogens.add_hydrogens(crambin)
    s = out.structure
    first, last = 0, s.n_residues - 1
    assert _on(s, first, "N") == ["H", "H2", "H3"] and s.charge[s.atoms_of(first)][0] == 1
    assert _on(s, last, "OXT") == [] and _on(s, last, "O") == []
    for r in range(1, s.n_residues):
        name = str(s.res_name[r])
        assert _on(s, r, "N") == ([] if name == "PRO" else ["H"]), s.residue_label(r)
        if name == "CYS":  # every Cys of crambin is in a disulfide
            assert _on(s, r, "SG") == []
        if name == "ASP":
            assert _on(s, r, "OD2") == []
        if name == "ARG":
            assert len(_on(s, r, "NH1") + _on(s, r, "NH2")) == 4
    names = s.res_name
    expected = (sum(names == "ARG") + sum(names == "LYS") - sum(names == "ASP")
                - sum(names == "GLU"))  # fmt: skip
    assert int(s.charge.sum()) == expected  # the termini cancel
    assert out.states == {}
    # the amide H lies in the peptide plane, opposite the carbonyl O
    r = 5
    at, prev = _positions(s, r), _positions(s, r - 1)
    dihedral = np.degrees(templates.dihedral_of(prev["O"], prev["C"], at["N"], at["H"]))
    assert abs(abs(dihedral) - 180) < 10, dihedral


def test_rotors_turn_polar_and_methyl_hydrogens(crambin) -> None:
    out = hydrogens.add_hydrogens(crambin)
    s = out.structure
    by_parent = {s.atom_label(int(s.indices(np.array([r.parent]))[0])): r for r in out.rotors}
    thr = [label for label in by_parent if "/THR " in label and not label.endswith("/N")]
    assert {label.rsplit("/", 1)[1] for label in thr} == {"OG1", "CG2"}
    n_term = by_parent[s.atom_label(int(s.atoms_of(0).start))]
    assert len(n_term.hydrogens) == 3  # NH3+ turns about CA–N
    assert not [label for label in by_parent if label.endswith("/NH1")]  # a planar NH2


def test_a_covalent_ligand_takes_a_hydrogens_place(thrombin) -> None:
    out = hydrogens.add_hydrogens(thrombin)
    s = out.structure
    his, ser = _residue(s, "H/HIS 57"), _residue(s, "H/SER 195")
    assert out.states[his] == "HID"  # NE2 is bonded to the ligand, so ND1 keeps the H
    assert _on(s, his, "ND1") == ["HD1"] and _on(s, his, "NE2") == []
    assert _on(s, ser, "OG") == []
    (ligand,) = np.flatnonzero(s.res_name == "0G6")
    span = s.atoms_of(int(ligand))
    assert (s.element[span] == 1).any()
    # no ligand hydrogen sits on top of the protein atom it is bonded to
    h = np.flatnonzero((s.element == 1) & (s.residue_of == ligand))
    protein = np.flatnonzero((s.element != 1) & (s.residue_of != ligand))
    d, _ = cKDTree(s.xyz()[protein]).query(s.xyz()[h])
    assert d.min() > 1.5


def test_metals_take_the_place_of_a_proton(kinase) -> None:
    s = _with_zinc(kinase, ["E/CYS 343/SG", "E/HIS 87/NE2"])
    out = hydrogens.add_hydrogens(s)
    p = out.structure
    cys, his = _residue(p, "E/CYS 343"), _residue(p, "E/HIS 87")
    assert _on(p, cys, "SG") == [] and out.states[cys] == "CYM"
    assert p.charge[p.atoms_of(cys)][p.name[p.atoms_of(cys)] == "SG"][0] == -1
    assert out.states[his] == "HID" and _on(p, his, "ND1") == ["HD1"]
    other = _residue(p, "E/CYS 199")
    assert _on(p, other, "SG") == ["HG"] and other not in out.states


def test_ligands_take_their_components_hydrogens(kinase) -> None:
    from rdkit import Chem

    from chemsuite.protein.chem import ligands

    out = hydrogens.add_hydrogens(kinase)
    s = out.structure
    (atp,) = np.flatnonzero(s.res_name == "ATP")
    mol = ligands.to_mol(kinase, int(atp)).mol
    expected = Chem.AddHs(mol).GetNumAtoms() - mol.GetNumAtoms()
    assert (s.element[s.atoms_of(int(atp))] == 1).sum() == expected
    names = s.name[s.atoms_of(int(atp))]
    assert len(set(names.tolist())) == len(names)
    # waters and ions get none
    for kind in ("HOH", "MN"):
        for r in np.flatnonzero(s.res_name == kind):
            assert not (s.element[s.atoms_of(int(r))] == 1).any()


def test_residues_with_hydrogens_are_left_alone() -> None:
    nmr = load("1d3z_2models.cif.gz")
    out = hydrogens.add_hydrogens(nmr)
    assert out.structure.n_atoms == nmr.n_atoms
    assert len(out.notes) == nmr.n_residues


def test_uids_are_kept_and_the_result_saves(crambin) -> None:
    s = hydrogens.add_hydrogens(crambin).structure
    assert s.uid[s.element != 1].tolist() == crambin.uid.tolist()
    again = io.parse(io.to_text(s, "pdb"), "pdb")
    assert (again.element == 1).sum() == (s.element == 1).sum()


def _with_zinc(s: Structure, atoms: list[str]) -> Structure:
    """*s* with a zinc ion 2.2 Å out from each named atom, along its bond from its parent."""
    text = io.to_text(s, "pdb")
    lines = [line for line in text.splitlines() if not line.startswith(("END", "MASTER"))]
    xyz = s.xyz()
    for k, label in enumerate(atoms):
        (i,) = [i for i in range(s.n_atoms) if s.atom_label(i) == label]
        partner = [b if a == i else a for a, b in s.bonds.tolist() if i in (a, b)]
        direction = xyz[i] - np.mean(xyz[partner], axis=0)
        x, y, z = xyz[i] + 2.2 * direction / np.linalg.norm(direction)
        lines.append(f"HETATM{9000 + k:5d} ZN    ZN Z{900 + k:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                     f"  1.00 20.00          ZN  ")  # fmt: skip
    return io.parse("\n".join([*lines, "END", ""]), "pdb")
