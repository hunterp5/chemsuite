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
"""Interactions in known pockets, against what the papers and PLIP report."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.protein.chem import interactions as it
from chemsuite.protein.structure import io
from chemsuite.protein.structure.model import Structure


def _ligand(s: Structure, name: str) -> int:
    return int(np.flatnonzero(s.res_name == name)[0])


def _found(s: Structure, name: str) -> it.Interactions:
    return it.find(s, it.ligand_query(s, _ligand(s, name)))


def _pairs(s: Structure, found: it.Interactions, kind: int) -> dict[tuple[str, str], float]:
    """{(query atom name, "RES123 NAME"): distance} for one kind."""
    out = {}
    for k in np.flatnonzero(found.kind == kind).tolist():
        a, b = int(found.a[k]), int(found.b[k])
        r = int(s.residue_of[b])
        partner = f"{s.res_name[r]}{s.res_seq[r]}{s.icode[r]} {s.name[b]}"
        out[(str(s.name[a]), partner)] = float(found.distance[k])
    return out


def test_ppack_in_thrombin(thrombin) -> None:
    found = _found(thrombin, "0G6")
    salt = _pairs(thrombin, found, it.SALT_BRIDGE)
    assert salt == {("NH2", "ASP189 OD2"): pytest.approx(2.68, abs=0.01)}  # the S1 pocket
    hbonds = _pairs(thrombin, found, it.HBOND)
    for pair in [("N", "GLY216 O"), ("O", "GLY216 N"),  # the antiparallel strand
                 ("O2", "GLY193 N"), ("O2", "SER195 N")]:  # the oxyanion hole  # fmt: skip
        assert pair in hbonds
    assert ("NH2", "ASP189 OD2") not in hbonds  # the salt bridge, drawn once
    assert all(it.HBOND_MIN <= d <= it.HBOND_MAX for d in hbonds.values())


def test_atp_in_protein_kinase_a(kinase) -> None:
    found = _found(kinase, "ATP")
    hbonds = _pairs(kinase, found, it.HBOND)
    assert hbonds[("N6", "GLU121 O")] == pytest.approx(3.01, abs=0.01)  # the hinge
    assert ("N1", "VAL123 N") in hbonds
    assert {b for _a, b in _pairs(kinase, found, it.SALT_BRIDGE)} == {"LYS72 NZ", "LYS168 NZ"}
    metals = _pairs(kinase, found, it.METAL)
    assert len(metals) == 4 and {b for _a, b in metals} == {"MN351 MN", "MN352 MN"}
    assert all(d < 2.3 for d in metals.values())
    # no hydrogen bond to a water: those are water bridges, each ending at the protein
    assert not any(b.startswith("HOH") for _a, b in hbonds)
    bridges = found.kind == it.WATER_BRIDGE
    waters = kinase.res_name[kinase.residue_of[found.b[bridges]]] == "HOH"
    ends = kinase.res_name[kinase.residue_of[found.b[bridges][~waters]]]
    assert waters.any() and "HOH" not in ends.tolist()


def test_a_covalent_inhibitor_in_egfr(egfr) -> None:
    found = _found(egfr, "80U")
    hbonds = _pairs(egfr, found, it.HBOND)
    assert hbonds[("N1", "MET793 N")] == pytest.approx(3.22, abs=0.01)  # the hinge
    assert hbonds[("NAB", "GLN791 O")] == pytest.approx(3.48, abs=0.01)
    halogen = _pairs(egfr, found, it.HALOGEN)
    assert halogen == {("CLA", "LEU788 O"): pytest.approx(3.21, abs=0.01)}
    # the warhead's own cysteine is bonded, not an interaction
    assert not any(b.startswith("CYS797") for _a, b in hbonds)


def test_a_hydrophobic_pocket(lysozyme) -> None:
    found = _found(lysozyme, "JZ4")
    contacts = _pairs(lysozyme, found, it.HYDROPHOBIC)
    residues = {b.split()[0] for _a, b in contacts}
    assert {"ALA99", "LEU118", "PHE153", "VAL103", "LEU84"} <= residues
    assert len(residues) == len(contacts)  # one line per residue
    assert all(d <= it.HYDROPHOBIC_MAX for d in contacts.values())
    assert found.counts().get(it.SALT_BRIDGE, 0) == 0


def test_lines_start_at_the_query(kinase) -> None:
    query = it.ligand_query(kinase, _ligand(kinase, "ATP"))
    found = it.find(kinase, query)
    assert np.allclose(found.start, kinase.xyz(0)[found.a], atol=1e-4)
    direct = found.kind != it.WATER_BRIDGE
    assert query[found.a[direct]].all() and not query[found.b[direct]].any()
    assert np.allclose(np.linalg.norm(found.end - found.start, axis=1), found.distance)


def test_a_peptide_can_be_the_query(kinase) -> None:
    peptide = kinase.atom_chain == 1
    found = it.find(kinase, peptide)
    assert found.counts()[it.SALT_BRIDGE] >= 3  # PKI's arginines on the kinase's acids
    assert not peptide[found.b[found.kind != it.WATER_BRIDGE]].any()


def test_nothing_to_find() -> None:
    s = io.parse(_PDB.replace("CL ", "C  ").replace("CL1", "C7 "), "pdb")
    assert len(it.find(s, np.zeros(s.n_atoms, bool))) == 0


_PDB = """\
HETATM    1  C1  CLB A 900       0.000   1.400   0.000  1.00 10.00           C
HETATM    2  C2  CLB A 900       1.212   0.700   0.000  1.00 10.00           C
HETATM    3  C3  CLB A 900       1.212  -0.700   0.000  1.00 10.00           C
HETATM    4  C4  CLB A 900       0.000  -1.400   0.000  1.00 10.00           C
HETATM    5  C5  CLB A 900      -1.212  -0.700   0.000  1.00 10.00           C
HETATM    6  C6  CLB A 900      -1.212   0.700   0.000  1.00 10.00           C
HETATM    7 CL1  CLB A 900       0.000   3.150   0.000  1.00 10.00          CL
ATOM      8  N   GLY B   1      -1.000   7.300   0.000  1.00 10.00           N
ATOM      9  CA  GLY B   1      -1.000   8.750   0.000  1.00 10.00           C
ATOM     10  C   GLY B   1       0.000   7.400   0.000  1.00 10.00           C
ATOM     11  O   GLY B   1       0.000   6.180   0.000  1.00 10.00           O
END
"""


def test_a_halogen_bond_needs_the_c_x_line_to_point_at_the_acceptor() -> None:
    s = io.parse(_PDB, "pdb")
    found = _found(s, "CLB")
    assert _pairs(s, found, it.HALOGEN) == {("CL1", "GLY1 O"): pytest.approx(3.03, abs=0.01)}
    bent = s.with_coords(s.coords + np.where(s.name == "O", 1, 0)[None, :, None]
                         * np.array([2.5, -0.8, 0.0], np.float32))  # fmt: skip
    assert it.HALOGEN not in _found(bent, "CLB").counts()


def test_describe(kinase) -> None:
    found = _found(kinase, "ATP")
    k = int(np.flatnonzero(found.kind == it.METAL)[0])
    text = it.describe(kinase, found, k)
    assert text.startswith("Metal  E/ATP 355/O") and text.endswith(" Å")


def test_a_display_shows_and_saves_its_interactions(kinase) -> None:
    from chemsuite.protein.look.display import Display, compose

    shown = Display("1ATP", "1ATP", kinase)
    assert len(compose([shown]).lines) == 0
    shown.show_interactions(it.ligand_query(kinase, _ligand(kinase, "ATP")))
    scene = compose([shown])
    found = scene.interactions[0]
    assert len(found) > 20 and shown.interactions() is found  # cached
    assert set(np.unique(scene.lines.owner)) == set(found.a.tolist())
    hbond = np.array(it.COLOURS[it.HBOND], np.float32)
    assert np.isclose(scene.lines_look[:, :3], hbond).all(axis=1).any()
    again = Display("1ATP", "1ATP", kinase)
    again.set_state(shown.state())
    assert np.array_equal(again.query, shown.query)
    shown.show_interactions(None)
    assert len(compose([shown]).lines) == 0
