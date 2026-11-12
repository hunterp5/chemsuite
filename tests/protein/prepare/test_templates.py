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

from chemsuite.protein.prepare import templates


@pytest.mark.parametrize("name", templates.BUNDLED)
def test_internal_coordinates_rebuild_every_ideal_hydrogen(name: str) -> None:
    t = templates.bundled()[name]
    x = t.xyz
    for p in t.hydrogens:
        i, par, a, d = (t.index[n] for n in (p.atom, p.parent, p.angle_atom, p.dihedral_atom))
        at = templates.place(x[d], x[a], x[par], p.length, p.angle, p.dihedral)
        assert np.linalg.norm(at - x[i]) < 0.01, (name, p.atom)
    placed = {p.atom for p in t.hydrogens}
    lone = {t.atoms[i] for i in range(len(t.atoms)) if t.is_hydrogen(i)} - placed
    # the caps' NH hydrogens are placed from the residue before them, not from the cap alone
    assert not lone or name in ("NH2", "NME"), lone


def test_rotors_are_the_hydrogens_that_can_turn() -> None:
    def rotors(name: str) -> set[str]:
        return {p.atom for p in templates.bundled()[name].hydrogens if p.rotor}

    assert rotors("SER") == {"HG", "HXT"}
    assert rotors("LYS") == {"HZ1", "HZ2", "HZ3", "HXT"}  # the CCD draws Lys charged
    assert {"HB1", "HB2", "HB3"} <= rotors("ALA")
    assert not {"HD21", "HD22"} & rotors("ASN")  # a planar amide
    assert not {"HH11", "HH12", "HH21", "HH22"} & rotors("ARG")


def test_a_template_gains_a_tetrahedral_hydrogen() -> None:
    ala = templates.bundled()["ALA"]  # the CCD draws the N-terminus as NH2
    charged = ala.with_hydrogen("N", "H3")
    n = charged.xyz[charged.index["N"]]
    around = [charged.xyz[charged.index[a]] for a in ("CA", "H", "H2", "H3")]
    angles = [np.degrees(templates.angle_of(a, n, b))
              for k, a in enumerate(around) for b in around[k + 1:]]  # fmt: skip
    assert all(100 < a < 118 for a in angles), angles
    assert "H3" in charged.hydrogens_on("N")
    assert len(charged.charges) == len(charged.atoms)
    assert charged.without(["H3"]).atoms == ala.atoms


def test_a_ccd_entry_reads_like_the_bundled_one(tmp_path) -> None:
    from chemsuite.protein.structure.components import CcdCache

    his = templates.bundled()["HIS"]
    # a residue not bundled is read from the CCD cache's file, which holds the same entry
    (tmp_path / "XHS.cif").write_text(_as_ccd(his, "XHS"), encoding="utf-8")
    library = templates.Library(CcdCache(tmp_path, offline=True))
    read = library.get("XHS")
    assert read is not None
    assert read.atoms == his.atoms and read.bonds == his.bonds
    assert np.allclose(read.xyz, his.xyz, atol=1e-3)
    assert library.get("ZZZ") is None  # offline and not in the folder


def _as_ccd(t: templates.Template, code: str) -> str:
    import gemmi

    rows = "\n".join(
        f"{code} {a} {gemmi.Element(z).name} {c} {x:.3f} {y:.3f} {w:.3f}"
        for a, z, c, (x, y, w) in zip(t.atoms, t.elements, t.charges, t.xyz, strict=True)
    )
    orders = {1: "SING", 2: "DOUB", 3: "TRIP"}
    bonds = "\n".join(f"{code} {t.atoms[i]} {t.atoms[j]} {orders[o]}" for i, j, o in t.bonds)
    return (
        f"data_{code}\n_chem_comp.id {code}\nloop_\n_chem_comp_atom.comp_id\n"
        "_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n"
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n"
        f"_chem_comp_atom.pdbx_model_Cartn_z_ideal\n{rows}\nloop_\n_chem_comp_bond.comp_id\n"
        "_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n"
        f"{bonds}\n"
    )
