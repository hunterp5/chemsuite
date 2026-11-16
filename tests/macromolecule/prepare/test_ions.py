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

from chemsuite.macromolecule.prepare import energy, hydrogens, ions
from chemsuite.macromolecule.structure import io
from chemsuite.macromolecule.structure.model import ION, Structure


def test_the_name_fixes_the_oxidation_state_before_the_element() -> None:
    assert ions.charge("FE", 26) == 3 and ions.charge("FE2", 26) == 2
    assert ions.charge("CU1", 29) == 1 and ions.charge("CU", 29) == 2
    assert ions.charge("NA", 11) == 1 and ions.charge("ZN", 30) == 2
    assert ions.charge("CL", 17) == -1 and ions.charge("MG", 12) == 2
    assert ions.charge("FE", 26, stated=2) == 2  # the file's own charge wins
    assert ions.charge("MNX", 25) is None  # Mn has no one common ion


def _ions(*names: tuple[str, str]) -> Structure:
    """One ion per (residue name, element), 10 Å apart, with blank charge columns."""
    lines = [
        f"HETATM{k + 1:5d} {res:<4} {res:>3} A{k + 1:4d}    {10.0 * k:8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"  1.00 20.00          {element:>2}  "
        for k, (res, element) in enumerate(names)
    ]
    return io.parse("\n".join([*lines, "END", ""]), "pdb")


def test_ions_with_no_charge_in_the_file_get_one() -> None:
    s = _ions(("ZN", "ZN"), ("FE2", "FE"), ("NA", "NA"), ("CL", "CL"))
    assert (s.atom_kind == ION).all() and not s.charge.any()
    out = hydrogens.add_hydrogens(s)
    assert out.structure.charge.tolist() == [2, 2, 1, -1]
    assert not out.notes


def test_an_ion_of_unknown_charge_is_reported_not_guessed() -> None:
    s = _ions(("ZN", "ZN"), ("MNX", "MN"))
    out = hydrogens.add_hydrogens(s)
    assert out.structure.charge.tolist() == [2, 0]
    assert [(n.residue, n.text) for n in out.notes] == [(1, hydrogens.UNKNOWN_ION)]


def test_the_energy_model_charges_ions_from_the_same_table() -> None:
    s = _ions(("FE", "FE"), ("NA", "NA"), ("MNX", "MN"))
    q = energy.group_charges(s)
    assert q.tolist() == [3.0, 1.0, float(energy.UNKNOWN_METAL)]
