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

from rdkit import Chem

from chemsuite.core.atom_names import AtomName, atoms_named, names_atoms, read_atom_names

#: paracetamol: C0 C1 O2 N3 C4…C9 O10 (atom 3 is the amide nitrogen)
PARACETAMOL = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")


def test_reads_one_atom_and_a_ranked_list_with_scores() -> None:
    assert read_atom_names("N3") == [AtomName("N", 3)]
    assert read_atom_names("N3 0.81, C9 0.05, O10 0.04") == [
        AtomName("N", 3, 0.81),
        AtomName("C", 9, 0.05),
        AtomName("O", 10, 0.04),
    ]
    assert read_atom_names("Cl12 1e-3") == [AtomName("Cl", 12, 0.001)]
    assert read_atom_names("") == []  # a row the tool left blank is not a failure


def test_refuses_text_that_is_not_atom_names() -> None:
    for text in ("not a site", "C8H10N4O2", "7", "C", "C7 strong", "C-3", "3A4"):
        assert read_atom_names(text) is None, text


def test_a_column_names_atoms_when_every_value_does() -> None:
    assert names_atoms(["N3 0.81, C9 0.05", "", "C7 0.40"])
    assert not names_atoms(["", "", ""])  # nothing to highlight
    assert not names_atoms(["N3 0.81", "see footnote"])
    assert not names_atoms(["C8H10N4O2", "C9H8O4"])  # formulas, not atom names


def test_names_become_atom_indices_checked_against_the_molecule() -> None:
    assert atoms_named(PARACETAMOL, "N3 0.81, C9 0.05") == ((3, 0.81), (9, 0.05))
    assert atoms_named(PARACETAMOL, "N3") == ((3, None),)
    assert atoms_named(PARACETAMOL, "C3 0.81") == ()  # atom 3 is nitrogen, not carbon
    assert atoms_named(PARACETAMOL, "C99 0.81") == ()  # past the end of the molecule
    assert atoms_named(PARACETAMOL, "N3 0.81, C3 0.05") == ((3, 0.81),)  # the wrong one dropped
    assert atoms_named(PARACETAMOL, "N3 0.81, N3 0.05") == ((3, 0.81),)  # named twice, kept once
    assert atoms_named(PARACETAMOL, "prose") == ()
