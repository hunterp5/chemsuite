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
"""The selection language."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.protein.structure.model import Structure
from chemsuite.protein.structure.select import SelectionError, residues_in, select


def count(s: Structure, text: str, **kw) -> int:
    return int(select(s, text, **kw).sum())


def test_basic_words(kinase: Structure) -> None:
    assert count(kinase, "all") == kinase.n_atoms
    assert count(kinase, "none") == count(kinase, "") == 0
    assert count(kinase, "chain I") == int((kinase.atom_chain == 1).sum())
    assert (
        count(kinase, "chain E+I")
        == count(kinase, "chain E,I")
        == count(kinase, "polymer or hetero")
    )
    assert count(kinase, "ligand") == 31
    assert count(kinase, "ion") == 2
    assert count(kinase, "resn MN") == 2
    assert count(kinase, "elem Mn") == count(kinase, "elem MN") == 2
    assert count(kinase, "water") == 103
    ca = count(kinase, "name CA")
    assert ca == count(kinase, "protein and name CA") == int((kinase.res_kind == 0).sum())


def test_case_and_wildcards(kinase: Structure) -> None:
    assert count(kinase, "NAME ca") == count(kinase, "name CA")
    c_atoms = select(kinase, "name C*")
    assert set(np.unique(kinase.name[c_atoms])) >= {"C", "CA", "CB", "CG"}
    assert not set(np.unique(kinase.name[c_atoms])) & {"N", "O"}
    assert count(kinase, "resn AT?") == 31  # ATP


def test_residue_numbers_ranges_and_insertion_codes(kinase, thrombin) -> None:
    residues = residues_in(kinase, select(kinase, "chain E and resi 10-50"))
    assert kinase.res_seq[residues].min() == 15 and kinase.res_seq[residues].max() == 50
    assert len(residues_in(kinase, select(kinase, "chain E and resi 10+20+30"))) == 2
    one = residues_in(thrombin, select(thrombin, "chain L and resi 1H"))
    assert [thrombin.residue_label(r) for r in one] == ["L/THR 1H"]
    plain = residues_in(thrombin, select(thrombin, "chain L and resi 1"))
    assert all(thrombin.icode[r] == "" for r in plain)


def test_logic_and_precedence(kinase: Structure) -> None:
    assert count(kinase, "not water") == kinase.n_atoms - 103
    assert count(kinase, "ion or ligand and name P*") == 2 + 3  # and binds tighter
    assert count(kinase, "(ion or ligand) and name P*") == 3
    assert count(kinase, "not not ion") == 2


def test_backbone_sidechain_and_hydrogens(crambin, nmr) -> None:
    backbone, side = count(crambin, "backbone"), count(crambin, "sidechain")
    assert backbone + side == count(crambin, "protein")
    assert backbone == 4 * crambin.n_residues + 1  # N CA C O, and the last OXT
    assert count(nmr, "hydrogen") == int((nmr.element == 1).sum())
    assert count(nmr, "helix") > 0 and count(nmr, "sheet") > 0


def test_byres_bychain_within_and_around(kinase: Structure) -> None:
    site = select(kinase, "byres (protein within 4 of ligand)")
    residues = residues_in(kinase, site)
    labels = {kinase.residue_label(r) for r in residues}
    assert "E/LYS 72" in labels  # the catalytic lysine binds the phosphates
    # every selected residue is complete
    for r in residues:
        span = kinase.atoms_of(r)
        assert site[span].all()
    assert count(kinase, "around 4 of ligand") == count(kinase, "within 4 of ligand") - 31
    assert count(kinase, "bychain resn ATP") == int((kinase.atom_chain == 0).sum())
    assert count(kinase, "within 3 of none") == 0
    assert count(kinase, "water within 3.5 of ion") == count(kinase, "water and within 3.5 of ion")


def test_saved_selections(kinase: Structure) -> None:
    pocket = select(kinase, "byres within 5 of ligand")
    assert count(kinase, "$pocket and not ligand", named={"pocket": pocket}) == int(
        pocket.sum() - 31
    )
    with pytest.raises(SelectionError, match="no saved selection"):
        select(kinase, "$missing")


@pytest.mark.parametrize(
    ("text", "position"),
    [
        ("chain", 5),
        ("protein and", 11),
        ("(protein", 8),
        ("protein )", 8),
        ("frobnicate", 0),
        ("within x of ligand", 7),
        ("resi 1-x", 5),
    ],
)
def test_mistakes_say_where(kinase: Structure, text: str, position: int) -> None:
    with pytest.raises(SelectionError) as error:
        select(kinase, text)
    assert error.value.position == position
