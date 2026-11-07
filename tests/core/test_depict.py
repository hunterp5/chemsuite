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

import time

from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.core.depict import depiction_mol, render_png

from .conftest import DRUGS


def test_keeps_file_2d_coordinates_and_replaces_3d() -> None:
    flat = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(flat)
    assert depiction_mol(flat) is flat
    three_d = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(three_d, randomSeed=1)
    drawn = depiction_mol(three_d)
    assert drawn is not three_d and not drawn.GetConformer().Is3D()
    assert three_d.GetConformer().Is3D()  # the stored molecule is untouched


def test_renders_png_quickly_in_light_and_dark() -> None:
    mol = Chem.MolFromSmiles(DRUGS["imatinib"])
    start = time.perf_counter()
    for dark in (False, True) * 10:
        png = render_png(mol, 220, 110, scale=2.0, dark=dark)
    per_mol = (time.perf_counter() - start) / 20
    assert png.startswith(b"\x89PNG")
    assert per_mol < 0.05, f"{per_mol * 1000:.1f} ms per structure"


def test_unkekulizable_structure_still_draws() -> None:
    mol = Chem.MolFromSmiles("c1cccc1", sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    assert render_png(mol, 100, 80).startswith(b"\x89PNG")


def test_aligned_drawings_put_the_core_where_the_core_is_drawn() -> None:
    from chemsuite.core.depict import aligned_mol, alignment_core

    core = alignment_core("c1ccc2ncccc2c1")
    assert core is not None and core.GetNumConformers() == 1
    want = core.GetConformer()
    for smiles in ("Cc1ccc2ncccc2c1", "OCCc1cccc2ncccc12"):
        mol = Chem.MolFromSmiles(smiles)
        drawn = aligned_mol(mol, core)
        assert drawn is not mol and Chem.MolToSmiles(drawn) == Chem.MolToSmiles(mol)
        got = drawn.GetConformer()
        for core_atom, atom in enumerate(drawn.GetSubstructMatch(core)):
            a, b = got.GetAtomPosition(atom), want.GetAtomPosition(core_atom)
            assert abs(a.x - b.x) < 0.01 and abs(a.y - b.y) < 0.01

    smarts = alignment_core("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1")  # SMARTS works as well
    assert smarts is not None
    other = aligned_mol(Chem.MolFromSmiles("CCO"), smarts)  # no core: drawn as usual
    assert other.GetNumConformers() == 1
    assert alignment_core("not [[ a structure") is None
