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


def test_attachment_points_draw_differently_from_dummy_atoms() -> None:
    fragment = Chem.MolFromSmiles("Cl[*:1]")
    plain = render_png(fragment, 120, 60)
    attached = render_png(fragment, 120, 60, attachments=True)
    assert attached.startswith(b"\x89PNG") and attached != plain


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


def test_graded_highlight_fades_with_the_ranking_and_ties_share_a_color() -> None:
    from chemsuite.core.depict import graded_colors

    def faded(color: tuple[float, float, float]) -> float:
        """How far the color has blended away from the full-strength highlight."""
        return sum(abs(a - b) for a, b in zip(color, colors[7], strict=True))

    colors = graded_colors(((7, 0.81), (3, 0.05), (9, 0.04)))
    assert list(colors) == [7, 3, 9]
    assert faded(colors[7]) == 0 and faded(colors[3]) < faded(colors[9])
    tied = graded_colors(((7, 0.40), (3, 0.40), (9, 0.02)))
    assert tied[7] == tied[3] != tied[9]  # symmetry-equivalent atoms score the same
    assert len(set(graded_colors(((7, None), (3, None))).values())) == 1  # no scores, no grading
    on_dark = graded_colors(((7, 0.81), (3, 0.05)), dark=True)
    assert on_dark[7] == colors[7] and on_dark[3] != colors[3]  # fades toward the dark base
    assert on_dark[3][2] < on_dark[7][2]  # toward the dark base, not toward white
    assert graded_colors(()) == {}


def test_highlighted_atoms_follow_the_hydrogens_a_drawing_leaves_out() -> None:
    from chemsuite.core.depict import drawn_highlight

    flat = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    atoms, colors = drawn_highlight(depiction_mol(flat), ((3, 0.81), (9, 0.05), (99, 0.01)))
    assert atoms == [3, 9] and set(colors) == {3, 9}  # the index past the end is left out

    three_d = Chem.AddHs(flat)
    AllChem.EmbedMolecule(three_d, randomSeed=1)
    shown = depiction_mol(three_d)  # the hydrogens are not drawn, so the rest renumber
    assert shown.GetNumAtoms() == flat.GetNumAtoms()
    hydrogen = flat.GetNumAtoms()
    atoms, colors = drawn_highlight(shown, ((3, 0.81), (9, 0.05), (hydrogen, 0.01)))
    assert [shown.GetAtomWithIdx(a).GetSymbol() for a in atoms] == ["N", "C"]
    assert set(colors) == set(atoms)  # a hydrogen is not a site and is left out

    # AddHs puts the hydrogens last, so nothing above renumbers. A 3D structure read from a
    # file can carry them between the heavy atoms, and then every later index moves.
    heavy = list(range(flat.GetNumAtoms()))
    hydrogens = list(range(flat.GetNumAtoms(), three_d.GetNumAtoms()))
    order = [i for pair in zip(hydrogens, heavy, strict=False) for i in pair]
    order += heavy[len(hydrogens) :] + hydrogens[len(heavy) :]
    mixed = Chem.RenumberAtoms(three_d, order)
    nitrogen = next(a.GetIdx() for a in mixed.GetAtoms() if a.GetSymbol() == "N")
    assert nitrogen != 3  # the site's index in the molecule is not its index in the drawing
    shown = depiction_mol(mixed)
    atoms, colors = drawn_highlight(shown, ((nitrogen, 0.81),))
    assert [shown.GetAtomWithIdx(a).GetSymbol() for a in atoms] == ["N"]
    assert set(colors) == set(atoms)
