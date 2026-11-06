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

"""Clean Structure (SKETCHER_PLAN §4.6): layout quality, place, turn and stereo kept.

The full benchmark (500 approved drugs, both layouts, stereo kept on all 500) is
``scripts/review/layout_review.py``; docs/sketcher/layout-review.md holds its results.
"""

from __future__ import annotations

import math
import random

import pytest
from rdkit import Chem

from chemsuite.sketch import convert, layout
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import Stereo

from .helpers import drawn

#: Drug-like structures with stereocentres, E/Z bonds, fused and bridged rings and a salt.
CASES = [
    "C[C@H](N)C(=O)O",  # L-alanine
    "C/C=C/C(=O)O",  # crotonic acid (E)
    r"C/C=C\[C@@H](F)Cl",
    "O=C1CCC[C@H]2CCCC[C@@H]12",  # trans-fused decalone
    "CN1[C@H]2CC[C@@H]1[C@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2",  # cocaine
    "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O",  # estradiol
    "CC(C)Cc1ccc(C(C)C(=O)[O-])cc1.[Na+]",  # ibuprofen sodium
    "c1ccc2c(c1)[nH]c1ccccc12",  # carbazole
]


def _clean(sketch, atoms=None):
    history = History(sketch)
    change = history.change("Clean")
    problem = layout.clean(change, atoms)
    history.commit(change)
    return history.sketch, problem


def _jittered(smiles: str, seed: int = 3):
    sketch = convert.from_text(smiles)
    rng = random.Random(seed)
    for aid, atom in list(sketch.atoms.items()):
        sketch.atoms[aid] = atom.moved(atom.x + rng.uniform(-0.35, 0.35),
                                       atom.y + rng.uniform(-0.35, 0.35))  # fmt: skip
    return sketch


@pytest.mark.parametrize("smiles", CASES)
def test_clean_keeps_the_molecule_and_its_stereo(smiles):
    sketch = _jittered(smiles)
    before = layout.stereo_key(sketch)
    assert before == Chem.MolToSmiles(convert.to_rdkit(convert.from_text(smiles)))
    cleaned, problem = _clean(sketch)
    assert problem == ""
    assert layout.stereo_key(cleaned) == before
    points = [a.pos for a in cleaned.atoms.values()]
    bonds = [(b.begin, b.end) for b in cleaned.bonds.values()]
    index = {aid: i for i, aid in enumerate(cleaned.atoms)}
    faults = layout.score(points, [(index[a], index[b]) for a, b in bonds])
    assert faults["overlaps"] == 0 and faults["crossings"] == 0
    assert faults["stretched"] <= (4 if "[C@H]2CC[C@@H]1" in smiles else 0)  # a drawn bridge


def test_clean_untangles_a_drawing_with_no_sensible_coordinates():
    sketch, _ = drawn("CC(C)Cc1ccc(cc1)C(C)C(=O)O")  # ibuprofen with atoms on a circle
    cleaned, problem = _clean(sketch)
    assert problem == ""
    assert convert.to_smiles(cleaned) == Chem.CanonSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    assert convert.bond_length_spread(cleaned) < 0.02


def test_clean_keeps_the_drawing_where_it_was_and_does_not_mirror_it():
    sketch = convert.from_text("CC(=O)Nc1ccc(O)cc1")  # paracetamol
    turned = sketch.copy()
    x0, y0, x1, y1 = sketch.bounds()
    centre = ((x0 + x1) / 2, (y0 + y1) / 2)
    c, s = math.cos(1.0), math.sin(1.0)
    for aid, a in sketch.atoms.items():
        dx, dy = a.x - centre[0], a.y - centre[1]
        turned.atoms[aid] = a.moved(10 + dx * c - dy * s, 5 + dx * s + dy * c)
    cleaned, _ = _clean(turned)
    for aid in turned.atoms:  # already clean: every atom stays within a small fraction
        assert math.dist(cleaned.atoms[aid].pos, turned.atoms[aid].pos) < 0.35


def test_clean_selection_moves_only_that_fragment():
    sketch = _jittered("CCO.c1ccccc1", seed=5)
    ethanol = next(a for a in sketch.atoms if sketch.atoms[a].symbol == "O")
    benzene = sketch.fragment(next(iter(set(sketch.atoms) - sketch.fragment(ethanol))))
    cleaned, problem = _clean(sketch, atoms={ethanol})
    assert problem == ""
    for aid in benzene:
        assert cleaned.atoms[aid] == sketch.atoms[aid]


def test_clean_rewedges_when_the_drawn_wedge_no_longer_reads_the_same():
    # Alanine with its atoms scrambled onto a circle: Clean keeps the configuration the
    # scrambled drawing has, whether through the drawn wedge or one RDKit chooses.
    sketch = convert.from_text("C[C@H](N)C(=O)O")
    n = len(sketch.atoms)
    for k, aid in enumerate(sorted(sketch.atoms)):
        t = 2 * math.pi * k / n
        sketch.atoms[aid] = sketch.atoms[aid].moved(2 * math.cos(t), 2 * math.sin(t))
    before = layout.stereo_key(sketch)
    assert "@" in before
    cleaned, problem = _clean(sketch)
    assert problem == ""
    assert layout.stereo_key(cleaned) == before
    assert any(b.stereo is not Stereo.NONE for b in cleaned.bonds.values())


def test_clean_refuses_a_drawing_rdkit_cannot_read():
    sketch, _ = drawn("C(C)(C)(C)(C)C")  # pentavalent carbon
    cleaned, problem = _clean(sketch)
    assert "RDKit cannot read" in problem


def test_fit_turns_but_never_mirrors():
    old = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]
    mirrored = [(0.0, 0.0), (-1.0, 0.0), (-1.0, 1.0)]
    out = layout.fit(mirrored, old)
    # the best rotation of a mirror image cannot match all three points
    assert max(math.dist(a, b) for a, b in zip(out, old, strict=True)) > 0.3
    turned = [(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0)]  # old turned by 90 degrees
    out = layout.fit(turned, old)
    assert max(math.dist(a, b) for a, b in zip(out, old, strict=True)) < 1e-9


@pytest.mark.parametrize("smiles", ["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", CASES[5]])
def test_rewedge_restores_each_centre_from_rdkits_stereo(smiles):
    from dataclasses import replace

    sketch = convert.from_text(smiles)
    before = layout.stereo_key(sketch)
    mol = convert.to_rdkit(sketch)
    bare = sketch.copy()
    bare.bonds = {bid: replace(b, stereo=Stereo.NONE) for bid, b in sketch.bonds.items()}
    assert layout.stereo_key(bare) != before
    layout._rewedge(bare, mol)
    assert layout.stereo_key(bare) == before


def test_rdkit_layout_is_rdkits_own_even_when_coordgen_is_preferred():
    from rdkit.Chem import rdDepictor

    mol = Chem.MolFromSmiles("FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F")
    prefer = rdDepictor.GetPreferCoordGen()
    try:
        rdDepictor.SetPreferCoordGen(False)
        plain = Chem.Mol(mol)
        layout.rdkit_coords(plain)
        rdDepictor.SetPreferCoordGen(True)  # as chemsuite.core.depict sets it
        again = Chem.Mol(mol)
        layout.rdkit_coords(again)
        assert rdDepictor.GetPreferCoordGen()  # restored
    finally:
        rdDepictor.SetPreferCoordGen(prefer)
    a, b = plain.GetConformer().GetPositions(), again.GetConformer().GetPositions()
    assert abs(a - b).max() < 1e-6


def test_text_without_a_drawing_keeps_bridged_stereo():
    # RDKit's wedges on its preferred layout of this quinuclidinium read back contradictory;
    # convert.from_mol then lays it out the other way (found by scripts/review/reaction_review.py)
    smiles = "C[N+]12C[C@H](OC(C)=O)C(CC1)CC2"
    assert convert.to_smiles(convert.from_text(smiles)) == Chem.CanonSmiles(smiles)
