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

"""Uni-pKa's microstates, conformers and thermodynamics (the network has its own tests)."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.ionisation.unipka import Ensemble, conformers, microstates, model, thermo

DATA = Path(__file__).resolve().parents[1] / "data" / "ionisation"


def _canon(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_microstates_match_upstream_on_sampl6() -> None:
    """Upstream's ``enumerator/main.py ensemble -u 2 -l -2 -t simple_smarts_pattern.tsv`` on
    SAMPL6 (its ``example_out.tsv``, Apache-2.0): the same microstates in every macrostate."""
    rows = list(csv.reader((DATA / "unipka_sampl6_ensemble.tsv").open(encoding="utf-8"),
                           delimiter="\t"))  # fmt: skip
    head = rows[0]
    microstates.use_sites(microstates.SIMPLE_SITES)
    try:
        for row in rows[1:]:
            want = {int(q): {_canon(s) for s in cell.split(",")}
                    for q, cell in zip(head[1:-1], row[1:-1], strict=True) if cell}  # fmt: skip
            found = microstates.ensemble(Chem.MolFromSmiles(row[-1]))
            got = {q: {_canon(m.smiles) for m in pool} for q, pool in found.items()}
            assert got == want, row[-1]
    finally:
        microstates.use_sites(None)


def test_full_templates_find_weak_nh_acids() -> None:
    """A hydantoin and a dibenzamide (Novartis acids, pKa 8.3 and 9.9) have no acid site in the
    simple templates; the bundled full set deprotonates their imide N–H."""
    for smiles in ("O=C1CN(c2ccccc2)C(=O)N1", "O=C(NC(=O)c1ccccc1)c1ccccc1"):
        assert -1 in microstates.ensemble(Chem.MolFromSmiles(smiles)), smiles
    microstates.use_sites(microstates.SIMPLE_SITES)
    try:
        assert -1 not in microstates.ensemble(Chem.MolFromSmiles("O=C1CN(c2ccccc2)C(=O)N1"))
    finally:
        microstates.use_sites(None)


def test_microstates_map_onto_the_input_atoms() -> None:
    """4-(pyridin-4-yl)butanoic acid drawn acid first: each microstate's hydrogens and charges
    sit on the input's own atom indices."""
    mol = Chem.MolFromSmiles("OC(=O)CCCc1ccncc1")  # O 0, the pyridine N 9
    found = microstates.ensemble(mol)
    assert set(found) == {-1, 0, 1}
    for pool in found.values():
        for m in pool:
            rebuilt = Chem.RWMol(mol)
            for atom, h, q in zip(rebuilt.GetAtoms(), m.hydrogens, m.charges, strict=True):
                atom.SetNumExplicitHs(h)
                atom.SetNoImplicit(True)
                atom.SetFormalCharge(q)
            assert _canon(Chem.MolToSmiles(rebuilt)) == _canon(m.smiles)
    (anion,) = found[-1]
    assert anion.charges[0] == -1 and anion.hydrogens[0] == 0 and anion.charges[9] == 0
    zwitterion = [m for m in found[0] if m.charges[9] == 1]
    assert len(zwitterion) == 1 and zwitterion[0].charges[0] == -1


def test_conformers_follow_upstream() -> None:
    found = conformers.atoms("C[NH3+]", n_3d=3)
    assert found.symbols == ("C", "N", "H", "H", "H", "H", "H", "H")
    assert found.charges == (0, 1, 0, 0, 0, 0, 0, 0)
    assert found.coords.shape == (4, 8, 3)
    assert not np.allclose(found.coords[-1, :, 2], found.coords[0, :, 2])  # 3D, then 2D
    assert np.allclose(found.coords[-1, :, 2], 0.0)


def test_seed_0_embeds_as_seed_1_so_it_is_made_once() -> None:
    """RDKit's engine takes seed 0 as 1, which is what lets :func:`.conformers.atoms` reuse it."""
    from rdkit.Chem import AllChem

    smiles = "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1"
    made = []
    for seed in (0, 1):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
        AllChem.MMFFOptimizeMolecule(mol)
        made.append(mol.GetConformer().GetPositions())
    assert np.array_equal(made[0], made[1])
    found = conformers.atoms(smiles, n_3d=3)
    assert np.allclose(found.coords[0], made[0], atol=1e-5)
    assert np.array_equal(found.coords[0], found.coords[1])
    assert not np.allclose(found.coords[1], found.coords[2])


def test_the_network_runs_in_float32() -> None:
    """A numpy float64 scalar in the encoder would promote it to float64, 4 × slower."""
    x = np.linspace(-3, 3, 7, dtype=np.float32)
    assert model._gelu(x).dtype == np.float32
    assert model._layer_norm(x, x, x).dtype == np.float32
    assert model._softmax(x).dtype == np.float32


def test_the_2d_layout_ignores_the_global_coordgen_preference() -> None:
    """chemsuite.core.depict prefers CoordGen; the model's 2D conformer must not change with it."""
    from rdkit.Chem import rdDepictor

    smiles = "OC(=O)CC(O)(C(=O)O)CC(=O)O"
    prefer = rdDepictor.GetPreferCoordGen()
    try:
        rdDepictor.SetPreferCoordGen(False)
        plain = conformers._layout(smiles)
        rdDepictor.SetPreferCoordGen(True)
        assert np.allclose(conformers._layout(smiles), plain)
        assert rdDepictor.GetPreferCoordGen()
    finally:
        rdDepictor.SetPreferCoordGen(prefer)


def test_tokens_frame_atoms_and_charges() -> None:
    atoms, charges = model.tokens(["C", "N", "Xe"], [0, 1, 3])
    assert atoms.tolist() == [1, 4, 5, 3, 2]  # [CLS] C N [UNK] [SEP]
    assert charges.tolist() == [1, 5, 4, 3, 2]  # [CLS] 0 1 [UNK] [SEP]


def test_one_site_titrates_at_its_macro_pka() -> None:
    """Acid HA and base A⁻: half and half at the pKa, and the pKa is the loss's formula."""
    g = np.array([-1.0, 2.0])  # HA, A⁻ (kT)
    protons = np.array([1, 0])
    mean = 6.5
    pka = thermo.macro_pka(g[:1], g[1:], mean)
    assert pka == pytest.approx(mean + 3.0 / np.log(10))
    assert thermo.populations(g, protons, pka, mean) == pytest.approx([0.5, 0.5])
    assert thermo.populations(g, protons, pka + 1, mean)[1] == pytest.approx(10 / 11)


def test_common_mean_keeps_each_folds_pka():
    """Folds trained on different means agree once moved onto one: a pair's pKa is unchanged."""
    means = np.array([6.6, 6.4])
    g_acid, g_base = np.array([1.0, 1.5]), np.array([-2.0, -1.0])  # charge 1 and 0, per fold
    each = [thermo.macro_pka(g_acid[k : k + 1], g_base[k : k + 1], means[k]) for k in range(2)]
    ga = thermo.common_mean(g_acid, 1, means, 6.5)
    gb = thermo.common_mean(g_base, 0, means, 6.5)
    assert thermo.macro_pka(np.array([ga]), np.array([gb]), 6.5) == pytest.approx(np.mean(each))


def test_tautomers_share_a_macrostate() -> None:
    """Two microstates of the acid macrostate: the pKa counts both, and each keeps its share."""
    states = (
        microstates.Microstate((1, 0), (0, 0), "a1"),
        microstates.Microstate((0, 1), (0, 0), "a2"),
        microstates.Microstate((0, 0), (-1, 0), "b"),
    )
    e = Ensemble(states, np.array([0.0, np.log(3.0), 0.0]), 7.0)
    ((q, pka),) = e.macro_pkas()
    assert q == 0 and pka == pytest.approx(7.0 + np.log10(4 / 3))
    share = e.populations(pka)
    assert share[0] + share[1] == pytest.approx(0.5) and share[0] == pytest.approx(3 * share[1])
