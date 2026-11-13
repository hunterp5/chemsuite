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
    for row in rows[1:]:
        want = {int(q): {_canon(s) for s in cell.split(",")}
                for q, cell in zip(head[1:-1], row[1:-1], strict=True) if cell}  # fmt: skip
        found = microstates.ensemble(Chem.MolFromSmiles(row[-1]))
        got = {q: {_canon(m.smiles) for m in pool} for q, pool in found.items()}
        assert got == want, row[-1]


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
