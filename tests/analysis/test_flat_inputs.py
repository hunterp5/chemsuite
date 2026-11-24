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

"""2D work on a 3D column reads each cell's SMILES instead of loading its every conformer
(``StructureColumn.flat_inputs``, ROADMAP §7.0f). These pin that every tool moved to it gives
the same answer from the text as from the blob it read before."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.analysis.fingerprints import KINDS, compute_counts, compute_fingerprints
from chemsuite.analysis.graphs.mhfp import compute_mhfp
from chemsuite.analysis.structure_search import StructureQuery, search
from chemsuite.core.mols import mol_to_blob, smiles_of
from chemsuite.core.table.columns import FLAT_MIN_BYTES, StructureColumn
from chemsuite.descriptors.catalog import CATALOG, PKA_KEYS, SHAPE_GROUP
from chemsuite.descriptors.compute import compute_descriptors
from chemsuite.sar.constellation import cores_chunk

SOURCES = [
    "C[C@H](N)C(=O)O",  # a defined stereocenter
    "CC(N)C(=O)O",  # an undefined one
    "C/C=C/C(=O)[O-]",  # E double bond and a charge
    "C[N+](C)(C)CCO",  # a quaternary nitrogen
    "[2H]C([2H])([2H])c1ccccc1",  # deuterium
    "[13CH3]c1ccccc1O",  # carbon-13
    "C/C=N/O",  # oxime geometry
    "[H]/N=C(\\C)c1ccccc1",  # imine geometry set by its hydrogen
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
]


def three_d(smiles: str, conformers: int = 60) -> tuple[bytes, str]:
    """A 3D cell as Generate Conformers writes one: the molecule with its hydrogens and
    conformers as the blob, canonical SMILES without the hydrogens as the text."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, conformers, randomSeed=7)
    return mol_to_blob(mol), smiles_of(Chem.RemoveHs(mol))


@pytest.fixture(scope="module")
def cells() -> tuple[list, list]:
    """``(blob inputs, flat inputs)`` of one 3D column, plus a blank cell."""
    blobs, texts = zip(*(three_d(s) for s in SOURCES), strict=True)
    column = StructureColumn.from_parsed([*blobs, b""], [*texts, ""])
    slots = np.arange(len(column))
    before, after = column.worker_inputs(slots), column.flat_inputs(slots)
    assert all(len(b) > FLAT_MIN_BYTES for b in blobs)  # every cell takes the text route
    assert after[:-1] == list(texts) and after[-1] is None
    return before, after


def test_small_blobs_molfiles_and_smarts_keep_their_blobs() -> None:
    small, small_text = three_d("CCO", conformers=1)
    big, _text = three_d("CCO")
    assert len(small) <= FLAT_MIN_BYTES < len(big)
    block = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
    texts = [small_text, block, "[#6]-[#8]", "[C;X4]O", "C~O", "CCO", ""]
    column = StructureColumn.from_parsed([small, big, big, big, big, b"", b""], texts)
    slots = np.arange(len(texts))
    assert column.flat_inputs(slots) == column.worker_inputs(slots)  # none handed as text
    assert StructureColumn.from_texts(["CCN"]).flat_inputs(np.arange(1)) == ["CCN"]  # unparsed


def test_fingerprints_are_the_same(cells) -> None:
    before, after = cells
    for key in KINDS:
        old, new = compute_fingerprints(before, key), compute_fingerprints(after, key)
        assert np.array_equal(old.valid, new.valid), key
        assert np.array_equal(old.packed, new.packed), key
    old, new = compute_counts(before, "morgan2"), compute_counts(after, "morgan2")
    assert np.array_equal(old[0], new[0])


def test_mhfp_is_the_same(cells) -> None:
    before, after = cells
    assert np.array_equal(compute_mhfp(before).hashes, compute_mhfp(after).hashes)


@pytest.mark.parametrize(
    "query",
    [
        StructureQuery("c1ccccc1"),
        StructureQuery("[2H]"),
        StructureQuery("[13C]"),
        StructureQuery("[#1]"),
        StructureQuery("[NH]"),
        StructureQuery("C=N"),
        StructureQuery("[O-]"),
        StructureQuery("[N+]"),
        StructureQuery("C/C=N/O", stereo=True),
        StructureQuery("C[C@@H](N)C(=O)O", stereo=True),
        *(StructureQuery(s, "exact", stereo=True) for s in SOURCES),
        *(StructureQuery(s, "exact") for s in SOURCES),
    ],
    ids=lambda q: f"{q.mode}:{q.text}:{q.stereo}",
)
def test_searches_find_the_same_rows(cells, query) -> None:
    before, after = cells
    hits = [
        search(query, c, compute_fingerprints(c, query.screen_kind)).hits for c in (before, after)
    ]
    assert np.array_equal(*hits)


def test_2d_descriptors_are_the_same(cells) -> None:
    before, after = cells
    keys = [k for k, d in CATALOG.items() if d.group != SHAPE_GROUP and k not in PKA_KEYS]
    old, new = compute_descriptors(before, keys), compute_descriptors(after, keys)
    for key in keys:
        a, b = old[key].values(), new[key].values()
        if CATALOG[key].is_text:
            assert list(a) == list(b), key
        else:  # sums over atoms in another order may differ in the last digit
            a, b = np.array(a, float), np.array(b, float)
            assert np.allclose(a, b, rtol=1e-9, equal_nan=True), key


def test_constellation_cores_are_the_same(cells) -> None:
    before, after = cells
    assert cores_chunk((0.67, before)) == cores_chunk((0.67, after))


def test_pka_is_the_same() -> None:
    """The pKa model makes its own conformers from SMILES, so a cell's do not reach it."""
    cells = [three_d(s) for s in ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "CN(C)CCCc1ccccc1")]
    keys = ["PKA_ACIDIC", "PKA_BASIC"]
    old = compute_descriptors([blob for blob, _ in cells], keys)
    new = compute_descriptors([text for _, text in cells], keys)
    for key in keys:
        a, b = np.array(old[key].values(), float), np.array(new[key].values(), float)
        assert np.allclose(a, b, equal_nan=True), key
