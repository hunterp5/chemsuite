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

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.analysis.structure_search import (
    QueryError,
    StructureQuery,
    match_atoms,
    parse_query,
    search,
)
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import mol_to_blob

SMILES = [
    "c1ccccc1O",  # 0 phenol
    "Cc1ccccc1O",  # 1 cresol
    "CCO",  # 2 ethanol
    "c1ccncc1",  # 3 pyridine
    "N[C@@H](C)C(=O)O",  # 4 L-alanine
    "N[C@H](C)C(=O)O",  # 5 D-alanine
    "",  # 6 blank
    "not a smiles",  # 7 invalid
    "Oc1ccccc1",  # 8 phenol again, written differently
]


def run(query: StructureQuery, cells=SMILES, **kwargs) -> list[int]:
    fps = compute_fingerprints([c or None for c in cells], query.screen_kind)
    result = search(query, cells, fps, **kwargs)
    return np.flatnonzero(result.hits).tolist()


def test_substructure_from_smiles_and_smarts() -> None:
    assert run(StructureQuery("c1ccccc1O")) == [0, 1, 8]
    assert run(StructureQuery("[OX2H][CX4]")) == [2]  # SMARTS: an aliphatic alcohol
    assert run(StructureQuery("n")) == [3]  # aromatic nitrogen (SMARTS, not valid SMILES)


def test_screen_leaves_only_plausible_candidates() -> None:
    query = StructureQuery("c1ccncc1")
    fps = compute_fingerprints([c or None for c in SMILES], query.screen_kind)
    result = search(query, SMILES, fps)
    assert result.count == 1
    assert result.candidates < len(SMILES) - 2  # the screen ruled rows out


def test_exact_match_with_and_without_stereo() -> None:
    assert run(StructureQuery("OC1=CC=CC=C1", "exact")) == [0, 8]
    assert run(StructureQuery("N[C@@H](C)C(=O)O", "exact")) == [4, 5]
    assert run(StructureQuery("N[C@@H](C)C(=O)O", "exact", stereo=True)) == [4]


def test_substructure_respecting_stereo() -> None:
    assert run(StructureQuery("N[C@@H](C)C=O")) == [4, 5]
    assert run(StructureQuery("N[C@@H](C)C=O", stereo=True)) == [4]


def test_similarity() -> None:
    query = StructureQuery("c1ccccc1O", "similarity", threshold=0.99)
    assert run(query) == [0, 8]
    fps = compute_fingerprints([c or None for c in SMILES], query.screen_kind)
    scores = search(query, SMILES, fps).scores
    assert np.isnan(scores[6]) and np.isnan(scores[7])


def test_query_errors() -> None:
    with pytest.raises(QueryError, match="Enter"):
        parse_query(StructureQuery("  "))
    with pytest.raises(QueryError, match="SMILES or SMARTS"):
        parse_query(StructureQuery("C(((("))
    with pytest.raises(QueryError, match="only for substructure"):
        parse_query(StructureQuery("[C,N]", "exact"))
    assert parse_query(StructureQuery("[C,N]")).is_smarts


def test_blobs_and_3d_cells_match_their_2d_source() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("Cc1ccccc1O"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    cells = [mol_to_blob(mol), mol_to_blob(Chem.MolFromSmiles("CCO"))]
    assert run(StructureQuery("c1ccccc1O"), cells) == [0]
    assert run(StructureQuery("Cc1ccccc1O", "exact"), cells) == [0]


def test_pool_and_chunks_give_the_same_hits() -> None:
    cells = SMILES * 50
    query = StructureQuery("c1ccccc1")
    expected = run(query, cells)
    with ThreadPoolExecutor(2) as pool:
        import chemsuite.analysis.structure_search as module

        old = module.POOL_MIN_ROWS
        module.POOL_MIN_ROWS = 1
        try:
            assert run(query, cells, executor=pool, chunk_rows=7) == expected
        finally:
            module.POOL_MIN_ROWS = old


def test_cancel() -> None:
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        run(StructureQuery("c1ccccc1"), cancel=token)


def test_match_atoms_covers_every_match() -> None:
    mol = Chem.MolFromSmiles("Oc1ccc(O)cc1")
    atoms, bonds = match_atoms(mol, Chem.MolFromSmarts("[OH]c"))
    assert atoms == [0, 1, 4, 5]
    assert len(bonds) == 2
