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

import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import Crippen

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import parse_structure
from chemsuite.descriptors import medchem
from chemsuite.descriptors.catalog import CATALOG, COMMON, GROUPS, in_group
from chemsuite.descriptors.compute import compute_chunk, compute_descriptors

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def mol(smiles: str) -> Chem.Mol:
    return Chem.MolFromSmiles(smiles)


def test_catalog_is_consistent() -> None:
    titles = [d.title for d in CATALOG.values()]
    assert len(titles) == len(set(titles)), "column titles must be unique"
    assert set(COMMON) <= set(CATALOG)
    assert sum(len(in_group(g)) for g in GROUPS) == len(CATALOG)


def test_every_descriptor_runs_on_a_drug() -> None:
    out = compute_chunk((tuple(CATALOG), [parse_structure(ASPIRIN)[0]]))
    for key, (value,) in out.items():
        spec = CATALOG[key]
        if spec.is_text:
            assert isinstance(value, str) and value, key
        else:
            assert math.isfinite(value), key


def test_known_values() -> None:
    aspirin = mol(ASPIRIN)
    assert medchem.lipinski_violations(aspirin) == 0
    assert medchem.ro5_pass(aspirin) == "Yes"
    assert medchem.lipinski_violations(mol("C" * 40)) == 2  # MW > 500 and cLogP > 5
    assert medchem.ro5_pass(mol("C" * 40)) == "No"
    assert medchem.esol_logs(mol("c1ccccc1")) == pytest.approx(-2.0, abs=0.3)  # exp. -1.64
    assert 0.0 <= medchem.cns_mpo(aspirin) <= 6.0
    assert 0.0 <= medchem.bbb_score(aspirin) <= 6.0
    assert medchem.net_formal_charge(mol("[NH4+].[Cl-]")) == 0
    assert medchem.net_formal_charge(mol("C[N+](C)(C)C")) == 1
    # A base with pKa above 7.4 has cLogD below cLogP.
    amine = mol("CCCCN(C)C")
    assert medchem.clogd74(amine) < Crippen.MolLogP(amine)


def test_bbb_score_terms() -> None:
    # Outside every polynomial window only the aromatic-ring term is left.
    assert medchem.bbb_score_from_descriptors(2, 0, 100, 0, 0, 200, 1) == pytest.approx(1.0)


def test_blank_invalid_and_unparsed_cells() -> None:
    cells = [parse_structure("CCO")[0], b"", None, "c1ccccc1", "not smiles"]
    columns = compute_descriptors(cells, ["MolWt", "FORMULA"])
    assert columns["MolWt"].displays(range(5)) == ["46.07", "", "", "78.11", ""]
    assert columns["FORMULA"].displays(range(5)) == ["C2H6O", "", "", "C6H6", ""]
    assert columns["MolWt"].kind == "number" and columns["FORMULA"].kind == "text"


def test_integer_counts_and_no_negative_zero() -> None:
    columns = compute_descriptors(["CCO"], ["NumHDonors", "MolLogP"])
    assert columns["NumHDonors"].display(0) == "1"
    assert columns["MolLogP"].display(0) == "0.00"  # -0.0014 rounds to 0, never "-0.00"


def test_chunks_progress_and_cancel() -> None:
    cells = ["CCO"] * 25
    seen: list[tuple[int, int]] = []
    columns = compute_descriptors(
        cells, ["MolWt"], chunk_rows=10, progress=lambda d, t: seen.append((d, t))
    )
    assert len(columns["MolWt"]) == 25
    assert seen == [(10, 25), (20, 25), (25, 25)]
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        compute_descriptors(cells, ["MolWt"], chunk_rows=10, cancel=token)


def test_unknown_key() -> None:
    with pytest.raises(KeyError, match="nope"):
        compute_descriptors(["C"], ["nope"])


def test_worker_processes_give_the_same_result() -> None:
    from chemsuite.core.jobs.parallel import process_pool

    cells = [parse_structure(s)[0] for s in (ASPIRIN, "CCN", "c1ccncc1")] * 4
    inline = compute_descriptors(cells, ["MolWt", "INCHIKEY"], chunk_rows=5)
    with process_pool(2) as pool:
        pooled = compute_descriptors(cells, ["MolWt", "INCHIKEY"], chunk_rows=5, executor=pool)
    for key in inline:
        rows = np.arange(len(cells))
        assert inline[key].displays(rows) == pooled[key].displays(rows)
