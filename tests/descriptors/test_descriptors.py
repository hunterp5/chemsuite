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
    from chemsuite.conformers.etkdg import generate
    from chemsuite.core.mols import mol_from_blob, mol_to_blob

    three_d = generate(mol_from_blob(parse_structure(ASPIRIN)[0])).mol  # for the 3D shapes
    out = compute_chunk((tuple(CATALOG), [mol_to_blob(three_d)]))
    blank = []
    for key, (value,) in out.items():
        if CATALOG[key].is_text:
            ok = isinstance(value, str) and bool(value)
        else:
            ok = math.isfinite(value) or key == "PKA_BASIC"  # aspirin has no basic site
        if not ok:
            blank.append(f"{key}: {value!r} ({_why(key, three_d)})")
    assert not blank, "\n".join(blank)


def _why(key: str, three_d: Chem.Mol) -> str:
    """What computing *key* on its own gives or raises, to say why a descriptor came out
    blank."""
    from chemsuite.descriptors.compute import SHAPE_GROUP, _as_drawn

    spec = CATALOG[key]
    try:
        return repr(spec.fn(three_d if spec.group == SHAPE_GROUP else _as_drawn(three_d)))
    except Exception as error:  # noqa: BLE001 (any error is the answer here)
        return f"{type(error).__name__}: {error}"


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


def test_sa_score() -> None:
    assert medchem.sa_score(mol(ASPIRIN)) == pytest.approx(1.58, abs=0.01)  # Ertl's paper: 1.58
    assert medchem.sa_score(mol("c1ccccc1")) == pytest.approx(1.0)
    erythromycin = (
        "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)"
        "[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)"
        "[C@@H](O)[C@]1(C)O"
    )
    assert medchem.sa_score(mol(erythromycin)) > 5  # a complex natural product
    assert CATALOG["SA_SCORE"].title == "SAScore" and CATALOG["SA_SCORE"].group == "Drug-likeness"


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


def test_cancel_keeps_the_rows_already_done() -> None:
    token = CancelToken()

    def progress(done: int, total: int) -> None:
        if done >= 20:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        compute_descriptors(["CCO"] * 50, ["MolWt", "SMILES"], chunk_rows=10,
                            progress=progress, cancel=token)  # fmt: skip
    partial = stopped.value.partial
    assert len(partial["MolWt"]) == len(partial["SMILES"]) == 20
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled) as stopped:
        compute_descriptors(["CCO"] * 5, ["MolWt"], cancel=token)
    assert stopped.value.partial is None  # nothing finished


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
