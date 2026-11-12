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
"""Cores, analog series and the constellation they make (the CCR method)."""

from __future__ import annotations

import time

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.sar.constellation import (
    constellation,
    find_cores,
    place,
    putative_cores,
    summarise,
)

#: a benzimidazole aryl ketone; the amides vary at the amide nitrogen
BENZIMIDAZOLE = "Cn1c(-c2ccc(C(={})cc2)nc2ccccc21"
ACYL = "Cn1c(-c2ccc(C=O)cc2)nc2ccccc21"  # the shared core, its attachment point an H
AMIDES = [BENZIMIDAZOLE.format(f"O){r}") for r in ("NC", "NCC", "NC3CC3", "N(C)C")]


def cores(smiles: str) -> list[str]:
    return putative_cores(Chem.MolFromSmiles(smiles))


def test_a_molecule_is_its_own_first_core() -> None:
    assert cores("c1ccccc1") == ["c1ccccc1"]
    found = cores(AMIDES[0])
    assert found[0] == Chem.CanonSmiles(AMIDES[0]) and ACYL in found


def test_cores_keep_two_thirds_of_the_heavy_atoms() -> None:
    # paracetamol: cutting the amide leaves 4-aminophenol (8 of 11 atoms), never the acetyl
    assert cores("CC(=O)Nc1ccc(O)cc1") == ["CC(=O)Nc1ccc(O)cc1", "Nc1ccc(O)cc1"]
    # a benzamide cut in half keeps neither half (8 of 16 atoms each)
    assert cores("O=C(NCc1ccncc1)c1ccccc1") == ["O=C(NCc1ccncc1)c1ccccc1"]
    loose = putative_cores(Chem.MolFromSmiles("O=C(NCc1ccncc1)c1ccccc1"), fraction=0.4)
    assert "NCc1ccncc1" in loose and "O=Cc1ccccc1" in loose


def test_salts_and_hydrogens_do_not_count() -> None:
    assert cores("[Na+].O=C([O-])c1ccccc1") == ["O=C([O-])c1ccccc1"]
    assert cores("[H]OC([H])([H])c1ccccc1") == ["OCc1ccccc1"]


def test_long_chains_stay_quick() -> None:
    start = time.perf_counter()
    found = cores("CC(N)C(=O)" + "NC(C)C(=O)" * 30 + "O")
    assert len(found) > 1 and time.perf_counter() - start < 10


def test_analogs_share_a_core_and_form_a_series() -> None:
    other = "O=C(NCc1ccncc1)c1ccccc1"
    lists = find_cores([*AMIDES, other, None, "not a molecule"])
    assert lists[-2:] == [(), ()]
    con = constellation(lists)
    assert con.cores == [ACYL] and con.n.tolist() == [4]
    assert con.members[0].tolist() == [0, 1, 2, 3] and con.series.tolist() == [1]
    assert con.series_size == {1: 4} and con.compounds == 5 and len(con.edges) == 0


def test_identical_member_sets_keep_the_largest_core() -> None:
    lists = [
        ("CCCCCCN", "CCCCN", "CCCCCCNC(C)=O"),
        ("CCCCCCN", "CCCCN", "CCCCCCNC(=O)CC"),
        ("c1ccccc1",),
    ]
    con = constellation(lists)
    assert con.cores == ["CCCCCCN"]  # CCCCN holds the same two compounds and is smaller


def test_series_edges_and_the_size_threshold() -> None:
    # compounds 0-2 share A, 2-3 share B: one series of four joined through compound 2
    lists = [("CCCCCCCO",), ("CCCCCCCO",), ("CCCCCCCO", "c1ccncc1"), ("c1ccncc1",), ("CCl",)]
    con = constellation(lists)
    assert sorted(con.cores) == ["CCCCCCCO", "c1ccncc1"]
    assert con.series.tolist() == [1, 1] and con.series_size == {1: 4}
    assert con.edges.tolist() == [[0, 1]]
    assert constellation(lists, min_compounds=3).cores == ["CCCCCCCO"]
    assert len(constellation([("C",), ("O",)])) == 0


def test_summaries_leave_cores_without_values_blank() -> None:
    con = constellation([("CCCCO",), ("CCCCO",), ("CCCCO",), ("CCCCN",), ("CCCCN",)])
    values = np.array([5.0, 7.0, 9.0, np.nan, np.nan])
    order = {core: i for i, core in enumerate(con.cores)}
    mean = summarise(con, values, "mean")
    assert mean[order["CCCCO"]] == 7.0 and np.isnan(mean[order["CCCCN"]])
    assert summarise(con, values, "count")[order["CCCCN"]] == 0
    assert summarise(con, values, "max")[order["CCCCO"]] == 9.0
    with pytest.raises(ValueError):
        summarise(con, values, "mode")


def test_cores_are_placed_on_a_map() -> None:
    few = constellation([("CCCCO",), ("CCCCO",)])
    assert place(few).shape == (1, 2)
    smiles = ["c1ccccc1O", "c1ccccc1N", "c1ccncc1", "C1CCCCC1", "CCOC(=O)C", "c1ccc2ccccc2c1"]
    con = constellation([(s,) for s in smiles for _ in range(2)])
    for method in ("pca", "tsne"):
        coords = place(con, method)
        assert coords.shape == (6, 2) and np.isfinite(coords).all()
    assert place(constellation([(s,) for s in smiles[:3] for _ in range(2)])).shape == (3, 2)


def test_finding_cores_can_be_cancelled() -> None:
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        find_cores(AMIDES * 100, cancel=token)
