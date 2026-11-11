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
"""Structure search at 200k rows (docs/ROADMAP.md §4.1). Run with ``pytest -m bench -s``.

The target is a substructure search under 2 s once the column's screening fingerprints are
cached (the first search computes them, like any fingerprint tool). On CI's shared runners
(``CI`` set) the limit has slack: their timing varies run to run by more than half a second.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest

from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.analysis.structure_search import StructureQuery, search
from chemsuite.core.jobs.parallel import process_pool
from chemsuite.core.mols import parse_structures

pytestmark = pytest.mark.bench

ROWS = 200_000
#: shared CI runners are slower and noisier than a developer's machine
SLACK = 1.75 if os.environ.get("CI") else 1.0
CORES = [
    "c1ccccc1",
    "c1ccncc1",
    "c1ccc2[nH]ccc2c1",
    "C1CCNCC1",
    "c1ccc2ncccc2c1",
    "c1cnc[nH]1",
    "C1CCOC1",
    "c1ccsc1",
    "O=C1CCCN1",
    "c1ccc2occc2c1",
]
GROUPS = ["C", "O", "N", "F", "Cl", "C(=O)O", "C(=O)N", "OC", "C#N", "S(=O)(=O)C", "CC", "N(C)C"]


def _library(rows: int) -> list[str]:
    """Drug-sized molecules: two cores joined by an amide, each with a substituent (some
    combinations are not valid structures; the fixture drops them)."""
    rng = np.random.default_rng(0)
    smiles = []
    for _ in range(rows):
        a, b = rng.choice(len(CORES), 2)
        g1, g2 = rng.choice(len(GROUPS), 2)
        smiles.append(f"{GROUPS[g1]}{_attach(CORES[b])}C(=O)N{_attach(CORES[a])}{GROUPS[g2]}")
    return smiles


def _attach(core: str) -> str:
    """The core with its ring-closure digits renumbered so two cores never clash."""
    return core.replace("1", "7").replace("2", "8")


@pytest.fixture(scope="module")
def library():
    smiles = _library(2 * ROWS)
    with process_pool() as pool:
        chunks = [smiles[i : i + 5000] for i in range(0, len(smiles), 5000)]
        blobs = [b for part in pool.map(parse_structures, chunks) for b in part if b][:ROWS]
        assert len(blobs) == ROWS
        start = time.perf_counter()
        fps = compute_fingerprints(blobs, "pattern", executor=pool)
        print(f"\n  pattern fingerprints for {ROWS:,} rows: {time.perf_counter() - start:.1f} s")
    assert fps.valid.all()
    return blobs, fps


@pytest.mark.parametrize(
    ("query", "limit"),
    [("c1ccc2[nH]ccc2c1", 2.0), ("C1CCOC1", 2.0), ("c1ccccc1C(=O)N", 2.0)],
)
def test_substructure_search_speed(library, query, limit) -> None:
    blobs, fps = library
    with process_pool() as pool:
        search(StructureQuery("C"), blobs[:10], fps.take(np.arange(10)), executor=pool)  # warm up
        start = time.perf_counter()
        result = search(StructureQuery(query), blobs, fps, executor=pool)
        seconds = time.perf_counter() - start
    print(
        f"\n  {query:<20} {seconds * 1000:8.0f} ms  {result.count:,} hits "
        f"of {result.candidates:,} candidates (limit {limit * SLACK * 1000:.0f} ms)"
    )
    assert result.count > 0
    assert seconds < limit * SLACK
