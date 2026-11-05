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

"""File and document benchmarks (``pytest -m bench -s``).

CI runners have 2–4 cores, so this uses 20k records and checks throughput. On an 8-core laptop a
200k-record, 589 MB SDF loads in about 9 s (target: 30 s); see docs/PLAN.md §2.6.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from chemsuite.core.document import load_document, save_document
from chemsuite.core.io.load import load_table, parse_pending_structures, save_table
from chemsuite.core.jobs.parallel import process_pool

from .conftest import DRUGS, molblock, sdf_record

RECORDS = 20_000
pytestmark = pytest.mark.bench


@pytest.fixture(scope="module")
def big_sdf(tmp_path_factory) -> Path:
    blocks = [molblock(smiles, title=name) for name, smiles in DRUGS.items()]
    path = tmp_path_factory.mktemp("bench") / "big.sdf"
    with path.open("w", encoding="utf-8") as handle:
        for i in range(RECORDS):
            handle.write(
                sdf_record(
                    blocks[i % len(blocks)],
                    {
                        "ID": f"CPD-{i:06d}",
                        "pIC50": f"{5 + (i % 400) / 100:.2f}",
                        "Series": "AB"[i % 2],
                    },
                )
            )
    return path


def timed(label: str, seconds: float, limit: float, rows: int) -> None:
    print(f"\n  {label:<28} {seconds:6.2f} s  ({rows / seconds:,.0f} rows/s; limit {limit:.0f} s)")
    assert seconds < limit


def test_load_sdf_parallel(big_sdf: Path) -> None:
    with process_pool() as pool:
        start = time.perf_counter()
        table = load_table(big_sdf, executor=pool)
        timed("load 20k-record SDF", time.perf_counter() - start, 10.0, len(table))
    assert len(table) == RECORDS


def test_csv_load_and_parse(big_sdf: Path, tmp_path: Path) -> None:
    csv_path = tmp_path / "big.csv"
    save_table(load_table(big_sdf), csv_path)
    start = time.perf_counter()
    table = load_table(csv_path, structure_column="Structure")
    timed("load 20k-row CSV (text)", time.perf_counter() - start, 2.0, len(table))
    with process_pool() as pool:
        start = time.perf_counter()
        parse_pending_structures(table, executor=pool)
        timed("parse 20k SMILES", time.perf_counter() - start, 10.0, len(table))


def test_document_save_and_open(big_sdf: Path, tmp_path: Path) -> None:
    table = load_table(big_sdf)
    path = tmp_path / "big.csd"
    start = time.perf_counter()
    save_document(table, path)
    timed("save .csd", time.perf_counter() - start, 3.0, len(table))
    start = time.perf_counter()
    loaded, _ = load_document(path)
    timed("open .csd", time.perf_counter() - start, 1.5, len(loaded))
