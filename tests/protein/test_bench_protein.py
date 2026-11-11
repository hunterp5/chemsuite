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
"""Protein benchmarks on a 100,000-atom structure (63 copies of 3HTB side by side, each with
its own chains). Run with ``pytest -m bench -s``.

Timings are the best of a few runs; the limits leave headroom over a laptop's (M-series:
reading 8 MB of mmCIF 0.8 s, a selection's new look 13 ms, the cartoon 0.13 s).
"""

from __future__ import annotations

import gc
import subprocess
import sys
import time
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from chemsuite.protein.chem import poses
from chemsuite.protein.look import cartoon
from chemsuite.protein.look import display as displays
from chemsuite.protein.structure import io
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.structure.select import select

from .conftest import load

ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.bench

ATOMS = 100_000
READ_S = 1.5  # the whole file, bonds and secondary structure included
SELECTION_S = 0.025  # a selection change: the new look of every atom and the cartoon
SELECT_S = 0.05  # a distance expression over every atom
CARTOON_S = 0.3
COMPOSE_S = 0.5


def best_of(fn: Callable[[], object], repeat: int = 3) -> float:
    """The best time of ``repeat`` runs, the garbage collector paused (in the full suite it
    would sweep the other tests' heap while a read allocates)."""
    best = float("inf")
    gc.collect()
    gc.disable()
    try:
        for _ in range(repeat):
            start = time.perf_counter()
            fn()
            best = min(best, time.perf_counter() - start)
    finally:
        gc.enable()
    return best


def tiled(base: Structure, atoms: int) -> Structure:
    """Copies of ``base`` on a grid 60 Å apart until there are ``atoms``; each copy's chains
    get their own IDs, as a large complex's have."""
    big, k = base, 1
    while big.n_atoms < atoms:
        shift = np.array([(k % 8) * 60.0, (k // 8) * 60.0, 0.0], np.float32)
        tile = replace(
            base.with_coords(base.coords + shift),
            chain_name=np.array([f"{c}{k}" for c in base.chain_name], dtype="<U4"),
            label_asym=np.array([f"{c}{k}" for c in base.label_asym], dtype="<U4"),
        )
        big = poses.merged(big, tile)
        k += 1
    return big


@pytest.fixture(scope="module")
def big() -> Structure:
    return tiled(load("3htb.cif.gz"), ATOMS)


READ = """
import sys, time
from chemsuite.protein.structure import io
best = float("inf")
for _ in range(2):
    start = time.perf_counter()
    io.read(sys.argv[1])
    best = min(best, time.perf_counter() - start)
print(best)
"""


def test_reading_a_large_mmcif_file(big, tmp_path) -> None:
    """In a process of its own: the rest of a test run (200,000-row tables held by the table
    benchmarks) slows allocation enough to triple the time."""
    path = tmp_path / "big.cif"
    path.write_text(io.to_text(big, io.MMCIF))
    out = subprocess.run([sys.executable, "-c", READ, str(path)], capture_output=True,
                         text=True, timeout=120, cwd=ROOT)  # fmt: skip
    assert out.returncode == 0, out.stderr
    seconds = float(out.stdout.split()[-1])
    print(f"\nread {path.stat().st_size >> 20} MB, {big.n_atoms:,} atoms: {seconds:.2f} s")
    assert seconds < READ_S


def test_a_selection_change_is_drawn_within_a_frame(big) -> None:
    shown = displays.Display("big", "big", big)
    scene = displays.compose([shown])
    shown.selected = select(big, "within 5 of ligand")

    def change() -> None:
        shown.selected = ~shown.selected
        displays.relook(scene)

    seconds = best_of(change, repeat=5)
    print(f"\nselection look, {big.n_atoms:,} atoms: {seconds * 1000:.1f} ms")
    assert seconds < SELECTION_S


def test_selecting_and_building(big) -> None:
    shown = displays.Display("big", "big", big)
    timings = {
        "select": (best_of(lambda: select(big, "byres (within 5 of ligand)")), SELECT_S),
        "cartoon": (
            best_of(
                lambda: cartoon.build(big, np.arange(big.n_residues), 0, shown.active), repeat=2
            ),
            CARTOON_S,
        ),  # fmt: skip
        "compose": (best_of(lambda: displays.compose([shown]), repeat=2), COMPOSE_S),
    }
    print("\n" + ", ".join(f"{k} {s * 1000:.0f} ms" for k, (s, _) in timings.items()))
    assert all(seconds < limit for seconds, limit in timings.values())
