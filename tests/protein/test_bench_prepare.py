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
"""Protein preparation against neutron structures (docs/prepare/p0-review.md).

Needs the benchmark data under ``samples/prepare`` (``scripts/review/prepare_review.py``
builds it); skipped without it. Run with ``pytest -m bench -s``.

P0 checks the harness itself: a method that flips nothing scores about half on orientation
(the inputs flip half the side chains at random), and Reduce keeps the scores recorded in the
review. Each later phase adds our own method here with the review's acceptance thresholds
(P1: chemsuite's flips, His states and hydrogens; P2: pKa values on PKAD-R's test half, at
PROPKA's accuracy, from ``scripts/review/pkad_review.py chemsuite``).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT / "scripts" / "review"
DATA = ROOT / "samples" / "prepare"

pytestmark = [
    pytest.mark.bench,
    pytest.mark.skipif(
        not (DATA / "neutron" / "ids.txt").exists(), reason="no samples/prepare benchmark data"
    ),
]

ORIENTATIONS = ("ASN orientation", "GLN orientation", "HIS orientation")
#: Reduce on the 262 entries (p0-review.md §3), less a point of headroom
REDUCE = {
    "ASN orientation": 0.84,
    "GLN orientation": 0.78,
    "HIS orientation": 0.87,
    "His tautomer (HID/HIE/HIP)": 0.61,
}


@pytest.fixture(scope="module")
def review():
    sys.path.insert(0, str(REVIEW))
    try:
        spec = importlib.util.spec_from_file_location(
            "prepare_review", REVIEW / "prepare_review.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(REVIEW))
    return module


def share(tally) -> float:
    return tally.right / tally.total


def test_a_method_that_never_flips_scores_half(review) -> None:
    covered, tallies = review.method_scores("naive")
    assert covered == len(review.entries())
    for question in ORIENTATIONS:
        print(f"naive {question}: {share(tallies[question]):.1%}")
        assert 0.45 <= share(tallies[question]) <= 0.55


def test_reduce_keeps_its_recorded_scores(review) -> None:
    if not (DATA / "out" / "reduce").is_dir():
        pytest.skip("no Reduce outputs")
    _, tallies = review.method_scores("reduce")
    for question, floor in REDUCE.items():
        print(f"reduce {question}: {share(tallies[question]):.1%}")
        assert share(tallies[question]) >= floor


#: P1's pass marks (p0-review.md §4): Reduce's scores, and the textbook states kept
P1 = {
    "ASN orientation": 0.85,
    "GLN orientation": 0.79,
    "HIS orientation": 0.88,
    "His tautomer (HID/HIE/HIP)": 0.62,
    "ASP state": 0.989,  # the textbook state: 3916 of 3956
    "GLU state": 0.98,
    "Ligand polar-atom H count": 0.955,
}


def test_chemsuite_meets_the_p1_marks(review) -> None:
    """Scores the outputs ``prepare_review.py chemsuite`` wrote (re-run it with ``--redo``
    after changing the pipeline)."""
    if not (DATA / "out" / "chemsuite").is_dir():
        pytest.skip("no chemsuite outputs")
    covered, tallies = review.method_scores("chemsuite")
    assert covered == len(review.entries())
    for question, floor in P1.items():
        print(f"chemsuite {question}: {share(tallies[question]):.1%}")
        assert share(tallies[question]) >= floor, question


def _load(name: str):
    sys.path.insert(0, str(REVIEW))
    try:
        spec = importlib.util.spec_from_file_location(name, REVIEW / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module  # its dataclasses look themselves up there
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(REVIEW))
    return module


def test_chemsuite_pka_values_match_propka() -> None:
    """P2's mark (docs/prepare/p2-review.md): on PKAD-R's test half, chemsuite's RMSE is no
    worse than PROPKA's on the same residues (the null model, every residue at its model pKa,
    is not yet beaten by either)."""
    pkad = _load("pkad_review")
    if not (pkad.DATA / "chemsuite.csv").is_file() or not (pkad.DATA / "propka").is_dir():
        pytest.skip("no PKAD-R outputs")
    found = {m: pkad.predictions(m) for m in ("null", "propka", "chemsuite")}
    shared = [e for e in pkad.entries() if all(e in f for f in found.values())]
    test = [e for e in shared if not pkad.training(e)]
    rmse = {m: float(np.sqrt(np.mean([(f[e] - e.pka) ** 2 for e in test])))
            for m, f in found.items()}  # fmt: skip
    print(f"PKAD-R test half (n={len(test)}):", {m: round(v, 2) for m, v in rmse.items()})
    assert len(test) > 250
    assert rmse["chemsuite"] <= rmse["propka"] + 0.005
