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
"""Simulation M5a: Prepare + Minimise on the P0 neutron set (docs/SIMULATE_PLAN.md, M5).

Reads what ``scripts/review/minimise_review.py run`` wrote under ``samples/prepare/out``
(about an hour on the OpenMM worker); skipped without it. Run with ``pytest -m bench -s``.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT / "scripts" / "review"
NUMBERS = ROOT / "samples" / "prepare" / "out" / "minimise"

pytestmark = [
    pytest.mark.bench,
    pytest.mark.skipif(not NUMBERS.is_dir(), reason="no minimise_review.py run to read"),
]

RMSD_LIMIT = 0.3  # Å, MinimiseSettings' default
#: entries that fail for a defect (``minimise_review.DEFECTS``; the rest of the failures are
#: what chemsuite does not simulate by design). A ratchet: lower it as defects are fixed
DEFECTS = 14
#: the share of minimised entries whose held heavy atoms end within the limit
WITHIN_LIMIT = 0.95
#: how far a neutron mark may fall from the prepared model's (minimisation keeps the states;
#: an Asn/Gln/His may still swing past the half-way point)
MARK_SLACK = 0.01
#: the share of OH/SH hydrogens turned away from the deposit's may not grow by more than this
TURN_SLACK = 0.01
#: surface lysines swing in implicit solvent: 0.9% of NH3 hydrogens turned after Prepare,
#: 8.9% after Minimise (2026-11-24); a ceiling, so it does not quietly get worse
LYS_TURNED = 0.10
#: Å a transition metal's contact may change by (held since 2026-11-24; a bare ion pulled them
#: in by up to 2.75 Å, median 0.34)
METAL_CONTACT_SLACK = 0.15


@pytest.fixture(scope="module")
def review():
    sys.path.insert(0, str(REVIEW))
    try:
        spec = importlib.util.spec_from_file_location(
            "minimise_review", REVIEW / "minimise_review.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(REVIEW))
    return module


@pytest.fixture(scope="module")
def rows(review) -> list[dict]:
    return review.load()


def _done(rows: list[dict]) -> list[dict]:
    return [r for r in rows if "result" in r]


def test_no_more_entries_fail_for_a_defect(review, rows) -> None:
    failed = [r for r in rows if "failed" in r]
    defects = [(r["id"], why) for r in failed for kind, why in [review.failure(r)]
               if kind == "defect"]  # fmt: skip
    print(f"{len(_done(rows))} of {len(rows)} minimised; {len(defects)} defects: {defects}")
    assert len(defects) <= DEFECTS


def test_the_held_heavy_atoms_stay_within_the_limit(rows) -> None:
    rmsd = np.array([r["result"]["rmsd"] for r in _done(rows)])
    print(f"heavy-atom RMSD median {np.median(rmsd):.2f}, max {rmsd.max():.2f} Å")
    assert np.mean(rmsd <= RMSD_LIMIT) >= WITHIN_LIMIT


def test_minimisation_lowers_the_clash_score(rows) -> None:
    done = _done(rows)
    before = np.array([r["clashes_before"] / r["atoms"] for r in done])
    after = np.array([r["clashes_after"] / r["atoms"] for r in done])
    print(f"clashes per 1000 atoms: {1000 * before.mean():.2f} -> {1000 * after.mean():.2f}")
    assert after.mean() < before.mean()
    assert np.mean(after <= before) >= WITHIN_LIMIT


def test_the_neutron_marks_hold(review, rows) -> None:
    ids = [r["id"] for r in _done(rows)]
    prepared = review.marks(review.PREP, ids)
    minimised = review.marks(review.MIN, ids)
    for label, tally in prepared.items():
        if not tally.total:
            continue
        before = tally.right / tally.total
        after_tally = minimised.get(label)
        assert after_tally is not None and after_tally.total, label
        after = after_tally.right / after_tally.total
        print(f"{label}: {before:.1%} -> {after:.1%}")
        assert after >= before - MARK_SLACK, label


def test_minimisation_does_not_turn_the_polar_hydrogens_away(review, rows) -> None:
    turned = {}
    for group, method, _distances, angles in review.h_table(_done(rows)):
        turned[group, method] = float(np.mean(angles > review.TURNED))
    for group in review.GROUPS:
        print(
            f"{group}: {turned[group, review.PREP]:.1%} -> {turned[group, review.MIN]:.1%} turned"
        )
    rotor = "OH/SH rotor"
    assert turned[rotor, review.MIN] <= turned[rotor, review.PREP] + TURN_SLACK
    assert turned["Lys NH3", review.MIN] <= LYS_TURNED


def test_the_metal_sites_hold(review, rows) -> None:
    sites = review.metal_sites(_done(rows))
    changes = np.array([change for _id, _metal, change, _came in sites])
    print(f"{len(sites)} metals; largest contact change {changes.max():.2f} Å")
    assert sites and changes.max() <= METAL_CONTACT_SLACK
    assert [(i, metal, came) for i, metal, _c, came in sites if came] == []
