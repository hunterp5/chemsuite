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
"""The preparation's finer settings (:class:`pipeline.Options`): each one changes what comes
out, and the defaults change nothing."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from chemsuite.protein.prepare import inspect, pipeline, repair, review, sites
from chemsuite.protein.structure.model import WATER
from tests.protein.prepare.test_hydrogens import _residue
from tests.protein.prepare.test_pipeline import _without

DRY = pipeline.Options(waters=pipeline.NONE)


def _replace(**changes) -> pipeline.Options:
    return replace(DRY, **changes)


def test_the_options_name_what_differs_from_the_defaults() -> None:
    assert pipeline.Options().changed() == {}
    assert pipeline.Options(caps=False, his="HID").changed() == {"caps": False, "his": "HID"}


def test_cut_ends_can_be_left_uncapped(lysozyme) -> None:
    done = pipeline.prepare(lysozyme, _replace(caps=False))
    s = done.structure
    assert not np.isin(s.res_name, list(inspect.CAPS)).any()
    assert not [f for f in done.report if f.kind == inspect.CAP]
    end = _residue(s, "A/ASN 163")
    assert "OXT" not in s.name[s.atoms_of(end)]  # left as the file has it, not made a terminus


def test_a_lone_missing_residue_can_be_left_out(egfr) -> None:
    built = pipeline.prepare(egfr, DRY)
    left = pipeline.prepare(egfr, _replace(fill_gaps=False))
    assert [f for f in built.report if f.kind == inspect.FILLED]
    assert not [f for f in left.report if f.kind == inspect.FILLED]
    caps = [sum(f.kind == inspect.CAP for f in d.report) for d in (built, left)]
    assert caps[1] > caps[0]  # its ends are capped instead
    s = left.structure
    assert "A/ILE 706" not in {s.residue_label(r) for r in range(s.n_residues)}


def test_side_chains_can_be_left_as_the_file_has_them(crambin, egfr) -> None:
    cut = _without(crambin, "A/ARG 10", {"CD", "NE", "CZ", "NH1", "NH2"})
    done = pipeline.prepare(cut, _replace(side_chains=False))
    s = done.structure
    assert "CZ" not in s.name[s.atoms_of(_residue(s, "A/ARG 10"))]
    assert "CZ" in pipeline.prepare(cut, DRY).structure.name  # the default builds it
    # a residue built from nothing still gets its side chain
    lone = pipeline.prepare(egfr, _replace(side_chains=False)).structure
    assert "CB" in lone.name[lone.atoms_of(_residue(lone, "A/ILE 706"))]


def test_the_first_alternate_location_can_be_kept(lysozyme) -> None:
    # one residue's B location made the more occupied, so the two rules part
    r = min(lysozyme.default_altlocs)
    occupancy = lysozyme.occupancy.copy()
    occupancy[(lysozyme.residue_of == r) & (lysozyme.altloc == "B")] = 0.9
    occupancy[(lysozyme.residue_of == r) & (lysozyme.altloc == "A")] = 0.1
    lysozyme = replace(lysozyme, occupancy=occupancy)
    first = pipeline.default_altlocs(lysozyme, pipeline.FIRST)
    highest = pipeline.default_altlocs(lysozyme)
    assert set(first) == set(highest) and set(first.values()) == {"A"}
    assert highest[r] == "B"
    kept = pipeline.prepare(lysozyme, _replace(altloc=pipeline.FIRST, optimise=False)).structure
    chosen = lysozyme.uid[(lysozyme.residue_of == r) & (lysozyme.altloc == "A")]
    assert set(chosen.tolist()) <= set(kept.uid.tolist())
    # and Inspect lists the first as the one kept
    found = review.survey("T4L", lysozyme, _replace(altloc=pipeline.FIRST, optimise=False))
    line = next(c for c in found.choices if c.kind == review.ALTLOC and c.residue == r)
    assert line.states[line.chosen].label == "A"


def test_waters_near_ligands_by_distance(kinase) -> None:
    def kept(reach: float) -> int:
        out = pipeline.keep_waters(kinase, pipeline.NEAR_LIGANDS, reach=reach)
        return len(np.unique(out.residue_of[out.atom_kind == WATER]))

    assert kept(3.0) < kept(pipeline.WATER_REACH) == 8 < kept(8.0)


def test_poorly_placed_waters_can_be_dropped(kinase) -> None:
    water = kinase.atom_kind == WATER
    limit = float(np.median(kinase.b[water]))
    out = pipeline.keep_waters(kinase, pipeline.ALL, b_max=limit)
    left = out.atom_kind == WATER
    assert 0 < left.sum() < water.sum()
    assert out.b[left].max() <= limit
    # the chemist's own keep still wins
    hot = int(kinase.residue_of[np.flatnonzero(water & (kinase.b > limit))[0]])
    keep = pipeline.Decisions(keep_waters=frozenset({kinase.residue_label(hot)}))
    again = pipeline.keep_waters(kinase, pipeline.ALL, keep, b_max=limit)
    assert kinase.residue_label(hot) in {again.residue_label(r) for r in range(again.n_residues)}


@pytest.mark.parametrize("titrate", [sites.OFF, sites.ALL])
def test_every_his_can_be_given_one_tautomer(lysozyme, titrate) -> None:
    done = pipeline.prepare(lysozyme, _replace(his="HID", titrate=titrate))
    his = [c for c in done.choices if c.site.kind == sites.HIS]
    assert his and all(c.site.states[c.state].tautomer == "HID" for c in his)
    assert all(st.tautomer == "HID" for c in his for st in c.site.states)  # titration too
    s = done.structure
    for c in his:
        names = set(s.name[s.atoms_of(c.site.residue)].tolist())
        assert "HD1" in names and "HE2" not in names


def test_flips_can_be_turned_off(lysozyme) -> None:
    done = pipeline.prepare(lysozyme, _replace(flips=False))
    turned = [c for c in done.choices if c.site.kind in (sites.AMIDE, sites.HIS)]
    assert turned and not any(st.flipped for c in turned for st in c.site.states)


def test_the_warning_threshold_moves_the_report_and_the_review_together(lysozyme) -> None:
    def warned(threshold: float) -> tuple[int, int]:
        options = _replace(ambiguous=threshold)
        found = review.survey("T4L", lysozyme, options)
        report = sum(f.warn for f in found.prepared.report if f.kind == pipeline.SITE)
        choices = sum(c.warn for c in found.choices
                      if c.kind not in (review.ALTLOC, review.KEEP))  # fmt: skip
        return report, choices

    low, high = warned(0.0), warned(5.0)
    assert high[0] > low[0] and high[1] > low[1]
    # the default threshold gives what the preparation always gave
    default = pipeline.prepare(lysozyme, DRY)
    assert [f.warn for f in default.report] == [
        f.warn for f in pipeline.prepare(lysozyme, _replace(ambiguous=pipeline.AMBIGUOUS)).report
    ]


def test_one_location_follows_the_rule(lysozyme) -> None:
    first = repair.one_location(lysozyme, pipeline._altlocs(lysozyme, {}, pipeline.FIRST))
    assert (first.altloc == "").all()
    assert pipeline._altlocs(lysozyme, {}) is None  # the defaults: nothing to say
