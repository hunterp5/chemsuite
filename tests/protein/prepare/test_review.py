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
"""A preparation reviewed and overruled: the survey, and the chemist's decisions applied."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.protein.prepare import inspect, pipeline, review, sites
from tests.protein.prepare.test_hydrogens import _residue
from tests.protein.prepare.test_pipeline import _without

NO_WATERS = pipeline.Options(waters=pipeline.NONE)


@pytest.fixture
def cut(crambin):
    """Crambin missing two side chains, so there is something to repair or leave alone."""
    arg_tip = {"CD", "NE", "CZ", "NH1", "NH2"}
    ring = {"CE1", "CZ", "OH"}
    return _without(_without(crambin, "A/ARG 10", arg_tip), "A/TYR 29", ring)


def _same(a, b) -> bool:
    return (
        a.n_atoms == b.n_atoms
        and np.array_equal(a.uid, b.uid)
        and np.allclose(a.coords, b.coords)
        and np.array_equal(a.charge, b.charge)
    )


def test_with_no_decisions_the_survey_proposes_exactly_what_prepare_does(lysozyme) -> None:
    """One code path: reviewing changes nothing unless the chemist overrules something."""
    plain = pipeline.prepare(lysozyme, NO_WATERS)
    found = review.survey("LYZ", lysozyme, NO_WATERS)
    assert _same(found.prepared.structure, plain.structure)
    assert found.net_charge == plain.net_charge
    again = pipeline.prepare(lysozyme, NO_WATERS, decisions=review.decisions_from(found))
    assert _same(again.structure, plain.structure)


def test_the_survey_lists_issues_and_every_state_of_every_choice(lysozyme) -> None:
    found = review.survey("LYZ", lysozyme, NO_WATERS)
    kinds = {issue.kind for issue in found.issues}
    assert {inspect.ALTLOCS, inspect.MISSING} <= kinds
    assert all(issue.fixable == (issue.kind in review.FIXABLE) for issue in found.issues
               if issue.residue >= 0)  # fmt: skip
    assert len({issue.id for issue in found.issues}) == len(found.issues)  # ids are unique
    assert found.choices and len({c.id for c in found.choices}) == len(found.choices)
    for choice in found.choices:
        assert len(choice.states) >= 1 and choice.state is not None and choice.state.chosen
        assert sum(state.chosen for state in choice.states) == 1
        # the chosen state is the cheapest one, and the gap is to the next
        costs = [state.cost for state in choice.states]
        assert min(costs) == pytest.approx(0.0, abs=1e-6)
        if len(costs) > 1:
            assert choice.state.cost == pytest.approx(0.0, abs=1e-6)
            others = sorted(c for k, c in enumerate(costs) if k != choice.chosen)
            assert choice.gap == pytest.approx(others[0], abs=1e-4)
    amides = [c for c in found.choices if c.kind == sites.AMIDE]
    assert amides and all(len(c.states) == 2 for c in amides)  # as found, flipped


def test_overruling_a_flip_holds_it_and_solves_the_rest_around_it(lysozyme) -> None:
    found = review.survey("LYZ", lysozyme, NO_WATERS)
    amide = next(c for c in found.choices if c.kind == sites.AMIDE)
    other = next(s.label for k, s in enumerate(amide.states) if k != amide.chosen)
    decisions = review.decisions_from(found, overrides={amide.id: other})
    assert decisions.states == {amide.id: other}
    again = review.survey("LYZ", lysozyme, NO_WATERS, decisions=decisions)
    held = again.choice(amide.id)
    assert held is not None and held.state.label == other
    # and the structure carries it: the amide's O and N have swapped places
    before = found.prepared.structure
    after = again.prepared.structure
    r = amide.residue
    o, n = ("OD1", "ND2") if before.res_name[r] == "ASN" else ("OE1", "NE2")

    def at(s, name):
        span = s.atoms_of(r)
        return s.xyz()[next(i for i in range(span.start, span.stop) if s.name[i] == name)]

    assert np.allclose(at(after, o), at(before, n), atol=0.05)
    assert np.allclose(at(after, n), at(before, o), atol=0.05)


def test_a_pinned_his_takes_the_tautomer_it_was_given(lysozyme) -> None:
    found = review.survey("LYZ", lysozyme, NO_WATERS)
    his = next((c for c in found.choices if c.kind == sites.HIS), None)
    if his is None:
        pytest.skip("no His in this structure")
    other = next(s for k, s in enumerate(his.states) if k != his.chosen)
    decisions = review.decisions_from(found, overrides={his.id: other.label})
    again = review.survey("LYZ", lysozyme, NO_WATERS, decisions=decisions)
    assert again.choice(his.id).state.label == other.label
    assert again.prepared.states.get(his.residue, "HIE") == (other.form or "HIE")


def test_an_override_naming_no_such_state_is_ignored(lysozyme) -> None:
    found = review.survey("LYZ", lysozyme, NO_WATERS)
    amide = next(c for c in found.choices if c.kind == sites.AMIDE)
    decisions = review.decisions_from(found, overrides={amide.id: "sideways", "no:such": "x"})
    assert decisions.states == {}


def test_an_unticked_repair_leaves_that_residue_as_the_file_has_it(cut) -> None:
    found = review.survey("CRN", cut, NO_WATERS)
    missing = {issue.label: issue for issue in found.issues if issue.kind == inspect.MISSING}
    assert set(missing) == {"A/ARG 10", "A/TYR 29"} and all(i.fixable for i in missing.values())
    decisions = review.decisions_from(found, skip_issues={missing["A/ARG 10"].id: True})
    assert decisions.repair == frozenset({missing["A/TYR 29"].residue})
    done = pipeline.prepare(cut, NO_WATERS, decisions=decisions)
    s = done.structure

    def names(label):
        span = s.atoms_of(_residue(s, label))
        return {str(s.name[i]) for i in range(span.start, span.stop)}

    assert "CZ" not in names("A/ARG 10")  # left alone
    assert {"CE1", "CZ", "OH"} <= names("A/TYR 29")  # built


def test_waters_can_be_kept_or_dropped_one_by_one(kinase) -> None:
    water = np.flatnonzero(kinase.res_name == "HOH")
    first, second = int(water[0]), int(water[1])
    keep_one = pipeline.Decisions(keep_waters=frozenset({first}))
    kept = pipeline.keep_waters(kinase, pipeline.NONE, keep_one)
    assert (kept.res_name == "HOH").sum() == 1
    drop_one = pipeline.Decisions(drop_waters=frozenset({second}))
    dropped = pipeline.keep_waters(kinase, pipeline.ALL, drop_one)
    assert (dropped.res_name == "HOH").sum() == len(water) - 1


def test_a_review_warns_exactly_where_the_report_does(kinase) -> None:
    """One ⚠ rule: a line the report does not flag is not flagged for review either, so
    waters and free rotors (near ties that mean nothing) stay out of "⚠ Only"."""
    found = review.survey("PKA", kinase, pipeline.Options(waters=pipeline.NEAR_LIGANDS))
    waters = [c for c in found.choices if c.kind == sites.WATER]
    assert waters and not any(c.warn for c in waters)
    flagged = {row.residue for row in found.prepared.report if row.warn}
    warned = {c.residue for c in found.choices if c.warn}
    assert warned and warned <= flagged
