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

"""Sites, the energy and the solver: flips, His states and rotors."""

from __future__ import annotations

import itertools
from dataclasses import replace

import numpy as np
import pytest

from chemsuite.protein.prepare import energy, hydrogens, sites, solve
from chemsuite.protein.structure.model import Structure
from tests.protein.conftest import load
from tests.protein.prepare.test_hydrogens import _on, _residue, _with_zinc


@pytest.fixture(scope="module")
def kinase_h() -> hydrogens.Protonated:
    return hydrogens.add_hydrogens(load("1atp.cif.gz"))


def _by_residue(found: list[sites.Site], kind: str) -> dict[int, sites.Site]:
    return {site.residue: site for site in found if site.kind == kind}


def test_every_site_has_its_states(kinase_h) -> None:
    s = kinase_h.structure
    found = sites.find_sites(kinase_h)
    his = _by_residue(found, sites.HIS)
    assert his and all(len(site.states) == 6 for site in his.values())
    one = next(iter(his.values()))
    hip = next(st for st in one.states if st.tautomer == "HIP" and not st.flipped)
    assert {"HD1", "HE2"} <= set(hip.names)
    flipped = next(st for st in one.states if st.tautomer == "HIP" and st.flipped)
    # a flip puts ND1 where CD2 was
    built = {str(s.name[i]): s.xyz()[i] for i in range(s.n_atoms) if s.residue_of[i] == one.residue}
    assert np.allclose(flipped.xyz[0, flipped.names.index("ND1")], built["CD2"], atol=1e-4)
    amides = _by_residue(found, sites.AMIDE)
    assert amides and all(len(site.states) == 2 for site in amides.values())
    rotors = {str(s.res_name[site.residue]): len(site.states)
              for site in found if site.kind == sites.ROTOR}  # fmt: skip
    assert rotors["SER"] == 24 and rotors["LYS"] == 8 and rotors["TYR"] == 2


def test_amide_hydrogens_follow_the_flip(kinase_h) -> None:
    site = next(iter(_by_residue(sites.find_sites(kinase_h), sites.AMIDE).values()))
    for state in site.states:
        names = list(state.names)
        n = state.xyz[0, names.index(names[1])]  # the amide N
        for h in names[2:]:
            assert 0.95 < np.linalg.norm(state.xyz[0, names.index(h)] - n) < 1.05


def test_a_hydrogen_bond_scores_and_a_clash_costs() -> None:
    def batch(xyz, hydrogen, polar, acceptor, donor):
        n = len(xyz)
        return energy.Batch(
            np.array(xyz, dtype=float)[None], np.array(hydrogen)[None], np.array(polar)[None],
            np.array(acceptor)[None], np.array(acceptor)[None], np.zeros((1, n), dtype=bool),
            np.array(donor, dtype=float)[None], np.ones((1, n), dtype=bool), np.zeros((1, n)),
            np.zeros((1, n), dtype=bool))  # fmt: skip

    nan = [np.nan] * 3
    donor = batch([[0, 0, 0], [1.0, 0, 0]], [False, True], [False, True], [False, False],
                  [nan, [0, 0, 0]])  # fmt: skip
    straight = batch([[2.9, 0, 0]], [False], [False], [True], [nan])
    bent = batch([[1.0, 1.9, 0]], [False], [False], [True], [nan])
    on_top = batch([[1.8, 0, 0]], [True], [False], [False], [nan])
    terms = energy.contact(donor, straight)[0, 0]
    assert terms[energy.HBOND] == pytest.approx(-1.0)
    assert energy.contact(donor, bent)[0, 0, energy.HBOND] == pytest.approx(0.0)
    assert energy.contact(donor, on_top)[0, 0, energy.CLASH] > 0
    w = energy.Weights()
    assert w.vector()[energy.CLASH] == w.clash and len(w.vector()) == len(energy.TERMS)


@pytest.mark.parametrize("seed", range(8))
def test_the_solver_finds_the_lowest_energy(seed: int) -> None:
    rng = np.random.default_rng(seed)
    sizes = rng.integers(2, 5, size=6).tolist()
    t = len(energy.TERMS)
    single = [np.zeros((n, t)) for n in sizes]
    for table in single:
        table[:, energy.HBOND] = rng.normal(size=len(table))
    pairs = {}
    for i, j in itertools.combinations(range(6), 2):
        if rng.random() < 0.6:
            pairs[(i, j)] = np.zeros((sizes[i], sizes[j], t))
            pairs[(i, j)][..., energy.HBOND] = rng.normal(size=(sizes[i], sizes[j]))
    picked = solve.solve(energy.Tables([], single, pairs), energy.Weights())

    def total(states) -> float:
        e = sum(single[i][s, energy.HBOND] for i, s in enumerate(states))
        return e + sum(m[states[i], states[j], energy.HBOND] for (i, j), m in pairs.items())

    best = min(itertools.product(*(range(n) for n in sizes)), key=total)
    assert total(picked) == pytest.approx(total(best))


def _flip_amide(s: Structure, residue: int) -> Structure:
    span = s.atoms_of(residue)
    names = [str(n) for n in s.name[span]]
    o, n = sites.AMIDES[str(s.res_name[residue])][1:3]
    i, j = span.start + names.index(o), span.start + names.index(n)
    coords = s.coords.copy()
    coords[:, [i, j]] = coords[:, [j, i]]
    return replace(s, coords=coords)


def test_a_wrongly_flipped_amide_is_put_back() -> None:
    kinase = load("1atp.cif.gz")
    _done, first = solve.optimise(hydrogens.add_hydrogens(kinase))
    clear = [c.site.residue for c in first
             if c.site.kind == sites.AMIDE and c.label == "as built" and c.gap > 3.0]  # fmt: skip
    assert len(clear) >= 3
    flipped = kinase
    for r in clear:
        flipped = _flip_amide(flipped, r)
    _done, again = solve.optimise(hydrogens.add_hydrogens(flipped))
    chosen = {c.site.residue: c.label for c in again if c.site.kind == sites.AMIDE}
    assert all(chosen[r] == "flipped" for r in clear)


def test_a_his_on_a_metal_turns_its_bare_n_to_it() -> None:
    kinase = load("1atp.cif.gz")
    zinc = _with_zinc(kinase, ["E/HIS 87/NE2"])
    his = _residue(zinc, "E/HIS 87")
    span = zinc.atoms_of(his)
    names = [str(n) for n in zinc.name[span]]
    coords = zinc.coords.copy()
    ring = [span.start + names.index(a) for a in sites.RING]
    flipped = [span.start + names.index(a) for a in sites.FLIPPED_RING]
    coords[:, ring] = zinc.coords[:, flipped]  # the input has the ring the wrong way round
    p = hydrogens.add_hydrogens(replace(zinc, coords=coords))
    done, picked = solve.optimise(p)
    s = done.structure
    (choice,) = [c for c in picked if c.site.residue == his]
    assert choice.site.states[choice.state].flipped
    (zn,) = np.flatnonzero(s.res_name[s.residue_of] == "ZN")
    ne2 = s.atoms_of(his).start + [str(n) for n in s.name[s.atoms_of(his)]].index("NE2")
    assert np.linalg.norm(s.xyz()[ne2] - s.xyz()[zn]) < 2.5
    assert _on(s, his, "NE2") == [] and _on(s, his, "ND1") == ["HD1"]
    assert done.states[his] == "HID"


def test_the_answer_is_written_back(kinase_h) -> None:
    done, picked = solve.optimise(kinase_h)
    s, before = done.structure, kinase_h.structure
    kept = np.isin(before.uid, s.uid)
    assert kept[before.element != 1].all()  # no heavy atom is lost or renumbered
    for choice in picked:
        if choice.site.kind != sites.HIS:
            continue
        r = choice.site.residue
        state = choice.site.states[choice.state]
        on = sorted(_on(s, r, "ND1") + _on(s, r, "NE2"))
        assert on == sorted(sites.HIS_HYDROGENS[state.tautomer])
        assert done.states[r] == state.tautomer
        charge = s.charge[s.atoms_of(r)].sum()
        assert charge == (1 if state.tautomer == "HIP" else 0)
    rotor_h = {u for r in done.rotors for u in r.hydrogens}
    assert rotor_h <= set(s.uid.tolist())
    assert all(np.isfinite(c.gap) for c in picked)
