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

"""Titration: the protonation forms, their free energies, populations and pKa values."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from chemsuite.protein.prepare import energy, hydrogens, populations, sites, solve
from chemsuite.protein.prepare.energy import RT_LN10, TERMS, Tables, Weights
from tests.protein.conftest import load
from tests.protein.prepare.test_hydrogens import _on, _residue


@pytest.fixture(scope="module")
def lysozyme_h() -> hydrogens.Protonated:
    from chemsuite.protein.prepare import repair, templates

    library = templates.Library()
    s = repair.build_missing(repair.one_location(load("3htb.cif.gz")), library)
    return hydrogens.add_hydrogens(s, library)


@pytest.fixture(scope="module")
def titrating(lysozyme_h) -> dict[str, sites.Site]:
    """One site of each kind: ASP, GLU, HIS, LYS, TYR, CYS and NTERM."""
    s = lysozyme_h.structure
    out: dict[str, sites.Site] = {}
    for site in sites.find_sites(lysozyme_h, titrate=True):
        if not site.titrates:
            continue
        nterm = any(st.form == "NH3+" for st in site.states)
        out.setdefault("NTERM" if nterm else str(s.res_name[site.residue]), site)
    return out


def _forms(site: sites.Site) -> dict[str, list[sites.State]]:
    out: dict[str, list[sites.State]] = {}
    for st in site.states:
        out.setdefault(st.form, []).append(st)
    return out


def test_every_titratable_group_has_its_forms(titrating) -> None:
    assert set(titrating) == {"ASP", "GLU", "HIS", "LYS", "TYR", "CYS", "NTERM"}
    counts = {name: {form: len(states) for form, states in _forms(site).items()}
              for name, site in titrating.items()}  # fmt: skip
    assert counts["ASP"] == {"ASP": 1, "ASH": 4} and counts["GLU"] == {"GLU": 1, "GLH": 4}
    assert counts["HIS"] == {"HID": 2, "HIE": 2, "HIP": 2}
    assert counts["LYS"] == {"LYS": 8, "LYN": 24} and counts["NTERM"] == {"NH3+": 8, "NH2": 24}
    assert counts["TYR"] == {"TYR": 2, "TYM": 1} and counts["CYS"] == {"CYS": 24, "CYM": 1}
    for site in titrating.values():
        for st in site.states:  # charges on the atoms, one proton between the forms
            assert sum(st.formal) == pytest.approx(sum(st.partial))
            assert st.protons in (0, 1)
        charged = {st.protons: sum(st.formal) for st in site.states}
        assert charged[1] - charged[0] == 1


def test_a_neutral_acid_has_its_hydrogen_syn_or_anti(titrating) -> None:
    from chemsuite.protein.prepare import templates

    for st in _forms(titrating["ASP"])["ASH"]:
        h = st.xyz[0, 2]
        oxygen = st.xyz[0, st.names.index(st.parents[2])]
        other = st.xyz[0, 1 - st.names.index(st.parents[2])]
        assert np.linalg.norm(h - oxygen) == pytest.approx(0.97, abs=1e-3)
        # syn: the H on the side of the other O, 0° from it about C–O; anti: 180°
        assert np.linalg.norm(h - other) < 2.4 if not st.anti else np.linalg.norm(h - other) > 2.8
        assert templates.angle_of(other, oxygen, h) > 0  # placed, not stacked


def test_titration_off_finds_p1s_sites(lysozyme_h) -> None:
    found = sites.find_sites(lysozyme_h)
    assert not [site for site in found if site.kind == sites.ACID]
    for site in found:
        if site.kind == sites.ROTOR:  # P1's rotors move only their hydrogens
            assert all(st.elements.tolist() == [hydrogens.HYDROGEN] * len(st.names)
                       for st in site.states)  # fmt: skip


def _alone(site: sites.Site, column: int | None = None, value: float = 0.0) -> Tables:
    """A table of one site with nothing around it: its intrinsic terms only, and *value* in
    *column* for its charged states."""
    single = energy.priors(site, 7.0, titrate=True)
    if column is not None:
        charged = np.array([abs(sum(st.partial)) > 0 for st in site.states])
        single[charged, column] = value
    return Tables([site], [single], {}, 7.0, True, [energy.protons(site)])


def test_alone_in_water_each_group_has_its_model_pka(titrating) -> None:
    for name, site in titrating.items():
        found = populations.titration(_alone(site), Weights()).pka[0]
        assert found == pytest.approx(sites.MODEL_PKA[name], abs=0.03), name


def test_burial_raises_an_acids_pka_and_lowers_a_bases(titrating) -> None:
    weights = Weights(desolv_acid=RT_LN10 * 2, desolv_base=RT_LN10 * 2)
    acid = _alone(titrating["ASP"], TERMS.index("desolv_acid"), 1.0)
    base = _alone(titrating["LYS"], TERMS.index("desolv_base"), 1.0)
    assert populations.titration(acid, weights).pka[0] == pytest.approx(3.8 + 2, abs=0.05)
    assert populations.titration(base, weights).pka[0] == pytest.approx(10.5 - 2, abs=0.05)


def _pair_of_acids(titrating, repulsion: float) -> Tables:
    """The Asp and the Glu, their charged forms *repulsion* apart (kcal/mol)."""
    a, b = titrating["ASP"], titrating["GLU"]
    pair = np.zeros((len(a.states), len(b.states), len(TERMS)))
    charged_a = np.array([st.protons == 0 for st in a.states])
    charged_b = np.array([st.protons == 0 for st in b.states])
    pair[np.ix_(charged_a, charged_b)] = 0.0
    pair[np.ix_(charged_a, charged_b, [TERMS.index("screened")])] = repulsion
    single = [energy.priors(a, 7.0, titrate=True), energy.priors(b, 7.0, titrate=True)]
    return Tables([a, b], single, {(0, 1): pair}, 7.0, True,
                  [energy.protons(a), energy.protons(b)])  # fmt: skip


def _exact(table: Tables, weights: Weights, grid: np.ndarray) -> list[float]:
    """Each site's pKa from the full Boltzmann sum over both sites' states."""
    single, pairs = table.energies(weights, 0.0)
    m = pairs[(0, 1)]
    curves = [np.zeros(len(grid)), np.zeros(len(grid))]
    for k, ph in enumerate(grid):
        total = np.zeros(2)
        z = 0.0
        for i, j in itertools.product(range(len(single[0])), range(len(single[1]))):
            e = single[0][i] + single[1][j] + m[i, j]
            e += RT_LN10 * ph * (table.protons[0][i] + table.protons[1][j])
            w = np.exp(-e / populations.KT)
            z += w
            total += w * np.array([table.protons[0][i], table.protons[1][j]])
        curves[0][k], curves[1][k] = total / z
    return [populations.pka(grid, c) for c in curves]


def test_two_acids_side_by_side_raise_each_other(titrating) -> None:
    grid = populations.GRID
    alone = _pair_of_acids(titrating, 0.0)
    near = _pair_of_acids(titrating, 4.0)
    weights = Weights(screened=1.0)  # the repulsion as given
    apart = populations.titration(alone, weights).pka
    assert apart == pytest.approx([3.8, 4.5], abs=0.03)
    exact = _exact(near, weights, grid)
    assert exact[0] > 3.8 + 0.1 and exact[1] > 4.5 + 0.1
    for method in (populations.MEAN_FIELD, populations.MONTE_CARLO):
        found = populations.titration(near, weights, method=method).pka
        assert found == pytest.approx(exact, abs=0.25), method


def test_apply_writes_each_form(lysozyme_h, titrating) -> None:
    chosen = {"ASP": "ASH", "LYS": "LYN", "TYR": "TYM", "CYS": "CYM", "HIS": "HIP"}
    found = [titrating[name] for name in chosen]
    picked = [next(k for k, st in enumerate(site.states) if st.form == form)
              for site, form in zip(found, chosen.values(), strict=True)]  # fmt: skip
    done = solve.apply(lysozyme_h, found, picked)
    s = done.structure
    labels = {name: lysozyme_h.structure.residue_label(titrating[name].residue) for name in chosen}
    at = {name: _residue(s, label) for name, label in labels.items()}

    def charge(name: str, atom: str) -> int:
        span = s.atoms_of(at[name])
        return int(next(s.charge[i] for i in range(span.start, span.stop) if s.name[i] == atom))

    ash = found[0].states[picked[0]]
    assert _on(s, at["ASP"], ash.parents[2]) == [ash.names[2]]
    assert charge("ASP", "OD1") == charge("ASP", "OD2") == 0
    assert len(_on(s, at["LYS"], "NZ")) == 2 and charge("LYS", "NZ") == 0
    assert _on(s, at["TYR"], "OH") == [] and charge("TYR", "OH") == -1
    assert _on(s, at["CYS"], "SG") == [] and charge("CYS", "SG") == -1
    assert charge("HIS", "ND1") == 1
    assert {done.states[at[name]] for name in chosen} == set(chosen.values())
    assert int(s.charge.sum()) == int(lysozyme_h.structure.charge.sum()) + 1 - 1 - 1 - 1 + 1


def test_the_ph_moves_the_answer(lysozyme_h) -> None:
    def forms(ph: float) -> set[str]:
        _, picked = solve.optimise(lysozyme_h, ph, titrate=True)
        return {c.site.states[c.state].form for c in picked if c.site.titrates}

    assert {"ASH", "GLH"} <= forms(1.0)
    assert not {"ASH", "GLH", "LYN", "TYM"} & forms(7.0)
    assert {"LYN", "TYM"} <= forms(13.0)
