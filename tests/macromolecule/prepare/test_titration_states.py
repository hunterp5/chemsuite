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

"""Titration decides states (docs/prepare/p5-review.md): a group leaves its form in water only
past its threshold, and a titrating His keeps P1's hard clash on its neutral states."""

from __future__ import annotations

import numpy as np

from chemsuite.macromolecule.prepare import energy, sites
from chemsuite.macromolecule.prepare.energy import TERMS
from tests.macromolecule.prepare.test_titration import lysozyme_h  # noqa: F401


def test_a_group_leaves_its_form_in_water_only_past_its_threshold() -> None:
    from chemsuite.macromolecule.prepare.pipeline import protonated

    strict = {"ASP": 0.9, "LYS": 0.9}
    # at pH 7 an Asp is charged and a Lys protonated in water: they switch only when the other
    # form reaches the threshold
    assert not protonated("ASP", 0.6, 7.0, strict) and protonated("ASP", 0.95, 7.0, strict)
    assert protonated("LYS", 0.4, 7.0, strict) and not protonated("LYS", 0.05, 7.0, strict)
    # at pH 3 the protonated Asp is the usual form, so the same rule holds it there
    assert protonated("ASP", 0.4, 3.0, strict) and not protonated("ASP", 0.05, 3.0, strict)
    # without a threshold, a plain majority
    assert protonated("GLU", 0.6, 7.0, strict) and not protonated("GLU", 0.4, 7.0, strict)


def test_every_titrating_protein_site_names_its_group(lysozyme_h) -> None:  # noqa: F811
    s = lysozyme_h.structure
    for site in sites.find_sites(lysozyme_h, titrate=sites.ALL):
        name = sites.group(site)
        if site.titrates:
            nterm = any(st.form == "NH3+" for st in site.states)
            assert name == ("NTERM" if nterm else str(s.res_name[site.residue]))
        else:
            assert name == ""


def test_a_titrating_his_keeps_hard_clash_on_its_neutral_states(lysozyme_h) -> None:  # noqa: F811
    s = lysozyme_h.structure
    found = sites.find_sites(lysozyme_h, titrate=sites.ALL)
    table = energy.tables(s, found, 7.0)
    seen_clash = False
    for i, site in enumerate(found):
        if site.kind != sites.HIS:
            continue
        assert site.titrating
        hip = np.array([st.tautomer == "HIP" for st in site.states])
        assert energy.soft_states(site).tolist() == hip.tolist()
        # its own terms, and its pairs with sites that do not titrate (a pair with a titrating
        # group is softened wherever that group's proton could be the one in the way)
        rows = [table.single[i]] + [
            np.moveaxis(t, 0 if a == i else 1, 0)
            for (a, b), t in table.pairs.items()
            if i in (a, b) and not found[b if a == i else a].titrating
        ]  # fmt: skip
        for t in rows:
            t = t.reshape(len(site.states), -1, len(TERMS))
            assert not t[~hip][..., energy.BUMP].any()  # neutral: hard clash only
            assert not t[hip][..., energy.CLASH].any()  # HIP: softened
            seen_clash |= bool(t[~hip][..., energy.CLASH].any())
    assert seen_clash
