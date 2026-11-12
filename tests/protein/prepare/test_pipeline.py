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

"""Inspect, repair and the whole preparation."""

from __future__ import annotations

import numpy as np

from chemsuite.protein.prepare import inspect, pipeline, repair, templates
from chemsuite.protein.structure import edit
from chemsuite.protein.structure.model import Structure
from tests.protein.prepare.test_hydrogens import _on, _residue


def _kinds(found: list[inspect.Finding]) -> dict[str, list[inspect.Finding]]:
    out: dict[str, list[inspect.Finding]] = {}
    for f in found:
        out.setdefault(f.kind, []).append(f)
    return out


def test_inspection_finds_what_needs_attention(lysozyme, egfr, thrombin, kinase) -> None:
    library = templates.Library()
    found = _kinds(inspect.inspect(lysozyme, library))
    assert "10 residues" in found[inspect.ALTLOCS][0].text
    assert [f.label for f in found[inspect.MISSING]] == ["A/ASN 163"]  # its OXT
    pocket = _kinds(inspect.inspect(egfr, library))
    assert len(pocket[inspect.GAP]) > 5 and all(f.warn for f in pocket[inspect.GAP])
    assert pocket[inspect.LINK][0].label == "A/CYS 797 – A/80U 1001"
    trypsin = _kinds(inspect.inspect(thrombin, library))
    assert len(trypsin[inspect.DISULFIDE]) == 4
    assert {f.label.split(" – ")[0] for f in trypsin[inspect.LINK]} == {"H/HIS 57", "H/SER 195"}
    metals = _kinds(inspect.inspect(kinase, library))[inspect.METAL]
    assert len(metals) == 2 and "E/ASP 184" in metals[0].text


def test_one_location_is_kept(lysozyme) -> None:
    one = repair.one_location(lysozyme)
    assert (one.altloc == "").all() and one.n_atoms == int(lysozyme.active().sum())
    assert set(one.uid.tolist()) <= set(lysozyme.uid.tolist())


def _without(s: Structure, label: str, names: set[str]) -> Structure:
    r = _residue(s, label)
    span = s.atoms_of(r)
    gone = np.zeros(s.n_atoms, dtype=bool)
    gone[[i for i in range(span.start, span.stop) if s.name[i] in names]] = True
    return edit.delete(s, gone)


def test_missing_side_chains_are_built(crambin) -> None:
    arg_tip = {"CD", "NE", "CZ", "NH1", "NH2"}
    ring = {"CE1", "CZ", "OH"}
    cut = _without(_without(crambin, "A/ARG 10", arg_tip), "A/TYR 29", ring)
    library = templates.Library()
    lacking = {f.label: set(f.atoms) for f in inspect.inspect(cut, library)
               if f.kind == inspect.MISSING}  # fmt: skip
    assert lacking == {"A/ARG 10": arg_tip, "A/TYR 29": ring}
    built = repair.build_missing(cut, library)
    assert not [f for f in inspect.inspect(built, library) if f.kind == inspect.MISSING]
    x = built.xyz()
    for label, name in (("A/ARG 10", "ARG"), ("A/TYR 29", "TYR")):
        r = _residue(built, label)
        span = built.atoms_of(r)
        at = {str(built.name[i]): x[i] for i in range(span.start, span.stop)}
        t = templates.bundled()[name]
        for i, j, _ in t.bonds:  # every bond the right length, the ring closed
            a, b = t.atoms[i], t.atoms[j]
            if a in at and b in at:
                ideal = np.linalg.norm(t.xyz[i] - t.xyz[j])
                assert abs(np.linalg.norm(at[a] - at[b]) - ideal) < 0.1, (label, a, b)
    assert set(cut.uid.tolist()) <= set(built.uid.tolist())


def test_preparation_in_one_call(lysozyme) -> None:
    done = pipeline.prepare(lysozyme, pipeline.Options(waters=pipeline.NONE))
    s = done.structure
    assert not (s.res_name == "HOH").any() and (s.altloc == "").all()
    end = _residue(s, "A/ASN 163")
    assert _on(s, end, "OXT") == [] and s.charge[s.atoms_of(end)].sum() == -1
    assert (s.element == 1).any() and isinstance(done.net_charge, int)
    kinds = {f.kind for f in done.report}
    assert {inspect.ALTLOCS, inspect.MISSING} <= kinds
    assert done.choices and all(c.gap >= 0 for c in done.choices)
    # every heavy atom shown in the input is still there, with its uid
    shown = lysozyme.uid[lysozyme.active() & (lysozyme.res_name[lysozyme.residue_of] != "HOH")]
    assert set(shown.tolist()) <= set(s.uid.tolist())


def test_waters_near_ligands(kinase) -> None:
    near = pipeline.keep_waters(kinase, pipeline.NEAR_LIGANDS)
    kept = (near.res_name == "HOH").sum()
    assert kept == 8 and kept < (kinase.res_name == "HOH").sum()  # 8 within 5 Å of ATP
    assert pipeline.keep_waters(kinase, pipeline.ALL) is kinase
