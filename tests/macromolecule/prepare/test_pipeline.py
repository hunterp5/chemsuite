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

from dataclasses import replace

import numpy as np
import pytest

from chemsuite.macromolecule.prepare import inspect, pipeline, repair, sites, templates
from chemsuite.macromolecule.structure import edit
from chemsuite.macromolecule.structure.model import Structure
from tests.macromolecule.prepare.test_hydrogens import _on, _residue


def _kinds(found: list[inspect.Finding]) -> dict[str, list[inspect.Finding]]:
    out: dict[str, list[inspect.Finding]] = {}
    for f in found:
        out.setdefault(f.kind, []).append(f)
    return out


def test_inspection_finds_what_needs_attention(lysozyme, egfr, thrombin, kinase) -> None:
    library = templates.Library()
    found = _kinds(inspect.inspect(lysozyme, library))
    assert "10 residues" in found[inspect.ALTLOCS][0].text
    assert inspect.MISSING not in found  # Asn 163 is not the end of the sequence: no OXT
    (cut,) = found[inspect.TRUNCATED]
    assert (cut.label, cut.text) == ("A/ASN 163", "1 C-terminal residues not modelled")
    pocket = _kinds(inspect.inspect(egfr, library))
    assert len(pocket[inspect.GAP]) > 5 and all(f.warn for f in pocket[inspect.GAP])
    assert pocket[inspect.LINK][0].label == "A/CYS 797 – A/80U 1001"
    trypsin = _kinds(inspect.inspect(thrombin, library))
    assert len(trypsin[inspect.DISULFIDE]) == 4
    assert {f.label.split(" – ")[0] for f in trypsin[inspect.LINK]} == {"H/HIS 57", "H/SER 195"}
    metals = _kinds(inspect.inspect(kinase, library))[inspect.METAL]
    assert len(metals) == 2 and "E/ASP 184" in metals[0].text


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u, w = a - b, c - b
    return float(np.degrees(np.arccos(np.dot(u, w) / np.linalg.norm(u) / np.linalg.norm(w))))


def test_a_residue_missing_alone_is_built_from_its_neighbours(egfr) -> None:
    egfr = repair.one_location(egfr)
    filled, built = repair.fill_single_gaps(egfr)
    assert built == ["A/ILE 706", "A/ILE 715", "A/ASN 826", "A/ARG 832", "A/LYS 860"]
    assert filled.n_residues == egfr.n_residues + len(built)
    assert set(egfr.uid.tolist()) <= set(filled.uid.tolist())
    x = filled.xyz()
    pairs = {tuple(sorted(p)) for p in filled.bonds.tolist()}
    for label in built:
        r = _residue(filled, label)
        assert filled.seqres[r] == filled.seqres[r - 1] + 1 == filled.seqres[r + 1] - 1
        assert (
            str(filled.sequence_of(int(filled.entity[r]))[filled.seqres[r]]) == filled.res_name[r]
        )
        at = {str(filled.name[i]): i for i in range(*filled.atoms_of(r).indices(filled.n_atoms))}
        assert set(at) == {"N", "CA", "C", "O"}  # the side chain comes from build_missing
        c0 = [i for i in range(*filled.atoms_of(r - 1).indices(filled.n_atoms))
              if filled.name[i] == "C"][0]  # fmt: skip
        n1 = filled.atoms_of(r + 1).start
        assert filled.name[n1] == "N"
        chain = [(c0, at["N"]), (at["N"], at["CA"]), (at["CA"], at["C"]), (at["C"], n1)]
        for (i, j), length in zip(chain, (1.335, 1.458, 1.525, 1.335), strict=True):
            assert tuple(sorted((i, j))) in pairs
            assert np.linalg.norm(x[i] - x[j]) == pytest.approx(length, abs=0.15), label
        assert _angle(x[at["N"]], x[at["CA"]], x[at["C"]]) == pytest.approx(111, abs=4)
    # in the preparation: the gaps it closed take no caps, and every report line still
    # selects its own residue although residues were put in
    done = pipeline.prepare(egfr, pipeline.Options(waters=pipeline.NONE))
    s = done.structure
    for label in built:
        r = _residue(s, label)
        assert s.res_name[r - 1] != "NME" and s.res_name[r + 1] != "ACE"
        assert "CB" in s.name[s.atoms_of(r)]  # its side chain built, turned clear of the rest
        heavy = np.flatnonzero(s.element != 1)
        far = heavy[np.abs(s.residue_of[heavy] - r) > 1]
        for i in range(*s.atoms_of(r).indices(s.n_atoms)):
            if s.element[i] != 1 and s.name[i] not in ("N", "CA", "C", "O"):
                assert np.linalg.norm(s.xyz()[far] - s.xyz()[i], axis=1).min() > 2.4, label
    rows = [f for f in done.report if f.kind == inspect.FILLED]
    assert [f.label for f in rows] == built and all(f.warn for f in rows)
    lines = [f for f in done.report if f.residue >= 0]
    assert all(f.label.startswith(s.residue_label(f.residue)) for f in lines)


def test_gaps_and_cut_ends_take_neutral_caps(egfr) -> None:
    library = templates.Library()
    egfr = repair.one_location(egfr)  # as the preparation caps them
    egfr, _built = repair.fill_single_gaps(egfr)  # and only once the lone residues are built
    ends = inspect.ends(egfr)
    gaps = [f for f in inspect.inspect(egfr, library) if f.kind == inspect.GAP]
    assert len(ends.open_c) >= len(gaps) and len(ends.open_n) >= len(gaps)
    capped, caps = repair.cap_ends(egfr)
    assert len(caps) == len(ends.open_n) + len(ends.open_c)
    assert capped.n_residues == egfr.n_residues + len(caps)
    assert set(egfr.uid.tolist()) <= set(capped.uid.tolist())
    assert not inspect.ends(capped).open_n and not inspect.ends(capped).open_c  # nothing left
    x = capped.xyz()
    pairs = {tuple(sorted(p)) for p in capped.bonds.tolist()}

    def atom(r: int, name: str) -> int:
        span = capped.atoms_of(r)
        (i,) = [i for i in range(span.start, span.stop) if capped.name[i] == name]
        return i

    for r in np.flatnonzero(capped.res_name == "NME").tolist():
        n, c = atom(r - 1, "C"), atom(r, "N")
        assert tuple(sorted((n, c))) in pairs
        assert np.linalg.norm(x[n] - x[c]) == pytest.approx(repair.PEPTIDE, abs=1e-3)
    for r in np.flatnonzero(capped.res_name == "ACE").tolist():
        c, n = atom(r, "C"), atom(r + 1, "N")
        assert tuple(sorted((n, c))) in pairs
        assert np.linalg.norm(x[n] - x[c]) == pytest.approx(repair.PEPTIDE, abs=1e-3)
        others = np.ones(capped.n_atoms, dtype=bool)
        others[capped.atoms_of(r)] = False
        others[[n, atom(r + 1, "CA")]] = False
        assert np.linalg.norm(x[others] - x[c], axis=1).min() > 2.0  # placed clear of the rest
    # with hydrogens: the capped ends are amides, not charged termini
    done = pipeline.prepare(egfr, pipeline.Options(waters=pipeline.NONE))
    s = done.structure
    for r in np.flatnonzero(s.res_name == "ACE").tolist():
        span = s.atoms_of(r + 1)
        n = [i for i in range(span.start, span.stop) if s.name[i] == "N"][0]
        assert s.charge[n] == 0 and len(_on(s, r + 1, "N")) == 1
        assert len(_on(s, r, "CH3")) == 3 and _on(s, r, "C") == []
    for r in np.flatnonzero(s.res_name == "NME").tolist():
        assert len(_on(s, r, "N")) == 1 and len(_on(s, r, "C")) == 3
        assert "OXT" not in s.name[s.atoms_of(r - 1)]
    assert sum(f.kind == inspect.CAP for f in done.report) == len(caps)


def test_a_chain_end_too_incomplete_to_cap_is_removed_and_the_one_before_capped(egfr) -> None:
    """4AGC's last residue is deposited as its N alone: no C to cap, so Minimize refused it
    ("The chain ends without an OXT or a cap"). It is removed and the residue before capped."""
    from chemsuite.simulate import build

    egfr = repair.one_location(egfr)
    tail = max(inspect.ends(egfr).open_c)
    label, before = egfr.residue_label(tail), egfr.residue_label(tail - 1)
    stub = _without(egfr, label, {"CA", "C", "O", "CB", "SG", "CG", "CD", "OG", "OG1", "CG1",
                                  "CG2", "CD1", "CD2", "CE", "CE1", "CE2", "CZ", "NZ", "OD1",
                                  "OD2", "OE1", "OE2", "ND1", "ND2", "NE", "NE1", "NE2", "NH1",
                                  "NH2", "OH", "SD", "CE3", "CZ2", "CZ3", "CH2"})  # fmt: skip
    assert stub.name[stub.atoms_of(_residue(stub, label))].tolist() == ["N"]
    done = pipeline.prepare(stub, pipeline.Options(waters=pipeline.NONE, optimise=False))
    s = done.structure
    assert label not in {s.residue_label(r) for r in range(s.n_residues)}
    removed = [f for f in done.report if f.kind == inspect.REMOVED]
    assert [f.label for f in removed] == [label] and removed[0].warn
    assert not any(f.kind == inspect.MISSING and f.label == label for f in done.report)
    assert any(f.kind == inspect.CAP and f.label == before for f in done.report)
    assert s.res_name[_residue(s, before) + 1] == "NME"
    blocking = build.problems(s, done.states)
    assert not any(p.label in (label, before) for p in blocking)


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


def test_a_metal_does_not_push_a_rebuilt_side_chain_away() -> None:
    """The turning that keeps rebuilt side chains clear of their neighbours leaves ions out: a
    Cys or His that binds a metal is meant to sit close to it. A zinc finger's Cys, cut back to
    CB with its zinc moved onto where SG belongs, is rebuilt as it would be with no zinc."""
    from tests.macromolecule.conftest import load

    s = repair.one_location(load("1zaa_pocket.cif.gz"))
    s = s.subset(s.element != 1)
    library = templates.Library()
    cys = next(r for r in range(s.n_residues) if s.res_name[r] == "CYS")
    label = s.residue_label(cys)
    cut = _without(s, label, {"SG"})
    free = repair.build_missing(cut, library)
    span = free.atoms_of(_residue(free, label))
    (sg,) = [i for i in range(span.start, span.stop) if free.name[i] == "SG"]
    zinc = np.flatnonzero(cut.element == 30)[0]
    coords = cut.coords.copy()
    coords[:, zinc] = free.coords[:, sg] + [0.0, 0.0, 1.0]  # 1 Å off: well inside CONTACT
    beside = repair.build_missing(replace(cut, coords=coords), library)
    span = beside.atoms_of(_residue(beside, label))
    (moved,) = [i for i in range(span.start, span.stop) if beside.name[i] == "SG"]
    assert np.allclose(beside.xyz()[moved], free.xyz()[sg], atol=1e-3)


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
    assert "OXT" not in s.name[s.atoms_of(end)] and s.charge[s.atoms_of(end)].sum() == 0
    assert s.residue_label(end + 1) == "A/NME 164"  # the sequence goes on: capped, not charged
    assert (s.element == 1).any() and isinstance(done.net_charge, int)
    kinds = {f.kind for f in done.report}
    assert {inspect.ALTLOCS, inspect.TRUNCATED, inspect.CAP} <= kinds
    (cap,) = [f for f in done.report if f.kind == inspect.CAP]
    assert (cap.label, cap.text, cap.residue) == ("A/ASN 163", "Capped with NME", end)
    assert done.choices and all(c.gap >= 0 for c in done.choices)
    # every heavy atom shown in the input is still there, with its uid
    shown = lysozyme.uid[lysozyme.active() & (lysozyme.res_name[lysozyme.residue_of] != "HOH")]
    assert set(shown.tolist()) <= set(s.uid.tolist())
    # and the heavy atoms it built are named, the cap's among them, no deposited one
    new = set(s.uid[s.element != 1].tolist()) - set(lysozyme.uid.tolist())
    assert done.built == new and not done.built & set(lysozyme.uid.tolist())
    assert set(s.uid[s.atoms_of(end + 1)][s.element[s.atoms_of(end + 1)] != 1].tolist()) <= new


def test_waters_near_ligands(kinase) -> None:
    near = pipeline.keep_waters(kinase, pipeline.NEAR_LIGANDS)
    kept = (near.res_name == "HOH").sum()
    assert kept == 8 and kept < (kinase.res_name == "HOH").sum()  # 8 within 5 Å of ATP
    assert pipeline.keep_waters(kinase, pipeline.ALL) is kinase


def test_titration_in_one_call(lysozyme) -> None:
    done = pipeline.prepare(lysozyme, pipeline.Options(waters=pipeline.NONE, titrate=sites.ALL))
    by_residue = {c.site.residue: c for c in done.choices
                  if c.site.titrates and c.site.kind != sites.LIGAND}  # fmt: skip
    ligand = {c.site.residue: c for c in done.choices if c.site.kind == sites.LIGAND}
    assert by_residue and set(done.pka) == set(by_residue)
    assert set(done.population) == set(by_residue) | set(ligand)
    for r, c in ligand.items():  # a ligand's share is that of the form it takes
        assert 1 / len({st.form for st in c.site.states}) <= done.population[r] <= 1.0
    for r, c in by_residue.items():
        share = done.population[r]
        assert 0.0 <= share <= 1.0 and 0.0 <= done.pka[r] <= 14.0
        n = [st.protons for st in c.site.states]
        # the structure carries the form most of the group is in
        assert n[c.state] == (max(n) if share >= 0.5 else min(n))
    mixed = [f for f in done.report if f.kind == pipeline.SITE and "protonated" in f.text]
    assert all(f.warn == (pipeline.MIXED < done.population[f.residue] < 1 - pipeline.MIXED)
               for f in mixed)  # fmt: skip


def test_titrating_the_pocket_titrates_only_the_residues_near_a_ligand(lysozyme) -> None:
    options = pipeline.Options(waters=pipeline.NONE, titrate=sites.POCKET)
    done = pipeline.prepare(lysozyme, options)
    near = pipeline._pocket(done.structure)
    assert 0 < int(near.sum()) < done.structure.n_residues
    titrated = {c.site.residue for c in done.choices if c.site.titrating}
    assert done.pka and set(done.pka) <= titrated
    # every titrated protein group is in the pocket; the rest keep their textbook form
    protein_groups = {c.site.residue for c in done.choices
                      if c.site.titrating and c.site.kind != sites.LIGAND}  # fmt: skip
    assert protein_groups and all(near[r] for r in protein_groups)
    everywhere = pipeline.prepare(lysozyme, replace(options, titrate=sites.ALL))
    assert len(everywhere.pka) > len(done.pka)


def test_only_rotors_with_something_to_decide_are_reported(lysozyme) -> None:
    done = pipeline.prepare(lysozyme, pipeline.Options(waters=pipeline.NONE))
    rotors = [c for c in done.choices if c.site.kind == sites.ROTOR]
    assert rotors
    rows = [f for f in done.report if f.kind == pipeline.SITE and f.text.startswith("Rotor:")]
    assert rows and len(rows) < len(rotors)
    for f in rows:  # each one chose between orientations that differ, and it bonds either way
        here = [c for c in rotors if c.site.residue == f.residue]
        assert any(c.energy <= pipeline.BONDED
                   and pipeline.DEGENERATE < c.gap < pipeline.AMBIGUOUS for c in here)  # fmt: skip
        assert f.warn
    free = sum(1 for c in rotors if c.energy > pipeline.BONDED)
    (line,) = [f for f in done.report if f.label == "Rotors"]
    assert line.text.endswith(f": {free}") and not line.warn  # free to spin: nothing to check
