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

"""Missing loops built by cyclic coordinate descent (docs/LOOPS_PLAN.md L1)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial import cKDTree

from chemsuite.macromolecule.prepare import loops, repair, templates
from chemsuite.macromolecule.structure import edit
from chemsuite.macromolecule.structure.model import Structure

BACKBONE = ("N", "CA", "C", "O")


def _cut(s: Structure, chain: str, first: int, last: int) -> Structure:
    """*s* (one location) without residues *first* to *last* of *chain*: a gap in the
    deposited sequence."""
    s = repair.one_location(s)
    gone = np.zeros(s.n_atoms, dtype=bool)
    for r in range(s.n_residues):
        if s.chain_label(int(s.chain_of[r])) == chain and first <= s.res_seq[r] <= last:
            gone[s.atoms_of(r)] = True
    return edit.delete(s, gone)


def _atom(s: Structure, label: str, name: str, model: int = 0) -> np.ndarray:
    r = next(r for r in range(s.n_residues) if s.residue_label(r) == label)
    span = s.atoms_of(r)
    (i,) = [i for i in range(span.start, span.stop) if s.name[i] == name]
    return s.coords[model, i].astype(np.float64)


@pytest.fixture
def library() -> templates.Library:
    return templates.Library()


def test_a_gap_is_found_with_its_sequence(lysozyme) -> None:
    s = _cut(lysozyme, "A", 53, 58)
    (gap,) = loops.gaps(s)
    assert s.res_seq[gap.before] == 52 and s.res_seq[gap.after] == 59
    original = repair.one_location(lysozyme)
    names = [str(original.res_name[r]) for r in range(original.n_residues)
             if 53 <= original.res_seq[r] <= 58 and original.chain_of[r] == 0]  # fmt: skip
    assert list(gap.names) == names
    assert loops.gaps(s, longest=5) == []


def test_a_loop_closes_onto_both_anchors_without_clashes(lysozyme, library) -> None:
    original = repair.one_location(lysozyme)
    s = _cut(lysozyme, "A", 53, 58)
    (gap,) = loops.gaps(s)
    before, after = s.residue_label(gap.before), s.residue_label(gap.after)
    out, built = loops.fill_gaps(s, library)
    labels = [original.residue_label(r) for r in range(original.n_residues)]
    assert built == labels[labels.index(before) + 1 : labels.index(after)]
    chain = [before, *built, after]
    for x, y in zip(chain[:-1], chain[1:], strict=True):
        bond = np.linalg.norm(_atom(out, x, "C") - _atom(out, y, "N"))
        assert abs(bond - loops.C_N) < 0.15, (x, y, bond)
    # every heavy atom of each new residue, as the template has it
    for label in built:
        r = next(r for r in range(out.n_residues) if out.residue_label(r) == label)
        t = library.get(str(out.res_name[r]))
        names = {str(out.name[i]) for i in range(*out.atoms_of(r).indices(out.n_atoms))}
        expected = {t.atoms[i] for i in range(len(t.atoms)) if t.elements[i] != 1} - {"OXT"}
        assert expected <= names, (label, expected - names)
    # no new heavy atom within 2.2 Å of one two or more residues away
    heavy = out.element != 1
    x = out.coords[0].astype(np.float64)
    new = np.isin(out.residue_of, [r for r in range(out.n_residues)
                                    if out.residue_label(r) in built])  # fmt: skip
    tree = cKDTree(x[heavy])
    owner = out.residue_of[heavy]
    for i in np.flatnonzero(new & heavy):
        for j in tree.query_ball_point(x[i], 2.2):
            assert abs(int(owner[j]) - int(out.residue_of[i])) < 2, out.atom_label(i)
    # the rest stays where it was (a's O aside, rebuilt in the new peptide plane)
    for label in (before, after, "A/THR 21"):
        for name in ("N", "CA", "C"):
            assert np.allclose(_atom(out, label, name), _atom(original, label, name), atol=1e-3)
    # near the deposit, generously: it is one plausible loop, not the crystal's
    rmsd = np.sqrt(np.mean([np.sum((_atom(out, b, n) - _atom(original, b, n)) ** 2)
                            for b in built for n in BACKBONE]))  # fmt: skip
    assert rmsd < 4.0


def test_new_residues_have_allowed_backbone_angles(lysozyme, library) -> None:
    out, built = loops.fill_gaps(_cut(lysozyme, "A", 53, 58), library)
    labels = [out.residue_label(r) for r in range(out.n_residues)]
    grids = loops.grids()
    for label in built:
        k = labels.index(label)
        prev, nxt = labels[k - 1], labels[k + 1]
        phi = loops._dihedral(_atom(out, prev, "C"), _atom(out, label, "N"),
                              _atom(out, label, "CA"), _atom(out, label, "C"))  # fmt: skip
        psi = loops._dihedral(_atom(out, label, "N"), _atom(out, label, "CA"),
                              _atom(out, label, "C"), _atom(out, nxt, "N"))  # fmt: skip
        name, following = label.split("/")[1].split()[0], nxt.split("/")[1].split()[0]
        grid = grids[loops.rama_class(name, following)]
        assert grid.at(phi, psi) >= grid.allowed, (label, np.degrees([phi, psi]))


def test_the_same_loop_every_time(lysozyme, library) -> None:
    s = _cut(lysozyme, "A", 53, 58)
    first, _ = loops.fill_gaps(s, library)
    again, _ = loops.fill_gaps(s, library)
    assert np.array_equal(first.coords, again.coords)


def test_a_gap_too_wide_for_its_residues_is_left_to_the_caps(lysozyme, library) -> None:
    s = _cut(lysozyme, "A", 53, 60)  # eight residues out of the structure ...
    seqres = s.seqres.copy()
    after = next(r for r in range(s.n_residues) if s.residue_label(r) == "A/ASP 61")
    chain_end = s.residues_of(int(s.chain_of[after])).stop
    seqres[after:chain_end] -= 6  # ... but the sequence says two are missing
    narrow = replace(s, seqres=seqres)
    (gap,) = loops.gaps(narrow)
    assert len(gap.names) == 2
    out, built = loops.fill_gaps(narrow, library)
    assert built == [] and out is narrow


def test_tracked_torsions_match_the_coordinates() -> None:
    """The CCD keeps each torsion's value as it turns it; they must stay the chain's own."""
    rng = np.random.default_rng(1)
    start = np.array([[0.0, 1.4, 0.0], [0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    angles = rng.uniform(-np.pi, np.pi, (5, 8))
    chain = loops._build(start, angles)
    target = chain[0, -3:] + 0.5
    classes = ["general"] * 5
    chain, angles, _rmsd = loops._ccd(chain, angles, target, classes, None, None, loops.grids())
    rebuilt = loops._build(start, angles)
    assert np.allclose(rebuilt, chain, atol=1e-6)


def test_a_peptide_can_be_built_cis() -> None:
    """A Pro's peptide may be cis (4RXV's Pro 199, 1A4I's Pro 102 in the loop benchmark)."""
    start = np.array([[0.0, 1.4, 0.0], [0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    angles = np.radians([[150.0, -70.0, 150.0, -70.0, 150.0, -70.0]])  # ψa, φ1, ψ1, φ2, ψ2, φb
    omega = np.array([[np.pi, 0.0, np.pi]])
    chain = loops._build(start, angles, omega=omega)
    for r, expected in ((1, 180.0), (2, 0.0), (3, 180.0)):
        ca_prev, c_prev, n, ca = (chain[0, i] for i in (3 * r - 2, 3 * r - 1, 3 * r, 3 * r + 1))
        w = np.degrees(loops._dihedral(ca_prev, c_prev, n, ca))
        assert abs(np.angle(np.exp(1j * np.radians(w - expected)))) < 1e-6


def test_each_model_closes_on_its_own_anchors(nmr, library) -> None:
    s = _cut(nmr, "A", 8, 10)
    s = s.subset(s.element != 1)
    out, built = loops.fill_gaps(s, library)
    assert len(built) == 3
    for model in range(out.n_models):
        c = _atom(out, built[-1], "C", model)
        n = _atom(out, "A/LYS 11", "N", model)
        assert abs(np.linalg.norm(c - n) - loops.C_N) < 0.15
    assert not np.allclose(_atom(out, built[0], "CA", 0), _atom(out, built[0], "CA", 1))


def test_close_gaps_names_each_gap_and_leaves_skipped_ones(lysozyme, library) -> None:
    s = _cut(lysozyme, "A", 55, 58)
    out, closed = loops.close_gaps(s, library, 8)
    (one,) = closed
    assert one.gap == "A/CYS 54 – A/THR 59"
    assert len(one.residues) == 4 and out.n_residues == s.n_residues + 4
    left, none = loops.close_gaps(s, library, 8, skip={one.gap})
    assert none == [] and left is s
    _same, too_short = loops.close_gaps(s, library, 3)  # 4 missing: longer than asked
    assert too_short == []


def test_a_remembered_loop_is_the_loop_that_would_be_built(lysozyme, library) -> None:
    s = _cut(lysozyme, "A", 55, 58)
    loops.forget()
    built, closed = loops.close_gaps(s, library, 8)
    assert len(loops._memo) == 1
    remembered, again_closed = loops.close_gaps(s, library, 8)  # from the memo
    assert again_closed == closed
    for a, b in ((remembered.coords, built.coords), (remembered.uid, built.uid),
                 (remembered.name, built.name), (remembered.charge, built.charge)):  # fmt: skip
        assert np.array_equal(a, b)
    moved = replace(s, coords=s.coords + np.float32(0.01))  # another structure: another key
    assert loops.fingerprint(moved) != loops.fingerprint(s)


def test_a_remembered_loop_takes_the_rest_from_the_structure_as_it_is_now(lysozyme, library):
    """A charge set by hand moves no atom, so the memo is hit, but the result keeps the new
    charge: only the backbone comes from the memo."""
    s = _cut(lysozyme, "A", 55, 58)
    loops.forget()
    loops.close_gaps(s, library, 8)
    charge = s.charge.copy()
    charge[0] = 3
    edited = replace(s, charge=charge)
    assert loops.fingerprint(edited) == loops.fingerprint(s)  # the same key
    out, _ = loops.close_gaps(edited, library, 8)
    first = out.uid == s.uid[0]
    assert out.charge[first].tolist() == [3]
