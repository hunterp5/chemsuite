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

"""Valence and implicit hydrogens (SKETCHER_PLAN §4.1; docs/sketcher/valence-review.md).

Reference data: the 61 species of the P0 review, each with a hydrogen count that is not in
doubt. RDKit's two known gaps (ClF₃ and IF₇) are expected failures.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.sketch import ctfile, valence
from chemsuite.sketch.convert import to_molblock
from tests.sketch.helpers import drawn, valence_review

REVIEW = valence_review()
KNOWN_GAPS = {"chlorine trifluoride", "iodine heptafluoride"}


def _cases():
    for smiles, index, want, species in REVIEW.REFERENCE:
        marks = (
            [pytest.mark.xfail(strict=True, reason="RDKit rejects")]
            if species in KNOWN_GAPS
            else []
        )
        yield pytest.param(smiles, index, want, id=species, marks=marks)


@pytest.mark.parametrize(("smiles", "index", "want"), list(_cases()))
def test_reference_species_hydrogens(smiles, index, want):
    sketch, ids = drawn(smiles)
    report = valence.check(sketch)
    assert report.valid, report.problems
    assert report.hydrogens[ids[index]] == want


def test_reference_species_score_matches_the_review():
    right = 0
    for smiles, index, want, _ in REVIEW.REFERENCE:
        sketch, ids = drawn(smiles)
        report = valence.check(sketch)
        right += report.valid and report.hydrogens[ids[index]] == want
    assert right == 59  # docs/sketcher/valence-review.md: RDKit 59/61


def _ctfile_hydrogens(molblock: str) -> list[int | None]:
    """Hydrogens a CTfile reader finds: the valence field if set, else the default model.

    The default model is the line-by-line port of Ketcher's ``Atom.calcValence`` in
    scripts/review/valence_review.py, standing in for other programs' molfile readers.
    """
    lines = molblock.splitlines()
    natoms, nbonds = int(lines[3][0:3]), int(lines[3][3:6])
    atoms = lines[4 : 4 + natoms]
    conn = [0] * natoms
    for line in lines[4 + natoms : 4 + natoms + nbonds]:
        a, b, order = int(line[0:3]) - 1, int(line[3:6]) - 1, int(line[6:9])
        conn[a] += order
        conn[b] += order
    charge, rad = [0] * natoms, [0] * natoms
    for line in lines:
        if line.startswith(("M  CHG", "M  RAD")):
            values = [int(v) for v in line[9:].split()]
            for i, v in zip(values[::2], values[1::2], strict=True):
                if line.startswith("M  CHG"):
                    charge[i - 1] = v
                else:
                    rad[i - 1] = {1: 2, 2: 1, 3: 2}[v]
    out: list[int | None] = []
    for i, line in enumerate(atoms):
        symbol, vvv = line[31:34].strip(), int(line[48:51])
        if vvv:
            out.append(0 if vvv == 15 else vvv - conn[i])
        else:
            out.append(ctfile.default_hydrogens(symbol, charge[i], conn[i], rad[i]))
    return out


@pytest.mark.parametrize(
    ("smiles", "index", "want"),
    [
        pytest.param(s, i, w, id=name)
        for s, i, w, name in REVIEW.REFERENCE
        if name not in KNOWN_GAPS
    ],
)
def test_exported_molfile_reads_the_same_in_ctfile_programs(smiles, index, want):
    """Where the CTfile model would differ, the export carries an explicit valence."""
    sketch, ids = drawn(smiles)
    block = to_molblock(sketch)
    report = valence.check(sketch)

    mol = Chem.MolFromMolBlock(block, removeHs=False)
    theirs = _ctfile_hydrogens(block)
    for atom in mol.GetAtoms():
        assert theirs[atom.GetIdx()] == atom.GetTotalNumHs(), (atom.GetIdx(), atom.GetSymbol())
    # The exported file keeps the reference count (nitro may be rewritten charge-separated).
    assert report.hydrogens[ids[index]] == want


def test_bad_valence_is_marked_with_a_reason_not_refused():
    sketch, ids = drawn("CC(C)(C)(C)C")
    report = valence.check(sketch)
    assert not report.valid
    assert set(report.problems) == {ids[1]}
    assert "Too many bonds" in report.problems[ids[1]]
    assert "C with 5" in report.problems[ids[1]]
    assert report.hydrogens[ids[0]] == 3  # the other atoms still get their hydrogens


def test_unkekulizable_aromatic_ring_is_marked():
    sketch, ids = drawn("c1cccc1")
    report = valence.check(sketch)
    assert not report.valid
    assert set(report.problems) <= set(ids) and report.problems


def test_empty_sketch_is_valid():
    from chemsuite.sketch.model import Sketch

    assert valence.check(Sketch()).valid
