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

"""Named cases from the literature, each on a pocket cut from its entry, with the defaults
(the ligand cases, which need chemsuite pKa, are in :mod:`.test_ligand_pka`)."""

from __future__ import annotations

import pytest

from chemsuite.macromolecule.prepare import pipeline
from chemsuite.macromolecule.structure.model import Structure
from tests.macromolecule.conftest import load


def _prepared(name: str, ph: float, **options) -> pipeline.Prepared:
    s = load(name)
    s = s.subset(s.element != 1)  # neutron entries: their H and D atoms off
    return pipeline.prepare(s, pipeline.Options(ph=ph, **options))


def _residue(s: Structure, name: str, seq: int | None = None) -> int:
    (r,) = [r for r in range(s.n_residues) if str(s.res_name[r]) == name
            and (seq is None or int(s.res_seq[r]) == seq)]  # fmt: skip
    return r


def _charge(s: Structure, r: int, atom: str) -> int:
    span = s.atoms_of(r)
    (i,) = [i for i in range(span.start, span.stop) if s.name[i] == atom]
    return int(s.charge[i])


@pytest.fixture(scope="module")
def zinc_finger() -> pipeline.Prepared:
    """Zif268's three Cys₂His₂ fingers (1ZAA)."""
    return _prepared("1zaa_pocket.cif.gz", 7.0)


@pytest.mark.parametrize("seq", [7, 12, 37, 40, 65, 68])
def test_zinc_finger_cysteines_are_thiolates(zinc_finger, seq) -> None:
    assert zinc_finger.states[_residue(zinc_finger.structure, "CYS", seq)] == "CYM"
