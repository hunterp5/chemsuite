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

"""Named cases from the literature, each on a pocket cut from its entry, with the defaults."""

from __future__ import annotations

import pytest

from chemsuite.protein.prepare import pipeline
from chemsuite.protein.structure.model import Structure
from tests.protein.conftest import load
from tests.protein.prepare.test_hydrogens import _on


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


def test_acetazolamide_binds_zinc_as_the_sulfonamide_anion() -> None:
    """Carbonic anhydrase II with acetazolamide (4G0C, neutron): N1 keeps one D, on the Zn."""
    done = _prepared("4g0c_pocket.cif.gz", 7.8)
    s = done.structure
    r = _residue(s, "AZM")
    assert len(_on(s, r, "N1")) == 1 and _charge(s, r, "N1") == -1
    assert done.states[r] == "−H N1"


def test_benzamidine_in_trypsin_is_the_amidinium() -> None:
    """Trypsin with benzamidine (5MO0, neutron): both N carry two D, the charge on Asp189."""
    done = _prepared("5mo0_pocket.cif.gz", 7.5)
    s = done.structure
    r = _residue(s, "BEN")
    assert len(_on(s, r, "N1")) == 2 and len(_on(s, r, "N2")) == 2
    assert _charge(s, r, "N1") + _charge(s, r, "N2") == 1


@pytest.fixture(scope="module")
def zinc_finger() -> pipeline.Prepared:
    """Zif268's three Cys₂His₂ fingers (1ZAA)."""
    return _prepared("1zaa_pocket.cif.gz", 7.0)


@pytest.mark.parametrize("seq", [7, 12, 37, 40, 65, 68])
def test_zinc_finger_cysteines_are_thiolates(zinc_finger, seq) -> None:
    assert zinc_finger.states[_residue(zinc_finger.structure, "CYS", seq)] == "CYM"
