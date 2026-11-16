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

"""A minimisation measured: RMSDs with built atoms left out, contacts, what moved most."""

from __future__ import annotations

import sys

import numpy as np
import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.simulate import build
from chemsuite.simulate.config import Interpreter
from chemsuite.simulate.minimise import (
    Done,
    closest_contact,
    free_atoms,
    measure,
    minimise,
)
from chemsuite.simulate.runner import Minimised, Simulator
from chemsuite.simulate.spec import Residue, SystemSpec

from .test_build import _without
from .test_runner import FAKE


def _index(system: SystemSpec, label: str) -> int:
    residue_of = system.residue_of()
    return next(i for i in range(system.n_atoms)
                if f"{system.residues[residue_of[i]].label}/{system.name[i]}" == label)  # fmt: skip


def _done(xyz: np.ndarray, worst: int = 0) -> Minimised:
    return Minimised(xyz, (), -100.0, -200.0, 9.0, worst, "CPU", False)


def test_built_atoms_are_free_and_left_out_of_what_moved(lysozyme) -> None:
    s, states = _without(lysozyme, ("PO4", "BME"))
    system = build.system_spec(s, states)
    free = free_atoms(system, lysozyme.built)
    heavy = np.array(system.element) != 1
    assert free.any() and not (free & ~heavy).any()  # the NME cap's; never a hydrogen
    cap = next(r for r in system.residues if r.name == "NME")
    assert free[cap.first : cap.first + cap.count][heavy[cap.first : cap.first + cap.count]].all()
    start = build.positions(s, system)
    end = start.copy()
    nz = _index(system, "A/LYS 124/NZ")
    end[nz, 0] += 1.0
    end[np.flatnonzero(free)[0], 0] += 5.0  # a built atom going far is not a deposited one moving
    found = measure(system, start, _done(end, nz), free, (("openmm", "8.4.0"),))
    held = heavy & ~free
    assert found.rmsd == pytest.approx(np.sqrt(1.0 / held.sum()))
    assert found.moved_most[0] == ("A/LYS 124", pytest.approx(1.0))
    assert all(label != cap.label for label, _d in found.moved_most)
    assert found.rmsd_ligand == 0.0 and found.rmsd_pocket is not None
    assert found.max_force_atom == "A/LYS 124/NZ" and found.digest == system.digest()
    assert found.versions == (("openmm", "8.4.0"),) and found.energy_after == -200.0


def test_no_ligand_no_pocket(crambin) -> None:
    system = build.system_spec(crambin.structure, crambin.states)
    start = build.positions(crambin.structure, system)
    found = measure(system, start, _done(start), free_atoms(system, crambin.built))
    assert found.rmsd == 0.0 and found.rmsd_pocket is None and found.rmsd_ligand is None
    assert found.moved_most == ()


def test_the_closest_contact_skips_bonded_neighbours(crambin) -> None:
    system = build.system_spec(crambin.structure, crambin.states)
    xyz = build.positions(crambin.structure, system)
    found = closest_contact(system, xyz)
    # an amide O and the next residue's N are 2.25 Å apart, and a disulfide's SGs bonded
    assert found is not None and found.distance > 2.5
    a, b = _index(system, "A/THR 1/OG1"), _index(system, "A/ASN 46/CB")
    planted = xyz.copy()
    planted[b] = planted[a] + [1.5, 0.0, 0.0]
    clash = closest_contact(system, planted)
    assert clash is not None and {clash.a, clash.b} == {"A/THR 1/OG1", "A/ASN 46/CB"}
    assert clash.distance == pytest.approx(1.5, abs=1e-6)


def test_an_ion_beside_a_water_is_no_contact() -> None:
    residues = (Residue("HOH", "A", 1, "", "water", "A/HOH 1", 0, 3),
                Residue("NA", "A", 2, "", "ion", "A/NA 2", 3, 1))  # fmt: skip
    system = SystemSpec((1, 2, 3, 4), (8, 1, 1, 11), ("O", "H1", "H2", "NA"), (0, 0, 0, 1),
                        residues, ((0, 1, 1), (0, 2, 1)))  # fmt: skip
    xyz = np.array([[0.0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0], [1.0, 1.0, 1.0]])
    assert closest_contact(system, xyz) is None


@pytest.fixture
def fake():
    found = Simulator(Interpreter(sys.executable))
    found.command = [sys.executable, "-I", "-c", FAKE]
    found.env = None
    yield found
    found.close()


def test_minimise_runs_the_worker_and_measures_what_it_gives_back(fake, crambin) -> None:
    s = crambin.structure
    system = build.system_spec(s, crambin.states)
    done = minimise(fake, s, system, built=crambin.built)
    assert isinstance(done, Done) and done.system == system  # no small molecule to charge
    result = done.result
    assert result.rmsd == pytest.approx(np.sqrt(3 * 0.1**2), abs=1e-6)  # the fake moves all
    assert result.platform == "Reference" and dict(result.versions)["openmm"] == "8.4"
    assert [name for name, _s in result.seconds] == ["Hydrogens", "Side chains and waters",
                                                    "Complex"]  # fmt: skip
    assert np.allclose(result.coordinates(s), s.xyz() + 0.1, atol=1e-3)


def test_a_cancelled_minimisation_is_measured_as_far_as_it_got(fake, crambin) -> None:
    s = crambin.structure
    system = build.system_spec(s, crambin.states)
    cancel = CancelToken()
    with pytest.raises(Cancelled) as raised:
        minimise(fake, s, system, progress=lambda *_e: cancel.cancel(), cancel=cancel)
    partial = raised.value.partial
    assert isinstance(partial, Done) and len(partial.result.seconds) == 1
