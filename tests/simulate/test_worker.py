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

"""The OpenMM worker in a real OpenMM environment.

These run only where ``CHEMSUITE_OPENMM_PYTHON`` names the interpreter of one (the env of
``envs/openmm-*.yml``); ``pytest -m openmm`` runs them alone.
"""

from __future__ import annotations

import os
import threading
from dataclasses import replace

import numpy as np
import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.simulate import build, spec
from chemsuite.simulate.config import Interpreter
from chemsuite.simulate.minimise import free_atoms, measure
from chemsuite.simulate.runner import Simulator
from chemsuite.simulate.spec import ForceFieldSpec, MinimiseSpec, Stage

from .test_build import _without

PYTHON = os.environ.get("CHEMSUITE_OPENMM_PYTHON", "")
pytestmark = [
    pytest.mark.openmm,
    pytest.mark.skipif(not PYTHON, reason="CHEMSUITE_OPENMM_PYTHON names no OpenMM interpreter"),
]
#: a short protocol, so the tests are quick: the stages of the real one, few iterations each
QUICK = MinimiseSpec(max_iterations=200)


@pytest.fixture(scope="module")
def worker():
    found = Simulator(Interpreter(PYTHON))
    yield found
    found.close()


@pytest.fixture
def complex_(lysozyme):
    """T4 lysozyme with its ligand, its phosphates and BME left out (their bond orders are
    guesses)."""
    s, states = _without(lysozyme, ("PO4", "BME"))
    return s, build.system_spec(s, states)


def test_hello_names_the_versions_and_a_working_platform(worker) -> None:
    found = worker.hello()
    assert found.versions["openmm"] and not found.missing
    assert any(p.works for p in found.platforms if p.name in ("CPU", "Reference"))


def test_a_ligand_is_charged_once_and_the_complex_built(worker, complex_) -> None:
    s, system = complex_
    forcefield = ForceFieldSpec()
    with pytest.raises(build.BuildError, match="no partial charges"):
        worker.build(system, forcefield)
    charged = worker.charges(system, forcefield)
    (molecule,) = charged.molecules
    assert molecule.charge_model.startswith("openff-gnn-am1bcc")
    formal = sum(charged.formal[charged.residues[molecule.residue].first + i]
                 for i in range(len(molecule.charges)))  # fmt: skip
    assert abs(sum(molecule.charges) - formal) < 1e-3
    assert charged.digest() == system.digest()
    assert worker.build(charged, forcefield) == system.n_atoms


def test_a_charge_the_force_field_disagrees_with_blocks_the_build(worker, complex_) -> None:
    s, system = complex_
    charged = worker.charges(system, ForceFieldSpec())
    lys = next(r for r in charged.residues if r.name == "LYS")
    formal = list(charged.formal)
    nz = next(i for i in range(lys.first, lys.first + lys.count) if charged.name[i] == "NZ")
    formal[nz] = 0  # Prepare would have said LYN, with a hydrogen fewer
    wrong = spec.SystemSpec(charged.uid, charged.element, charged.name, tuple(formal),
                            charged.residues, charged.bonds, charged.disulfides,
                            charged.molecules)  # fmt: skip
    with pytest.raises(build.BuildError) as raised:
        worker.build(wrong, ForceFieldSpec())
    assert [p.label for p in raised.value.problems] == [lys.label]


def test_minimising_moves_hydrogens_first_and_holds_heavy_atoms_near(worker, complex_) -> None:
    s, system = complex_
    charged = worker.charges(system, ForceFieldSpec())
    xyz = build.positions(s, charged)
    events = []
    done = worker.minimise(charged, ForceFieldSpec(), xyz, QUICK, platform="CPU",
                           progress=lambda *e: events.append(e))  # fmt: skip
    assert [r.name for r in done.stages] == [s.name for s in QUICK.stages] and not done.stopped
    assert done.energy_after < done.energy_before
    assert {e[0] for e in events} == {s.name for s in QUICK.stages}
    heavy = np.array(charged.element) != 1
    moved = np.linalg.norm(done.xyz - xyz, axis=1)
    # the plan's mark is the heavy-atom RMSD; single surface side-chain tips swing further in
    # the side-chain stage (Glu 108's OE1 about 1.1 Å), which M3's protocol tuning looks at
    assert np.sqrt((moved[heavy] ** 2).mean()) < 0.3 and moved[heavy].max() < 1.5
    assert moved[~heavy].max() > 0.05  # the hydrogens settled


def test_a_hydrogens_only_stage_leaves_every_heavy_atom_where_it_was(worker, complex_) -> None:
    s, system = complex_
    charged = worker.charges(system, ForceFieldSpec())
    xyz = build.positions(s, charged)
    only = MinimiseSpec(stages=(Stage("Hydrogens", spec.HYDROGENS),), max_iterations=100)
    done = worker.minimise(charged, ForceFieldSpec(), xyz, only, platform="CPU")
    heavy = np.array(charged.element) != 1
    assert np.abs(done.xyz[heavy] - xyz[heavy]).max() < 1e-4


def test_a_cancelled_minimisation_stops_and_says_where_it_got(worker, complex_) -> None:
    s, system = complex_
    charged = worker.charges(system, ForceFieldSpec())
    cancel = CancelToken()
    long = MinimiseSpec(tolerance=0.001, max_iterations=0)  # would run a long while
    seen = threading.Event()

    def progress(stage: str, iteration: int, energy: float) -> None:
        if iteration >= 50:
            seen.set()
            cancel.cancel()

    with pytest.raises(Cancelled) as raised:
        worker.minimise(charged, ForceFieldSpec(), build.positions(s, charged), long,
                        platform="CPU", progress=progress, cancel=cancel)  # fmt: skip
    assert seen.is_set() and worker.alive  # stopped, not killed: the worker stays warm
    partial = raised.value.partial
    assert partial is not None and partial.stopped and len(partial.stages) == 1
    assert worker.hello().versions["openmm"]


def test_built_atoms_go_free_and_the_worker_measures_as_the_app_does(worker, lysozyme) -> None:
    s, states = _without(lysozyme, ("PO4", "BME"))
    charged = worker.charges(build.system_spec(s, states), ForceFieldSpec())
    xyz = build.positions(s, charged)
    free = free_atoms(charged, lysozyme.built)
    done = worker.minimise(charged, ForceFieldSpec(), xyz, QUICK, platform="CPU",
                           free=np.flatnonzero(free).tolist())  # fmt: skip
    found = measure(charged, xyz, done, free)
    assert done.rmsd == pytest.approx(found.rmsd, abs=1e-3) and found.rmsd < 0.3
    assert found.contact_after is not None and found.contact_after.distance > 2.2


def test_the_last_stage_is_tightened_while_the_heavy_atoms_end_too_far(worker, complex_) -> None:
    s, system = complex_
    charged = worker.charges(system, ForceFieldSpec())
    xyz = build.positions(s, charged)
    once = worker.minimise(charged, ForceFieldSpec(), xyz,
                           replace(QUICK, rmsd_limit=0.01, tighten=0), platform="CPU")  # fmt: skip
    tight = worker.minimise(charged, ForceFieldSpec(), xyz,
                            replace(QUICK, rmsd_limit=0.01, tighten=2), platform="CPU")  # fmt: skip
    assert [r.name for r in tight.stages][-2:] == ["Complex, restraints x2",
                                                   "Complex, restraints x4"]  # fmt: skip
    assert tight.rmsd < once.rmsd
    # (within the limit nothing is tightened: the default protocol above ran its three stages)
