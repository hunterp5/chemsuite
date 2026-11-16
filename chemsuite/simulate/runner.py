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

"""Driving the OpenMM worker (:mod:`.worker`) from the app.

One warm worker per interpreter (:func:`worker_for`), so OpenMM, OpenFF and the force fields
load once per session and a system built for one request is reused by the next. A cancelled
minimisation is asked to stop between iterations and comes back with where it got to.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.lines import LineWorker, WarmWorkers, WorkerError, child_environment

from . import worker as worker_script
from .config import Interpreter
from .spec import ForceFieldSpec, MinimiseSpec, SystemSpec

#: seconds the worker may stay silent: starting (OpenMM and its plugins load), answering hello
#: (a context is made on every platform, which compiles kernels on a GPU the first time), and
#: computing (a build of a large system, or a stage's context, sends nothing meanwhile)
START_TIMEOUT_S = 180.0
HELLO_TIMEOUT_S = 300.0
WORK_TIMEOUT_S = 1800.0

#: what minimisation reports as it goes: (stage, iteration, energy in kJ/mol)
Progress = Callable[[str, int, float], None]


class SimulateError(WorkerError):
    """The OpenMM worker could not be started, died, or could not do what was asked."""


@dataclass(frozen=True)
class Platform:
    name: str
    works: bool
    detail: str  # the device, or why a context could not be made


@dataclass(frozen=True)
class Hello:
    python: str
    versions: dict[str, str]
    platforms: tuple[Platform, ...]

    @property
    def missing(self) -> list[str]:
        """The packages a minimisation needs that the environment lacks."""
        need = ("openmmforcefields", "openff-toolkit", "openff-nagl", "openff-nagl-models", "rdkit")
        return [name for name in need if not self.versions.get(name)]


@dataclass(frozen=True)
class StageRun:
    name: str
    energy_before: float  # kJ/mol, restraints left out
    energy_after: float
    iterations: int
    seconds: float


@dataclass(frozen=True)
class Minimised:
    xyz: np.ndarray  # (n, 3) Å, in the spec's atom order
    stages: tuple[StageRun, ...]
    energy_before: float
    energy_after: float
    max_force: float  # kJ/mol/nm
    max_force_atom: int  # an index into the spec's atoms
    platform: str
    stopped: bool  # asked to stop before the last stage ended
    rmsd: float = 0.0  # Å, the held heavy atoms from where they started


class Simulator(LineWorker):
    """One OpenMM worker process."""

    error = SimulateError

    def __init__(self, where: Interpreter) -> None:
        self.where = where
        super().__init__(
            worker_script.command(where.python, where.wsl),
            name="The OpenMM worker",
            env=None if where.wsl else child_environment(where.python),
            start_timeout=START_TIMEOUT_S,
        )

    def hello(self, cancel: CancelToken | None = None) -> Hello:
        """The versions and the platforms that work. Builds nothing."""
        reply = self._ok(self.request("hello", timeout=HELLO_TIMEOUT_S, cancel=cancel))
        platforms = tuple(Platform(str(p[0]), bool(p[1]), str(p[2])) for p in reply["platforms"])
        versions = {str(k): str(v) for k, v in reply["versions"].items()}
        found = Hello(str(reply["python"]), versions, platforms)
        _SEEN[self.where] = found
        return found

    def charges(
        self, system: SystemSpec, forcefield: ForceFieldSpec, cancel: CancelToken | None = None
    ) -> SystemSpec:
        """*system* with every small molecule given its partial charges (those it has kept)."""
        if system.charged:
            return system
        body = {"system": system.to_json(), "small_molecule": forcefield.small_molecule}
        reply = self._ok(self.request("charges", body, timeout=WORK_TIMEOUT_S, cancel=cancel))
        for found in reply["charges"]:
            system = system.with_charges(
                int(found["residue"]), tuple(float(q) for q in found["charges"]),
                str(found["model"]),
            )  # fmt: skip
        return system

    def build(
        self, system: SystemSpec, forcefield: ForceFieldSpec, cancel: CancelToken | None = None
    ) -> int:
        """Build the OpenMM system (kept by the worker); its number of particles. Raises
        :class:`~.build.BuildError` with what the force field could not cover."""
        body = {"system": system.to_json(), "forcefield": forcefield.to_json()}
        reply = self._ok(self.request("build", body, timeout=WORK_TIMEOUT_S, cancel=cancel))
        return int(reply["particles"])

    def minimise(
        self,
        system: SystemSpec,
        forcefield: ForceFieldSpec,
        xyz: np.ndarray,
        protocol: MinimiseSpec | None = None,
        *,
        platform: str = "",
        free: Sequence[int] = (),
        progress: Progress | None = None,
        cancel: CancelToken | None = None,
    ) -> Minimised:
        """Minimise *system* from *xyz* (Å, the spec's order), stage by stage, every heavy atom
        held to where it started but the *free* ones (indices: the atoms Prepare built).
        Cancelling stops it between iterations: :class:`Cancelled` carries what it got to (or
        None)."""
        body = {
            "system": system.to_json(),
            "forcefield": forcefield.to_json(),
            "positions": np.asarray(xyz, dtype=np.float64).reshape(-1).tolist(),
            "protocol": (protocol or MinimiseSpec()).to_json(),
            "platform": platform,
            "free": sorted(int(i) for i in free),
        }

        def event(found: dict) -> None:
            if progress is not None:
                progress(str(found["stage"]), int(found["iteration"]), float(found["energy"]))

        try:
            reply = self.request("minimise", body, timeout=WORK_TIMEOUT_S, cancel=cancel,
                                 events=event, stoppable=True)  # fmt: skip
        except Cancelled as stopped:
            partial = stopped.partial
            raise Cancelled(
                _minimised(partial) if isinstance(partial, dict) and partial.get("ok") else None
            ) from None
        return _minimised(self._ok(reply))

    def _ok(self, reply: dict) -> dict:
        if reply.get("ok"):
            return reply
        problems = reply.get("problems")
        if isinstance(problems, list) and problems:
            from .build import BuildError, Problem  # (the protein chemistry: not at start-up)

            raise BuildError([Problem(str(p["label"]), str(p["text"])) for p in problems])
        raise SimulateError(str(reply.get("error") or "The OpenMM worker gave no reason."))


def _minimised(reply: dict) -> Minimised:
    stages = tuple(
        StageRun(str(s["name"]), float(s["energy_before"]), float(s["energy_after"]),
                 int(s["iterations"]), float(s["seconds"]))
        for s in reply["stages"]
    )  # fmt: skip
    return Minimised(
        xyz=np.asarray(reply["positions"], dtype=np.float64).reshape(-1, 3),
        stages=stages,
        energy_before=float(reply["energy_before"]),
        energy_after=float(reply["energy_after"]),
        max_force=float(reply["max_force"]),
        max_force_atom=int(reply["max_force_atom"]),
        platform=str(reply["platform"]),
        stopped=bool(reply["stopped"]),
        rmsd=float(reply.get("rmsd", 0.0)),
    )


# -- the warm workers ---------------------------------------------------------------------------
_WORKERS: WarmWorkers[Simulator] = WarmWorkers(Simulator)
#: what each interpreter's worker last said it has (Settings' Test, or a run), this session
_SEEN: dict[Interpreter, Hello] = {}


def platforms_seen(where: Interpreter) -> tuple[str, ...]:
    """The platforms a worker in *where* has said work, without starting one (empty until a
    worker there has said hello this session)."""
    found = _SEEN.get(where)
    return tuple(p.name for p in found.platforms if p.works) if found is not None else ()


def worker_for(where: Interpreter) -> Simulator:
    """The warm worker for this interpreter, started on first use."""
    return _WORKERS.get(where)


def close_workers() -> None:
    """Stop every warm worker (the app closing, or the paths changing in Settings)."""
    _WORKERS.close_all()


def probe(where: Interpreter) -> tuple[bool, str]:
    """``(ok, message)`` for the Settings ▸ External Tools **Test** button: a throwaway worker
    reports the versions and the platforms that work."""
    from chemsuite.core.jobs.lines import is_app_python

    if not where.python:
        return False, "Set the interpreter of the OpenMM environment."
    if not where.wsl and is_app_python(where.python):
        return False, "That is chemsuite's own interpreter; OpenMM needs its own."
    worker = Simulator(where)
    try:
        found = worker.hello()
    except SimulateError as err:
        return False, str(err).splitlines()[0]
    finally:
        worker.close()
    if found.missing:
        return False, f"The environment lacks {', '.join(found.missing)}."
    working = [p.name + (f" ({p.detail})" if p.detail else "") for p in found.platforms if p.works]
    return True, (
        f"OpenMM {found.versions['openmm']}, openmmforcefields "
        f"{found.versions['openmmforcefields']}, OpenFF {found.versions['openff-toolkit']} "
        f"on Python {found.python}: {', '.join(working) or 'no working platform'}."
    )
