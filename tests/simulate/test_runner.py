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

"""Where OpenMM runs, and what the app makes of the worker's replies, without OpenMM: the
replies come from a stand-in worker (the real one is tested in ``test_worker.py``)."""

from __future__ import annotations

import sys

import numpy as np
import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.settings import Settings
from chemsuite.simulate import config, runner, worker
from chemsuite.simulate.build import BuildError
from chemsuite.simulate.config import Interpreter
from chemsuite.simulate.runner import SimulateError, Simulator
from chemsuite.simulate.spec import ForceFieldSpec, Molecule, Residue, SystemSpec

#: answers as the OpenMM worker does, with made-up values
FAKE = r"""
import json, sys, threading, time
stop = threading.Event()
work = []

def send(payload):
    sys.stdout.write(json.dumps(payload) + "\n")
    sys.stdout.flush()

def read():
    for line in sys.stdin:
        request = json.loads(line)
        if request["op"] == "stop":
            stop.set()
        else:
            work.append(request)
    work.append(None)

threading.Thread(target=read, daemon=True).start()
send({"id": 0, "ok": True, "ready": True})
while True:
    while not work:
        time.sleep(0.01)
    request = work.pop(0)
    if request is None or request["op"] == "quit":
        sys.exit(0)
    ident, op = request["id"], request["op"]
    stop.clear()
    if op == "hello":
        send({"id": ident, "ok": True, "python": "3.12.16", "versions": {"openmm": "8.4",
              "openmmforcefields": "0.15.1", "openff-toolkit": "0.19.0", "openff-nagl": "0.6.1",
              "openff-nagl-models": "", "rdkit": "2026.09.2", "numpy": "2"},
              "platforms": [["Reference", True, ""], ["CUDA", False, "No compatible CUDA device"]]})
    elif op == "charges":
        model = "openff-gnn-am1bcc-1.0.0"
        found = [{"residue": m["residue"], "charges": [-0.5, 0.5], "model": model}
                 for m in request["system"]["molecules"] if not m["charges"]]
        send({"id": ident, "ok": True, "charges": found})
    elif op == "build":
        send({"id": ident, "ok": False, "error": "Problems: ...",
              "problems": [{"label": "A/TYR 29", "text": "ff14SB has no TYM template"}]})
    elif op == "minimise":
        n = len(request["positions"]) // 3
        stages = []
        for stage in request["protocol"]["stages"]:
            for i in range(0, 100, 25):
                send({"id": ident, "event": {"stage": stage["name"], "iteration": i, "energy": -i}})
                time.sleep(0.02)
                if stop.is_set():
                    break
            stages.append({"name": stage["name"], "energy_before": 0.0, "energy_after": -1.0,
                           "iterations": 100, "seconds": 0.1})
            if stop.is_set():
                break
        send({"id": ident, "ok": True, "positions": [p + 0.1 for p in request["positions"]],
              "stages": stages, "energy_before": 0.0, "energy_after": -3.0, "max_force": 9.0,
              "max_force_atom": n - 1, "platform": request["platform"] or "Reference",
              "stopped": stop.is_set()})
    else:
        send({"id": ident, "ok": False, "error": "unknown " + op})
"""


@pytest.fixture
def fake():
    found = Simulator(Interpreter(sys.executable))
    found.command = [sys.executable, "-I", "-c", FAKE]
    found.env = None
    yield found
    found.close()


def _system() -> SystemSpec:
    """Water and a two-atom 'molecule': enough for the runner, which reads no chemistry."""
    residues = (Residue("HOH", "A", 1, "", "water", "A/HOH 1", 0, 3),
                Residue("LIG", "A", 2, "", "molecule", "A/LIG 2", 3, 2))  # fmt: skip
    return SystemSpec((1, 2, 3, 4, 5), (8, 1, 1, 6, 8), ("O", "H1", "H2", "C1", "O1"),
                      (0, 0, 0, 0, 0), residues, ((0, 1, 1), (0, 2, 1), (3, 4, 2)),
                      molecules=(Molecule(1, "molblock"),))  # fmt: skip


def test_the_native_interpreter_wins_and_none_is_looked_up_on_path(tmp_path) -> None:
    python = tmp_path / "python.exe"
    python.write_text("", encoding="utf-8")
    assert config.interpreter(Settings()) is None
    assert "not set up" in config.ready_problem(Settings())
    both = Settings(tool_paths={config.PYTHON_TOOL: str(python), config.WSL_TOOL: "/usr/bin/py"})
    assert config.interpreter(both) == Interpreter(str(python))
    assert config.ready_problem(both) == ""
    own = Settings(tool_paths={config.PYTHON_TOOL: sys.executable})
    assert "running chemsuite" in config.ready_problem(own)


@pytest.mark.skipif(sys.platform != "win32", reason="WSL is a Windows path")
def test_on_windows_a_wsl_interpreter_is_used_when_no_native_one_is_set() -> None:
    wsl = Settings(tool_paths={config.WSL_TOOL: "/home/me/miniconda3/envs/sim/bin/python"})
    found = config.interpreter(wsl)
    assert found == Interpreter("/home/me/miniconda3/envs/sim/bin/python", wsl=True)
    assert worker.command(found.python, found.wsl)[:3] == ["wsl.exe", "-e", found.python]
    assert Simulator(found).env is None  # WSL has its own environment
    browsed = r"\\wsl.localhost\Ubuntu\home\me\miniconda3\envs\sim\bin\python"
    assert config.wsl_python(Settings(tool_paths={config.WSL_TOOL: browsed})) == (
        "/home/me/miniconda3/envs/sim/bin/python"
    )
    assert config.linux_path(r"\\wsl$\Ubuntu\opt\py") == "/opt/py"


def test_the_wsl_set_in_preferences_runs_a_wsl_interpreter(monkeypatch) -> None:
    monkeypatch.setattr(config.sys, "platform", "win32")
    chosen = Settings(tool_paths={config.WSL_TOOL: "/opt/py", config.WSL_EXE_TOOL: "D:/wsl.exe"})
    found = config.interpreter(chosen)
    assert found == Interpreter("/opt/py", wsl=True, launcher="D:/wsl.exe")
    argv = worker.command(found.python, found.wsl, found.launcher)
    assert argv[:3] == ["D:/wsl.exe", "-e", "/opt/py"]
    plain = config.interpreter(Settings(tool_paths={config.WSL_TOOL: "/opt/py"}))
    assert plain is not None and plain.launcher == "wsl.exe"


def test_the_worker_script_is_python_this_interpreter_reads() -> None:
    compile(worker.SCRIPT, "<openmm worker>", "exec")
    assert worker.command("python")[:3] == ["python", "-I", "-c"]


def test_the_worker_script_fits_on_a_windows_command_line() -> None:
    """The script travels in ``-c``, and Windows caps a whole command line at 32,767
    characters (``WinError 206`` beyond). When this fails, start the worker with a short
    bootstrap in ``-c`` that reads the script from its first stdin line instead."""
    import subprocess

    line = subprocess.list2cmdline(worker.command(sys.executable, wsl=True))
    assert len(line) < 24_000, len(line)


def test_hello_reports_what_is_missing_and_which_platforms_work(fake) -> None:
    found = fake.hello()
    assert found.missing == ["openff-nagl-models"]
    assert [(p.name, p.works) for p in found.platforms] == [("Reference", True), ("CUDA", False)]


def test_charges_are_written_into_the_spec(fake) -> None:
    system = _system()
    charged = fake.charges(system, ForceFieldSpec())
    assert charged.molecules[0].charges == (-0.5, 0.5) and charged.charged
    assert fake.charges(charged, ForceFieldSpec()) is charged  # nothing left to ask for


def test_what_the_force_field_cannot_cover_comes_back_as_problems(fake) -> None:
    with pytest.raises(BuildError) as raised:
        fake.build(_system(), ForceFieldSpec())
    assert [str(p) for p in raised.value.problems] == ["A/TYR 29: ff14SB has no TYM template"]
    with pytest.raises(SimulateError, match="unknown"):
        fake._ok(fake.request("nonsense", timeout=10))


def test_a_minimisation_reports_its_stages_and_comes_back_in_angstroms(fake) -> None:
    xyz = np.arange(15, dtype=float).reshape(5, 3)
    seen = []
    done = fake.minimise(_system(), ForceFieldSpec(), xyz, progress=lambda *e: seen.append(e))
    assert np.allclose(done.xyz, xyz + 0.1) and not done.stopped
    assert [s.name for s in done.stages] == ["Hydrogens", "Side chains and waters", "Complex"]
    assert seen[0] == ("Hydrogens", 0, 0.0) and len(seen) == 12
    assert done.max_force_atom == 4 and done.platform == "Reference"


def test_a_cancelled_minimisation_brings_back_what_it_did(fake) -> None:
    cancel = CancelToken()
    xyz = np.zeros((5, 3))

    def progress(stage: str, iteration: int, energy: float) -> None:
        if iteration == 25:
            cancel.cancel()

    with pytest.raises(Cancelled) as raised:
        fake.minimise(_system(), ForceFieldSpec(), xyz, progress=progress, cancel=cancel)
    partial = raised.value.partial
    assert isinstance(partial, runner.Minimised) and partial.stopped
    assert len(partial.stages) == 1 and fake.alive


def test_workers_are_kept_warm_per_interpreter_and_closed_together() -> None:
    a, b = Interpreter("/a/python"), Interpreter("/b/python", wsl=True)
    try:
        assert runner.worker_for(a) is runner.worker_for(a)
        assert runner.worker_for(b) is not runner.worker_for(a)
    finally:
        runner.close_workers()
    assert runner.worker_for(a) is not None
    runner.close_workers()


def test_the_platforms_a_worker_named_are_remembered_without_starting_one() -> None:
    from chemsuite.simulate import runner

    where = Interpreter(sys.executable)
    runner._SEEN.pop(where, None)
    assert runner.platforms_seen(where) == ()
    worker = Simulator(where)
    worker.command = [sys.executable, "-I", "-c", FAKE]
    worker.env = None
    try:
        worker.hello()
    finally:
        worker.close()
    assert "Reference" in runner.platforms_seen(where)
