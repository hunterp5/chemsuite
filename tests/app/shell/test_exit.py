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

"""Closing the last window ends everything: the app, its worker processes and the programs
it started, even with a job that never looks at its cancel token."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

from tests.gui import gui
from tests.processes import alive, kill

APP = r"""
import json, os, subprocess, sys, time
from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QLabel

from chemsuite.app.document.main_window import MainWindow
from chemsuite.app.shell import quitting
from chemsuite.app.shell.jobs import QUEUE, JobRunner
from chemsuite.app.shell.main import run_app
from chemsuite.app.shell.windows import ToolWindow
from chemsuite.app.table.cells import StructureImages
from chemsuite.core.jobs import check, children

PARENT = (
    "import subprocess, sys, time\n"
    "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(600)'])\n"
    "print(child.pid, flush=True)\n"
    "time.sleep(600)\n"
)
pids = {}


def deaf(job):  # stuck in a call that never sees its cancel token
    time.sleep(600)


def pooled(job):  # waits on a worker process busy for ten minutes
    pool = job.parent().process_pool()
    pids["worker"] = pool.submit(os.getpid).result(timeout=60)  # a worker is up
    busy = pool.submit(time.sleep, 600)
    pids["workers"] = [p.pid for p in pool._processes.values()]
    while not busy.done():
        check(job.cancel_token)
        time.sleep(0.05)


def background(job):
    while True:
        check(job.cancel_token)
        time.sleep(0.05)


def main():  # pool workers are spawned: they import this file without running it
    out = Path(sys.argv[1])
    app = QApplication([])
    runner = JobRunner(app)
    window = MainWindow(runner, StructureImages())
    quitting.confirm = lambda parent, question: True  # "Cancel them and quit?" answered Yes
    pids["app"] = os.getpid()
    program = children.start([sys.executable, "-c", PARENT], stdout=subprocess.PIPE, text=True)
    pids["program"] = program.pid
    pids["grandchild"] = int(program.stdout.readline())
    runner.submit("deaf", deaf, owner=window)
    runner.submit("pooled", pooled, lane=QUEUE)
    runner.submit("view", background).background = True
    window.show_log()
    tool = ToolWindow(QLabel("tool"), "Tool", name=lambda: "")
    tool.show()

    def ready():
        if "workers" not in pids:
            QTimer.singleShot(50, ready)
            return
        time.sleep(0.2)  # the worker is inside time.sleep(600)
        out.write_text(json.dumps(pids))
        window.close()

    QTimer.singleShot(200, ready)
    return run_app(app, runner, window)


if __name__ == "__main__":
    sys.exit(main())
"""


@gui
def test_closing_the_last_window_ends_everything(tmp_path) -> None:
    script = tmp_path / "app.py"
    script.write_text(APP)
    out = tmp_path / "pids.json"
    env = dict(
        os.environ,
        QT_QPA_PLATFORM="offscreen",
        CHEMSUITE_CONFIG_DIR=str(tmp_path / "config"),
        CHEMSUITE_CACHE_DIR=str(tmp_path / "cache"),
        PYTHONPATH=os.pathsep.join(
            [str(Path(__file__).parents[3]), os.environ.get("PYTHONPATH", "")]
        ),
    )
    with (tmp_path / "stderr.txt").open("w") as errors:
        app = subprocess.Popen(
            [sys.executable, str(script), str(out)],
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,  # not pipes: a lingering child would hold a pipe open
            stderr=errors,
        )
    deadline = time.monotonic() + 60  # start-up: imports and the first pool worker
    while not out.exists() and app.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    timed_out = False
    try:
        app.wait(timeout=20)  # the window has been asked to close: the exit itself
    except subprocess.TimeoutExpired:
        timed_out = True
        app.kill()
        app.wait(timeout=30)
    stderr = (tmp_path / "stderr.txt").read_text(errors="replace")
    assert not timed_out, "the app did not exit within 20 s\n" + stderr
    assert out.exists(), stderr
    assert "failed:" not in stderr, stderr  # every job ran until it was stopped
    pids = json.loads(out.read_text())
    started = [pids["program"], pids["grandchild"], *pids["workers"]]
    deadline = time.monotonic() + 10
    while any(alive(pid) for pid in started) and time.monotonic() < deadline:
        time.sleep(0.05)
    lingering = [pid for pid in [pids["app"], *started] if alive(pid)]
    for pid in lingering:
        kill(pid)
    assert lingering == [], stderr
