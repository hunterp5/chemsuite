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

"""The JSON-lines worker client, against a worker written in plain Python (so it runs on every
system, Windows included)."""

from __future__ import annotations

import sys
import threading
import time

import pytest

from chemsuite.core.jobs import Cancelled, CancelToken, lines
from chemsuite.core.jobs.lines import LineWorker, WarmWorkers, WorkerError

#: a worker speaking the protocol: ``echo``, ``events`` (n events, *gap* s apart), ``slow``
#: (runs until stopped), ``deaf`` (ignores a stop), ``noisy`` (prints to stdout), ``crash``
FAKE = r"""
import json, os, sys, threading, time
out = os.fdopen(os.dup(1), "w", encoding="utf-8", newline="\n")
os.dup2(2, 1)
sys.stdout = sys.stderr
stop = threading.Event()
work = []
lock = threading.Lock()

def send(payload):
    with lock:
        out.write(json.dumps(payload) + "\n")
        out.flush()

def read():
    for line in sys.stdin:
        request = json.loads(line)
        if request["op"] == "stop":
            stop.set()
        else:
            work.append(request)
    stop.set()
    work.append(None)

threading.Thread(target=read, daemon=True).start()
send({"id": 0, "ok": True, "ready": True})
while True:
    while not work:
        time.sleep(0.01)
    request = work.pop(0)
    if request is None:
        sys.exit(0)
    ident, op = request["id"], request["op"]
    stop.clear()
    if op == "quit":
        send({"id": ident, "ok": True})
        sys.exit(0)
    elif op == "echo":
        send({"id": ident, "ok": True, "body": request.get("body")})
    elif op == "events":
        for i in range(request["n"]):
            time.sleep(request["gap"])
            send({"id": ident, "event": {"i": i}})
        send({"id": ident, "ok": True})
    elif op == "slow":
        steps = 0
        while not stop.is_set():
            steps += 1
            time.sleep(0.02)
        send({"id": ident, "ok": True, "steps": steps})
    elif op == "deaf":
        time.sleep(60)
    elif op == "noisy":
        print("a library talking")
        os.write(1, b"straight to fd 1\n")
        send({"id": ident, "ok": True})
    elif op == "crash":
        sys.stderr.write("something broke\n")
        sys.stderr.flush()
        os._exit(3)
    else:
        send({"id": ident, "ok": False, "error": "unknown " + op})
"""


@pytest.fixture
def worker():
    found = LineWorker([sys.executable, "-I", "-c", FAKE], name="The fake", start_timeout=30)
    yield found
    found.close()


def test_a_request_comes_back_by_its_id(worker) -> None:
    assert worker.request("echo", {"body": "é"}, timeout=10)["body"] == "é"
    assert worker.alive
    assert worker.request("other", timeout=10) == {"id": 2, "ok": False, "error": "unknown other"}


def test_events_arrive_in_order_on_the_callers_thread(worker) -> None:
    seen = []
    caller = threading.get_ident()
    worker.request("events", {"n": 4, "gap": 0.01}, timeout=10,
                   events=lambda e: seen.append((e["i"], threading.get_ident())))  # fmt: skip
    assert seen == [(i, caller) for i in range(4)]


def test_the_timeout_is_for_silence_so_events_keep_a_request_alive(worker) -> None:
    worker.request("events", {"n": 5, "gap": 0.3}, timeout=1.0)  # 1.5 s in all
    with pytest.raises(WorkerError, match="did not answer within 1 s"):
        worker.request("events", {"n": 1, "gap": 3.0}, timeout=1.0)


def test_a_stoppable_request_is_stopped_and_the_worker_stays_warm(worker) -> None:
    cancel = CancelToken()
    threading.Timer(0.3, cancel.cancel).start()
    with pytest.raises(Cancelled) as raised:
        worker.request("slow", timeout=10, cancel=cancel, stoppable=True)
    assert raised.value.partial["steps"] > 0  # the reply to the stopped request
    assert worker.alive
    assert worker.request("echo", {"body": 1}, timeout=10)["body"] == 1


def test_a_request_that_cannot_stop_is_killed_at_once(worker) -> None:
    cancel = CancelToken()
    threading.Timer(0.3, cancel.cancel).start()
    began = time.monotonic()
    with pytest.raises(Cancelled) as raised:
        worker.request("deaf", timeout=10, cancel=cancel)
    assert time.monotonic() - began < 3 and raised.value.partial is None
    assert not worker.alive
    # started again at once: the killed worker's end of stream does not reach the new one
    for body in range(5):
        assert worker.request("echo", {"body": body}, timeout=10)["body"] == body


def test_a_worker_deaf_to_a_stop_is_killed_after_the_grace(worker, monkeypatch) -> None:
    monkeypatch.setattr(lines, "STOP_GRACE_S", 0.5)
    cancel = CancelToken()
    threading.Timer(0.2, cancel.cancel).start()
    with pytest.raises(Cancelled):
        worker.request("deaf", timeout=10, cancel=cancel, stoppable=True)
    assert not worker.alive


def test_what_a_library_prints_does_not_corrupt_the_replies(worker) -> None:
    assert worker.request("noisy", timeout=10)["ok"]
    assert worker.request("echo", {"body": 2}, timeout=10)["body"] == 2


def test_a_worker_that_dies_is_reported_with_its_exit_code_and_output(worker) -> None:
    with pytest.raises(WorkerError) as raised:
        worker.request("crash", timeout=10)
    assert "exit code 3" in str(raised.value) and "something broke" in str(raised.value)
    assert not worker.alive
    assert worker.request("echo", {"body": 3}, timeout=10)["body"] == 3  # started again


def test_closing_ends_the_worker(worker) -> None:
    worker.request("echo", timeout=10)
    process = worker._process
    worker.close()
    assert process is not None and process.poll() == 0 and not worker.alive


def test_a_program_that_is_not_there_is_reported() -> None:
    with pytest.raises(WorkerError, match="Could not start"):
        LineWorker(["no-such-program-anywhere"]).request("echo", timeout=1)


def test_one_warm_worker_at_a_time_when_asked() -> None:
    made = []

    class Quiet(LineWorker):
        closed = False

        def close(self) -> None:
            self.closed = True

    def make(key):
        made.append(Quiet([key]))
        return made[-1]

    single = WarmWorkers(make, single=True)
    a = single.get("a")
    assert single.get("a") is a
    b = single.get("b")
    assert a.closed and not b.closed and len(made) == 2
    several = WarmWorkers(make)
    c, d = several.get("c"), several.get("d")
    assert not c.closed and several.get("c") is c
    several.close_all()
    assert c.closed and d.closed


def test_a_conda_interpreter_finds_its_own_libraries(tmp_path, monkeypatch) -> None:
    root = tmp_path / "env"
    (root / "conda-meta").mkdir(parents=True)
    (root / "Library" / "bin").mkdir(parents=True)
    python = root / "python.exe"
    python.write_text("", encoding="utf-8")
    monkeypatch.setenv("PYTHONPATH", "chemsuite's own")
    env = lines.child_environment(str(python))
    assert "PYTHONPATH" not in env and env["PYTHONNOUSERSITE"] == "1"
    first = env["PATH"].split(lines.os.pathsep)[:2]
    assert first == [str(root), str(root / "Library" / "bin")]
    assert env["CONDA_PREFIX"] == str(root)
    assert lines.is_app_python(sys.executable) and not lines.is_app_python(str(python))
