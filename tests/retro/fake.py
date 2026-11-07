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

"""A stand-in for an AiZynthFinder interpreter, so the worker protocol can be tested.

:func:`fake_interpreter` writes an executable that chemsuite starts exactly as it would start
the real thing. It honours ``-c``: the script chemsuite passes is the real
:data:`chemsuite.retro.worker.SCRIPT`, executed against a stub ``aizynthfinder`` package. So
the protocol, the reply stream, the per-target error handling and the settings the worker
applies are all exercised, with recorded output in place of a search.

What a target asks the stub to do:

====================  ==========================================================
target                behaviour
====================  ==========================================================
``FAIL``              ``prepare_tree`` raises ``ValueError`` (one target fails)
``CRASH``             the process exits at once (the worker dies mid-job)
``SLOW``              the search takes 30 s (so a test can cancel it)
``NOISY``             writes to stdout and stderr before answering
anything else         the recorded paracetamol routes
====================  ==========================================================

Every search appends what the worker asked for to ``FAKE_AZF_LOG``, one JSON object per line,
so a test can check that the dialog's options reached AiZynthFinder's own settings.
"""

from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

from chemsuite.retro import worker as worker_module

#: the stub package the worker's ``import aizynthfinder.aizynthfinder`` finds
_STUB = '''
"""A stub of AiZynthFinder's public API, answering with recorded output."""
import json
import os
import sys
import time

RECORDED = json.loads(os.environ["FAKE_AZF_RECORDED"])


class _Collection:
    def __init__(self, items):
        self.items = list(items)
        self.selection = []

    def select(self, value, append=False):
        keys = [value] if isinstance(value, str) else list(value)
        for key in keys:
            if key not in self.items:
                raise KeyError(key)
        self.selection = keys

    def select_all(self):
        self.selection = list(self.items)

    def reset_cache(self):
        pass


class _Search:
    algorithm = "mcts"
    max_transforms = 6
    iteration_limit = 100
    time_limit = 120
    return_first = False
    exclude_target_from_stock = True

    def __init__(self):
        self.algorithm_config = {"C": 1.4, "default_prior": 0.5}


class _Post:
    min_routes = 5
    max_routes = 25


class _Config:
    def __init__(self):
        self.search = _Search()
        self.post_processing = _Post()


class _Routes:
    def __init__(self, trees):
        self._trees = trees

    def compute_scores(self, *scorers):
        pass

    def dict_with_extra(self, include_metadata=False, include_scores=False):
        return self._trees

    def cluster(self, n_clusters=0):
        return [0, 1, 0][: len(self._trees)]


class _Scorers:
    def objects(self):
        return ["state score"]


class AiZynthFinder:
    def __init__(self, configfile=None, configdict=None):
        if configfile and not os.path.isfile(configfile):
            raise FileNotFoundError(configfile)
        self.config = _Config()
        self.expansion_policy = _Collection(["uspto", "ringbreaker"])
        self.filter_policy = _Collection(["uspto"])
        self.stock = _Collection(["zinc", "bench"])
        self.scorers = _Scorers()
        self.target_smiles = ""
        self.routes = _Routes([])

    def prepare_tree(self):
        if self.target_smiles == "FAIL":
            raise ValueError("Target molecule unsanitizable")
        if self.target_smiles == "CRASH":
            os._exit(3)
        self.routes = _Routes([])

    def tree_search(self, show_progress=False):
        if self.target_smiles == "SLOW":
            time.sleep(30)
        if self.target_smiles == "NOISY":
            sys.stdout.write("a library printed this\\n")
            sys.stdout.flush()
            sys.stderr.write("loading models\\n")
        self._log()
        return 1.0

    def build_routes(self, selection=None, scorer=None):
        self.routes = _Routes(RECORDED["trees"])

    def extract_statistics(self):
        return dict(RECORDED["statistics"], target=self.target_smiles)

    def stock_info(self):
        return dict(RECORDED["stock_info"])

    def _log(self):
        path = os.environ.get("FAKE_AZF_LOG")
        if not path:
            return
        entry = {
            "target": self.target_smiles,
            "algorithm": self.config.search.algorithm,
            "time_limit": self.config.search.time_limit,
            "iteration_limit": self.config.search.iteration_limit,
            "max_transforms": self.config.search.max_transforms,
            "return_first": self.config.search.return_first,
            "exclude_target_from_stock": self.config.search.exclude_target_from_stock,
            "algorithm_config": self.config.search.algorithm_config,
            "max_routes": self.config.post_processing.max_routes,
            "policies": self.expansion_policy.selection,
            "filters": self.filter_policy.selection,
            "stocks": self.stock.selection,
        }
        with open(path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry) + "\\n")
'''

#: run by the shim: puts the stub on the path, then runs the script chemsuite passed to ``-c``
_RUNNER = """
import importlib.metadata as metadata
import sys

stub, args = sys.argv[1], sys.argv[2:]
sys.path.insert(0, stub)
_version = metadata.version
metadata.version = lambda name: "4.4.1" if name == "aizynthfinder" else _version(name)
script = args[args.index("-c") + 1]
sys.argv = ["aizynth-worker", *args[args.index("-c") + 2 :]]
exec(compile(script, "<worker>", "exec"), {"__name__": "__main__", "__file__": "<worker>"})
"""

CONFIG_YML = """\
expansion:
  uspto:
    - uspto_model.onnx
    - uspto_templates.csv.gz
  ringbreaker:
    - uspto_ringbreaker_model.onnx
    - uspto_ringbreaker_templates.csv.gz
filter:
  uspto: uspto_filter_model.onnx
stock:
  zinc: zinc_stock.hdf5
  bench:
    type: inchiset
    path: bench.csv
"""


def fake_interpreter(folder: Path) -> Path:
    """Write the fake interpreter under *folder* and return its path.

    Only POSIX: the shim is a shell script. The tests that use it skip on Windows.
    """
    folder = Path(folder)
    stub = folder / "stub"
    (stub / "aizynthfinder").mkdir(parents=True, exist_ok=True)
    (stub / "aizynthfinder" / "__init__.py").write_text("", encoding="utf-8")
    (stub / "aizynthfinder" / "aizynthfinder.py").write_text(_STUB, encoding="utf-8")
    runner = folder / "runner.py"
    runner.write_text(_RUNNER, encoding="utf-8")
    shim = folder / "fake-python"
    shim.write_text(
        f'#!/bin/sh\nexec "{sys.executable}" "{runner}" "{stub}" "$@"\n', encoding="utf-8"
    )
    shim.chmod(shim.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return shim


def write_config(folder: Path, body: str = CONFIG_YML) -> Path:
    """Write a ``config.yml`` under *folder* and return its path."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "config.yml"
    path.write_text(body, encoding="utf-8")
    return path


def recorded_environment(monkeypatch, trees, statistics, stock_info, log: Path | None = None):
    """Put the recorded output (and the request log) where the stub will find them."""
    payload = {"trees": trees, "statistics": statistics, "stock_info": stock_info}
    monkeypatch.setenv("FAKE_AZF_RECORDED", json.dumps(payload))
    if log is not None:
        monkeypatch.setenv("FAKE_AZF_LOG", str(log))
    else:
        monkeypatch.delenv("FAKE_AZF_LOG", raising=False)


def requests_logged(log: Path) -> list[dict]:
    """What the stub recorded, one entry per search."""
    if not Path(log).is_file():
        return []
    return [json.loads(line) for line in Path(log).read_text(encoding="utf-8").splitlines() if line]


#: the fake runs chemsuite's own worker script, never a copy of it
assert worker_module.SCRIPT
