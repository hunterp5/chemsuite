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

"""The script that runs inside AiZynthFinder's own interpreter, and the protocol it speaks.

The script is passed with ``python -I -c``, so no file is written and the interpreter ignores
``PYTHONPATH`` and the user site directory: it sees only its own environment. It loads the
models once and then answers one request per line, which is what keeps a second target (and
the next job) from paying the model load again.

**Protocol.** One JSON object per line in each direction.

Requests
    ``{"id": 1, "op": "hello"}``
        the AiZynthFinder version. Imports the package but loads no models.
    ``{"id": 2, "op": "search", "target": "...", ...}``
        search one target. The rest of the body is :meth:`~.config.SearchOptions.request`.
    ``{"id": 3, "op": "quit"}``
        stop.

Replies
    ``{"id": 2, "ok": true, "target": "...", "stats": {...}, "trees": [...]}``
    ``{"id": 2, "ok": false, "error": "..."}``
        the request failed; the worker stays up for the next one (a target that cannot be
        sanitized must not end the job).

Replies go to a duplicate of the original stdout, and file descriptor 1 is pointed at stderr
before anything is imported, so that AiZynthFinder's logging, a stray ``print`` in a model
library or a progress bar cannot corrupt the stream.
"""

from __future__ import annotations

from pathlib import Path

#: argv for the worker: ``python -I -c SCRIPT <config.yml>``
#:
#: ``-I`` (isolated) keeps the child's imports to its own environment, whatever chemsuite's
#: own PYTHONPATH says.
ARGS = ("-I", "-c")

SCRIPT = r'''
"""chemsuite's AiZynthFinder worker. One JSON request per line; see chemsuite/retro/worker.py."""
import json
import os
import sys
import traceback

# Replies go here; fd 1 becomes stderr so no library can write into the stream.
_replies = os.fdopen(os.dup(1), "w", encoding="utf-8", newline="\n")
os.dup2(2, 1)
sys.stdout = sys.stderr


def _plain(value):
    """JSON-safe form of an AiZynthFinder field (numpy scalars, sets, arrays, objects)."""
    if isinstance(value, (set, frozenset, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        try:
            return _plain(value.tolist())
        except Exception:
            return str(value)
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_plain(item) for item in value]
    if isinstance(value, (str, bool, int, float)) or value is None:
        return value
    return str(value)


def _reply(payload):
    _replies.write(json.dumps(payload, default=str) + "\n")
    _replies.flush()


class Session:
    """The loaded finder, created on the first search."""

    def __init__(self, config):
        self.config = config
        self.finder = None

    def finder_for(self, request):
        if self.finder is None:
            from aizynthfinder.aizynthfinder import AiZynthFinder

            self.finder = AiZynthFinder(configfile=self.config)
        finder = self.finder
        search = request.get("search") or {}
        for key, value in search.items():
            if key == "algorithm_config":
                finder.config.search.algorithm_config.update(value or {})
            elif hasattr(finder.config.search, key):
                setattr(finder.config.search, key, value)
        for key, value in (request.get("post_processing") or {}).items():
            if hasattr(finder.config.post_processing, key):
                setattr(finder.config.post_processing, key, value)
        policies = list(request.get("policies") or []) or list(finder.expansion_policy.items)
        finder.expansion_policy.select(policies)
        chosen_filters = list(request.get("filters") or [])
        if chosen_filters:
            finder.filter_policy.select(chosen_filters)
        elif finder.filter_policy.items:
            finder.filter_policy.select_all()
        finder.stock.select(list(request.get("stocks") or []) or finder.stock.items)
        return finder

    def search(self, request):
        target = str(request.get("target") or "")
        finder = self.finder_for(request)
        finder.target_smiles = target
        finder.prepare_tree()  # ValueError for an unsanitizable target: reported per target
        finder.tree_search()
        finder.build_routes()
        finder.routes.compute_scores(*finder.scorers.objects())
        stats = dict(finder.extract_statistics())
        if request.get("cluster"):
            try:
                stats["cluster_labels"] = finder.routes.cluster(n_clusters=0)
            except Exception as err:  # clustering needs >2 routes and optional deps
                stats["cluster_error"] = str(err)
        stats["stock_info"] = finder.stock_info()
        trees = finder.routes.dict_with_extra(include_metadata=True, include_scores=True)
        return {"target": target, "stats": _plain(stats), "trees": _plain(trees)}


def _hello():
    import importlib.metadata as metadata

    from aizynthfinder.aizynthfinder import AiZynthFinder  # noqa: F401  (import must work)

    return {"version": metadata.version("aizynthfinder"), "python": sys.version.split()[0]}


def main():
    session = Session(sys.argv[1] if len(sys.argv) > 1 else "")
    _reply({"id": 0, "ok": True, "ready": True})
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            request = json.loads(line)
        except ValueError:
            _reply({"id": None, "ok": False, "error": "could not read the request"})
            continue
        ident = request.get("id")
        op = request.get("op")
        if op == "quit":
            _reply({"id": ident, "ok": True})
            return
        try:
            if op == "hello":
                _reply({"id": ident, "ok": True, **_hello()})
            elif op == "search":
                _reply({"id": ident, "ok": True, **session.search(request)})
            else:
                _reply({"id": ident, "ok": False, "error": "unknown request " + repr(op)})
        except Exception as err:
            _reply(
                {
                    "id": ident,
                    "ok": False,
                    "error": "{}: {}".format(type(err).__name__, err),
                    "traceback": traceback.format_exc(limit=8),
                }
            )


main()
'''


def command(python: str, config: Path | str) -> list[str]:
    """Argv that starts the worker in *python* against *config*."""
    return [str(python), *ARGS, SCRIPT, str(config)]
