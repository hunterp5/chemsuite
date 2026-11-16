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
"""Session files (``.cspx``): a macromolecule workspace on its own, to open again anywhere.

A session is the SQLite container every chemsuite file uses (:mod:`chemsuite.core.container`):
the ``assets`` table holds the structures, and ``meta`` holds the format, its version and the
workspace's saved state (``workspace``, the JSON of
:meth:`~.workspace.MacromoleculeWorkspace.capture`: looks, named selections, focus, labels,
measurements, the view, reports, poses, reviews and preparation records).

Opening a session beside structures already open may meet ids that are taken: :func:`remap`
renames them everywhere the session names a structure.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from chemsuite import __version__
from chemsuite.core import container
from chemsuite.core.assets import Asset
from chemsuite.core.jobs import CancelToken

from .session import DERIVED_FROM

FORMAT = "chemsuite-macromolecule-session"
VERSION = 1  # the newest session format this version reads
SUFFIX = ".cspx"
WHAT = "chemsuite session"


class SessionFileError(container.ContainerError):
    """Not a session file, a damaged one, or one from a newer chemsuite."""


@dataclass(frozen=True)
class Session:
    """What a session file holds: the structures, and the workspace state (None: none saved)."""

    assets: tuple[Asset, ...]
    state: dict[str, Any] | None


def save(
    path: str | Path,
    assets: Iterable[Asset],
    state: Any,
    cancel: CancelToken | None = None,
) -> None:
    """Write *assets* (the structures) and the workspace *state* to a session file at *path*,
    replacing it in one step."""
    meta = {
        "format": FORMAT.encode(),
        "version": str(VERSION).encode(),
        "app_version": __version__.encode(),
        "workspace": json.dumps(state).encode(),
    }

    def write(connection) -> None:
        connection.executescript(container.META)
        container.write_meta(connection, meta)
        container.write_assets(connection, assets, cancel=cancel)

    container.atomic_save(path, write)


def read(path: str | Path) -> Session:
    """The session file at *path*. Raises :class:`SessionFileError` for a file that is not
    one, is damaged, or was saved by a newer chemsuite."""

    def load(connection) -> Session:
        meta = container.read_meta(connection)
        container.check_version(meta, FORMAT, VERSION, error=SessionFileError, what=WHAT)
        state = json.loads(meta["workspace"]) if "workspace" in meta else None
        assets = container.read_assets(connection, error=SessionFileError)
        return Session(tuple(assets), state if isinstance(state, dict) else None)

    return container.read_container(path, load, error=SessionFileError, what=WHAT)


def is_session(path: str | Path) -> bool:
    return Path(path).suffix.lower() == SUFFIX


# -- ids that clash ------------------------------------------------------------------------------
def new_ids(incoming: Iterable[str], unique: Callable[[str], str]) -> dict[str, str]:
    """Each incoming id that *unique* changes (taken already), old → new. *unique* is the
    store's ``unique_asset_id``; ids given out here are not taken yet, so each is checked
    against them too."""
    out: dict[str, str] = {}
    given: set[str] = set()
    for old in incoming:
        new = unique(old)
        n = 2
        while new in given:
            new = unique(f"{old}-{n}")
            n += 1
        given.add(new)
        if new != old:
            out[old] = new
    return out


def remap_assets(assets: Iterable[Asset], ids: Mapping[str, str]) -> list[Asset]:
    """*assets* with their ids, and the structure each was made from, renamed by *ids*."""
    out = []
    for asset in assets:
        meta = dict(asset.meta)
        source = meta.get(DERIVED_FROM)
        if source in ids:
            meta[DERIVED_FROM] = ids[source]
        out.append(replace(asset, id=ids.get(asset.id, asset.id), meta=meta))
    return out


def remap(state: Mapping[str, Any], ids: Mapping[str, str]) -> dict[str, Any]:
    """A saved workspace *state* with every structure id renamed by *ids*: displays, named
    selections, focus, labels, measurements; reports, poses, reviews and preparation records
    (their keys, the structure each names, and the one a minimised structure came from)."""
    if not ids:
        return dict(state)
    out = dict(state)
    if isinstance(state.get("session"), dict):
        out["session"] = _remap_session(state["session"], ids)
    for key in ("reports", "reviews", "prepared"):
        if key in out:
            out[key] = _keyed(out[key], ids, lambda v: _record(v, ids))
    if "poses" in out:
        out["poses"] = _keyed(out["poses"], ids, lambda runs: _runs(runs, ids))
    return out


def _to(name: Any, ids: Mapping[str, str]) -> Any:
    return ids.get(name, name) if isinstance(name, str) else name


def _first(item: Any, ids: Mapping[str, str]) -> Any:
    """[structure id, …] with the id renamed."""
    return [_to(item[0], ids), *item[1:]] if isinstance(item, list) and item else item


def _keyed(saved: Any, ids: Mapping[str, str], each: Callable[[Any], Any] = lambda v: v) -> Any:
    if not isinstance(saved, dict):
        return saved
    return {_to(k, ids): each(v) for k, v in saved.items()}


def _record(saved: Any, ids: Mapping[str, str]) -> Any:
    if not isinstance(saved, dict):
        return saved
    out = dict(saved)
    if "protein" in out:
        out["protein"] = _to(out["protein"], ids)
    minimised = out.get("minimised")
    if isinstance(minimised, dict) and "source" in minimised:
        out["minimised"] = {**minimised, "source": _to(minimised["source"], ids)}
    return out


def _runs(runs: Any, ids: Mapping[str, str]) -> Any:
    return [_record(r, ids) for r in runs] if isinstance(runs, list) else runs


def _remap_session(session: Mapping[str, Any], ids: Mapping[str, str]) -> dict[str, Any]:
    out = dict(session)
    out["proteins"] = _keyed(session.get("proteins"), ids)
    named = session.get("named")
    if isinstance(named, dict):
        out["named"] = {n: _keyed(by, ids) for n, by in named.items()}
    out["focused"] = _first(session.get("focused"), ids)
    measurements = session.get("measurements")
    if isinstance(measurements, list):
        out["measurements"] = [
            [_first(p, ids) for p in m] if isinstance(m, list) else m for m in measurements
        ]
    labels = session.get("labels")
    if isinstance(labels, list):
        out["labels"] = [_first(label, ids) for label in labels]
    return out
