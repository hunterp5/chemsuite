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
"""Session files (.cspx): written and read back whole, refused when not one, and ids renamed
everywhere a session names them."""

from __future__ import annotations

import json
import sqlite3

import pytest

from chemsuite.core.assets import AssetStore
from chemsuite.macromolecule.structure import store as structure_store
from chemsuite.macromolecule.workspace import file
from chemsuite.macromolecule.workspace.session import DERIVED_FROM
from chemsuite.macromolecule.workspace.workspace import MacromoleculeWorkspace

from .test_session import _add


@pytest.fixture
def workspace() -> MacromoleculeWorkspace:
    store = AssetStore()
    _add(store)
    _add(store, "1CRN", "1crn.cif.gz")
    found = MacromoleculeWorkspace(store, deferred=False)
    session = found.session()
    session.select_text("chain E and resi 10-20")
    session.save_selection("site")
    session.display("1CRN").visible = False
    session.measurements.append((("1ATP", 1), ("1ATP", 2)))
    session.labels[("1CRN", 1, "atom")] = None
    found.set_view({"background": "black"})
    return found


def test_a_session_comes_back_whole(workspace, tmp_path) -> None:
    path = tmp_path / "kinase.cspx"
    assets = list(workspace.store.assets.values())
    state = workspace.capture()
    file.save(path, assets, state)
    found = file.read(path)
    assert [a.id for a in found.assets] == ["1ATP", "1CRN"]
    back = structure_store.as_structure(found.assets[0].payload)
    assert back.n_atoms == structure_store.as_structure(assets[0].payload).n_atoms
    assert found.state == json.loads(json.dumps(state))
    assert not list(tmp_path.glob(".*.tmp"))  # written whole, then swapped in


def test_a_file_that_is_not_a_session_is_refused(tmp_path) -> None:
    junk = tmp_path / "junk.cspx"
    junk.write_bytes(b"not a database at all, just bytes" * 10)
    with pytest.raises(file.SessionFileError):
        file.read(junk)
    other = tmp_path / "other.cspx"
    connection = sqlite3.connect(other)
    connection.executescript("CREATE TABLE meta (key TEXT PRIMARY KEY, value BLOB NOT NULL);")
    connection.execute("INSERT INTO meta VALUES ('format', ?)", (b"chemsuite-document",))
    connection.execute("INSERT INTO meta VALUES ('version', ?)", (b"1",))
    connection.commit()
    connection.close()
    with pytest.raises(file.SessionFileError, match="not a chemsuite session"):
        file.read(other)


def test_a_session_from_a_newer_chemsuite_is_refused(tmp_path, monkeypatch) -> None:
    path = tmp_path / "new.cspx"
    monkeypatch.setattr(file, "VERSION", file.VERSION + 1)
    file.save(path, [], {})
    monkeypatch.undo()
    with pytest.raises(file.SessionFileError, match="newer chemsuite"):
        file.read(path)


def test_ids_that_clash_are_renamed_everywhere(workspace) -> None:
    taken = {"1ATP", "1ATP-2"}

    def unique(wanted: str) -> str:
        found, n = wanted, 2
        while found in taken:
            found, n = f"{wanted}-{n}", n + 1
        return found

    ids = file.new_ids(["1ATP", "1CRN"], unique)
    assert ids == {"1ATP": "1ATP-3"}
    state = workspace.capture()
    state["session"]["focused"] = ["1ATP", 5]
    state["prepared"] = {"1ATP": {"protein": "1ATP", "minimised": {"source": "1ATP"}}}
    state["poses"] = {"1ATP": [{"protein": "1ATP", "number": 1}]}
    out = file.remap(state, ids)
    session = out["session"]
    assert set(session["proteins"]) == {"1ATP-3", "1CRN"}
    assert session["named"]["site"].keys() == {"1ATP-3"}
    assert session["focused"] == ["1ATP-3", 5]
    assert session["measurements"] == [[["1ATP-3", 1], ["1ATP-3", 2]]]
    assert session["labels"] == [["1CRN", 1, "atom"]]
    assert out["prepared"] == {"1ATP-3": {"protein": "1ATP-3", "minimised": {"source": "1ATP-3"}}}
    assert out["poses"] == {"1ATP-3": [{"protein": "1ATP-3", "number": 1}]}
    assets = list(workspace.store.assets.values())
    copy = file.remap_assets([assets[1]], {})
    assert copy[0].id == "1CRN"
    derived = file.remap_assets(
        [type(assets[1])("1CRN", "protein", "1CRN", b"", {DERIVED_FROM: "1ATP"})], ids
    )
    assert derived[0].meta[DERIVED_FROM] == "1ATP-3"
