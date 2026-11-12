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
"""The workbench session without Qt or a table: an asset store of its own drives it."""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pytest

from chemsuite.core.assets import Asset, AssetStore, PutAsset
from chemsuite.protein.look import scene as scenes
from chemsuite.protein.structure import io
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.workspace.events import Event
from chemsuite.protein.workspace.session import RESIDUE_LABEL, WorkbenchSession

DATA = Path(__file__).resolve().parents[2] / "data" / "protein"


@cache
def structure(name: str) -> Structure:
    return io.read(DATA / name)


def _add(store: AssetStore, name: str = "1ATP", file: str = "1atp.cif.gz") -> Asset:
    asset = Asset(name, "protein", name, structure(file), {})
    store.edit(f"Open {name}", [PutAsset(name, None, asset)])
    return asset


@pytest.fixture
def store() -> AssetStore:
    store = AssetStore()
    _add(store)
    return store


def test_events_call_back_until_disconnected() -> None:
    event, heard = Event(), []
    stop = event.connect(heard.append)
    event.emit(1)
    stop()
    event.emit(2)
    assert heard == [1] and len(event) == 0


def test_a_store_alone_drives_the_session_and_its_undo(store) -> None:
    session = WorkbenchSession(store)
    told = []
    session.proteinsChanged.connect(lambda: told.append("proteins"))
    assert [d.id for d in session.displays] == ["1ATP"]
    _add(store, "1CRN", "1crn.cif.gz")
    assert [d.id for d in session.displays] == ["1ATP", "1CRN"] and told == ["proteins"]
    store.undo()
    assert [d.id for d in session.displays] == ["1ATP"]
    session.close()
    store.redo()
    assert [d.id for d in session.displays] == ["1ATP"]  # closed: no longer follows


def test_an_edit_is_one_step_of_the_stores_history_and_keeps_the_look(store) -> None:
    session = WorkbenchSession(store)
    session.select_text("resn ATP")
    session.show_as(scenes.SPACEFILL)
    edited = []
    session.edited.connect(lambda: edited.append(True))
    session.select_text("water")
    removed = session.delete_selected()
    assert removed > 0 and store.history.undo_label is not None and edited == [True]
    shown = session.display("1ATP")
    atp = shown.structure.res_name[shown.structure.residue_of] == "ATP"
    assert (shown.styles[atp] == scenes.SPACEFILL).all()  # the look came across by uid
    store.undo()
    assert session.display("1ATP").structure.n_atoms == structure("1atp.cif.gz").n_atoms


def test_the_revision_counts_what_is_worth_saving_only(store) -> None:
    session = WorkbenchSession(store)
    start = session.revision
    session.select_text("resn ATP")  # selecting is not saved
    session.hovered.emit("", -1)
    assert session.revision == start
    session.colour((1.0, 0.0, 0.0))
    session.label(RESIDUE_LABEL)
    session.save_selection("atp")
    assert session.revision == start + 3


def test_saved_state_comes_back_in_a_new_session(store) -> None:
    session = WorkbenchSession(store)
    session.select_text("resn ATP")
    session.save_selection("atp")
    session.label(RESIDUE_LABEL)
    residue = int(np.flatnonzero(session.display("1ATP").structure.res_name == "ATP")[0])
    session.focus("1ATP", residue)
    state = session.state()
    again = WorkbenchSession(store)
    again.set_state(state)
    assert again.named == session.named and again.labels == session.labels
    assert again.focused == session.focused


def test_slow_parts_wait_for_the_schedule_when_deferred(store) -> None:
    session = WorkbenchSession(store, deferred=True)
    asked = []
    session.schedule = asked.append
    session.select_text("polymer")
    session.show_surface(True)
    _scene = session.scene
    assert asked and asked[-1] == session.scene.displays
