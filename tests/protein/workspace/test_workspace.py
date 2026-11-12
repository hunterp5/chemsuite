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
"""A document's protein workspace: kept as saved until used, versioned, counting changes."""

from __future__ import annotations

from chemsuite.core.assets import AssetStore
from chemsuite.protein.workspace.workspace import VERSION, ProteinWorkspace

from .test_session import _add


def test_saved_state_is_kept_untouched_until_the_session_is_made() -> None:
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {"named": {"x": {"1ATP": [1, 2]}}}, "view": {}}
    workspace = ProteinWorkspace(store, saved, deferred=False)
    assert not workspace.started and workspace.capture() == saved and workspace.revision == 0
    assert workspace.session().named == {"x": {"1ATP": [1, 2]}}
    assert workspace.revision == 0  # putting the state back is not a change


def test_a_newer_state_is_saved_back_as_it_was() -> None:
    newer = {"version": VERSION + 1, "something": "new"}
    workspace = ProteinWorkspace(AssetStore(), newer)
    assert workspace.capture() == newer


def test_nothing_saved_makes_nothing_to_save_and_may_adopt_an_older_viewer() -> None:
    workspace = ProteinWorkspace(AssetStore())
    assert workspace.capture() is None and workspace.fresh
    assert workspace.adopt({"session": {"named": {}}, "background": "black"})
    assert workspace.view == {"background": "black"} and not workspace.fresh
    assert not workspace.adopt({"background": "white"})  # only once, and only when fresh


def test_the_views_background_counts_but_the_camera_does_not() -> None:
    workspace = ProteinWorkspace(AssetStore())
    workspace.set_view({"camera": {"distance": 3}})
    assert workspace.revision == 0
    workspace.set_view({"background": "white"})
    assert workspace.revision == 1
