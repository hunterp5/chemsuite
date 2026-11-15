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
"""Docked poses: their saved form, and the workspace keeping them."""

from __future__ import annotations

from chemsuite.core.assets import AssetStore
from chemsuite.protein.workspace.poses import DOCKED, RESCORED, DockedPoses, Pose, runs_from_json
from chemsuite.protein.workspace.workspace import VERSION, ProteinWorkspace

from .test_session import _add

MOLBLOCK = "ethanol\n  made\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
RUN = DockedPoses(
    "1ATP",
    (
        Pose("ethanol", 1, -8.5, 0.0, 0.0, MOLBLOCK),
        Pose("ethanol", 2, -7.25, 1.2, 2.4, MOLBLOCK),
        Pose("phenol", 1, -9.1, 0.0, 0.0, MOLBLOCK),
    ),
    failed=1,
    skipped=2,
)


def test_a_run_survives_its_saved_form() -> None:
    saved = RUN.to_json()
    assert DockedPoses.from_json(saved) == RUN
    assert RUN.ligands == ("ethanol", "phenol") and RUN.best == -9.1
    assert RUN.kind == DOCKED and saved["poses"][1][2] == -7.25


def test_a_rescore_keeps_what_kind_of_run_it_was_and_how_it_was_scored() -> None:
    rescored = DockedPoses("1ATP", RUN.poses, RESCORED, scoring="vinardo", number=3)
    assert DockedPoses.from_json(rescored.to_json()) == rescored
    assert rescored.label == "Run 3 · Rescored (vinardo)"


def test_damaged_saved_poses_are_ignored_rather_than_half_read() -> None:
    assert DockedPoses.from_json({"protein": "1ATP"}) is None  # no poses
    assert DockedPoses.from_json({"poses": []}) is None  # no protein
    assert DockedPoses.from_json({"protein": "1ATP", "poses": [["a", "b"]]}) is None
    assert DockedPoses.from_json("1ATP") is None
    bare = DockedPoses.from_json({"protein": "1ATP", "poses": []})
    assert bare is not None and bare.kind == DOCKED and bare.best is None
    assert bare.failed == 0 and bare.skipped == 0
    odd = DockedPoses.from_json({"protein": "1ATP", "kind": "guessed", "poses": []})
    assert odd is not None and odd.kind == DOCKED  # an unknown kind reads as a dock


def test_the_workspace_keeps_poses_and_saves_them_with_the_document() -> None:
    store = AssetStore()
    _add(store)
    workspace = ProteinWorkspace(store, deferred=False)
    assert workspace.docked("1ATP") is None and workspace.runs("1ATP") == ()
    first = workspace.add_docked(RUN)
    assert first.number == 1 and workspace.docked("1ATP") == first and workspace.revision == 1
    saved = workspace.capture()
    assert saved["poses"] == {"1ATP": [first.to_json()]}

    again = ProteinWorkspace(store, saved, deferred=False)
    assert again.docked("1ATP") == first and again.revision == 0


def test_every_run_is_kept_and_any_one_can_be_removed() -> None:
    store = AssetStore()
    _add(store)
    workspace = ProteinWorkspace(store, deferred=False)
    first = workspace.add_docked(RUN)
    second = workspace.add_docked(DockedPoses("1ATP", RUN.poses[:1], RESCORED, scoring="vina"))
    assert [run.number for run in workspace.runs("1ATP")] == [1, 2]
    assert workspace.docked("1ATP") == second and workspace.docked("1ATP", 1) == first
    assert workspace.drop_docked("1ATP", 1) and not workspace.drop_docked("1ATP", 1)
    third = workspace.add_docked(RUN)
    assert third.number == 3  # after the newest, not in the gap: a run keeps its name
    again = ProteinWorkspace(store, workspace.capture(), deferred=False)
    assert again.runs("1ATP") == (second, third)
    assert again.drop_docked("1ATP") and again.runs("1ATP") == ()


def test_a_document_that_kept_one_run_reads_it_as_run_one() -> None:
    """Before runs were kept side by side a protein saved one run, not a list."""
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {}, "view": {}, "poses": {"1ATP": RUN.to_json()}}
    workspace = ProteinWorkspace(store, saved, deferred=False)
    (run,) = workspace.runs("1ATP")
    assert run.number == 1 and run.poses == RUN.poses
    twins = runs_from_json([RUN.to_json(), RUN.to_json(), "damaged"])
    assert [run.number for run in twins] == [1, 2]


def test_poses_without_their_protein_are_not_saved() -> None:
    workspace = ProteinWorkspace(AssetStore(), deferred=False)  # no 1ATP in this store
    workspace.add_docked(RUN)
    assert workspace.capture() is None
    assert workspace.drop_docked("1ATP") and not workspace.drop_docked("1ATP")


def test_an_older_document_without_poses_still_loads() -> None:
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {"named": {}}, "view": {}}
    workspace = ProteinWorkspace(store, saved, deferred=False)
    assert workspace.poses == {} and workspace.capture() == saved


def test_damaged_poses_leave_the_rest_of_the_state_alone() -> None:
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {"named": {}}, "view": {}, "poses": {"1ATP": "no"}}
    workspace = ProteinWorkspace(store, saved, deferred=False)
    assert workspace.poses == {} and workspace.session() is not None
