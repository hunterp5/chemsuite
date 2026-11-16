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
"""Preparations being reviewed: their saved form, and the workspace keeping them."""

from __future__ import annotations

from chemsuite.core.assets import AssetStore
from chemsuite.protein.workspace.reviews import SETTINGS, Review
from chemsuite.protein.workspace.workspace import VERSION, ProteinWorkspace

from .test_session import _add

REVIEW = Review(
    "1ATP",
    {**SETTINGS, "ph": 7.4, "waters": "near ligands"},
    ("missing:12:CZ,NH1",),
    {"amide:24:CB/OD1": "flipped"},
)


def test_a_review_survives_its_saved_form() -> None:
    assert Review.from_json(REVIEW.to_json()) == REVIEW and REVIEW.changed == 2


def test_damaged_settings_fall_back_to_the_defaults() -> None:
    saved = {"protein": "1ATP", "settings": {"ph": "seven", "repair": 1, "waters": "none"}}
    review = Review.from_json(saved)
    assert review is not None
    assert review.settings["ph"] == 7.0 and review.settings["repair"] is True
    assert review.settings["waters"] == "none"  # what was valid is kept


def test_damaged_saved_reviews_are_ignored_rather_than_half_read() -> None:
    assert Review.from_json({"skipped": []}) is None  # no protein
    assert Review.from_json({"protein": "1ATP", "skipped": "x"}) is None
    assert Review.from_json("1ATP") is None
    odd = Review.from_json({"protein": "1ATP", "skipped": ["a", 3], "overrides": {"b": 4}})
    assert odd is not None and odd.skipped == ("a",) and odd.overrides == {}


def test_a_neutral_n_terminus_chosen_under_its_old_name_reads_as_the_new_one() -> None:
    from chemsuite.protein.prepare import sites

    saved = {"protein": "1ATP", "overrides": {"nterm:1:N": "NH2 120°", "lys:7:NZ": "LYN 15°"}}
    review = Review.from_json(saved)
    assert review is not None
    assert review.overrides == {"nterm:1:N": f"{sites.NEUTRAL_N} 120°", "lys:7:NZ": "LYN 15°"}


def test_the_workspace_keeps_a_review_until_it_is_applied() -> None:
    store = AssetStore()
    _add(store)
    workspace = ProteinWorkspace(store, deferred=False)
    workspace.set_review(REVIEW)
    assert workspace.review("1ATP") == REVIEW and workspace.revision == 1
    workspace.set_review(REVIEW)  # the same again is not a change
    assert workspace.revision == 1
    saved = workspace.capture()
    assert saved["reviews"] == {"1ATP": REVIEW.to_json()}
    again = ProteinWorkspace(store, saved, deferred=False)
    assert again.review("1ATP") == REVIEW and again.revision == 0
    assert again.drop_review("1ATP") and not again.drop_review("1ATP")
    assert "reviews" not in (again.capture() or {})


def test_an_older_document_without_reviews_still_loads() -> None:
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {"named": {}}, "view": {}}
    workspace = ProteinWorkspace(store, saved, deferred=False)
    assert workspace.reviews == {} and workspace.capture() == saved
