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

"""Preparation reports: their saved form, and the workspace keeping them."""

from __future__ import annotations

from chemsuite.core.assets import AssetStore
from chemsuite.protein.workspace.reports import PrepareReport, Row
from chemsuite.protein.workspace.workspace import VERSION, ProteinWorkspace

from .test_session import _add

REPORT = PrepareReport(
    "1ATP",
    (Row("A/HIS 87", "His: HIE", False, 12), Row("A/ASN 99", "Amide: flipped", True, 24)),
    ph=7.4,
    net_charge=-3,
)


def test_a_report_survives_its_saved_form() -> None:
    saved = REPORT.to_json()
    assert PrepareReport.from_json(saved) == REPORT and REPORT.warnings == 1
    assert saved["rows"][1] == ["A/ASN 99", "Amide: flipped", True, 24]


def test_damaged_saved_reports_are_ignored_rather_than_half_read() -> None:
    assert PrepareReport.from_json({"protein": "1ATP"}) is None  # no rows
    assert PrepareReport.from_json({"rows": []}) is None  # no protein
    assert PrepareReport.from_json({"protein": "1ATP", "rows": [["a", "b"]]}) is None
    assert PrepareReport.from_json("1ATP") is None
    bare = PrepareReport.from_json({"protein": "1ATP", "rows": []})
    assert bare is not None and bare.ph == 7.0 and bare.net_charge == 0


def test_the_workspace_keeps_a_report_and_saves_it_with_the_document() -> None:
    store = AssetStore()
    _add(store)
    workspace = ProteinWorkspace(store, deferred=False)
    assert workspace.report("1ATP") is None
    workspace.set_report(REPORT)
    assert workspace.report("1ATP") == REPORT and workspace.revision == 1
    saved = workspace.capture()
    assert saved["reports"] == {"1ATP": REPORT.to_json()}

    again = ProteinWorkspace(store, saved, deferred=False)
    assert again.report("1ATP") == REPORT and again.revision == 0


def test_a_report_without_its_protein_is_not_saved() -> None:
    workspace = ProteinWorkspace(AssetStore(), deferred=False)  # no 1ATP in this store
    workspace.set_report(REPORT)
    assert workspace.capture() is None
    assert workspace.drop_report("1ATP") and not workspace.drop_report("1ATP")


def test_an_older_document_without_reports_still_loads() -> None:
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {"named": {}}, "view": {}}
    workspace = ProteinWorkspace(store, saved, deferred=False)
    assert workspace.reports == {} and workspace.capture() == saved


def test_a_damaged_report_leaves_the_rest_of_the_state_alone() -> None:
    store = AssetStore()
    _add(store)
    saved = {
        "version": VERSION,
        "session": {"named": {"x": {"1ATP": [1]}}},
        "view": {"background": "white"},
        "reports": {"1ATP": {"protein": "1ATP"}, "2XYZ": REPORT.to_json()},
    }
    workspace = ProteinWorkspace(store, saved, deferred=False)
    assert list(workspace.reports) == ["2XYZ"]  # the damaged one dropped, the sound one kept
    assert workspace.view == {"background": "white"}
    assert "reports" not in workspace.capture()  # 2XYZ is not a protein of this document
