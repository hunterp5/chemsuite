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
"""The Inspect panel's lines that are chosen, not solved for: locations and crystal waters."""

from __future__ import annotations

import pytest
from PySide6.QtWidgets import QComboBox

from chemsuite.macromolecule.prepare import pipeline, review
from chemsuite.macromolecule.ui.inspector import STATE_COLUMN, InspectorPanel
from tests.gui import gui

pytestmark = gui


@pytest.fixture
def panel(lysozyme, qtbot) -> InspectorPanel:
    found = review.survey("LYZ", lysozyme, pipeline.Options())
    made = InspectorPanel(found)
    qtbot.addWidget(made)
    return made


def _editor(panel: InspectorPanel, id: str) -> QComboBox:
    item = panel.line(id)
    index = panel.rows.indexFromItem(item, STATE_COLUMN)
    delegate = panel.rows.itemDelegateForColumn(STATE_COLUMN)
    editor = delegate.createEditor(panel.rows.viewport(), None, index)
    assert isinstance(editor, QComboBox)
    return editor


def test_a_location_is_chosen_by_letter_with_its_occupancy(panel) -> None:
    altloc = next(c for c in panel.survey.choices if c.kind == review.ALTLOC)
    editor = _editor(panel, altloc.id)
    texts = [editor.itemText(i) for i in range(editor.count())]
    assert all("occupancy" in text for text in texts) and not any("best" in t for t in texts)
    other = next(s.label for k, s in enumerate(altloc.states) if k != altloc.chosen)
    assert panel.set_state(altloc.id, other)
    assert panel.decisions() == ((), {altloc.id: other})
    assert panel.line(altloc.id).text(STATE_COLUMN) == f"Alternate location: {other}"


def test_a_crystal_water_is_kept_or_dropped(panel) -> None:
    water = next(c for c in panel.survey.choices
                 if c.kind == review.KEEP and c.state.label == review.KEPT)  # fmt: skip
    editor = _editor(panel, water.id)
    assert [editor.itemText(i) for i in range(editor.count())] == [review.KEPT, review.DROPPED]
    assert panel.set_state(water.id, review.DROPPED)
    assert panel.decisions() == ((), {water.id: review.DROPPED})


def test_the_new_lines_stay_out_of_warnings_only(panel) -> None:
    panel.only_warnings.setChecked(True)
    for line in panel.lines():
        choice = line.data(0, 0x0100)  # LINE role
        if isinstance(choice, review.Choice) and not choice.scored:
            assert line.isHidden()


def test_a_gap_box_says_it_is_capped_when_unticked(lysozyme, qtbot) -> None:
    from PySide6.QtCore import Qt

    from chemsuite.macromolecule.prepare import inspect as inspection
    from tests.macromolecule.prepare.test_pipeline import LOOPED, _looped

    found = review.survey("LYZ", _looped(lysozyme), LOOPED)
    made = InspectorPanel(found)
    qtbot.addWidget(made)
    (gap,) = [i for i in found.issues if i.kind == inspection.GAP]
    item = made.line(gap.id)
    assert item.checkState(2) == Qt.CheckState.Checked
    assert item.toolTip(2) == "Ticked: built when applied. Unticked: its ends are capped"
