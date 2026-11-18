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
"""What large conformer ensembles cost in the app (ROADMAP §7.0f; ``pytest -m bench -s``).

Uses the documents ``scripts/review/conformer_scale.py data`` makes (1,000 molecules with
CONFORGE conformers, local under ``samples/``); skipped without them. The targets are in
``docs/conformers/scale.md``.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QElapsedTimer, QTimer
from PySide6.QtWidgets import QApplication

from tests.gui import gui

DATA = Path(__file__).resolve().parents[3] / "samples" / "conformers" / "scale"
COUNTS = [10, 100, 300]
THREE_D = "SMILES 3D"
pytestmark = [pytest.mark.bench, gui]


def report(count: int, label: str, seconds: float) -> None:
    print(f"\n  {count:>3} conformers  {label:<44} {seconds * 1e3:8.1f} ms")


@pytest.fixture(params=COUNTS, ids=[f"{c}x1k" for c in COUNTS])
def scaled(request, qtbot, window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.document import load_document
    from chemsuite.core.table.columns import StructureColumn

    path = DATA / f"scale_{request.param}.csd"
    if not path.exists():
        pytest.skip(f"run scripts/review/conformer_scale.py data first ({path.name})")
    table, _extras = load_document(path)
    StructureColumn._cache.clear()
    window._install(Workspace(table))
    window.resize(1400, 900)
    window.show()
    qtbot.waitExposed(window)
    QApplication.processEvents()
    return request.param, window


def _cold(window) -> None:
    from chemsuite.core.table.columns import StructureColumn

    StructureColumn._cache.clear()
    window._images.clear()


def _timed(fn) -> float:
    start = time.perf_counter()
    fn()
    return time.perf_counter() - start


def test_painting_and_scrolling_the_3d_column(scaled) -> None:
    count, window = scaled
    view = window.view
    model = window.workspace.model
    column = next(i for i in range(model.columnCount()) if model.column_name(i) == THREE_D)
    view.horizontalScrollBar().setValue(view.columnViewportPosition(column))
    _cold(window)
    report(count, "first screen of the 3D column (cold)", _timed(view.viewport().repaint))
    bar = view.verticalScrollBar()
    pages = []
    for _ in range(5):
        bar.setValue(bar.value() + bar.pageStep())
        pages.append(_timed(view.viewport().repaint))
    report(count, "a page down (pictures not cached, mean of 5)", float(np.mean(pages)))
    bar.setValue(0)
    report(count, "back to the top (cached)", _timed(view.viewport().repaint))


def test_sorting_another_column(scaled) -> None:
    count, window = scaled
    table = window.workspace.table

    def sort() -> None:
        table.set_sort([("ID", False)])
        QApplication.processEvents()

    report(count, "sort by ID and repaint", _timed(sort))


def test_applying_a_new_3d_column_and_undoing_it(scaled) -> None:
    from chemsuite.app.document.tools import ToolContext

    count, window = scaled
    table = window.workspace.table
    rows = np.arange(len(table))
    copy = table.column(THREE_D).take(rows)
    ctx = ToolContext(window)

    def apply() -> None:
        ctx.apply_results({"Copy 3D": (table.oids.copy(), copy)}, label="Bench")
        QApplication.processEvents()

    report(count, "apply a 1k-row 3D column", _timed(apply))
    report(count, "undo it", _timed(lambda: (table.undo(), QApplication.processEvents())))
    report(count, "redo it", _timed(lambda: (table.redo(), QApplication.processEvents())))


def test_the_browser_steps_and_moves(scaled, monkeypatch) -> None:
    from chemsuite.app.document.tools import ToolContext
    from chemsuite.app.viewer.panel import StructureViewer
    from chemsuite.tools.viewer.tool import open_viewer

    monkeypatch.setattr(StructureViewer, "_failed", lambda self, *_: None)  # offscreen: no GPU
    count, window = scaled
    ctx = ToolContext(window)
    ctx.set_current_oid(int(window.workspace.table.oids[0]))
    panel = open_viewer(ctx, "3d")
    panel.column.setCurrentText(THREE_D)
    panel.refresh()
    steps = [_timed(lambda: panel.step(1)) for _ in range(10)]
    report(count, "step a conformer (mean of 10)", float(np.mean(steps)))
    moves = []
    for _ in range(10):
        panel.go_by(1)
        moves.append(_timed(panel.refresh))
    report(count, "move to the next row (mean of 10)", float(np.mean(moves)))
    panel.overlay_action.setChecked(True)
    report(count, "overlay every conformer", _timed(panel.refresh))


def test_the_ui_stays_responsive_while_saving(scaled, tmp_path) -> None:
    count, window = scaled
    clock = QElapsedTimer()
    gaps: list[int] = []

    def tick() -> None:
        if clock.isValid():
            gaps.append(clock.restart())
        else:
            clock.start()

    timer = QTimer()
    timer.setInterval(10)
    timer.timeout.connect(tick)
    timer.start()
    took = _timed(lambda: window._save_to(tmp_path / "copy.csd"))
    timer.stop()
    report(count, "save (whole job)", took)
    report(count, "longest UI-thread stall while saving", max(gaps, default=0) / 1e3)
