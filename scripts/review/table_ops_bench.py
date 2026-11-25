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

"""Time table operations on a large document in a real (offscreen) main window.

Each step is timed until the event queue is empty again (repaints, listeners, deferred
work), so the number is how long the window would be frozen.

usage: python scripts/review/table_ops_bench.py DOC.csd [--profile STEP] [--onscreen]
"""

from __future__ import annotations

import argparse
import cProfile
import os
import pstats
import sys
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("doc", type=Path)
    parser.add_argument("--profile", default="", help="step name to profile (or 'all')")
    parser.add_argument("--onscreen", action="store_true")
    args = parser.parse_args()
    if not args.onscreen:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    import numpy as np
    from PySide6.QtCore import QItemSelection, QItemSelectionModel
    from PySide6.QtWidgets import QApplication

    from chemsuite.app.document.main_window import MainWindow
    from chemsuite.app.document.workspace import Workspace
    from chemsuite.app.shell.jobs import JobRunner
    from chemsuite.app.table.cells import StructureImages
    from chemsuite.core.document import load_document

    app = QApplication.instance() or QApplication(sys.argv)
    t = time.perf_counter()
    table, extras = load_document(args.doc)
    print(
        f"load_document {time.perf_counter() - t:.2f}s",
        f"rows={len(table)} cols={len(table.column_names)}",
    )
    runner = JobRunner()
    window = MainWindow(runner, StructureImages(), Workspace(table, args.doc, extras))
    window.resize(1400, 900)
    window.show()
    settle(app)

    def step(name: str, fn) -> None:
        profiler = cProfile.Profile() if args.profile in (name, "all") else None
        t = time.perf_counter()
        if profiler:
            profiler.enable()
        fn()
        t_call = time.perf_counter() - t
        settle(app)
        if profiler:
            profiler.disable()
        total = time.perf_counter() - t
        sel = len(window.workspace.selection)
        print(
            f"{name:28s} call {t_call:6.2f}s  until idle {total:6.2f}s",
            f"rows={len(table)} selected={sel}",
        )
        if profiler:
            pstats.Stats(profiler).sort_stats("cumulative").print_stats(25)

    view = window.view
    model = view.model()
    selection = window.workspace.selection

    def select_alternate() -> None:  # worst case for Qt's ranges: one per row
        sel = QItemSelection()
        last = model.columnCount() - 1
        for r in range(0, model.rowCount(), 2):
            sel.select(model.index(r, 0), model.index(r, last))
        view.selectionModel().select(sel, QItemSelectionModel.SelectionFlag.ClearAndSelect)

    oids = table.oids
    step("select all (view)", view.selectAll)
    step("clear selection", window.clear_selection)
    step(
        "select 50k (store, block)", lambda: selection.set(oids[: len(oids) // 2].tolist(), "bench")
    )
    step("invert selection", lambda: selection.invert("menu"))
    step("select every other (store)", lambda: selection.set(oids[::2].tolist(), "bench"))
    step("invert (alternating)", lambda: selection.invert("menu"))
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest

    def click(row: int, col: int, mods=Qt.KeyboardModifier.NoModifier) -> None:
        rect = view.visualRect(model.index(row, col))
        QTest.mouseClick(view.viewport(), Qt.MouseButton.LeftButton, mods, rect.center())

    def header_click(col: int, mods=Qt.KeyboardModifier.NoModifier) -> None:
        header = view.horizontalHeader()
        x = header.sectionViewportPosition(col) + header.sectionSize(col) // 2
        QTest.mouseClick(header.viewport(), Qt.MouseButton.LeftButton, mods, QPoint(x, 5))

    def scroll() -> None:
        bar = view.verticalScrollBar()
        for _ in range(5):
            bar.setValue(bar.value() + bar.pageStep())
            view.viewport().repaint()
            view.verticalHeader().viewport().repaint()
            view.horizontalHeader().viewport().repaint()

    from PySide6.QtCore import QPoint

    ctrl = Qt.KeyboardModifier.ControlModifier
    step("scroll 5 pages (scattered)", scroll)
    step("ctrl+click a cell (scattered)", lambda: click(3, 2, ctrl))
    step("shift+click a cell (scattered)", lambda: click(5, 2, Qt.KeyboardModifier.ShiftModifier))
    step("store: every other again", lambda: selection.set(oids[::2].tolist(), "bench"))
    step("plain click a cell (scattered)", lambda: click(4, 2))
    step("header click column 2", lambda: header_click(2))
    step("scroll 5 pages (1 column)", scroll)
    step("ctrl+click header column 3", lambda: header_click(3, ctrl))
    step("column names (header menu)", view.selected_column_names)
    step("select all, column names", lambda: (view.selectAll(), view.selected_column_names()))
    step("scroll 5 pages (all)", scroll)
    step("store: every other once more", lambda: selection.set(oids[::2].tolist(), "bench"))
    step("delete selection (50k)", window.delete_rows)
    step("undo delete", table.undo)
    step("select all (view) again", view.selectAll)
    step("delete all", window.delete_rows)
    step("undo delete all", table.undo)
    window.workspace.mark_saved(table.history.checkpoint(), window.workspace.lists.revision)
    window.close()
    runner.shutdown()
    _ = np


def settle(app) -> None:
    """Process events until nothing is left (deferred timers included, briefly)."""
    from PySide6.QtCore import QCoreApplication, QEventLoop

    end = time.perf_counter() + 0.05
    while True:
        app.processEvents(QEventLoop.ProcessEventsFlag.AllEvents)
        QCoreApplication.sendPostedEvents()
        if time.perf_counter() > end:
            break


if __name__ == "__main__":
    main()
