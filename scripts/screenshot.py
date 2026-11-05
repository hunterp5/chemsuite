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

"""Render chemsuite windows offscreen to PNG files, for checking layout without a display.

    .venv/bin/python scripts/screenshot.py OUT_DIR [FILE]

With FILE, that file is opened (any supported format); otherwise a small sample table is used.
Writes ``main.png`` plus one PNG per top-level menu.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

SAMPLE = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    "broken": "C1CC",
    "blank": "",
}


def sample_table():
    from chemsuite.core.table import Table
    from chemsuite.core.table.columns import StructureColumn

    table = Table()
    names = list(SAMPLE)
    table.add_rows(
        {
            "Structure": StructureColumn.from_texts(list(SAMPLE.values())),
            "ID": [f"CPD-{i + 1}" for i in range(len(names))],
            "Name": names,
            "pIC50": [f"{5 + i * 0.37:.2f}" for i in range(len(names))],
            "Notes": ["", "a long note " * 6, "", "two\nlines", "", "", ""],
        },
        undoable=False,
    )
    return table


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    if argv[1] == "--check":  # can Qt list its own plugin directory here? (sandboxes break it)
        from PySide6.QtCore import QDir, QLibraryInfo

        plugins = QLibraryInfo.path(QLibraryInfo.LibraryPath.PluginsPath) + "/platforms"
        print("python sees:", sorted(os.listdir(plugins)))
        print("qt sees:    ", QDir(plugins).entryList())
        return 0
    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)

    from PySide6.QtCore import QCoreApplication, QEventLoop
    from PySide6.QtWidgets import QApplication

    from chemsuite.app.jobs import JobRunner
    from chemsuite.app.main_window import MainWindow, Workspace
    from chemsuite.app.table_view import StructureImages
    from chemsuite.core.io.load import load_table, parse_pending_structures

    app = QApplication([argv[0]])
    runner = JobRunner(app)
    if len(argv) > 2:
        source = Path(argv[2])
        table = load_table(source)
    else:
        source = None
        table = sample_table()
    parse_pending_structures(table)  # synchronously, so drawings are ready for the shot

    window = MainWindow(runner, StructureImages(), Workspace(table, None))
    window.resize(1280, 800)
    window.show()
    window.view.selectRow(1)

    def settle() -> None:
        for _ in range(5):
            QCoreApplication.processEvents(QEventLoop.ProcessEventsFlag.AllEvents, 50)

    settle()
    window.grab().save(str(out / "main.png"))
    for action in window.menuBar().actions():
        menu = action.menu()
        if menu is None:
            continue
        menu.popup(window.mapToGlobal(window.rect().topLeft()))
        settle()
        name = action.text().replace("&", "").lower()
        menu.grab().save(str(out / f"menu-{name}.png"))
        menu.hide()
    table.set_sort([("pIC50", False)])
    settle()
    window.grab().save(str(out / "main-sorted.png"))
    column = window.workspace.model.column_index("ID")
    menu = window.columns.build_menu("ID", [], column)
    menu.popup(window.mapToGlobal(window.rect().center()))
    settle()
    menu.grab().save(str(out / "menu-header.png"))
    menu.hide()
    print("\n".join(str(p) for p in sorted(out.glob("*.png"))))
    runner.shutdown()
    window.workspace.table.history.mark_clean()  # skip the save prompt
    window.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
