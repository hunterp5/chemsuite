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
# Never read or write the user's own settings while rendering.
os.environ.setdefault(
    "CHEMSUITE_CONFIG_DIR",
    os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), ".screenshot-config"),
)

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


def main(argv: list[str]) -> int:  # noqa: C901 - a linear script
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
    from chemsuite.core import logbook
    from chemsuite.core.io.load import load_table, parse_pending_structures

    logbook.setup_logging(to_file=False, to_stderr=False)  # the Log window shows the session
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
    table.set_sort([])
    window.filter_column("pIC50").set_bounds(5.5, 7.0)
    category = window.filter_column("Name")
    category.set_checked({n for n in category.listed_values() if n != "broken"})
    window.filters.apply_now()
    settle()
    window.grab().save(str(out / "filters.png"))
    window.filters.remove_all()
    window.filter_dock.hide()
    window.search.edit.setText("in")
    window.search.apply_now()
    dialog = window.show_find_replace()
    dialog.find_edit.setText("CPD")
    dialog.find_next()
    settle()
    window.grab().save(str(out / "search.png"))
    dialog.grab().save(str(out / "find-replace.png"))
    dialog.close_dialog()
    window.search.edit.clear()
    window.search.apply_now()
    window.workspace.lists.add("Actives", table.oids[:4])
    window.workspace.selection.set(table.oids[1:5].tolist(), "script")
    viewer = window.show_browser()
    viewer.go(1)
    settle()
    window.grab().save(str(out / "browser.png"))
    next(d for d in window._panels if d.widget() is viewer).close()
    lists_menu = window._lists_menu
    lists_menu.popup(window.mapToGlobal(window.rect().topLeft()))
    settle()
    lists_menu.grab().save(str(out / "menu-lists.png"))
    lists_menu.hide()
    window.view.setCurrentIndex(window.workspace.model.index(1, 0))
    editor = window.edit_current_structure()
    settle()
    editor.grab().save(str(out / "structure-editor.png"))
    editor.set_text("c1cccc1C(=O)N")
    settle()
    editor.grab().save(str(out / "structure-editor-error.png"))
    editor.reject()
    split = window.split_column("Notes")
    split.update_preview()
    settle()
    split.grab().save(str(out / "split-column.png"))
    split.reject()
    join = window.join_columns(["ID", "Name"])
    join.update_preview()
    settle()
    join.grab().save(str(out / "join-columns.png"))
    join.reject()
    descriptors = window.open_tool("descriptors")
    settle()
    descriptors.grab().save(str(out / "descriptors.png"))
    descriptors.reject()
    menu = next(a.menu() for a in window.menuBar().actions() if a.text() == "&Chemistry")
    menu.popup(window.mapToGlobal(window.rect().topLeft()))
    settle()
    menu.grab().save(str(out / "menu-chemistry.png"))
    menu.hide()
    conformers = window.open_tool("conformers")
    settle()
    conformers.grab().save(str(out / "conformers.png"))
    conformers.count.setValue(5)
    job = conformers.generate(conformers.settings())
    while job.running:
        settle()
    conformers.reject()
    window.view.scrollTo(window.workspace.model.index(0, len(table.column_names) - 1))
    settle()
    window.grab().save(str(out / "conformers-table.png"))
    scatter = window.open_tool("scatter")
    scatter.set_spec(x="pIC50", y="pIC50")
    settle()
    window.grab().save(str(out / "plot.png"))
    stats = window.open_tool("statistics")
    stats.scope.setCurrentIndex(stats.scope.findData("all"))
    for tab in (stats.summary, stats.correlation, stats.outliers):
        stats.tabs.setCurrentWidget(tab)
        stats.refresh()
        settle()
        name = stats.tabs.tabText(stats.tabs.currentIndex()).lower()
        window.grab().save(str(out / f"statistics-{name}.png"))
    for tool in (
        *("mpo", "binning", "selectivity", "similarity", "diverse_subset", "cluster"),
        *("rgroups", "free_wilson", "hansch", "mmp", "sali"),
        *("boiled_egg", "space_map", "som"),
    ):
        dialog = window.open_tool(tool)
        if tool == "mpo":
            dialog.add_property()
            dialog.add_property()
        settle()
        dialog.grab().save(str(out / f"{tool}.png"))
        dialog.reject()
    import threading

    from chemsuite.app.tools import ToolContext
    from chemsuite.core.jobs import check

    release = threading.Event()

    def slow(job):
        for i in range(10**6):
            if release.wait(0.01):
                return None
            job.report(min(i, 60), 100)
            check(job.cancel_token)
        return None

    def broken(job):
        raise ValueError("column 'pIC50' has no numbers")

    ctx = ToolContext(window)
    ctx.run("Calculating 12 descriptor(s)…", lambda job: None, lambda r: ctx.show_message(
        "Added 12 descriptor column(s) for 7 row(s)"))  # fmt: skip
    failed = ctx.run("Hansch analysis…", broken, lambda r: None)
    while failed.running:
        settle()
    window._failure_box.close()
    running = ctx.run("Generating conformers for 7 row(s)…", slow, lambda r: None)
    ctx.run("Finding matched pairs…", lambda job: None, lambda r: None)
    ctx.run("Making a UMAP map…", lambda job: None, lambda r: None)
    for _ in range(10):
        settle()
    window.grab().save(str(out / "jobs-status.png"))
    jobs = window.show_jobs()
    settle()
    jobs.grab().save(str(out / "jobs.png"))
    jobs.close()
    release.set()
    while running.running or window._runner.queued_jobs:
        settle()
    log_window = window.show_log()
    log_window.level.setCurrentIndex(0)
    log_window.view.selectRow(
        next(i for i, e in enumerate(log_window.model.rows) if e.details) if any(
            e.details for e in log_window.model.rows) else 0
    )  # fmt: skip
    settle()
    log_window.grab().save(str(out / "log.png"))
    log_window.close()
    settings = window.open_settings()
    for row in range(settings.pages.count()):
        settings.pages.setCurrentRow(row)
        settle()
        name = settings.pages.item(row).text().lower().replace(" ", "-")
        settings.grab().save(str(out / f"settings-{name}.png"))
    page = settings.appearance
    page.theme.setCurrentIndex(page.theme.findData("dark"))
    settings.pages.setCurrentRow(0)
    settle()
    settings.grab().save(str(out / "settings-dark.png"))
    window.grab().save(str(out / "main-dark.png"))
    settings.reject()
    from chemsuite.app.custom_theme_dialog import CustomThemeDialog

    editor = CustomThemeDialog({}, None, window)
    editor.show()
    settle()
    editor.grab().save(str(out / "custom-theme.png"))
    editor.reject()
    print("\n".join(str(p) for p in sorted(out.glob("*.png"))))
    runner.shutdown()
    workspace = window.workspace
    workspace.mark_saved(workspace.table.history.checkpoint(), workspace.lists.revision)
    window.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
