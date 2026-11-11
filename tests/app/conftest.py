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

"""Qt test setup.

Set ``CHEMSUITE_NO_GUI=1`` to run the model-level tests under a ``QCoreApplication`` where no Qt
platform plugin can load (some sandboxes); tests marked ``gui`` are then skipped. CI and
developer machines run everything offscreen (set here unless ``QT_QPA_PLATFORM`` is already
set), so tests never open windows on screen.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
NO_GUI = bool(os.environ.get("CHEMSUITE_NO_GUI"))
gui = pytest.mark.skipif(NO_GUI, reason="needs a Qt GUI platform (CHEMSUITE_NO_GUI is set)")


@pytest.fixture(scope="session")
def qapp_cls():
    if NO_GUI:
        from PySide6.QtCore import QCoreApplication

        return QCoreApplication
    from PySide6.QtWidgets import QApplication

    return QApplication


@pytest.fixture(autouse=True)
def no_modal_dialogs(monkeypatch):
    """A modal dialog nobody answers would hang the run; make it fail the test instead."""
    if NO_GUI:
        return

    from PySide6.QtWidgets import QFileDialog, QInputDialog, QMessageBox

    def refuse(owner: str, name: str):
        def raiser(*args, **kwargs):
            raise AssertionError(f"unexpected modal dialog: {owner}.{name}{args[1:3]!r}")

        return raiser

    for name in ("question", "warning", "critical", "information"):
        monkeypatch.setattr(QMessageBox, name, refuse("QMessageBox", name))
    for name in ("getInt", "getText", "getItem", "getDouble"):
        monkeypatch.setattr(QInputDialog, name, refuse("QInputDialog", name))
    for name in ("getOpenFileName", "getOpenFileNames", "getSaveFileName"):
        monkeypatch.setattr(QFileDialog, name, refuse("QFileDialog", name))


@pytest.fixture(autouse=True)
def private_settings(tmp_path, monkeypatch):
    """Each test gets empty settings and no machine policy; the user's own are never touched."""
    monkeypatch.setenv("CHEMSUITE_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.setenv("CHEMSUITE_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("CHEMSUITE_POLICY", str(tmp_path / "no-policy.json"))
    if NO_GUI:
        yield
        return
    from chemsuite.app.shell import theme
    from chemsuite.app.shell.settings_hub import reset_hub

    reset_hub()
    yield
    reset_hub()
    theme.reset()


@pytest.fixture(autouse=True)
def clear_clipboard():
    """Qt's offscreen platform crashes at exit if the clipboard still owns mime data."""
    yield
    if NO_GUI:
        return
    from PySide6.QtGui import QGuiApplication

    if QGuiApplication.instance() is not None:
        QGuiApplication.clipboard().clear()


@pytest.fixture(autouse=True)
def graphics_state():
    """Each test starts with 3D as available as the session found it: a 3D view's ``failed``
    signal (they cannot draw offscreen) arrives through a timer, sometimes in a later test,
    and would turn 3D off for every test after it."""
    if NO_GUI:
        yield
        return
    from chemsuite.app.shell import graphics

    before = graphics.problem()
    yield
    graphics.set_problem(before)


@pytest.fixture
def window(qtbot):
    """A main window on a small table (ID, pIC50), closed without a save prompt."""
    from chemsuite.app.document.main_window import MainWindow, Workspace
    from chemsuite.app.shell.jobs import JobRunner
    from chemsuite.app.table.table_view import StructureImages

    runner = JobRunner()
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": ["A", "B", "C"], "pIC50": ["5", "6", "7"]}, undoable=False)
    w = MainWindow(runner, StructureImages(), Workspace(table))
    # Not qtbot.addWidget: pytest-qt would close it before this fixture can mark it clean.
    yield w
    w.workspace.mark_saved(w.workspace.table.history.checkpoint(), w.workspace.lists.revision)
    w.close()
    runner.shutdown()
    destroy_now(w)


def destroy_now(widget) -> None:
    """Delete *widget* (and its docks, popups, panels) while the test's QApplication runs.

    Left for the interpreter's exit, dozens of windows are destroyed in no particular order,
    after Python has half torn down their wrappers: an intermittent crash at exit.
    """
    from PySide6.QtCore import QCoreApplication, QEvent

    widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


@pytest.fixture(scope="session", autouse=True)
def delete_leftover_windows():
    """At the end of the run, delete any top-level widget a test left behind (see above)."""
    yield
    if NO_GUI:
        return
    from PySide6.QtWidgets import QApplication

    if QApplication.instance() is None:
        return
    for widget in QApplication.topLevelWidgets():
        widget.hide()
        widget.deleteLater()
    from PySide6.QtCore import QCoreApplication, QEvent

    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
