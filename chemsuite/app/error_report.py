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
"""Showing errors: a short message, with the details one click away and never a raw traceback.

:func:`show_error` is used for failed jobs and for uncaught exceptions (from any thread, via
:func:`install_crash_handler`). It never blocks: other work carries on while it is open.
"""

from __future__ import annotations

import traceback

import shiboken6
from PySide6.QtCore import QObject, Qt, QUrl, Signal, Slot
from PySide6.QtGui import QDesktopServices, QGuiApplication
from PySide6.QtWidgets import QApplication, QMessageBox, QWidget

from chemsuite.core import logbook

MAX_OPEN = 3  # more errors at once than this go to the log only


def show_log_folder() -> None:
    folder = logbook.log_folder()
    folder.mkdir(parents=True, exist_ok=True)
    QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))


def details_of(error: BaseException) -> str:
    return "".join(traceback.format_exception(error))


_open: list[QMessageBox] = []


def show_error(
    parent: QWidget | None,
    summary: str,
    reason: str = "",
    details: str = "",
    *,
    icon: QMessageBox.Icon = QMessageBox.Icon.Warning,
) -> QMessageBox | None:
    """A non-blocking error box: *summary*, *reason*, and *details* behind "Show Details…"."""
    _open[:] = [b for b in _open if shiboken6.isValid(b) and b.isVisible()]
    if len(_open) >= MAX_OPEN:
        return None
    box = QMessageBox(parent)
    box.setIcon(icon)
    box.setWindowTitle("chemsuite")
    box.setText(summary)
    if reason:
        box.setInformativeText(reason)
    if details:
        box.setDetailedText(details)
        copy = box.addButton("Copy Details", QMessageBox.ButtonRole.ActionRole)
        copy.clicked.disconnect()  # an action button would otherwise close the box
        copy.clicked.connect(lambda: QGuiApplication.clipboard().setText(f"{summary}\n{details}"))
    folder = box.addButton("Show Log Folder", QMessageBox.ButtonRole.ActionRole)
    folder.clicked.disconnect()
    folder.clicked.connect(show_log_folder)
    box.addButton(QMessageBox.StandardButton.Close)
    box.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
    box.open()
    _open.append(box)
    # Closed boxes delete themselves (and die with their window): forget them at once.
    box.destroyed.connect(lambda *_: _open.remove(box) if box in _open else None)
    if parent is not None and not parent.isActiveWindow():
        QApplication.alert(parent)
    return box


class _Bridge(QObject):
    """Carries an error from any thread to the GUI thread.

    ``error`` is connected to a slot of this object (which lives on the GUI thread), so Qt
    queues the call; a plain function would run on the emitting thread.
    """

    error = Signal(str, str)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.error.connect(self.show)

    @Slot(str, str)
    def show(self, summary: str, details: str) -> None:
        show_error(
            QApplication.activeWindow(),
            "Something went wrong in chemsuite.",
            f"{summary}\n\nYour work is still open; saving a copy now is a good idea. "
            "The full details are in the log.",
            details,
            icon=QMessageBox.Icon.Critical,
        )


_bridge: _Bridge | None = None


def install_crash_handler() -> _Bridge:
    """Log every uncaught exception in full and show a short message about it."""
    global _bridge
    _bridge = _Bridge(QApplication.instance())
    bridge = _bridge
    logbook.install_excepthooks(lambda summary, details: bridge.error.emit(summary, details))
    return bridge
