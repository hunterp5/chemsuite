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
"""The window's status bar: everything on one left-aligned line.

Left to right: the document's row counts, what the running job is doing (with its progress,
"Generating conformers… (120/500)"), then the latest message. ``QStatusBar`` would put messages
in front of everything and hide its other widgets while a message shows. Here
:meth:`StatusBar.showMessage` writes into its own label instead, so the counts and the job's
progress stay visible.
"""

from __future__ import annotations

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QHBoxLayout, QLabel, QSizePolicy, QStatusBar, QWidget

SPACING = 18


class StatusBar(QStatusBar):
    def __init__(self, jobs: QWidget, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.document = QLabel()
        self.jobs = jobs
        self.message = QLabel()
        # A long message must not widen the window: it is cut off instead.
        self.message.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        holder = QWidget()
        row = QHBoxLayout(holder)
        row.setContentsMargins(4, 0, 4, 0)
        row.setSpacing(SPACING)
        row.addWidget(self.document)
        row.addWidget(self.jobs)
        row.addWidget(self.message, 1)
        self.addWidget(holder, 1)
        self._clear_timer = QTimer(self)
        self._clear_timer.setSingleShot(True)
        self._clear_timer.timeout.connect(self.clearMessage)

    # QStatusBar's message API, redirected to our own label (Python callers only; the window
    # routes status tips here too, see MainWindow.event).
    def showMessage(self, text: str, timeout: int = 0) -> None:
        self.message.setText(text)
        self.message.setToolTip(text)
        if timeout > 0:
            self._clear_timer.start(timeout)
        else:
            self._clear_timer.stop()
        self.messageChanged.emit(text)

    def clearMessage(self) -> None:
        self._clear_timer.stop()
        if self.message.text():
            self.message.clear()
            self.message.setToolTip("")
            self.messageChanged.emit("")

    def currentMessage(self) -> str:
        return self.message.text()
