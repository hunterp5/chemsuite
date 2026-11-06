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
"""The window's status bar, on one line.

Left to right: what the running job is doing (with its progress, "Generating conformers…
(120/500)"), then the latest message; the document's row and selection counts sit at the far
right. ``QStatusBar`` would put messages
in front of everything and hide its other widgets while a message shows. Here
:meth:`StatusBar.showMessage` writes into its own label instead, so the counts and the job's
progress stay visible.
"""

from __future__ import annotations

from PySide6.QtCore import QEvent, QTimer
from PySide6.QtGui import QResizeEvent
from PySide6.QtWidgets import QHBoxLayout, QLabel, QSizePolicy, QStatusBar, QWidget

SPACING = 18
#: pixels above and below the text: a slim bar that does not grow while a job shows
PADDING = 2


class StatusBar(QStatusBar):
    def __init__(self, jobs: QWidget, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.document = QLabel()
        self.jobs = jobs
        self.message = QLabel()
        # A long message must not widen the window: it is cut off instead.
        self.message.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        # Not added with addWidget: QStatusBar's own layout puts 3 px above its items, which
        # pushes the text off centre in a slim bar. The row fills the bar instead (resizeEvent).
        holder = self.holder = QWidget(self)
        row = QHBoxLayout(holder)
        row.setContentsMargins(6, 0, 6, 0)
        row.setSpacing(SPACING)
        row.addWidget(self.jobs)
        row.addWidget(self.message, 1)
        row.addWidget(self.document)
        self.setSizeGripEnabled(False)
        self._fit_height()
        self._clear_timer = QTimer(self)
        self._clear_timer.setSingleShot(True)
        self._clear_timer.timeout.connect(self.clearMessage)

    def _fit_height(self) -> None:
        height = self.fontMetrics().height()
        self.setFixedHeight(height + 2 * PADDING)
        self.jobs.setMaximumHeight(height + 2)

    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        self.holder.setGeometry(self.rect())

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.FontChange:  # Settings ▸ Appearance ▸ font size
            self._fit_height()

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
