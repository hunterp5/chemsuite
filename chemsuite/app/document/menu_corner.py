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
"""Layout and Log, at the right end of the menu bar.

Windows and Linux show them as two buttons in the menu bar's corner. macOS has one menu bar
for the whole app, so they go in a Window menu there instead.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QAction, QColor, QIcon, QPainter, QPalette, QPixmap
from PySide6.QtWidgets import QHBoxLayout, QMenu, QToolButton, QWidget

from ..shell.jobs import QUEUED, RUNNING
from ..shell.jobs_ui import FAILED_COLOUR, JobsIndicator
from . import dock_layout

if TYPE_CHECKING:
    from .main_window import MainWindow

#: menu-bar buttons look like the menu titles next to them
BUTTON_STYLE = (
    "QToolButton { border: none; border-radius: 4px; padding: 2px 8px; }"
    "QToolButton:hover, QToolButton:open { background: palette(midlight); }"
    "QToolButton::menu-indicator { image: none; width: 0px; }"
)


def layout_menu(window: MainWindow) -> QMenu:
    menu = QMenu("&Layout", window)
    menu.aboutToShow.connect(lambda: dock_layout.fill_layout_menu(menu, window))
    dock_layout.fill_layout_menu(menu, window)
    return menu


def install(window: MainWindow, layout: QMenu, log: QAction) -> None:
    """Put Layout and Log in the corner of *window*'s menu bar, or a Window menu on macOS."""
    bar = window.menuBar()
    if bar.isNativeMenuBar():
        menu = bar.addMenu("&Window")
        menu.addMenu(layout)
        menu.addAction(log)
        return
    corner = QWidget(bar)
    row = QHBoxLayout(corner)
    row.setContentsMargins(0, 0, 8, 0)
    row.setSpacing(2)
    layout_button = QToolButton()
    layout_button.setText("Layout")
    layout_button.setMenu(layout)
    layout_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
    layout_button.setToolTip("Arrange the panes around the table")
    log_button = LogButton(log, window.jobs_indicator)
    for button in (layout_button, log_button):
        button.setAutoRaise(True)
        button.setStyleSheet(BUTTON_STYLE)
        button.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        row.addWidget(button)
    bar.setCornerWidget(corner, Qt.Corner.TopRightCorner)
    window.layout_button, window.log_button = layout_button, log_button


class LogButton(QToolButton):
    """Opens the Log window; a dot shows jobs running (accent) or a failure not yet seen
    (red)."""

    def __init__(self, action: QAction, indicator: JobsIndicator) -> None:
        super().__init__()
        self.indicator = indicator
        self.setDefaultAction(action)
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._shown: str | None = None
        indicator.refreshed.connect(self.refresh)
        self.refresh()

    def state(self) -> str:
        """``failed``, ``busy`` or ``""``."""
        if self.indicator.unseen_failures():
            return "failed"
        if any(j.state in (RUNNING, QUEUED) for j in self.indicator.runner.jobs):
            return "busy"
        return ""

    def refresh(self) -> None:
        state = self.state()
        if state == self._shown:  # the indicator refreshes often while jobs run
            return
        self._shown = state
        colour = {
            "failed": FAILED_COLOUR,
            "busy": self.palette().color(QPalette.ColorRole.Highlight),
        }
        self.setIcon(_dot(colour[state], self) if state else QIcon())
        tips = {"failed": "A job failed. Show the jobs and messages.",
                "busy": "Jobs are running. Show the jobs and messages."}  # fmt: skip
        self.setToolTip(tips.get(state, "Show the jobs and messages"))


def _dot(colour: QColor, widget: QWidget) -> QIcon:
    size, ratio = 8, widget.devicePixelRatioF()
    pixmap = QPixmap(round(size * ratio), round(size * ratio))
    pixmap.setDevicePixelRatio(ratio)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(colour)
    painter.drawEllipse(QRectF(0.5, 0.5, size - 1, size - 1))
    painter.end()
    return QIcon(pixmap)
