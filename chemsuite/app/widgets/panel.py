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
"""The base of every tool panel, docked beside the table or in a window of its own.

A host (``PanelDock``, ``ToolWindow``, the document's saved session) asks a panel only what
this class declares. A panel overrides what it has; the defaults do nothing, so the host never
has to ask whether a panel has something.
"""

from __future__ import annotations

from PySide6.QtCore import QPoint, Qt, Signal
from PySide6.QtWidgets import QFormLayout, QFrame, QWidget


class OptionsPopup(QFrame):
    """A panel's options, in a popup that closes when you click elsewhere (opened by the
    header's settings button: a panel exposes one as ``ToolPanel.options``)."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent, Qt.WindowType.Popup)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.form = QFormLayout(self)
        self.form.setContentsMargins(10, 10, 10, 10)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)

    def popup(self, anchor: QWidget) -> None:
        """Show below *anchor* (the settings button), kept on its screen."""
        self.adjustSize()
        below = anchor.mapToGlobal(QPoint(anchor.width() - self.width(), anchor.height()))
        screen = anchor.screen().availableGeometry()
        x = min(max(below.x(), screen.left()), screen.right() - self.width())
        y = below.y()
        if y + self.height() > screen.bottom():  # no room below: open above
            y = anchor.mapToGlobal(QPoint(0, 0)).y() - self.height()
        self.move(x, max(y, screen.top()))
        self.show()


class ToolPanel(QWidget):
    #: the panel's title changed: the dock or window follows it
    titleChanged = Signal(str)
    #: saved with the document and reopened by the tool of this id; "" is never saved
    tool_id: str = ""
    #: the popup behind the header's gear button, if the panel has one
    options: OptionsPopup | None = None
    #: buttons shown in the dock's header, if any
    toolbar: QWidget | None = None
    #: the dock's title can be renamed by the user (``set_custom_title``)
    renamable: bool = False

    def state(self) -> dict:
        """What reopens the panel with the document (``Tool.restore`` takes it back)."""
        return {}

    def set_custom_title(self, text: str | None) -> None:
        """The user renamed the dock (None: back to the panel's own title)."""

    def attach_header(self, attached: bool) -> None:
        """The panel's header moved into the dock's title bar (or back)."""

    def close_panel(self) -> None:
        """The dock or window is closing for good: stop work, let go of the document."""
