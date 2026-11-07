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

"""Shared look of side panes (Filters, Advanced Search, Browser …).

Native macOS inspector style: 12 px margins, 8 px between controls, flat sections divided
by hairlines rather than boxed cards, and pop-up buttons that draw their own menu arrow.
"""

from __future__ import annotations

from PySide6.QtWidgets import QFrame, QLayout, QMenu, QPushButton, QWidget

MARGINS = (12, 10, 12, 12)
SPACING = 8


def pane_layout(layout: QLayout) -> QLayout:
    """Give a pane's outer layout the standard margins and spacing."""
    layout.setContentsMargins(*MARGINS)
    layout.setSpacing(SPACING)
    return layout


def section(frame: QFrame, name: str) -> None:
    """A flat section: no box, a hairline under it (the last one's rule ends the list)."""
    frame.setObjectName(name)
    frame.setFrameShape(QFrame.Shape.NoFrame)
    frame.setStyleSheet(f"QFrame#{name} {{ border: none; border-bottom: 1px solid palette(mid); }}")


def menu_button(text: str, menu: QMenu | None = None, parent: QWidget | None = None):
    """A push button that opens *menu* (native pop-up arrow); returns (button, menu)."""
    button = QPushButton(text, parent)
    menu = menu or QMenu(button)
    button.setMenu(menu)
    return button, menu
