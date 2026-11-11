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

"""Desktop notifications, for a long background job that ends while chemsuite is not in front.

Each platform has its own most reliable route:

- **macOS:** ``osascript``'s ``display notification``. Qt's tray-icon messages need a menu-bar
  icon and an app bundle, which a chemsuite started from Python does not have. The title and
  message go in as arguments (``on run argv``), never pasted into the script, so a job title
  cannot inject AppleScript.
- **Linux:** ``notify-send`` when it is installed, else a tray-icon message.
- **Windows:** a tray-icon balloon; the icon is shown only while the message is.

A notification is best effort: when none of these is available nothing happens, and the dock or
taskbar alert the window already raises still stands.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QSystemTrayIcon

log = logging.getLogger(__name__)

APP = "chemsuite"
TRAY_MS = 8000

_APPLESCRIPT = (
    "on run argv",
    "display notification (item 2 of argv) with title (item 1 of argv)",
    "end run",
)
_tray: QSystemTrayIcon | None = None


def command(title: str, message: str, platform: str | None = None) -> list[str] | None:
    """The program that posts the notification on *platform*, or None to use the tray."""
    platform = platform or sys.platform
    if platform == "darwin":
        script = [part for line in _APPLESCRIPT for part in ("-e", line)]
        return ["osascript", *script, title, message]
    if platform.startswith("linux"):
        found = shutil.which("notify-send")
        if found:
            return [found, "--app-name", APP, title, message]
    return None


def notify(title: str, message: str) -> bool:
    """Post a desktop notification; True when one was handed to the system."""
    argv = command(title, message)
    if argv is not None:
        try:
            subprocess.Popen(
                argv,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return True
        except OSError as err:
            log.info("could not post a notification with %s: %s", argv[0], err)
    return _tray_message(title, message)


def _tray_message(title: str, message: str) -> bool:
    global _tray
    if not (QSystemTrayIcon.isSystemTrayAvailable() and QSystemTrayIcon.supportsMessages()):
        return False
    app = QApplication.instance()
    if app is None:
        return False
    if _tray is None:
        _tray = QSystemTrayIcon(app.windowIcon(), app)
        _tray.setToolTip(APP)
    _tray.show()
    _tray.showMessage(title, message, QSystemTrayIcon.MessageIcon.Information, TRAY_MS)
    QTimer.singleShot(TRAY_MS + 1000, _tray.hide)
    return True
