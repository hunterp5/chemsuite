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
"""Keep the computer from sleeping while long jobs run.

macOS runs ``caffeinate``, Linux ``systemd-inhibit`` where it exists, and Windows sets the
thread's execution state. The display may still turn off; only system sleep is held off.
Every failure is logged and ignored: staying awake is a convenience, never a requirement.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys

log = logging.getLogger(__name__)

_ES_CONTINUOUS = 0x80000000
_ES_SYSTEM_REQUIRED = 0x00000001


class KeepAwake:
    """``hold(True)`` keeps the system awake until ``hold(False)``. Call from one thread."""

    def __init__(self) -> None:
        self._held = False
        self._process: subprocess.Popen | None = None

    @property
    def held(self) -> bool:
        return self._held

    def hold(self, on: bool) -> None:
        if on == self._held:
            return
        self._held = on
        try:
            self._start() if on else self._stop()
        except Exception as err:  # noqa: BLE001 - optional nicety; never break a job over it
            log.warning("could not %s sleep prevention: %s", "start" if on else "stop", err)

    def _start(self) -> None:
        if sys.platform == "win32":
            import ctypes

            ctypes.windll.kernel32.SetThreadExecutionState(_ES_CONTINUOUS | _ES_SYSTEM_REQUIRED)
            return
        command = _command()
        if command:
            self._process = subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    def _stop(self) -> None:
        if sys.platform == "win32":
            import ctypes

            ctypes.windll.kernel32.SetThreadExecutionState(_ES_CONTINUOUS)
            return
        process, self._process = self._process, None
        if process is not None and process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                process.kill()


def _command() -> list[str] | None:
    pid = str(os.getpid())
    if sys.platform == "darwin" and shutil.which("caffeinate"):
        return ["caffeinate", "-i", "-w", pid]  # -w: exits by itself if chemsuite dies
    if shutil.which("systemd-inhibit"):
        return [
            "systemd-inhibit",
            "--what=sleep:idle",
            "--who=chemsuite",
            "--why=Running jobs",
            "--mode=block",
            "tail",
            f"--pid={pid}",
            "-f",
            os.devnull,
        ]
    return None
