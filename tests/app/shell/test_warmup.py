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
from __future__ import annotations

import subprocess
import sys
import textwrap

# In a fresh interpreter, so the modules are really imported on the warm-up thread.
SCRIPT = textwrap.dedent(
    """
    import sys
    from PySide6.QtCore import QObject, qInstallMessageHandler
    from PySide6.QtWidgets import QApplication

    messages = []
    qInstallMessageHandler(lambda _mode, _ctx, text: messages.append(text))
    app = QApplication([])
    before = set(sys.modules)
    from chemsuite.app.shell import warmup

    warmup.start().join(120)
    missing = [name for name in warmup.MODULES if name not in sys.modules]
    objects = [
        f"{name}.{key}"
        for name in set(sys.modules) - before
        if not name.startswith("__")  # this script, as multiprocessing names it
        for key, value in list(vars(sys.modules[name]).items())
        if isinstance(value, QObject)
    ]
    print("missing", missing)
    print("objects", objects)
    print("messages", messages)
    """
)


def test_the_warm_up_imports_everything_and_makes_no_qt_object_on_its_thread() -> None:
    done = subprocess.run(
        [sys.executable, "-c", SCRIPT],
        capture_output=True,
        text=True,
        timeout=180,
        check=True,
    )
    assert done.stdout.splitlines() == ["missing []", "objects []", "messages []"]
