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

"""Application entry point: ``chemsuite [files…]`` or ``python -m chemsuite.app.main``."""

from __future__ import annotations

import logging
import multiprocessing
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    multiprocessing.freeze_support()  # required for the process pool in frozen builds
    argv = list(sys.argv if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    from PySide6.QtWidgets import QApplication

    from .jobs import JobRunner
    from .main_window import MainWindow
    from .table_view import StructureImages

    app = QApplication(argv)
    app.setApplicationName("chemsuite")
    app.setOrganizationName("chemsuite")
    runner = JobRunner(app, keep_awake=True)
    images = StructureImages()
    app.aboutToQuit.connect(runner.shutdown)

    window = MainWindow(runner, images)
    window.show()
    for arg in argv[1:]:
        if not arg.startswith("-"):
            window.open_path(Path(arg))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
