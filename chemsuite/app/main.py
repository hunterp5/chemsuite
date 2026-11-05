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
import platform
import sys
from importlib import metadata
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    multiprocessing.freeze_support()  # required for the process pool in frozen builds
    argv = list(sys.argv if argv is None else argv)
    from chemsuite.core import logbook

    logbook.setup_logging()  # at once, so start-up problems are logged too

    from PySide6.QtCore import qVersion
    from PySide6.QtWidgets import QApplication

    from .error_report import install_crash_handler
    from .jobs import JobRunner
    from .main_window import MainWindow
    from .settings_hub import settings_hub
    from .table_view import StructureImages

    app = QApplication(argv)
    app.setApplicationName("chemsuite")
    app.setOrganizationName("chemsuite")
    policy = settings_hub().store.policy
    log_file = logbook.setup_logging(level=logbook.level_from(policy.get("log_level")))
    audit_file = logbook.setup_audit(policy.get("audit_log"))
    install_crash_handler()
    log = logging.getLogger("chemsuite")
    log.info(
        "chemsuite %s starting: Python %s, Qt %s, %s",
        _version(),
        platform.python_version(),
        qVersion(),
        platform.platform(),
    )
    log.info("log file: %s%s", log_file, f"; audit log: {audit_file}" if audit_file else "")
    logbook.audit("start", version=_version())
    runner = JobRunner(app, keep_awake=True)
    images = StructureImages()
    app.aboutToQuit.connect(runner.shutdown)

    window = MainWindow(runner, images)
    window.show()
    for arg in argv[1:]:
        if not arg.startswith("-"):
            window.open_path(Path(arg))
    return app.exec()


def _version() -> str:
    try:
        return metadata.version("chemsuite")
    except metadata.PackageNotFoundError:
        return "(development)"


if __name__ == "__main__":
    sys.exit(main())
