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

"""Application entry point: ``chemsuite [files…]`` or ``python -m chemsuite.app.shell.main``."""

from __future__ import annotations

import logging
import multiprocessing
import os
import platform
import sys
from collections.abc import Iterable
from importlib import metadata
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from PySide6.QtWidgets import QApplication

    from ..document.main_window import MainWindow
    from .jobs import JobRunner


def main(argv: list[str] | None = None) -> int:
    multiprocessing.freeze_support()  # required for the process pool in frozen builds
    from chemsuite.core.jobs.parallel import limit_app_threads

    limit_app_threads()  # before numpy and friends load: leave a core for the window
    argv = list(sys.argv if argv is None else argv)
    from chemsuite.core import logbook
    from chemsuite.core.settings import Policy, SettingsStore

    logbook.setup_logging()  # at once, so start-up problems are logged too

    from PySide6.QtCore import QCoreApplication, qVersion
    from PySide6.QtWidgets import QApplication

    from ..document.main_window import MainWindow
    from ..table.table_view import StructureImages
    from . import graphics
    from .error_report import install_crash_handler
    from .jobs import JobRunner
    from .settings_hub import config_folder, settings_hub

    QCoreApplication.setApplicationName("chemsuite")  # before config_folder()
    QCoreApplication.setOrganizationName("chemsuite")
    mode = SettingsStore(config_folder(), Policy.load()).effective().performance.graphics
    graphics.prepare(mode)  # before the QApplication: the graphics are chosen once
    app = QApplication(argv)
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
    graphics.check()  # 3D views and plots read graphics.problem()
    logbook.audit("start", version=_version())
    runner = JobRunner(app, keep_awake=True)
    window = MainWindow(runner, StructureImages())
    paths = [Path(arg) for arg in argv[1:] if not arg.startswith("-")]
    return run_app(app, runner, window, paths)


def run_app(
    app: QApplication, runner: JobRunner, window: MainWindow, paths: Iterable[Path] = ()
) -> int:
    """Show *window*, run until the last document window closes, then stop everything.

    Everything ends with the app: jobs are cancelled, worker processes and external programs
    (Vina, the AiZynthFinder worker) are ended whatever they are doing. If a job's thread is
    stuck in a call that never looks at its cancel token, the process ends without waiting for
    it: documents and settings are always written atomically, so nothing is half-saved.
    """
    from chemsuite.core import logbook
    from chemsuite.retro.runner import close_workers

    from . import quitting

    quitting.state.with_last_window = True
    window.show()
    for path in paths:
        window.open_path(path)
    code = app.exec()
    close_workers()  # asks the warm AiZynthFinder worker to quit before it is ended
    clean = runner.shutdown()
    logbook.audit("exit", code=code, forced=not clean)
    if clean:
        return code
    logging.getLogger("chemsuite").warning("forced exit: a job did not stop in time")
    logging.shutdown()
    os._exit(code)


def _version() -> str:
    try:
        return metadata.version("chemsuite")
    except metadata.PackageNotFoundError:
        return "(development)"


if __name__ == "__main__":
    sys.exit(main())
