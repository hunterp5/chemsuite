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

"""The ligand pKa model a preparation runs with: chemsuite pKa, fetched once when the chemist
agrees, or MolGpKa.

:func:`ligand_pka` runs on the GUI thread before the job starts: with the weights already in
the cache it uses them; without them it asks (once per session) whether to download them, unless
IT policy turns the network off. What it returns is called inside the job, where the download
happens with the job's progress and cancel; a download that fails leaves MolGpKa for that run
(the report's "Protonation states by" line says which model ran).
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any

from PySide6.QtWidgets import QMessageBox

from chemsuite.app.document.tools import ToolContext
from chemsuite.core.download import DownloadError
from chemsuite.ionisation.unipka import weights
from chemsuite.macromolecule.prepare.ligand_pka import LigandPka, ModelStates

from .tool import dialog_parent

log = logging.getLogger(__name__)

FOLDER = "pka"
_declined = False  # the download was declined in this session: MolGpKa until restart

Loader = Callable[[Any], LigandPka | None]  # called with the job


def _none(_job: Any) -> None:
    return None


def ligand_pka(ctx: ToolContext, ask: bool = True) -> Loader:
    """For a background job: the ligand pKa model to prepare with (None: MolGpKa)."""
    global _declined
    folder = ctx.cache_folder(FOLDER)
    path = weights.cached(folder)
    if path is None:
        if _declined or not ask or not ctx.network_allowed():
            return _none
        answer = QMessageBox.question(
            dialog_parent(ctx),
            "chemsuite pKa",
            f"Download the chemsuite pKa model ({weights.SIZE / 1e6:.0f} MB) for ligand "
            "protonation states?",
        )
        if answer != QMessageBox.StandardButton.Yes:
            _declined = True
            return _none

    def load(job: Any) -> LigandPka | None:
        found = path
        if found is None:
            try:
                found = weights.fetch(folder, cancel=job.cancel_token,
                                      report=lambda got: job.report(got, weights.SIZE))  # fmt: skip
            except DownloadError as error:
                log.warning("chemsuite pKa not downloaded, MolGpKa used: %s", error)
                return None
        return ModelStates(weights.predictor(found))

    return load
