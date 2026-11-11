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
"""Fetch from PDB…: download an entry by its ID and add it to the document."""

from __future__ import annotations

from PySide6.QtWidgets import QLineEdit

from chemsuite.app.dialogs import ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.download import DownloadError
from chemsuite.protein import fetch as pdb
from chemsuite.protein.io import ProteinFileError


class FetchDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Fetch from PDB", guide="proteins")
        self.ctx = ctx
        entry = self.section("Entry")
        self.code = entry.row("PDB ID", QLineEdit(),
                              help="Four characters, such as 1ATP, or pdb_00001atp")  # fmt: skip
        self.code.setPlaceholderText("1ATP")
        self.code.textChanged.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        if not self.ctx.network_allowed():
            return "Downloads are turned off on this computer (IT policy)"
        return pdb.problem(self.code.text()) or ""

    def verb(self) -> str:
        code = pdb.normalise(self.code.text())
        return f"Fetch {code}" if code else "Fetch"

    def run(self):
        code = pdb.normalise(self.code.text()) or ""
        return fetch(self.ctx, code)


def fetch(ctx: ToolContext, code: str):
    """Download entry *code* in the background and add it (a message says what went wrong)."""
    from .tool import add_protein, ligand_work

    rows_of = ligand_work(ctx)

    def work(job):
        try:
            structure = pdb.fetch(code, cancel=job.cancel_token)
        except (DownloadError, ProteinFileError) as err:
            return None, None, str(err)
        return structure, rows_of(structure, job), ""

    def done(result) -> None:
        structure, rows, error = result
        if structure is None:
            ctx.show_message(f"Could not fetch {code}: {error}", 10_000)
        else:
            add_protein(ctx, structure, pdb.name_of(structure, code), f"PDB {code}", rows)

    return ctx.run(f"Fetching {code}…", work, done)
