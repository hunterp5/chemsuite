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
"""Render the Protein Viewer to PNG files in light and dark, to compare with the table and plots.

    .venv/bin/python scripts/screenshot_protein.py OUT_DIR [FILE]

FILE is a PDB or mmCIF file (default: 1ATP from the test data). The 3D view needs the GPU, so
this opens real windows on the screen (the offscreen platform cannot draw it) and closes them.
Writes ``protein-light.png``, ``protein-dark.png`` and the same with the first ligand focused.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("CHEMSUITE_CONFIG_DIR", tempfile.mkdtemp(prefix="chemsuite-shot-"))


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    from PySide6.QtWidgets import QApplication

    app = QApplication([argv[0]])
    from chemsuite.app.shell import theme
    from chemsuite.app.shell.windows import ToolWindow
    from chemsuite.core.settings import DARK, LIGHT, Appearance
    from chemsuite.core.table import Asset, PutAsset, Table
    from chemsuite.protein.chem import ligands
    from chemsuite.protein.structure import io
    from chemsuite.protein.ui.panel import ProteinPanel

    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    source = Path(argv[2]) if len(argv) > 2 else ROOT / "tests/data/protein/1atp.cif.gz"
    structure = io.read(source)
    name = str(structure.meta.get("id") or source.stem).upper()
    table = Table()
    table.edit("Open", [PutAsset(name, None, Asset(name, "protein", name, structure, {}))])
    panel = ProteinPanel(table)
    window = ToolWindow(panel, "Protein Viewer", name=lambda: name, menu_bar=panel.menu_bar)
    window.resize(1400, 900)
    window.show()

    def settle() -> None:
        for _ in range(40):
            app.processEvents()

    found = ligands.row_residues(structure)
    for look in (LIGHT, DARK):
        theme.apply_appearance(app, Appearance(theme=look))  # Fusion, as the app's own themes
        panel.session.clear_focus()
        settle()
        window.grab().save(str(out / f"protein-{look}.png"))
        if found:
            panel.session.focus(name, found[0])
            settle()
            window.grab().save(str(out / f"protein-{look}-focus.png"))
    for path in sorted(out.glob("protein-*.png")):
        print(path)
    window.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
