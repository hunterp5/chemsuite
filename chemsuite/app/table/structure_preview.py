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

"""The large structure drawing shown when the pointer rests on a structure cell."""

from __future__ import annotations

import base64
from collections import OrderedDict
from typing import TYPE_CHECKING

from PySide6.QtGui import QGuiApplication, QPalette
from PySide6.QtWidgets import QApplication
from rdkit import Chem

from chemsuite.core.depict import aligned_mol, fitted_size, render_png

from .highlights import AtomHighlight, Highlight, marked

if TYPE_CHECKING:  # the cells draw with it; importing it here would close a circle
    from .cells import Alignment

#: largest logical size of the drawing in the tooltip; smaller molecules get a smaller box
PREVIEW_WIDTH = 360
PREVIEW_HEIGHT = 270
#: smallest box, so a tiny molecule still reads as a picture
_MIN_WIDTH = 120
_MIN_HEIGHT = 90
#: drawn a little larger than the table's thumbnails
_ZOOM = 1.3
_KEEP = 32

_cache: OrderedDict[tuple[str, bool, float, str | None, str | None], str] = OrderedDict()


def _dark() -> bool:
    """Whether tooltips have a dark background (the drawing goes on it).

    From the application's palette, not ``QToolTip.palette()``: Qt freezes the tooltip class
    palette the first time anything reads it — which this did, on the first structure hovered —
    and does not revisit it, so after a theme change in Preferences the tooltip the user sees
    is dark while its class palette still reads light.
    """
    app = QApplication.instance()
    if not isinstance(app, QApplication):
        return False
    return app.palette().color(QPalette.ColorRole.ToolTipBase).lightness() < 128


def preview_size(mol: Chem.Mol) -> tuple[int, int]:
    """The tooltip drawing's logical size: fitted to *mol* at the preview's bond length."""
    width, height = fitted_size(
        mol,
        round(PREVIEW_WIDTH / _ZOOM),
        round(PREVIEW_HEIGHT / _ZOOM),
        min_width=round(_MIN_WIDTH / _ZOOM),
        min_height=round(_MIN_HEIGHT / _ZOOM),
    )
    return round(width * _ZOOM), round(height * _ZOOM)


def preview_html(
    mol: Chem.Mol,
    key: str,
    highlight: Highlight | AtomHighlight | None = None,
    alignment: Alignment | None = None,
) -> str:
    """Tooltip HTML: *mol* alone, centered in a box fitted to it (at most 360 x 270 points).
    *key* (its SMILES) identifies the drawing in the cache.

    *highlight* and *alignment* are the cell's own (see :mod:`.highlights` and
    :class:`.cells.Alignment`), so resting the pointer on a cell enlarges the drawing it
    already shows: the same atoms marked, laid out on the same core.
    """
    screen = QGuiApplication.primaryScreen() if QGuiApplication.instance() else None
    dpr = float(screen.devicePixelRatio()) if screen is not None else 1.0
    dark = _dark()
    width, height = preview_size(mol)
    cache_key = (
        key,
        dark,
        dpr,
        highlight.key if highlight else None,
        alignment.text if alignment else None,
    )
    image = _cache.get(cache_key)
    if image is None:
        atoms: list[int] = []
        bonds: list[int] = []
        colors: dict[int, tuple[float, float, float]] = {}
        if alignment is not None:
            mol = aligned_mol(mol, alignment.core)
        if highlight is not None:
            mol, atoms, bonds, colors = marked(mol, highlight, dark)
        png = render_png(
            mol,
            round(width * dpr),
            round(height * dpr),
            scale=dpr * _ZOOM,
            dark=dark,
            highlight_atoms=atoms,
            highlight_bonds=bonds,
            highlight_colors=colors,
        )
        image = base64.b64encode(png).decode("ascii")
        _cache[cache_key] = image
        while len(_cache) > _KEEP:
            _cache.popitem(last=False)
    else:
        _cache.move_to_end(cache_key)
    return (
        f'<div align="center"><img src="data:image/png;base64,{image}" '
        f'width="{width}" height="{height}"></div>'
    )
