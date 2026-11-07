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
import html
from collections import OrderedDict

from PySide6.QtGui import QGuiApplication
from rdkit import Chem

from chemsuite.core.depict import render_png

from .theme import is_dark

#: logical size of the drawing in the tooltip
PREVIEW_WIDTH = 360
PREVIEW_HEIGHT = 270
#: drawn a little larger than the table's thumbnails
_ZOOM = 1.3
_KEEP = 32

_cache: OrderedDict[tuple[str, bool, float], str] = OrderedDict()


def preview_html(mol: Chem.Mol, text: str) -> str:
    """Tooltip HTML: *mol* drawn at about 360 x 270 points, with *text* (its SMILES) below."""
    app = QGuiApplication.instance()
    screen = QGuiApplication.primaryScreen() if app is not None else None
    dpr = float(screen.devicePixelRatio()) if screen is not None else 1.0
    dark = is_dark(QGuiApplication.palette()) if app is not None else False
    key = (text, dark, dpr)
    image = _cache.get(key)
    if image is None:
        png = render_png(
            mol,
            round(PREVIEW_WIDTH * dpr),
            round(PREVIEW_HEIGHT * dpr),
            scale=dpr * _ZOOM,
            dark=dark,
        )
        image = base64.b64encode(png).decode("ascii")
        _cache[key] = image
        while len(_cache) > _KEEP:
            _cache.popitem(last=False)
    else:
        _cache.move_to_end(key)
    shown = html.escape(text if len(text) <= 200 else text[:200] + "…")
    return (
        f'<img src="data:image/png;base64,{image}" '
        f'width="{PREVIEW_WIDTH}" height="{PREVIEW_HEIGHT}"><br>'
        f'<span style="font-family: monospace;">{shown}</span>'
    )
