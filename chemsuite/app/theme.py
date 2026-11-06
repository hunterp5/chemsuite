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
"""Application themes: System, Light, Dark, Groovy (a random palette) and named custom themes.

System keeps the platform's own style and follows the OS light/dark switch. The others use the
Fusion style with an explicit palette, so they look the same on every OS. Widgets that draw
their own colours read them from the palette, and repaint on ``ApplicationPaletteChange``.
"""

from __future__ import annotations

import colorsys
import random

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QProxyStyle,
    QStyle,
    QStyleFactory,
    QStyleOptionHeader,
)

from chemsuite.core.settings import DARK, GROOVY, LIGHT, SYSTEM, Appearance, custom_theme_name

Role = QPalette.ColorRole

# Colours a custom theme can set (key, label, role).
CUSTOM_ROLES: tuple[tuple[str, str, QPalette.ColorRole], ...] = (
    ("window", "Window", Role.Window),
    ("window_text", "Window text", Role.WindowText),
    ("base", "Base (tables, fields)", Role.Base),
    ("alternate_base", "Alternate rows", Role.AlternateBase),
    ("text", "Text", Role.Text),
    ("button", "Button", Role.Button),
    ("button_text", "Button text", Role.ButtonText),
    ("highlight", "Selection", Role.Highlight),
    ("highlighted_text", "Selected text", Role.HighlightedText),
    ("mid", "Borders", Role.Mid),
    ("light", "Light", Role.Light),
    ("dark", "Dark", Role.Dark),
    ("link", "Link", Role.Link),
    ("tooltip_base", "Tooltip background", Role.ToolTipBase),
    ("tooltip_text", "Tooltip text", Role.ToolTipText),
)

_original: dict[str, object] = {}  # the platform's style, palette and font, saved on first use


def _remember(app: QApplication) -> None:
    if not _original:
        _original["style"] = app.style().name()
        _original["palette"] = QPalette(app.palette())
        _original["font"] = QFont(app.font())


def fusion_palette() -> QPalette:
    style = QStyleFactory.create("Fusion")
    return style.standardPalette() if style is not None else QPalette()


def light_palette() -> QPalette:
    """A true light palette, fixed colours.

    Not Fusion's ``standardPalette()``: since Qt 6.5 that follows the OS colour scheme, so with
    macOS or Windows in dark mode "Light" came out as a grey dark theme.
    """
    text = QColor(28, 28, 30)
    disabled = QColor(160, 160, 165)
    p = QPalette()
    for role, colour in (
        (Role.Window, QColor(246, 246, 247)),
        (Role.WindowText, text),
        (Role.Base, QColor(255, 255, 255)),
        (Role.AlternateBase, QColor(245, 247, 250)),
        (Role.ToolTipBase, QColor(255, 255, 240)),
        (Role.ToolTipText, text),
        (Role.PlaceholderText, QColor(128, 128, 134)),
        (Role.Text, text),
        (Role.Button, QColor(252, 252, 252)),
        (Role.ButtonText, text),
        (Role.BrightText, QColor(220, 30, 30)),
        (Role.Link, QColor(0, 102, 204)),
        (Role.Highlight, QColor(48, 140, 198)),
        (Role.HighlightedText, QColor(255, 255, 255)),
        (Role.Mid, QColor(196, 196, 200)),
        (Role.Midlight, QColor(230, 230, 233)),
        (Role.Dark, QColor(160, 160, 165)),
        (Role.Light, QColor(255, 255, 255)),
        (Role.Shadow, QColor(110, 110, 115)),
    ):
        p.setColor(role, colour)
    for role in (Role.Text, Role.ButtonText, Role.WindowText):
        p.setColor(QPalette.ColorGroup.Disabled, role, disabled)
    p.setColor(QPalette.ColorGroup.Disabled, Role.Base, QColor(240, 240, 242))
    return p


def dark_palette() -> QPalette:
    """Fusion's layout with dark surfaces; selection and link colours as in light mode."""
    light = light_palette()
    text = QColor(240, 240, 240)
    disabled = QColor(128, 128, 128)
    p = QPalette()
    for role, colour in (
        (Role.Window, QColor(53, 53, 53)),
        (Role.WindowText, text),
        (Role.Base, QColor(35, 35, 35)),
        (Role.AlternateBase, QColor(45, 45, 45)),
        (Role.ToolTipBase, QColor(53, 53, 53)),
        (Role.ToolTipText, text),
        (Role.PlaceholderText, QColor(150, 150, 150)),
        (Role.Text, text),
        (Role.Button, QColor(68, 68, 68)),
        (Role.ButtonText, text),
        (Role.BrightText, light.color(Role.BrightText)),
        (Role.Link, QColor(110, 170, 255)),
        (Role.Highlight, light.color(Role.Highlight)),
        (Role.HighlightedText, light.color(Role.HighlightedText)),
        (Role.Mid, QColor(128, 128, 128)),
        (Role.Dark, QColor(30, 30, 30)),
        (Role.Light, QColor(75, 75, 75)),
        (Role.Shadow, QColor(15, 15, 15)),
    ):
        p.setColor(role, colour)
    for role in (Role.Text, Role.ButtonText, Role.WindowText):
        p.setColor(QPalette.ColorGroup.Disabled, role, disabled)
    return p


def _hsv(h: float, s: float, v: float) -> QColor:
    r, g, b = colorsys.hsv_to_rgb(h % 1.0, min(1.0, max(0.0, s)), min(1.0, max(0.0, v)))
    return QColor(int(r * 255), int(g * 255), int(b * 255))


def contrasting_text(background: QColor) -> QColor:
    luminance = (
        0.2126 * background.red() + 0.7152 * background.green() + 0.0722 * background.blue()
    ) / 255.0
    return QColor(18, 18, 22) if luminance > 0.55 else QColor(250, 248, 255)


def groovy_palette(rng: random.Random | None = None) -> QPalette:
    """A random, saturated palette around one hue, with readable text on every surface."""
    rng = rng or random.Random()
    hue = rng.random()

    def shift(delta: float, s: tuple[float, float], v: tuple[float, float]) -> QColor:
        return _hsv(hue + delta, rng.uniform(*s), rng.uniform(*v))

    window = shift(0.00, (0.55, 0.95), (0.35, 0.75))
    base = shift(0.12, (0.45, 0.90), (0.22, 0.55))
    alternate = shift(0.22, (0.50, 0.95), (0.30, 0.65))
    button = shift(0.35, (0.60, 1.00), (0.40, 0.85))
    highlight = shift(0.55, (0.75, 1.00), (0.55, 0.95))
    link = shift(0.70, (0.70, 1.00), (0.55, 0.95))
    mid = shift(0.08, (0.25, 0.55), (0.35, 0.60))
    light = shift(0.05, (0.35, 0.70), (0.70, 0.95))
    dark = shift(0.02, (0.40, 0.80), (0.12, 0.30))
    tip = shift(0.40, (0.50, 0.90), (0.35, 0.70))
    disabled = QColor(mid.red(), mid.green(), mid.blue(), 160)
    p = QPalette()
    for role, colour in (
        (Role.Window, window),
        (Role.WindowText, contrasting_text(window)),
        (Role.Base, base),
        (Role.AlternateBase, alternate),
        (Role.ToolTipBase, tip),
        (Role.ToolTipText, contrasting_text(tip)),
        (Role.Text, contrasting_text(base)),
        (Role.Button, button),
        (Role.ButtonText, contrasting_text(button)),
        (Role.BrightText, QColor(255, 255, 80)),
        (Role.Link, link),
        (Role.Highlight, highlight),
        (Role.HighlightedText, contrasting_text(highlight)),
        (Role.Mid, mid),
        (Role.Dark, dark),
        (Role.Light, light),
        (Role.Shadow, QColor(10, 5, 20)),
    ):
        p.setColor(role, colour)
    for role in (Role.Text, Role.ButtonText, Role.WindowText):
        p.setColor(QPalette.ColorGroup.Disabled, role, disabled)
    return p


def palette_colours(palette: QPalette) -> dict[str, str]:
    """The custom-theme colours of *palette* (to start a custom theme from it)."""
    return {key: palette.color(role).name() for key, _label, role in CUSTOM_ROLES}


def default_custom_colours() -> dict[str, str]:
    return palette_colours(light_palette())


def custom_palette(colours: dict[str, str]) -> QPalette:
    merged = default_custom_colours()
    merged.update({k: v for k, v in colours.items() if k in merged and QColor(v).isValid()})
    p = QPalette()
    for key, _label, role in CUSTOM_ROLES:
        p.setColor(role, QColor(merged[key]))
    mid = QColor(merged["mid"])
    disabled = QColor(mid.red(), mid.green(), mid.blue(), 160)
    p.setColor(Role.BrightText, QColor(255, 255, 80))
    p.setColor(Role.Shadow, QColor(20, 20, 20))
    for role in (Role.Text, Role.ButtonText, Role.WindowText):
        p.setColor(QPalette.ColorGroup.Disabled, role, disabled)
    return p


def palette_for(appearance: Appearance, rng: random.Random | None = None) -> QPalette | None:
    """The palette for the chosen theme (None: the platform's own)."""
    theme = appearance.theme
    if theme == LIGHT:
        return light_palette()
    if theme == DARK:
        return dark_palette()
    if theme == GROOVY:
        return groovy_palette(rng)
    name = custom_theme_name(theme)
    if name is not None:
        colours = appearance.custom_themes.get(name)
        return custom_palette(colours) if colours is not None else None
    return None


_applied: tuple | None = None  # what is on screen now, so unchanged settings do nothing


class PolishStyle(QProxyStyle):
    """App-wide details on top of the platform or Fusion style.

    Horizontal headers whose text Qt would centre by default are left-aligned, as in native
    macOS tables; a header that asks for an alignment of its own keeps it.
    """

    def drawControl(self, element, option, painter, widget=None) -> None:
        if (
            element == QStyle.ControlElement.CE_HeaderLabel
            and isinstance(option, QStyleOptionHeader)
            and option.orientation == Qt.Orientation.Horizontal
            and option.textAlignment & Qt.AlignmentFlag.AlignHCenter
        ):
            option = QStyleOptionHeader(option)
            option.textAlignment = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        super().drawControl(element, option, painter, widget)


def _set_style(app: QApplication, name: str) -> None:
    app.setStyle(PolishStyle(name))


def apply_appearance(app: QApplication, appearance: Appearance, *, reroll: bool = False) -> None:
    """Style, palette and font for the whole application.

    Groovy draws a new palette only when first chosen or when *reroll* is set, so changing an
    unrelated setting does not repaint everything in new colours.
    """
    global _applied
    _remember(app)
    colours = appearance.custom_themes.get(custom_theme_name(appearance.theme) or "")
    key = (appearance.theme, repr(sorted(colours.items())) if colours else None)
    if reroll or _applied is None or _applied[0] != key:
        palette = palette_for(appearance)
        if palette is None or appearance.theme == SYSTEM:
            _set_style(app, str(_original["style"]))
            app.setPalette(QPalette(_original["palette"]))
        else:
            _set_style(app, "Fusion")
            app.setPalette(palette)
        _set_colour_scheme(app, None if palette is None else palette)
    font = QFont(_original["font"])
    if appearance.app_font_pt:
        font.setPointSize(appearance.app_font_pt)
    if app.font() != font:
        app.setFont(font)
    _applied = (key,)


def _set_colour_scheme(app: QApplication, palette: QPalette | None) -> None:
    """Tell the platform whether the theme is light or dark, so what Qt does not paint itself
    (macOS menus and title bars, native dialogs) matches. None: follow the OS again."""
    hints = app.styleHints()
    if not hasattr(hints, "setColorScheme"):  # Qt before 6.8
        return
    if palette is None:
        hints.unsetColorScheme()
    else:
        hints.setColorScheme(Qt.ColorScheme.Dark if is_dark(palette) else Qt.ColorScheme.Light)


def is_dark(palette: QPalette) -> bool:
    return palette.color(Role.Base).lightness() < 128


def reset() -> None:
    """Back to the platform's look (tests)."""
    global _applied
    app = QApplication.instance()
    if _original and app is not None:
        app.setStyle(str(_original["style"]))
        app.setPalette(QPalette(_original["palette"]))
        app.setFont(QFont(_original["font"]))
        _set_colour_scheme(app, None)
    _applied = None
