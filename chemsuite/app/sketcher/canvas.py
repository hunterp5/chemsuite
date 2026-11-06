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

"""The drawing surface: paints a sketch and passes mouse and key input to the active tool.

It holds no chemistry. What to draw comes from ``chemsuite.sketch.render``, hydrogens and
valence problems from ``chemsuite.sketch.valence``, and every change is an edit from
``chemsuite.sketch.edits`` so undo is exact. Keys follow docs/sketcher/interaction-spec.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from PySide6.QtCore import QEvent, QMimeData, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontMetricsF,
    QGuiApplication,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QPainter,
    QPainterPath,
    QPen,
    QPolygonF,
)
from PySide6.QtWidgets import QToolTip, QWidget
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from chemsuite.app.theme import is_dark
from chemsuite.sketch import convert, edits, reaction, render, stereo, valence
from chemsuite.sketch.edits import Change, History
from chemsuite.sketch.model import (
    ARROW_TITLES,
    ATOM_LIST,
    ATTACH_POINTS,
    BOND_QUERIES,
    MAX_RGROUP,
    QUERY_ATOMS,
    RGROUP,
    TOPOLOGY,
    Order,
    Sketch,
    enhanced_text,
    next_enhanced,
    next_rgroup,
    parse_enhanced,
)

Point = tuple[float, float]

#: Pixels per bond length at zoom 1.
BOND_PX = 40.0
#: How close (in bond lengths) the pointer must be to pick an atom or a bond.
ATOM_PICK = 0.3
BOND_PICK = 0.22
OBJECT_PICK = 0.3  # reaction arrows and plus signs
#: A press that moves less than this many pixels is a click.
CLICK_SLOP = 4.0

#: Hover keys that set an element (interaction-spec §1). Shifted letters are upper case.
ELEMENT_KEYS = {
    "c": "C",
    "n": "N",
    "o": "O",
    "s": "S",
    "p": "P",
    "f": "F",
    "h": "H",
    "i": "I",
    "l": "Cl",
    "b": "Br",
    "d": "D",
    "S": "Si",
    "N": "Na",
    "B": "B",
}
BOND_KEYS = {
    "1": Order.SINGLE,
    "2": Order.DOUBLE,
    "3": Order.TRIPLE,
    "4": Order.AROMATIC,
    "0": Order.ANY,  # Ketcher: any-bond query
}


@dataclass(frozen=True)
class Hit:
    kind: str  # "atom", "bond", "arrow", "plus" or "" for empty canvas
    id: int = 0

    def __bool__(self) -> bool:
        return bool(self.kind)


NOTHING = Hit("")


def _next_ring(tool):
    """``t`` cycles the ring bar, starting from benzene (Ketcher)."""
    from chemsuite.app.sketcher.tools import TemplateTool
    from chemsuite.sketch.templates import RING_BAR

    if isinstance(tool, TemplateTool) and tool.template in RING_BAR:
        return RING_BAR[(RING_BAR.index(tool.template) + 1) % len(RING_BAR)]
    return RING_BAR[0]


@lru_cache(maxsize=256)
def _element_rgb(symbol: str, dark: bool) -> tuple[float, float, float]:
    options = rdMolDraw2D.MolDrawOptions()
    if dark:
        rdMolDraw2D.SetDarkMode(options)
    palette = options.getAtomPalette()
    try:
        number = Chem.GetPeriodicTable().GetAtomicNumber(symbol)
    except (RuntimeError, ValueError):
        number = 0
    r, g, b, _a = palette.get(number, palette.get(-1, (0.0, 0.0, 0.0, 1.0)))
    return r, g, b


#: CIP labelling is skipped on drawings larger than this (it can take seconds on huge ones)
STEREO_ATOM_LIMIT = 250


class SketchCanvas(QWidget):
    """Draws a ``History``'s sketch and edits it through the current tool."""

    #: show R/S and E/Z labels (shared by every sketcher in the session)
    show_stereo = True

    #: The drawing changed (an edit, undo or redo).
    changed = Signal()
    #: The tool changed (its name), so a toolbar can follow keyboard tool picks.
    tool_changed = Signal(str)
    #: A short message for the user (why a label could not be used).
    message = Signal(str)

    def __init__(self, history: History | None = None, parent: QWidget | None = None):
        super().__init__(parent)
        from chemsuite.app.sketcher.tools import SelectTool

        self.history = history if history is not None else History()
        self.history.listeners.append(self._history_changed)
        self.zoom = 1.0
        self.origin = QPointF(0.0, 0.0)  # screen position of model (0, 0)
        self.hover: Hit = NOTHING
        self.mouse: QPointF | None = None
        self.selected_atoms: set[int] = set()
        self.selected_bonds: set[int] = set()
        #: selected reaction arrows and plus signs
        self.selected_objects: set[int] = set()
        self._preview: Change | None = None
        self._label_editor = None
        self.tool = SelectTool()
        self.report = valence.Report()
        self.scene = render.Scene()
        self.stereo = stereo.StereoLabels()
        #: R/S and E/Z labels the user hid one by one: ("atom" | "bond", id)
        self.hidden_stereo: set[tuple[str, int]] = set()
        #: wedges that do not say what they seem to: bond id -> why (stereo.wedge_problems)
        self.wedge_problems: dict[int, str] = {}
        self._fitted = False
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumSize(360, 260)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent)
        self._recompute()

    # -- model -------------------------------------------------------------------------------

    @property
    def sketch(self) -> Sketch:
        return self.history.sketch

    def set_tool(self, tool) -> None:
        self.cancel_preview()
        self.tool = tool
        self.setCursor(tool.cursor)
        self.tool_changed.emit(tool.name)
        self.update()

    def notify(self, text: str) -> None:
        """A short message for the sketcher window (why an edit was not made)."""
        self.message.emit(text)

    def edit(self, name: str, apply) -> None:
        """Make one undoable edit: ``apply(change)`` records it."""
        change = self.history.change(name)
        apply(change)
        if self.history.commit(change) is None:
            self._recompute()

    def preview(self, name: str, apply) -> None:
        """Show an edit while a drag is in progress, replacing the previous preview."""
        if self._preview is not None:
            self.history.discard(self._preview)
        self._preview = self.history.change(name)
        apply(self._preview)
        self._recompute()

    def commit_preview(self) -> None:
        change, self._preview = self._preview, None
        if change is not None and self.history.commit(change) is None:
            self._recompute()

    def cancel_preview(self) -> None:
        change, self._preview = self._preview, None
        if change is not None:
            self.history.discard(change)
            self._recompute()

    def set_sketch(self, sketch: Sketch, name: str = "Replace drawing") -> None:
        self.history.replace(sketch, name)
        self.clear_selection()
        self.fit()

    def _history_changed(self) -> None:
        self.selected_atoms &= set(self.sketch.atoms)
        self.selected_bonds &= set(self.sketch.bonds)
        self.selected_objects &= set(self.sketch.arrows) | set(self.sketch.pluses)
        self._recompute()

    def _recompute(self) -> None:
        self.report = valence.check(self.sketch)
        self.scene = render.build(self.sketch, self.report)
        big = len(self.sketch.atoms) > STEREO_ATOM_LIMIT
        show = SketchCanvas.show_stereo and not big and self.report.valid
        self.stereo = stereo.labels(self.sketch) if show else stereo.StereoLabels()
        checks = not big and self.report.valid
        self.wedge_problems = stereo.wedge_problems(self.sketch) if checks else {}
        if self.hover.kind == "atom" and self.hover.id not in self.sketch.atoms:
            self.hover = NOTHING
        if self.hover.kind == "bond" and self.hover.id not in self.sketch.bonds:
            self.hover = NOTHING
        self.changed.emit()
        self.update()

    # -- coordinates -------------------------------------------------------------------------

    @property
    def scale(self) -> float:
        return BOND_PX * self.zoom

    def to_model(self, p: QPointF) -> Point:
        return ((p.x() - self.origin.x()) / self.scale, (self.origin.y() - p.y()) / self.scale)

    def to_screen(self, p: Point) -> QPointF:
        return QPointF(self.origin.x() + p[0] * self.scale, self.origin.y() - p[1] * self.scale)

    def hit_at(self, p: Point, exclude_atoms: frozenset[int] = frozenset()) -> Hit:
        """The atom or bond under a point. Atoms inside a contracted abbreviation cannot be
        hit; its label is hit through the group's attachment atom."""
        hidden = self.sketch.hidden_atoms()
        exclude = frozenset(exclude_atoms) | hidden
        atom = self.sketch.atom_at(*p, within=ATOM_PICK, exclude=exclude)
        if atom is not None:
            return Hit("atom", atom)
        bond = self.sketch.bond_at(*p, within=BOND_PICK)
        if bond is not None and not (
            {self.sketch.bonds[bond].begin, self.sketch.bonds[bond].end} & exclude
        ):
            return Hit("bond", bond)
        plus = self.sketch.plus_at(*p, within=OBJECT_PICK)
        if plus is not None:
            return Hit("plus", plus)
        arrow = self.sketch.arrow_at(*p, within=OBJECT_PICK)
        if arrow is not None:
            return Hit("arrow", arrow)
        return NOTHING

    def with_groups(self, atoms) -> set[int]:
        """``atoms`` plus every atom of a contracted abbreviation whose label is among them."""
        out = set(atoms)
        for group in self.sketch.groups.values():
            if not group.expanded and group.attach in out:
                out |= group.atoms
        return out

    def contracted_group(self, atom: int):
        """(id, group) when ``atom`` is the label of a contracted abbreviation."""
        found = self.sketch.group_of(atom)
        if found and not found[1].expanded and found[1].attach == atom:
            return found
        return None

    def fit(self) -> None:
        """Centre the drawing; shrink it if it does not fit (never enlarge past zoom 1)."""
        w, h = max(self.width(), 1), max(self.height(), 1)
        if self.sketch.is_empty():
            self.zoom = 1.0
            self.origin = QPointF(w / 2, h / 2)
        else:
            x0, y0, x1, y1 = self.sketch.bounds()
            span_x, span_y = (x1 - x0) + 2, (y1 - y0) + 2
            self.zoom = min(1.0, w / (span_x * BOND_PX), h / (span_y * BOND_PX))
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            self.origin = QPointF(w / 2 - cx * self.scale, h / 2 + cy * self.scale)
        self._fitted = True
        self.update()

    def zoom_by(self, factor: float, about: QPointF | None = None) -> None:
        about = about or QPointF(self.width() / 2, self.height() / 2)
        model = self.to_model(about)
        self.zoom = min(4.0, max(0.2, self.zoom * factor))
        self.origin = QPointF(about.x() - model[0] * self.scale, about.y() + model[1] * self.scale)
        self.update()

    # -- selection ---------------------------------------------------------------------------

    def clear_selection(self) -> None:
        self.selected_atoms.clear()
        self.selected_bonds.clear()
        self.selected_objects.clear()
        self.update()

    def select_objects(self, objects, add: bool = True) -> None:
        """Select reaction arrows and plus signs (added to the selection by default)."""
        if not add:
            self.clear_selection()
        self.selected_objects |= set(objects)
        self.update()

    def objects_in(self, inside) -> set[int]:
        """Arrows (both ends) and plus signs for which ``inside(point)`` holds."""
        found = {pid for pid, plus in self.sketch.pluses.items() if inside(plus.pos)}
        found |= {
            rid for rid, arrow in self.sketch.arrows.items()
            if inside(arrow.tail) and inside(arrow.head)
        }  # fmt: skip
        return found

    def select_atoms(self, atoms, add: bool = False) -> None:
        if not add:
            self.selected_atoms.clear()
            self.selected_bonds.clear()
            self.selected_objects.clear()
        self.selected_atoms |= self.with_groups(atoms)
        self.selected_bonds |= {
            bid
            for bid, b in self.sketch.bonds.items()
            if b.begin in self.selected_atoms and b.end in self.selected_atoms
        }
        self.update()

    def select_hit(self, hit: Hit, add: bool = False) -> None:
        if hit.kind == "atom":
            self.select_atoms([hit.id], add)
        elif hit.kind == "bond":
            bond = self.sketch.bonds[hit.id]
            self.select_atoms([bond.begin, bond.end], add)
            self.selected_bonds.add(hit.id)
        elif hit.kind in ("arrow", "plus"):
            self.select_objects([hit.id], add)

    def selection_atoms(self) -> set[int]:
        return set(self.selected_atoms)

    # -- clipboard ---------------------------------------------------------------------------

    def copy(self) -> bool:
        """Copy the selection (or the whole drawing): a molfile and SMILES, or for a
        reaction an RXN file and reaction SMILES."""
        selected = self.selected_atoms or self.selected_objects
        part = (
            edits.subset(self.sketch, self.selected_atoms, self.selected_objects)
            if selected
            else self.sketch
        )
        if part.is_empty():
            return False
        mime = QMimeData()
        if part.is_reaction:
            try:
                block = reaction.to_rxn_block(part)
                mime.setData("chemical/x-mdl-rxnfile", block.encode())
                mime.setText(reaction.to_smiles(part))
                QGuiApplication.clipboard().setMimeData(mime)
                return True
            except (reaction.ReactionError, convert.SketchError):
                pass  # not one reaction (two arrows): copy the molecules alone
        if not part.atoms:
            return False
        block = convert.to_molblock(part)
        mime.setData("chemical/x-mdl-molfile", block.encode())
        try:
            mime.setText(convert.to_smiles(part))
        except convert.SketchError:
            mime.setText(block)
        QGuiApplication.clipboard().setMimeData(mime)
        return True

    def cut(self) -> None:
        if (self.selected_atoms or self.selected_objects) and self.copy():
            self.erase_selection("Cut")

    def paste(self) -> bool:
        from chemsuite.app.structure_editor import clipboard_structure_text

        text = clipboard_structure_text()
        incoming = reaction.from_text(text) if reaction.looks_like_reaction(text) else None
        incoming = incoming or convert.from_text(text)
        if incoming is None or incoming.is_empty():
            return False
        if self.mouse is not None and self.rect().contains(self.mouse.toPoint()):
            at = self.to_model(self.mouse)
        elif self.sketch.is_empty():
            at = self.to_model(QPointF(self.width() / 2, self.height() / 2))
        else:
            x0, y0, x1, y1 = self.sketch.bounds()
            px0, _, px1, _ = incoming.bounds()
            at = (x1 + 1.5 + (px1 - px0) / 2, (y0 + y1) / 2)
        placed = convert.placed(incoming, at)
        added: list[int] = []
        before = set(self.sketch.arrows) | set(self.sketch.pluses)
        self.edit("Paste", lambda ch: added.extend(edits.insert(ch, placed)))
        self.select_atoms(added)
        self.select_objects((set(self.sketch.arrows) | set(self.sketch.pluses)) - before)
        return True

    def erase_selection(self, name: str = "Erase") -> None:
        atoms, bonds = set(self.selected_atoms), set(self.selected_bonds)
        objects = set(self.selected_objects)
        if atoms or bonds or objects:
            self.edit(name, lambda ch: edits.erase(ch, atoms=atoms, bonds=bonds, objects=objects))
            self.clear_selection()

    # -- painting ----------------------------------------------------------------------------

    def _atom_colour(self, symbol: str, dark: bool) -> QColor:
        """RDKit's element colour (its dark-mode palette on dark themes); C and H use the text
        colour so they read on any background, R-groups and query atoms the link colour."""
        if symbol in ("C", "H"):
            return self.palette().text().color()
        if symbol in QUERY_ATOMS or symbol in (RGROUP, ATOM_LIST):
            return self.palette().link().color()
        r, g, b = _element_rgb(symbol, dark)
        return QColor.fromRgbF(r, g, b)

    def paintEvent(self, _event) -> None:
        if not self._fitted:
            self.fit()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        dark = is_dark(pal)
        p.fillRect(self.rect(), pal.base())
        ink = pal.text().color()
        highlight = QColor(pal.highlight().color())
        s = self.scale

        self._paint_highlights(p, highlight)

        pen = QPen(ink, max(1.2, 0.04 * s))
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        for line in self.scene.lines:
            pen.setStyle(Qt.PenStyle.DashLine if line.dashed else Qt.PenStyle.SolidLine)
            p.setPen(pen)
            p.drawLine(self.to_screen(line.a), self.to_screen(line.b))
        pen.setStyle(Qt.PenStyle.SolidLine)
        p.setPen(pen)
        p.setBrush(ink)
        for wedge in self.scene.wedges:
            p.drawPolygon(QPolygonF([self.to_screen(wedge.tip), *map(self.to_screen, wedge.base)]))
        p.setBrush(Qt.BrushStyle.NoBrush)
        for _oid, a, b in self.scene.strokes:
            p.drawLine(self.to_screen(a), self.to_screen(b))
        for wavy in [*self.scene.wavy, *self.scene.stubs]:
            path = QPainterPath(self.to_screen(wavy.points[0]))
            for point in wavy.points[1:]:
                path.lineTo(self.to_screen(point))
            p.drawPath(path)
        self._paint_tags(p, pal.link().color(), pal.base().color())
        self._paint_maps(p, pal.link().color())

        font = QFont(self.font())
        font.setPixelSize(max(6, round(0.42 * s)))
        small = QFont(font)
        small.setPixelSize(max(5, round(0.3 * s)))
        problem = QColor("#d32f2f")
        for label in self.scene.labels.values():
            colour = problem if label.problem else self._atom_colour(label.symbol, dark)
            self._paint_label(p, label, font, small, colour, pal.base().color())
        self._paint_stereo(p, font, ink, pal.base().color())
        for aid in self.scene.problems:
            if aid in self.sketch.atoms:
                centre = self.to_screen(self.sketch.atoms[aid].pos)
                p.setPen(QPen(problem, max(1.0, 0.03 * s), Qt.PenStyle.DotLine))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawEllipse(centre, 0.36 * s, 0.36 * s)
        warn = QColor("#e65100")  # amber: a warning, not an error
        for bid in self.wedge_problems:
            bond = self.sketch.bonds.get(bid)
            if bond is None:
                continue
            a, b = self.sketch.atoms[bond.begin].pos, self.sketch.atoms[bond.end].pos
            mid = self.to_screen(((a[0] + b[0]) / 2, (a[1] + b[1]) / 2))
            p.setPen(QPen(warn, max(1.0, 0.03 * s), Qt.PenStyle.DotLine))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawEllipse(mid, 0.42 * s, 0.42 * s)
        self.tool.paint(self, p)
        p.end()

    def _paint_tags(self, p: QPainter, colour: QColor, base: QColor) -> None:
        """Query-bond tags (S/D, Rng …) and attachment-point numbers, small, in the link
        colour that query atoms use."""
        if not self.scene.tags:
            return
        font = QFont(self.font())
        font.setPixelSize(max(5, round(0.26 * self.scale)))
        metrics = QFontMetricsF(font)
        p.setFont(font)
        for at, text in self.scene.tags:
            centre = self.to_screen(at)
            width = metrics.horizontalAdvance(text)
            p.fillRect(QRectF(centre.x() - width / 2 - 1, centre.y() - metrics.ascent() / 2 - 1,
                              width + 2, metrics.height()), base)  # fmt: skip
            p.setPen(colour)
            p.drawText(QPointF(centre.x() - width / 2, centre.y() + metrics.capHeight() / 2), text)

    def _paint_maps(self, p: QPainter, colour: QColor) -> None:
        """Atom-atom map numbers, small, beside their atoms."""
        if not self.scene.maps:
            return
        font = QFont(self.font())
        font.setPixelSize(max(5, round(0.28 * self.scale)))
        metrics = QFontMetricsF(font)
        p.setFont(font)
        p.setPen(colour)
        for _aid, at, text in self.scene.maps:
            centre = self.to_screen(at)
            width = metrics.horizontalAdvance(text)
            p.drawText(QPointF(centre.x() - width / 2, centre.y() + metrics.capHeight() / 2), text)

    def stereo_label_positions(self) -> list[tuple[str, QPointF]]:
        """(text, screen position) of each R/S and E/Z label shown."""
        sketch, out = self.sketch, []
        shown = {aid: f"({code})" for aid, code in self.stereo.atoms.items()}
        for aid, atom in sketch.atoms.items():
            group = enhanced_text(atom.enhanced)
            if group and SketchCanvas.show_stereo:
                shown[aid] = f"{shown[aid]} {group}" if aid in shown else group
        for aid, text in sorted(shown.items()):
            if ("atom", aid) in self.hidden_stereo or aid not in sketch.atoms:
                continue
            here = sketch.atoms[aid].pos
            spots = [sketch.atoms[n].pos for n in sketch.neighbours(aid)]
            out.append((text, self.to_screen(_away(here, spots, 0.62))))
        for bid, code in sorted(self.stereo.bonds.items()):
            if ("bond", bid) in self.hidden_stereo or bid not in sketch.bonds:
                continue
            bond = sketch.bonds[bid]
            a, b = sketch.atoms[bond.begin].pos, sketch.atoms[bond.end].pos
            mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            dx, dy = b[0] - a[0], b[1] - a[1]
            length = math.hypot(dx, dy) or 1.0
            normal = (-dy / length, dx / length)
            others = [
                sketch.atoms[n].pos
                for end in (bond.begin, bond.end)
                for n in sketch.neighbours(end)
                if n not in (bond.begin, bond.end)
            ]
            side = sum((o[0] - mid[0]) * normal[0] + (o[1] - mid[1]) * normal[1] for o in others)
            sign = -1.0 if side > 0 else 1.0  # the side with fewer substituents
            spot = (mid[0] + sign * normal[0] * 0.42, mid[1] + sign * normal[1] * 0.42)
            out.append((f"({code})", self.to_screen(spot)))
        return out

    def _paint_stereo(self, p: QPainter, font: QFont, ink: QColor, base: QColor) -> None:
        labels = self.stereo_label_positions()
        if not labels:
            return
        italic = QFont(font)
        italic.setItalic(True)
        italic.setPixelSize(max(5, round(font.pixelSize() * 0.72)))
        metrics = QFontMetricsF(italic)
        colour = QColor(ink)
        colour.setAlphaF(0.75)
        p.setFont(italic)
        for text, at in labels:
            width = metrics.horizontalAdvance(text)
            box = QRectF(at.x() - width / 2 - 1, at.y() - metrics.ascent() / 2 - 2,
                         width + 2, metrics.height())  # fmt: skip
            p.fillRect(box, base)
            p.setPen(colour)
            p.drawText(QPointF(at.x() - width / 2, at.y() + metrics.capHeight() / 2), text)

    def _paint_highlights(self, p: QPainter, highlight: QColor) -> None:
        s = self.scale
        selected = QColor(highlight)
        selected.setAlpha(110)
        hovered = QColor(highlight)
        hovered.setAlpha(70)
        hover_object = {self.hover.id} if self.hover.kind in ("arrow", "plus") else set()
        for colour, atoms, bonds, objects in (
            (selected, self.selected_atoms, self.selected_bonds, self.selected_objects),
            (
                hovered,
                {self.hover.id} if self.hover.kind == "atom" else set(),
                {self.hover.id} if self.hover.kind == "bond" else set(),
                hover_object,
            ),
        ):
            pen = QPen(colour, 0.3 * s)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            p.setPen(pen)
            for bid in bonds:
                bond = self.sketch.bonds.get(bid)
                if bond is not None:
                    p.drawLine(
                        self.to_screen(self.sketch.atoms[bond.begin].pos),
                        self.to_screen(self.sketch.atoms[bond.end].pos),
                    )
            for _oid, a, b in self.scene.strokes:
                if _oid in objects:
                    p.drawLine(self.to_screen(a), self.to_screen(b))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(colour)
            for aid in atoms:
                atom = self.sketch.atoms.get(aid)
                if atom is not None:
                    p.drawEllipse(self.to_screen(atom.pos), 0.28 * s, 0.28 * s)
            p.setBrush(Qt.BrushStyle.NoBrush)

    def _paint_label(
        self, p: QPainter, label: render.Label, font: QFont, small: QFont, colour, base
    ) -> None:
        big, little = QFontMetricsF(font), QFontMetricsF(small)
        widths = [
            (little if kind != "normal" else big).horizontalAdvance(text)
            for text, kind in label.runs
        ]
        centre = self.to_screen(label.pos)
        left = centre.x() - widths[label.centre] / 2 - sum(widths[: label.centre])
        baseline = centre.y() + big.capHeight() / 2
        box = QRectF(left - 1, baseline - big.ascent(), sum(widths) + 2, big.height())
        hovered = self.hover.kind == "atom" and self.hover.id == label.atom
        if not hovered and label.atom not in self.selected_atoms:
            p.fillRect(box, base)  # keep bond ends from showing through the label
        p.setPen(colour)
        x = left
        for (text, kind), width in zip(label.runs, widths, strict=True):
            if kind == "normal":
                p.setFont(font)
                p.drawText(QPointF(x, baseline), text)
            else:
                p.setFont(small)
                shift = -big.capHeight() * 0.55 if kind == "sup" else big.capHeight() * 0.35
                p.drawText(QPointF(x, baseline + shift), text)
            x += width

    # -- input -------------------------------------------------------------------------------

    def mousePressEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MouseButton.MiddleButton:
            self._pan_from = event.position()
            return
        if event.button() != Qt.MouseButton.LeftButton:
            return
        self.setFocus()
        pos = self.to_model(event.position())
        self.tool.press(self, pos, self.hit_at(pos), event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self.mouse = event.position()
        if event.buttons() & Qt.MouseButton.MiddleButton and hasattr(self, "_pan_from"):
            delta = event.position() - self._pan_from
            self._pan_from = event.position()
            self.origin += delta
            self.update()
            return
        pos = self.to_model(event.position())
        if event.buttons() & Qt.MouseButton.LeftButton:
            self.tool.drag(self, pos, event)
        else:
            hit = self.hit_at(pos)
            if hit != self.hover:
                self.hover = hit
                self.update()

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            pos = self.to_model(event.position())
            self.tool.release(self, pos, event)
            self.hover = self.hit_at(pos)
            self.update()

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:
        pos = self.to_model(event.position())
        self.tool.double_click(self, pos, self.hit_at(pos))

    def leaveEvent(self, _event) -> None:
        self.mouse = None
        if self.hover:
            self.hover = NOTHING
            self.update()

    def wheelEvent(self, event) -> None:
        delta = event.angleDelta()
        if event.modifiers() & (
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
        ):
            self.zoom_by(math.pow(1.0015, delta.y()), event.position())
        else:
            pixel = event.pixelDelta()
            step = pixel if not pixel.isNull() else delta / 4
            self.origin += QPointF(step.x(), step.y())
            self.update()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if not self._fitted:
            self.fit()

    def event(self, event) -> bool:
        if event.type() == QEvent.Type.ToolTip:
            pos = self.to_model(event.position() if hasattr(event, "position") else event.pos())
            hit = self.hit_at(pos)
            if hit.kind == "atom":
                message = self.report.problems.get(hit.id, "")
            else:
                message = self.wedge_problems.get(hit.id, "") if hit.kind == "bond" else ""
            if message:
                QToolTip.showText(event.globalPos(), message, self)
            else:
                QToolTip.hideText()
                event.ignore()
            return True
        return super().event(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if self._command(event) or self._hover_key(event) or self._tool_key(event):
            event.accept()
            return
        super().keyPressEvent(event)

    def _command(self, event: QKeyEvent) -> bool:
        key = event.key()
        mods = event.modifiers()
        mod = mods & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
        shift = bool(mods & Qt.KeyboardModifier.ShiftModifier)
        alt = bool(mods & Qt.KeyboardModifier.AltModifier)
        if event.matches(QKeySequence.StandardKey.Undo):
            self._undo()
        elif event.matches(QKeySequence.StandardKey.Redo) or (mod and key == Qt.Key.Key_Y):
            self._redo()
        elif mod and shift and key == Qt.Key.Key_M:
            part = edits.subset(self.sketch, self.selected_atoms) if self.selected_atoms else None
            block = convert.to_molblock(part or self.sketch)
            QGuiApplication.clipboard().setText(block)
        elif event.matches(QKeySequence.StandardKey.Copy):
            self.copy()
        elif event.matches(QKeySequence.StandardKey.Cut):
            self.cut()
        elif event.matches(QKeySequence.StandardKey.Paste):
            self.paste()
        elif mod and shift and key == Qt.Key.Key_K:
            self.clean()
        elif mod and shift and key == Qt.Key.Key_A:
            self.clear_selection()
        elif event.matches(QKeySequence.StandardKey.SelectAll):
            self.select_atoms(self.sketch.atoms)
        elif alt and key == Qt.Key.Key_E:
            self.popup_enhanced()
        elif alt and key in (Qt.Key.Key_H, Qt.Key.Key_V):
            atoms = self.selected_atoms or set(self.sketch.atoms)
            horizontal = key == Qt.Key.Key_H
            self.edit("Flip", lambda ch: edits.flip(ch, atoms, horizontal))
        elif key == Qt.Key.Key_Escape or (key == Qt.Key.Key_Backtab):
            from chemsuite.app.sketcher.tools import SelectTool

            self.cancel_preview()
            self.set_tool(SelectTool())
        elif key in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace) and not self.hover:
            self.erase_selection()
        else:
            return False
        return True

    def _hover_key(self, event: QKeyEvent) -> bool:
        """Keys over the atom or bond under the pointer (interaction-spec §1)."""
        if not self.hover or self._preview is not None:
            return False
        mods = event.modifiers() & ~Qt.KeyboardModifier.ShiftModifier
        if mods & ~Qt.KeyboardModifier.KeypadModifier:
            return False
        text, key, hit = event.text(), event.key(), self.hover
        if key in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace) and hit.kind in ("arrow", "plus"):
            self.edit("Erase", lambda ch: edits.erase(ch, objects=[hit.id]))
        elif key in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            if hit.kind == "atom":
                atoms = self.with_groups([hit.id])
                self.edit("Erase", lambda ch: edits.erase(ch, atoms=atoms))
            else:
                self.edit("Erase", lambda ch: edits.erase(ch, bonds=[hit.id]))
        elif hit.kind == "atom" and key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self.edit_label(hit.id)
        elif hit.kind == "atom" and text in ELEMENT_KEYS:
            label = ELEMENT_KEYS[text]

            def apply(ch):
                self._drop_group(ch, hit.id)
                edits.set_element(ch, hit.id, label)

            self.edit(f"Set {label}", apply)
        elif hit.kind == "atom" and text in ("1", "2", "3"):
            order = BOND_KEYS[text]
            self.edit("Add bond", lambda ch: edits.bond_from_atom(ch, hit.id, order))
        elif hit.kind == "bond" and text in BOND_KEYS:
            order = BOND_KEYS[text]
            self.edit("Bond order", lambda ch: edits.set_bond_order(ch, hit.id, order))
        elif hit.kind == "atom" and text in ("=", "+", "-", "−"):
            delta = 1 if text in ("=", "+") else -1
            self.edit("Charge", lambda ch: edits.add_charge(ch, hit.id, delta))
        else:
            return False
        if self.mouse is not None:
            self.hover = self.hit_at(self.to_model(self.mouse))
        return True

    def _drop_group(self, change: Change, atom: int) -> None:
        """Typing an element on an abbreviation's label replaces the abbreviation."""
        found = self.contracted_group(atom)
        if found is not None:
            gid, group = found
            change.set_group(gid, None)
            edits.erase(change, atoms=group.atoms - {atom})

    # -- labels and abbreviations ------------------------------------------------------------

    def edit_label(self, atom: int, text: str | None = None) -> None:
        """Type an atom's label: an element (N, Cl, 13C, NH2, O-) or an abbreviation (CO2Me)."""
        from PySide6.QtWidgets import QCompleter, QLineEdit

        from chemsuite.sketch.abbrev import TABLE

        self.cancel_preview()
        if self._label_editor is not None:
            self._label_editor.deleteLater()
        found = self.contracted_group(atom)
        current = found[1].label if found else self.sketch.atoms[atom].label
        editor = QLineEdit(current if text is None else text, self)
        editor.setObjectName("sketch-label-editor")
        names = sorted({e.label for e in TABLE} | {a for e in TABLE for a in e.aliases})
        completer = QCompleter(names, editor)
        completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        editor.setCompleter(completer)
        editor.setFixedWidth(110)
        at = self.to_screen(self.sketch.atoms[atom].pos)
        editor.move(int(at.x() - 20), int(at.y() - editor.sizeHint().height() / 2))
        editor.selectAll()
        editor.returnPressed.connect(lambda: self._apply_label(atom, editor.text()))
        editor.editingFinished.connect(self._close_label_editor)
        self._label_editor = editor
        editor.show()
        editor.setFocus()

    def _apply_label(self, atom: int, text: str) -> None:
        self._close_label_editor()
        if atom not in self.sketch.atoms or not text.strip():
            return
        from chemsuite.sketch.abbrev import set_label

        problem = ""

        def apply(ch):
            nonlocal problem
            self._drop_group(ch, atom)
            problem = set_label(ch, atom, text)

        change = self.history.change(f"Label {text.strip()}")
        apply(change)
        if problem:
            self.history.discard(change)
            self.message.emit(problem)
        elif self.history.commit(change) is None:
            self._recompute()
        self.setFocus()

    def _close_label_editor(self) -> None:
        editor, self._label_editor = self._label_editor, None
        if editor is not None:
            editor.hide()
            editor.deleteLater()

    def _arrow_menu(self, menu, rid: int) -> None:
        from dataclasses import replace

        arrow = self.sketch.arrows[rid]
        kinds = menu.addMenu("Arrow")
        for kind, title in ARROW_TITLES.items():
            action = kinds.addAction(
                title,
                lambda k=kind: self.edit(
                    "Arrow style", lambda ch: ch.set_arrow(rid, replace(arrow, kind=k))
                ),
            )
            action.setCheckable(True)
            action.setChecked(kind is arrow.kind)
        menu.addAction(
            "Reverse Arrow",
            lambda: self.edit(
                "Reverse arrow",
                lambda ch: ch.set_arrow(
                    rid, replace(arrow, x0=arrow.x1, y0=arrow.y1, x1=arrow.x0, y1=arrow.y0)
                ),
            ),
        )
        menu.addAction(
            "Delete Arrow", lambda: self.edit("Erase", lambda ch: edits.erase(ch, objects=[rid]))
        )

    def contextMenuEvent(self, event) -> None:
        hit = self.hit_at(self.to_model(QPointF(event.pos())))
        self.context_menu(hit).exec(event.globalPos())

    def context_menu(self, hit: Hit):
        """The right-click menu for what is under the pointer (an atom, a bond or nothing)."""
        from PySide6.QtWidgets import QMenu

        from chemsuite.sketch import abbrev

        found = self.sketch.group_of(hit.id) if hit.kind == "atom" else None
        menu = QMenu(self)
        if hit.kind == "atom":
            menu.addAction("Edit Label…", lambda: self.edit_label(hit.id))
            if found is None:
                self._query_menus(menu, hit.id)
                self._attach_menu(menu, hit.id)
                self._hydrogen_actions(menu, hit.id)
                self._centre_menu(menu, hit.id)
        if hit.kind == "bond":
            self._double_bond_menu(menu, hit.id)
            self._bond_query_menu(menu, hit.id)
        if hit.kind == "arrow":
            self._arrow_menu(menu, hit.id)
        if hit.kind == "atom" and self.sketch.atoms[hit.id].mapping:
            menu.addAction(
                "Clear Map Number",
                lambda: self.edit("Clear map", lambda ch: edits.set_mapping(ch, [hit.id], 0)),
            )
        if any(a.mapping for a in self.sketch.atoms.values()):
            mapped = [aid for aid, a in self.sketch.atoms.items() if a.mapping]
            menu.addAction(
                "Clear All Mapping",
                lambda: self.edit("Clear mapping", lambda ch: edits.set_mapping(ch, mapped, 0)),
            )
        if found is not None:
            gid, group = found
            if group.expanded:
                menu.addAction(
                    f"Contract to {group.label}",
                    lambda: self.edit("Contract", lambda ch: abbrev.set_expanded(ch, gid, False)),
                )
            else:
                menu.addAction(
                    f"Expand {group.label}",
                    lambda: self.edit("Expand", lambda ch: abbrev.set_expanded(ch, gid, True)),
                )
            menu.addAction(
                f"Remove the {group.label} label (keep the atoms)",
                lambda: self.edit("Ungroup", lambda ch: abbrev.ungroup(ch, gid)),
            )
        if not menu.isEmpty():
            menu.addSeparator()
        history = self.history
        menu.addAction("Undo", self._undo).setEnabled(history.can_undo)
        menu.addAction("Redo", self._redo).setEnabled(history.can_redo)
        menu.addSeparator()
        has_atoms = not self.sketch.is_empty()
        menu.addAction("Cut", self.cut).setEnabled(bool(self.selected_atoms))
        menu.addAction("Copy", self.copy).setEnabled(has_atoms)
        menu.addAction("Paste", self.paste)
        menu.addAction("Select All", lambda: self.select_atoms(self.sketch.atoms)).setEnabled(
            has_atoms
        )
        menu.addSeparator()
        self._salts_menu(menu.addMenu("Add Salt or Solvent"))
        menu.addSeparator()
        labels = menu.addAction("Show Stereo Labels (R/S, E/Z)", self.toggle_stereo_labels)
        labels.setCheckable(True)
        labels.setChecked(SketchCanvas.show_stereo)
        menu.addAction("Show All Hydrogens", lambda: self.show_hydrogens(None)).setEnabled(
            any(self.report.hydrogens.values())
        )
        menu.addAction("Hide All Hydrogens", lambda: self.hide_hydrogens(None)).setEnabled(
            any(a.symbol == "H" for a in self.sketch.atoms.values())
        )
        menu.addSeparator()
        what = "Clean Selection" if self.selected_atoms else "Clean Structure"
        clean = menu.addAction(f"{what}\t{_shortcut('Ctrl+Shift+K')}", self.clean)
        clean.setEnabled(has_atoms and self.report.valid)
        menu.addAction("Clear", lambda: self.set_sketch(Sketch(), "Clear")).setEnabled(has_atoms)
        return menu

    def _salts_menu(self, menu) -> None:
        from chemsuite.sketch import salts

        for group in salts.GROUPS:
            sub = menu.addMenu(group)
            for entry in salts.TABLE:
                if entry.group == group:
                    sub.addAction(entry.title, lambda e=entry: self.add_salt(e))

    def add_salt(self, entry) -> None:
        """Drop a salt or solvent to the right of the drawing and select it."""
        from chemsuite.sketch import salts

        added: list[int] = []
        self.edit(f"Add {entry.name}", lambda ch: added.extend(salts.place(ch, entry)))
        self.select_atoms(added)
        self._ensure_visible()

    def _ensure_visible(self) -> None:
        """Refit when the drawing has grown past the window's edges."""
        if self.sketch.is_empty():
            return
        x0, y0, x1, y1 = self.sketch.bounds()
        a, b = self.to_screen((x0, y1)), self.to_screen((x1, y0))
        if not self.rect().adjusted(8, 8, -8, -8).contains(QRectF(a, b).toRect()):
            self.fit()

    def clean(self) -> None:
        """Clean Structure (ChemDraw's Mod+Shift+K): a fresh layout for the fragments of the
        selection, or for the whole drawing, keeping place, turn and stereo (``layout.py``)."""
        from chemsuite.sketch import layout

        atoms = set(self.selected_atoms) or None
        name = "Clean Selection" if atoms else "Clean Structure"
        self._set_stereo(name, lambda ch: layout.clean(ch, atoms))

    # -- stereo and hydrogens ------------------------------------------------------------------
    def toggle_stereo_labels(self) -> None:
        SketchCanvas.show_stereo = not SketchCanvas.show_stereo
        self.hidden_stereo.clear()
        self._recompute()

    def _targets(self, atom: int) -> list[int]:
        """The selected atoms when *atom* is one of them, else just *atom*."""
        return sorted(self.selected_atoms) if atom in self.selected_atoms else [atom]

    def show_hydrogens(self, atom: int | None) -> None:
        atoms = list(self.sketch.atoms) if atom is None else self._targets(atom)
        counts = dict(self.report.hydrogens)
        self.edit("Show Hydrogens", lambda ch: stereo.show_hydrogens(ch, atoms, counts))

    def hide_hydrogens(self, atom: int | None) -> None:
        atoms = list(self.sketch.atoms) if atom is None else self._targets(atom)
        kept: list[int] = []

        def apply(ch) -> None:
            kept.append(stereo.hide_hydrogens(ch, atoms)[1])

        self.edit("Hide Hydrogens", apply)
        if kept and kept[0]:
            self.notify(
                f"{kept[0]} hydrogen{'s' * (kept[0] != 1)} on wedges kept: they define "
                "stereocentres."
            )

    def _hydrogen_actions(self, menu, atom: int) -> None:
        atoms = self._targets(atom)
        implicit = sum(self.report.hydrogens.get(a, 0) for a in atoms)
        drawn = sum(len(edits.drawn_hydrogens(self.sketch, a)[0]) for a in atoms)
        show = menu.addAction("Show Explicit Hydrogens", lambda: self.show_hydrogens(atom))
        show.setEnabled(implicit > 0)
        hide = menu.addAction("Hide Explicit Hydrogens", lambda: self.hide_hydrogens(atom))
        hide.setEnabled(drawn > 0)

    def _set_stereo(self, name: str, apply) -> None:
        problem: list[str] = []
        change = self.history.change(name)
        problem.append(apply(change))
        if problem[0]:
            self.history.discard(change)
            self._recompute()
            self.notify(problem[0])
            return
        if self.history.commit(change) is None:
            self._recompute()

    def _centre_menu(self, menu, atom: int) -> None:
        code = self.stereo.atoms.get(atom)
        if code is None and atom not in self.stereo.undefined_atoms:
            return
        sub = menu.addMenu("Stereocentre (R/S)")
        for want in ("R", "S"):
            action = sub.addAction(
                want,
                lambda w=want: self._set_stereo(
                    f"Make {w}", lambda ch: stereo.set_centre(ch, atom, w)
                ),
            )
            action.setCheckable(True)
            action.setChecked(code is not None and code.upper() == want)
        sub.addSeparator()
        shown = sub.addAction("Show Label", lambda: self._toggle_label("atom", atom))
        shown.setCheckable(True)
        shown.setChecked(code is not None and ("atom", atom) not in self.hidden_stereo)
        shown.setEnabled(code is not None)
        self.enhanced_menu(menu.addMenu("Enhanced Stereo"), self._centres(self._targets(atom)))

    def popup_enhanced(self):
        """Alt+E (Ketcher): the Enhanced Stereo menu for the selected stereocentres, or all."""
        from PySide6.QtGui import QCursor
        from PySide6.QtWidgets import QMenu

        atoms = self._centres(sorted(self.selected_atoms) or sorted(self.sketch.atoms))
        if not atoms:
            self.notify("Select one or more stereocentres for enhanced stereo.")
            return None
        menu = QMenu(self)
        self.enhanced_menu(menu, atoms)
        menu.popup(QCursor.pos())
        return menu

    def _centres(self, atoms) -> list[int]:
        """The stereocentres among *atoms* (labelled or not yet defined)."""
        known = set(self.stereo.atoms) | set(self.stereo.undefined_atoms)
        return [a for a in atoms if a in known]

    def enhanced_menu(self, menu, atoms: list[int]) -> None:
        """Enhanced stereo for *atoms* (V3000 ABS / AND / OR groups): absolute, "and" (a
        racemic mixture of this and its mirror image) or "or" (one of the two, unknown)."""
        values = {self.sketch.atoms[a].enhanced for a in atoms if a in self.sketch.atoms}
        current = values.pop() if len(values) == 1 else None

        def setter(value: str, name: str):
            return lambda: self.edit(name, lambda ch: edits.set_enhanced(ch, atoms, value))

        plain = menu.addAction("Absolute (as drawn)", setter("", "Enhanced stereo: absolute"))
        plain.setCheckable(True)
        plain.setChecked(current == "")
        tagged = menu.addAction("Absolute, labelled abs", setter("abs", "Enhanced stereo: abs"))
        tagged.setCheckable(True)
        tagged.setChecked(current == "abs")
        for kind, title, meaning in (
            ("and", "AND", "racemic: this and its mirror image"),
            ("or", "OR", "one of this or its mirror image, not known which"),
        ):
            sub = menu.addMenu(f"{title} Group ({meaning})")
            used = sorted(
                {p[1] for a in self.sketch.atoms.values()
                 if (p := parse_enhanced(a.enhanced)) and p[0] == kind}
            )  # fmt: skip
            for n in used:
                value = f"{kind}{n}"
                action = sub.addAction(
                    enhanced_text(value), setter(value, f"Enhanced stereo: {value}")
                )
                action.setCheckable(True)
                action.setChecked(current == value)
            new = f"{kind}{next_enhanced(self.sketch, kind)}"
            sub.addAction(
                f"New Group ({enhanced_text(new)})", setter(new, f"Enhanced stereo: {new}")
            )
        menu.setEnabled(bool(atoms))

    def _double_bond_menu(self, menu, bond: int) -> None:
        code = self.stereo.bonds.get(bond)
        unknown = bond in self.stereo.unknown_bonds
        if code is None and not unknown and bond not in self.stereo.undefined_bonds:
            return
        sub = menu.addMenu("Double Bond (E/Z)")
        for want, text in (("E", "E"), ("Z", "Z"), (stereo.UNKNOWN, "Unknown (crossed bond)")):
            action = sub.addAction(
                text,
                lambda w=want, t=text: self._set_stereo(
                    "Make E/Z unknown" if w == stereo.UNKNOWN else f"Make {t}",
                    lambda ch: stereo.set_double_bond(ch, bond, w),
                ),
            )
            action.setCheckable(True)
            action.setChecked(code == want or (unknown and want == stereo.UNKNOWN))
        sub.addSeparator()
        shown = sub.addAction("Show Label", lambda: self._toggle_label("bond", bond))
        shown.setCheckable(True)
        shown.setChecked(code is not None and ("bond", bond) not in self.hidden_stereo)
        shown.setEnabled(code is not None)
        if not menu.isEmpty():
            menu.addSeparator()

    def _bond_query_menu(self, menu, bond: int) -> None:
        """Bond Query: CTfile query orders and ring/chain topology, for structure searches."""
        here = self.sketch.bonds[bond]
        sub = menu.addMenu("Bond Query")
        for order, text in BOND_QUERIES.items():
            action = sub.addAction(
                text,
                lambda o=order, t=text: self.edit(
                    f"Bond: {t.lower()}", lambda ch: edits.set_bond_query(ch, bond, o)
                ),
            )
            action.setCheckable(True)
            action.setChecked(here.order is order)
        if here.order.is_query or here.order is Order.AROMATIC:
            sub.addAction(
                "Back to Single",
                lambda: self.edit(
                    "Bond: single", lambda ch: edits.set_bond_query(ch, bond, Order.SINGLE)
                ),  # fmt: skip
            )
        sub.addSeparator()
        for value, text in TOPOLOGY.items():
            action = sub.addAction(
                text,
                lambda v=value, t=text: self.edit(
                    f"Bond: {t.lower()}", lambda ch: edits.set_bond_query(ch, bond, topology=v)
                ),
            )
            action.setCheckable(True)
            action.setChecked(here.topology == value)

    def _attach_menu(self, menu, atom: int) -> None:
        """Attachment Point for R-group members (CTfile ``M  APO``)."""
        atoms = self._targets(atom)
        now = {self.sketch.atoms[a].attach for a in atoms}
        sub = menu.addMenu("Attachment Point")
        for value, text in ATTACH_POINTS.items():
            action = sub.addAction(
                text,
                lambda v=value, t=text: self.edit(
                    f"Attachment point: {t.lower()}", lambda ch: edits.set_attach(ch, atoms, v)
                ),
            )
            action.setCheckable(True)
            action.setChecked(now == {value})

    def _toggle_label(self, kind: str, id_: int) -> None:
        key = (kind, id_)
        if key in self.hidden_stereo:
            self.hidden_stereo.discard(key)
            if not SketchCanvas.show_stereo:
                SketchCanvas.show_stereo = True
                self._recompute()
        else:
            self.hidden_stereo.add(key)
        self.update()

    def _query_menus(self, menu, atom: int) -> None:
        """R-Group and Query Atom submenus for an atom."""
        here = self.sketch.atoms[atom]

        def setter(label: str):
            return lambda: self.edit(f"Set {label}", lambda ch: edits.set_element(ch, atom, label))

        rgroups = menu.addMenu("R-Group")
        top = min(next_rgroup(self.sketch, skip=atom), MAX_RGROUP)
        for n in range(1, top + 1):
            action = rgroups.addAction(f"R{n}", setter(f"R{n}"))
            action.setCheckable(True)
            action.setChecked(here.symbol == RGROUP and here.rgroup == n)
        queries = menu.addMenu("Query Atom")
        for symbol, meaning in QUERY_ATOMS.items():
            action = queries.addAction(f"{symbol}\t{meaning}", setter(symbol))
            action.setCheckable(True)
            action.setChecked(here.symbol == symbol)
        queries.addSeparator()
        queries.addAction("Atom List…", lambda: self.edit_label(atom, "[C,N,O]"))
        if here.is_query:
            queries.addAction("Back to Carbon", setter("C"))

    def _undo(self) -> None:
        self.cancel_preview()
        self.history.undo()

    def _redo(self) -> None:
        self.cancel_preview()
        self.history.redo()

    def _tool_key(self, event: QKeyEvent) -> bool:
        """Without hovering, the same keys pick the matching tool (as in Ketcher)."""
        from chemsuite.app.sketcher import tools

        if event.modifiers() & ~(
            Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.KeypadModifier
        ):
            return False
        text = event.text()
        if text in ELEMENT_KEYS:
            self.set_tool(tools.AtomTool(ELEMENT_KEYS[text]))
        elif text in ("1", "2", "3"):
            self.set_tool(tools.BondTool(BOND_KEYS[text]))
        elif text in ("=", "+"):
            self.set_tool(tools.ChargeTool(+1))
        elif text in ("-", "−"):
            self.set_tool(tools.ChargeTool(-1))
        elif event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.set_tool(tools.EraseTool())
        elif text == "t":
            self.set_tool(tools.TemplateTool(_next_ring(self.tool)))
        else:
            return False
        return True


def _shortcut(portable: str) -> str:
    """A shortcut as the platform writes it (⇧⌘K on macOS)."""
    return QKeySequence(portable).toString(QKeySequence.SequenceFormat.NativeText)


def _away(here, spots, distance: float):
    """A point *distance* from *here* in the middle of the widest gap between its bonds to
    *spots* (below-right when it has none)."""
    if not spots:
        return (here[0] + distance * 0.7, here[1] - distance * 0.7)
    angles = sorted(math.atan2(y - here[1], x - here[0]) for x, y in spots)
    if len(angles) == 1:
        best = angles[0] + math.pi
    else:
        gaps = [
            ((angles[(k + 1) % len(angles)] - a) % (2 * math.pi) or 2 * math.pi, a)
            for k, a in enumerate(angles)
        ]
        width, start = max(gaps)
        best = start + width / 2
    return (here[0] + math.cos(best) * distance, here[1] + math.sin(best) * distance)
