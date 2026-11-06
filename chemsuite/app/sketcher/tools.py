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

"""The sketcher's tools: each turns mouse input into edits from ``chemsuite.sketch.edits``.

Mouse behaviour follows docs/sketcher/interaction-spec.md §3, which cites the Ketcher tool
each rule comes from. A drag shows its result as a preview edit, which is committed on
release (one undo step) or discarded on Escape.
"""

from __future__ import annotations

import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from chemsuite.sketch import edits, geometry
from chemsuite.sketch.model import MAX_RGROUP, RGROUP, Order, Stereo, next_rgroup

Point = tuple[float, float]


def _free_angle(event) -> bool:
    """Ctrl (⌘ on macOS) frees the angle from the 15 degree steps."""
    return bool(
        event.modifiers() & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
    )


class Tool:
    name = "tool"
    cursor = Qt.CursorShape.CrossCursor

    def __init__(self):
        self.start: Point | None = None
        self.start_hit = None
        self.dragging = False
        self._press_screen: QPointF | None = None

    def press(self, canvas, pos: Point, hit, event) -> None:
        self.start, self.start_hit, self.dragging = pos, hit, False
        self._press_screen = event.position()

    def drag(self, canvas, pos: Point, event) -> None:
        if self.start is None:
            return
        if not self.dragging:
            moved = event.position() - self._press_screen
            if math.hypot(moved.x(), moved.y()) < 4.0:
                return
            self.dragging = True
        self.dragged(canvas, pos, event)

    def release(self, canvas, pos: Point, event) -> None:
        if self.start is None:
            return
        if self.dragging:
            self.finish_drag(canvas, pos, event)
        else:
            self.click(canvas, self.start, self.start_hit, event)
        self.start, self.start_hit, self.dragging = None, None, False

    def double_click(self, canvas, pos: Point, hit) -> None:
        pass

    # Subclasses override these.
    def click(self, canvas, pos: Point, hit, event) -> None:
        pass

    def dragged(self, canvas, pos: Point, event) -> None:
        pass

    def finish_drag(self, canvas, pos: Point, event) -> None:
        canvas.commit_preview()

    def paint(self, canvas, painter: QPainter) -> None:
        pass


#: The rotate handle sits this many pixels above the selection's box.
HANDLE_GAP = 22.0
HANDLE_RADIUS = 6.0


def inside_polygon(point: Point, polygon: list[Point]) -> bool:
    """Even-odd test: is *point* inside the closed *polygon*?"""
    x, y = point
    inside = False
    j = len(polygon) - 1
    for i, (xi, yi) in enumerate(polygon):
        xj, yj = polygon[j]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi) + xi:
            inside = not inside
        j = i
    return inside


class SelectTool(Tool):
    """Click to select, Shift-click to add, drag empty canvas for a box (or a lasso), drag
    the selection to move it, drag the handle above it to rotate it.

    Double-click selects a whole fragment. Atoms dropped onto other atoms merge with them
    (Ketcher ``atomMerge``). Rotation snaps to 15 degree steps; Ctrl (⌘) frees it (Ketcher
    ``tool/rotate.ts``).
    """

    cursor = Qt.CursorShape.ArrowCursor

    def __init__(self, lasso: bool = False):
        super().__init__()
        self.lasso = lasso
        self.name = "lasso" if lasso else "select"
        self.moving: set[int] = set()
        self.box_end: Point | None = None
        self.path: list[Point] = []
        self.rotating: tuple[set[int], Point] | None = None
        self.angle = 0.0
        self.tip: Point | None = None

    @staticmethod
    def handle(canvas) -> tuple[QPointF, Point] | None:
        """The rotate handle (screen point) and the centre it turns about (model point), when
        two or more atoms are selected."""
        atoms = [canvas.sketch.atoms[a] for a in canvas.selected_atoms if a in canvas.sketch.atoms]
        if len(atoms) < 2:
            return None
        xs, ys = [a.x for a in atoms], [a.y for a in atoms]
        centre = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
        top = canvas.to_screen((centre[0], max(ys)))
        return QPointF(top.x(), top.y() - 0.35 * canvas.scale - HANDLE_GAP), centre

    def press(self, canvas, pos, hit, event) -> None:
        super().press(canvas, pos, hit, event)
        self.rotating, self.path, self.angle = None, [], 0.0
        found = self.handle(canvas)
        if found is not None:
            at, centre = found
            gap = event.position() - at
            if math.hypot(gap.x(), gap.y()) <= HANDLE_RADIUS + 4:
                self.rotating = (set(canvas.selected_atoms), centre)
                self.moving = set()
                return
        add = bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
        self.moving = set()
        already = (hit.kind == "atom" and hit.id in canvas.selected_atoms) or (
            hit.kind == "bond" and hit.id in canvas.selected_bonds
        )
        if already:
            self.moving = set(canvas.selected_atoms)
        elif hit:
            canvas.select_hit(hit, add)
            self.moving = set(canvas.selected_atoms)
        elif not add:
            canvas.clear_selection()

    def _rotation(self, pos, event) -> float:
        _, centre = self.rotating
        start = geometry.angle_of(centre, self.start)
        angle = geometry.angle_of(centre, pos) - start
        return angle if _free_angle(event) else geometry.snap_angle(angle)

    def dragged(self, canvas, pos, event) -> None:
        if self.rotating is not None:
            atoms, centre = self.rotating
            angle = self.angle = self._rotation(pos, event)
            self.box_end = None
            self.tip = pos
            canvas.preview("Rotate", lambda ch: edits.rotate(ch, atoms, centre, angle))
        elif self.moving:
            dx, dy = pos[0] - self.start[0], pos[1] - self.start[1]
            atoms = set(self.moving)
            canvas.preview("Move", lambda ch: edits.move(ch, atoms, dx, dy, merge=False))
        else:
            add = bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
            if self.lasso:
                if not self.path:
                    self.path = [self.start]
                self.path.append(pos)
                inside = [
                    aid
                    for aid, a in canvas.sketch.atoms.items()
                    if inside_polygon(a.pos, self.path)
                ]
            else:
                self.box_end = pos
                x0, x1 = sorted((self.start[0], pos[0]))
                y0, y1 = sorted((self.start[1], pos[1]))
                inside = [
                    aid
                    for aid, a in canvas.sketch.atoms.items()
                    if x0 <= a.x <= x1 and y0 <= a.y <= y1
                ]
            canvas.select_atoms(inside, add)

    def finish_drag(self, canvas, pos, event) -> None:
        if self.rotating is not None:
            atoms, centre = self.rotating
            angle = self._rotation(pos, event)
            canvas.cancel_preview()
            canvas.edit("Rotate", lambda ch: edits.rotate(ch, atoms, centre, angle))
            canvas.select_atoms(atoms & set(canvas.sketch.atoms))
            self.rotating, self.tip = None, None
        elif self.moving:
            dx, dy = pos[0] - self.start[0], pos[1] - self.start[1]
            atoms = set(self.moving)
            canvas.cancel_preview()
            canvas.edit("Move", lambda ch: edits.move(ch, atoms, dx, dy, merge=True))
            canvas.select_atoms(atoms & set(canvas.sketch.atoms))
        self.box_end = None
        self.path = []
        canvas.update()

    def double_click(self, canvas, pos, hit) -> None:
        if hit.kind == "atom":
            canvas.edit_label(hit.id)
        elif hit.kind == "bond":
            canvas.select_atoms(canvas.sketch.fragment(canvas.sketch.bonds[hit.id].begin))

    def paint(self, canvas, painter) -> None:
        colour = QColor(canvas.palette().highlight().color())
        fill = QColor(colour)
        fill.setAlpha(30)
        if self.box_end is not None and self.start is not None:
            painter.setPen(QPen(colour, 1, Qt.PenStyle.DashLine))
            painter.setBrush(fill)
            painter.drawRect(QRectF(canvas.to_screen(self.start), canvas.to_screen(self.box_end)))
        if len(self.path) > 1:
            painter.setPen(QPen(colour, 1, Qt.PenStyle.DashLine))
            painter.setBrush(fill)
            painter.drawPolygon(QPolygonF([canvas.to_screen(p) for p in self.path]))
        if self.rotating is not None and self.dragging and self.tip is not None:
            _, centre = self.rotating  # while turning: a spoke from the centre to the pointer
            anchor, tip = canvas.to_screen(centre), canvas.to_screen(self.tip)
            painter.setPen(QPen(colour, 1.4, Qt.PenStyle.DashLine))
            painter.drawLine(anchor, tip)
            painter.setPen(QPen(colour, 1.4))
            painter.setBrush(canvas.palette().base())
            painter.drawEllipse(anchor, 2.5, 2.5)
            painter.drawEllipse(tip, HANDLE_RADIUS, HANDLE_RADIUS)
            painter.setPen(canvas.palette().placeholderText().color())
            painter.drawText(tip + QPointF(10, -8), f"{math.degrees(self.angle):+.0f}°")
            return
        found = self.handle(canvas)
        if found is None or self.dragging:
            return
        at, _centre = found
        painter.setPen(QPen(colour, 1.4))
        painter.drawLine(QPointF(at.x(), at.y() + HANDLE_RADIUS), QPointF(at.x(), at.y() + 12))
        painter.setBrush(canvas.palette().base())
        painter.drawEllipse(at, HANDLE_RADIUS, HANDLE_RADIUS)


class BondTool(Tool):
    """Bonds of one order and wedge style (interaction-spec §3, Ketcher ``tool/bond.ts``)."""

    def __init__(self, order: Order = Order.SINGLE, stereo: Stereo = Stereo.NONE):
        super().__init__()
        self.order, self.stereo = order, stereo
        self.name = f"bond:{order.name.lower()}:{stereo.name.lower()}"

    def click(self, canvas, pos, hit, event) -> None:
        order, stereo = self.order, self.stereo
        if hit.kind == "atom":
            canvas.edit("Add bond", lambda ch: edits.bond_from_atom(ch, hit.id, order, stereo))
        elif hit.kind == "bond":
            canvas.edit("Change bond", lambda ch: edits.click_bond(ch, hit.id, order, stereo))
        else:
            canvas.edit("Add bond", lambda ch: edits.bond_on_canvas(ch, pos, order, stereo))

    def dragged(self, canvas, pos, event) -> None:
        if self.start_hit.kind == "bond":
            return
        order, stereo = self.order, self.stereo
        start_atom = self.start_hit.id if self.start_hit.kind == "atom" else None
        origin = canvas.sketch.atoms[start_atom].pos if start_atom is not None else self.start
        angle = geometry.angle_of(origin, pos)
        if not _free_angle(event):
            angle = geometry.snap_angle(angle)
        end = geometry.along(origin, angle)
        exclude = frozenset({start_atom}) if start_atom is not None else frozenset()

        def apply(ch):
            first = start_atom
            if first is None:
                first = edits.atom_on_canvas(ch, origin)
            target = canvas.sketch.atom_at(*pos, within=0.3, exclude=exclude | {first})
            point = canvas.sketch.atoms[target].pos if target is not None else end
            edits.bond_to_point(ch, first, point, order, stereo, snap_to=0.01)

        canvas.preview("Add bond", apply)


class ChainTool(Tool):
    """Drag a zig-zag carbon chain, one bond per bond length (Ketcher ``tool/chain.ts``)."""

    name = "chain"

    def __init__(self):
        super().__init__()
        self.count = 0
        self.tip: Point | None = None

    def click(self, canvas, pos, hit, event) -> None:
        if hit.kind == "atom":
            canvas.edit("Add bond", lambda ch: edits.bond_from_atom(ch, hit.id))
        elif not hit:
            canvas.edit("Add bond", lambda ch: edits.bond_on_canvas(ch, pos))

    def dragged(self, canvas, pos, event) -> None:
        start_atom = self.start_hit.id if self.start_hit.kind == "atom" else None
        origin = canvas.sketch.atoms[start_atom].pos if start_atom is not None else self.start
        angle = geometry.angle_of(origin, pos)
        if not _free_angle(event):
            angle = geometry.snap_angle(angle)
        self.count = geometry.chain_length(math.dist(origin, pos))
        start = start_atom if start_atom is not None else origin
        count = self.count
        canvas.preview("Chain", lambda ch: edits.chain(ch, start, angle, count))
        self.tip = pos

    def finish_drag(self, canvas, pos, event) -> None:
        canvas.commit_preview()
        self.tip = None

    def paint(self, canvas, painter) -> None:
        if self.tip is None or not self.dragging:
            return
        painter.setPen(canvas.palette().placeholderText().color())
        at = canvas.to_screen(self.tip) + QPointF(12, -12)
        painter.drawText(at, f"{self.count} bond{'s' if self.count != 1 else ''}")


class AtomTool(Tool):
    """Set an atom's element, place a lone atom, or drag a bond to a new atom of this element."""

    def __init__(self, label: str = "C"):
        super().__init__()
        self.label = label
        self.name = f"atom:{label}"

    def resolved(self, canvas, atom: int | None = None) -> str:
        """The label to apply: the R tool numbers each new R-group with the next free number."""
        if self.label != "R":
            return self.label
        if atom is not None and canvas.sketch.atoms[atom].symbol == RGROUP:
            return canvas.sketch.atoms[atom].label  # already an R-group: keep its number
        return f"R{min(next_rgroup(canvas.sketch, skip=atom), MAX_RGROUP)}"

    def click(self, canvas, pos, hit, event) -> None:
        if hit.kind == "atom":
            label = self.resolved(canvas, hit.id)
            canvas.edit(f"Set {label}", lambda ch: edits.set_element(ch, hit.id, label))
        elif not hit:
            label = self.resolved(canvas)
            canvas.edit(f"Add {label}", lambda ch: edits.atom_on_canvas(ch, pos, label))

    def double_click(self, canvas, pos, hit) -> None:
        if hit.kind == "atom":
            canvas.edit_label(hit.id)

    def dragged(self, canvas, pos, event) -> None:
        if self.start_hit.kind != "atom":
            return
        start = self.start_hit.id
        origin = canvas.sketch.atoms[start].pos
        angle = geometry.angle_of(origin, pos)
        if not _free_angle(event):
            angle = geometry.snap_angle(angle)
        end = geometry.along(origin, angle)
        label = self.resolved(canvas)

        def apply(ch):
            bid = edits.bond_to_point(ch, start, end)
            edits.set_element(ch, ch.sketch.bonds[bid].other(start), label)

        canvas.preview(f"Add {label}", apply)


class ChargeTool(Tool):
    def __init__(self, delta: int = 1):
        super().__init__()
        self.delta = delta
        self.name = "charge:plus" if delta > 0 else "charge:minus"

    def click(self, canvas, pos, hit, event) -> None:
        if hit.kind == "atom":
            delta = self.delta
            canvas.edit("Charge", lambda ch: edits.add_charge(ch, hit.id, delta))


class EraseTool(Tool):
    """Click or drag over atoms and bonds to delete them (one undo step per stroke)."""

    name = "erase"

    def __init__(self):
        super().__init__()
        self.atoms: set[int] = set()
        self.bonds: set[int] = set()

    def press(self, canvas, pos, hit, event) -> None:
        super().press(canvas, pos, hit, event)
        self.atoms, self.bonds = set(), set()
        self._collect(canvas, hit)

    def drag(self, canvas, pos, event) -> None:
        self.dragging = True
        self._collect(canvas, canvas.hit_at(pos, frozenset()))

    def _collect(self, canvas, hit) -> None:
        if hit.kind == "atom":
            self.atoms |= canvas.with_groups([hit.id])
        elif hit.kind == "bond":
            self.bonds.add(hit.id)
        else:
            return
        atoms, bonds = set(self.atoms), set(self.bonds)
        canvas.preview("Erase", lambda ch: edits.erase(ch, atoms=atoms, bonds=bonds))

    def release(self, canvas, pos, event) -> None:
        canvas.commit_preview()
        self.start, self.start_hit, self.dragging = None, None, False


class TemplateTool(Tool):
    """Rings and templates (interaction-spec §3, Ketcher ``tool/template.ts``).

    Click empty canvas to place, an atom to attach there, a bond to fuse. Dragging rotates on
    empty canvas, aims from an atom (past one bond length it attaches through a new bond), and
    picks the side on a bond.
    """

    def __init__(self, template=None):
        super().__init__()
        from chemsuite.sketch.templates import BENZENE

        self.template = template or BENZENE
        self.name = f"template:{self.template.name}"

    def click(self, canvas, pos, hit, event) -> None:
        from chemsuite.sketch import templates

        t = self.template
        if hit.kind == "atom":
            canvas.edit(f"Add {t.name}", lambda ch: templates.on_atom(ch, t, hit.id))
        elif hit.kind == "bond":
            canvas.edit(f"Fuse {t.name}", lambda ch: templates.on_bond(ch, t, hit.id))
        else:
            canvas.edit(f"Add {t.name}", lambda ch: templates.on_canvas(ch, t, pos))

    def dragged(self, canvas, pos, event) -> None:
        from chemsuite.sketch import templates

        t, hit = self.template, self.start_hit
        if hit.kind == "bond":
            bond = canvas.sketch.bonds[hit.id]
            a, b = canvas.sketch.atoms[bond.begin].pos, canvas.sketch.atoms[bond.end].pos
            side = 1 if geometry.side_of(a, b, pos) > 0 else -1
            canvas.preview(f"Fuse {t.name}", lambda ch: templates.on_bond(ch, t, hit.id, side))
            return
        origin = canvas.sketch.atoms[hit.id].pos if hit.kind == "atom" else self.start
        angle = geometry.angle_of(origin, pos)
        if not _free_angle(event):
            angle = geometry.snap_angle(angle)
        if hit.kind == "atom":
            extra = math.dist(origin, pos) > 1.0
            canvas.preview(
                f"Add {t.name}", lambda ch: templates.on_atom(ch, t, hit.id, angle, extra)
            )
        else:
            start = self.start
            canvas.preview(f"Add {t.name}", lambda ch: templates.on_canvas(ch, t, start, angle))
