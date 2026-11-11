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

"""The Routes pane: the retrosynthetic routes found for each target, drawn as a tree.

The target sits on the right and its precursors branch to the left, the way a retrosynthetic
scheme is read. A molecule's card is framed green when it is in stock and amber when it still
needs a route of its own; each branch carries a numbered step button, and clicking one shows
that reaction's policy, probability, template and other metadata underneath.

The whole pane is drawn on one :class:`QGraphicsScene` (cards, lines and step buttons as
items), so panning and zooming are the view's own and a wide route costs nothing extra. A card's
structure is a picture rendered for the zoom it is shown at (in steps of √2, at the screen's
pixel ratio), so it stays sharp however far you zoom in; until the sharper picture is ready the
nearest one already made is drawn, so zooming never stalls on RDKit. The routes themselves
are plain JSON, so :meth:`RoutesPanel.state` saves them with the document and the pane comes
back exactly as it was left.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from PySide6.QtCore import QEvent, QPointF, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QGuiApplication,
    QImage,
    QPainter,
    QPen,
    QPixmap,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsSimpleTextItem,
    QGraphicsView,
    QHeaderView,
    QLabel,
    QMenu,
    QSplitter,
    QStyleOptionGraphicsItem,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)
from rdkit import Chem, rdBase

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.table.structure_copy import FORMATS, copy_structures
from chemsuite.app.widgets.flow_layout import FlowRow, compact, field
from chemsuite.core.depict import render_png
from chemsuite.retro.routes import leaf_stock_counts, route_tree, sort_routes

CARD_W = 150
CARD_H = 104
CARD_GAP_Y = 14
CARD_PAD = 6  # between a card's frame and its structure
STEP_W = 86
STEP_H = 22
BRANCH = 34
IN_STOCK = QColor("#2e7d32")
NOT_IN_STOCK = QColor("#ef6c00")
ROUTE_COLUMNS = ("#", "Status", "Score", "Steps", "In stock")
ZOOM_MIN = 0.2
ZOOM_MAX = 4.0
SCENE_MARGIN = 24
#: structure pictures are made at these scales (device pixels per scene unit)
PICTURE_MIN, PICTURE_MAX = 1.0, 16.0
#: pictures of one molecule kept at once; the ones farthest from the current zoom go first
PICTURES_PER_MOLECULE = 3
#: how long zooming must pause before sharper pictures are made
RENDER_DELAY_MS = 90


@dataclass
class _Placed:
    """A laid-out subtree: its size and where its product card sits inside it."""

    width: float
    height: float
    anchor: float  # the vertical middle of the product card, from the top of the subtree


class _StepButton(QGraphicsItem):
    """A numbered step on a branch line. Clicking it shows the reaction's metadata."""

    def __init__(self, number: int, step: dict, panel: RoutesPanel) -> None:
        super().__init__()
        self.number = number
        self.step = step
        self.panel = panel
        self.selected = False
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip(_step_tooltip(step))

    def boundingRect(self) -> QRectF:
        return QRectF(0, 0, STEP_W, STEP_H)

    def paint(self, painter: QPainter, _option, _widget=None) -> None:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        palette = self.panel.palette()
        fill = palette.highlight().color() if self.selected else palette.button().color()
        text = palette.highlightedText().color() if self.selected else palette.buttonText().color()
        painter.setPen(QPen(palette.mid().color(), 1))
        painter.setBrush(QBrush(fill))
        painter.drawRoundedRect(self.boundingRect().adjusted(0.5, 0.5, -0.5, -0.5), 4, 4)
        painter.setPen(QPen(text))
        font = painter.font()
        font.setPointSizeF(max(8.0, font.pointSizeF() - 1))
        painter.setFont(font)
        painter.drawText(self.boundingRect(), Qt.AlignmentFlag.AlignCenter, f"Step {self.number}")

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            self.panel.show_step(self)
            event.accept()
            return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event) -> None:
        self.panel.step_menu(self, event.screenPos())
        event.accept()


class _MoleculeCard(QGraphicsItem):
    """A structure in its own well, framed green (in stock) or amber (still to be made)."""

    def __init__(self, node: dict, panel: RoutesPanel) -> None:
        super().__init__()
        self.node = node
        self.smiles = str(node.get("smiles") or "")
        self.panel = panel
        self.in_stock = bool(node.get("in_stock"))
        self.setToolTip(_card_tooltip(node))

    def boundingRect(self) -> QRectF:
        return QRectF(0, 0, CARD_W, CARD_H)

    def paint(self, painter: QPainter, _option, _widget=None) -> None:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        rect = self.boundingRect().adjusted(1, 1, -1, -1)
        painter.setPen(QPen(IN_STOCK if self.in_stock else NOT_IN_STOCK, 2))
        painter.setBrush(QBrush(self.panel.palette().base().color()))
        painter.drawRoundedRect(rect, 5, 5)
        device = painter.device()
        ratio = device.devicePixelRatioF() if device is not None else 1.0
        zoom = QStyleOptionGraphicsItem.levelOfDetailFromTransform(painter.worldTransform())
        pixmap = self.panel.picture(self.smiles, zoom * ratio)
        if pixmap is not None and not pixmap.isNull():
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
            painter.drawPixmap(_inner_rect(), pixmap, QRectF(pixmap.rect()))

    def contextMenuEvent(self, event) -> None:
        self.panel.card_menu(self, event.screenPos())
        event.accept()


class _Canvas(QGraphicsView):
    """Pans with a drag or a scroll; zooms around the pointer with ⌘/Ctrl-scroll or a pinch.

    Until the user pans or zooms, the route is fitted again whenever the pane is resized, so a
    route drawn before the pane had its real size is not left tiny.
    """

    blankClicked = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.setScene(QGraphicsScene(self))
        self.setRenderHints(
            QPainter.RenderHint.Antialiasing | QPainter.RenderHint.SmoothPixmapTransform
        )
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setFrameShape(QGraphicsView.Shape.NoFrame)
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.moved = False  # the user panned or zoomed since the route was fitted

    @property
    def zoom(self) -> float:
        return self.transform().m11()

    def zoom_by(self, factor: float, at: QPointF | None = None) -> None:
        """Zoom by *factor*, kept within the limits, holding the scene point under *at* still."""
        target = min(ZOOM_MAX, max(ZOOM_MIN, self.zoom * factor))
        factor = target / self.zoom
        if abs(factor - 1.0) < 1e-6:
            return
        self.moved = True
        if at is None:
            self.scale(factor, factor)
            return
        anchor = self.transformationAnchor()
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.NoAnchor)
        before = self.mapToScene(at.toPoint())
        self.scale(factor, factor)
        after = self.mapToScene(at.toPoint())
        shift = after - before
        self.translate(shift.x(), shift.y())
        self.setTransformationAnchor(anchor)

    def wheelEvent(self, event) -> None:
        delta = event.angleDelta()
        if event.modifiers() & (
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
        ):
            if delta.y():
                self.zoom_by(math.pow(1.0015, delta.y()), event.position())
            event.accept()
            return
        self.moved = True
        super().wheelEvent(event)  # scrolls: pans with a trackpad or a mouse wheel

    def viewportEvent(self, event) -> bool:
        if (  # a trackpad pinch (macOS, Wayland)
            event.type() == QEvent.Type.NativeGesture
            and event.gestureType() == Qt.NativeGestureType.ZoomNativeGesture
        ):
            self.zoom_by(1.0 + event.value(), event.position())
            return True
        return super().viewportEvent(event)

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton and self.itemAt(event.pos()) is None:
            self.blankClicked.emit()
            self.moved = True  # a drag on empty space pans
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton and self.itemAt(event.pos()) is None:
            self.fit()
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if not self.moved:
            self.fit()

    def fit(self) -> None:
        """Show the whole route, never magnified past its natural size."""
        self.moved = False
        scene = self.scene()
        if scene is None or scene.sceneRect().isEmpty():
            return
        self.resetTransform()
        rect = scene.sceneRect()
        view = self.viewport().rect()
        if rect.width() > view.width() or rect.height() > view.height():
            self.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)
            if self.zoom < ZOOM_MIN:
                self.resetTransform()
                self.scale(ZOOM_MIN, ZOOM_MIN)
        self.centerOn(rect.center())


class RoutesPanel(ToolPanel):
    """The pane: a target, its routes, and the tree of the chosen one."""

    tool_id = "retro_routes"

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__()
        self.ctx = ctx
        self._targets: list[dict] = []
        self._routes: list[dict] = []
        #: (smiles, dark) -> {picture scale: pixmap}; None for a structure RDKit cannot draw
        self._pictures: dict[tuple[str, bool], dict[float, QPixmap | None]] = {}
        self._wanted: set[tuple[str, bool, float]] = set()
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(RENDER_DELAY_MS)
        self._render_timer.timeout.connect(self._render_wanted)
        self._step: _StepButton | None = None

        self.target = QComboBox()
        self.target.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.target.setMinimumContentsLength(28)
        compact(self.target)
        self.target.currentIndexChanged.connect(self._target_chosen)
        controls = FlowRow()
        controls.addWidget(field("Target:", self.target))

        self.routes = QTableWidget()
        self.routes.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.routes.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.routes.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.routes.verticalHeader().hide()
        self.routes.verticalHeader().setDefaultSectionSize(22)
        self.routes.horizontalHeader().setStretchLastSection(True)
        self.routes.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.routes.currentCellChanged.connect(
            lambda row, *_a: self._show_route(row) if row >= 0 else None
        )

        self.canvas = _Canvas()
        self.canvas.blankClicked.connect(self.hide_step)

        self.steps = QTableWidget(0, 2)
        self.steps.setHorizontalHeaderLabels(["Field", "Value"])
        self.steps.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.steps.verticalHeader().hide()
        self.steps.verticalHeader().setDefaultSectionSize(20)
        self.steps.horizontalHeader().setStretchLastSection(True)
        self.steps.hide()  # until a step is opened

        self.empty = QLabel("Chemistry ▸ Retrosynthesis fills this pane.")
        self.empty.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.empty.setEnabled(False)

        self.splitter = QSplitter(Qt.Orientation.Vertical)
        self.splitter.addWidget(self.routes)
        self.splitter.addWidget(self.canvas)
        self.splitter.addWidget(self.steps)
        self.splitter.setStretchFactor(1, 1)
        self.splitter.setChildrenCollapsible(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        layout.addWidget(controls)
        layout.addWidget(self.empty)
        layout.addWidget(self.splitter, 1)
        self._show_empty(True)

    # -- filling ------------------------------------------------------------------------------
    def title(self) -> str:
        if not self._targets:
            return "Routes"
        return f"Routes ({len(self._targets)})"

    def add_targets(self, targets: list[dict], *, focus: str = "") -> None:
        """Add or replace searched targets. Each is ``{"smiles": …, "routes": [...]}``."""
        for entry in targets:
            smiles = str(entry.get("smiles") or "").strip()
            if not smiles:
                continue
            routes = sort_routes([dict(route) for route in entry.get("routes") or []])
            found = next(
                (item for item in self._targets if item["smiles"] == smiles),
                None,
            )
            if found is not None:
                found["routes"] = routes
            else:
                self._targets.append({"smiles": smiles, "routes": routes})
        self._fill_targets(focus or (targets[-1].get("smiles", "") if targets else ""))
        self.titleChanged.emit(self.title())

    def _fill_targets(self, focus: str = "") -> None:
        self.target.blockSignals(True)
        self.target.clear()
        index = 0
        for position, entry in enumerate(self._targets):
            smiles = entry["smiles"]
            count = len(entry.get("routes") or [])
            label = smiles if len(smiles) <= 44 else smiles[:41] + "…"
            self.target.addItem(f"{label}  ({count})", smiles)
            if focus and smiles == focus:
                index = position
        self.target.blockSignals(False)
        self._show_empty(not self._targets)
        if self._targets:
            self.target.setCurrentIndex(index)
            self._target_chosen(index)

    def _show_empty(self, empty: bool) -> None:
        self.empty.setVisible(empty)
        self.splitter.setVisible(not empty)

    def _target_chosen(self, index: int) -> None:
        entry = self._targets[index] if 0 <= index < len(self._targets) else None
        self._routes = list(entry.get("routes") or []) if entry else []
        self._fill_routes()
        if self._routes:
            self.routes.setCurrentCell(0, 0)
            self._show_route(0)
        else:
            self.canvas.scene().clear()
            self.hide_step()

    def _fill_routes(self) -> None:
        extras = _extra_scores(self._routes)
        headers = [*ROUTE_COLUMNS, *extras]
        table = self.routes
        table.blockSignals(True)
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setRowCount(len(self._routes))
        for row, route in enumerate(self._routes):
            in_stock, total = leaf_stock_counts(route)
            solved = bool(route.get("solved"))
            cells = [
                str(row + 1),
                "Solved" if solved else "Unsolved",
                str(route.get("score") or "—"),
                str(len(route.get("steps") or [])),
                f"{in_stock}/{total}",
                *[_score_cell(route, name) for name in extras],
            ]
            for column, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if column == 1:
                    item.setForeground(QBrush(IN_STOCK if solved else NOT_IN_STOCK))
                table.setItem(row, column, item)
        table.setMaximumHeight(min(160, 26 + 22 * max(1, len(self._routes))))
        table.blockSignals(False)

    # -- drawing ------------------------------------------------------------------------------
    def _show_route(self, row: int) -> None:
        if not (0 <= row < len(self._routes)):
            return
        self.hide_step()
        scene = self.canvas.scene()
        scene.clear()
        tree = route_tree(self._routes[row])
        if tree is None:
            scene.setSceneRect(QRectF())
            return
        steps = list(self._routes[row].get("steps") or [])
        numbers = {
            str(step.get("reaction") or ""): index for index, step in enumerate(steps, start=1)
        }
        self._trim_pictures()
        placed = self._place(tree, scene, 0.0, 0.0, numbers)
        margin = SCENE_MARGIN
        scene.setSceneRect(-margin, -margin, placed.width + 2 * margin, placed.height + 2 * margin)
        self.canvas.fit()

    def _place(
        self, node: dict, scene: QGraphicsScene, x: float, y: float, numbers: dict[str, int]
    ) -> _Placed:
        """Lay *node* out with its precursors to the left, and return the space it took.

        The product card is placed at the right edge of its own subtree, level with the middle
        of the precursors that make it.
        """
        reactions = [item for item in node.get("reactions") or [] if isinstance(item, dict)]
        reaction = reactions[0] if reactions else None
        reactants = (
            [item for item in reaction.get("reactants") or [] if isinstance(item, dict)]
            if reaction
            else []
        )
        if not reactants:
            card = _MoleculeCard(node, self)
            card.setPos(x, y)
            scene.addItem(card)
            return _Placed(CARD_W, CARD_H, CARD_H / 2)

        children: list[_Placed] = []
        offset = y
        branch_x = x
        for item in reactants:
            placed = self._place(item, scene, x, offset, numbers)
            children.append(placed)
            offset += placed.height + CARD_GAP_Y
            branch_x = max(branch_x, x + placed.width)
        height = offset - y - CARD_GAP_Y
        anchors = [
            sum(child.height + CARD_GAP_Y for child in children[:index]) + child.anchor
            for index, child in enumerate(children)
        ]
        middle = (anchors[0] + anchors[-1]) / 2
        step_x = branch_x + BRANCH
        card_x = step_x + STEP_W + BRANCH
        card = _MoleculeCard(node, self)
        card.setPos(card_x, y + middle - CARD_H / 2)
        scene.addItem(card)
        pen = QPen(self.palette().mid().color(), 1.5)
        for child, anchor in zip(children, anchors, strict=True):
            scene.addLine(x + child.width, y + anchor, branch_x + BRANCH / 2, y + anchor, pen)
        if len(children) > 1:  # the vertical that gathers the precursors
            scene.addLine(
                branch_x + BRANCH / 2,
                y + anchors[0],
                branch_x + BRANCH / 2,
                y + anchors[-1],
                pen,
            )
        scene.addLine(branch_x + BRANCH / 2, y + middle, step_x, y + middle, pen)
        scene.addLine(step_x + STEP_W, y + middle, card_x, y + middle, pen)
        step = _StepButton(
            numbers.get(str(reaction.get("smiles") or ""), 0),
            _step_of(reaction, node),
            self,
        )
        step.setPos(step_x, y + middle - STEP_H / 2)
        scene.addItem(step)
        if len(reactions) > 1:  # AiZynthFinder found more than one way to this molecule
            note = QGraphicsSimpleTextItem(f"+{len(reactions) - 1} more")
            note.setBrush(QBrush(self.palette().mid().color()))
            font = QFont(note.font())
            font.setPointSizeF(max(7.0, font.pointSizeF() - 2))
            note.setFont(font)
            note.setPos(step_x, y + middle + STEP_H)
            scene.addItem(note)
        width = max(card_x + CARD_W, branch_x)
        return _Placed(width, max(height, CARD_H), middle)

    # -- structure pictures -------------------------------------------------------------------
    def picture(self, smiles: str, scale: float) -> QPixmap | None:
        """The structure for a card shown at *scale* device pixels per scene unit.

        Returns the picture made for that scale if there is one. Otherwise the nearest one
        already made is returned and the right one is queued, to be made once zooming pauses;
        only a molecule with no picture at all is drawn on the spot.
        """
        dark = self.palette().base().color().lightness() < 128
        wanted = picture_scale(scale)
        made = self._pictures.setdefault((smiles, dark), {})
        if wanted in made:
            return made[wanted]
        if not made:
            made[wanted] = _draw(smiles, dark, wanted)
            return made[wanted]
        if any(pixmap is None for pixmap in made.values()):
            return None  # RDKit cannot draw it at any size
        self._wanted.add((smiles, dark, wanted))
        self._render_timer.start()
        sharper = [have for have in made if have >= wanted]
        return made[min(sharper) if sharper else max(made)]

    def _render_wanted(self) -> None:
        wanted, self._wanted = self._wanted, set()
        for smiles, dark, scale in wanted:
            made = self._pictures.setdefault((smiles, dark), {})
            if scale in made:
                continue
            made[scale] = _draw(smiles, dark, scale)
            for extra in sorted(made, key=lambda have: abs(have - scale))[PICTURES_PER_MOLECULE:]:
                del made[extra]
        if wanted:
            self.canvas.viewport().update()

    def _trim_pictures(self) -> None:
        """Keep one small picture per molecule when a new route is drawn."""
        self._wanted.clear()
        for made in self._pictures.values():
            for extra in sorted(made)[1:]:
                del made[extra]

    # -- steps --------------------------------------------------------------------------------
    def show_step(self, button: _StepButton) -> None:
        """Show one step's fields under the canvas."""
        if self._step is not None and self._step is not button:
            self._step.selected = False
            self._step.update()
        self._step = button
        button.selected = True
        button.update()
        rows = _step_rows(button.step)
        self.steps.setRowCount(len(rows))
        for row, (name, value) in enumerate(rows):
            self.steps.setItem(row, 0, QTableWidgetItem(name))
            self.steps.setItem(row, 1, QTableWidgetItem(value))
        self.steps.setVisible(bool(rows))
        self.steps.resizeColumnToContents(0)

    def hide_step(self) -> None:
        if self._step is not None:
            self._step.selected = False
            self._step.update()
            self._step = None
        self.steps.setRowCount(0)
        self.steps.hide()

    def step_menu(self, button: _StepButton, at) -> None:
        menu = QMenu(self)
        smiles = str(button.step.get("reaction") or "")
        copy = menu.addAction("Copy Reaction SMILES")
        draw = menu.addAction("Draw Reaction in the Sketcher")
        copy.setEnabled(bool(smiles))
        draw.setEnabled(bool(smiles))
        chosen = menu.exec(at.toPoint() if hasattr(at, "toPoint") else at)
        if chosen is copy and smiles:
            QGuiApplication.clipboard().setText(smiles)
            self.ctx.show_message("Copied the reaction as reaction SMILES")
        elif chosen is draw and smiles:
            self._draw_reaction(smiles)

    def card_menu(self, card: _MoleculeCard, at) -> None:
        menu = QMenu(self)
        actions = {}
        for key, title in FORMATS:
            actions[menu.addAction(f"Copy as {title}")] = key
        menu.addSeparator()
        add = menu.addAction("Add to Table")
        chosen = menu.exec(at.toPoint() if hasattr(at, "toPoint") else at)
        if chosen in actions:
            mol = _parse(card.smiles)
            self.ctx.show_message(copy_structures([mol], actions[chosen]))
        elif chosen is add and card.smiles and self.ctx.add_compound(card.smiles) is not None:
            self.ctx.show_message(f"Added {card.smiles} as a new row")

    def _draw_reaction(self, smiles: str) -> None:
        from chemsuite.app.sketcher import SketcherDialog

        dialog = SketcherDialog(smiles, "Reaction Step", self, reaction_only=True)
        self.ctx.show_dialog(dialog)

    # -- the document -------------------------------------------------------------------------
    def state(self) -> dict:
        """What reopens this pane: the targets, their routes and which one is shown."""
        return {
            "targets": self._targets,
            "target": self.target.currentData() or "",
            "route": self.routes.currentRow(),
        }

    def set_state(self, state: dict) -> None:
        targets = state.get("targets")
        if isinstance(targets, list):
            self.add_targets([item for item in targets if isinstance(item, dict)])
        focus = str(state.get("target") or "")
        if focus:
            index = self.target.findData(focus)
            if index >= 0:
                self.target.setCurrentIndex(index)
        row = state.get("route")
        if isinstance(row, int) and 0 <= row < len(self._routes):
            self.routes.setCurrentCell(row, 0)
        self.titleChanged.emit(self.title())


# -- helpers ------------------------------------------------------------------------------------
def _parse(smiles: str) -> Chem.Mol | None:
    if not smiles:
        return None
    with rdBase.BlockLogs():
        return Chem.MolFromSmiles(smiles)


def _inner_rect() -> QRectF:
    """Where a card's structure goes, in the card's own coordinates."""
    return QRectF(CARD_PAD, CARD_PAD, CARD_W - 2 * CARD_PAD, CARD_H - 2 * CARD_PAD)


def picture_scale(scale: float) -> float:
    """*scale* rounded up to the next step of √2, within the picture limits.

    Rounding up means a picture is never stretched; the steps mean a zoom makes a new picture
    only when it crosses one.
    """
    if not math.isfinite(scale) or scale <= PICTURE_MIN:
        return PICTURE_MIN
    step = math.ceil(round(2 * math.log2(scale), 6)) / 2
    return min(PICTURE_MAX, 2.0**step)


def _draw(smiles: str, dark: bool, scale: float) -> QPixmap | None:
    """A molecule drawn to fill a card's inner box at *scale* device pixels per scene unit.

    The box and the bond length grow together, so every scale draws the same layout, just
    with more pixels: swapping one picture for a sharper one never shifts the structure.
    """
    mol = _parse(smiles)
    if mol is None:
        return None
    inner = _inner_rect()
    width = max(1, round(inner.width() * scale))
    height = max(1, round(inner.height() * scale))
    png = render_png(mol, width, height, scale=scale, dark=dark)
    image = QImage.fromData(png, "PNG")
    if image.isNull():
        return None
    return QPixmap.fromImage(image)


def _step_of(reaction: dict, product: dict) -> dict:
    """The step a branch stands for: its reaction, product and metadata."""
    return {
        "product": str(product.get("smiles") or ""),
        "reaction": str(reaction.get("smiles") or ""),
        "metadata": reaction.get("metadata") or {},
        "reactants": [
            {"smiles": str(item.get("smiles") or ""), "in_stock": bool(item.get("in_stock"))}
            for item in reaction.get("reactants") or []
            if isinstance(item, dict)
        ],
    }


def _step_rows(step: dict) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    if step.get("reaction"):
        rows.append(("reaction", str(step["reaction"])))
    metadata = step.get("metadata")
    if isinstance(metadata, dict):
        for key, value in metadata.items():
            if value in (None, ""):
                continue
            rows.append((str(key).replace("_", " "), _value_text(value)))
    return rows


def _step_tooltip(step: dict) -> str:
    metadata = step.get("metadata")
    if not isinstance(metadata, dict) or not metadata:
        return "Click for this reaction's details"
    parts = []
    name = str(metadata.get("policy_name") or "").strip()
    if name:
        parts.append(name)
    if metadata.get("policy_probability") not in (None, ""):
        parts.append(f"p = {_value_text(metadata['policy_probability'])}")
    classification = str(metadata.get("classification") or "").strip()
    if classification:
        parts.append(classification)
    return " · ".join(parts) or "Click for this reaction's details"


def _card_tooltip(node: dict) -> str:
    smiles = str(node.get("smiles") or "")
    stock = node.get("stock")
    names = [str(item) for item in stock] if isinstance(stock, list) else []
    where = (
        f"In stock: {', '.join(names)}"
        if names
        else ("In stock" if node.get("in_stock") else "Not in stock")
    )
    return f"{smiles}\n{where}"


def _extra_scores(routes: list[dict]) -> list[str]:
    """Score names other than the state score, in first-seen order."""
    names: list[str] = []
    for route in routes:
        scores = route.get("scores")
        if not isinstance(scores, dict):
            continue
        for key in scores:
            if str(key).lower() != "state score" and str(key) not in names:
                names.append(str(key))
    return names


def _score_cell(route: dict, name: str) -> str:
    scores = route.get("scores")
    if not isinstance(scores, dict) or name not in scores:
        return "—"
    return _value_text(scores[name])


def _value_text(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, (list, tuple)):
        return ", ".join(_value_text(item) for item in value)
    if isinstance(value, dict):
        return ", ".join(f"{key}: {_value_text(item)}" for key, item in value.items())
    return str(value)


def restore(ctx: ToolContext, state: dict) -> RoutesPanel:
    """Reopen the pane from a saved document (see ``Tool.restore``)."""
    panel = RoutesPanel(ctx)
    panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel
