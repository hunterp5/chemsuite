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
"""A plot's settings: the controls in its options popup, and the :class:`PlotSpec` they make.

:class:`PlotControls` fills the popup's form, shows only the rows the current kind and columns
use, and says when anything changed; the panel draws.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from PySide6.QtCore import QObject, QSize, Qt, Signal
from PySide6.QtGui import QColor, QIcon, QLinearGradient, QPainter, QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QPushButton,
    QSpinBox,
    QWidget,
)

from chemsuite.analysis.medchem_space import GUIDES
from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER, STRUCTURE, TEXT
from chemsuite.plotting.data import colorable_columns, number_columns
from chemsuite.plotting.stats import AGGREGATES, CURVES
from chemsuite.plotting.style import COLOR_MAPS, DEFAULT_MAP, parse_range
from chemsuite.sar.constellation import MAPS, SUMMARIES

from ..widgets.flow_layout import compact
from .spec import PlotSpec, Range

NONE = "(none)"
KINDS = {
    "scatter": "Scatter",
    "line": "Line",
    "histogram": "Histogram",
    "bar": "Bar Chart",
    "box": "Box Plot",
    "density": "Density",
    "radar": "Radar",
    "constellation": "Constellation",
    "scatter3d": "3D Scatter",
    "landscape": "3D Landscape",
}
#: drawn by ``scene3d`` on the GPU, imported only when one is first shown
KINDS_3D = ("scatter3d", "landscape")
#: kinds plotting one number column against another
XY_KINDS = ("scatter", "line", "density", *KINDS_3D)
#: kinds colouring their points by a column
COLORED = ("scatter", "line", "scatter3d")
#: kinds sizing their points by a column
SIZED = ("scatter", "scatter3d")
CONSTELLATION = "constellation"
#: kinds always coloured along a colour map (counts, heights)
MAPPED = ("density", "landscape")
POINT_SIZES = (4, 18)
#: every setting saved with the document, besides the camera, title and cards
STATE_KEYS = (
    "kind",
    "x",
    "y",
    "color",
    "size",
    "bins",
    "fit",
    "curve",
    "guide",
    "aggregate",
    "errors",
    "columns",
    "z",
    "drop",
    "cmap",
    "reverse",
    "color_range",
    "size_range",
    "point_sizes",
    "structure",
    "summary",
    "layout",
    "min_compounds",
    "core_percent",
    "labels",
)


def gradient_icon(key: str, width: int = 48, height: int = 12) -> QIcon:
    """A swatch of colour map *key*, for its entry in the Colors list."""
    pixmap = QPixmap(width, height)
    gradient = QLinearGradient(0, 0, width, 0)
    anchors = COLOR_MAPS[key][1]
    for i, rgb in enumerate(anchors):
        gradient.setColorAt(i / max(len(anchors) - 1, 1), QColor(*rgb))
    painter = QPainter(pixmap)
    painter.fillRect(pixmap.rect(), gradient)
    painter.end()
    return QIcon(pixmap)


class RangeEdit(QWidget):
    """Two ends of a value range, each a number or empty (Auto: the data's own end)."""

    changed = Signal()

    def __init__(self, tip: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(8)
        self.low, self.high = QLineEdit(), QLineEdit()
        for edit, end in ((self.low, "Lowest"), (self.high, "Highest")):
            edit.setPlaceholderText("Auto")
            edit.setToolTip(f"{end} value of the range ({tip}); empty: the data's own")
            edit.setMaximumWidth(80)
            edit.editingFinished.connect(self._finished)
        row.addWidget(self.low)
        row.addWidget(QLabel("to"))
        row.addWidget(self.high)
        row.addStretch(1)
        self._last = self.value()

    def value(self) -> Range:
        return parse_range(self.low.text()), parse_range(self.high.text())

    def set_value(self, value: object) -> None:
        ends = value if isinstance(value, Sequence) and len(value) == 2 else (None, None)
        for edit, end in zip((self.low, self.high), ends, strict=True):
            number = end if isinstance(end, int | float) and not isinstance(end, bool) else None
            edit.setText("" if number is None else f"{number:g}")
        self._last = self.value()

    def _finished(self) -> None:
        for edit in (self.low, self.high):  # what is not a number goes back to Auto
            if parse_range(edit.text()) is None:
                edit.clear()
        if self.value() != self._last:
            self._last = self.value()
            self.changed.emit()


class PlotControls(QObject):
    """The plot's controls, in the rows of *form*; :attr:`changed` when one changes."""

    changed = Signal()

    def __init__(
        self, table: Table, form: QFormLayout, kind: str, parent: QWidget | None = None
    ) -> None:
        super().__init__(parent)
        self.table = table
        self.form = form
        self.kind = QComboBox()
        for key, title in KINDS.items():
            if key == KINDS_3D[0]:
                self.kind.insertSeparator(self.kind.count())
            self.kind.addItem(title, key)
        self.kind.setCurrentIndex(max(0, self.kind.findData(kind)))  # unknown kind: scatter
        self.x = QComboBox()
        self.y = QComboBox()
        self.z = QComboBox()
        self.color = QComboBox()
        self.size = QComboBox()
        self.bins = QSpinBox()
        self.bins.setRange(0, 200)
        self.bins.setSpecialValueText("Auto")
        self.fit = QCheckBox("Fit line")
        self.fit.setToolTip("Least-squares line; one per category when coloured by text")
        self.curve = QComboBox()
        self.curve.addItem("None", None)
        for key, title in CURVES.items():
            self.curve.addItem(title, key)
        self.guide = QComboBox()
        self.guide.addItem("None", None)
        for key, guide in GUIDES.items():
            self.guide.addItem(guide.title, key)
            self.guide.setItemData(
                self.guide.count() - 1,
                f"Expects X: {guide.x_hint}, Y: {guide.y_hint}",
                Qt.ItemDataRole.ToolTipRole,
            )
        self.guide.setToolTip("Reference regions or lines drawn behind the points")
        self.aggregate = QComboBox()
        for key, title in AGGREGATES.items():
            self.aggregate.addItem(title, key)
        self.aggregate.setToolTip("How each bar sums up its rows' values")
        self.errors = QCheckBox("Error bars")
        self.errors.setToolTip("± one standard deviation (with Mean)")
        self.drop = QCheckBox("Drop lines")
        self.drop.setToolTip("A line from each point down to the floor, so heights read easily")
        self.axes: list[str] = []  # radar spokes
        self.axes_button = QPushButton("None")
        self.axes_button.setToolTip("The number columns drawn as spokes (at least three)")
        self.axes_menu = QMenu(self.axes_button)
        self.axes_menu.aboutToShow.connect(self._fill_axes_menu)
        self.axes_button.setMenu(self.axes_menu)
        self.axes_button.setStyleSheet("text-align: left; padding: 3px 8px;")
        self._make_scale_controls()
        self._make_constellation_controls()
        compact(
            self.x,
            self.y,
            self.z,
            self.color,
            self.size,
            self.curve,
            self.guide,
            self.aggregate,
            chars=16,
        )
        form.addRow("Plot:", self.kind)
        self.labels: dict[QWidget, QLabel] = {}
        for title, widget in (
            ("Structure", self.structure),
            ("X", self.x),
            ("Y", self.y),
            ("Z", self.z),
            ("Color", self.color),
            ("Size", self.size),
            ("Bins", self.bins),
            ("Curve", self.curve),
            ("Guide", self.guide),
            ("Summary", self.aggregate),
            ("Axes", self.axes_button),
            ("Summary", self.core_summary),
            ("Map", self.layout),
            ("Min compounds", self.min_compounds),
            ("Core size", self.core_percent),
            ("Colors", self.colors),
            ("Color range", self.color_range),
            ("Size range", self.size_range),
            ("Point size", self.point_sizes),
        ):
            label = QLabel(title + ":")
            self.labels[widget] = label
            form.addRow(label, widget)
        form.addRow("", self.fit)
        form.addRow("", self.errors)
        form.addRow("", self.drop)
        form.addRow("", self.series_labels)
        self.fill_columns(initial=True)
        self.kind.currentIndexChanged.connect(self._kind_changed)
        for combo in (self.x, self.y, self.z, self.curve, self.guide, self.aggregate, self.cmap):
            combo.currentIndexChanged.connect(self.changed)
        for combo in (self.structure, self.core_summary, self.layout):
            combo.currentIndexChanged.connect(self.changed)
        for combo in (self.color, self.size):
            combo.currentIndexChanged.connect(self._column_changed)
        for spin in (self.bins, self.small, self.large, self.min_compounds, self.core_percent):
            spin.valueChanged.connect(self.changed)
        for box in (self.fit, self.errors, self.drop, self.reverse, self.series_labels):
            box.toggled.connect(self.changed)
        for edit in (self.color_range, self.size_range):
            edit.changed.connect(self.changed)

    def _make_scale_controls(self) -> None:
        self.cmap = QComboBox()
        for key, (title, _anchors) in COLOR_MAPS.items():
            self.cmap.addItem(gradient_icon(key), title, key)
        self.cmap.setIconSize(QSize(48, 12))
        self.cmap.setCurrentIndex(self.cmap.findData(DEFAULT_MAP))
        self.cmap.setToolTip("The colours number values run through, low to high")
        self.reverse = QCheckBox("Reverse")
        self.reverse.setToolTip("Run the colours high to low")
        self.colors = _row(self.cmap, self.reverse)
        self.color_range = RangeEdit("colours")
        self.size_range = RangeEdit("point sizes")
        self.small, self.large = QSpinBox(), QSpinBox()
        for spin, end in ((self.small, "Smallest"), (self.large, "Largest")):
            spin.setRange(1, 60)
            spin.setSuffix(" px")
            spin.setToolTip(f"{end} point, for the low (high) end of the size range")
        self.small.setValue(POINT_SIZES[0])
        self.large.setValue(POINT_SIZES[1])
        self.point_sizes = _row(self.small, QLabel("to"), self.large)

    def _make_constellation_controls(self) -> None:
        self.structure = QComboBox()
        self.structure.setToolTip("The structures whose cores are found")
        self.core_summary = QComboBox()
        for key, title in SUMMARIES.items():
            self.core_summary.addItem(title, key)
        self.core_summary.setToolTip(
            "How each core sums up its compounds' Color values (Mean of a 0/1 column: the "
            "share of actives)"
        )
        self.layout = QComboBox()
        for key, title in MAPS.items():
            self.layout.addItem(title, key)
        self.layout.setToolTip(
            "How the cores are placed, by their Morgan fingerprints (under five cores: PCA)"
        )
        self.min_compounds = QSpinBox()
        self.min_compounds.setRange(2, 10_000)
        self.min_compounds.setToolTip("Cores shared by fewer compounds are left out")
        self.core_percent = QSpinBox()
        self.core_percent.setRange(10, 100)
        self.core_percent.setValue(67)
        self.core_percent.setSuffix(" %")
        self.core_percent.setToolTip(
            "The share of a compound's heavy atoms a core must keep (two thirds in the paper)"
        )
        self.series_labels = QCheckBox("Series labels")
        self.series_labels.setChecked(True)
        self.series_labels.setToolTip("Number each analog series beside its largest core")
        compact(self.structure, self.core_summary, self.layout, chars=16)

    @property
    def controls(self) -> tuple[QObject, ...]:
        """Everything that says when it changed (blocked while settings are put back)."""
        return (
            self.kind,
            self.x,
            self.y,
            self.z,
            self.color,
            self.size,
            self.bins,
            self.fit,
            self.curve,
            self.guide,
            self.aggregate,
            self.errors,
            self.drop,
            self.cmap,
            self.reverse,
            self.color_range,
            self.size_range,
            self.small,
            self.large,
            self.structure,
            self.core_summary,
            self.layout,
            self.min_compounds,
            self.core_percent,
            self.series_labels,
        )

    # -- the spec -------------------------------------------------------------------------------
    def spec(self) -> PlotSpec:
        def name(combo: QComboBox) -> str | None:
            text = combo.currentText()
            return text if text and text != NONE else None

        kind = self.kind.currentData()
        return PlotSpec(
            kind=kind,
            x=name(self.x) if kind not in ("radar", CONSTELLATION) else None,
            y=name(self.y) if kind in (*XY_KINDS, "box", "bar") else None,
            color=name(self.color) if kind in (*COLORED, CONSTELLATION) else None,
            size=name(self.size) if kind in SIZED else None,
            bins=self.bins.value() or None
            if kind in ("histogram", "density", "landscape")
            else None,
            fit=self.fit.isChecked() and kind == "scatter",
            curve=self.curve.currentData() if kind == "histogram" else None,
            guide=self.guide.currentData() if kind == "scatter" else None,
            aggregate=self.aggregate.currentData() if kind == "bar" else "mean",
            errors=self.errors.isChecked() and kind == "bar",
            columns=tuple(self.axes) if kind == "radar" else (),
            z=name(self.z) if kind in KINDS_3D else None,
            drop=self.drop.isChecked() and kind == "scatter3d",
            cmap=self.cmap.currentData(),
            reverse=self.reverse.isChecked(),
            color_range=self.color_range.value(),
            size_range=self.size_range.value(),
            point_sizes=self._point_sizes(),
            structure=name(self.structure),
            summary=self.core_summary.currentData(),
            layout=self.layout.currentData(),
            min_compounds=self.min_compounds.value(),
            core_percent=self.core_percent.value(),
            labels=self.series_labels.isChecked(),
        )

    def _point_sizes(self) -> tuple[float, float]:
        small, large = self.small.value(), self.large.value()
        return float(min(small, large)), float(max(small, large))

    def state(self) -> dict:
        """The settings, JSON-friendly (a hidden control keeps its value)."""
        spec = self.spec()
        return {
            "kind": spec.kind,
            "x": spec.x,
            "y": spec.y,
            "color": spec.color,
            "size": spec.size,
            "bins": spec.bins,
            "fit": spec.fit,
            "curve": spec.curve,
            "guide": spec.guide,
            "aggregate": self.aggregate.currentData(),
            "errors": self.errors.isChecked(),
            "columns": list(self.axes),
            "z": spec.z,
            "drop": self.drop.isChecked(),
            "cmap": spec.cmap,
            "reverse": spec.reverse,
            "color_range": list(spec.color_range),
            "size_range": list(spec.size_range),
            "point_sizes": [self.small.value(), self.large.value()],
            "structure": spec.structure,
            "summary": spec.summary,
            "layout": spec.layout,
            "min_compounds": spec.min_compounds,
            "core_percent": spec.core_percent,
            "labels": spec.labels,
        }

    def set_values(self, values: Mapping[str, object]) -> None:
        """Set controls by their :data:`STATE_KEYS` names, telling no one."""
        for control in self.controls:
            control.blockSignals(True)
        try:
            if "kind" in values:
                self.kind.setCurrentIndex(max(0, self.kind.findData(values["kind"])))
                self.fill_columns()
            combos = {"x": self.x, "y": self.y, "z": self.z, "color": self.color}
            combos["size"], combos["structure"] = self.size, self.structure
            pick = {"curve": self.curve, "guide": self.guide, "aggregate": self.aggregate}
            pick["summary"], pick["layout"] = self.core_summary, self.layout
            checks = {"fit": self.fit, "errors": self.errors, "drop": self.drop}
            checks["reverse"], checks["labels"] = self.reverse, self.series_labels
            spins = {"min_compounds": self.min_compounds, "core_percent": self.core_percent}
            setters = {
                "bins": lambda v: self.bins.setValue(v if isinstance(v, int) else 0),
                "columns": self.set_axes,
                "cmap": lambda v: self.cmap.setCurrentIndex(
                    self.cmap.findData(v if v in COLOR_MAPS else DEFAULT_MAP)
                ),
                "color_range": self.color_range.set_value,
                "size_range": self.size_range.set_value,
                "point_sizes": self._set_point_sizes,
            }
            for key, value in values.items():
                if key in combos:
                    combo = combos[key]
                    if value is None or (isinstance(value, str) and combo.findText(value) >= 0):
                        combo.setCurrentText(value or NONE)
                elif key in pick:
                    pick[key].setCurrentIndex(max(0, pick[key].findData(value)))
                elif key in checks:
                    checks[key].setChecked(bool(value))
                elif key in spins:
                    if isinstance(value, int) and not isinstance(value, bool):
                        spins[key].setValue(value)
                elif key in setters:
                    setters[key](value)
        finally:
            for control in self.controls:
                control.blockSignals(False)
        self._show_rows()

    def _set_point_sizes(self, value: object) -> None:
        ends = value if isinstance(value, Sequence) and len(value) == 2 else POINT_SIZES
        for spin, end, default in zip((self.small, self.large), ends, POINT_SIZES, strict=True):
            spin.setValue(int(end) if isinstance(end, int | float) else default)

    # -- the rows each kind uses ----------------------------------------------------------------
    def fill_columns(self, *, initial: bool = False) -> None:
        """The column lists (after the kind or the table's columns changed), and the rows."""
        table = self.table
        numbers = number_columns(table)
        kind = self.kind.currentData()
        texts = [n for n in table.column_names if table.kind(n) == TEXT]
        groupable = [n for n in table.column_names if table.kind(n) in (NUMBER, TEXT)]
        x_choices = {
            "scatter": numbers,
            "line": numbers,
            "density": numbers,
            "histogram": groupable,
            "bar": texts or groupable,  # the grouping column
            "box": [NONE, *texts],  # the grouping column (optional)
            "radar": [],
            "scatter3d": numbers,
            "landscape": numbers,
            CONSTELLATION: [],
        }[kind]
        lists = {
            self.x: x_choices,
            self.y: [NONE, *numbers] if kind == "bar" else numbers,
            self.z: numbers,
            self.color: [NONE, *(numbers if kind == CONSTELLATION else colorable_columns(table))],
            self.structure: [n for n in table.column_names if table.kind(n) == STRUCTURE],
            self.size: [NONE, *numbers],
        }
        for combo, names in lists.items():
            current = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(names)
            if current in names:
                combo.setCurrentText(current)
            elif initial and combo is self.y and len(names) > 1:
                combo.setCurrentIndex(1)
            elif combo is self.z:  # a column not already on X or Y
                taken = {self.x.currentText(), self.y.currentText()}
                spare = [n for n in names if n not in taken]
                if spare:
                    combo.setCurrentText(spare[0])
            elif combo is self.y and kind == "bar" and current not in names:
                combo.setCurrentIndex(min(1, len(names) - 1))  # a number column, not a count
            combo.blockSignals(False)
        self.axes = [a for a in self.axes if a in numbers]
        if initial or (kind == "radar" and len(self.axes) < 3):
            self.axes = self.axes if len(self.axes) >= 3 else numbers[:6]
        self._show_axes()
        self._show_rows()

    def _show_rows(self) -> None:
        kind = self.kind.currentData()
        shown = {
            "scatter": (self.x, self.y, self.color, self.size, self.fit, self.guide),
            "line": (self.x, self.y, self.color),
            "density": (self.x, self.y, self.bins),
            "histogram": (self.x, self.bins, self.curve),
            "bar": (self.x, self.y, self.aggregate, self.errors),
            "box": (self.x, self.y),
            "radar": (self.axes_button,),
            "scatter3d": (self.x, self.y, self.z, self.color, self.size, self.drop),
            "landscape": (self.x, self.y, self.z, self.bins),
            CONSTELLATION: (
                self.structure,
                self.color,
                self.layout,
                self.min_compounds,
                self.core_percent,
                self.series_labels,
                self.size_range,
                self.point_sizes,
            ),
        }[kind]
        numbered = self._is_number(self.color) and kind in (*COLORED, CONSTELLATION)
        sized = self._is_number(self.size) and kind in SIZED
        shown += (self.colors,) if numbered or kind in MAPPED else ()
        shown += (self.color_range,) if numbered else ()
        shown += (self.size_range, self.point_sizes) if sized else ()
        shown += (self.core_summary,) if numbered and kind == CONSTELLATION else ()
        for widget in self.labels:
            self.form.setRowVisible(widget, widget in shown)
        for widget in (self.fit, self.errors, self.drop, self.series_labels):
            self.form.setRowVisible(widget, widget in shown)
        self.labels[self.x].setText("Group by:" if kind in ("box", "bar") else "X:")
        self.labels[self.y].setText("Values:" if kind in ("box", "bar") else "Y:")
        self.labels[self.z].setText("Height:" if kind == "landscape" else "Z:")
        self.labels[self.bins].setText("Grid:" if kind == "landscape" else "Bins:")
        self.labels[self.size_range].setText(
            "Compounds:" if kind == CONSTELLATION else "Size range:"
        )
        self.bins.setToolTip(
            "Grid cells along each side (Auto picks from the number of points)"
            if kind == "landscape"
            else "Number of bins (Auto picks a width from the data)"
        )

    def _is_number(self, combo: QComboBox) -> bool:
        name = combo.currentText()
        return name in self.table and self.table.kind(name) == NUMBER

    def _kind_changed(self) -> None:
        self.fill_columns()
        self.changed.emit()

    def _column_changed(self) -> None:
        self._show_rows()  # a number column brings its scale rows
        self.changed.emit()

    # -- radar spokes ---------------------------------------------------------------------------
    def _fill_axes_menu(self) -> None:
        menu = self.axes_menu
        menu.clear()
        for name in number_columns(self.table):
            action = menu.addAction(name)
            action.setCheckable(True)
            action.setChecked(name in self.axes)
            action.toggled.connect(lambda on, n=name: self._toggle_axis(n, on))

    def _toggle_axis(self, name: str, on: bool) -> None:
        axes = [a for a in self.axes if a != name]
        if on:
            axes.append(name)
        self.set_axes(axes)
        self.changed.emit()

    def set_axes(self, names: object) -> None:
        numbers = set(number_columns(self.table))
        listed = names if isinstance(names, Sequence) and not isinstance(names, str) else ()
        self.axes = [n for n in listed if isinstance(n, str) and n in numbers]
        self._show_axes()

    def _show_axes(self) -> None:
        text = ", ".join(self.axes) if self.axes else "None"
        metrics = self.axes_button.fontMetrics()
        self.axes_button.setText(metrics.elidedText(text, Qt.TextElideMode.ElideRight, 180))


def _row(*widgets: QWidget) -> QWidget:
    """Widgets side by side in one form row."""
    box = QWidget()
    row = QHBoxLayout(box)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(8)
    for widget in widgets:
        row.addWidget(widget)
    row.addStretch(1)
    return box
