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
"""The Bin Values dialog, with a live preview of the bins and their counts."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QHeaderView,
    QLineEdit,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)

from chemsuite.analysis.binning import (
    METHODS,
    Binned,
    BinSpec,
    bin_values,
    parse_names,
    parse_numbers,
)
from chemsuite.app.dialogs import FIELD_WIDTH, RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.plotting.data import number_columns, numbers


class BinningDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Bin Values", guide="bin-values")
        self.ctx = ctx
        self.result: Binned | None = None
        self._problem_text = ""

        inputs = self.section("Input")
        self.column = QComboBox()
        self.column.addItems(number_columns(ctx.table))
        current = ctx.current_column()
        if current and self.column.findText(current) >= 0:
            self.column.setCurrentText(current)
        inputs.row("Column", self.column)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.method = QComboBox()
        for key, title in METHODS.items():
            self.method.addItem(title, key)
        method.row("Method", self.method)
        self.count = QSpinBox()
        self.count.setRange(2, 256)
        self.count.setValue(3)
        method.row("Bins", self.count)
        self.edges = QLineEdit()
        self.edges.setPlaceholderText("5, 6, 7")
        method.row("Edges", self.edges, help="The bin edges, low to high, separated by commas.")
        self.open_ends = QCheckBox("Bins below the first and above the last edge")
        self.open_ends.setChecked(True)
        method.row("", self.open_ends)
        self.names = QLineEdit()
        self.names.setPlaceholderText("Low, Medium, High")
        method.row(
            "Class names",
            self.names,
            help="One name per bin, low to high. Empty: the classes are numbered.",
        )
        self._method_section = method

        more = self.section("More Options", collapsed=True)
        self.width = QDoubleSpinBox()
        self.width.setRange(0.0, 1e12)
        self.width.setDecimals(4)
        self.width.setSpecialValueText("Automatic")
        more.row(
            "Bin width",
            self.width,
            help="Automatic: the range divided by the number of bins. A width sets the "
            "number of bins itself.",
        )
        self.minimum = QLineEdit()
        self.minimum.setPlaceholderText("Lowest value")
        more.row("From", self.minimum, help="Where the first bin starts (empty: the data's).")
        self.maximum = QLineEdit()
        self.maximum.setPlaceholderText("Highest value")
        more.row("To", self.maximum, help="Where the last bin ends (empty: the data's).")
        self._more = more

        output = self.section("Output")
        self.write_range, self.range_name = self._named_output(
            output, "Range column", "A text column with each row's range, such as “[5, 6)”"
        )
        self.write_class, self.class_name = self._named_output(
            output, "Class column", "A text column with each row's class name"
        )
        self._output = output

        preview = self.section("Preview")
        self.preview = QTableWidget(0, 3)
        self.preview.setHorizontalHeaderLabels(["Range", "Class", "Rows"])
        self.preview.verticalHeader().hide()
        self.preview.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.preview.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.preview.setFixedSize(FIELD_WIDTH, 150)
        self.preview.setShowGrid(False)
        self.preview.setAlternatingRowColors(True)
        preview.row("Bins", self.preview)

        self.column.currentIndexChanged.connect(self._column_changed)
        for signal in (
            self.method.currentIndexChanged,
            self.count.valueChanged,
            self.width.valueChanged,
            self.minimum.textChanged,
            self.maximum.textChanged,
            self.edges.textChanged,
            self.open_ends.toggled,
            self.names.textChanged,
            self.rows.changed,
            self.write_range.toggled,
            self.write_class.toggled,
            self.range_name.textChanged,
            self.class_name.textChanged,
        ):
            signal.connect(lambda *_: self.update_preview())
        self._column_changed()

    @staticmethod
    def _named_output(section, label: str, help: str) -> tuple[QCheckBox, QLineEdit]:
        box = QCheckBox()
        box.setChecked(True)
        name = QLineEdit()
        field = QWidget()
        line = QHBoxLayout(field)
        line.setContentsMargins(0, 0, 0, 0)
        line.setSpacing(6)
        line.addWidget(box)
        line.addWidget(name, 1)
        field.setFixedWidth(FIELD_WIDTH)
        box.toggled.connect(name.setEnabled)
        section.row(label, field, help=help)
        return box, name

    def _column_changed(self, *_args) -> None:
        name = self.column.currentText()
        self.range_name.setText(f"{name} Range" if name else "")
        self.class_name.setText(f"{name} Class" if name else "")
        self.update_preview()

    def _values(self) -> np.ndarray:
        name = self.column.currentText()
        if not name:
            return np.empty(0)
        return numbers(self.ctx.table, name, self.rows.oids())

    @staticmethod
    def _limit(edit: QLineEdit) -> float | None:
        text = edit.text().strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            raise ValueError(f"“{text}” is not a number.") from None

    def spec(self) -> BinSpec:
        method = self.method.currentData()
        names = parse_names(self.names.text())
        if method == "custom":
            return BinSpec(
                "custom",
                edges=parse_numbers(self.edges.text()),
                open_ends=self.open_ends.isChecked(),
                names=names,
            )
        if method == "equal_count":
            return BinSpec("equal_count", self.count.value(), names=names)
        return BinSpec(
            "equal_width",
            self.count.value(),
            width=self.width.value() or None,
            minimum=self._limit(self.minimum),
            maximum=self._limit(self.maximum),
            names=names,
        )

    def update_preview(self) -> None:
        method = self.method.currentData()
        width = method == "equal_width"
        section = self._method_section
        section.set_row_visible(self.count, method != "custom")
        section.set_row_visible(self.edges, method == "custom")
        section.set_row_visible(self.open_ends, method == "custom")
        self._more.setVisible(width)
        self.count.setEnabled(not (width and self.width.value()))
        self.result = None
        self.preview.setRowCount(0)
        self._problem_text = ""
        try:
            result = bin_values(self._values(), self.spec())
        except ValueError as err:
            self._problem_text = str(err)
            self.refresh()
            return
        self.result = result
        counts = result.counts()
        self.preview.setRowCount(len(result.intervals))
        for row, (rng, cls, count) in enumerate(
            zip(result.range_names, result.class_names, counts.tolist(), strict=True)
        ):
            for col, text in enumerate((rng, cls, f"{count:,}")):
                self.preview.setItem(row, col, QTableWidgetItem(text))
        self.refresh()

    # -- what will happen ---------------------------------------------------------------------
    def problem(self) -> str:
        if not self.column.currentText():
            return "There is no number column to bin."
        if self._problem_text:
            return self._problem_text
        if not (self.write_range.isChecked() or self.write_class.isChecked()):
            return "Choose at least one column to write."
        return ""

    def verb(self) -> str:
        return f"Bin {len(self.rows.oids()):,} Rows"

    def run(self):
        if self.result is not None:
            self.apply()
        return None

    def apply(self) -> dict[str, str]:
        result = self.result
        assert result is not None
        oids = self.rows.oids()
        source = self.column.currentText()
        results = {}
        if self.write_range.isChecked():
            name = self.range_name.text().strip() or f"{source} Range"
            results[name] = (oids, result.labels(result.range_names))
        if self.write_class.isChecked():
            name = self.class_name.text().strip() or f"{source} Class"
            results[name] = (oids, result.labels(result.class_names))
        names = self.ctx.apply_results(results, label=f"Bin {source}")
        self.ctx.show_message(f"Binned {source} into {len(result.intervals)} bin(s)")
        return names
