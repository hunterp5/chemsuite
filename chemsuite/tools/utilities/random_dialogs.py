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

"""Tools ▸ Table ▸ Random Numbers… and Random Rows…"""

from __future__ import annotations

from PySide6.QtGui import QIntValidator
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox

from chemsuite.analysis.random_values import (
    DISTRIBUTIONS,
    INTEGER,
    LOGNORMAL,
    NORMAL,
    ORDER,
    UNIFORM,
    RandomSpec,
    random_numbers,
    random_sample,
    sample_size,
)
from chemsuite.app.dialogs import OutputChoice, RowsChoice, Segmented, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.lists import ListError


def _number(value: float, decimals: int = 4) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(-1e12, 1e12)
    spin.setDecimals(decimals)
    spin.setValue(value)
    return spin


def _seed_field() -> QLineEdit:
    field = QLineEdit()
    field.setValidator(QIntValidator(0, 2_147_483_647))
    field.setPlaceholderText("None: a new draw each time")
    return field


def _seed(field: QLineEdit) -> int | None:
    text = field.text().strip()
    return int(text) if text else None


class RandomNumbersDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Random Numbers", guide="random-numbers")
        self.ctx = ctx

        method = self.section("Method")
        self.distribution = QComboBox()
        for key, title in DISTRIBUTIONS.items():
            self.distribution.addItem(title, key)
        method.row("Distribution", self.distribution)
        self.low = method.row("Low", _number(0.0), "The smallest value.")
        self.high = method.row("High", _number(1.0), "The largest value.")
        self.mean = method.row("Mean", _number(0.0))
        self.sd = method.row("Standard deviation", _number(1.0))
        self.sd.setMinimum(0)
        self.clip = QCheckBox("Keep within low and high")
        method.row("", self.clip, "Values drawn outside the range are set to low or high.")
        self.decimals = QSpinBox()
        self.decimals.setRange(0, 12)
        self.decimals.setValue(3)
        method.row("Decimals", self.decimals, "Values are rounded to this many decimals.")
        self.seed = method.row(
            "Seed",
            _seed_field(),
            "The same seed and settings give the same numbers again (for example to repeat a "
            "training / test split).",
        )

        rows = self.section("Rows")
        self.rows = RowsChoice(ctx, "all")
        rows.row("Fill", self.rows, "Rows outside the choice are left blank.")
        output = self.section("Output")
        self.name = QLineEdit("Random")
        output.row("Column", self.name)
        self.output = OutputChoice()
        output.row("If it exists", self.output)
        self.method = method

        self.distribution.currentIndexChanged.connect(self.refresh)
        for widget in (self.low, self.high, self.mean, self.sd, self.decimals):
            widget.valueChanged.connect(self.refresh)
        self.clip.toggled.connect(self.refresh)
        self.rows.changed.connect(lambda *_: self.refresh())
        self.name.textChanged.connect(self.refresh)
        self.refresh()

    def spec(self) -> RandomSpec:
        kind = self.distribution.currentData()
        return RandomSpec(
            kind,
            low=self.low.value(),
            high=self.high.value(),
            mean=self.mean.value(),
            sd=self.sd.value(),
            clip=self.clip.isChecked(),
            decimals=None if kind in (INTEGER, ORDER) else self.decimals.value(),
            seed=_seed(self.seed),
        )

    def refresh(self) -> None:
        if not hasattr(self, "method"):
            return
        kind = self.distribution.currentData()
        shaped = kind in (NORMAL, LOGNORMAL)
        ranged = kind in (UNIFORM, INTEGER) or (shaped and self.clip.isChecked())
        form = self.method
        form.set_row_visible(self.low, ranged)
        form.set_row_visible(self.high, ranged)
        form.set_row_visible(self.mean, shaped)
        form.set_row_visible(self.sd, shaped)
        form.set_row_visible(self.clip, shaped)
        form.set_row_visible(self.decimals, kind in (UNIFORM, NORMAL, LOGNORMAL))
        lognormal = kind == LOGNORMAL
        self.method.form.labelForField(self.mean).setText("Median:" if lognormal else "Mean:")
        self.method.form.labelForField(self.sd).setText(
            "Spread factor:" if lognormal else "Standard deviation:"
        )
        form.set_help(
            self.sd,
            "About two thirds of values fall between median ÷ factor and median × factor."
            if lognormal
            else "About two thirds of values fall within one standard deviation of the mean.",
        )
        super().refresh()

    def problem(self) -> str:
        if not self.name.text().strip():
            return "Name the new column."
        return self.spec().problem()

    def verb(self) -> str:
        return f"Fill {len(self.rows.oids()):,} Rows"

    def run(self):
        oids = self.rows.oids()
        values = random_numbers(self.spec(), len(oids))
        name = self.name.text().strip()
        names = self.ctx.apply_results(
            {name: (oids, [float(v) for v in values])},
            label=f"Random numbers in {name}",
            overwrite=self.output.overwrite(),
        )
        self.ctx.show_message(f"Filled {names.get(name, name)} with random numbers")
        return None


SELECT, ADD_SELECT, LIST = "select", "add", "list"


class RandomRowsDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Random Rows", guide="random-rows")
        self.ctx = ctx

        source = self.section("Rows")
        self.rows = RowsChoice(ctx, "shown")
        source.row("Pick from", self.rows)

        sample = self.section("Sample")
        self.unit = Segmented((("count", "Number of rows"), ("percent", "Percent")), "count")
        sample.row("Size as", self.unit)
        self.amount = QDoubleSpinBox()
        self.amount.setDecimals(1)
        self.amount.setRange(0, 1e9)
        self.amount.setValue(10)
        sample.row("Size", self.amount)
        self.seed = sample.row(
            "Seed", _seed_field(), "The same seed picks the same rows again from the same rows."
        )

        output = self.section("Output")
        self.action = Segmented(
            ((SELECT, "Select them"), (ADD_SELECT, "Add to selection"), (LIST, "New list")),
            SELECT,
        )
        output.row("Then", self.action)
        self.list_name = QLineEdit("Random sample")
        output.row("List name", self.list_name, "The rows are saved as a named list.")
        self.output_form = output

        self.unit.changed.connect(self._unit_changed)
        self.amount.valueChanged.connect(self.refresh)
        self.rows.changed.connect(lambda *_: self.refresh())
        self.action.changed.connect(lambda *_: self.refresh())
        self.list_name.textChanged.connect(self.refresh)
        self.refresh()

    def _unit_changed(self, unit: str) -> None:
        self.amount.setSuffix(" %" if unit == "percent" else "")
        self.amount.setMaximum(100 if unit == "percent" else 1e9)
        self.refresh()

    def count(self) -> int:
        return sample_size(
            len(self.rows.oids()), self.amount.value(), self.unit.value() == "percent"
        )

    def refresh(self) -> None:
        if not hasattr(self, "output_form"):
            return
        self.output_form.set_row_visible(self.list_name, self.action.value() == LIST)
        super().refresh()

    def problem(self) -> str:
        if not len(self.rows.oids()):
            return "There are no rows to pick from."
        if self.count() == 0:
            return "The sample is empty; choose a larger size."
        if self.action.value() == LIST:
            name = self.list_name.text().strip()
            if not name:
                return "Name the list."
            if name in self.ctx.lists:
                return f"A list named “{name}” already exists."
        return ""

    def verb(self) -> str:
        return f"Pick {self.count():,} Rows"

    def run(self):
        chosen = random_sample(self.rows.oids(), self.count(), _seed(self.seed))
        action = self.action.value()
        if action == LIST:
            name = self.list_name.text().strip()
            try:
                self.ctx.lists.add(name, chosen)
            except ListError as err:
                self.ctx.show_message(str(err))
                return None
            self.ctx.show_message(f"Saved {len(chosen):,} random rows as the list “{name}”")
        elif action == ADD_SELECT:
            self.ctx.selection.add(chosen.tolist(), "random rows")
            self.ctx.show_message(f"Added {len(chosen):,} random rows to the selection")
        else:
            self.ctx.selection.set(chosen.tolist(), "random rows")
            self.ctx.show_message(f"Selected {len(chosen):,} random rows")
        return None
