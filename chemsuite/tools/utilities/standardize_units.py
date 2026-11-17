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
"""Tools ▸ Standardize…: a value column and its unit column (nM, µM, mg/mL…) as one new
column, every row in the unit chosen."""

from __future__ import annotations

from PySide6.QtWidgets import QComboBox

from chemsuite.analysis.units import (
    MASS,
    MOLAR,
    convert,
    molecular_weights,
    needs_weight,
    parse_unit,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.table.columns import NUMBER, STRUCTURE, TEXT, NumberColumn
from chemsuite.core.wording import plural

#: words in a column name that mark it as the value or the unit column
VALUE_WORDS = ("ic50", "ec50", "ki", "kd", "potency", "value", "activity", "conc")
UNIT_WORDS = ("unit",)
DEFAULT_UNIT = "nM"

_last: str = DEFAULT_UNIT  # the unit chosen, remembered for the session


def _guess(names: list[str], words: tuple[str, ...]) -> str | None:
    return next((n for n in names if any(w in n.lower() for w in words)), None)


def output_name(values: str, unit: str) -> str:
    return f"{values} ({unit})"


class StandardizeUnitsDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Standardize", guide="standardize-values")
        self.ctx = ctx
        table = ctx.table
        numbers = [n for n in table.column_names if table.kind(n) in (NUMBER, TEXT)]
        texts = [n for n in table.column_names if table.kind(n) == TEXT]

        inputs = self.section("Input")
        self.values = QComboBox()
        self.values.addItems(numbers)
        guess = _guess([n for n in numbers if table.kind(n) == NUMBER], VALUE_WORDS) or _guess(
            numbers, VALUE_WORDS
        )
        if guess:
            self.values.setCurrentText(guess)
        inputs.row(
            "Values",
            self.values,
            help="Numbers, or text such as “>10”: the qualifier is kept in a column of its own.",
        )
        self.units = QComboBox()
        self.units.addItems(texts)
        guess = _guess(texts, UNIT_WORDS)
        if guess:
            self.units.setCurrentText(guess)
        inputs.row(
            "Units",
            self.units,
            help="Each row's unit: M, mM, µM (or uM), nM, pM, fM, g/L, mg/mL, mg/L, µg/mL, "
            "µg/L, ng/mL, ng/L or pg/mL, in any case.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        method = self.section("Method")
        self.target = QComboBox()
        for unit in MOLAR:
            self.target.addItem(unit, unit)
        self.target.insertSeparator(self.target.count())
        for unit in MASS:
            self.target.addItem(unit, unit)
        self.target.setCurrentIndex(max(0, self.target.findData(_last)))
        method.row("Convert to", self.target)
        self.structure = QComboBox()
        self.structure.addItems([n for n in table.column_names if table.kind(n) == STRUCTURE])
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        method.row(
            "Molecular weight from",
            self.structure,
            help="Molar and mass units convert through each row's molecular weight, worked out "
            "from its structure as drawn (a salt's counterion included).",
        )

        for signal in (
            self.values.currentIndexChanged,
            self.units.currentIndexChanged,
            self.target.currentIndexChanged,
            self.structure.currentIndexChanged,
            self.rows.changed,
        ):
            signal.connect(lambda *_: self._update())
        self._update()

    def _units_used(self) -> list[str | None]:
        column = self.units.currentText()
        if not column:
            return []
        return [parse_unit(t) for t in self.ctx.distinct_values(column, self.rows.oids())]

    def _needs_weight(self) -> bool:
        return needs_weight(self._units_used(), self.target.currentData())

    def _update(self) -> None:
        self.structure.setEnabled(self._needs_weight())
        self.refresh()

    def problem(self) -> str:
        if not self.values.currentText():
            return "There is no column of values."
        if not self.units.currentText():
            return "There is no text column of units."
        if self.rows.problem():
            return self.rows.problem()
        if not len(self.rows.oids()):
            return "There are no rows to work on."
        if not any(self._units_used()):
            return f"No units are read in “{self.units.currentText()}”."
        if self._needs_weight() and not self.structure.currentText():
            return "Molar and mass units are mixed: a structure column is needed."
        return ""

    def verb(self) -> str:
        return f"Standardize {len(self.rows.oids()):,} Rows"

    def run(self):
        global _last
        ctx = self.ctx
        _last = target = self.target.currentData()
        values_name, units_name = self.values.currentText(), self.units.currentText()
        oids = self.rows.oids()
        read = ctx.table.read([values_name, units_name], oids)
        values, units = read[values_name].tolist(), read[units_name].tolist()
        cells = ctx.structure_inputs(self.structure.currentText(), oids) if (
            self._needs_weight()) else None  # fmt: skip

        def work(_job):
            weights = molecular_weights(cells) if cells is not None else None
            return convert(values, units, target, weights)

        def done(result) -> None:
            name = output_name(values_name, target)
            results: dict = {name: (oids, NumberColumn.from_floats(result.values))}
            if any(result.qualifiers):
                results[f"{name} Qualifier"] = (oids, result.qualifiers)
            names = ctx.apply_results(results, label="Standardize values")
            made = int((result.values == result.values).sum())
            message = f"Wrote {plural(made, 'value')} in {target} to “{names[name]}”"
            if result.unknown_units:
                message += f"; {plural(result.unknown_units, 'unit')} not read"
            if result.no_weight:
                message += f"; {plural(result.no_weight, 'row')} with no molecular weight"
            ctx.show_message(message)

        return ctx.run("Standardizing values…", work, done)
