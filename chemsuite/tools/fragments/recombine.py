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

"""Chemistry ▸ Enumerate ▸ Recombine Fragments: new molecules from a BRICS fragment pool."""

from __future__ import annotations

from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox

from chemsuite.analysis.recombine import (
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_PRODUCTS,
    PRESETS,
    PROPERTIES,
    Limits,
    Recombined,
    recombine,
    usable_fragments,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog, inline
from chemsuite.core.mols import mol_from_text, parse_structure
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.wording import plural
from chemsuite.tools.reactions.enumerate import new_table

NO_LIMIT = -1.0  # a spin box at its minimum means "not checked"
PRODUCT = "Product"


def fragment_columns(ctx: ToolContext) -> list[str]:
    """The columns a BRICS run wrote (``BRICS 1``, ``BRICS 2``…)."""
    return [
        name
        for name in ctx.table.column_names
        if name.startswith("BRICS ") and name.split(" ")[-1].isdigit()
    ]


class RecombineDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Recombine Fragments", guide="recombine-fragments")
        self.ctx = ctx

        inputs = self.section("Input")
        self.columns = fragment_columns(ctx)
        self.rows = inputs.row(
            "Rows",
            RowsChoice(ctx),
            help="The fragment pool is taken from the BRICS columns of these rows.",
        )

        method = self.section("Method")
        self.max_products = QSpinBox()
        self.max_products.setRange(1, 100_000)
        self.max_products.setValue(DEFAULT_MAX_PRODUCTS)
        method.row("Most products", self.max_products)
        self.max_depth = QSpinBox()
        self.max_depth.setRange(1, 8)
        self.max_depth.setValue(DEFAULT_MAX_DEPTH)
        method.row(
            "Most fragments joined",
            self.max_depth,
            help="How many fragments one product may be built from.",
        )

        more = self.section("More Options", collapsed=True)
        self.preset = QComboBox()
        self.preset.addItems(list(PRESETS))
        self.preset.currentTextChanged.connect(self._preset_chosen)
        more.row(
            "Drug-likeness",
            self.preset,
            help="Rule of 5 (Lipinski 1997) or Lead-like (Teague 1999) fills the ranges "
            "below; changing a range makes it your own.",
        )
        self.limits: dict[str, tuple[QDoubleSpinBox, QDoubleSpinBox]] = {}
        for name, (_measure, label) in PROPERTIES.items():
            low, high = QDoubleSpinBox(), QDoubleSpinBox()
            for box in (low, high):
                box.setRange(NO_LIMIT, 10_000.0)
                box.setDecimals(1)
                box.setValue(NO_LIMIT)
                box.setSpecialValueText("—")  # the minimum means "not checked"
                box.valueChanged.connect(lambda *_: self.refresh())
            self.limits[name] = (low, high)
            more.row(name, inline(low, high, stretch=0), help=f"Least and most {label}.")
        self.novel = QCheckBox("Skip compounds already in the table")
        self.novel.setChecked(True)
        self.novel.toggled.connect(lambda *_: self.refresh())
        more.row("Novelty", self.novel)

        self.rows.changed.connect(lambda *_: self.refresh())
        for box in (self.max_products, self.max_depth):
            box.valueChanged.connect(lambda *_: self.refresh())
        self.refresh()

    # -- the pool and the limits ------------------------------------------------------------
    def _preset_chosen(self, name: str) -> None:
        chosen = PRESETS.get(name)
        if chosen is None:
            return
        for prop, (low, high) in self.limits.items():
            want = chosen.ranges.get(prop, (None, None))
            low.setValue(want[0] if want[0] is not None else NO_LIMIT)
            high.setValue(want[1] if want[1] is not None else NO_LIMIT)
        self.refresh()

    def chosen_limits(self) -> Limits:
        ranges: dict[str, tuple[float | None, float | None]] = {}
        for prop, (low, high) in self.limits.items():
            start = None if low.value() == NO_LIMIT else low.value()
            stop = None if high.value() == NO_LIMIT else high.value()
            if start is not None or stop is not None:
                ranges[prop] = (start, stop)
        return Limits(ranges)

    def pool(self) -> list[str]:
        """The distinct, usable fragments in the BRICS columns of the chosen rows."""
        oids = self.rows.oids()
        found: list[str] = []
        for name in self.columns:
            found.extend(self.ctx.distinct_values(name, oids))
        return usable_fragments(found)

    def existing(self) -> list[str]:
        """The structures already in the table, as written (canonicalised in the worker)."""
        if not self.novel.isChecked():
            return []
        column = self.ctx.default_structure_column()
        if not column:
            return []
        return self.ctx.distinct_values(column, self.ctx.scope_oids("all"))

    # -- ToolDialog ---------------------------------------------------------------------------
    def problem(self) -> str:
        if not self.columns:
            return "Run Fragments with BRICS first: its columns are the pool to build from."
        if not len(self.rows.oids()):
            return "There are no rows to take fragments from."
        found = self.pool()
        if len(found) < 2:
            return (
                "Fewer than two fragments have numbered attachment points. Run Fragments "
                "with BRICS and “Plain attachment points” off."
            )
        for prop, (low, high) in self.limits.items():
            start, stop = low.value(), high.value()
            if start != NO_LIMIT and stop != NO_LIMIT and start > stop:
                return f"{prop}: the least is more than the most."
        return ""

    def verb(self) -> str:
        count = len(self.pool()) if self.columns else 0
        return f"Recombine {plural(count, 'Fragment')}" if count else "Recombine"

    def run(self):
        ctx = self.ctx
        fragments = self.pool()
        limits = self.chosen_limits()
        existing = self.existing()
        cap = self.max_products.value()
        depth = self.max_depth.value()

        def work(job):
            known = canonical_smiles(existing)
            return recombine(
                fragments,
                max_depth=depth,
                max_products=cap,
                limits=limits,
                known=known,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(made: Recombined) -> None:
            finish(ctx, made, len(fragments))

        return ctx.run("Recombining fragments…", work, done, keep_partial=True)


def canonical_smiles(texts: list[str]) -> list[str]:
    """RDKit canonical SMILES of the readable *texts*, to compare with products."""
    from rdkit import Chem

    out = []
    for text in texts:
        mol = mol_from_text(text)
        if mol is not None:
            out.append(Chem.MolToSmiles(mol))
    return out


def products_table(made: Recombined) -> Table:
    """One row per product."""
    blobs = [parse_structure(text)[0] for text in made.products]
    return new_table({PRODUCT: StructureColumn.from_parsed(blobs, made.products)})


def finish(ctx: ToolContext, made: Recombined, pool: int) -> None:
    """Open the products in a new window, saying what became of the rest."""
    if not made.products:
        ctx.show_message(
            f"Recombining {plural(pool, 'fragment')} made nothing: "
            f"{made.filtered:,} candidates were outside the property limits"
            if made.filtered
            else f"Recombining {plural(pool, 'fragment')} made nothing"
        )
        return
    ctx.open_table(
        products_table(made),
        f"{len(made.products):,} molecules from {pool:,} fragments",
    )
    message = f"Recombined {plural(pool, 'fragment')} into {plural(len(made.products), 'molecule')}"
    extra = []
    if made.filtered:
        extra.append(f"{made.filtered:,} outside the limits")
    if made.known:
        extra.append(f"{made.known:,} already in the table")
    if extra:
        message += f"; {', '.join(extra)}"
    if made.stopped:
        message += "; the search gave up before the cap"
    ctx.show_message(message)
