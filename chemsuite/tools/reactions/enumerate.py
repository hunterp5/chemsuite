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

"""Chemistry ▸ Enumerate Products."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QComboBox, QFileDialog

from chemsuite.analysis.reactions import (
    POOL_MIN_ROWS,
    Library,
    TemplateError,
    combinations_count,
    enumerate_library,
    enumerate_rows,
    parse_template,
)
from chemsuite.app.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.jobs import Cancelled
from chemsuite.core.mols import parse_structure
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn

from .common import ReactionText, cells, current_reaction

MAX_REACTANTS = 4
#: every-combination runs above this are refused: the result would not fit a table
MAX_COMBINATIONS = 2_000_000
FILE = "\x00file"  # the combo item that asks for a file
PRODUCT = "Product"
FOUND = "Products Found"
FILE_FILTER = (
    "Structure files (*.sdf *.sd *.mol *.csv *.tsv *.txt *.smi *.smiles *.xlsx *.gz);;All files (*)"
)


class EnumerateDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Enumerate Products", guide="enumerate-products")
        self.ctx = ctx
        self.files: dict[int, Path] = {}  # reactant slot -> building block file

        template = self.section("Template")
        reactions = ctx.reaction_columns()
        current = current_reaction(ctx, lambda: reactions[0]) if reactions else None
        self.template = template.row(
            "Reaction",
            ReactionText("Reaction SMARTS, e.g. [C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]", current),
            help="Atom-mapped reaction SMARTS, or Draw… and map the atoms with the Atom Map "
            "tool. Each reactant of the template takes one building block.",
        )

        blocks = self.section("Building Blocks")
        self._blocks_section = blocks
        self.sources: list[QComboBox] = []
        for i in range(MAX_REACTANTS):
            combo = QComboBox()
            for name in ctx.structure_columns():
                combo.addItem(name, name)
            combo.addItem("From a File…", FILE)
            if i < combo.count() - 1:
                combo.setCurrentIndex(i)
            combo.activated.connect(lambda _index, slot=i: self._source_chosen(slot))
            combo.currentIndexChanged.connect(lambda _index: self.refresh())
            self.sources.append(blocks.row(f"Reactant {i + 1}", combo))
        self.rows = blocks.row("Rows", RowsChoice(ctx))
        self.mode = blocks.row(
            "Combine",
            Segmented((("rows", "Row by row"), ("all", "Every combination")), "all"),
            help="Row by row: reactant 1 from one column and reactant 2 from another in the "
            "same row; the product goes in a new column. Every combination: each distinct "
            "building block with every other (a virtual library), in a new window.",
        )

        self.template.changed.connect(self.refresh)
        for signal in (self.rows.changed, self.mode.changed):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    # -- the template and its building blocks -------------------------------------------------
    def _template(self):
        try:
            return parse_template(self.template.text()), ""
        except TemplateError as err:
            return None, str(err)

    def reactant_count(self) -> int:
        rxn, _ = self._template()
        return rxn.GetNumReactantTemplates() if rxn is not None else 1

    def _source_chosen(self, slot: int) -> None:
        combo = self.sources[slot]
        if combo.currentData() != FILE:
            self.files.pop(slot, None)
            self.refresh()
            return
        path, _ = QFileDialog.getOpenFileName(self, "Building Blocks", "", FILE_FILTER)
        if path:
            self.set_file(slot, Path(path))
        else:
            combo.setCurrentIndex(0)
            self.files.pop(slot, None)
        self.refresh()

    def set_file(self, slot: int, path: Path) -> None:
        """Take reactant *slot*'s building blocks from the structures in *path*."""
        self.files[slot] = path
        combo = self.sources[slot]
        combo.setItemText(combo.count() - 1, f"File: {path.name}")
        combo.setCurrentIndex(combo.count() - 1)

    def _slots(self) -> range:
        return range(min(self.reactant_count(), MAX_REACTANTS))

    def _uses_files(self) -> bool:
        return any(self.sources[i].currentData() == FILE for i in self._slots())

    def _row_mode(self) -> bool:
        return self.mode.value() == "rows"

    def _distinct(self, slot: int) -> list:
        """The distinct structures of reactant *slot*'s column among the chosen rows."""
        name = self.sources[slot].currentData()
        oids = self.rows.oids()
        seen: dict[str, None] = {}
        for text in self.ctx.table.column(name).displays(self.ctx.table.slots_of(oids)):
            if text:
                seen.setdefault(text, None)
        return list(seen)

    def combinations(self) -> int | None:
        """How many combinations an every-combination run tries (None when a file is used)."""
        if self._uses_files():
            return None
        return combinations_count([len(self._distinct(i)) for i in self._slots()])

    # -- ToolDialog ---------------------------------------------------------------------------
    def refresh(self) -> None:
        if hasattr(self, "_blocks_section"):
            count = self.reactant_count()
            for i, combo in enumerate(self.sources):
                self._blocks_section.set_row_visible(combo, i < count)
        super().refresh()

    def problem(self) -> str:
        rxn, error = self._template()
        if rxn is None:
            return error
        if rxn.GetNumReactantTemplates() > MAX_REACTANTS:
            return f"Templates with more than {MAX_REACTANTS} reactants are not supported."
        if not len(self.rows.oids()) and not all(
            self.sources[i].currentData() == FILE for i in self._slots()
        ):
            return "There are no rows to take building blocks from."
        if self._row_mode() and self._uses_files():
            return "Row by row takes every building block from the table: choose columns."
        count = None if self._row_mode() else self.combinations()
        if count is not None and count > MAX_COMBINATIONS:
            return (
                f"{count:,} combinations is more than {MAX_COMBINATIONS:,}: choose fewer "
                "rows (Selected or Shown)."
            )
        return ""

    def verb(self) -> str:
        if self._row_mode():
            return f"Enumerate {len(self.rows.oids()):,} Rows"
        count = self.combinations()
        return f"Enumerate {count:,} Combinations" if count is not None else "Enumerate"

    def run(self):
        template = self.template.text()
        if self._row_mode():
            return self._run_rows(template)
        return self._run_library(template)

    def _run_rows(self, template: str):
        ctx = self.ctx
        oids = self.rows.oids()
        columns = [cells(ctx.table, self.sources[i].currentData(), oids) for i in self._slots()]
        use_pool = len(oids) >= POOL_MIN_ROWS

        def work(job):
            return enumerate_rows(
                template,
                columns,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(products: list[list[str]]) -> None:
            finish_rows(ctx, oids[: len(products)], products)

        return ctx.run("Enumerating products…", work, done, keep_partial=True)

    def _run_library(self, template: str):
        ctx = self.ctx
        slots = list(self._slots())
        sources: list[list | Path] = [
            self.files[i] if self.sources[i].currentData() == FILE else self._distinct(i)
            for i in slots
        ]

        def work(job):
            blocks = [load_blocks(s) if isinstance(s, Path) else s for s in sources]
            total = combinations_count([len(b) for b in blocks])
            if total > MAX_COMBINATIONS:
                raise TemplateError(
                    f"{total:,} combinations is more than {MAX_COMBINATIONS:,}: use fewer "
                    "building blocks."
                )
            try:
                library = enumerate_library(
                    template,
                    blocks,
                    executor=ctx.process_pool() if total >= POOL_MIN_ROWS else None,
                    progress=job.report,
                    cancel=job.cancel_token,
                )
            except Cancelled as stopped:  # keep the products of the combinations finished
                raise Cancelled((blocks, stopped.partial) if stopped.partial else None) from None
            return blocks, library

        def done(result: tuple[list[list], Library]) -> None:
            finish_library(ctx, *result)

        return ctx.run("Enumerating a library…", work, done, keep_partial=True)


def load_blocks(path: Path) -> list[str]:
    """The distinct structures of the first structure column in *path*, as text."""
    from chemsuite.core.io.load import load_table

    table = load_table(path)
    name = next(
        (n for n in table.column_names if isinstance(table.column(n), StructureColumn)), None
    )
    if name is None:
        raise TemplateError(f"{path.name} has no structure column.")
    seen: dict[str, None] = {}
    for text in table.column(name).displays(np.arange(len(table))):
        if text:
            seen.setdefault(text, None)
    return list(seen)


def finish_rows(ctx: ToolContext, oids: np.ndarray, products: list[list[str]]) -> dict[str, str]:
    """The first product of each row in a structure column, and how many there were."""
    texts = [p[0] if p else "" for p in products]
    blobs = [parse_structure(t)[0] if t else b"" for t in texts]
    written = ctx.apply_results(
        {
            PRODUCT: (oids, StructureColumn.from_parsed(blobs, texts)),
            FOUND: (oids, [float(len(p)) for p in products]),
        },
        label="Enumerate Products",
    )
    made = sum(1 for p in products if p)
    several = sum(1 for p in products if len(p) > 1)
    message = f"Products for {made:,} of {len(oids):,} row(s)"
    if several:
        message += f"; {several:,} with more than one (the first is shown)"
    ctx.show_message(message, 8000)
    return written


def library_table(blocks: list[list], library: Library) -> Table:
    """One row per product: the building blocks it came from and the product."""
    data: dict[str, object] = {}
    for i, block in enumerate(blocks):
        texts = [block[k] if isinstance(block[k], str) else "" for k in library.blocks[i]]
        data[f"Reactant {i + 1}"] = StructureColumn.from_texts(texts)
    blobs = [parse_structure(t)[0] for t in library.products]
    data[PRODUCT] = StructureColumn.from_parsed(blobs, library.products)
    table = Table()
    table.add_rows(data, undoable=False)
    table.history.mark_clean()
    return table


def finish_library(ctx: ToolContext, blocks: list[list], library: Library) -> None:
    if not library.products:
        ctx.show_message(
            f"None of the {library.combinations:,} combinations gave a product: check that "
            "each building block has the group its reactant needs",
            10000,
        )
        return
    ctx.open_table(
        library_table(blocks, library),
        f"{len(library.products):,} products from {library.combinations:,} combinations",
    )
    ctx.show_message(
        f"Enumerated {len(library.products):,} product(s) from {library.combinations:,} "
        f"combination(s); {library.failed:,} gave none",
        10000,
    )
