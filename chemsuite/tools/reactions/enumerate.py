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

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QComboBox, QFileDialog, QPushButton

from chemsuite.analysis.reaction_templates import any_fits
from chemsuite.analysis.reactions import (
    POOL_MIN_ROWS,
    Library,
    TemplateError,
    combinations_count,
    enumerate_library,
    enumerate_rows,
    parse_template,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog
from chemsuite.core.jobs import Cancelled
from chemsuite.core.mols import mol_from_text, parse_structure
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.wording import plural

from .common import ReactionText, TemplateChoice, ask_structure, cells, current_reaction

MAX_REACTANTS = 4
#: every-combination runs above this are refused: the result would not fit a table
MAX_COMBINATIONS = 2_000_000
FILE = "\x00file"  # the combo item that asks for a file
SINGLE = "\x00single"  # the combo item that asks for one typed or drawn structure
#: building blocks parsed when checking that a component's source fits its role
FIT_SAMPLE = 100
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
        self.singles: dict[int, str] = {}  # reactant slot -> one typed or drawn structure
        self._filling = False  # a named template is writing the text, so leave the choice

        template = self.section("Template")
        self.choice = template.row(
            "Reaction",
            TemplateChoice(),
            help="A named reaction, or Custom to type or draw your own. Each component of "
            "the reaction takes one building block.",
        )
        reactions = ctx.reaction_columns()
        current = current_reaction(ctx, lambda: reactions[0]) if reactions else None
        self.template = template.row(
            "Template",
            ReactionText("Reaction SMARTS, e.g. [C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]", current),
            help="Atom-mapped reaction SMARTS, or Draw… and map the atoms with the Atom Map "
            "tool. Editing a named reaction's template makes it Custom.",
        )

        blocks = self.section("Building Blocks")
        self._blocks_section = blocks
        self.sources: list[QComboBox] = []
        for i in range(MAX_REACTANTS):
            combo = QComboBox()
            for name in ctx.structure_columns():
                combo.addItem(name, name)
            combo.addItem("From a File…", FILE)
            combo.addItem("A Single Structure…", SINGLE)
            if i < combo.count() - 2:
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

        more = self.section("More Options", collapsed=True)
        self.swap = QPushButton("Swap Components")
        self.swap.setAutoDefault(False)
        self.swap.clicked.connect(self._swap)
        more.row(
            "Order",
            self.swap,
            help="Exchange the first two components, when the building blocks are the other "
            "way round.",
        )

        self.choice.picked.connect(self._template_picked)
        self.template.changed.connect(self._text_edited)
        for signal in (self.rows.changed, self.mode.changed):
            signal.connect(lambda *_: self.refresh())
        self.refresh()

    # -- the named template ---------------------------------------------------------------------
    def _template_picked(self) -> None:
        """A named reaction was chosen: write its SMARTS into the template field."""
        template = self.choice.template()
        if template is None:
            return
        self._filling = True
        try:
            self.template.set_text(template.smarts)
        finally:
            self._filling = False
        self.refresh()

    def _text_edited(self) -> None:
        """The template text changed; editing a named reaction's SMARTS makes it Custom."""
        if not self._filling:
            template = self.choice.template()
            if template is not None and self.template.text() != template.smarts:
                self.choice.show_custom()
        self.refresh()

    def role(self, slot: int) -> str:
        """What this component is called: the template's role, else "Reactant n"."""
        template = self.choice.template()
        if template is not None and slot < len(template.roles):
            return template.roles[slot]
        return f"Reactant {slot + 1}"

    def _swap(self) -> None:
        """Exchange the first two components' sources."""
        first, second = self.sources[0], self.sources[1]
        index = first.currentIndex()
        first.setCurrentIndex(second.currentIndex())
        second.setCurrentIndex(index)
        for store in (self.files, self.singles):
            held = (store.pop(0, None), store.pop(1, None))
            if held[1] is not None:
                store[0] = held[1]
            if held[0] is not None:
                store[1] = held[0]
        for slot in (0, 1):
            self._relabel_source(slot)
        self.refresh()

    def _relabel_source(self, slot: int) -> None:
        """Keep the File and Single items of *slot*'s combo showing what they hold."""
        combo = self.sources[slot]
        path = self.files.get(slot)
        single = self.singles.get(slot)
        combo.setItemText(combo.findData(FILE), f"File: {path.name}" if path else "From a File…")
        combo.setItemText(
            combo.findData(SINGLE),
            f"Structure: {single}" if single else "A Single Structure…",
        )

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
        chosen = combo.currentData()
        if chosen == FILE:
            path, _ = QFileDialog.getOpenFileName(self, self.role(slot), "", FILE_FILTER)
            if path:
                self.singles.pop(slot, None)
                self.set_file(slot, Path(path))
            else:
                combo.setCurrentIndex(0)
                self.files.pop(slot, None)
        elif chosen == SINGLE:
            text = ask_structure(self, self.role(slot), self.singles.get(slot, ""))
            if text:
                self.files.pop(slot, None)
                self.set_single(slot, text)
            else:
                combo.setCurrentIndex(0)
                self.singles.pop(slot, None)
        else:
            self.files.pop(slot, None)
            self.singles.pop(slot, None)
        self.refresh()

    def set_file(self, slot: int, path: Path) -> None:
        """Take reactant *slot*'s building blocks from the structures in *path*."""
        self.files[slot] = path
        self._relabel_source(slot)
        self.sources[slot].setCurrentIndex(self.sources[slot].findData(FILE))

    def set_single(self, slot: int, text: str) -> None:
        """Use one structure for reactant *slot*, against every other component."""
        self.singles[slot] = text
        self._relabel_source(slot)
        self.sources[slot].setCurrentIndex(self.sources[slot].findData(SINGLE))

    def _slots(self) -> range:
        return range(min(self.reactant_count(), MAX_REACTANTS))

    def _uses_files(self) -> bool:
        return any(self.sources[i].currentData() == FILE for i in self._slots())

    def _row_mode(self) -> bool:
        return self.mode.value() == "rows"

    def _distinct(self, slot: int) -> list:
        """The distinct building blocks of reactant *slot*: one typed structure, or the
        distinct structures of its column among the chosen rows."""
        chosen = self.sources[slot].currentData()
        if chosen == SINGLE:
            single = self.singles.get(slot)
            return [single] if single else []
        oids = self.rows.oids()
        seen: dict[str, None] = {}
        for text in self.ctx.table.column(chosen).displays(self.ctx.table.slots_of(oids)):
            if text:
                seen.setdefault(text, None)
        return list(seen)

    def combinations(self) -> int | None:
        """How many combinations an every-combination run tries (None when a file is used)."""
        if self._uses_files():
            return None
        return combinations_count([len(self._distinct(i)) for i in self._slots()])

    def _misfit(self) -> str:
        """The first component whose building blocks cannot fill its role, if any. Only a
        named template knows its roles; a file is not read until the run."""
        template = self.choice.template()
        if template is None:
            return ""
        for slot in self._slots():
            if self.sources[slot].currentData() == FILE:
                continue
            texts = self._distinct(slot)[:FIT_SAMPLE]
            mols = [mol_from_text(text) for text in texts]
            if mols and not any_fits(template.id, slot, mols):
                where = self.sources[slot].currentText().strip()
                return f"Nothing in {where} is a {self.role(slot).lower()}."
        return ""

    # -- ToolDialog ---------------------------------------------------------------------------
    def refresh(self) -> None:
        if hasattr(self, "_blocks_section"):
            count = self.reactant_count()
            for i, combo in enumerate(self.sources):
                self._blocks_section.set_row_visible(combo, i < count)
                if i < count:
                    self._blocks_section.caption(combo).setText(self.role(i))
            self.swap.setEnabled(count > 1)
        super().refresh()

    def problem(self) -> str:
        rxn, error = self._template()
        if rxn is None:
            return error
        if rxn.GetNumReactantTemplates() > MAX_REACTANTS:
            return f"Templates with more than {MAX_REACTANTS} reactants are not supported."
        if not len(self.rows.oids()) and not all(
            self.sources[i].currentData() in (FILE, SINGLE) for i in self._slots()
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
        return self._misfit()

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

    def _row_cells(self, slot: int, oids) -> list:
        """Reactant *slot* for every chosen row: its column, or one structure in every row."""
        chosen = self.sources[slot].currentData()
        if chosen == SINGLE:
            return [self.singles.get(slot) or None] * len(oids)
        return cells(self.ctx.table, chosen, oids)

    def _run_rows(self, template: str):
        ctx = self.ctx
        oids = self.rows.oids()
        columns = [self._row_cells(i, oids) for i in self._slots()]
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
        roles = [self.role(i) for i in slots]

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
            finish_library(ctx, *result, names=roles)

        return ctx.run("Enumerating a library…", work, done, keep_partial=True)


def load_blocks(path: Path) -> list[str]:
    """The distinct structures of the first structure column in *path*, as text."""
    from chemsuite.io.load import load_table

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
    message = f"Enumerated products for {made:,} of {plural(len(oids), 'row')}"
    if several:
        message += f"; {several:,} with more than one (the first is shown)"
    ctx.show_message(message)
    return written


def library_table(
    blocks: list[list], library: Library, names: Sequence[str] | None = None
) -> Table:
    """One row per product: the building blocks it came from and the product. *names* titles
    the building-block columns (a named template's roles); without it they are numbered."""
    data: dict[str, object] = {}
    for i, block in enumerate(blocks):
        texts = [block[k] if isinstance(block[k], str) else "" for k in library.blocks[i]]
        title = names[i] if names is not None and i < len(names) else f"Reactant {i + 1}"
        data[title if title not in data else f"{title} {i + 1}"] = StructureColumn.from_texts(texts)
    blobs = [parse_structure(t)[0] for t in library.products]
    data[PRODUCT] = StructureColumn.from_parsed(blobs, library.products)
    table = Table()
    table.add_rows(data, undoable=False)
    table.history.mark_clean()
    return table


def finish_library(
    ctx: ToolContext,
    blocks: list[list],
    library: Library,
    names: Sequence[str] | None = None,
) -> None:
    if not library.products:
        ctx.show_message(
            f"Could not enumerate a product from {plural(library.combinations, 'combination')}:"
            " check that each building block has the group its reactant needs"
        )
        return
    ctx.open_table(
        library_table(blocks, library, names),
        f"{len(library.products):,} products from {library.combinations:,} combinations",
    )
    ctx.show_message(
        f"Enumerated {plural(len(library.products), 'product')} from "
        f"{plural(library.combinations, 'combination')}; {library.failed:,} gave none"
    )
