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
"""Edit ▸ Structure Search…: substructure, exact and similarity search over a structure column.

:class:`StructureSearcher` runs searches as background jobs for both this dialog and the
structure card in the Filters panel. The matches can be shown alone (as a filter card), selected
or added to a list, and substructure matches can be highlighted in the table's drawings.
"""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import Executor
from typing import TYPE_CHECKING, Any

import numpy as np
from PySide6.QtCore import QEvent, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QToolButton,
    QWidget,
)

from chemsuite.analysis.fingerprints import DEFAULT_KIND, KINDS, fingerprint_task
from chemsuite.analysis.structure_search import (
    MODES,
    QueryError,
    SearchResult,
    StructureQuery,
    drawn_query_text,
    parse_query,
    search,
)
from chemsuite.analysis.structure_search import (
    POOL_MIN_ROWS as SEARCH_POOL_MIN_ROWS,
)
from chemsuite.core.lists import ListError
from chemsuite.core.table import Table
from chemsuite.core.table.columns import STRUCTURE

from ..shell.jobs import QUICK, Job
from ..widgets.dialogs import FIELD_WIDTH, RowsChoice, Segmented, ToolDialog
from .table_view import Highlight

if TYPE_CHECKING:
    from ..document.main_window import MainWindow

Run = Callable[..., Job]
#: ``on_done(oids searched, result)``
SearchDone = Callable[[np.ndarray, SearchResult], None]


def structure_columns(table: Table) -> list[str]:
    return [name for name in table.column_names if table.kind(name) == STRUCTURE]


def highlight_for(query: StructureQuery) -> Highlight | None:
    """What to highlight for *query*: substructure matches only (an exact or similar hit has
    no atoms to single out)."""
    if query.mode != "substructure":
        return None
    try:
        parsed = parse_query(query)
    except QueryError:
        return None
    return Highlight(f"{query.text}|{query.stereo}", parsed.mol, query.stereo)


class StructureSearcher:
    """Runs structure searches over one table as background jobs.

    *run* is the window's ``(label, work, on_done, lane=…, guard=…) -> Job``; *pool* returns the
    shared process pool; *highlight* shows atoms matching a query in a column's drawings.
    """

    def __init__(
        self,
        table: Table,
        run: Run,
        pool: Callable[[], Executor],
        highlight: Callable[[str, Highlight | None], None] = lambda _c, _h: None,
    ) -> None:
        self.table = table
        self._run = run
        self._pool = pool
        self.highlight = highlight

    def start(
        self,
        column: str,
        query: StructureQuery,
        on_done: SearchDone,
        *,
        oids: np.ndarray | None = None,
    ) -> Job:
        """Search *oids* (default: every row). Raises :class:`QueryError` before starting."""
        parse_query(query)
        table = self.table
        oids = table.oids.copy() if oids is None else np.asarray(oids, dtype=np.int64)
        compute, remember = fingerprint_task(table, column, query.screen_kind, oids, self._pool)
        structures = table.column(column)
        slots = table.slots_of(oids)
        blobs = structures.blobs(slots)
        texts = structures.displays(slots)

        def cells_at(rows: np.ndarray) -> list:
            # A cell not parsed yet has no blob: the worker parses its text.
            return [blobs[i] if blobs[i] is not None else (texts[i] or None) for i in rows.tolist()]

        # Too few rows to ever need workers: do not start the process pool for them.
        pool = self._pool if len(oids) >= SEARCH_POOL_MIN_ROWS else (lambda: None)

        def work(job: Job) -> tuple[Any, SearchResult]:
            fps = compute(job)
            result = search(
                query,
                cells_at,
                fps,
                executor=pool(),
                progress=job.report,
                cancel=job.cancel_token,
            )
            return fps, result

        def done(value: tuple[Any, SearchResult]) -> None:
            fps, result = value
            remember(fps)
            on_done(oids, result)

        return self._run(
            f"Structure search: {query.describe()}", work, done, lane=QUICK, guard=True
        )

    def start_reactions(
        self,
        column: str,
        query,
        on_done: Callable[[np.ndarray, np.ndarray], None],
    ) -> Job:
        """Search reaction column *column* for a
        :class:`~chemsuite.analysis.reactions.ReactionQuery`; *on_done* gets ``(oids, hits)``.
        Raises :class:`QueryError` before starting."""
        from chemsuite.analysis.reactions import (
            POOL_MIN_ROWS,
            parse_reaction_query,
            search_reactions,
        )

        parse_reaction_query(query)
        table = self.table
        oids = table.oids.copy()
        data = table.column(column)
        slots = table.slots_of(oids)
        cells = [
            blob if blob is not None else (text or None)
            for blob, text in zip(data.blobs(slots), data.displays(slots), strict=True)
        ]
        pool = self._pool if len(oids) >= POOL_MIN_ROWS else (lambda: None)

        def work(job: Job) -> np.ndarray:
            return search_reactions(
                query, cells, executor=pool(), progress=job.report, cancel=job.cancel_token
            )

        return self._run(
            f"Reaction search: {query.describe()}",
            work,
            lambda hits: on_done(oids, hits),
            lane=QUICK,
            guard=True,
        )


class QueryInput(QWidget):
    """A query structure: typed SMILES or SMARTS, drawn, or taken from the current row."""

    changed = Signal()

    def __init__(
        self,
        current_structure: Callable[[], str] | None = None,
        parent: QWidget | None = None,
        *,
        compact: bool = False,
    ) -> None:
        super().__init__(parent)
        self._current_structure = current_structure
        self.edit = QLineEdit()
        self.edit.setPlaceholderText("SMILES or SMARTS")
        self.edit.setClearButtonEnabled(True)
        if compact:  # a pencil glyph, leaving the text box room in a narrow field
            self.draw_button = QToolButton()
            self.draw_button.setAutoRaise(True)
            self.draw_button.setFixedSize(22, 22)
            self.draw_button.setToolTip("Draw the query in the sketcher")
            self._compact = True
            self._paint_pencil()
        else:
            self.draw_button = QPushButton("Draw…")
            self._compact = False
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4 if compact else 6)
        row.addWidget(self.edit, 1)
        row.addWidget(self.draw_button)
        if current_structure is not None:
            self.current_button = QPushButton("Current Row")
            self.current_button.setToolTip("Use the structure in the table's current row")
            self.current_button.clicked.connect(self.use_current_row)
            row.addWidget(self.current_button)
        self.edit.textChanged.connect(lambda _t: self.changed.emit())
        self.draw_button.clicked.connect(self.draw)

    def _paint_pencil(self) -> None:
        from ..shell.glyphs import glyph

        self.draw_button.setIcon(glyph("pencil", self.palette().text().color()))

    def changeEvent(self, event) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and getattr(self, "_compact", False):
            self._paint_pencil()

    def text(self) -> str:
        return self.edit.text()

    def set_text(self, text: str) -> None:
        self.edit.setText(text)

    #: the molfile last drawn, and the query text it became (to redraw a SMARTS query)
    _drawn: tuple[str, str] | None = None

    def use_current_row(self) -> None:
        text = self._current_structure() if self._current_structure else ""
        if text:
            self.edit.setText(text)

    def draw(self) -> None:
        """Draw the query in the sketcher. A drawing with query features comes back as SMARTS
        (``drawn_query_text``); drawing again reopens the drawing, not the SMARTS."""
        from chemsuite.app.sketcher import SketcherDialog
        from chemsuite.core.structure_input import read_structure

        start = self.edit.text()
        if self._drawn is not None and self._drawn[1] == start:
            start = self._drawn[0]  # SMARTS cannot be drawn back; the molfile can
        dialog = SketcherDialog(start, "Query Structure", self)
        if not dialog.exec():
            return
        result = dialog.molblock
        parsed = read_structure(result) if result else None
        if parsed is None or parsed.mol is None:
            return
        text = drawn_query_text(parsed.mol) if parsed.kind == "molfile" else result
        self._drawn = (result, text)
        self.edit.setText(text)


def describe_query(query: StructureQuery) -> tuple[str, bool]:
    """``(a line about the query, whether it can be searched)``."""
    if not query.text.strip():
        return "Enter a structure, draw one, or use the current row.", False
    try:
        parsed = parse_query(query)
    except QueryError as err:
        return str(err), False
    if parsed.is_smarts:
        return f"SMARTS pattern, {parsed.mol.GetNumAtoms()} atom(s)", True
    from rdkit.Chem import rdMolDescriptors

    return f"Structure {rdMolDescriptors.CalcMolFormula(parsed.mol)}", True


ACTIONS = {"show": "Filter", "select": "Select", "list": "Add to List"}


class StructureSearchDialog(ToolDialog):
    """Asks for a query and what to do with the matches; the search runs after OK."""

    def __init__(self, window: MainWindow, last: StructureQuery | None = None) -> None:
        super().__init__(window, "Structure Search", guide="structure-search-ctrlshiftm")
        self.window_ = window
        workspace = window.workspace
        self.table = workspace.table

        inputs = self.section("Input")
        self.column = QComboBox()
        self.column.addItems(structure_columns(self.table))
        current = window.current_column_name()
        if current in structure_columns(self.table):
            self.column.setCurrentText(current)
        inputs.row("Structures", self.column)
        self.query = QueryInput(self._current_structure, compact=True)
        self.query.setFixedWidth(FIELD_WIDTH)
        inputs.row(
            "Query",
            self.query,
            help="SMILES or SMARTS: type it, draw it, or use the structure in the current row.",
        )
        self.rows = inputs.row(
            "Rows",
            RowsChoice(workspace),
            help="A filter applies to every row (with the other filters), so it has no rows "
            "choice.",
        )
        self._inputs = inputs

        method = self.section("Method")
        self.mode = QComboBox()
        for key, title in MODES.items():
            self.mode.addItem(title, key)
        method.row("Search for", self.mode)
        self.stereo = QCheckBox("Match stereochemistry")
        method.row("Stereo", self.stereo)
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.05, 1.0)
        self.threshold.setSingleStep(0.05)
        self.threshold.setValue(0.7)
        method.row("At least", self.threshold, help="The lowest similarity that counts as a match.")
        self.fingerprint = QComboBox()
        for kind in KINDS.values():
            if kind.listed:
                self.fingerprint.addItem(kind.title, kind.key)
        self.fingerprint.setCurrentIndex(self.fingerprint.findData(DEFAULT_KIND))
        method.row("Fingerprint", self.fingerprint)
        self._method = method

        output = self.section("Output")
        self.then = Segmented(tuple(ACTIONS.items()))
        self.actions = self.then.buttons
        output.row(
            "Matches",
            self.then,
            help="Filter: show only the matches (adds a structure filter). Select: select "
            "them. Add to List: add them to a named list.",
        )
        self.list_name = QComboBox()
        self.list_name.setEditable(True)
        self.list_name.addItems(workspace.lists.names)
        self.list_name.setCurrentText(workspace.lists.names[0] if len(workspace.lists) else "")
        self.list_name.lineEdit().setPlaceholderText("List name")
        output.row("List", self.list_name)
        self.highlight = QCheckBox("Matched atoms in the table")
        self.highlight.setChecked(True)
        output.row("Highlight", self.highlight)
        self._output = output

        if last is not None:
            self.set_query(last)
        elif window.current_structure_text():
            self.query.use_current_row()
        self.query.changed.connect(self._update)
        self.mode.currentIndexChanged.connect(lambda _i: self._update())
        self.then.changed.connect(lambda _k: self._update())
        self.rows.changed.connect(lambda _k: self._update())
        self.column.currentIndexChanged.connect(lambda _i: self._update())
        self.list_name.editTextChanged.connect(lambda _t: self._update())
        self._update()

    def _current_structure(self) -> str:
        return self.window_.current_structure_text(self.column.currentText())

    # -- state --------------------------------------------------------------------------------
    def current_query(self) -> StructureQuery:
        return StructureQuery(
            self.query.text(),
            self.mode.currentData(),
            stereo=self.stereo.isChecked(),
            threshold=self.threshold.value(),
            fingerprint=self.fingerprint.currentData(),
        )

    def set_query(self, query: StructureQuery) -> None:
        self.query.set_text(query.text)
        self.mode.setCurrentIndex(self.mode.findData(query.mode))
        self.stereo.setChecked(query.stereo)
        self.threshold.setValue(query.threshold)
        index = self.fingerprint.findData(query.fingerprint)
        if index >= 0:
            self.fingerprint.setCurrentIndex(index)

    def action(self) -> str:
        return self.then.value()

    def _update(self) -> None:
        mode = self.mode.currentData()
        similar = mode == "similarity"
        self._method.set_row_visible(self.stereo, not similar)
        self._method.set_row_visible(self.threshold, similar)
        self._method.set_row_visible(self.fingerprint, similar)
        # A filter applies to every row (ANDed with the other filters), so no row choice.
        self._inputs.set_row_visible(self.rows, self.action() != "show")
        self._output.set_row_visible(self.list_name, self.action() == "list")
        self._output.set_row_visible(self.highlight, mode == "substructure")
        info, ok = describe_query(self.current_query())
        self.query.edit.setToolTip(info if ok else "")
        self.refresh()

    def problem(self) -> str:
        if not self.column.count():
            return "The table has no structure column."
        info, ok = describe_query(self.current_query())
        if not ok:
            return info
        if self.action() == "list" and not self.list_name.currentText().strip():
            return "Name the list to add the matches to."
        return ""

    def verb(self) -> str:
        if self.action() == "show":
            return "Search"
        return f"Search {len(self.rows.oids()):,} Rows"

    # -- running ------------------------------------------------------------------------------
    def accept(self) -> None:
        if self._problem():
            self.refresh()
            return
        if self.start():
            QDialog.accept(self)

    def start(self) -> bool:
        """Start the search; False (dialog stays open) when the query cannot be searched."""
        query = self.current_query()
        column = self.column.currentText()
        window = self.window_
        try:
            parse_query(query)
        except QueryError as err:
            self.report(str(err))
            return False
        window.remember_structure_query(query)
        highlight = self.highlight.isChecked() and query.mode == "substructure"
        action = self.action()
        if action == "show":
            card = window.filter_column(column, "structure")
            card.set_query(query, highlight=highlight)
            self.job = card.job
            return True
        workspace = window.workspace
        oids = None if self.rows.value() == "all" else self.rows.oids()
        list_name = self.list_name.currentText().strip()
        searcher = window.structure_searcher

        def done(searched: np.ndarray, result: SearchResult) -> None:
            hits = searched[result.hits]
            if highlight:
                searcher.highlight(column, highlight_for(query))
            if action == "select":
                workspace.selection.set(hits.tolist(), "structure search")
                message = f"Selected {len(hits):,} match(es) for {query.describe()}"
            else:
                try:
                    workspace.lists.add(list_name, hits)
                except ListError as err:
                    window.show_message(f"Could not add to “{list_name}”: {err}")
                    return
                message = f"Added {len(hits):,} match(es) for {query.describe()} to “{list_name}”"
            window.show_message(message, 8000)
            if self.job is not None:
                self.job.set_summary(message)

        self.job = searcher.start(column, query, done, oids=oids)
        return True
