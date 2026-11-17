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
"""Similarity Search: the table's rows, or a library file, scored against one or more
references by fingerprints, pharmacophores, substructure and 3D shape (``chemsuite.similarity``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QLineEdit,
    QPushButton,
    QSpinBox,
)

from chemsuite.analysis.fingerprints import cache_for
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, Segmented, ToolDialog, inline
from chemsuite.core.table.columns import NumberColumn, StructureColumn
from chemsuite.core.wording import plural
from chemsuite.io.formats import UnsupportedFormat, read_filter, split_suffix
from chemsuite.similarity import model as sm
from chemsuite.similarity.fusion import unit
from chemsuite.similarity.run import RANKING
from chemsuite.tools.superpose.poses import (
    DELTA_E,
    new_reference_id,
    output_name,
    put,
    reference_asset,
)

from .methods import Methods
from .references import Entry, References

TABLE, FILE = "table", "file"
#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 20
FILE_FILTER = read_filter(reactions=False) + ";;All files (*)"

_last: dict | None = None  # the settings of the last run, remembered for the session


class SimilarityDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Similarity Search", guide="similarity-search")
        self.ctx = ctx

        inputs = self.inputs = self.section("Input")
        self.where = Segmented(((TABLE, "This Table"), (FILE, "A File")))
        inputs.row(
            "Search in",
            self.where,
            help="This Table: score its rows and write the scores as columns. A File: screen a "
            "library (SDF, SMILES, CSV…) without opening it; the best hits open as a new "
            "document.",
        )
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row(
            "Structures",
            self.structure,
            help="The rows' structures; the current row and selected rows are read from it too.",
        )
        self.rows = inputs.row("Rows", RowsChoice(ctx))
        self.path = QLineEdit()
        self.path.setPlaceholderText("Library file")
        browse = QPushButton("Browse…")
        browse.setAutoDefault(False)
        browse.clicked.connect(self.browse)
        self.file_field = inline(self.path, browse, stretch=0)
        inputs.row("File", self.file_field)

        self.refs = References(ctx, self.section("References"), self.structure.currentText)
        method = self.section("Method")
        output = self.output = self.section("Output")
        more = self.more = self.section("More Options", collapsed=True)
        self.methods = Methods(method, more)

        self.sort = QCheckBox("Most similar first")
        self.sort.setChecked(True)
        output.row(
            "Sort",
            self.sort,
            help="Sort the table by Consensus when there is one, else by the first score.",
        )
        self.select = QCheckBox("At least")
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.0, 1.0)
        self.threshold.setSingleStep(0.05)
        self.threshold.setValue(0.7)
        self.threshold.setEnabled(False)
        self.select.toggled.connect(self.threshold.setEnabled)
        self.select_field = inline(self.select, self.threshold)
        output.row(
            "Select",
            self.select_field,
            help="Select the rows scoring at least this much, by the score the table is sorted "
            "by, on a 0–1 scale (Shape + Features, 0–2, is halved).",
        )
        self.poses = QCheckBox()
        output.row(
            "Keep aligned poses",
            self.poses,
            help="Write each row's best 3D pose, fitted onto the reference it is most like, in "
            "“<column> superposed”; the Browser shows it with that reference.",
        )
        self.top = QSpinBox()
        self.top.setRange(0, 1_000_000)
        self.top.setSpecialValueText("Every hit")
        self.top.setValue(1000)
        output.row(
            "Keep",
            self.top,
            help="The best this many compounds of the file, ranked by the fingerprint as it is "
            "read; every method then scores them.",
        )
        self.minimum = QDoubleSpinBox()
        self.minimum.setRange(0.0, 1.0)
        self.minimum.setSingleStep(0.05)
        output.row(
            "At least",
            self.minimum,
            help="Keep only compounds whose fingerprint similarity reaches this.",
        )

        if _last is not None:
            self._restore(_last)
        self.where.changed.connect(lambda _k: self._update())
        self.structure.currentIndexChanged.connect(lambda _i: self._update())
        self.rows.changed.connect(lambda *_a: self.refresh())
        self.path.textChanged.connect(lambda _t: self.refresh())
        self.refs.changed.connect(self._update)
        self.methods.changed.connect(self._update)
        if not ctx.structure_columns():
            self.where.set_value(FILE)
            self.where.set_enabled(TABLE, False, "The table has no structure column.")
        self._update()

    # -- settings -----------------------------------------------------------------------------
    def _restore(self, last: dict) -> None:
        self.methods.set_search(last["search"])
        self.refs.fuse.setCurrentIndex(max(0, self.refs.fuse.findData(last["search"].fuse)))
        self.sort.setChecked(last["sort"])
        self.poses.setChecked(last["poses"])
        self.top.setValue(last["top"])
        self.minimum.setValue(last["minimum"])

    def in_file(self) -> bool:
        return self.where.value() == FILE

    def _update(self) -> None:
        file, table = self.in_file(), not self.in_file()
        inputs, output = self.inputs, self.output
        inputs.set_row_visible(self.structure, table)
        inputs.set_row_visible(self.rows, table)
        inputs.set_row_visible(self.file_field, file)
        for field in (self.sort, self.select_field):
            output.set_row_visible(field, table)
        output.set_row_visible(self.poses, table and self.methods.has_3d())
        output.set_row_visible(self.top, file)
        output.set_row_visible(self.minimum, file)
        self.refs.set_fuse_visible(self.refs.from_selection() or len(self.refs.entries()) > 1)
        self.refresh()

    def browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Library", self.path.text(), FILE_FILTER)
        if path:
            self.path.setText(path)

    # -- checks -------------------------------------------------------------------------------
    def targets(self) -> np.ndarray:
        """The rows to score (none when screening a file)."""
        return np.empty(0, dtype=np.int64) if self.in_file() else self.rows.oids()

    def search(self) -> sm.Search:
        return self.methods.search(self.refs.fusion())

    def problem(self) -> str:
        if self.in_file():
            problem = self._file_problem()
        elif not self.structure.currentText():
            problem = "There is no structure column."
        else:
            problem = "" if len(self.rows.oids()) else "There are no rows to work on."
        problem = problem or self.refs.problem()
        if problem:
            return problem
        if not self.methods.count():
            return "Choose at least one method."
        try:
            search = self.search()
            sm.check_references(self.refs.references(self.targets()), search)
        except ValueError as err:
            text = str(err)
            return text[:1].upper() + text[1:] + "."
        return ""

    def _file_problem(self) -> str:
        text = self.path.text().strip()
        if not text:
            return "Choose a library file to screen."
        if not Path(text).is_file():
            return "The file was not found."
        try:
            split_suffix(text)
        except UnsupportedFormat:
            return "That file type cannot be read."
        return ""

    def verb(self) -> str:
        if self.in_file():
            return "Screen File"
        return f"Compare {len(self.rows.oids()):,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        global _last
        search = self.search()
        _last = {
            "search": search,
            "sort": self.sort.isChecked(),
            "poses": self.poses.isChecked(),
            "top": self.top.value(),
            "minimum": self.minimum.value(),
        }
        return self.screen(search) if self.in_file() else self.calculate(search)

    def calculate(self, search: sm.Search):
        ctx = self.ctx
        column = self.structure.currentText()
        targets = self.targets()
        entries = self.refs.entries()
        references = self.refs.references(targets, entries)
        cells = ctx.structure_inputs(column, targets)
        cache = cache_for(ctx.table)
        generation = cache.generation(column)
        kinds = {m.fingerprint for m in search.methods if m.kind == sm.FINGERPRINT}
        kinds.add(RANKING.fingerprint)
        known = {k: fps for k in kinds if (fps := cache.get(column, k, targets)) is not None}
        keep_poses = self.poses.isChecked() and self.methods.has_3d()
        pose_ids = [new_reference_id(ctx) for _ in references] if keep_poses else None

        def work(job):
            from chemsuite.similarity.run import score_cells

            return score_cells(
                cells,
                references,
                search,
                fingerprints=known,
                pose_ids=pose_ids,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(scores: sm.Scores) -> None:
            for kind, fps in scores.fingerprints.items():
                if kind not in known and len(fps) == len(targets):
                    cache.put(column, kind, targets, fps, generation)
            self.finish(scores, search, targets[: len(scores)], entries, column, pose_ids)

        return ctx.run("Searching for similar compounds…", work, done, keep_partial=True)

    def finish(
        self,
        scores: sm.Scores,
        search: sm.Search,
        oids: np.ndarray,
        entries: list[Entry],
        column: str,
        pose_ids: list[str] | None,
    ) -> None:
        ctx = self.ctx
        results: dict[str, tuple] = {}
        for method in search.methods:
            values = scores.columns.get(method.column)
            if values is not None:
                results[method.column] = (oids, NumberColumn.from_floats(values, 3))
        if scores.consensus is not None:
            results[sm.CONSENSUS_COLUMN] = (oids, NumberColumn.from_floats(scores.consensus, 3))
        if len(entries) > 1:
            names = [entries[i].name if i >= 0 else "" for i in scores.nearest.tolist()]
            results[sm.NEAREST_COLUMN] = (oids, names)
        extra = []
        if pose_ids is not None and scores.poses is not None:
            poses = scores.poses
            blobs = [b"" if p is None else p[0] for p in poses]
            displays = ["" if p is None else p[1] for p in poses]
            results[output_name(column)] = (oids, StructureColumn.from_parsed(blobs, displays))
            delta_e = np.array([np.nan if p is None else p[3] for p in poses], dtype=np.float64)
            if np.isfinite(delta_e).any():
                results[DELTA_E] = (oids, NumberColumn.from_floats(delta_e, 2))
            method = next(m for m in search.methods if m.kind in sm.THREE_D)
            for index in sorted({p[2] for p in poses if p is not None}):
                entry = entries[index]
                meta = {**entry.meta, "method": method.kind}
                asset = reference_asset(pose_ids[index], entry.name, entry.mol, meta)
                extra.append(put(asset))
        if not results:
            ctx.show_message("Nothing was scored")
            return
        written = ctx.apply_results(results, label="Similarity search", extra_ops=extra)
        ranked, values = self._ranking(scores, search)
        message = f"Scored {plural(int(np.isfinite(values).sum()), 'row')}"
        if ranked in written:
            message += f" in “{written[ranked]}”"
            if self.sort.isChecked():
                ctx.table.set_sort([(written[ranked], False)])
        if self.select.isChecked():
            hits = oids[np.nan_to_num(values, nan=-1.0) >= self.threshold.value()]
            ctx.selection.set(hits.tolist(), "similarity")
            message += f"; selected {len(hits):,} at ≥ {self.threshold.value():.2f}"
        ctx.show_message(message)

    @staticmethod
    def _ranking(scores: sm.Scores, search: sm.Search) -> tuple[str, np.ndarray]:
        """The column the table is sorted and selected by, and its values on a 0–1 scale."""
        if scores.consensus is not None:
            return sm.CONSENSUS_COLUMN, scores.consensus
        method = next(m for m in search.methods if m.column in scores.columns)
        return method.column, unit(method, scores.columns[method.column])

    def screen(self, search: sm.Search):
        ctx = self.ctx
        path = Path(self.path.text().strip())
        references = self.refs.references(np.empty(0, dtype=np.int64))  # none is a file row
        top, minimum = self.top.value(), self.minimum.value()

        def work(job):
            from chemsuite.similarity.screen import screen_file

            return screen_file(
                path,
                references,
                search,
                top_k=top,
                threshold=minimum,
                executor=ctx.process_pool(),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(hits) -> None:
            read = f"{hits.scanned:,} compounds of {path.name}"
            if hits.cancelled:
                read = f"the first {read} (stopped)"
            if not len(hits):
                ctx.show_message(f"No hits in {read}")
                return
            ctx.open_table(hits.table, f"{plural(len(hits), 'hit')} from {read}")

        return ctx.run(f"Screening {path.name}…", work, done, keep_partial=True)
