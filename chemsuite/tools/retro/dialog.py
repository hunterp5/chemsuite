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

"""The Retrosynthesis dialog: pick rows and limits, run the search, write the columns."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.dialogs import RowsChoice, ToolDialog
from chemsuite.app.settings_hub import settings_hub
from chemsuite.app.tools import ToolContext
from chemsuite.retro import (
    OUTPUT_COLUMNS,
    SearchOptions,
    catalogue,
    ready_problem,
    resolve_config,
    resolve_python,
)
from chemsuite.retro.runner import SearchReply, search, unique_targets

from .tool import ROUTES_TOOL_ID

#: the options last used, so a second run starts where the first left off (this session only)
_last: SearchOptions | None = None

ALGORITHM_HELP = (
    "Monte Carlo tree search is AiZynthFinder's default. Retro*, breadth-first and proof "
    "number search explore the same reaction templates in a different order; a hard target "
    "sometimes falls to one of them when Monte Carlo runs out of time."
)
EXPLORATION_HELP = (
    "How widely Monte Carlo tries less likely reactions (the C constant). Higher explores "
    "more broadly, lower follows the policy's favourite first."
)


class RetrosynthesisDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Retrosynthesis", guide="retrosynthesis")
        self.ctx = ctx
        settings = settings_hub().settings
        self.python = resolve_python(settings)
        self.config = resolve_config(settings)
        self.catalogue = catalogue(self.config)
        saved = _last or SearchOptions()

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        # a search costs minutes per structure, so it starts on the selection when there is
        # one rather than on the whole table
        start = "selected" if len(ctx.selection) else "all"
        self.rows = inputs.row("Rows", RowsChoice(ctx, start))

        search_section = self.section("Search")
        self.policies = _checklist(self.catalogue.policies, saved.policies)
        search_section.row(
            "Policies",
            self.policies,
            help="The expansion policies in your config.yml. Every one is used unless you "
            "clear some. Ringbreaker is trained on ring-forming reactions.",
        )
        self.stocks = _checklist(self.catalogue.stocks, saved.stocks)
        search_section.row(
            "Stocks",
            self.stocks,
            help="Which stock files count as buyable. A route is solved when every precursor "
            "is in one of them.",
        )
        self.time_limit = _spin(5, 3600, saved.time_limit, suffix=" s")
        search_section.row(
            "Time limit",
            self.time_limit,
            help="Each target gets this long. A molecule often takes a minute or more.",
        )
        self.iterations = _spin(1, 10000, saved.iteration_limit)
        search_section.row(
            "Iterations",
            self.iterations,
            help="The search also stops after this many iterations, whichever comes first.",
        )
        self.max_routes = _spin(1, 500, saved.max_routes)
        search_section.row(
            "Routes to keep",
            self.max_routes,
            help="How many of the best routes to keep for the Routes pane.",
        )

        more = self.section("More Options", collapsed=True)
        self.algorithm = QComboBox()
        for label, value in _algorithms():
            self.algorithm.addItem(label, value)
        index = self.algorithm.findData(saved.algorithm)
        self.algorithm.setCurrentIndex(max(0, index))
        more.row("Algorithm", self.algorithm, help=ALGORITHM_HELP)
        self.exploration = QDoubleSpinBox()
        self.exploration.setRange(0.0, 10.0)
        self.exploration.setDecimals(2)
        self.exploration.setSingleStep(0.1)
        self.exploration.setValue(saved.exploration)
        more.row("Exploration (C)", self.exploration, help=EXPLORATION_HELP)
        self.filters = _checklist(self.catalogue.filters, saved.filters)
        more.row(
            "Filters",
            self.filters,
            help="Filter policies reject reactions the model believes would not work. Every "
            "one in the config is used unless you clear some.",
        )
        self.max_transforms = _spin(1, 20, saved.max_transforms)
        more.row(
            "Longest route",
            self.max_transforms,
            help="The most reactions a route may have (AiZynthFinder's max_transforms).",
        )
        self.return_first = QCheckBox("Stop at the first solved route")
        self.return_first.setChecked(saved.return_first)
        more.row("", self.return_first)
        self.exclude_target = QCheckBox("Ignore the target's own stock entry")
        self.exclude_target.setChecked(saved.exclude_target_from_stock)
        more.row(
            "",
            self.exclude_target,
            help="On, a target that is itself buyable is still taken apart. Off, it counts as "
            "solved on the spot.",
        )
        self.cluster = QCheckBox("Cluster similar routes")
        self.cluster.setChecked(saved.cluster)
        more.row("", self.cluster, help="Group routes that share most of their reactions.")

        output = self.section("Output")
        self.skip_done = QCheckBox("Skip rows already searched")
        self.skip_done.setChecked(True)
        output.row(
            "",
            self.skip_done,
            help="Rows that already have AiZynth values are left alone. Clear those cells, or "
            "turn this off, to search them again.",
        )

        self.structure.currentTextChanged.connect(lambda *_: self.refresh())
        self.rows.changed.connect(lambda *_: self.refresh())
        for box in (self.skip_done, self.return_first):
            box.toggled.connect(lambda *_: self.refresh())
        self.refresh()

    # -- what will happen ---------------------------------------------------------------------
    def options(self) -> SearchOptions:
        return SearchOptions(
            policies=_chosen(self.policies),
            filters=_chosen(self.filters),
            stocks=_chosen(self.stocks),
            algorithm=str(self.algorithm.currentData() or "mcts"),
            exploration=self.exploration.value(),
            max_transforms=self.max_transforms.value(),
            iteration_limit=self.iterations.value(),
            time_limit=self.time_limit.value(),
            return_first=self.return_first.isChecked(),
            exclude_target_from_stock=self.exclude_target.isChecked(),
            max_routes=self.max_routes.value(),
            cluster=self.cluster.isChecked(),
        )

    def oids(self, *, skip: bool | None = None) -> np.ndarray:
        """The rows to search: the chosen rows, less the ones already searched."""
        oids = self.rows.oids()
        if not (self.skip_done.isChecked() if skip is None else skip):
            return oids
        table = self.ctx.table
        done = next((name for name in OUTPUT_COLUMNS if name in table), "")
        if not done:
            return oids
        slots = table.slots_of(oids)
        searched = table.column(done).displays(slots)
        keep = [bool(str(text or "").strip()) is False for text in searched]
        return oids[np.asarray(keep, dtype=bool)]

    def targets(self, *, skip: bool | None = None) -> list[str]:
        """The distinct structures to search, in row order (blank cells are not targets)."""
        table = self.ctx.table
        name = self.structure.currentText()
        if name not in table:
            return []
        slots = table.slots_of(self.oids(skip=skip))
        return unique_targets(table.column(name).displays(slots))

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to search."
        ready = ready_problem(settings_hub().settings)
        if ready:
            return ready
        if self.catalogue.problem:
            return self.catalogue.problem
        if not len(self.rows.oids()):
            return "There are no rows to search."
        if not self.targets():
            if self.skip_done.isChecked() and self.targets(skip=False):
                return "Every structure in the chosen rows has been searched already."
            return "The chosen rows have no structures."
        return ""

    def verb(self) -> str:
        targets = len(self.targets())
        if targets == 1:
            return "Search 1 Structure"
        return f"Search {targets:,} Structures"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        global _last
        ctx = self.ctx
        options = self.options()
        _last = options
        targets = self.targets()
        oids = self.oids()
        slots = ctx.table.slots_of(oids)
        of_row = [
            str(text or "").strip()
            for text in ctx.table.column(self.structure.currentText()).displays(slots)
        ]
        python, config = self.python, self.config
        resume = self.skip_done.isChecked()

        def work(job):
            return search(
                python,
                str(config),
                targets,
                options,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(replies: list[SearchReply]) -> None:
            finish(ctx, oids, of_row, replies, resume=resume)

        label = f"Planning routes to {len(targets):,} structure(s)…"
        return ctx.run(label, work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    of_row: list[str],
    replies: list[SearchReply],
    *,
    resume: bool,
) -> dict[str, str]:
    """Write the ``AiZynth …`` columns for the rows the search covered, and show the routes."""
    by_target = {reply.target: reply for reply in replies}
    names = _column_order(replies)
    rows = [
        (oid, by_target[target])
        for oid, target in zip(oids.tolist(), of_row, strict=True)
        if target in by_target
    ]
    written: dict[str, str] = {}
    if rows and names:
        covered = np.asarray([oid for oid, _reply in rows], dtype=np.int64)
        columns = {
            name: (covered, [reply.columns.get(name, "") for _oid, reply in rows]) for name in names
        }
        written = ctx.apply_results(columns, label="Retrosynthesis", overwrite=resume)
    show_routes(ctx, replies)
    ctx.show_message(_summary(replies, len(rows)), 8000)
    return written


def show_routes(ctx: ToolContext, replies: list[SearchReply]) -> None:
    """Add the routes to the Routes pane, opening it if it is not already there."""
    targets = [
        {"smiles": reply.target, "routes": reply.routes} for reply in replies if reply.routes
    ]
    if not targets:
        return
    panel = ctx.find_panel(ROUTES_TOOL_ID)
    if panel is None:
        panel = ctx.open_tool(ROUTES_TOOL_ID)
    if panel is not None:
        panel.add_targets(targets)


def _summary(replies: list[SearchReply], rows: int) -> str:
    solved = sum(1 for reply in replies if reply.columns.get("AiZynth Solved") == "Yes")
    failed = [reply for reply in replies if not reply.ok]
    text = (
        f"Searched {len(replies):,} structure(s) in {rows:,} row(s); {solved:,} with a solved route"
    )
    if failed:
        text += f"; {len(failed):,} could not be searched ({failed[0].error.splitlines()[0]})"
    return text


def _column_order(replies: list[SearchReply]) -> list[str]:
    """The named columns first, then anything else AiZynthFinder reported, in first-seen order."""
    seen: list[str] = []
    for reply in replies:
        for name in reply.columns:
            if name not in seen:
                seen.append(name)
    named = [name for name in OUTPUT_COLUMNS if name in seen]
    return named + [name for name in seen if name not in named]


# -- small widgets ------------------------------------------------------------------------------
def _algorithms() -> tuple[tuple[str, str], ...]:
    from chemsuite.retro.config import ALGORITHMS

    return ALGORITHMS


def _spin(low: int, high: int, value: int, *, suffix: str = "") -> QSpinBox:
    spin = QSpinBox()
    spin.setRange(low, high)
    spin.setValue(max(low, min(high, int(value))))
    if suffix:
        spin.setSuffix(suffix)
    return spin


def _checklist(names: tuple[str, ...], chosen: tuple[str, ...]) -> QWidget:
    """A checkbox per name; nothing checked means "everything in the config"."""
    holder = QWidget()
    layout = QVBoxLayout(holder)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(4)
    boxes: dict[str, QCheckBox] = {}
    for name in names:
        box = QCheckBox(name)
        box.setChecked(name in chosen if chosen else True)
        layout.addWidget(box)
        boxes[name] = box
    holder.boxes = boxes  # type: ignore[attr-defined]
    return holder


def _chosen(holder: QWidget) -> tuple[str, ...]:
    """The ticked names, or ``()`` when every one is ticked (the config's own choice)."""
    boxes: dict[str, QCheckBox] = getattr(holder, "boxes", {})
    picked = tuple(name for name, box in boxes.items() if box.isChecked())
    return () if len(picked) == len(boxes) else picked
