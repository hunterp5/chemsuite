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
"""Data ▸ Dimensionality Reduction ▸ PCA…, UMAP…, t-SNE… or PaCMAP…: the map's coordinates
as columns, plotted.

The dialog is the same for every method: Input (features, rows, and the rows the map is
fitted on), Method and More Options (each method's own, from :mod:`.map_options`, plus how
the features are made), and Output (what the plot is colored and sized by). It opens with the
settings it last ran with (``ToolDialog.remember``).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox

from chemsuite.analysis import embedding as em
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import FitChoice, RowsChoice, Segmented, ToolDialog
from chemsuite.core.table.columns import NumberColumn
from chemsuite.core.wording import plural
from chemsuite.plotting.data import colorable_columns, number_columns
from chemsuite.tools.sar.common import activity_combo
from chemsuite.tools.similarity.common import FeatureChooser

from .common import features_task
from .map_options import NONE, OPTIONS, MethodOptions, PcaOptions

#: methods that hold Python's GIL for seconds (numba): run in a worker process, or the window
#: freezes while they work (PCA is NumPy, which lets go of it, and reports its progress)
IN_PROCESS = frozenset(["umap", "tsne", "pacmap"])


class SpaceMapDialog(ToolDialog):
    """One method's map (Data ▸ Dimensionality Reduction has an entry for each)."""

    def __init__(self, ctx: ToolContext, method: str = "umap") -> None:
        super().__init__(ctx.parent, em.METHODS[method], guide="dimensionality-reduction")
        self.ctx = ctx
        self.method = method
        self.plot = None
        self.loadings_panel = None

        inputs = self.section("Input")
        section = self.section("Method")
        more = self.section("More Options", collapsed=False)
        output = self.section("Output")
        self._inputs, self._method_section = inputs, section
        self._more_section, self._output = more, output
        self.features = FeatureChooser(ctx, metric=False, section=inputs, both=True, more=more)
        self.rows = inputs.row("Rows", RowsChoice(ctx))
        self.fit = inputs.row("Fit on", FitChoice(ctx),
                              help="The rows the map is made from; the other rows are then "
                              "placed onto it. Fit on a reference set to see where new "
                              "compounds fall in its space.")  # fmt: skip
        self.dims = Segmented((("2", "2D"), ("3", "3D")), "2")
        section.row(
            "Dimensions",
            self.dims,
            help="3D writes three columns and opens a 3D Scatter you can turn.",
        )
        self.settings: MethodOptions = OPTIONS[method](ctx, section, more)

        self.color = QComboBox()
        self.size = QComboBox()
        self.color.addItem(NONE, None)
        self.size.addItem(NONE, None)
        for name in colorable_columns(ctx.table):
            self.color.addItem(name, name)
        for name in number_columns(ctx.table):
            self.size.addItem(name, name)
        activity = activity_combo(ctx).currentText()
        self.color.setCurrentIndex(max(0, self.color.findData(activity or None)))
        output.row("Color by", self.color, help="The column the plot colors its points by.")
        output.row("Size by", self.size, help="The number column the plot sizes points by.")
        self.open_plot = QCheckBox("Open a plot")
        self.open_plot.setChecked(True)
        output.row("Plot", self.open_plot, help="Off: only write the map's columns to the table.")

        for signal in (
            self.rows.changed,
            self.features.changed,
            self.dims.changed,
            self.fit.currentIndexChanged,
            self.open_plot.toggled,
            *self.settings.signals(),
        ):
            signal.connect(lambda *_: self._update())
        self._update()
        self.remember(ctx, f"space_map_{method}", self._fields())

    def _fields(self) -> dict:
        fields = {"dims": self.dims, "fit": self.fit, "open_plot": self.open_plot}
        fields |= self.features.remembered()
        return fields | self.settings.fields()

    def _count(self) -> int:
        return len(self.rows.oids())

    def _fit_rows(self, oids: np.ndarray) -> np.ndarray | None:
        return self.fit.rows(oids) if self.settings.fits else None

    def problem(self) -> str:
        problem = self.features.problem()
        if problem:
            return problem
        oids = self.rows.oids()
        if len(oids) < 3:
            return "A map needs at least 3 rows."
        fit = self._fit_rows(oids)
        if fit is not None and len(fit) < 5:
            return "Fit on at least 5 of the rows mapped."
        return self.settings.problem(
            rows=len(oids), fitting=fit is not None, bits=self.features.bits()
        )

    def verb(self) -> str:
        return f"Map {self._count():,} Rows"

    def _update(self) -> None:
        three = self.dims.value() == "3"
        columns_only = self.features.kind.currentData() is None
        self.settings.update(bits=self.features.bits(), three=three, columns_only=columns_only)
        self._inputs.set_row_visible(self.fit, self.settings.fits)
        for widget in (self.color, self.size):
            self._output.set_row_visible(widget, self.open_plot.isChecked())
        self.refresh()

    def run(self):
        return self.make_map()

    def options(self) -> dict:
        """What ``embed`` is called with (bar the rows to fit on)."""
        return {"dims": int(self.dims.value()), **self.settings.options(self.rows.oids())}

    def make_map(self):
        ctx = self.ctx
        oids = self.rows.oids()
        method = self.method
        options = {"dims": int(self.dims.value()), **self.settings.options(oids)}
        fit = self._fit_rows(oids)
        if fit is not None:
            options["fit"] = fit
        build, remember = features_task(ctx, self.features, oids)
        plot = self.open_plot.isChecked()
        color, size = self.color.currentData(), self.size.currentData()
        columns = self.features.number_columns()
        settings = self.settings
        loadings = isinstance(settings, PcaOptions) and settings.loadings.isChecked()
        loadings = loadings and self.features.kind.currentData() is None

        def work(job):
            features = build(job)
            if method in IN_PROCESS:
                return ctx.in_process(job, em.embed, features, method, **options)
            return em.embed(
                features, method, progress=job.report, cancel=job.cancel_token, **options
            )

        def done(result: em.Embedding) -> None:
            remember()
            self.finish(oids, result, color, size, plot=plot)
            if loadings and result.loadings is not None:
                self.show_loadings(result, columns)

        return ctx.run(f"Calculating {em.METHODS[method]}…", work, done)

    def finish(
        self,
        oids: np.ndarray,
        result: em.Embedding,
        color: str | None = None,
        size: str | None = None,
        *,
        plot: bool = True,
    ):
        ctx = self.ctx
        dims = result.coords.shape[1]
        title = em.METHODS[result.method]
        if result.components is not None:  # PC1…PCk; the plot takes the chosen ones
            k = result.components.shape[1]
            columns = fresh_names(ctx.table.column_names, [f"PC{i + 1}" for i in range(k)])
            values = result.components
            settings = self.settings
            pair = settings.pair() if isinstance(settings, PcaOptions) and dims == 2 else (1, 2, 3)
            axes = [columns[c - 1] for c in pair[:dims]]
        else:
            columns = map_columns(ctx.table.column_names, title, dims)
            values = result.coords
            axes = columns
        names = ctx.apply_results(
            {
                name: (oids, NumberColumn.from_floats(values[:, i], 3))
                for i, name in enumerate(columns)
            },
            label=f"{title} map",
        )
        axes = [names[c] for c in axes]
        message = f"Placed {result.placed:,} of {plural(len(oids), 'row')} on the map"
        if result.notes:
            message += "; " + "; ".join(result.notes)
        ctx.show_message(message)
        if not plot:
            return None
        spec = {"color": color if color is not None and color in ctx.table else None}
        if size is not None and size in ctx.table:
            spec["size"] = size
        if dims == 3:
            self.plot = ctx.open_tool("scatter3d")
            self.plot.set_spec(x=axes[0], y=axes[1], z=axes[2], **spec)
        else:
            self.plot = ctx.open_tool("scatter")
            self.plot.set_spec(x=axes[0], y=axes[1], **spec)
        return self.plot

    def show_loadings(self, result: em.Embedding, columns: list[str]):
        from .loadings import open_loadings

        kept = [name for name, keep in zip(columns, result.kept, strict=True) if keep]
        self.loadings_panel = open_loadings(self.ctx, kept, result.loadings)
        return self.loadings_panel


def fresh_names(existing, names: list[str]) -> list[str]:
    """*names* as new columns, never an existing one: "PC1", "PC2"; then "PC1 (1)", … — one
    suffix for the whole set, so they stay recognizably together."""
    taken = set(existing)
    run = 0
    while True:
        suffix = f" ({run})" if run else ""
        out = [f"{name}{suffix}" for name in names]
        if not taken.intersection(out):
            return out
        run += 1


def map_columns(existing, method: str, dims: int) -> list[str]:
    """New column names for one map, never an existing column: "UMAP 1", "UMAP 2"; the next
    map "UMAP 1 (1)", "UMAP 2 (1)"; then "(2)" and so on. A map's columns share one number,
    so they stay recognizably together (rename them in the table afterwards)."""
    return fresh_names(existing, [f"{method} {i + 1}" for i in range(dims)])
