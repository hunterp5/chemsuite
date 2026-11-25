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

"""The Permeability dialog."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QCheckBox, QComboBox, QVBoxLayout, QWidget

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.wording import plural
from chemsuite.permeability import COMMON, ENDPOINTS, Endpoint
from chemsuite.permeability import models as shipped
from chemsuite.permeability.predict import Results, predict_permeability

# The columns each endpoint can write, as tools/qsar/columns.py names them. The prediction
# itself is not optional: a run that writes nothing has nothing to say.
from chemsuite.tools.qsar.columns import names_for, values_for
from chemsuite.tools.qsar.work import PARALLEL_FROM

#: output key -> (label in the dialog, tooltip)
OUTPUTS = {
    "interval": (
        "Interval",
        "The range the measured value should fall in, at the model's confidence. Wide where "
        "the model has little to go on, narrow where it has company.",
    ),
    "domain": (
        "Applicability",
        "How far the compound sits from the model's training set, and whether that is further "
        "than the training compounds sit from each other. A prediction outside the domain is "
        "one the interval was never calibrated for.",
    ),
}
DEFAULT_OUTPUTS = ("interval",)


class PermeabilityDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Permeability", guide="permeability")
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        self.endpoints: dict[str, QCheckBox] = {}
        method = self.section("Method")
        method.row(
            "Endpoints",
            self._endpoint_boxes(),
            help="Each endpoint is a model of its own, fitted on its own public assay data, "
            "and writes its own columns.",
        )

        output = self.section("Output")
        boxes = QWidget()
        layout = QVBoxLayout(boxes)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        self.outputs: dict[str, QCheckBox] = {}
        for key, (label, tip) in OUTPUTS.items():
            box = QCheckBox(label)
            box.setToolTip(tip)
            box.setChecked(key in DEFAULT_OUTPUTS)
            box.toggled.connect(lambda *_: self.refresh())
            layout.addWidget(box)
            self.outputs[key] = box
        output.row(
            "Also write",
            boxes,
            help="The prediction is always written. When a name is already taken, every "
            "column of the run gets the same suffix: (2), (3)…",
        )

        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def _endpoint_boxes(self) -> QWidget:
        """One checkbox per endpoint, each with what its assay measures in its tooltip."""
        holder = QWidget()
        layout = QVBoxLayout(holder)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        for endpoint in ENDPOINTS:
            box = QCheckBox(endpoint.title)
            box.setToolTip(endpoint.tip)
            box.setChecked(endpoint.id in COMMON and shipped.installed(endpoint))
            box.setEnabled(shipped.installed(endpoint))
            box.toggled.connect(lambda *_: self.refresh())
            layout.addWidget(box)
            self.endpoints[endpoint.id] = box
        return holder

    # -- what will happen ---------------------------------------------------------------------
    def chosen_endpoints(self) -> list[Endpoint]:
        """In the package's own order, not the order of the checkboxes."""
        return [e for e in ENDPOINTS if self.endpoints[e.id].isChecked()]

    def chosen(self) -> list[str]:
        """The optional column groups, plus the prediction, as columns.py names them."""
        return ["prediction", *(key for key, box in self.outputs.items() if box.isChecked())]

    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to predict from."
        if not self._rows():
            return "There are no rows to predict."
        if not any(shipped.installed(e) for e in ENDPOINTS):
            return "No permeability model is installed."
        if not self.chosen_endpoints():
            return "Choose at least one endpoint."
        return ""

    def verb(self) -> str:
        return f"Predict for {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        ctx = self.ctx
        oids = self.rows.oids()
        structures = self.structure.currentText()
        cells = ctx.flat_structure_inputs(structures, oids)
        endpoints = self.chosen_endpoints()
        chosen = self.chosen()
        # the descriptor pool is the only heavy part, and it is worth starting the shared
        # workers only once there are enough rows to pay for them
        executor = ctx.workers_for(len(cells), PARALLEL_FROM)

        def work(job):
            return predict_permeability(
                cells,
                endpoints,
                executor=executor,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(results: Results) -> None:
            finish(ctx, oids[: len(results)], results, endpoints, chosen)

        return ctx.run("Predicting permeability…", work, done, keep_partial=True)


def finish(
    ctx: ToolContext,
    oids: np.ndarray,
    results: Results,
    endpoints: list[Endpoint],
    chosen: list[str],
) -> dict[str, str]:
    """Write the chosen columns for *oids* (the rows *results* covers), always as new columns."""
    requested: list[str] = []
    flat: list[object] = []
    for endpoint in endpoints:
        model = shipped.model(endpoint)
        said = results.per_endpoint[endpoint.id].as_prediction()
        requested += names_for(endpoint.column, chosen, endpoint.task)
        flat += values_for(model, said, chosen)

    fresh = ctx.table.unique_names(requested)
    columns = {name: (oids, value) for name, value in zip(fresh, flat, strict=True)}
    written = ctx.apply_results(columns, label="Permeability")
    message = (
        f"Predicted permeability for {plural(len(oids), 'row')}, "
        f"{plural(len(endpoints), 'endpoint')}"
    )
    if results.failed:
        message += f"; {plural(results.failed, 'structure')} could not be processed"
    ctx.show_message(message)
    return written
