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
"""The QSAR Models panel: what this document's models are, and how good they are.

Four tabs, because four questions get asked of a model in roughly this order. *Scores* — how
well does it do, and on what? *Leaderboard* — what else was tried, and by how much did this
one win? (A narrow win is worth knowing: it means the choice of model hardly matters, and the
features are doing the work.) *Features* — what is it actually using? *History* — what has
been done to it since.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
)

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.qsar.features import title_of
from chemsuite.qsar.metrics import TITLES
from chemsuite.qsar.model import SUFFIX, Model

from .part import PART_KEY, QsarModels


def _table(headers: list[str]) -> QTableWidget:
    table = QTableWidget(0, len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setAlternatingRowColors(True)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    table.horizontalHeader().setStretchLastSection(True)
    return table


def _fill(table: QTableWidget, rows: list[list[str]]) -> None:
    table.setRowCount(len(rows))
    for r, row in enumerate(rows):
        for c, text in enumerate(row):
            item = QTableWidgetItem(text)
            if c and _is_number(text):
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            table.setItem(r, c, item)


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def _number(value: float) -> str:
    return "—" if value != value else f"{value:.3g}"


def _setting(value: object) -> str:
    """A searched hyperparameter as a person reads it.

    The search samples from continuous ranges, so a setting arrives as a full-precision float
    (``learning_rate=0.02771473968728253``), which tells a reader nothing the third digit has
    not already told them.
    """
    return _number(value) if isinstance(value, float) else str(value)


class ModelsPanel(ToolPanel):
    """The models this document holds."""

    tool_id = "qsar_models"

    def __init__(self, ctx: ToolContext) -> None:
        super().__init__()
        self.ctx = ctx
        self.part: QsarModels = ctx.document_part(
            PART_KEY, lambda saved: QsarModels(ctx.table, saved)
        )

        self.which = QComboBox()
        self.which.setMinimumWidth(220)
        self.which.currentIndexChanged.connect(lambda _i: self._show())
        self.headline = QLabel()
        self.headline.setTextFormat(Qt.TextFormat.PlainText)

        self.scores = _table(["", "Out of fold"])
        self.leaderboard = _table(["Model", "Score", "Settings"])
        self.features = _table(["Feature", "Kind"])
        self.history = _table(["Date", "What", "Rows", "Score"])
        self.tabs = QTabWidget()
        for tab, title in (
            (self.scores, "Scores"),
            (self.leaderboard, "Leaderboard"),
            (self.features, "Features"),
            (self.history, "History"),
        ):
            self.tabs.addTab(tab, title)

        self.apply_button = QPushButton("Apply to Table…")
        self.apply_button.setToolTip("Predict with this model for rows in this document.")
        self.apply_button.clicked.connect(lambda: self.ctx.open_tool("qsar_apply"))
        self.refine_button = QPushButton("Refine…")
        self.refine_button.setToolTip("Fit this model again on its rows plus newer data.")
        self.refine_button.clicked.connect(lambda: self.ctx.open_tool("qsar_refine"))
        self.save_button = QPushButton("Save Model…")
        self.save_button.setToolTip(
            f"Write this model to a {SUFFIX} file, to use in another document or send to a "
            "colleague. The file carries the compounds it was trained on."
        )
        self.save_button.clicked.connect(self.save_model)
        self.export_button = QPushButton("Export for Python…")
        self.export_button.setToolTip(
            "Write the fitted scikit-learn object as a pickle, for use in a script outside "
            "chemsuite. chemsuite never reads a pickle back."
        )
        self.export_button.clicked.connect(self.export_for_python)
        self.forget_button = QPushButton("Remove")
        self.forget_button.setToolTip("Take this model out of the document.")
        self.forget_button.clicked.connect(self.forget_model)

        top = QHBoxLayout()
        top.setSpacing(8)
        top.addWidget(self.which, 1)
        top.addWidget(self.forget_button)
        actions = QHBoxLayout()
        actions.setSpacing(8)
        actions.addWidget(self.apply_button)
        actions.addWidget(self.refine_button)
        actions.addWidget(self.save_button)
        actions.addWidget(self.export_button)
        actions.addStretch(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        layout.addLayout(top)
        layout.addWidget(self.headline)
        layout.addWidget(self.tabs, 1)
        layout.addLayout(actions)
        self.reload()

    def title(self) -> str:
        return "QSAR Models"

    # -- filling ---------------------------------------------------------------------------
    def reload(self) -> None:
        """Bring the list of models up to date (a fit finished, or an undo took one away)."""
        kept = self.part.kept()
        wanted = self.which.currentData()
        self.which.blockSignals(True)
        self.which.clear()
        for entry in reversed(kept):
            self.which.addItem(entry.title, entry.asset)
        if wanted is not None:
            found = self.which.findData(wanted)
            if found >= 0:
                self.which.setCurrentIndex(found)
        self.which.blockSignals(False)
        self._show()

    def current(self) -> Model | None:
        asset = self.which.currentData()
        return None if asset is None else self.part.model(str(asset))

    def _show(self) -> None:
        model = self.current()
        for button in (
            self.apply_button,
            self.refine_button,
            self.save_button,
            self.export_button,
            self.forget_button,
        ):
            button.setEnabled(model is not None)
        if model is None:
            self.headline.setText("This document holds no model yet.")
            for table in (self.scores, self.leaderboard, self.features, self.history):
                table.setRowCount(0)
            return
        task = "class" if model.recipe.task == "classification" else "number"
        self.headline.setText(
            f"{model.candidate} · {len(model.training):,} compounds · "
            f"{model.features} features · predicts a {task}"
        )
        _fill(
            self.scores,
            [[TITLES.get(key, key), _number(value)] for key, value in model.scores.items()],
        )
        _fill(
            self.leaderboard,
            [
                [
                    str(entry.get("title", entry.get("candidate", ""))),
                    _number(float(entry.get("scores", {}).get(model.objective, float("nan")))),
                    ", ".join(
                        f"{k}={_setting(v)}" for k, v in dict(entry.get("params", {})).items()
                    ),
                ]
                for entry in model.leaderboard
            ],
        )
        _fill(self.features, [[title_of(n), _kind_of(n)] for n in model.kept_names()])
        _fill(self.history, [_history_row(entry, model.objective) for entry in model.lineage])

    # -- the actions -----------------------------------------------------------------------
    def save_model(self) -> None:
        model = self.current()
        if model is None:
            return
        suggested = f"{model.recipe.target}{SUFFIX}".replace("/", "-")
        chosen, _filter = QFileDialog.getSaveFileName(
            self, "Save Model", suggested, f"chemsuite model (*{SUFFIX})"
        )
        if not chosen:
            return
        path = model.save(Path(chosen))
        self.ctx.show_message(f"Saved {path.name}")

    def export_for_python(self) -> None:
        """Write the scikit-learn object as a pickle, for use outside chemsuite.

        Write-only on purpose: nothing in chemsuite reads it back, because unpickling runs
        whatever the file says. The warning is a real one — a pickle from someone else is
        code, not data — so the user sees it before the file dialog, not after.
        """
        model = self.current()
        if model is None:
            return
        from .export_python import export_for_python, exportable

        why = exportable(self.part, str(self.which.currentData()))
        if why:
            QMessageBox.information(self, "Export for Python", why)
            return
        suggested = f"{model.recipe.target}.joblib".replace("/", "-")
        chosen, _filter = QFileDialog.getSaveFileName(
            self, "Export for Python", suggested, "joblib pickle (*.joblib)"
        )
        if not chosen:
            return
        written = export_for_python(self.part, str(self.which.currentData()), Path(chosen))
        self.ctx.show_message(f"Wrote {written.name} and its manifest")

    def forget_model(self) -> None:
        asset = self.which.currentData()
        model = self.current()
        if asset is None or model is None:
            return
        answer = QMessageBox.question(
            self,
            "Remove Model",
            f"Remove {model.title} from this document? Its columns stay in the table.",
            QMessageBox.StandardButton.Cancel | QMessageBox.StandardButton.Yes,
            QMessageBox.StandardButton.Cancel,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        op = self.part.forget(str(asset))
        if op is not None:
            self.ctx.edit(f"Remove {model.title}", [op])
        self.reload()

    # -- persistence -----------------------------------------------------------------------
    def state(self) -> dict:
        return {"model": self.which.currentData(), "tab": self.tabs.currentIndex()}

    def set_state(self, state: dict) -> None:
        asset = state.get("model")
        if asset is not None:
            found = self.which.findData(asset)
            if found >= 0:
                self.which.setCurrentIndex(found)
        tab = state.get("tab")
        if isinstance(tab, int) and 0 <= tab < self.tabs.count():
            self.tabs.setCurrentIndex(tab)


def show_panel(ctx: ToolContext) -> ModelsPanel | None:
    """The document's models panel, opened if it is not already, and brought up to date.

    The one place a panel is opened, so a fit and a refit cannot drift apart, and so the
    ``add_panel`` call sits beside the class that carries ``tool_id`` (see
    ``tests/app/document/test_session.py``).
    """
    panel = ctx.find_panel(ModelsPanel.tool_id)
    if panel is None:
        panel = ModelsPanel(ctx)
        ctx.add_panel(panel, panel.title())
    if isinstance(panel, ModelsPanel):
        panel.reload()
        return panel
    return None


def _history_row(entry: dict, objective: str) -> list[str]:
    """One line of a model's history, as the History tab shows it."""
    scores = dict(entry.get("scores", {}))
    return [
        str(entry.get("date", "")),
        str(entry.get("note", "")),
        f"{int(entry.get('rows', 0)):,}",
        _number(float(scores.get(objective, float("nan")))),
    ]


def _kind_of(name: str) -> str:
    """What a pool column's name says it came from (see ``qsar.features.title_of``)."""
    prefix, _, _rest = name.partition(":")
    return {"d": "Descriptor", "f": "Fingerprint", "c": "Column"}.get(prefix, "")
