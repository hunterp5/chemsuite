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
"""Stepping through a protein's models: the NMR models of a file now, trajectory frames later.

A slim bar under the 3D view, shown only while a visible protein has more than one model:
previous and next, the model shown (type a number to jump to it) of how many, play (a loop),
and a slider to scrub. When several shown proteins have models, a list picks which one the bar
steps. It reads and changes the session only, so the tree, sequence and view follow.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QEvent, QSignalBlocker, QSize, Qt, QTimer
from PySide6.QtGui import QAction
from PySide6.QtWidgets import QComboBox, QLabel, QSizePolicy, QSlider, QSpinBox, QToolBar

from chemsuite.app.shell.glyphs import glyph
from chemsuite.protein.display import Display

from .session import ProteinSession

#: models shown per second while playing
FPS = 10


class ModelsBar(QToolBar):
    def __init__(self, session: ProteinSession, parent=None) -> None:
        super().__init__(parent)
        self.session = session
        self.setIconSize(QSize(15, 15))
        self.setMovable(False)
        self._glyphs: dict[QAction, str] = {}
        self._ids: list[str] = []

        self.protein = QComboBox()
        self.protein.setToolTip("Which protein's models to step through")
        self.protein.currentIndexChanged.connect(self.refresh)
        self._protein_action = self.addWidget(self.protein)
        self.previous = self._add(
            "previous", "Previous Model", "one back ([)", lambda: self.step(-1)
        )
        self.position = QSpinBox()
        self.position.setKeyboardTracking(False)
        self.position.setButtonSymbols(QSpinBox.ButtonSymbols.NoButtons)
        self.position.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.position.setToolTip("Model shown; type a number to jump to it")
        self.position.valueChanged.connect(lambda value: self.go(value - 1))
        self.addWidget(self.position)
        self.count = QLabel()
        self.count.setContentsMargins(4, 0, 4, 0)
        self.addWidget(self.count)
        self.next = self._add("next", "Next Model", "one on (])", lambda: self.step(1))
        self.play = self._add("play", "Play", "step through the models in a loop",
                              self.set_playing, checkable=True)  # fmt: skip
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setToolTip("Drag to scrub through the models")
        self.slider.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        self.slider.valueChanged.connect(self.go)
        self.addWidget(self.slider)

        self._timer = QTimer(self)
        self._timer.setInterval(1000 // FPS)
        self._timer.timeout.connect(self._tick)
        session.proteinsChanged.connect(self.refresh)
        session.sceneChanged.connect(self.refresh)
        self._restyle()
        self.refresh()

    def _add(self, name: str, text: str, tip: str, slot: Callable, *, checkable=False) -> QAction:
        action = self.addAction(text)
        action.setToolTip(f"{text}: {tip}")
        action.setCheckable(checkable)
        (action.toggled if checkable else action.triggered).connect(slot)
        self._glyphs[action] = name
        return action

    # which protein
    def candidates(self) -> list[Display]:
        """The shown proteins with more than one model."""
        return [d for d in self.session.displays if d.visible and d.n_models > 1]

    def current(self) -> Display | None:
        shown = self.session.display(self.protein.currentData() or "")
        return shown if shown is not None and shown in self.candidates() else None

    def refresh(self) -> None:
        candidates = self.candidates()
        ids = [d.id for d in candidates]
        if ids != self._ids:
            with QSignalBlocker(self.protein):
                keep = self.protein.currentData()
                self.protein.clear()
                for shown in candidates:
                    self.protein.addItem(shown.name, shown.id)
                self.protein.setCurrentIndex(max(self.protein.findData(keep), 0))
            self._ids = ids
        self._protein_action.setVisible(len(candidates) > 1)
        self.setVisible(bool(candidates))
        shown = self.current()
        if shown is None:
            self.play.setChecked(False)
            return
        n, k = shown.n_models, shown.model
        with QSignalBlocker(self.position), QSignalBlocker(self.slider):
            self.position.setRange(1, n)
            self.position.setValue(k + 1)
            self.slider.setRange(0, n - 1)
            self.slider.setValue(k)
        self.count.setText(f"of {n}")
        self.previous.setEnabled(k > 0)
        self.next.setEnabled(k < n - 1)

    # stepping
    def go(self, model: int) -> None:
        shown = self.current()
        if shown is not None:
            self.session.set_model(shown.id, model)

    def step(self, delta: int) -> None:
        shown = self.current()
        if shown is not None:
            self.go(min(max(shown.model + delta, 0), shown.n_models - 1))

    def set_playing(self, on: bool) -> None:
        if on and self.current() is not None:
            self._timer.start()
        else:
            self._timer.stop()
        self.play.setText("Pause" if self._timer.isActive() else "Play")
        self._glyphs[self.play] = "pause" if self._timer.isActive() else "play"
        self._restyle()

    @property
    def playing(self) -> bool:
        return self._timer.isActive()

    def _tick(self) -> None:
        shown = self.current()
        if shown is None:
            self.play.setChecked(False)
            return
        self.go((shown.model + 1) % shown.n_models)

    # the theme
    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange:
            self._restyle()

    def _restyle(self) -> None:
        ink = self.palette().text().color()
        for action, name in self._glyphs.items():
            action.setIcon(glyph(name, ink))
