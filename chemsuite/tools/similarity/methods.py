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
"""The Method section of Similarity Search: one row per family of methods, each switched on
by its box, with its settings in More Options; and what they make, an engine ``Search``."""

from __future__ import annotations

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox

from chemsuite.analysis.fingerprints import ASYMMETRIC_METRICS, DEFAULT_KIND, KINDS, METRICS
from chemsuite.app.widgets.dialogs import inline
from chemsuite.similarity import model as sm

SUBSTRUCTURES = {
    sm.MCS: "Common Substructure (MCS)",
    sm.SCAFFOLD: "Same Scaffold",
    sm.GENERIC_SCAFFOLD: "Same Generic Scaffold",
}
SHAPES = {sm.SHAPE: "Shape", sm.COLOR: "Shape + Features", sm.FIELD: "Field + Shape"}
PHARMACOPHORE = "pharm2d"


def _combo(choices: dict[str, str]) -> QComboBox:
    box = QComboBox()
    for key, title in choices.items():
        box.addItem(title, key)
    return box


def _weight() -> QDoubleSpinBox:
    box = QDoubleSpinBox()
    box.setRange(0.0, 1.0)
    box.setSingleStep(0.1)
    box.setDecimals(2)
    return box


class Methods(QObject):
    """The method rows (in *section*) and their settings (in *more*, folded)."""

    changed = Signal()

    def __init__(self, section, more) -> None:
        super().__init__(section)
        self.section, self.more = section, more

        self.use_fingerprint = QCheckBox()
        self.use_fingerprint.setChecked(True)
        self.kind = QComboBox()
        for kind in (k for k in KINDS.values() if k.listed):
            self.kind.addItem(kind.title, kind.key)
            if kind.help:
                self.kind.setItemData(self.kind.count() - 1, kind.help, Qt.ItemDataRole.ToolTipRole)
        self.kind.setCurrentIndex(self.kind.findData(DEFAULT_KIND))
        section.row("Fingerprint", inline(self.use_fingerprint, self.kind, stretch=1))
        self.metric = _combo({**METRICS, **ASYMMETRIC_METRICS})
        section.row(
            "Similarity",
            self.metric,
            help="Tanimoto, Dice and Cosine treat the row and the reference alike. Tversky "
            "weighs the bits only the reference has (α) and those only the row has (β), in "
            "More Options: α 1, β 0 asks how much of the reference the row contains.",
        )
        self.use_pharmacophore = QCheckBox("Gobbi 2D pharmacophore")
        section.row(
            "Pharmacophore",
            self.use_pharmacophore,
            help="Pairs and triplets of donors, acceptors, charges, rings and hydrophobes and the "
            "bonds between them, compared by Tanimoto.",
        )
        self.use_substructure = QCheckBox()
        self.substructure = _combo(SUBSTRUCTURES)
        section.row(
            "Substructure",
            inline(self.use_substructure, self.substructure, stretch=1),
            help="Common Substructure: the largest substructure the two share, as a share of "
            "their heavy atoms (an atom Tanimoto). Same Scaffold: 1 when the Bemis–Murcko "
            "scaffolds are the same, else 0; the generic one ignores elements and bond orders.",
        )
        self.use_3d = QCheckBox()
        self.shape = _combo(SHAPES)
        section.row(
            "3D",
            inline(self.use_3d, self.shape, stretch=1),
            help="Each row is superposed onto the reference (as Superpose does) and scored by "
            "the fit: shape overlap; shape plus donors, acceptors, rings and charges (0–2); or "
            "electrostatic field plus shape. The references must have 3D coordinates; a row "
            "with only 2D gets conformers first.",
        )
        self.consensus = _combo(sm.CONSENSUS)
        section.row(
            "Combine",
            self.consensus,
            help="One Consensus score from the methods: the mean of their scores (each on a 0–1 "
            "scale), the mean of their rank among the rows, or the mean after stretching each "
            "method's scores to fill 0–1.",
        )

        self.prefilter = QSpinBox()
        self.prefilter.setRange(0, 1_000_000)
        self.prefilter.setSpecialValueText("Every row")
        self.prefilter.setValue(500)
        more.row(
            "Prefilter",
            self.prefilter,
            help="The common substructure and 3D methods take seconds a row, so they score only "
            "this many rows, the ones the fingerprint ranks best; the rest are left blank.",
        )
        self.timeout = QSpinBox()
        self.timeout.setRange(1, 600)
        self.timeout.setSuffix(" s")
        self.timeout.setValue(2)
        more.row("MCS time limit", self.timeout, help="Per row: the best substructure so far.")
        self.conformers = QSpinBox()
        self.conformers.setRange(1, 300)
        self.conformers.setValue(30)
        more.row(
            "Conformers for 2D rows",
            self.conformers,
            help="Generated as Generate Conformers does (CONFORGE where installed, else ETKDG); "
            "a row with 3D is fitted as it is stored.",
        )
        self.color_weight = _weight()
        self.color_weight.setValue(0.5)
        more.row(
            "Feature weight",
            self.color_weight,
            help="The share of the features in the fit: 0 is shape only, 1 features only.",
        )
        self.field_weight = _weight()
        self.field_weight.setValue(0.5)
        more.row(
            "Field weight",
            self.field_weight,
            help="The share of the field in the fit: 0 is shape only, 1 field only.",
        )
        self.alpha = _weight()
        self.alpha.setValue(0.5)
        more.row("Tversky α", self.alpha, help="The weight of the bits only the reference has.")
        self.beta = _weight()
        self.beta.setValue(0.5)
        more.row("Tversky β", self.beta, help="The weight of the bits only the row has.")

        for box in (self.use_fingerprint, self.use_pharmacophore, self.use_substructure,
                    self.use_3d):  # fmt: skip
            box.toggled.connect(lambda _on: self._update())
        for combo in (self.kind, self.metric, self.substructure, self.shape, self.consensus):
            combo.currentIndexChanged.connect(lambda _i: self._update())
        for spin in (self.prefilter, self.timeout, self.conformers, self.color_weight,
                     self.field_weight, self.alpha, self.beta):  # fmt: skip
            spin.valueChanged.connect(lambda _v: self.changed.emit())
        self._update()

    # -- what is chosen -----------------------------------------------------------------------
    def kinds(self) -> list[str]:
        """The chosen methods' kinds, in the order their columns come (the first is primary)."""
        out = []
        if self.use_fingerprint.isChecked() or self.use_pharmacophore.isChecked():
            out.append(sm.FINGERPRINT)
        if self.use_substructure.isChecked():
            out.append(self.substructure.currentData())
        if self.use_3d.isChecked():
            out.append(self.shape.currentData())
        return out

    def methods(self) -> tuple[sm.Method, ...]:
        """The chosen methods as the engine takes them (raises ``ValueError`` on bad settings)."""
        out = []
        if self.use_fingerprint.isChecked():
            out.append(self._method(sm.FINGERPRINT, self.kind.currentData(),
                                    self.metric.currentData()))  # fmt: skip
        if self.use_pharmacophore.isChecked():
            out.append(self._method(sm.FINGERPRINT, PHARMACOPHORE, "tanimoto"))
        if self.use_substructure.isChecked():
            out.append(self._method(self.substructure.currentData()))
        if self.use_3d.isChecked():
            out.append(self._method(self.shape.currentData()))
        return tuple(out)

    def _method(self, kind: str, fingerprint: str = DEFAULT_KIND, metric: str = "tanimoto"):
        return sm.Method(
            kind=kind,
            fingerprint=fingerprint,
            metric=metric,
            alpha=self.alpha.value(),
            beta=self.beta.value(),
            mcs_timeout=self.timeout.value(),
            conformers=self.conformers.value(),
            color_weight=self.color_weight.value(),
            field_weight=self.field_weight.value(),
        )

    def search(self, fuse: str) -> sm.Search:
        """The search the rows describe (raises ``ValueError`` when it cannot be run)."""
        return sm.Search(self.methods(), fuse=fuse, consensus=self.consensus.currentData(),
                         prefilter=self.prefilter.value())  # fmt: skip

    def has_3d(self) -> bool:
        return self.use_3d.isChecked()

    def count(self) -> int:
        boxes = (self.use_fingerprint, self.use_pharmacophore, self.use_substructure, self.use_3d)
        return sum(b.isChecked() for b in boxes)

    def set_search(self, search: sm.Search) -> None:
        """Show *search*'s methods and settings (a search remembered from the last run)."""
        kinds = {m.kind for m in search.methods}
        prints = [m for m in search.methods if m.kind == sm.FINGERPRINT]
        plain = [m for m in prints if m.fingerprint != PHARMACOPHORE]
        self.use_fingerprint.setChecked(bool(plain))
        if plain:
            self.kind.setCurrentIndex(max(0, self.kind.findData(plain[0].fingerprint)))
            self.metric.setCurrentIndex(max(0, self.metric.findData(plain[0].metric)))
        self.use_pharmacophore.setChecked(any(m.fingerprint == PHARMACOPHORE for m in prints))
        found = kinds & set(SUBSTRUCTURES)
        self.use_substructure.setChecked(bool(found))
        if found:
            self.substructure.setCurrentIndex(self.substructure.findData(found.pop()))
        found = kinds & set(SHAPES)
        self.use_3d.setChecked(bool(found))
        if found:
            self.shape.setCurrentIndex(self.shape.findData(found.pop()))
        self.consensus.setCurrentIndex(max(0, self.consensus.findData(search.consensus)))
        self.prefilter.setValue(search.prefilter)
        first = search.methods[0]
        self.timeout.setValue(first.mcs_timeout)
        self.conformers.setValue(first.conformers)
        self.color_weight.setValue(first.color_weight)
        self.field_weight.setValue(first.field_weight)
        self.alpha.setValue(first.alpha)
        self.beta.setValue(first.beta)

    def _update(self) -> None:
        fingerprint = self.use_fingerprint.isChecked()
        kinds = set(self.kinds())
        tversky = fingerprint and self.metric.currentData() in ASYMMETRIC_METRICS
        self.kind.setEnabled(fingerprint)
        self.substructure.setEnabled(self.use_substructure.isChecked())
        self.shape.setEnabled(self.use_3d.isChecked())
        self.section.set_row_visible(self.metric, fingerprint)
        self.section.set_row_visible(self.consensus, self.count() > 1)
        more = self.more
        more.set_row_visible(self.prefilter, bool(kinds & sm.EXPENSIVE))
        more.set_row_visible(self.timeout, sm.MCS in kinds)
        more.set_row_visible(self.conformers, bool(kinds & sm.THREE_D))
        more.set_row_visible(self.color_weight, sm.COLOR in kinds)
        more.set_row_visible(self.field_weight, sm.FIELD in kinds)
        more.set_row_visible(self.alpha, tversky)
        more.set_row_visible(self.beta, tversky)
        self.changed.emit()
