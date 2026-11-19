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
"""Each map method's own settings: the rows it adds to the Method and More Options sections,
what it hands to ``embedding.embed``, and what stops it running.

One dialog shows one method, so each class builds only its own rows. Every row explains
itself in its tooltip.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox, QWidget

from chemsuite.analysis import embedding as em
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import Section, inline
from chemsuite.plotting.data import BLANK_LABEL, categories, colorable_columns, numbers

#: the PCA axes offered: (first, second) component, 1-based
COMPONENT_PAIRS = ((1, 2), (1, 3), (2, 3), (3, 4))
#: default neighbors per method
NEIGHBOURS = {"umap": 15, "pacmap": 10}
NONE = "(none)"


def spin(low: int, high: int, value: int, step: int = 1, special: str = "") -> QSpinBox:
    box = QSpinBox()
    box.setRange(low, high)
    box.setSingleStep(step)
    box.setValue(value)
    if special:
        box.setSpecialValueText(special)
    box.setGroupSeparatorShown(True)
    return box


def number(low: float, high: float, value: float, step: float, decimals: int = 2):
    box = QDoubleSpinBox()
    box.setDecimals(decimals)
    box.setRange(low, high)
    box.setSingleStep(step)
    box.setValue(value)
    return box


def choice(items: dict[str, str], current: str) -> QComboBox:
    box = QComboBox()
    for key, title in items.items():
        box.addItem(title, key)
    box.setCurrentIndex(max(0, box.findData(current)))
    return box


def fit_bit_distance(box: QComboBox, bits: bool, memory: dict) -> None:
    """Bit distances only where they apply: switch to Euclidean when they stop fitting, and
    back to the bit distance once they fit again (*memory* holds what was switched away)."""
    model = box.model()
    for row in range(box.count()):
        if box.itemData(row) in em.BIT_METRICS:
            model.item(row).setEnabled(bits)
    current = box.currentData()
    blocked = box.blockSignals(True)
    if not bits and current in em.BIT_METRICS:
        memory[box] = current
        box.setCurrentIndex(box.findData("euclidean"))
    elif bits and current == "euclidean" and box in memory:
        box.setCurrentIndex(box.findData(memory.pop(box)))
    elif bits:
        memory.pop(box, None)
    box.blockSignals(blocked)


class MethodOptions:
    """The settings of one method. *method* and *more* are the dialog's sections."""

    key = ""
    fits = True  # can fit on some rows and place the rest

    def __init__(self, ctx: ToolContext, method: Section, more: Section) -> None:
        self.ctx = ctx
        self.method = method
        self.more = more
        self._bit_choice: dict[QComboBox, str] = {}

    def signals(self) -> list[Signal]:
        """What changes the dialog must react to."""
        return []

    def fields(self) -> dict[str, QWidget]:
        """The settings remembered between runs."""
        return {}

    def update(self, *, bits: bool, three: bool, columns_only: bool) -> None:
        """Show the rows that apply now (*bits*: fingerprint bits alone)."""

    def problem(self, *, rows: int, fitting: bool, bits: bool) -> str:
        return ""

    def options(self, oids: np.ndarray) -> dict:
        """What ``embed`` takes, for the rows *oids*."""
        return {}

    def seed_row(self) -> QSpinBox:
        self.seed = spin(-1, 2_147_483_647, 42, special="Random")
        self.more.row("Seed", self.seed,
                      help="The same seed gives the same map every run. Random is faster but "
                      "differs each time.")  # fmt: skip
        return self.seed

    def seed_value(self) -> int | None:
        return None if self.seed.value() < 0 else self.seed.value()


class PcaOptions(MethodOptions):
    key = "pca"

    def __init__(self, ctx, method, more) -> None:
        super().__init__(ctx, method, more)
        self.components = QComboBox()
        for first, second in COMPONENT_PAIRS:
            self.components.addItem(f"PC{first} and PC{second}", f"{first},{second}")
        method.row("Components", self.components,
                   help="Which principal components make the two axes. Later components show "
                   "differences the first two miss.")  # fmt: skip
        self.write = spin(2, 20, 2, special="Plotted only")
        more.row("Write components", self.write,
                 help="Also write the scores of PC1 to this many components as columns "
                 "(PC1, PC2, …), for other plots and models.")  # fmt: skip
        self.loadings = QCheckBox("Show what each column contributes")
        more.row("Loadings", self.loadings,
                 help="A table of each number column's weight on each component. Only for "
                 "number columns without a fingerprint: single bits mean little.")  # fmt: skip

    def signals(self):
        return [self.components.currentIndexChanged, self.write.valueChanged,
                self.loadings.toggled]  # fmt: skip

    def fields(self):
        return {"components": self.components, "write": self.write, "loadings": self.loadings}

    def update(self, *, bits, three, columns_only):
        self.method.set_row_visible(self.components, not three)
        self.loadings.setEnabled(columns_only)

    def pair(self) -> tuple[int, int]:
        first, second = self.components.currentData().split(",")
        return int(first), int(second)

    def options(self, oids):
        written = self.write.value()
        return {"components": self.pair(), "write": written if written > 2 else 0}


class UmapOptions(MethodOptions):
    key = "umap"

    def __init__(self, ctx, method, more) -> None:
        super().__init__(ctx, method, more)
        self.neighbours = spin(2, 200, NEIGHBOURS["umap"])
        method.row("Neighbors", self.neighbours,
                   help="Fewer keeps local detail; more shows the overall layout.")  # fmt: skip
        self.min_dist = number(0.0, 0.99, 0.1, 0.05)
        method.row("Min. distance", self.min_dist, help="How tightly points may pack together.")
        self.guide = QComboBox()
        self.guide.addItem(NONE, None)
        for name in colorable_columns(ctx.table):
            self.guide.addItem(name, name)
        method.row("Guide by", self.guide,
                   help="Supervised UMAP: compounds with the same value of this column are "
                   "pulled together. A number column is cut into groups; rows without a value "
                   "are placed by their features alone.")  # fmt: skip
        self.groups = spin(2, 20, em.GUIDE_GROUPS)
        method.row(
            "Groups",
            self.groups,
            help="How many groups of about equal size a number column is cut into.",
        )
        self.guide_weight = number(0.0, 1.0, 0.5, 0.1)
        method.row(
            "Guide weight",
            self.guide_weight,
            help="0: the features alone shape the map; 1: the guide column does.",
        )
        self.distance = choice(em.UMAP_METRICS, "jaccard")
        more.row("Distance", self.distance,
                 help="How unlike two compounds are. Jaccard, Dice, Rogers-Tanimoto and "
                 "Hamming work on fingerprint bits alone; anything else is measured on the "
                 "numbers.")  # fmt: skip
        self.spread = number(0.1, 10.0, 1.0, 0.1)
        more.row("Spread", self.spread,
                 help="The overall scale of the map; with Min. distance, how clumped the "
                 "islands are.")  # fmt: skip
        self.epochs = spin(0, 5000, 0, 50, special="Automatic")
        more.row(
            "Epochs",
            self.epochs,
            help="Optimization rounds. Automatic: 500 for small sets, 200 for large ones.",
        )
        self.fit_max = spin(1_000, 1_000_000, em.UMAP_FIT_MAX, 10_000)
        more.row("Fit on up to", self.fit_max,
                 help="Larger sets are fitted on a random sample of this many rows and the "
                 "rest placed onto that map: faster, and nearly the same picture.")  # fmt: skip
        self.init = choice(em.UMAP_INITS, "spectral")
        more.row(
            "Start from",
            self.init,
            help="The layout optimization starts from. Spectral follows the neighbor "
            "graph; PCA keeps the overall layout steadier; Random is the plainest.",
        )
        self.densmap = QCheckBox("Keep how dense each region is")
        more.row("densMAP", self.densmap,
                 help="densMAP: tight groups stay small and loose ones spread out, as in the "
                 "data. It cannot place rows, so it fits every row.")  # fmt: skip
        self.density_weight = number(0.0, 10.0, 2.0, 0.5)
        more.row(
            "Density weight",
            self.density_weight,
            help="How strongly densMAP holds the densities against the layout.",
        )
        self.local_connectivity = number(0.1, 20.0, 1.0, 0.5)
        more.row("Local connectivity", self.local_connectivity,
                 help="How many nearest neighbors each point is surely joined to; raise it "
                 "for noisy data that falls apart into specks.")  # fmt: skip
        self.repulsion = number(0.1, 20.0, 1.0, 0.5)
        more.row("Repulsion", self.repulsion,
                 help="How hard unrelated points push apart; higher separates the islands "
                 "further.")  # fmt: skip
        self.negative_samples = spin(1, 50, 5)
        more.row("Negative samples", self.negative_samples,
                 help="Unrelated points pushed away per step: more is slower and a cleaner "
                 "separation.")  # fmt: skip
        self.learning_rate = number(0.01, 10.0, 1.0, 0.1)
        more.row("Learning rate", self.learning_rate, help="Step size of the optimization.")
        self.seed_row()
        self.seed.setToolTip(self.seed.toolTip() + " UMAP then uses one processor core.")

    def signals(self):
        return [self.distance.currentIndexChanged, self.min_dist.valueChanged,
                self.spread.valueChanged, self.guide.currentIndexChanged,
                self.densmap.toggled, self.fit_max.valueChanged]  # fmt: skip

    def fields(self):
        return {"neighbours": self.neighbours, "min_dist": self.min_dist,
                "groups": self.groups, "guide_weight": self.guide_weight,
                "distance": self.distance, "spread": self.spread, "epochs": self.epochs,
                "fit_max": self.fit_max, "init": self.init, "densmap": self.densmap,
                "density_weight": self.density_weight,
                "local_connectivity": self.local_connectivity, "repulsion": self.repulsion,
                "negative_samples": self.negative_samples,
                "learning_rate": self.learning_rate, "seed": self.seed}  # fmt: skip

    def _numeric_guide(self) -> bool:
        name = self.guide.currentData()
        return name is not None and self.ctx.table.kind(name) == "number"

    def update(self, *, bits, three, columns_only):
        fit_bit_distance(self.distance, bits, self._bit_choice)
        guided = self.guide.currentData() is not None
        self.method.set_row_visible(self.groups, self._numeric_guide())
        self.method.set_row_visible(self.guide_weight, guided)
        self.more.set_row_visible(self.density_weight, self.densmap.isChecked())

    def problem(self, *, rows, fitting, bits):
        metric = self.distance.currentData()
        if metric in em.BIT_METRICS and not bits:
            return f"{em.UMAP_METRICS[metric]} distance needs fingerprint bits alone."
        if self.min_dist.value() > self.spread.value():
            return "Min. distance must not be larger than the spread."
        if self.densmap.isChecked() and fitting:
            return "densMAP cannot place rows onto a map: fit on all rows."
        if self.densmap.isChecked() and rows > self.fit_max.value():
            return f"densMAP fits every row: map at most {self.fit_max.value():,} rows."
        return ""

    def guide_labels(self, oids: np.ndarray) -> np.ndarray | None:
        name = self.guide.currentData()
        if name is None or name not in self.ctx.table:
            return None
        if self._numeric_guide():
            return em.guide_labels(numbers(self.ctx.table, name, oids), self.groups.value())
        found = categories(self.ctx.table, name, oids, limit=10_000)
        blank = [i for i, label in enumerate(found.labels) if label == BLANK_LABEL]
        return np.where(np.isin(found.codes, blank), -1, found.codes)

    def options(self, oids):
        return {
            "neighbours": self.neighbours.value(),
            "min_dist": self.min_dist.value(),
            "spread": self.spread.value(),
            "metric": self.distance.currentData(),
            "epochs": self.epochs.value() or None,
            "fit_max": self.fit_max.value(),
            "seed": self.seed_value(),
            "guide": self.guide_labels(oids),
            "guide_weight": self.guide_weight.value(),
            "densmap": self.densmap.isChecked(),
            "density_weight": self.density_weight.value(),
            "init": self.init.currentData(),
            "local_connectivity": self.local_connectivity.value(),
            "repulsion": self.repulsion.value(),
            "negative_samples": self.negative_samples.value(),
            "learning_rate": self.learning_rate.value(),
        }


def reduce_row(more: Section, what: str) -> QSpinBox:
    box = spin(0, 500, 0, special="Off")
    more.row("PCA dimensions", box,
             help=f"Off: the map is made from the features themselves ({what}). A number: "
             "from that many principal components; faster on large sets, but bit distances "
             "are then not available.")  # fmt: skip
    return box


class TsneOptions(MethodOptions):
    key = "tsne"
    fits = False

    def __init__(self, ctx, method, more) -> None:
        super().__init__(ctx, method, more)
        self.perplexity = number(2, 200, 30, 1)
        method.row("Perplexity", self.perplexity,
                   help="Roughly how many neighbors each point attends to.")  # fmt: skip
        self.learning_rate = number(0, 5000, 0, 50)
        self.learning_rate.setSpecialValueText("Automatic")
        more.row("Learning rate", self.learning_rate,
                 help="Step size. Automatic scales it with the number of rows; raise it if "
                 "the map looks like a compressed ball.")  # fmt: skip
        self.iterations = spin(250, 10_000, 1000, 250)
        more.row("Iterations", self.iterations, help="More can settle a large map further.")
        self.exaggeration = number(1.0, 100.0, 12.0, 1.0)
        more.row("Exaggeration", self.exaggeration,
                 help="Early exaggeration: higher pushes the groups further apart at the "
                 "start.")  # fmt: skip
        self.reduce_to = reduce_row(more, "fingerprint bits by Jaccard")
        self.distance = choice(em.TSNE_METRICS, "euclidean")
        more.row("Distance", self.distance,
                 help="How unlike two compounds are, for number columns, counts or PCA "
                 "dimensions. Fingerprint bits are always compared by Jaccard.")  # fmt: skip
        self.exact = QCheckBox("Exact")
        more.row("Algorithm", self.exact,
                 help="Exact measures every pair: slower and only up to "
                 f"{em.TSNE_EXACT_MAX:,} rows. Otherwise Barnes-Hut approximates far "
                 "points.")  # fmt: skip
        self.angle = number(0.1, 1.0, 0.5, 0.05)
        more.row(
            "Angle", self.angle, help="Barnes-Hut's trade-off: lower is more accurate and slower."
        )
        self.init = choice({"random": "Random", "pca": "PCA"}, "random")
        more.row("Start from", self.init,
                 help="PCA keeps the overall layout more stable; Random is the classic "
                 "t-SNE.")  # fmt: skip
        self.seed_row()

    def signals(self):
        return [self.reduce_to.valueChanged, self.exact.toggled]

    def fields(self):
        return {"perplexity": self.perplexity, "learning_rate": self.learning_rate,
                "iterations": self.iterations, "exaggeration": self.exaggeration,
                "reduce_to": self.reduce_to, "distance": self.distance, "exact": self.exact,
                "angle": self.angle, "init": self.init, "seed": self.seed}  # fmt: skip

    def update(self, *, bits, three, columns_only):
        on_bits = bits and self.reduce_to.value() == 0
        self.more.set_row_visible(self.distance, not on_bits)
        self.more.set_row_visible(self.angle, not self.exact.isChecked())

    def problem(self, *, rows, fitting, bits):
        if self.exact.isChecked() and rows > em.TSNE_EXACT_MAX:
            return f"Exact t-SNE takes at most {em.TSNE_EXACT_MAX:,} rows."
        return ""

    def options(self, oids):
        return {
            "perplexity": self.perplexity.value(),
            "learning_rate": self.learning_rate.value() or None,
            "iterations": self.iterations.value(),
            "exaggeration": self.exaggeration.value(),
            "reduce_to": self.reduce_to.value(),
            "init": self.init.currentData(),
            "seed": self.seed_value(),
            "metric": self.distance.currentData(),
            "exact": self.exact.isChecked(),
            "angle": self.angle.value(),
        }


class PacmapOptions(MethodOptions):
    key = "pacmap"

    def __init__(self, ctx, method, more) -> None:
        super().__init__(ctx, method, more)
        self.neighbours = spin(2, 200, NEIGHBOURS["pacmap"])
        method.row("Neighbors", self.neighbours,
                   help="Fewer keeps local detail; more shows the overall layout.")  # fmt: skip
        self.variant = choice(em.PACMAP_VARIANTS, "pacmap")
        method.row("Variant", self.variant,
                   help="LocalMAP (Wang et al., 2025) sharpens local groups by letting go of "
                   "neighbor pairs that end up far apart.")  # fmt: skip
        self.mn_ratio = number(0.0, 10.0, 0.5, 0.1)
        more.row("Mid-near pairs", self.mn_ratio,
                 help="Mid-near pairs per neighbor pair; more keeps the overall layout "
                 "(global structure) better.")  # fmt: skip
        self.fp_ratio = number(0.0, 20.0, 2.0, 0.5)
        more.row("Further pairs", self.fp_ratio,
                 help="Further pairs per neighbor pair; more pushes separate groups further "
                 "apart.")  # fmt: skip
        self.distance = choice(em.PACMAP_METRICS, "jaccard")
        more.row("Distance", self.distance,
                 help="Jaccard and Hamming work on fingerprint bits alone, with PCA "
                 "dimensions Off. Angular is the cosine distance.")  # fmt: skip
        self.reduce_to = reduce_row(more, "fingerprint bits")
        self.local_threshold = number(0.1, 100.0, 10.0, 1.0)
        more.row(
            "Local distance limit",
            self.local_threshold,
            help="LocalMAP lets go of neighbor pairs farther apart than this on the map.",
        )
        self.learning_rate = number(0.01, 10.0, 1.0, 0.1)
        more.row("Learning rate", self.learning_rate, help="Step size of the optimization.")
        self.phases = [spin(0, 5000, n, 50) for n in (100, 100, 250)]
        more.row("Iterations", inline(*self.phases),
                 help="The three phases: settling the overall layout, refining it, then "
                 "tightening local groups.")  # fmt: skip
        self.init = choice({"pca": "PCA", "random": "Random"}, "pca")
        more.row("Start from", self.init, help="Where the layout starts.")
        self.seed_row()

    def signals(self):
        return [self.distance.currentIndexChanged, self.reduce_to.valueChanged,
                self.variant.currentIndexChanged]  # fmt: skip

    def fields(self):
        return {"neighbours": self.neighbours, "variant": self.variant,
                "mn_ratio": self.mn_ratio, "fp_ratio": self.fp_ratio,
                "distance": self.distance, "reduce_to": self.reduce_to,
                "local_threshold": self.local_threshold, "learning_rate": self.learning_rate,
                "phase1": self.phases[0], "phase2": self.phases[1], "phase3": self.phases[2],
                "init": self.init, "seed": self.seed}  # fmt: skip

    def bits_fit(self, bits: bool) -> bool:
        return bits and self.reduce_to.value() == 0

    def update(self, *, bits, three, columns_only):
        fit_bit_distance(self.distance, self.bits_fit(bits), self._bit_choice)
        local = self.variant.currentData() == "localmap"
        self.more.set_row_visible(self.local_threshold, local)

    def problem(self, *, rows, fitting, bits):
        distance = self.distance.currentData()
        if distance in em.BIT_METRICS and not self.bits_fit(bits):
            return (f"{em.PACMAP_METRICS[distance]} distance needs fingerprint bits alone, "
                    "with PCA dimensions Off.")  # fmt: skip
        if distance in em.BIT_METRICS and fitting:
            return (f"PaCMAP cannot place rows by {em.PACMAP_METRICS[distance]} distance: "
                    "choose another distance, or fit on all rows.")  # fmt: skip
        return ""

    def options(self, oids):
        return {
            "neighbours": self.neighbours.value(),
            "mn_ratio": self.mn_ratio.value(),
            "fp_ratio": self.fp_ratio.value(),
            "distance": self.distance.currentData(),
            "reduce_to": self.reduce_to.value(),
            "seed": self.seed_value(),
            "variant": self.variant.currentData(),
            "learning_rate": self.learning_rate.value(),
            "iterations": tuple(p.value() for p in self.phases),
            "init": self.init.currentData(),
            "local_threshold": self.local_threshold.value(),
        }


OPTIONS = {"pca": PcaOptions, "umap": UmapOptions, "tsne": TsneOptions, "pacmap": PacmapOptions}
