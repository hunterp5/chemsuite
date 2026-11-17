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
"""The minimisation's settings as dialog fields: Tools ▸ Minimise… and Prepare's
Minimisation section share them."""

from __future__ import annotations

from PySide6.QtWidgets import QComboBox, QDoubleSpinBox

from chemsuite.app.widgets.dialogs import Section, Segmented
from chemsuite.simulate import spec
from chemsuite.simulate.spec import MinimiseSettings

from . import minimise

FORCE_FIELDS = ((spec.FF14SB, "ff14SB"), (spec.FF19SB, "ff19SB"))
SOLVENTS = ((spec.GBN2, "GBn2"), (spec.VACUUM, "Vacuum"))
RESTRAINT_NAMES = ("Loose", "Normal", "Firm")
RESTRAINTS = tuple(
    (str(x), name) for x, name in zip(spec.RESTRAINT_SCALES, RESTRAINT_NAMES, strict=True)
)
AUTO = ""  # the platform: the fastest that works


class MinimiseFields:
    """The minimisation's settings in one section of a dialog."""

    def __init__(self, section: Section, settings: MinimiseSettings | None = None) -> None:
        settings = settings or MinimiseSettings()
        self.section = section
        self.forcefield = section.row(
            "Force Field", Segmented(FORCE_FIELDS, settings.protein),
            help="ff14SB with TIP3P water, or ff19SB with OPC water; small molecules take "
            "OpenFF Sage either way.",
        )  # fmt: skip
        self.solvent = section.row(
            "Solvent", Segmented(SOLVENTS, settings.solvent),
            help="GBn2 implicit solvent, or none (vacuum). Kept crystal waters stay as "
            "explicit waters in both.",
        )  # fmt: skip
        self.restraints = section.row(
            "Restraints", Segmented(RESTRAINTS, str(settings.restraints)),
            help="How hard the deposited heavy atoms are held where they started: half, as "
            "planned, or twice as hard. Atoms Prepare built are never held.",
        )  # fmt: skip
        self.rmsd_limit = QDoubleSpinBox()
        self.rmsd_limit.setRange(0.1, 2.0)
        self.rmsd_limit.setSingleStep(0.05)
        self.rmsd_limit.setDecimals(2)
        self.rmsd_limit.setSuffix(" Å")
        self.rmsd_limit.setValue(settings.rmsd_limit)
        section.row("RMSD Limit", self.rmsd_limit,
                    help="The last stage is run again, held harder, while the deposited heavy "
                    "atoms end further than this; a run that still ends further is marked "
                    "⚠.")  # fmt: skip
        self.platform = QComboBox()
        self.platform.addItem("Auto", AUTO)
        for name in minimise.platforms():
            self.platform.addItem(name, name)
        found = self.platform.findData(settings.platform)
        self.platform.setCurrentIndex(max(found, 0))
        section.row("Platform", self.platform,
                    help="Where OpenMM runs: Auto takes the fastest that works. The others "
                    "are listed once OpenMM has run (or Settings ▸ External ▸ "
                    "Test).")  # fmt: skip

    def settings(self) -> MinimiseSettings:
        return MinimiseSettings(
            protein=self.forcefield.value(),
            solvent=self.solvent.value(),
            restraints=float(self.restraints.value()),
            rmsd_limit=float(self.rmsd_limit.value()),
            platform=str(self.platform.currentData() or AUTO),
        )

    def allow_waters(self, waters: bool) -> None:
        """Offer only the force fields that can take crystal waters, when there are some."""
        why = MinimiseSettings(protein=spec.FF19SB).problem(has_waters=True)
        self.forcefield.set_enabled(spec.FF19SB, not waters, why if waters else "")

    def set_enabled(self, on: bool) -> None:
        for field in (self.forcefield, self.solvent, self.restraints, self.rmsd_limit,
                      self.platform):  # fmt: skip
            field.setEnabled(on)
