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
"""The preparation's settings: their fields in the Prepare and Inspect dialogs, the form the
document keeps them in, and the line a report gives the ones that were changed.

Every field of :class:`~chemsuite.protein.prepare.pipeline.Options` but the fitted weights is
here. The everyday ones go in the dialog's Method section; the rest in More Options, each
defaulting to what the preparation has always done."""

from __future__ import annotations

from collections.abc import Mapping

from PySide6.QtWidgets import QCheckBox, QDoubleSpinBox

from chemsuite.app.widgets.dialogs import Section, Segmented
from chemsuite.protein.prepare import pipeline, sites

WATERS = ((pipeline.ALL, "All"), (pipeline.NEAR_LIGANDS, "Near Ligands"), (pipeline.NONE, "None"))
TITRATE = ((sites.OFF, "Off"), (sites.POCKET, "Pocket"), (sites.ALL, "All"))
ALTLOCS = ((pipeline.HIGHEST, "Most Occupied"), (pipeline.FIRST, "First"))
HIS = ((pipeline.AUTO, "Auto"), *((form, form) for form in pipeline.HIS_FORMS))
#: what a report calls each setting, in the dialog's order
NAMES = {
    "ph": "pH",
    "waters": "Waters",
    "repair": "Repair",
    "optimise": "Optimise",
    "titrate": "Titrate",
    "side_chains": "Side Chains",
    "fill_gaps": "Lone Residues",
    "caps": "Cap Ends",
    "altloc": "Alternate Locations",
    "water_reach": "Water Distance",
    "water_b_max": "Water B-factor Limit",
    "orient_waters": "Orient Waters",
    "ligand_states": "Ligand Forms",
    "his": "His",
    "flips": "Flips",
    "ambiguous": "⚠ Threshold",
}
CHOICES = {"waters": WATERS, "titrate": TITRATE, "altloc": ALTLOCS, "his": HIS}
B_OFF = 0.0  # the B-factor limit's spin box at this value means no limit


def _check(on: bool = True) -> QCheckBox:
    box = QCheckBox()
    box.setChecked(on)
    return box


def _spin(low: float, high: float, step: float, value: float, suffix: str = "") -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(low, high)
    spin.setSingleStep(step)
    spin.setDecimals(1)
    spin.setValue(value)
    spin.setSuffix(suffix)
    return spin


class OptionsFields:
    """The settings' fields, laid out in a dialog's Method and More Options sections."""

    def __init__(self, method: Section, more: Section) -> None:
        d = pipeline.Options()
        self.ph = _spin(0.0, 14.0, 0.5, d.ph)
        method.row("pH", self.ph, help="Which form each acid, base and ligand takes.")
        self.waters = Segmented(WATERS, d.waters)
        method.row("Waters", self.waters,
                   help="Keep every water, only those near a ligand (More Options ▸ Water "
                   "Distance), or none.")  # fmt: skip
        self.repair = method.row("Repair", _check(d.repair),
                                 help="Put the structure in order before hydrogens go on: one "
                                 "location per atom, missing side chains, a residue missing "
                                 "on its own, caps on cut ends (each under More Options). "
                                 "Longer loops are not built.")  # fmt: skip
        self.optimise = method.row("Optimise", _check(d.optimise),
                                   help="Choose each Asn/Gln/His flip, each His tautomer and "
                                   "each hydroxyl and amine rotor by the hydrogen bonds they "
                                   "make.")  # fmt: skip

        self.titrate = more.row("Titrate", Segmented(TITRATE, d.titrate),
                                help="Work out each group's pKa and take the form it mostly "
                                "holds at this pH, instead of the textbook form. Less accurate "
                                "on average than the textbook forms.")  # fmt: skip
        self.side_chains = more.row("Side Chains", _check(d.side_chains),
                                    help="Build missing side-chain atoms from the component "
                                    "templates, turned clear of their neighbours.")  # fmt: skip
        self.fill_gaps = more.row("Lone Residues", _check(d.fill_gaps),
                                  help="Build a residue the file is missing on its own, from "
                                  "its neighbours. Off: the ends either side are "
                                  "capped.")  # fmt: skip
        self.caps = more.row("Cap Ends", _check(d.caps),
                             help="Cap the ends the file cuts short of the sequence with ACE "
                             "and NME, as a force field needs. Off: they are left as the file "
                             "has them.")  # fmt: skip
        self.altloc = more.row("Alternate Locations", Segmented(ALTLOCS, d.altloc),
                               help="Which location each atom keeps when the file models it "
                               "more than once.")  # fmt: skip
        self.water_reach = more.row("Water Distance", _spin(2.0, 12.0, 0.5, d.water_reach, " Å"),
                                    help="With Waters ▸ Near Ligands, the waters kept are those "
                                    "this close to a ligand atom.")  # fmt: skip
        self.water_b_max = _spin(B_OFF, 200.0, 5.0, B_OFF, " Å²")
        self.water_b_max.setSpecialValueText("Off")
        more.row("Water B-factor Limit", self.water_b_max,
                 help="Drop waters whose B-factor is above this: placed with little "
                 "confidence. Waters kept by hand in Inspect stay.")  # fmt: skip
        self.orient_waters = more.row("Orient Waters", _check(d.orient_waters),
                                      help="Turn the waters near ligands to the hydrogen bonds "
                                      "they make, with the rest of the network.")  # fmt: skip
        self.ligand_states = more.row("Ligand Forms", _check(d.ligand_states),
                                      help="Let each ligand take another protomer or tautomer "
                                      "when the pocket favours it. Off: as drawn.")  # fmt: skip
        self.his = more.row("His", Segmented(HIS, d.his),
                            help="Auto chooses each His tautomer by its hydrogen bonds; or "
                            "give every His the same one. A His bound to a metal keeps its "
                            "own.")  # fmt: skip
        self.flips = more.row("Flips", _check(d.flips),
                              help="Let Asn, Gln and His flip (swap the atoms X-ray cannot "
                              "tell apart). Off: as the file has them.")  # fmt: skip
        self.ambiguous = more.row("⚠ Threshold", _spin(0.0, 5.0, 0.1, d.ambiguous, " kcal/mol"),
                                  help="A choice whose next best is closer than this is marked "
                                  "⚠ in the report and in Inspect.")  # fmt: skip
        for box in (self.repair, self.optimise):
            box.toggled.connect(lambda *_: self.follow())
        self.waters.changed.connect(lambda *_: self.follow())
        self.follow()

    def follow(self) -> None:
        """Settings that only mean something under another follow it."""
        repair, optimise = self.repair.isChecked(), self.optimise.isChecked()
        for field in (self.side_chains, self.fill_gaps, self.caps, self.altloc):
            field.setEnabled(repair)
        for field in (self.titrate, self.orient_waters, self.ligand_states, self.his,
                      self.flips, self.ambiguous):  # fmt: skip
            field.setEnabled(optimise)
        self.water_reach.setEnabled(self.waters.value() == pipeline.NEAR_LIGANDS)

    def options(self) -> pipeline.Options:
        b_max = float(self.water_b_max.value())
        return pipeline.Options(
            ph=float(self.ph.value()),
            waters=self.waters.value(),
            repair=self.repair.isChecked(),
            optimise=self.optimise.isChecked(),
            titrate=self.titrate.value(),
            ligand_states=self.ligand_states.isChecked(),
            orient_waters=self.orient_waters.isChecked(),
            side_chains=self.side_chains.isChecked(),
            fill_gaps=self.fill_gaps.isChecked(),
            caps=self.caps.isChecked(),
            altloc=self.altloc.value(),
            water_reach=float(self.water_reach.value()),
            water_b_max=None if b_max == B_OFF else b_max,
            his=self.his.value(),
            flips=self.flips.isChecked(),
            ambiguous=float(self.ambiguous.value()),
        )

    def set_options(self, options: pipeline.Options) -> None:
        self.ph.setValue(options.ph)
        self.waters.set_value(options.waters)
        self.repair.setChecked(options.repair)
        self.optimise.setChecked(options.optimise)
        self.titrate.set_value(options.titrate)
        self.ligand_states.setChecked(options.ligand_states)
        self.orient_waters.setChecked(options.orient_waters)
        self.side_chains.setChecked(options.side_chains)
        self.fill_gaps.setChecked(options.fill_gaps)
        self.caps.setChecked(options.caps)
        self.altloc.set_value(options.altloc)
        self.water_reach.setValue(options.water_reach)
        self.water_b_max.setValue(B_OFF if options.water_b_max is None else options.water_b_max)
        self.his.set_value(options.his)
        self.flips.setChecked(options.flips)
        self.ambiguous.setValue(options.ambiguous)
        self.follow()


# -- as the document keeps them ------------------------------------------------------------------
def settings_of(options: pipeline.Options) -> dict[str, object]:
    """*options* as plain values (every setting but the fitted weights)."""
    return {
        "ph": options.ph,
        "waters": options.waters,
        "repair": options.repair,
        "optimise": options.optimise,
        "titrate": options.titrate,
        "ligand_states": options.ligand_states,
        "orient_waters": options.orient_waters,
        "side_chains": options.side_chains,
        "fill_gaps": options.fill_gaps,
        "caps": options.caps,
        "altloc": options.altloc,
        "water_reach": options.water_reach,
        "water_b_max": options.water_b_max,
        "his": options.his,
        "flips": options.flips,
        "ambiguous": options.ambiguous,
    }


def options_from(settings: Mapping) -> pipeline.Options:
    """The options *settings* holds; a setting it lacks, or holds in a form not understood (an
    older document), takes its default."""
    d = pipeline.Options()

    def choice(name: str, default: str) -> str:
        value = settings.get(name)
        return value if value in {key for key, _ in CHOICES[name]} else default

    def flag(name: str, default: bool) -> bool:
        value = settings.get(name)
        return value if isinstance(value, bool) else default

    def number(name: str, default: float, low: float, high: float) -> float:
        value = settings.get(name)
        ok = isinstance(value, (int, float)) and not isinstance(value, bool)
        return min(max(float(value), low), high) if ok else default

    b_max = number("water_b_max", B_OFF, B_OFF, 200.0)
    return pipeline.Options(
        ph=number("ph", d.ph, 0.0, 14.0),
        waters=choice("waters", d.waters),
        repair=flag("repair", d.repair),
        optimise=flag("optimise", d.optimise),
        titrate=choice("titrate", d.titrate),
        ligand_states=flag("ligand_states", d.ligand_states),
        orient_waters=flag("orient_waters", d.orient_waters),
        side_chains=flag("side_chains", d.side_chains),
        fill_gaps=flag("fill_gaps", d.fill_gaps),
        caps=flag("caps", d.caps),
        altloc=choice("altloc", d.altloc),
        water_reach=number("water_reach", d.water_reach, 2.0, 12.0),
        water_b_max=None if b_max == B_OFF else b_max,
        his=choice("his", d.his),
        flips=flag("flips", d.flips),
        ambiguous=number("ambiguous", d.ambiguous, 0.0, 5.0),
    )


def describe(options: pipeline.Options) -> str:
    """The settings changed from the defaults, as a report line gives them ("" for none)."""
    parts = []
    for name, value in options.changed().items():
        if name not in NAMES or name == "ph":  # the report gives the pH already
            continue
        if isinstance(value, bool):
            shown = "on" if value else "off"
        elif name in CHOICES:
            shown = dict(CHOICES[name]).get(str(value), str(value))
        elif name == "water_b_max":
            shown = "off" if value is None else f"{value:g} Å²"
        elif name == "water_reach":
            shown = f"{value:g} Å"
        elif name == "ambiguous":
            shown = f"{value:g} kcal/mol"
        else:
            shown = f"{value:g}" if isinstance(value, float) else str(value)
        parts.append(f"{NAMES[name]} {shown}")
    return ", ".join(parts)
