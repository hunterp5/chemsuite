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

"""What a molecule's microstates mean at a given pH, from chemsuite pKa (:mod:`.unipka`).

The model gives every protonation microstate a free energy, so the sites are not treated as
independent: populations at a pH are Boltzmann weights over whole microstates, and the pKa
values are macroscopic, one for each step in net charge. Each step is marked acid or base by
the atom that loses the most hydrogen across it (each macrostate's microstates weighted by
their populations): a base when that atom is cationic in the protonated form, an acid
otherwise.

Microstates that protonate an oxygen or sulfur are left out (:func:`.unipka.ranges.plausible`).

The model is loaded once per process, from the weights shipped with chemsuite
(:mod:`.unipka.weights`); :func:`use` puts another in its place (tests use a stub). The charges
enumerated widen (:mod:`.unipka.ranges`) while a form at an end holds more than
:data:`EDGE` of the molecules at the pH asked for, so the pKa values listed are those of the
charges that matter there (a dication's pKa may appear only at a lower pH). The microstates
start from the range the model was evaluated on (charges −1 to +1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen
from rdkit.Chem.MolStandardize import rdMolStandardize

from .unipka import Ensemble, microstates, ranges, thermo
from .unipka.microstates import Microstate

ACID = "acid"
BASE = "base"
EDGE = 0.01  # a charge at the end of the range with this share at the pH widens it


class FreeEnergies(Protocol):
    """What :func:`profile` needs from the model (a :class:`.unipka.Predictor`)."""

    mean: float

    def free_energy(self, smiles: str) -> float: ...


@dataclass(frozen=True)
class Step:
    """One macroscopic pKa: between the forms of net charge ``charge`` and ``charge − 1``."""

    pka: float
    kind: str  # ACID or BASE
    charge: int


@dataclass(frozen=True)
class Profile:
    """A molecule's neutral parent (see :func:`prepare`), its microstates and their free
    energies, and the pKa steps between them, lowest first."""

    mol: Chem.Mol
    ensemble: Ensemble
    steps: tuple[Step, ...]

    @property
    def most_acidic(self) -> float | None:
        """The lowest acid pKa, or None without an acid step."""
        values = [s.pka for s in self.steps if s.kind == ACID]
        return min(values) if values else None

    @property
    def most_basic(self) -> float | None:
        """The highest base pKa (of the conjugate acid), or None without a base step."""
        values = [s.pka for s in self.steps if s.kind == BASE]
        return max(values) if values else None

    def populations(self, ph: float) -> np.ndarray:
        return self.ensemble.populations(ph)

    def charge(self, ph: float) -> float:
        """Mean net charge at *ph*, counting permanent charges such as quaternary nitrogen."""
        charges = np.array([m.charge for m in self.ensemble.states])
        return float(self.populations(ph) @ charges)

    def fraction_neutral(self, ph: float) -> float:
        """Share at *ph* of the forms with no site ionized: every atom's charge as in the
        neutral parent (zwitterions count as ionized; a permanent charge does not)."""
        parent = tuple(a.GetFormalCharge() for a in self.mol.GetAtoms())
        neutral = np.array([m.charges == parent for m in self.ensemble.states])
        return float(self.populations(ph)[neutral].sum())

    def logd(self, logp: float, ph: float) -> float:
        """*logp* corrected for ionization: only the neutral forms are assumed to partition."""
        fraction = self.fraction_neutral(ph)
        return logp + math.log10(fraction) if fraction > 0.0 else -math.inf

    def major_species(self, ph: float) -> Chem.Mol:
        """The most populated microstate at *ph*, on the parent's atoms. Carbanions are never
        chosen; the parent comes back when no microstate can be built."""
        populations = self.populations(ph)
        order = np.argsort(-populations, kind="stable")
        for i in order.tolist():
            state = self.ensemble.states[i]
            if any(q < 0 and a.GetAtomicNum() == 6
                   for a, q in zip(self.mol.GetAtoms(), state.charges, strict=True)):  # fmt: skip
                continue
            built = with_state(self.mol, state)
            if built is not None:
                return built
        return Chem.Mol(self.mol)


def with_state(mol: Chem.Mol, state: Microstate) -> Chem.Mol | None:
    """*mol* (heavy atoms) with *state*'s hydrogens and charges, or None if RDKit refuses it."""
    out = Chem.RWMol(mol)
    for atom, h, q in zip(out.GetAtoms(), state.hydrogens, state.charges, strict=True):
        atom.SetNumExplicitHs(h)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(q)
    made = out.GetMol()
    try:
        Chem.SanitizeMol(made)
    except (ValueError, RuntimeError):  # RDKit's MolSanitizeException derives from ValueError
        return None
    return made


def prepare(mol: Chem.Mol) -> Chem.Mol | None:
    """The neutral parent the model sees: the largest (organic) fragment, uncharged, heavy
    atoms only; None when RDKit cannot make it."""
    try:
        with rdBase.BlockLogs():
            parent = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
            parent = rdMolStandardize.Uncharger(canonicalOrder=True).uncharge(parent)
            parent = Chem.RemoveHs(parent)
            Chem.SanitizeMol(parent)
    except (ValueError, RuntimeError):
        return None
    return parent if parent.GetNumAtoms() else None


def steps(ensemble: Ensemble) -> tuple[Step, ...]:
    """The macroscopic pKa between each pair of adjacent charges, lowest first, each marked by
    the atom that loses the most hydrogen across it (see the module notes)."""
    g = ensemble.g
    charges = np.array([m.charge for m in ensemble.states])
    hydrogens = np.array([m.hydrogens for m in ensemble.states], dtype=np.float64)
    formal = np.array([m.charges for m in ensemble.states], dtype=np.float64)
    out = []
    for q in sorted(set(charges.tolist()))[1:]:
        acid, base = charges == q, charges == q - 1
        if not base.any():
            continue
        pka = thermo.macro_pka(g[acid], g[base], ensemble.mean)
        w_acid, w_base = _weights(g[acid]), _weights(g[base])
        lost = w_acid @ hydrogens[acid] - w_base @ hydrogens[base]
        site = int(np.argmax(lost))
        kind = BASE if (w_acid @ formal[acid])[site] > 0.5 else ACID
        out.append(Step(pka, kind, q))
    return tuple(sorted(out, key=lambda s: s.pka))


def _weights(g: np.ndarray) -> np.ndarray:
    """Boltzmann weights of free energies *g* (kT) within one macrostate."""
    e = np.exp(-(g - g.min()))
    return e / e.sum()


# -- the model, per process ------------------------------------------------------------------
_model: FreeEnergies | None = None
_last: tuple[Chem.Mol, float, Profile | None] | None = None


def use(model: Path | FreeEnergies | None) -> None:
    """Set the model :func:`profile` runs: the weights at a path, a model object, or None for
    the weights shipped with chemsuite (the default)."""
    global _model, _last
    if isinstance(model, Path):
        from .unipka import weights

        model = weights.predictor(model)
    if model is not _model:
        _model, _last = model, None


def _current() -> FreeEnergies:
    model = _model
    if model is None:
        from .unipka import weights

        model = weights.predictor()
    return model


def profile(mol: Chem.Mol, ph: float = 7.4) -> Profile | None:
    """The pKa profile of *mol* for use at *ph* (which sets the charges enumerated), or None
    when it cannot be prepared or has more microstates than the model is given.

    The last result is remembered, so several descriptors of one molecule (cLogD, CNS MPO…)
    run the model once.
    """
    global _last
    model = _current()
    last = _last  # one read: another thread may replace it
    if last is not None and last[0] is mol and last[1] == ph:
        return last[2]
    parent = prepare(mol)
    result = None
    if parent is not None:
        try:
            ensemble = _ensemble(parent, model, ph)
        except ranges.TooManyMicrostates:
            ensemble = None
        if ensemble is not None and ensemble.states:
            result = Profile(parent, ensemble, steps(ensemble))
    _last = (mol, ph, result)
    return result


def _ensemble(parent: Chem.Mol, model: FreeEnergies, ph: float) -> Ensemble:
    def pools_of(lower: int, upper: int) -> ranges.Pools:
        return ranges.plausible(microstates.ensemble(parent, lower, upper), parent)

    def accept(pools: ranges.Pools) -> list[Ensemble]:
        states = tuple(m for pool in pools.values() for m in pool)
        g = np.array([model.free_energy(m.smiles) for m in states])
        return [Ensemble(states, g, model.mean)] if states else []

    def ends(found: list[Ensemble], low: int, high: int) -> tuple[bool, bool]:
        ensemble = found[0]
        charges = np.array([m.charge for m in ensemble.states])
        populations = ensemble.populations(ph)
        below = populations[charges == low].sum() > EDGE
        above = populations[charges == high].sum() > EDGE
        return bool(below), bool(above)

    found = ranges.widen(pools_of, accept, ends)
    return found[0] if found else Ensemble((), np.zeros(0), model.mean)


def clogd(mol: Chem.Mol, ph: float = 7.4) -> float:
    """cLogD at *ph*: Crippen cLogP of the neutral parent, corrected with chemsuite pKa's
    populations. The parent is the largest fragment, uncharged, so a salt gives its parent's
    value."""
    found = profile(mol, ph)
    if found is None:
        raise ValueError("chemsuite pKa could not prepare the structure")
    return found.logd(float(Crippen.MolLogP(found.mol)), ph)
