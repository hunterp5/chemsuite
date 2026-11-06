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

"""What a molecule's predicted pKa values mean at a given pH.

MolGpKa predicts each site on its own, so the sites are treated as independent
(Henderson-Hasselbalch per site). That is a good approximation when sites are far apart and a
rough one when they interact, as in zwitterions or polyamines.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Crippen

from .molgpka import ACID, BASE, Site, predict_sites, prepare


@dataclass(frozen=True)
class Profile:
    """A molecule as MolGpKa sees it (see :func:`~.molgpka.prepare`), and its sites."""

    mol: Chem.Mol
    sites: tuple[Site, ...]

    @property
    def most_acidic(self) -> float | None:
        """The lowest acid pKa, or None without an acidic site."""
        values = [s.pka for s in self.sites if s.kind == ACID]
        return min(values) if values else None

    @property
    def most_basic(self) -> float | None:
        """The highest pKa of a conjugate acid (the strongest base), or None without a base."""
        values = [s.pka for s in self.sites if s.kind == BASE]
        return max(values) if values else None

    def charge(self, ph: float) -> float:
        """Mean net charge at *ph*, counting permanent charges such as quaternary nitrogen."""
        fixed = sum(atom.GetFormalCharge() for atom in self.mol.GetAtoms())
        return fixed + sum(
            _charged(site, ph) * (1 if site.kind == BASE else -1) for site in self.sites
        )

    def fraction_neutral(self, ph: float) -> float:
        """Fraction of molecules with no site ionised at *ph* (zwitterions count as ionised)."""
        fraction = 1.0
        for site in self.sites:
            fraction *= 1.0 - _charged(site, ph)
        return fraction

    def logd(self, logp: float, ph: float) -> float:
        """*logp* corrected for ionisation: only the neutral form is assumed to partition."""
        fraction = self.fraction_neutral(ph)
        return logp + math.log10(fraction) if fraction > 0.0 else -math.inf

    def major_species(self, ph: float) -> Chem.Mol:
        """The dominant protonation state at *ph*: each site in its majority form."""
        ionised: dict[int, Site] = {}
        for site in self.sites:
            if _charged(site, ph) <= 0.5:
                continue
            other = ionised.get(site.atom)  # both a base and an acid: the more ionised wins
            if other is None or _charged(site, ph) > _charged(other, ph):
                ionised[site.atom] = site
        mol = Chem.RWMol(self.mol)
        for atom_index, site in ionised.items():
            atom = mol.GetAtomWithIdx(atom_index)
            hydrogens = atom.GetTotalNumHs()
            if site.kind == BASE:
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                atom.SetNumExplicitHs(hydrogens + 1)
            else:
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                atom.SetNumExplicitHs(max(0, hydrogens - 1))
            atom.SetNoImplicit(True)
        out = mol.GetMol()
        try:
            Chem.SanitizeMol(out)
        except (ValueError, RuntimeError):  # RDKit's MolSanitizeException derives from ValueError
            return Chem.Mol(self.mol)
        return out


def _charged(site: Site, ph: float) -> float:
    """Fraction of the site in its charged form at *ph*."""
    exponent = (site.pka - ph) if site.kind == BASE else (ph - site.pka)
    exponent = max(-50.0, min(50.0, exponent))
    return 1.0 / (1.0 + 10.0 ** (-exponent))


_last: tuple[Chem.Mol, Profile | None] | None = None


def profile(mol: Chem.Mol) -> Profile | None:
    """The pKa profile of *mol*, or None when it cannot be prepared.

    The last result is remembered, so several descriptors of one molecule (cLogD, CNS MPO…)
    run the model once.
    """
    global _last
    last = _last  # one read: another thread may replace it
    if last is not None and last[0] is mol:
        return last[1]
    prepared = prepare(mol)
    result = Profile(prepared, tuple(predict_sites(prepared))) if prepared is not None else None
    _last = (mol, result)
    return result


def clogd(mol: Chem.Mol, ph: float = 7.4) -> float:
    """cLogD at *ph*: Crippen cLogP of the neutral parent, corrected with MolGpKa's pKa values.

    The parent is the largest fragment, uncharged, so a salt gives its parent's value.
    """
    found = profile(mol)
    if found is None:
        raise ValueError("MolGpKa could not prepare the structure")
    return found.logd(float(Crippen.MolLogP(found.mol)), ph)
