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
"""Conformers with CONFORGE, CDPKit's conformer generator (optional: the ``cdpkit`` package).

CONFORGE (Seidel, Permann, Wieder, Kohlbacher and Langer, J. Chem. Inf. Model. 2023, 63,
5549–5570, doi:10.1021/acs.jcim.3c00563) builds conformers from a fragment library and a
torsion library and scores them with MMFF94s. CDPKit is LGPL-2.1 and is imported, not copied.

The molecule goes to CDPKit with its hydrogens explicit, so CDPKit adds no atoms and the atom
order is RDKit's; the conformers come back onto the same RDKit molecule, lowest energy first.
CONFORGE is deterministic: :attr:`ConformerSettings.seed` and ``optimize`` do not apply.
"""

from __future__ import annotations

import importlib.util

from rdkit import Chem, rdBase
from rdkit.Geometry import Point3D

from chemsuite.conformers.etkdg import Conformers, ConformerSettings

FORCE_FIELD = "MMFF94s"
#: an energy window wide enough to keep every conformer (``energy_window=None``)
_NO_WINDOW = 1.0e6
#: conformers sampled for a macrocycle (CONFORGE's stochastic mode). Its default, 2,000, took
#: 40–70 s on plain 13- and 18-membered rings (ETKDG: under 1 s); 200 took 7 s. Provisional
#: until the benchmark (docs/conformers/bench.md) tunes it.
MACROCYCLE_SAMPLES = 200
#: the MOL block format: V2000 unless the molecule is too big for it
_V2000_ATOMS = 999


def available() -> bool:
    """Whether CDPKit (the ``cdpkit`` package, imported as ``CDPL``) is installed."""
    return importlib.util.find_spec("CDPL") is not None


def _messages() -> dict[int, str]:
    import CDPL.ConfGen as ConfGen

    code = ConfGen.ReturnCode
    return {
        code.TIMEOUT: "timed out",
        code.ABORTED: "aborted",
        code.FORCEFIELD_SETUP_FAILED: "no MMFF94 parameters",
        code.FORCEFIELD_MINIMIZATION_FAILED: "force-field minimisation failed",
        code.FRAGMENT_CONF_GEN_FAILED: "fragment conformers failed",
        code.FRAGMENT_CONF_GEN_TIMEOUT: "fragment conformers timed out",
        code.TORSION_DRIVING_FAILED: "torsion driving failed",
        code.CONF_GEN_FAILED: "conformer generation failed",
        code.MAX_ROT_BOND_COUNT_EXCEEDED: "too many rotatable bonds",
    }


def _read(block: str) -> object | None:
    """The CDPKit molecule of the MOL *block*, or None."""
    import CDPL.Base as Base
    import CDPL.Chem as CChem

    target = CChem.BasicMolecule()
    return target if CChem.MOLMoleculeReader(Base.StringIOStream(block)).read(target) else None


def generate(mol: Chem.Mol, settings: ConformerSettings) -> Conformers:
    """Conformers of *mol* (which is not changed). Failures come back with ``mol=None``."""
    import CDPL.Base.Exceptions as CDPLErrors
    import CDPL.ConfGen as ConfGen

    with rdBase.BlockLogs():
        work = Chem.AddHs(Chem.Mol(mol))
        work.RemoveAllConformers()
        block = Chem.MolToMolBlock(work, forceV3000=work.GetNumAtoms() > _V2000_ATOMS)
    try:
        target = _read(block)
        if target is None:
            return Conformers(None, error="CDPKit could not read the structure")
        ConfGen.prepareForConformerGeneration(target)
        if target.numAtoms != work.GetNumAtoms():
            return Conformers(None, error="CDPKit changed the atoms")
        generator = ConfGen.ConformerGenerator()
        options = generator.settings
        options.setMaxNumOutputConformers(settings.count)
        options.setMinRMSD(settings.prune_rms)
        window = settings.energy_window
        options.setEnergyWindow(_NO_WINDOW if window is None else window)
        options.setTimeout(settings.timeout * 1000)
        options.setMaxNumSampledConformers(MACROCYCLE_SAMPLES)
        status = generator.generate(target)
    except CDPLErrors.Exception as err:  # CDPKit's own errors (bad input, I/O)
        return Conformers(None, error=str(err) or type(err).__name__)
    code = ConfGen.ReturnCode
    if status not in (code.SUCCESS, code.TOO_MUCH_SYMMETRY) or generator.numConformers == 0:
        return Conformers(None, error=_messages().get(status, f"CONFORGE failed ({status})"))

    found = [generator.getConformer(i) for i in range(generator.numConformers)]
    found.sort(key=lambda c: c.energy)
    for new_id, coords in enumerate(found):
        conformer = Chem.Conformer(work.GetNumAtoms())
        for atom in range(work.GetNumAtoms()):
            v = coords[atom]
            conformer.SetAtomPosition(atom, Point3D(v[0], v[1], v[2]))
        conformer.Set3D(True)
        conformer.SetId(new_id)
        work.AddConformer(conformer, assignId=False)
    return Conformers(work, [float(c.energy) for c in found], FORCE_FIELD)
