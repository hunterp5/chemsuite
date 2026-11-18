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

"""The ligand pKa model Prepare and Inspect run with."""

from __future__ import annotations

import importlib

from chemsuite.ionisation.unipka import weights
from chemsuite.macromolecule.prepare.ligand_pka import ModelStates
from chemsuite.macromolecule.prepare.pipeline import Options


def test_ligand_forms_take_the_shipped_model() -> None:
    pka_model = importlib.import_module("chemsuite.tools.pka.model")  # not the conftest's stand-in
    importlib.reload(pka_model)
    model = pka_model.ligand_pka(Options())
    assert isinstance(model, ModelStates) and model.predictor is weights.predictor()
    assert pka_model.ligand_pka(Options(ligand_states=False)) is None
