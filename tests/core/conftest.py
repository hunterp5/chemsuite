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

from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

DRUGS = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    "alanine": "C[C@H](N)C(=O)O",
}


def molblock(smiles: str, title: str = "", three_d: bool = False) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if three_d:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=7)
    else:
        AllChem.Compute2DCoords(mol)
    mol.SetProp("_Name", title)
    return Chem.MolToMolBlock(mol)


def sdf_record(block: str, fields: dict[str, str]) -> str:
    out = block.rstrip("\n") + "\n"
    for name, value in fields.items():
        out += f"> <{name}>\n{value}\n\n"
    return out + "$$$$\n"


@pytest.fixture
def drug_sdf(tmp_path: Path) -> Path:
    """Six drugs with fields, one record with a broken molecule, and one with no title."""
    records = []
    for i, (name, smiles) in enumerate(DRUGS.items()):
        fields = {"ID": f"CPD-{i + 1}", "pIC50": f"{5 + i * 0.25:.2f}", "Name": name}
        if name == "caffeine":
            fields["Note"] = "line one\nline two"
        records.append(sdf_record(molblock(smiles, title=name if i else ""), fields))
    broken = molblock("CC").replace("    0.0000 C ", "    0.0000 Xx", 1)  # unknown element
    records.append(sdf_record(broken, {"ID": "CPD-BAD", "pIC50": "4.00", "Name": "broken"}))
    path = tmp_path / "drugs.sdf"
    path.write_text("".join(records), encoding="utf-8")
    return path
