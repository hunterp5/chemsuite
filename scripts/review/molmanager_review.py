#!/usr/bin/env python3
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

"""Sketcher review: test MolManager's sketcher chemistry against independent references.

    python scripts/review/molmanager_review.py MOLMANAGER_SKETCHER_DIR KETCHER_FG_SDF

MOLMANAGER_SKETCHER_DIR is ``molmanager/ui/sketcher`` of a MolManager checkout and
KETCHER_FG_SDF is Ketcher's ``packages/ketcher-react/src/templates/fg.sdf``. The modules are
loaded read-only for testing; nothing is copied into chemsuite. Results print as Markdown and
are discussed in docs/sketcher/molmanager-review.md.

1. ``alkene_stereo.infer_alkene_ez_for_sketch_mol`` against RDKit's CIP labeller
   (``rdCIPLabeler``, Hanson et al., J. Chem. Inf. Model. 2018) on alkenes of known
   configuration.
2. ``contracted_labels`` against Ketcher's functional-group definitions and against the
   isomer each label's name or alias states.
3. ``valence.target_valence_for_implicit_h`` against the reference species of
   ``valence_review.py``.
"""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdCIPLabeler, rdDepictor

sys.path.insert(0, str(Path(__file__).parent))
from valence_review import REFERENCE, _drawn

# -- loading MolManager modules without their Qt dependencies --------------------------------


def load_sketcher(directory: Path) -> types.ModuleType:
    """A stand-in package for MolManager's sketcher so its modules import on their own."""
    package = types.ModuleType("mm_sketcher")
    package.__path__ = [str(directory)]
    sys.modules["mm_sketcher"] = package
    wildcards = types.ModuleType("mm_sketcher.wildcards")  # imports PyQt5; not under test
    wildcards._is_wildcard_node = lambda node: False
    wildcards._normalize_wildcard_elements = lambda node: [node.get("element", "C")]
    sys.modules["mm_sketcher.wildcards"] = wildcards
    return package


def module(name: str) -> types.ModuleType:
    directory = Path(sys.modules["mm_sketcher"].__path__[0])
    spec = importlib.util.spec_from_file_location(f"mm_sketcher.{name}", directory / f"{name}.py")
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


# -- 1. alkene E/Z ---------------------------------------------------------------------------

#: Alkenes with stereo in the SMILES; RDKit's CIP labeller supplies the reference label.
ALKENES = [
    "F/C=C/Cl",
    "F/C=C\\Cl",
    "C/C=C/C",
    "C/C=C\\C",
    "OC/C=C/C",
    "OC/C=C\\CC",
    "Cl/C=C/Br",
    "BrC(/Cl)=C/F",
    "C/C(CO)=C/C",
    "C/C(CC)=C(/C)CO",
    "CC/C(C)=C(/C)CCO",
    "N#C/C=C/C(=O)O",
    "OC(=O)/C=C\\C(=O)O",
    "OC(=O)/C=C/C(=O)O",
    "C/C=C/c1ccccc1",
    "CC(C)/C=C/C(C)(C)C",
    "C/C(=C\\c1ccccc1)C(=O)O",
    "C\\C(F)=C(\\Cl)C",
    "CC/C(=C(/CC)c1ccccc1)c1ccc(OCCN(C)C)cc1",  # tamoxifen (Z)
    "C/C(=C\\CO)CCC=C(C)C",  # geraniol (E)
    "C/C(=C/CO)CCC=C(C)C",  # nerol (Z)
    "O=C(O)/C=C/c1ccc(O)cc1",  # p-coumaric acid (E)
    "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O",  # linoleic acid (Z,Z)
    "C/C=C/C=C/C(=O)O",  # sorbic acid (E,E)
    "CO/N=C/c1ccccc1",  # oxime ether C=N
    "OC/C(Cl)=C(/Br)C",
    "C/C(I)=C(/Br)Cl",
    "[2H]/C=C/C",
    "C/C(=C(/C)C(C)C)CC",
    "OC(=O)C(/C)=C(\\C)C(C)(C)C",
]


def review_alkenes(alkene_stereo: types.ModuleType) -> list[str]:
    rows = ["| SMILES | Bond | CIP (RDKit) | MolManager |", "|---|---|---|---|"]
    agree = total = 0
    for smiles in ALKENES:
        mol = Chem.MolFromSmiles(smiles)
        rdDepictor.Compute2DCoords(mol)  # keeps the double-bond geometry of the SMILES
        rdCIPLabeler.AssignCIPLabels(mol)
        theirs = alkene_stereo.infer_alkene_ez_for_sketch_mol(Chem.Mol(mol))
        for bond in mol.GetBonds():
            if not bond.HasProp("_CIPCode"):
                continue
            key = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
            want = bond.GetProp("_CIPCode")
            got = theirs.get(key, "none")
            total += 1
            agree += got == want
            mark = "" if got == want else " ✗"
            rows.append(f"| `{smiles}` | {key[0]}={key[1]} | {want} | {got}{mark} |")
    return [f"**{agree}/{total} double bonds agree with CIP.**", "", *rows]


# -- 2. contracted labels --------------------------------------------------------------------

#: Labels whose name or alias states one isomer: label -> SMILES of that group on [*].
NAMED_ISOMERS = {
    "Naph": ("1-NAPH", "*c1cccc2ccccc12"),
    "Tol": ("P-TOL", "*c1ccc(C)cc1"),
    "Mes": ("MESITYL", "*c1c(C)cc(C)cc1C"),
    "Ad": ("ADAMANTYL (1-)", "*C12CC3CC(CC(C3)C1)C2"),
    "sBu": ("SBU", "*C(C)CC"),
    "iBu": ("IBU", "*CC(C)C"),
    "PMB": ("4-MEOC6H4CH2", "*Cc1ccc(OC)cc1"),
}


def _on_star(smiles: str) -> str | None:
    """Canonical SMILES of a group given with its attachment atom first, bonded to [*]."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    star = rw.AddAtom(Chem.Atom(0))
    rw.AddBond(0, star, Chem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(rw)
    except (ValueError, RuntimeError):
        return None
    return Chem.MolToSmiles(rw)


def ketcher_groups(path: Path) -> dict[str, str]:
    """Ketcher's functional groups: name (upper case) -> canonical SMILES on [*]."""
    out: dict[str, str] = {}
    for block in path.read_text(encoding="utf-8").split("$$$$\n"):
        block = block.lstrip("\n")  # Ketcher's files put a blank line after each $$$$
        if "M  END" not in block:
            continue
        name = block.splitlines()[0].strip()
        molblock = block[: block.index("M  END") + len("M  END")]
        mol = Chem.MolFromMolBlock(molblock)
        sap = [line for line in molblock.splitlines() if line.startswith("M  SAP")]
        if mol is None or not sap:
            continue
        attach = int(sap[0].split()[4]) - 1
        rw = Chem.RWMol(mol)
        star = rw.AddAtom(Chem.Atom(0))
        rw.AddBond(attach, star, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except (ValueError, RuntimeError):
            continue
        out[name.upper()] = Chem.MolToSmiles(rw)
    return out


def review_labels(labels: types.ModuleType, fg_sdf: Path) -> list[str]:
    ketcher = ketcher_groups(fg_sdf)
    unique = {lab.display: lab for lab in labels._CONTRACTED.values()}
    aliases: dict[str, list[str]] = {}
    for key, lab in labels._CONTRACTED.items():
        aliases.setdefault(lab.display, []).append(key)
    rows = ["| Label | MolManager on [*] | Ketcher on [*] | Verdict |", "|---|---|---|---|"]
    counts = {"agree": 0, "differ": 0, "not in Ketcher": 0, "invalid": 0}
    for display, lab in sorted(unique.items()):
        mine = _on_star(lab.smiles)
        theirs = ketcher.get(display.upper())
        if mine is None:
            verdict = "invalid SMILES"
            counts["invalid"] += 1
        elif theirs is None:
            verdict = "not in Ketcher"
            counts["not in Ketcher"] += 1
        elif theirs == mine:
            verdict = "agree"
            counts["agree"] += 1
        else:
            verdict = "**differ**"
            counts["differ"] += 1
        rows.append(f"| {display} | `{mine}` | `{theirs or ''}` | {verdict} |")
    named = ["| Label | Name or alias says | Expected | MolManager | |", "|---|---|---|---|---|"]
    for display, (alias, expected) in NAMED_ISOMERS.items():
        lab = unique.get(display)
        if lab is None:
            continue
        mine = _on_star(lab.smiles)
        want = Chem.MolToSmiles(Chem.MolFromSmiles(expected))
        named.append(
            f"| {display} | {alias} | `{want}` | `{mine}` | {'ok' if mine == want else '✗'} |"
        )
    shared: dict[str, set[str]] = {}
    for key, lab in labels._CONTRACTED.items():
        shared.setdefault(key, set()).add(lab.display)
    summary = ", ".join(f"{k}: {v}" for k, v in counts.items())
    return [
        f"**{len(unique)} labels: {summary}.**",
        "",
        *rows,
        "",
        "Labels whose name states an isomer:",
        "",
        *named,
    ]


# -- 3. valence ------------------------------------------------------------------------------


def review_valence(valence: types.ModuleType) -> list[str]:
    rows = ["| Species | Expected H | MolManager |", "|---|---|---|"]
    ok = 0
    for smiles, index, want, species in REFERENCE:
        atom = _drawn(smiles).GetAtomWithIdx(index)
        bonds = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
        target = valence.target_valence_for_implicit_h(
            atom.GetSymbol(), bonds, atom.GetFormalCharge()
        )
        got = target - bonds
        ok += got == want
        if got != want:
            rows.append(f"| {species} (`{smiles}`) | {want} | {got} ✗ |")
    return [f"**{ok}/{len(REFERENCE)} reference species right.** Failures:", "", *rows]


def main(argv: list[str]) -> int:
    RDLogger.DisableLog("rdApp.*")
    load_sketcher(Path(argv[1]))
    out = sys.stdout
    out.write("## 1. alkene_stereo vs CIP\n\n")
    out.write("\n".join(review_alkenes(module("alkene_stereo"))) + "\n\n")
    out.write("## 2. contracted_labels vs Ketcher fg.sdf\n\n")
    out.write("\n".join(review_labels(module("contracted_labels"), Path(argv[2]))) + "\n\n")
    out.write("## 3. valence.py vs reference species\n\n")
    module("rgroups")
    out.write("\n".join(review_valence(module("valence"))) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
