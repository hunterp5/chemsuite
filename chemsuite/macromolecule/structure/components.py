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
"""Chemical components: a residue's atoms and bonds as the PDB's Chemical Component Dictionary
(CCD) defines them.

They come from the structure file's own ``_chem_comp_atom``/``_chem_comp_bond`` tables (RCSB
mmCIF files carry them for every component in the entry), or from the CCD itself, one small
file per component (:class:`CcdCache`), fetched only when the network is allowed and kept on
disk so each is downloaded once.

A component is kept as plain JSON-able lists, the shape of the bundled ``residues.json``: atoms
``[name, element, charge]`` and bonds ``[atom, atom, order, aromatic]`` with the Kekulé order.
Hydrogens are listed only to work out formal charges when the table has none, and then dropped.
"""

from __future__ import annotations

import json
import urllib.error
from collections.abc import Mapping
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import gemmi
from rdkit import Chem, rdBase

from chemsuite.core.download import DownloadError, Opener, fetch_bytes
from chemsuite.core.jobs import CancelToken

URL = "https://files.rcsb.org/ligands/download/{name}.cif"
META_KEY = "components"  # Structure.meta: the file's components for its ligands, as JSON
_ORDERS = {"SING": 1, "DOUB": 2, "TRIP": 3, "QUAD": 4}
#: the usual valence of elements whose formal charge follows from their bonds
_VALENCE = {"N": 3, "O": 2, "S": 2, "P": 3, "B": 3}
_HYPERVALENT = {"S": (4, 6), "P": (5,)}


@dataclass(frozen=True)
class Component:
    name: str
    atoms: tuple[tuple[str, str, int], ...]  # heavy atoms: (name, element, formal charge)
    bonds: tuple[tuple[str, str, int, bool], ...]  # between heavy atoms: (a, b, order, aromatic)
    leaving: frozenset[str] = frozenset()  # atoms lost when it bonds to another residue

    def as_json(self) -> dict:
        return {"atoms": [list(a) for a in self.atoms], "bonds": [list(b) for b in self.bonds],
                "leaving": sorted(self.leaving)}  # fmt: skip

    @classmethod
    def from_json(cls, name: str, entry: Mapping) -> Component:
        return cls(
            name,
            tuple((str(n), str(e), int(c)) for n, e, c in entry["atoms"]),
            tuple((str(a), str(b), int(o), bool(r)) for a, b, o, r in entry["bonds"]),
            frozenset(map(str, entry.get("leaving", ()))),
        )

    @property
    def template_bonds(self) -> list[tuple[str, str, int]]:
        """The bonds as :func:`.bonds.perceive` takes them (aromatic ones as AROMATIC)."""
        from .model import AROMATIC

        return [(a, b, AROMATIC if aromatic else order) for a, b, order, aromatic in self.bonds]

    def to_mol(self) -> Chem.Mol | None:
        """The component as a sanitised RDKit molecule (heavy atoms, hydrogens implicit), or
        None when RDKit cannot make sense of it."""
        mol = Chem.RWMol()
        index = {}
        table = Chem.GetPeriodicTable()
        for name, element, charge in self.atoms:
            try:
                number = table.GetAtomicNumber(element.capitalize())
            except RuntimeError:  # the CCD's "X" for an unknown atom, say
                return None
            atom = Chem.Atom(number)
            atom.SetFormalCharge(charge)
            index[name] = mol.AddAtom(atom)
        for a, b, order, _aromatic in self.bonds:
            if a in index and b in index:
                mol.AddBond(index[a], index[b], bond_type(order))
        try:
            with rdBase.BlockLogs():
                Chem.SanitizeMol(mol)
        except (Chem.MolSanitizeException, RuntimeError):
            return None
        return mol.GetMol()


def bond_type(order: int) -> Chem.BondType:
    return {2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}.get(order, Chem.BondType.SINGLE)


def from_block(block: gemmi.cif.Block, names: set[str] | None = None) -> dict[str, Component]:
    """The components a CIF block defines (only those in ``names``, when given)."""
    atoms = block.find("_chem_comp_atom.", ["comp_id", "atom_id", "type_symbol", "?charge",
                                            "?pdbx_leaving_atom_flag"])  # fmt: skip
    bonds = block.find("_chem_comp_bond.", ["comp_id", "atom_id_1", "atom_id_2", "value_order",
                                            "?pdbx_aromatic_flag"])  # fmt: skip
    found: dict[str, list[tuple[str, str, int | None]]] = {}
    leaving: dict[str, set[str]] = {}
    for row in atoms:
        name = row.str(0)
        if names is None or name in names:
            charge = _int(row.str(3)) if row.has(3) else None
            found.setdefault(name, []).append((row.str(1), row.str(2).upper(), charge))
            if row.has(4) and row.str(4).upper() == "Y":
                leaving.setdefault(name, set()).add(row.str(1))
    links: dict[str, list[tuple[str, str, int, bool]]] = {}
    for row in bonds:
        name = row.str(0)
        if name in found:
            order = _ORDERS.get(row.str(3).upper(), 1)
            aromatic = row.has(4) and row.str(4).upper() == "Y"
            links.setdefault(name, []).append((row.str(1), row.str(2), order, aromatic))
    return {name: _component(name, found[name], links.get(name, []), leaving.get(name, set()))
            for name in found}  # fmt: skip


def _int(text: str) -> int | None:
    try:
        return int(text)
    except ValueError:
        return None


def _component(
    name: str,
    atoms: list[tuple[str, str, int | None]],
    bonds: list[tuple[str, str, int, bool]],
    leaving: set[str],
) -> Component:
    element = {atom: symbol for atom, symbol, _ in atoms}
    valence = dict.fromkeys(element, 0)
    for a, b, order, _ in bonds:
        if a in valence and b in valence:
            valence[a] += order
            valence[b] += order
    heavy = []
    for atom, symbol, charge in atoms:
        if symbol in ("H", "D"):
            continue
        if charge is None:
            charge = _charge(symbol, valence[atom])
        heavy.append((atom, symbol, charge))
    kept = {atom for atom, _, _ in heavy}
    return Component(name, tuple(heavy), tuple(b for b in bonds if b[0] in kept and b[1] in kept),
                     frozenset(leaving & kept))  # fmt: skip


def _charge(symbol: str, valence: int) -> int:
    """The formal charge of an atom with ``valence`` bonds (hydrogens counted): the CCD's own
    files give it, an entry's tables often do not."""
    usual = _VALENCE.get(symbol)
    if usual is None or valence == 0 or valence in _HYPERVALENT.get(symbol, ()):
        return 0
    if valence == usual + 1:
        return 1
    if valence == usual - 1:
        return -1
    return 0


def to_meta(components: Mapping[str, Component]) -> str:
    return json.dumps({name: c.as_json() for name, c in sorted(components.items())})


@lru_cache(maxsize=32)
def from_meta(text: str) -> dict[str, Component]:
    """The components :func:`to_meta` wrote ({} for anything else)."""
    try:
        entries = json.loads(text) if text else {}
        return {name: Component.from_json(name, e) for name, e in entries.items()}
    except (ValueError, KeyError, TypeError):
        return {}


class CcdCache:
    """Components from the CCD, one ``<NAME>.cif`` each in ``folder``. A component the CCD does
    not have is remembered as ``<NAME>.missing`` and not asked for again; any other failure
    (no network, downloads turned off, a server error) just gives None. ``offline`` reads only
    what is in the folder already (for work on the GUI thread)."""

    def __init__(
        self, folder: Path | str, opener: Opener | None = None, offline: bool = False
    ) -> None:
        self.folder = Path(folder)
        self.opener = opener
        self.offline = offline

    def get(self, name: str, cancel: CancelToken | None = None) -> Component | None:
        block = self.block(name, cancel)
        return None if block is None else from_block(block, {name.upper()}).get(name.upper())

    def block(self, name: str, cancel: CancelToken | None = None) -> gemmi.cif.Block | None:
        """The whole CCD entry for *name* (with its hydrogens and ideal coordinates)."""
        name = name.upper()
        if not name.isalnum() or len(name) > 5:
            return None
        path = self.folder / f"{name}.cif"
        if not path.is_file() and not self._fetched(name, path, cancel):
            return None
        try:
            return gemmi.cif.read(str(path)).sole_block()
        except (RuntimeError, ValueError, OSError):
            return None

    def _fetched(self, name: str, path: Path, cancel: CancelToken | None) -> bool:
        if self.offline or (self.folder / f"{name}.missing").exists():
            return False
        try:
            data = fetch_bytes(URL.format(name=name), purpose=f"Chemical component {name}",
                               opener=self.opener, cancel=cancel)  # fmt: skip
        except DownloadError as error:
            if isinstance(error.__cause__, urllib.error.HTTPError) and error.__cause__.code == 404:
                self._write(self.folder / f"{name}.missing", b"")
            return False
        return self._write(path, data)

    def _write(self, path: Path, data: bytes) -> bool:
        part = path.with_name(path.name + ".part")
        try:
            self.folder.mkdir(parents=True, exist_ok=True)
            part.write_bytes(data)
            part.replace(path)
        except OSError:
            return False
        return True
