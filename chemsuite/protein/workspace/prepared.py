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
"""What a preparation leaves behind beyond the structure, kept with the document.

A minimisation, and later MD and MM-GBSA, need more than the prepared structure. They need
the form Prepare gave each titratable residue (``HID``, ``CYM``…), and which heavy atoms
Prepare built rather than read from the file. Those atoms are never restrained to where they
were built, since nothing deposited says where they belong. :class:`PreparedRecord` keeps
both, with the pKa values and populations of a titrated run. A minimised protein's record
also holds its :class:`Minimisation`.

A record knows the structure it describes by a digest of its atoms (:func:`atoms_digest`). An
edit after the preparation (an atom deleted, a residue added) leaves it describing another
structure: :meth:`PreparedRecord.matches` says so. It is checked when used, so nothing has to
listen for edits.

Residues are named by their labels (``A/HIS 57``), not their indices, as a review's decisions
are. Built atoms are named by uid, which an edit never reuses.

Saved as plain JSON under ``prepared`` of the protein workspace's state. This module never
imports :mod:`chemsuite.simulate`, which depends on the protein packages, so the minimised
system and its result are kept as the JSON they arrive as.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..structure.model import Structure


def atoms_digest(s: Structure) -> str:
    """A fingerprint of *s*'s atoms (their uids and elements, in order). An edit that adds,
    removes or changes an atom's element changes it; moving atoms does not."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(s.uid, dtype=np.int64).tobytes())
    h.update(np.ascontiguousarray(s.element, dtype=np.uint8).tobytes())
    return h.hexdigest()[:32]


@dataclass(frozen=True)
class Minimisation:
    """How a minimised protein was made: from which prepared protein, with what, and how it
    came out. *system* is the :class:`~chemsuite.simulate.spec.SystemSpec`, ligand charges
    included, as its JSON; *result* the measurements, as the report's numbers."""

    source: str
    forcefield: dict[str, Any] = field(default_factory=dict)
    system: dict[str, Any] = field(default_factory=dict)
    result: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {"source": self.source, "forcefield": self.forcefield, "system": self.system,
                "result": self.result}  # fmt: skip

    @classmethod
    def from_json(cls, saved: Any) -> Minimisation | None:
        if not isinstance(saved, dict) or not isinstance(saved.get("source"), str):
            return None
        parts = [saved.get(k) for k in ("forcefield", "system", "result")]
        if not all(isinstance(p, dict) for p in parts):
            return None
        return cls(saved["source"], *parts)


@dataclass(frozen=True)
class PreparedRecord:
    protein: str
    atoms: str  # :func:`atoms_digest` of the structure it describes
    ph: float = 7.0
    settings: Mapping[str, Any] = field(default_factory=dict)  # the preparation's settings
    states: Mapping[str, str] = field(default_factory=dict)  # residue label → form
    pka: Mapping[str, float] = field(default_factory=dict)  # residue label → pKa
    population: Mapping[str, float] = field(default_factory=dict)  # label → protonated share
    built: frozenset[int] = frozenset()  # uids of the heavy atoms Prepare built
    minimised: Minimisation | None = None

    @classmethod
    def of(
        cls,
        protein: str,
        s: Structure,
        *,
        ph: float = 7.0,
        settings: Mapping[str, Any] | None = None,
        states: Mapping[int, str] | None = None,
        pka: Mapping[int, float] | None = None,
        population: Mapping[int, float] | None = None,
        built: Iterable[int] = (),
    ) -> PreparedRecord:
        """The record of *s* as prepared, its per-residue values given by residue index."""

        def named(values: Mapping[int, Any] | None) -> dict[str, Any]:
            return {s.residue_label(int(r)): v for r, v in (values or {}).items()}

        return cls(protein, atoms_digest(s), float(ph), dict(settings or {}), named(states),
                   named(pka), named(population), frozenset(int(u) for u in built))  # fmt: skip

    def matches(self, s: Structure) -> bool:
        """Whether this record still describes *s* (it has not been edited since)."""
        return self.atoms == atoms_digest(s)

    def states_on(self, s: Structure) -> dict[int, str]:
        """The forms, by residue index in *s* (a label *s* no longer has is passed over)."""
        where = {s.residue_label(r): r for r in range(s.n_residues)}
        return {where[label]: form for label, form in self.states.items() if label in where}

    def built_in(self, s: Structure) -> frozenset[int]:
        """The built atoms *s* still has."""
        return self.built & frozenset(s.uid.tolist())

    def to_json(self) -> dict[str, Any]:
        saved: dict[str, Any] = {
            "protein": self.protein,
            "atoms": self.atoms,
            "ph": self.ph,
            "settings": dict(self.settings),
            "states": dict(self.states),
            "pka": dict(self.pka),
            "population": dict(self.population),
            "built": sorted(self.built),
        }
        if self.minimised is not None:
            saved["minimised"] = self.minimised.to_json()
        return saved

    @classmethod
    def from_json(cls, saved: Any) -> PreparedRecord | None:
        """The record *saved* held, or None when it is not one."""
        if not isinstance(saved, dict):
            return None
        protein, atoms, ph = saved.get("protein"), saved.get("atoms"), saved.get("ph", 7.0)
        if not isinstance(protein, str) or not isinstance(atoms, str):
            return None
        if not isinstance(ph, (int, float)) or isinstance(ph, bool):
            return None
        settings = saved.get("settings", {})
        states = _strings(saved.get("states", {}), str)
        pka = _strings(saved.get("pka", {}), float)
        population = _strings(saved.get("population", {}), float)
        built = saved.get("built", [])
        if not isinstance(settings, dict) or states is None or pka is None or population is None:
            return None
        if not isinstance(built, list) or not all(_is_int(u) for u in built):
            return None
        minimised = None
        if "minimised" in saved:
            minimised = Minimisation.from_json(saved["minimised"])
            if minimised is None:
                return None
        return cls(protein, atoms, float(ph), settings, states, pka, population,
                   frozenset(built), minimised)  # fmt: skip


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _strings(saved: Any, kind: type) -> dict[str, Any] | None:
    """*saved* as a mapping of residue labels to values of *kind*, or None when it is not."""
    if not isinstance(saved, dict):
        return None
    out: dict[str, Any] = {}
    for label, value in saved.items():
        if not isinstance(label, str) or isinstance(value, bool):
            return None
        if kind is str and isinstance(value, str):
            out[label] = value
        elif kind is float and isinstance(value, (int, float)):
            out[label] = float(value)
        else:
            return None
    return out


def carried_built(
    new: Iterable[int], before: PreparedRecord | None, s: Structure
) -> frozenset[int]:
    """The atoms Prepare built in *s*: those it built this time (*new*), and those an earlier
    preparation of the same protein built (*before*) that *s* still has.

    A second preparation finds the first one's caps and side chains already there, so it
    would not call them built on its own, and minimisation would restrain them where they
    were guessed."""
    old = before.built_in(s) if before is not None else frozenset()
    return frozenset(int(u) for u in new) | old
