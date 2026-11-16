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
"""A preparation being reviewed (Tools ▸ Inspect…), kept with the document until applied.

A review is the settings it was surveyed with and what the chemist has said so far: the
issues they unticked and the states they chose by hand. That is all that is saved. The survey
itself — every state and what it costs — is worked out again from the structure when the
document reopens, so the numbers shown are always the structure's own, never a stale copy.

Once applied, the review becomes the protein's preparation report (:mod:`.reports`), which
records what was overruled, and the pending review is dropped.

Plain JSON under ``reviews`` of the protein workspace's state, holding no
:mod:`~chemsuite.macromolecule.prepare` types (as :mod:`.reports` does not).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

#: the neutral N-terminus's state name (``prepare.sites.NEUTRAL_N``)
NEUTRAL_N = "N-neutral"
#: the settings a survey is run with, and their defaults (``pipeline.Options``'s)
SETTINGS = {"ph": 7.0, "waters": "all", "repair": True, "optimise": True, "titrate": "off"}


def _settings(saved: Any) -> dict[str, Any]:
    out = dict(SETTINGS)
    if not isinstance(saved, dict):
        return out
    for key, default in SETTINGS.items():
        value = saved.get(key)
        if isinstance(default, bool):
            if isinstance(value, bool):
                out[key] = value
        elif isinstance(default, float):
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                out[key] = float(value)
        elif isinstance(value, str):
            out[key] = value
    return out


@dataclass(frozen=True)
class Review:
    """What the chemist has decided so far about one protein's preparation."""

    protein: str
    settings: dict[str, Any] = field(default_factory=lambda: dict(SETTINGS))
    #: issue ids left alone (an unticked repair)
    skipped: tuple[str, ...] = ()
    #: site id → the state chosen by hand
    overrides: dict[str, str] = field(default_factory=dict)

    @property
    def changed(self) -> int:
        """How many lines the chemist overruled."""
        return len(self.skipped) + len(self.overrides)

    def to_json(self) -> dict[str, Any]:
        return {
            "protein": self.protein,
            "settings": dict(self.settings),
            "skipped": list(self.skipped),
            "overrides": dict(self.overrides),
        }

    @classmethod
    def from_json(cls, saved: Any) -> Review | None:
        """The review *saved* held, or None when it is not one."""
        if not isinstance(saved, dict) or not isinstance(saved.get("protein"), str):
            return None
        skipped = saved.get("skipped") or []
        overrides = saved.get("overrides") or {}
        if not isinstance(skipped, list) or not isinstance(overrides, dict):
            return None
        return cls(
            saved["protein"],
            _settings(saved.get("settings")),
            tuple(s for s in skipped if isinstance(s, str)),
            {
                k: _renamed(v)
                for k, v in overrides.items()
                if isinstance(k, str) and isinstance(v, str)
            },  # fmt: skip
        )


def _renamed(label: str) -> str:
    """A state *label* as this chemsuite names it: the neutral N-terminus's states were "NH2 …°"
    until 2026-11-15 (``NH2`` is also the amide cap's residue name), now "N-neutral …°"."""
    old = "NH2 "
    return NEUTRAL_N + label[len(old) - 1 :] if label.startswith(old) else label
