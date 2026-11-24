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
"""The models a document keeps, whether or not the panel is open.

A fitted model belongs to the document, not to a window: a chemist fits a model, closes the
panel, saves, comes back next week and applies it to new rows. So the model is a **document
part** (``docs/ARCHITECTURE.md``) whose own state is small — a list of asset ids and the
panel's layout — while each model's bytes live in the document's **asset store**, the same
place proteins do.

Assets are the right home for two reasons beyond size. Adding a model is an op, so it undoes
and redoes with everything else, and the store already knows how to save a payload with the
document and hand back the decoded object — :class:`chemsuite.qsar.model.Model` satisfies
``Packable``, so a model is unpacked once per session and not once per prediction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from chemsuite.core.assets import Asset, PutAsset
from chemsuite.core.history import Op
from chemsuite.qsar.model import Model, unpack

#: The part's key in the document, and the kind of asset a model is stored as.
PART_KEY = "qsar_models"
MODEL_KIND = "qsar_model"


def _asset_name(target: str) -> str:
    """A plain asset id from a target column's name, which may hold anything a chemist types.

    An asset id becomes a name inside the saved document, so it keeps to plain ASCII letters,
    digits and dashes rather than carrying a space, a slash or a micro sign into the archive.
    """
    stem = "".join(c if (c.isascii() and c.isalnum()) else "-" for c in target.lower()).strip("-")
    while "--" in stem:
        stem = stem.replace("--", "-")
    return f"qsar-{stem}" if stem else "qsar"


@dataclass(frozen=True)
class Kept:
    """One model the document holds: where it is, and what to call it in a list."""

    asset: str
    title: str
    target: str
    task: str

    def as_json(self) -> dict[str, str]:
        return {"asset": self.asset, "title": self.title, "target": self.target, "task": self.task}

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Kept:
        return cls(
            asset=str(data["asset"]),
            title=str(data.get("title", "")),
            target=str(data.get("target", "")),
            task=str(data.get("task", "regression")),
        )

    @classmethod
    def of(cls, asset: str, model: Model) -> Kept:
        return cls(asset, model.title, model.recipe.target, model.recipe.task)


class QsarModels:
    """The document's fitted models, and the panel's state.

    The part never holds a model's bytes itself: :meth:`put` returns the op that adds the
    asset, so the caller can write its result columns and store the model in one undo step.
    """

    def __init__(self, host: Any, saved: Any = None) -> None:
        #: the table (an ``AssetHost``), which owns the asset store
        self.host = host
        self._kept: list[Kept] = []
        self._panel: dict[str, Any] = {}
        self._revision = 0
        self._decoded: dict[str, Model] = {}
        if isinstance(saved, dict):
            for entry in saved.get("models", ()):
                if isinstance(entry, dict) and "asset" in entry:
                    self._kept.append(Kept.from_json(entry))
            panel = saved.get("panel")
            if isinstance(panel, dict):
                self._panel = dict(panel)

    # -- the part --------------------------------------------------------------------------
    @property
    def revision(self) -> int:
        return self._revision

    def capture(self) -> dict[str, Any]:
        """What the document saves: the asset ids, not the models."""
        if not self._kept and not self._panel:
            return {}
        return {"models": [k.as_json() for k in self.kept()], "panel": dict(self._panel)}

    def close(self) -> None:
        self._decoded.clear()

    # -- the models ------------------------------------------------------------------------
    def kept(self) -> list[Kept]:
        """The models still in the document, newest last. An entry whose asset has gone —
        undone, or the row it came with deleted — is dropped rather than reported broken."""
        assets = self.host.assets
        return [k for k in self._kept if k.asset in assets]

    def model(self, asset: str) -> Model:
        """The model stored under *asset*, unpacked once per session."""
        kept = self._decoded.get(asset)
        if kept is not None:
            return kept
        payload = self.host.assets[asset].payload
        model = payload if isinstance(payload, Model) else unpack(payload)
        self._decoded[asset] = model
        return model

    def put(self, model: Model, *, name: str = "") -> tuple[str, Op]:
        """An asset id for *model* and the op that stores it.

        The caller passes the op to ``ctx.apply_results(..., extra_ops=[op])`` so the columns
        and the model land as one undoable step.
        """
        asset_id = self.host.unique_asset_id(_asset_name(model.recipe.target))
        title = name or model.title
        asset = Asset(
            asset_id,
            MODEL_KIND,
            title,
            model,  # a Model is Packable, so the store keeps it decoded
            {
                "target": model.recipe.target,
                "task": model.recipe.task,
                "candidate": model.candidate,
            },
        )
        self._kept.append(Kept.of(asset_id, model))
        self._decoded[asset_id] = model
        self._revision += 1
        return asset_id, PutAsset(asset_id, None, asset)

    def forget(self, asset: str) -> Op | None:
        """Take a model out of the document. Returns the op that removes its asset."""
        existing = self.host.assets.get(asset)
        self._kept = [k for k in self._kept if k.asset != asset]
        self._decoded.pop(asset, None)
        self._revision += 1
        return None if existing is None else PutAsset(asset, existing, None)

    def rename(self, asset: str, title: str) -> None:
        self._kept = [
            Kept(k.asset, title, k.target, k.task) if k.asset == asset else k for k in self._kept
        ]
        self._revision += 1

    # -- the panel -------------------------------------------------------------------------
    def panel_state(self) -> dict[str, Any]:
        return dict(self._panel)

    def set_panel_state(self, state: dict[str, Any]) -> None:
        if state != self._panel:
            self._panel = dict(state)
            self._revision += 1
