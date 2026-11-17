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
"""Writing superposed structures, shared by Superpose and Similarity Search.

A superposed structure names the reference it was fitted onto (``blob_with_reference``), and
that reference is kept with the document as an asset, so the Browser can draw the two together.
The asset's id is chosen before the work starts (the workers write it into each structure) and
the asset is added in the same undo step as the column (``PutAsset`` in ``extra_ops``).
"""

from __future__ import annotations

import uuid
from collections.abc import Mapping

from rdkit import Chem

from chemsuite.app.document.tools import ToolContext
from chemsuite.core.assets import Asset, PutAsset
from chemsuite.superpose import model as sp

SUFFIX = " superposed"


def output_name(structure: str) -> str:
    """The column a superposed copy of *structure* goes in."""
    return structure if structure.endswith(SUFFIX) else structure + SUFFIX


def new_reference_id(ctx: ToolContext) -> str:
    """An id for a reference asset not yet in the document."""
    return ctx.table.unique_asset_id(f"superpose-{uuid.uuid4().hex[:8]}")


def reference_asset(asset_id: str, name: str, mol: Chem.Mol, meta: Mapping[str, str]) -> Asset:
    """The reference asset: *mol*'s first conformer, shown as *name*."""
    return Asset(asset_id, sp.REFERENCE_KIND, name, sp.reference_block(mol), dict(meta))


def put(asset: Asset) -> PutAsset:
    """The step that adds *asset* to the document (for ``apply_results``' ``extra_ops``)."""
    return PutAsset(asset.id, None, asset)
