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

"""The compound table core: typed columns, view (filter + sort), undoable edits, events."""

from chemsuite.core.assets import Asset, AssetStore, Packable, PutAsset
from chemsuite.core.history import History

from .columns import NUMBER, TEXT, Column, NumberColumn, TextColumn
from .events import AboutToChange, AssetsChanged, Changed
from .filters import CategoryFilter, Filter, OidFilter, RangeFilter, TextFilter
from .style import ColourRule, ColumnStyle
from .table import SortKey, Table

__all__ = [
    "NUMBER",
    "TEXT",
    "AboutToChange",
    "Asset",
    "AssetStore",
    "AssetsChanged",
    "CategoryFilter",
    "Changed",
    "ColourRule",
    "Column",
    "ColumnStyle",
    "Filter",
    "History",
    "NumberColumn",
    "OidFilter",
    "Packable",
    "PutAsset",
    "RangeFilter",
    "SortKey",
    "Table",
    "TextColumn",
    "TextFilter",
]
