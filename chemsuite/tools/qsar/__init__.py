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
"""Data ▸ QSAR: fit a model on the table's data, apply a saved one, refine it with new rows.

The bridge between the table and :mod:`chemsuite.qsar`, which knows nothing about either Qt or
the document. This package reads the columns, hands the engine a feature pool, writes the
predictions back as columns, and keeps the fitted model with the document (:mod:`.part`).
"""
