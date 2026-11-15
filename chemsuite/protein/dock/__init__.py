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
"""Docking, without Qt and without the table: the search box, PDBQT files and AutoDock Vina.

- :mod:`.box`: where a ligand may go (from a bound ligand, a residue or a selection).
- :mod:`.pdbqt`: the Meeko boundary — the receptor and each ligand as PDBQT, and poses read
  back as molecules with their bond orders.
- :mod:`.vina`: the bundled AutoDock Vina binary (``bin/``, Apache-2.0), one ligand at a time.
- :mod:`.run`: many ligands, with progress and cancel, keeping what finished.

The receptor must be prepared first (:mod:`chemsuite.protein.prepare`): Vina needs hydrogens
to assign its atom types.
"""
