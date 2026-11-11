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
"""The Protein Workbench: proteins, their ligands, and how they are drawn.

Everything outside :mod:`.ui` is Qt-free:

- :mod:`.structure`: a protein as a :class:`~.structure.model.Structure`, immutable numpy
  arrays per atom, residue and chain. It is read and written with gemmi (:mod:`.structure.io`),
  bonded from residue templates (:mod:`.structure.bonds`), sorted into polymer, ligand, water
  and the rest (:mod:`.structure.classify`), stored in a document as packed arrays
  (:mod:`.structure.store`), and queried with a PyMOL-like language (:mod:`.structure.select`).
- :mod:`.chem`: ligands as RDKit molecules, poses, contacts and interaction profiles.
- :mod:`.look`: how each protein is shown, and the scene, cartoon and surfaces drawn for it.
- :mod:`.ui`: the Qt workbench (window, view, menus, panels). It never reaches the table: the
  bridge to the table is :mod:`chemsuite.tools.protein`.
"""
