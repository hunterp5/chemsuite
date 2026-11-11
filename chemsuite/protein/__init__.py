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
"""Proteins without Qt: structures, their parts, and what can be asked of them.

A protein is a :class:`~.model.Structure`, immutable numpy arrays per atom, residue and chain
(:mod:`.model`). It is read and written with gemmi (:mod:`.io`), bonded from residue templates
(:mod:`.bonds`), sorted into polymer, ligand, water and the rest (:mod:`.classify`), stored in
a document as packed arrays (:mod:`.store`), and queried with a PyMOL-like language
(:mod:`.select`). Secondary structure (:mod:`.secondary`), the sequence against SEQRES
(:mod:`.sequence`) and ligands as RDKit molecules (:mod:`.ligands`) come from the same arrays.
``chemsuite.app.protein`` draws a structure.
"""
