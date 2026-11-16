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

"""Simulation: a prepared complex described for a force field, and run in OpenMM's own
interpreter (``docs/SIMULATE_PLAN.md``).

- :mod:`.spec`: what is simulated (:class:`~.spec.SystemSpec`), with which force field and how
  it is minimised, as plain JSON the worker reads without chemsuite;
- :mod:`.build`: a prepared structure → a :class:`~.spec.SystemSpec`, refusing what no force
  field here can describe;
- :mod:`.result`: what a minimisation gives back.

Qt-free; OpenMM and OpenFF are never imported here, only in the worker's interpreter.
"""
