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

"""Sites of metabolism: which atom of a molecule a cytochrome P450 attacks (Qt-free).

An atom-level gradient boosting machine over cheap RDKit descriptors
(:mod:`.features`), evaluated in numpy over exported tree arrays (:mod:`.trees`), with
:mod:`.sites` as the per-molecule interface and :mod:`.predict` the batch driver.
docs/SOM_PLAN.md holds the design, the licence chain and the measured accuracy.
"""

from __future__ import annotations

#: The nine CYP isoforms the model is trained for, in the column order of the one-hot
#: block of :mod:`.features` — which is also ``som-benchmark``'s and ATTNSOM's order.
#: Changing it invalidates every exported model.
ISOFORMS = ("1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4")

#: The four a dialog offers first: the isoforms that clear most marketed drugs.
COMMON_ISOFORMS = ("3A4", "2D6", "2C9", "1A2")
