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

from chemsuite.core.wording import plural


def test_plural_counts_with_separators_and_real_plurals():
    assert plural(1, "row") == "1 row"
    assert plural(0, "row") == "0 rows"
    assert plural(1234, "conformer") == "1,234 conformers"
    assert plural(2, "matched pair") == "2 matched pairs"
    assert plural(3, "analysis", "analyses") == "3 analyses"
