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

"""The chemsuite pKa weights ship with chemsuite, whole, with their licenses."""

from __future__ import annotations

import hashlib

from chemsuite.ionisation.unipka import weights


def test_the_shipped_weights_are_the_released_ones() -> None:
    """A truncated file or a Git LFS pointer would load wrongly or not at all."""
    assert weights.PATH.stat().st_size == weights.SIZE
    digest = hashlib.sha256(weights.PATH.read_bytes()).hexdigest()
    assert digest == weights.SHA256


def test_the_weights_ship_with_their_notice_and_licenses() -> None:
    data = weights.PATH.parent
    notice = (data / "NOTICE-chemsuite-pka.txt").read_text(encoding="utf-8")
    assert "GNU General Public License, version 3" in notice and weights.SHA256 in notice
    assert "MIT License" in (data / "LICENSE-Uni-Mol.txt").read_text(encoding="utf-8")


def test_one_predictor_per_process() -> None:
    assert weights.predictor() is weights.predictor(weights.PATH)
    assert weights.predictor().nets
