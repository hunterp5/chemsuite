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

"""The chemsuite pKa weights: fetched once, kept only when they match their checksum."""

from __future__ import annotations

import hashlib
import io

import pytest

from chemsuite.core.download import DownloadError
from chemsuite.ionisation.unipka import weights

PAYLOAD = b"chemsuite pKa weights, for the test" * 100


class _Response(io.BytesIO):
    status = 200

    def __init__(self, data: bytes) -> None:
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data))}


@pytest.fixture
def pinned(monkeypatch):
    monkeypatch.setattr(weights, "SIZE", len(PAYLOAD))
    monkeypatch.setattr(weights, "SHA256", hashlib.sha256(PAYLOAD).hexdigest())


def _opener(data: bytes, calls: list[str]):
    def opener(request, timeout):
        calls.append(request.full_url)
        return _Response(data)

    return opener


def test_fetch_downloads_once_and_then_uses_the_cache(tmp_path, pinned) -> None:
    calls: list[str] = []
    assert weights.cached(tmp_path) is None
    path = weights.fetch(tmp_path, opener=_opener(PAYLOAD, calls))
    assert path.read_bytes() == PAYLOAD and calls == [weights.URL]
    assert weights.cached(tmp_path) == path
    assert weights.fetch(tmp_path, opener=_opener(PAYLOAD, calls)) == path
    assert len(calls) == 1


def test_a_file_that_does_not_match_its_checksum_is_removed(tmp_path, pinned) -> None:
    wrong = PAYLOAD[::-1]
    with pytest.raises(DownloadError, match="checksum"):
        weights.fetch(tmp_path, opener=_opener(wrong, []))
    assert not (tmp_path / weights.NAME).exists()


def test_a_damaged_file_in_the_cache_is_not_used(tmp_path, pinned) -> None:
    (tmp_path / weights.NAME).write_bytes(PAYLOAD[::-1])
    assert weights.cached(tmp_path) is None


def test_the_weights_come_from_chemsuites_release_over_https() -> None:
    assert weights.URL.startswith("https://github.com/hunterp5/chemsuite/releases/download/")
    assert weights.URL.endswith(weights.NAME)
