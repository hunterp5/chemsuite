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
"""Shared downloads: the IT policy, HTTPS only, size limits (a stub opener: no network)."""

from __future__ import annotations

import io
import json
import ssl
import urllib.error

import pytest

from chemsuite.core import download
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.settings import POLICY_ENV


class Opener:
    def __init__(self, body: bytes = b"data", error: int = 0) -> None:
        self.body = body
        self.error = error
        self.urls: list[str] = []

    def __call__(self, request, timeout=None):
        self.urls.append(request.full_url)
        if self.error:
            raise urllib.error.HTTPError(request.full_url, self.error, "Not Found", {}, None)
        response = io.BytesIO(self.body)
        response.headers = {"Content-Length": str(len(self.body))}
        return response


@pytest.fixture
def blocked(tmp_path, monkeypatch):
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"network_allowed": False}))
    monkeypatch.setenv(POLICY_ENV, str(policy))


def test_reads_a_small_file(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv(POLICY_ENV, str(tmp_path / "none.json"))
    opener = Opener(b"x" * 3_000_000)
    assert download.network_allowed()
    data = download.fetch_bytes("https://example.org/a", purpose="test", opener=opener)
    assert len(data) == 3_000_000 and opener.urls == ["https://example.org/a"]


def test_the_policy_can_forbid_the_network(blocked, tmp_path) -> None:
    opener = Opener()
    assert not download.network_allowed()
    with pytest.raises(download.NetworkBlocked, match="turned off"):
        download.fetch_bytes("https://example.org/a", purpose="test", opener=opener)
    with pytest.raises(download.NetworkBlocked):
        download.download_file("https://example.org/a", tmp_path / "a", purpose="test",
                               opener=opener)  # fmt: skip
    assert opener.urls == []  # nothing was asked of the network


def test_only_https(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv(POLICY_ENV, str(tmp_path / "none.json"))
    with pytest.raises(download.DownloadError, match="not HTTPS"):
        download.fetch_bytes("http://example.org/a", purpose="test", opener=Opener())


def test_errors_and_limits_are_said_plainly(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv(POLICY_ENV, str(tmp_path / "none.json"))
    with pytest.raises(download.DownloadError, match="HTTP 404 Not Found"):
        download.fetch_bytes("https://example.org/a", purpose="test", opener=Opener(error=404))
    with pytest.raises(download.DownloadError, match="larger than"):
        download.fetch_bytes("https://example.org/a", purpose="test",
                             opener=Opener(b"x" * (3 << 20)), max_bytes=2 << 20)  # fmt: skip
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        download.fetch_bytes("https://example.org/a", purpose="test", opener=Opener(),
                             cancel=token)  # fmt: skip


def test_https_has_certificates() -> None:
    context = download.ssl_context()
    assert isinstance(context, ssl.SSLContext) and context.verify_mode == ssl.CERT_REQUIRED
