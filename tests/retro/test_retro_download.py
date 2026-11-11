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

"""Downloading AiZynthFinder's public data (RETRO_PLAN A4), with a stub opener: no network."""

from __future__ import annotations

import io
import json
import urllib.error

import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.retro import download
from chemsuite.retro.config import catalogue, config_problem


class _Response(io.BytesIO):
    def __init__(self, body: bytes, *, status: int = 200, length: int | None = None) -> None:
        super().__init__(body)
        self.status = status
        self.headers = {"Content-Length": str(len(body) if length is None else length)}


class Opener:
    """Serves each file's name as its body, and records every request."""

    def __init__(self, *, honour_range: bool = True, fail: str = "", short: str = "") -> None:
        self.requests: list = []
        self.honour_range = honour_range
        self.fail = fail
        self.short = short

    def body(self, url: str) -> bytes:
        name = next(item.name for item in download.FILES if item.url == url)
        return (name * 50).encode()

    def __call__(self, request, timeout=None):
        assert timeout
        self.requests.append(request)
        url = request.full_url
        if self.fail and self.fail in url:
            raise urllib.error.HTTPError(url, 503, "Service Unavailable", {}, None)
        body = self.body(url)
        if self.short and self.short in url:
            return _Response(body[:10], length=len(body))
        wanted = request.get_header("Range")
        if wanted and self.honour_range:
            start = int(wanted.split("=")[1].rstrip("-"))
            return _Response(body[start:], status=206)
        return _Response(body)


def test_every_file_is_fetched_and_the_config_written(tmp_path) -> None:
    opener = Opener()
    seen = []
    path = download.download_public_data(
        tmp_path / "data", opener=opener, progress=lambda d, t: seen.append((d, t))
    )
    folder = (tmp_path / "data").resolve()
    assert path == folder / "config.yml"
    for item in download.FILES:
        assert (folder / item.name).read_bytes() == opener.body(item.url)
        assert not (folder / (item.name + ".part")).exists()
    assert all(r.full_url.startswith("https://") for r in opener.requests)
    assert seen  # progress was reported
    body = path.read_text(encoding="utf-8")
    assert config_problem(body) == ""
    found = catalogue(path)
    assert found.policies == ("uspto", "ringbreaker")
    assert found.filters == ("uspto",)
    assert found.stocks == ("zinc",)
    assert json.dumps(str(folder / "zinc_stock.hdf5")) in body  # backslashes escaped


def test_files_already_there_are_not_fetched_again(tmp_path) -> None:
    for item in download.FILES[:-1]:
        (tmp_path / item.name).write_bytes(b"done")
    opener = Opener()
    download.download_public_data(tmp_path, opener=opener)
    assert [r.full_url for r in opener.requests] == [download.FILES[-1].url]
    assert (tmp_path / download.FILES[0].name).read_bytes() == b"done"


def test_a_partial_file_resumes_with_a_range_request(tmp_path) -> None:
    item = download.FILES[0]
    opener = Opener()
    full = opener.body(item.url)
    (tmp_path / (item.name + ".part")).write_bytes(full[:100])
    download.download_public_data(tmp_path, opener=opener)
    first = opener.requests[0]
    assert first.get_header("Range") == "bytes=100-"
    assert (tmp_path / item.name).read_bytes() == full


def test_a_server_that_ignores_the_range_restarts_the_file(tmp_path) -> None:
    item = download.FILES[0]
    opener = Opener(honour_range=False)
    full = opener.body(item.url)
    (tmp_path / (item.name + ".part")).write_bytes(b"stale bytes")
    download.download_public_data(tmp_path, opener=opener)
    assert (tmp_path / item.name).read_bytes() == full


def test_cancel_keeps_the_part_file(tmp_path) -> None:
    token = CancelToken()

    def progress(_done, _total) -> None:
        token.cancel()  # after the first chunk of the first file

    with pytest.raises(Cancelled):
        download.download_public_data(tmp_path, opener=Opener(), cancel=token, progress=progress)
    first = download.FILES[0].name
    assert (tmp_path / (first + ".part")).is_file()
    assert not (tmp_path / first).exists()
    assert not (tmp_path / "config.yml").exists()


def test_a_short_download_is_not_kept_as_complete(tmp_path) -> None:
    with pytest.raises(download.DownloadError, match="stopped after 10"):
        download.download_public_data(tmp_path, opener=Opener(short="uspto_model"))
    assert not (tmp_path / "uspto_model.onnx").exists()


def test_an_http_error_says_which_file(tmp_path) -> None:
    with pytest.raises(download.DownloadError, match=r"zinc_stock.hdf5: HTTP 503"):
        download.download_public_data(tmp_path, opener=Opener(fail="figshare"))


def test_a_users_own_config_is_left_alone(tmp_path) -> None:
    own = tmp_path / "config.yml"
    own.write_text("expansion:\n  mine: x.onnx\n")
    path = download.write_config(tmp_path)
    assert path.name == "chemsuite_config.yml"
    assert own.read_text() == "expansion:\n  mine: x.onnx\n"
    # chemsuite's own config is rewritten in place
    (tmp_path / "b").mkdir()
    first = download.write_config(tmp_path / "b")
    assert download.write_config(tmp_path / "b") == first
