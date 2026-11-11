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
"""Every download chemsuite makes: public files over HTTPS, nothing about the user's data sent.

- The machine policy can forbid the network (``"network_allowed": false``, see
  :mod:`chemsuite.core.settings`): then nothing is fetched and :class:`NetworkBlocked` says so.
- Only HTTPS. Each request is written to the audit log with its purpose.
- :func:`download_file` writes ``<name>.part`` and renames it only when complete, so a finished
  file is never a truncated one; an interrupted ``.part`` resumes with an HTTP range request
  where the server allows it; a cancel stops within one chunk and keeps the ``.part``.
- :func:`fetch_bytes` reads a small file into memory (a PDB entry), with a size limit.
- Certificates: Python's own store, or ``certifi``'s when installed; where Python has none (the
  python.org macOS build before its "Install Certificates" step) the system's
  ``/etc/ssl/cert.pem``, so HTTPS works without the user fixing their Python.

The *opener* is injectable so tests run without a network.
"""

from __future__ import annotations

import contextlib
import logging
import os
import ssl
import sys
import urllib.request
from collections.abc import Callable
from functools import cache
from pathlib import Path
from typing import Any

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.logbook import audit

log = logging.getLogger(__name__)

CHUNK = 1 << 20  # 1 MiB
TIMEOUT_S = 60.0
MAX_BYTES = 256 << 20  # fetch_bytes refuses more than this
SYSTEM_CERTS = Path("/etc/ssl/cert.pem")

Opener = Callable[..., Any]  # urllib.request.urlopen(request, timeout=…) -> response


class DownloadError(RuntimeError):
    """A file could not be fetched or written."""


class NetworkBlocked(DownloadError):
    """The machine policy does not allow chemsuite to use the network."""

    def __init__(self) -> None:
        super().__init__("Downloads are turned off on this computer (IT policy).")


def network_allowed() -> bool:
    """False when the machine policy forbids the network."""
    from chemsuite.core.settings import Policy

    return Policy.load().get("network_allowed", True) is not False


@cache
def ssl_context() -> ssl.SSLContext:
    try:
        import certifi
    except ImportError:
        certifi = None
    if certifi is not None:
        return ssl.create_default_context(cafile=certifi.where())
    context = ssl.create_default_context()
    if sys.platform == "darwin" and not context.get_ca_certs() and SYSTEM_CERTS.is_file():
        context.load_verify_locations(cafile=str(SYSTEM_CERTS))
    return context


def _open(request: urllib.request.Request, timeout: float) -> Any:
    return urllib.request.urlopen(request, timeout=timeout, context=ssl_context())


def _start(url: str, purpose: str) -> urllib.request.Request:
    if not network_allowed():
        raise NetworkBlocked
    if not url.startswith("https://"):
        raise DownloadError(f"Refusing a download that is not HTTPS: {url}")
    audit("network", url=url, purpose=purpose)
    return urllib.request.Request(url, headers={"User-Agent": "chemsuite"})


def fetch_bytes(
    url: str,
    *,
    purpose: str,
    opener: Opener | None = None,
    cancel: CancelToken | None = None,
    max_bytes: int = MAX_BYTES,
) -> bytes:
    """The whole file at *url*, read into memory."""
    request = _start(url, purpose)
    log.info("downloading %s", url)
    try:
        response = (opener or _open)(request, timeout=TIMEOUT_S)
    except OSError as err:  # URLError and HTTPError are OSErrors
        raise DownloadError(f"Could not download {url}: {reason(err)}") from err
    chunks: list[bytes] = []
    got = 0
    with contextlib.closing(response):
        try:
            while chunk := response.read(CHUNK):
                if cancel is not None and cancel.cancelled:
                    raise Cancelled
                got += len(chunk)
                if got > max_bytes:
                    raise DownloadError(f"{url} is larger than {max_bytes >> 20} MB.")
                chunks.append(chunk)
        except OSError as err:
            raise DownloadError(f"Could not download {url}: {reason(err)}") from err
    return b"".join(chunks)


def download_file(
    url: str,
    target: Path,
    *,
    purpose: str,
    opener: Opener | None = None,
    cancel: CancelToken | None = None,
    report: Callable[[int], None] = lambda _bytes: None,
) -> int:
    """Download *url* to *target* (through ``.part``) and return the bytes it has. *report*
    is told the bytes so far."""
    request = _start(url, purpose)
    part = target.with_name(target.name + ".part")
    have = part.stat().st_size if part.is_file() else 0
    if have:
        request.add_header("Range", f"bytes={have}-")
    log.info("downloading %s%s", url, f" from byte {have}" if have else "")
    try:
        response = (opener or _open)(request, timeout=TIMEOUT_S)
    except OSError as err:
        raise DownloadError(f"Could not download {target.name}: {reason(err)}") from err
    with contextlib.closing(response):
        status = int(getattr(response, "status", 200) or 200)
        if have and status != 206:
            have = 0  # the server ignored the range: start the file again
        length = _length(response)
        expected = have + length if length is not None else None
        try:
            with open(part, "ab" if have else "wb") as out:
                got = have
                while True:
                    if cancel is not None and cancel.cancelled:
                        raise Cancelled
                    chunk = response.read(CHUNK)
                    if not chunk:
                        break
                    out.write(chunk)
                    got += len(chunk)
                    report(got)
        except OSError as err:
            raise DownloadError(f"Could not download {target.name}: {reason(err)}") from err
    if expected is not None and got != expected:
        raise DownloadError(
            f"{target.name} stopped after {got:,} of {expected:,} bytes. Download again to resume."
        )
    if got == 0:
        raise DownloadError(f"{target.name} came back empty.")
    os.replace(part, target)
    return got


def _length(response: Any) -> int | None:
    headers = getattr(response, "headers", None)
    raw = headers.get("Content-Length") if headers is not None else None
    try:
        return int(raw) if raw is not None else None
    except (TypeError, ValueError):
        return None


def reason(err: OSError) -> str:
    """What went wrong, in words: ``HTTP 404 Not Found`` or the network's reason."""
    why = getattr(err, "reason", None)
    code = getattr(err, "code", None)
    if code is not None:
        return f"HTTP {code} {why or ''}".strip()
    return str(why or err)
