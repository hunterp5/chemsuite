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

"""What the machine has: physical cores and memory (no psutil; best effort, never raises).

Used to size the worker processes so that background work leaves the rest of the machine
usable: one worker per physical core but one (a second worker on a core's other hardware
thread competes with every other program for that core), and no more workers than a share of
the memory holds.
"""

from __future__ import annotations

import ctypes
import os
import sys
from ctypes import wintypes
from dataclasses import dataclass
from pathlib import Path

_MB = 1 << 20


def logical_cores() -> int:
    return os.cpu_count() or 1


def physical_cores() -> int:
    """Physical cores (hardware threads when they cannot be told apart)."""
    try:
        if sys.platform == "win32":
            found = _windows_cores()
        elif sys.platform == "darwin":
            found = _sysctl("hw.physicalcpu")
        else:
            found = _linux_cores()
    except (OSError, ValueError, AttributeError):
        found = 0
    return found if 0 < found <= logical_cores() else logical_cores()


@dataclass(frozen=True)
class Memory:
    total_mb: int
    available_mb: int | None  # free or reclaimable now (None: not known)


def memory() -> Memory | None:
    """The machine's physical memory, or None when it cannot be read."""
    try:
        if sys.platform == "win32":
            return _windows_memory()
        if sys.platform == "darwin":
            total = _sysctl("hw.memsize")
            return Memory(total // _MB, None) if total > 0 else None
        return _linux_memory()
    except (OSError, ValueError, AttributeError):
        return None


# -- Windows ------------------------------------------------------------------------------
def _windows_cores() -> int:
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    size = wintypes.DWORD(0)
    relation_core = 0
    kernel.GetLogicalProcessorInformationEx(relation_core, None, ctypes.byref(size))
    if size.value == 0:
        return 0
    buffer = ctypes.create_string_buffer(size.value)
    if not kernel.GetLogicalProcessorInformationEx(relation_core, buffer, ctypes.byref(size)):
        return 0
    count, offset = 0, 0
    while offset < size.value:  # one variable-length record per core: Relationship, Size, …
        record = int.from_bytes(buffer.raw[offset + 4 : offset + 8], "little")
        if record <= 0:
            break
        count += 1
        offset += record
    return count


class _MemoryStatus(ctypes.Structure):
    _fields_ = [
        ("dwLength", wintypes.DWORD),
        ("dwMemoryLoad", wintypes.DWORD),
        ("ullTotalPhys", ctypes.c_uint64),
        ("ullAvailPhys", ctypes.c_uint64),
        ("ullTotalPageFile", ctypes.c_uint64),
        ("ullAvailPageFile", ctypes.c_uint64),
        ("ullTotalVirtual", ctypes.c_uint64),
        ("ullAvailVirtual", ctypes.c_uint64),
        ("ullAvailExtendedVirtual", ctypes.c_uint64),
    ]


def _windows_memory() -> Memory | None:
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    status = _MemoryStatus()
    status.dwLength = ctypes.sizeof(status)
    if not kernel.GlobalMemoryStatusEx(ctypes.byref(status)):
        return None
    return Memory(status.ullTotalPhys // _MB, status.ullAvailPhys // _MB)


# -- macOS --------------------------------------------------------------------------------
def _sysctl(name: str) -> int:
    libc = ctypes.CDLL("/usr/lib/libc.dylib")
    value = ctypes.c_uint64(0)
    size = ctypes.c_size_t(ctypes.sizeof(value))
    if libc.sysctlbyname(name.encode(), ctypes.byref(value), ctypes.byref(size), None, 0):
        return 0
    return int(value.value)


# -- Linux --------------------------------------------------------------------------------
def _linux_cores() -> int:
    cores = set()
    physical = core = ""
    for line in Path("/proc/cpuinfo").read_text().splitlines() + [""]:
        key, _, value = line.partition(":")
        key = key.strip()
        if key == "physical id":
            physical = value.strip()
        elif key == "core id":
            core = value.strip()
        elif not line.strip():
            if core:
                cores.add((physical, core))
            physical = core = ""
    return len(cores)


def _linux_memory() -> Memory | None:
    fields = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, _, value = line.partition(":")
        parts = value.split()
        if parts and parts[0].isdigit():
            fields[key] = int(parts[0]) // 1024  # kB
    if "MemTotal" not in fields:
        return None
    return Memory(fields["MemTotal"], fields.get("MemAvailable"))
