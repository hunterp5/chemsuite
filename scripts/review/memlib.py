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

"""Windows process memory counters via ctypes (no psutil), for any PID."""

import ctypes
import threading
import time
from ctypes import wintypes


class PMC(ctypes.Structure):
    _fields_ = [
        ("cb", wintypes.DWORD),
        ("PageFaultCount", wintypes.DWORD),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
        ("PrivateUsage", ctypes.c_size_t),
        ("PrivateWorkingSetSize", ctypes.c_size_t),  # Task Manager's "Memory" column
        ("SharedCommitUsage", ctypes.c_uint64),
    ]


_k = ctypes.windll.kernel32
_psapi = ctypes.windll.psapi
_k.OpenProcess.restype = wintypes.HANDLE
_k.GetCurrentProcess.restype = wintypes.HANDLE
_psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
_k.CloseHandle.argtypes = [wintypes.HANDLE]
MB = 1 << 20


def counters(pid=None):
    """(working set, peak working set, private commit, private working set) MB of *pid*
    (default: this process). The last is Task Manager's default "Memory" column."""
    if pid is None:
        h, close = _k.GetCurrentProcess(), False
    else:
        h, close = _k.OpenProcess(0x1000 | 0x0010, False, pid), True
        if not h:
            return (0, 0, 0, 0)
    c = PMC()
    c.cb = ctypes.sizeof(c)
    _psapi.GetProcessMemoryInfo(h, ctypes.byref(c), c.cb)
    if close:
        _k.CloseHandle(h)
    return (
        round(c.WorkingSetSize / MB),
        round(c.PeakWorkingSetSize / MB),
        round(c.PrivateUsage / MB),
        round(c.PrivateWorkingSetSize / MB),
    )


class Sampler:
    """Samples main + worker private working set MB every *dt* seconds; tracks peaks."""

    def __init__(self, pids=lambda: [], dt=0.2):
        self.pids, self.dt = pids, dt
        self.peak_main = self.peak_total = self.peak_workers = 0
        self._stop = threading.Event()
        self._t = threading.Thread(target=self._run, daemon=True)

    def now(self):
        main = counters()[3]
        workers = [counters(p)[3] for p in self.pids()]
        return main, sum(workers), workers

    def _run(self):
        while not self._stop.is_set():
            main, w, _ = self.now()
            self.peak_main = max(self.peak_main, main)
            self.peak_workers = max(self.peak_workers, w)
            self.peak_total = max(self.peak_total, main + w)
            time.sleep(self.dt)

    def __enter__(self):
        self._t.start()
        return self

    def __exit__(self, *a):
        self._stop.set()
        self._t.join()
