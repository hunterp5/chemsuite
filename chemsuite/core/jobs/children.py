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

"""External programs chemsuite starts (AutoDock Vina, the AiZynthFinder worker, …).

Every one is started through :func:`start`, so that :func:`stop` ends it together with any
programs it started in turn, and :func:`stop_all` ends them all when chemsuite quits:

- on macOS and Linux each program leads a process group of its own, and stopping it signals
  the whole group (SIGTERM, then SIGKILL after a grace period);
- on Windows each program runs in a Job Object of its own that kills everything in it when
  its last handle closes. :func:`stop` terminates the job; and when chemsuite ends, even by a
  crash or End Task, Windows closes the handles and the programs end with it. Only the
  programs started here are put in a job: files and links chemsuite opens in other
  applications are left alone.
"""

from __future__ import annotations

import contextlib
import logging
import os
import signal
import subprocess
import sys
import threading
from collections.abc import Sequence

log = logging.getLogger(__name__)

#: seconds a program gets to end after being asked, before it is killed
GRACE_S = 2.0

_lock = threading.Lock()
_live: dict[subprocess.Popen, object] = {}  # process -> its Windows job (None elsewhere)


def start(command: Sequence[str], **kwargs) -> subprocess.Popen:
    """``subprocess.Popen(command, **kwargs)``, registered so chemsuite can stop it.

    No console window opens on Windows; elsewhere the program gets a process group of its own.
    """
    if sys.platform == "win32":
        kwargs["creationflags"] = kwargs.get("creationflags", 0) | subprocess.CREATE_NO_WINDOW
    else:
        kwargs.setdefault("start_new_session", True)
    process = subprocess.Popen(list(command), **kwargs)
    job = _windows_job(process.pid) if sys.platform == "win32" else None
    with _lock:
        for ended in [p for p in _live if p.poll() is not None and p is not process]:
            _close_job(_live.pop(ended))
        _live[process] = job
    return process


def running() -> list[subprocess.Popen]:
    """The registered programs still running."""
    with _lock:
        return [p for p in _live if p.poll() is None]


def stop(process: subprocess.Popen, grace: float = GRACE_S) -> None:
    """End *process* and every program it started, then forget it. Never raises."""
    with _lock:
        job = _live.pop(process, None)
    try:
        if sys.platform == "win32":
            _stop_windows(process, job, grace)
        else:
            _stop_group(process, grace)
    except (OSError, subprocess.SubprocessError) as err:
        log.warning("could not stop %s (pid %s): %s", _name(process), process.pid, err)
    finally:
        _close_job(job)


def stop_all(grace: float = GRACE_S) -> None:
    """End every registered program (on the way out)."""
    with _lock:
        processes = list(_live)
    for process in processes:
        stop(process, grace)


def _name(process: subprocess.Popen) -> str:
    args = process.args
    first = args[0] if isinstance(args, list | tuple) and args else args
    return os.path.basename(str(first))


# -- macOS and Linux: process groups ----------------------------------------------------------


def _stop_group(process: subprocess.Popen, grace: float) -> None:
    group = process.pid  # start_new_session: the program leads its own group
    if process.poll() is None:
        _signal_group(group, signal.SIGTERM)
        try:
            process.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            _signal_group(group, signal.SIGKILL)
            process.wait(timeout=grace)
    _signal_group(group, signal.SIGKILL)  # whatever it started and left behind


def _signal_group(group: int, sig: int) -> None:
    with contextlib.suppress(ProcessLookupError, PermissionError):  # gone, or never ours
        os.killpg(group, sig)


# -- Windows: Job Objects ---------------------------------------------------------------------

_KILL_ON_JOB_CLOSE = 0x2000
_EXTENDED_LIMIT_INFORMATION = 9
_PROCESS_SET_QUOTA = 0x0100
_PROCESS_TERMINATE = 0x0001


def _windows_job(pid: int) -> object:
    """A kill-on-close Job Object holding *pid*, or None (logged) when one cannot be made."""
    if sys.platform != "win32":
        return None
    import ctypes
    from ctypes import wintypes

    class Basic(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_int64),
            ("PerJobUserTimeLimit", ctypes.c_int64),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class Counters(ctypes.Structure):
        _fields_ = [(name, ctypes.c_uint64) for name in ("r", "w", "o", "rb", "wb", "ob")]

    class Extended(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", Basic),
            ("IoInfo", Counters),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    kernel = _kernel32()
    job = kernel.CreateJobObjectW(None, None)
    if not job:
        log.warning("could not create a job object: error %s", ctypes.get_last_error())
        return None
    info = Extended()
    info.BasicLimitInformation.LimitFlags = _KILL_ON_JOB_CLOSE
    handle = kernel.OpenProcess(_PROCESS_SET_QUOTA | _PROCESS_TERMINATE, False, pid)
    ok = bool(handle) and bool(
        kernel.SetInformationJobObject(
            job, _EXTENDED_LIMIT_INFORMATION, ctypes.byref(info), ctypes.sizeof(info)
        )
    )
    ok = ok and bool(kernel.AssignProcessToJobObject(job, handle))
    error = ctypes.get_last_error()
    if handle:
        kernel.CloseHandle(handle)
    if not ok:
        kernel.CloseHandle(job)
        log.warning("could not put pid %s in a job object: error %s", pid, error)
        return None
    return job


def _stop_windows(process: subprocess.Popen, job: object, grace: float) -> None:
    if job is not None and sys.platform == "win32":
        _kernel32().TerminateJobObject(job, 1)  # the program and everything it started
    elif process.poll() is None:
        process.kill()
    process.wait(timeout=grace)


def _close_job(job: object) -> None:
    if job is not None and sys.platform == "win32":
        _kernel32().CloseHandle(job)


def _kernel32():
    import ctypes
    from ctypes import wintypes

    if sys.platform != "win32":
        raise OSError("kernel32 exists only on Windows")
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateJobObjectW.restype = wintypes.HANDLE
    kernel.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel.SetInformationJobObject.argtypes = [
        wintypes.HANDLE,
        ctypes.c_int,
        ctypes.c_void_p,
        wintypes.DWORD,
    ]
    kernel.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    kernel.TerminateJobObject.argtypes = [wintypes.HANDLE, wintypes.UINT]
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    return kernel
