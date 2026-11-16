"""M0 gate 6 for docs/SIMULATE_PLAN.md: a worker inside WSL, started from Windows.

Run with chemsuite's own interpreter on the Windows PC:  python docs/simulate/m0_wsl.py

It answers: does text survive the pipe, does `children.stop` (a Job Object kill of wsl.exe)
end the Linux process, does the Linux process end when its stdin closes, and how a Windows
path is handed over.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time

from chemsuite.core.jobs import children

#: the Linux side: echoes JSON lines, exits at end of input; MARK makes it findable by pgrep
WORKER = (
    "import sys, json, os\n"
    "MARK = 'cs_m0_worker'\n"
    "print(json.dumps({'ready': True, 'pid': os.getpid()}), flush=True)\n"
    "for line in sys.stdin:\n"
    "    print(json.dumps({'echo': json.loads(line)}), flush=True)\n"
)


def start() -> subprocess.Popen:
    return children.start(
        ["wsl.exe", "-e", "python3", "-I", "-c", WORKER],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        bufsize=1,
    )


def linux_alive(pid: int) -> bool:
    done = subprocess.run(
        ["wsl.exe", "-e", "sh", "-c", f"kill -0 {pid} 2>/dev/null && echo yes || echo no"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    return done.stdout.strip() == "yes"


def wait_gone(pid: int, seconds: float = 6.0) -> float | None:
    """Seconds until the Linux process is gone, or None when it outlives the wait."""
    begun = time.perf_counter()
    while time.perf_counter() - begun < seconds:
        if not linux_alive(pid):
            return round(time.perf_counter() - begun, 2)
        time.sleep(0.25)
    return None


def main() -> int:
    out: dict = {}

    process = start()
    assert process.stdin and process.stdout
    t = time.perf_counter()
    pid = json.loads(process.stdout.readline())["pid"]
    out["start_seconds"] = round(time.perf_counter() - t, 2)
    process.stdin.write(json.dumps({"text": "Å µ → 日本"}) + "\n")
    process.stdin.flush()
    out["echo"] = json.loads(process.stdout.readline())["echo"]["text"]
    out["alive_before_stop"] = linux_alive(pid)
    children.stop(process)
    out["gone_after_job_kill_s"] = wait_gone(pid)
    if out["gone_after_job_kill_s"] is None:
        subprocess.run(["wsl.exe", "-e", "kill", "-9", str(pid)], check=False)

    process = start()
    assert process.stdin and process.stdout
    pid = json.loads(process.stdout.readline())["pid"]
    process.stdin.close()
    out["gone_after_stdin_close_s"] = wait_gone(pid)
    out["wsl_exit_code"] = process.wait(timeout=10)
    children.stop(process)

    done = subprocess.run(
        ["wsl.exe", "-e", "wslpath", "-a", sys.executable],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    out["wslpath"] = {sys.executable: done.stdout.strip()}
    print("GATE", json.dumps({"wsl": out}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
