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

"""User settings and the machine policy that IT can use to override them.

Settings are one JSON file in the user's config folder, written atomically. Loading never
fails: unknown keys are ignored, a bad value falls back to its default (and is logged), and an
unreadable file is set aside as ``settings.json.bad`` so the user starts from defaults.

The **policy file** is JSON that IT places on the machine (or points ``CHEMSUITE_POLICY`` at)::

    {"settings": {"appearance.theme": "light", "tool_paths.obabel": "C:/Tools/obabel.exe"},
     "network_allowed": false}

Every key under ``settings`` replaces the user's value and is *locked* (shown greyed with a
padlock in the Settings window). A key naming a whole group (``"tool_paths"``) locks the group.
"""

from __future__ import annotations

import contextlib
import copy
import json
import logging
import os
import sys
import tempfile
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

SCHEMA_VERSION = 1
SETTINGS_FILE = "settings.json"
POLICY_ENV = "CHEMSUITE_POLICY"

SYSTEM, LIGHT, DARK, GROOVY = "system", "light", "dark", "groovy"
BUILTIN_THEMES = (SYSTEM, LIGHT, DARK, GROOVY)
CUSTOM_PREFIX = "custom:"
ALIGN_H = ("left", "center", "right")
ALIGN_V = ("top", "center", "bottom")
MIN_FONT_PT, MAX_FONT_PT = 6, 32
MIN_ROW_HEIGHT, MAX_ROW_HEIGHT = 40, 400
DEFAULT_ROW_HEIGHT = 110
RECENT_LIMIT = 10


def custom_theme_id(name: str) -> str:
    return CUSTOM_PREFIX + name


def custom_theme_name(theme: str) -> str | None:
    return theme[len(CUSTOM_PREFIX) :] if theme.startswith(CUSTOM_PREFIX) else None


@dataclass
class Appearance:
    theme: str = SYSTEM
    custom_themes: dict[str, dict[str, str]] = field(default_factory=dict)  # name → role → #hex
    app_font_pt: int = 0  # 0: the system's size
    table_font_pt: int = 0  # 0: same as the application
    table_align_h: str = "left"  # text cells (numbers stay right-aligned)
    table_align_v: str = "center"
    structure_row_height: int = DEFAULT_ROW_HEIGHT
    status_bar: bool = True


@dataclass
class Files:
    recent: list[str] = field(default_factory=list)
    recent_limit: int = RECENT_LIMIT
    last_directory: str = ""


@dataclass
class Performance:
    max_workers: int = 0  # 0: all cores but one


@dataclass
class Jobs:
    notify: bool = True  # a desktop notification when a long background job ends unseen


@dataclass
class Settings:
    version: int = SCHEMA_VERSION
    appearance: Appearance = field(default_factory=Appearance)
    shortcuts: dict[str, list[str]] = field(default_factory=dict)  # command id → key sequences
    tool_paths: dict[str, str] = field(default_factory=dict)  # external tool id → path
    files: Files = field(default_factory=Files)
    performance: Performance = field(default_factory=Performance)
    jobs: Jobs = field(default_factory=Jobs)

    def copy(self) -> Settings:
        return copy.deepcopy(self)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    # -- dotted access (policy keys, the Settings window) --------------------------------------
    def get(self, key: str) -> Any:
        target: Any = self
        for part in key.split("."):
            target = target[part] if isinstance(target, dict) else getattr(target, part)
        return target

    def set(self, key: str, value: Any) -> None:
        *parents, last = key.split(".")
        target: Any = self
        for part in parents:
            target = target[part] if isinstance(target, dict) else getattr(target, part)
        if isinstance(target, dict):
            target[last] = value
        else:
            setattr(target, last, value)


# -- validation ---------------------------------------------------------------------------------
def _font(v: Any) -> bool:
    return (
        isinstance(v, int)
        and not isinstance(v, bool)
        and (v == 0 or MIN_FONT_PT <= v <= MAX_FONT_PT)
    )


def _theme(v: Any) -> bool:
    return isinstance(v, str) and (
        v in BUILTIN_THEMES or (v.startswith(CUSTOM_PREFIX) and len(v) > 7)
    )


def _colour(v: Any) -> bool:
    return (
        isinstance(v, str)
        and len(v) == 7
        and v[0] == "#"
        and all(c in "0123456789abcdefABCDEF" for c in v[1:])
    )


def _str_list(v: Any) -> bool:
    return isinstance(v, list) and all(isinstance(s, str) for s in v)


def _int_between(lo: int, hi: int):
    return lambda v: isinstance(v, int) and not isinstance(v, bool) and lo <= v <= hi


CHECKS = {
    "appearance.theme": _theme,
    "appearance.custom_themes": lambda v: (
        isinstance(v, dict)
        and all(
            isinstance(k, str)
            and k
            and isinstance(c, dict)
            and all(isinstance(r, str) and _colour(x) for r, x in c.items())
            for k, c in v.items()
        )
    ),
    "appearance.app_font_pt": _font,
    "appearance.table_font_pt": _font,
    "appearance.table_align_h": lambda v: v in ALIGN_H,
    "appearance.table_align_v": lambda v: v in ALIGN_V,
    "appearance.structure_row_height": _int_between(MIN_ROW_HEIGHT, MAX_ROW_HEIGHT),
    "appearance.status_bar": lambda v: isinstance(v, bool),
    "shortcuts": lambda v: (
        isinstance(v, dict) and all(isinstance(k, str) and _str_list(s) for k, s in v.items())
    ),
    "tool_paths": lambda v: (
        isinstance(v, dict) and all(isinstance(k, str) and isinstance(p, str) for k, p in v.items())
    ),
    "files.recent": _str_list,
    "files.recent_limit": _int_between(0, 50),
    "files.last_directory": lambda v: isinstance(v, str),
    "performance.max_workers": _int_between(0, 1024),
    "jobs.notify": lambda v: isinstance(v, bool),
}


def valid(key: str, value: Any) -> bool:
    """Whether *value* is acceptable for *key* (a dotted name, or an entry of a dict group)."""
    check = CHECKS.get(key)
    if check is not None:
        return bool(check(value))
    group, _, _entry = key.partition(".")
    if group == "tool_paths":
        return isinstance(value, str)
    if group == "shortcuts":
        return _str_list(value)
    return False


def _fill(target: Any, raw: dict, prefix: str, problems: list[str]) -> None:
    for f in fields(target):
        if f.name not in raw:
            continue
        key = f"{prefix}{f.name}"
        value = raw[f.name]
        current = getattr(target, f.name)
        if is_dataclass(current):
            if isinstance(value, dict):
                _fill(current, value, key + ".", problems)
            else:
                problems.append(key)
        elif key == "version":
            continue
        elif valid(key, value):
            setattr(target, f.name, copy.deepcopy(value))
        else:
            problems.append(key)


def from_dict(raw: Any) -> tuple[Settings, list[str]]:
    """Settings from parsed JSON, and the keys whose values were refused."""
    settings = Settings()
    problems: list[str] = []
    if not isinstance(raw, dict):
        return settings, ["(the file is not a JSON object)"]
    _fill(settings, raw, "", problems)
    settings.files.recent = settings.files.recent[: max(settings.files.recent_limit, 0)]
    return settings, problems


# -- the machine policy -------------------------------------------------------------------------
def machine_policy_path() -> Path:
    override = os.environ.get(POLICY_ENV, "").strip()
    if override:
        return Path(override)
    if sys.platform == "win32":
        return Path(os.environ.get("PROGRAMDATA", r"C:\ProgramData")) / "chemsuite" / "policy.json"
    if sys.platform == "darwin":
        return Path("/Library/Application Support/chemsuite/policy.json")
    return Path("/etc/chemsuite/policy.json")


@dataclass
class Policy:
    settings: dict[str, Any] = field(default_factory=dict)  # dotted key → forced value
    values: dict[str, Any] = field(default_factory=dict)  # everything else in the file
    path: Path | None = None
    problems: list[str] = field(default_factory=list)

    @classmethod
    def load(cls, path: Path | None = None) -> Policy:
        path = path or machine_policy_path()
        if not path.is_file():
            return cls(path=path)
        try:
            raw = json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError) as err:
            log.error("ignoring unreadable policy file %s: %s", path, err)
            return cls(path=path, problems=[f"unreadable: {err}"])
        if not isinstance(raw, dict):
            return cls(path=path, problems=["not a JSON object"])
        forced = raw.get("settings") or {}
        policy = cls(path=path, values={k: v for k, v in raw.items() if k != "settings"})
        if not isinstance(forced, dict):
            policy.problems.append("settings: not an object")
            return policy
        for key, value in forced.items():
            if key in ("tool_paths", "shortcuts") and isinstance(value, dict):
                for entry, v in value.items():  # a whole group: every entry, and lock it all
                    if valid(f"{key}.{entry}", v):
                        policy.settings[f"{key}.{entry}"] = v
                policy.settings[key] = None
            elif valid(key, value):
                policy.settings[key] = value
            else:
                policy.problems.append(key)
                log.error("policy file %s: refusing %s = %r", path, key, value)
        return policy

    def locked(self, key: str) -> bool:
        if key in self.settings:
            return True
        group = key.partition(".")[0]
        return group in self.settings and self.settings[group] is None

    def apply(self, settings: Settings) -> Settings:
        """A copy of *settings* with the policy's values in place."""
        out = settings.copy()
        for key, value in self.settings.items():
            if value is not None:
                out.set(key, copy.deepcopy(value))
        return out

    def get(self, name: str, default: Any = None) -> Any:
        return self.values.get(name, default)


# -- storage ------------------------------------------------------------------------------------
class SettingsStore:
    """Reads and writes the user's settings file; :meth:`effective` adds the policy."""

    def __init__(self, folder: Path, policy: Policy | None = None) -> None:
        self.path = Path(folder) / SETTINGS_FILE
        self.policy = policy if policy is not None else Policy.load()
        self.user = self.load()

    def load(self) -> Settings:
        if not self.path.is_file():
            return Settings()
        try:
            raw = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as err:
            log.error("settings file %s is unreadable (%s); starting from defaults", self.path, err)
            with contextlib.suppress(OSError):
                os.replace(self.path, self.path.with_name(SETTINGS_FILE + ".bad"))
            return Settings()
        settings, problems = from_dict(raw)
        for key in problems:
            log.warning("settings file %s: ignoring a bad value for %s", self.path, key)
        return settings

    def save(self, settings: Settings) -> None:
        """Write atomically: a crash mid-write leaves the previous file intact."""
        self.user = settings.copy()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(settings.to_dict(), indent=2, sort_keys=True)
        handle, temp = tempfile.mkstemp(prefix=".settings-", dir=self.path.parent)
        try:
            with os.fdopen(handle, "w", encoding="utf-8") as out:
                out.write(text)
            os.replace(temp, self.path)
        except OSError:
            Path(temp).unlink(missing_ok=True)
            raise

    def effective(self) -> Settings:
        return self.policy.apply(self.user)

    def locked(self, key: str) -> bool:
        return self.policy.locked(key)


def add_recent(settings: Settings, path: str) -> None:
    """Put *path* first in the recent-files list (once), trimmed to the limit."""
    recent = [p for p in settings.files.recent if os.path.normcase(p) != os.path.normcase(path)]
    settings.files.recent = [path, *recent][: settings.files.recent_limit]
