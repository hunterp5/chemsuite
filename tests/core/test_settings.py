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

from __future__ import annotations

import json
import sys

import pytest

from chemsuite.core import external_tools as et
from chemsuite.core.settings import (
    DARK,
    GROOVY,
    Policy,
    Settings,
    SettingsStore,
    add_recent,
    custom_theme_id,
    from_dict,
)


def test_defaults_round_trip(tmp_path) -> None:
    store = SettingsStore(tmp_path, Policy())
    assert store.user == Settings()  # no file yet
    settings = store.user.copy()
    settings.appearance.theme = GROOVY
    settings.appearance.custom_themes["Ocean"] = {"window": "#102030", "text": "#ffffff"}
    settings.shortcuts["tool.descriptors"] = ["Ctrl+Alt+D"]
    settings.tool_paths["obabel"] = "/opt/obabel"
    store.save(settings)
    again = SettingsStore(tmp_path, Policy())
    assert again.user == settings
    assert not list(tmp_path.glob(".settings-*"))  # no temporary file left behind


def test_bad_values_fall_back_to_defaults(tmp_path) -> None:
    raw = {
        "appearance": {
            "theme": "neon",
            "app_font_pt": 12,
            "table_font_pt": 99,
            "status_bar": "yes",
        },
        "files": {"recent": ["a", "b", "c"], "recent_limit": 2},
        "shortcuts": {"edit.undo": "Ctrl+Z"},  # must be a list
        "unknown": 1,
    }
    settings, problems = from_dict(raw)
    assert settings.appearance.theme == "system" and settings.appearance.app_font_pt == 12
    assert settings.appearance.table_font_pt == 0 and settings.appearance.status_bar is True
    assert settings.files.recent == ["a", "b"] and settings.shortcuts == {}
    assert set(problems) == {
        "appearance.theme",
        "appearance.table_font_pt",
        "appearance.status_bar",
        "shortcuts",
    }


def test_unreadable_file_is_set_aside(tmp_path) -> None:
    (tmp_path / "settings.json").write_text("{not json")
    store = SettingsStore(tmp_path, Policy())
    assert store.user == Settings()
    assert (tmp_path / "settings.json.bad").exists()


def test_policy_overrides_and_locks(tmp_path) -> None:
    policy_file = tmp_path / "policy.json"
    policy_file.write_text(
        json.dumps(
            {
                "settings": {
                    "appearance.theme": DARK,
                    "performance.max_workers": 2,
                    "tool_paths": {"gnina": "/shared/gnina"},
                    "appearance.app_font_pt": "big",  # refused
                },
                "network_allowed": False,
            }
        )
    )
    policy = Policy.load(policy_file)
    assert policy.problems == ["appearance.app_font_pt"]
    store = SettingsStore(tmp_path / "user", policy)
    user = store.user.copy()
    user.appearance.theme = GROOVY
    user.tool_paths["obabel"] = "/mine/obabel"
    store.save(user)
    effective = store.effective()
    assert effective.appearance.theme == DARK and effective.performance.max_workers == 2
    assert effective.tool_paths == {"obabel": "/mine/obabel", "gnina": "/shared/gnina"}
    assert store.user.appearance.theme == GROOVY  # the user's own choice is kept underneath
    assert store.locked("appearance.theme") and store.locked("tool_paths.obabel")
    assert not store.locked("appearance.table_font_pt")
    assert policy.get("network_allowed") is False


def test_policy_path_from_environment(tmp_path, monkeypatch) -> None:
    from chemsuite.core.settings import machine_policy_path

    monkeypatch.setenv("CHEMSUITE_POLICY", str(tmp_path / "p.json"))
    assert machine_policy_path() == tmp_path / "p.json"
    assert Policy.load().settings == {}  # a missing file is no policy


def test_dotted_access_and_recent() -> None:
    settings = Settings()
    settings.set("appearance.theme", custom_theme_id("Ocean"))
    assert settings.get("appearance.theme") == "custom:Ocean"
    settings.set("tool_paths.java", "/usr/bin/java")
    assert settings.get("tool_paths.java") == "/usr/bin/java"
    for name in ("a", "b", "c", "a"):
        add_recent(settings, name)
    assert settings.files.recent == ["a", "c", "b"]


def test_external_tools(tmp_path) -> None:
    settings = Settings()
    settings.tool_paths["obabel"] = "/configured/obabel"
    assert et.resolve("obabel", settings) == "/configured/obabel"
    assert et.resolve("biotransformer", settings) is None  # a file: nothing to look up
    ok, message = et.test_tool("java", str(tmp_path / "missing"))
    assert not ok and "No such file" in message
    jar = tmp_path / "bt.jar"
    jar.write_text("")
    assert et.test_tool("biotransformer", str(jar)) == (True, "File found.")
    ok, message = et.test_tool("aizynthfinder_python", sys.executable)
    assert ok and "Python" in message
    assert all(t.id in et.BY_ID for t in et.visible_tools())


@pytest.mark.parametrize("theme", ["system", "light", "dark", "groovy", "custom:Mine"])
def test_theme_names(theme) -> None:
    assert from_dict({"appearance": {"theme": theme}})[0].appearance.theme == theme
