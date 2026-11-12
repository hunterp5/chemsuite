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
"""Settings window, custom themes, shortcuts, Open Recent and the machine policy."""

from __future__ import annotations

import json

from PySide6.QtGui import QKeySequence, QPalette
from PySide6.QtWidgets import QApplication

from chemsuite.app.shell.settings_hub import settings_hub
from tests.gui import gui

pytestmark = gui


def _base() -> str:
    return QApplication.instance().palette().color(QPalette.ColorRole.Base).name()


def test_theme_previews_then_cancel_restores(window) -> None:
    before = _base()
    dialog = window.open_settings()
    page = dialog.appearance
    page.theme.setCurrentIndex(page.theme.findData("dark"))
    assert _base() != before  # live preview
    assert settings_hub().settings.appearance.theme == "dark"
    dialog.reject()
    assert _base() == before
    assert settings_hub().settings.appearance.theme == "system"


def test_ok_saves_and_survives_reload(window, tmp_path) -> None:
    dialog = window.open_settings()
    page = dialog.appearance
    page.theme.setCurrentIndex(page.theme.findData("dark"))
    page.row_height.setValue(150)
    page.status_bar.setChecked(False)
    assert not window.statusBar().isVisible()
    dialog.accept()
    saved = json.loads((tmp_path / "config" / "settings.json").read_text())
    assert saved["appearance"]["theme"] == "dark"
    assert saved["appearance"]["structure_row_height"] == 150
    assert window.view._structure_row_height == 150


def test_groovy_reroll(window) -> None:
    dialog = window.open_settings()
    page = dialog.appearance
    page.theme.setCurrentIndex(page.theme.findData("groovy"))
    assert not page.reroll.isHidden()
    first = _base()
    for _ in range(5):  # random: a repeat is possible but not five times running
        page.reroll.click()
        if _base() != first:
            break
    assert _base() != first
    page.theme.setCurrentIndex(page.theme.findData("light"))
    assert page.reroll.isHidden()
    dialog.reject()


def test_custom_theme_editor(qtbot, monkeypatch) -> None:
    from chemsuite.app.shell.custom_theme_dialog import CustomThemeDialog

    editor = CustomThemeDialog({}, None)
    qtbot.addWidget(editor)
    assert not editor.save_button.isEnabled()  # no name yet
    editor.set_colour("base", "#123456")
    assert editor.preview.palette().color(QPalette.ColorRole.Base).name() == "#123456"
    editor.name.setText("  Night   Shift ")
    editor.save()
    assert editor.chosen == "Night Shift"
    assert editor.themes["Night Shift"]["base"] == "#123456"

    editor = CustomThemeDialog({"A": {"base": "#000000"}, "B": {}}, "A", None)
    qtbot.addWidget(editor)
    assert editor.name.text() == "A"
    assert editor.colours["base"] == "#000000"
    monkeypatch.setattr(editor, "confirm_delete", lambda name: True)
    editor.delete_current()
    assert set(editor.themes) == {"B"}
    editor.close_editor()
    assert editor.result() == editor.DialogCode.Accepted  # a deletion counts as a change
    assert editor.chosen is None


def test_custom_theme_applies_from_settings(window, monkeypatch) -> None:
    from chemsuite.app.shell import custom_theme_dialog

    class Fake:
        themes = {"Mine": {"base": "#204060"}}
        chosen = "Mine"

        def __init__(self, *args) -> None:
            pass

        def exec(self):
            return custom_theme_dialog.QDialog.DialogCode.Accepted

    monkeypatch.setattr(custom_theme_dialog, "CustomThemeDialog", Fake)
    dialog = window.open_settings()
    dialog.appearance.edit_custom()
    assert dialog.appearance.theme.currentData() == "custom:Mine"
    assert _base() == "#204060"
    dialog.accept()
    assert settings_hub().settings.appearance.custom_themes == {"Mine": {"base": "#204060"}}


def test_shortcut_change_and_conflict(window) -> None:
    dialog = window.open_settings()
    keys = dialog.keyboard
    keys.select("data.add_column")
    keys.add_recorded(QKeySequence("Ctrl+Shift+K"))
    assert window.command_action("data.add_column").shortcut() == QKeySequence("Ctrl+Shift+K")
    # The same keys on a second command: flagged, and OK refuses until it is fixed.
    keys.select("data.add_rows")
    keys.add_recorded(QKeySequence("Ctrl+Shift+K"))
    problem = dialog.problem()
    assert "Add Column" in problem and "Add Rows" in problem, problem
    assert "Add Column" in keys.status.text()
    dialog.accept()
    assert dialog.isVisible() and dialog.pages.currentRow() == 1
    keys.set_shortcuts("data.add_rows", None)  # back to its default
    assert dialog.problem() == ""
    dialog.accept()
    assert settings_hub().settings.shortcuts == {"data.add_column": ["Ctrl+Shift+K"]}


def test_search_filters_commands(window) -> None:
    dialog = window.open_settings()
    dialog.keyboard.search.setText("redo")
    titles = [
        dialog.keyboard.table.item(r, 0).text() for r in range(dialog.keyboard.table.rowCount())
    ]
    assert titles and all("redo" in t.lower() for t in titles)
    dialog.reject()


def test_policy_locks_settings(tmp_path, monkeypatch, qtbot) -> None:
    policy = tmp_path / "policy.json"
    policy.write_text(
        json.dumps({"settings": {"appearance.theme": "dark", "tool_paths": {"obabel": "/x"}}})
    )
    monkeypatch.setenv("CHEMSUITE_POLICY", str(policy))
    from chemsuite.app.shell.settings_dialog import SettingsDialog
    from chemsuite.app.shell.settings_hub import reset_hub

    reset_hub()
    dialog = SettingsDialog()
    qtbot.addWidget(dialog)
    assert not dialog.appearance.theme.isEnabled()
    assert dialog.appearance.row_height.isEnabled()
    assert not dialog.tools.edits["obabel"].isEnabled()
    assert settings_hub().settings.appearance.theme == "dark"
    assert settings_hub().settings.tool_paths["obabel"] == "/x"
    dialog.reject()


def test_external_tool_test_button(window, tmp_path) -> None:
    dialog = window.open_settings()
    tools = dialog.tools
    tools.edits["mafft"].setText(str(tmp_path / "nope"))
    ok, _message = tools.test("mafft")
    assert not ok
    assert tools.results["mafft"].text().startswith("✗")
    dialog.accept()
    assert settings_hub().settings.tool_paths["mafft"] == str(tmp_path / "nope")


def test_open_recent_and_last_folder(window, tmp_path) -> None:
    path = tmp_path / "a.csv"
    path.write_text("ID\nX\n")
    window._remember_path(path)
    window._remember_path(tmp_path / "gone.csv")
    window._fill_recent_menu()
    texts = [a.text() for a in window._recent_menu.actions()]
    assert texts[:2] == ["gone.csv (missing)", "a.csv"]
    assert window._start_folder() == str(tmp_path)
    window._recent_menu.actions()[-1].trigger()  # Clear Recent
    assert settings_hub().settings.files.recent == []


def test_general_page_workers(window) -> None:
    dialog = window.open_settings()
    dialog.general.workers.setValue(2)
    dialog.accept()
    assert settings_hub().settings.performance.max_workers == 2


def test_light_and_dark_are_fixed_whatever_the_os_uses(qapp) -> None:
    from chemsuite.app.shell import theme
    from chemsuite.core.settings import Appearance

    light = theme.light_palette()
    assert light.color(theme.Role.Base).lightness() == 255
    assert light.color(theme.Role.Window).lightness() > 230
    assert light.color(theme.Role.Text).lightness() < 40
    assert theme.is_dark(theme.dark_palette()) and not theme.is_dark(light)
    appearance = Appearance()
    appearance.theme = "light"
    theme.apply_appearance(qapp, appearance)
    assert not theme.is_dark(qapp.palette())


def test_status_bar_is_slim_while_a_job_runs(window, qtbot) -> None:
    import threading

    window.show()
    bar = window.statusBar()
    idle = bar.height()
    release = threading.Event()
    window.jobs.run("Calculating…", lambda job: release.wait(5), lambda _r: None)
    qtbot.waitUntil(window.jobs_indicator.isVisible)
    assert bar.height() == idle <= bar.fontMetrics().height() + 8
    release.set()
