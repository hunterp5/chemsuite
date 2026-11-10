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
"""Settings…: appearance, keyboard shortcuts, external tools and general options.

Changes preview live in every window; OK saves them, Cancel puts everything back. Settings the
IT policy fixes are shown disabled with a padlock.
"""

from __future__ import annotations

import os
from pathlib import Path

import shiboken6
from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QKeySequenceEdit,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core import external_tools as et
from chemsuite.core.jobs.parallel import default_workers
from chemsuite.core.settings import (
    ALIGN_H,
    ALIGN_V,
    BUILTIN_THEMES,
    GROOVY,
    MAX_FONT_PT,
    MAX_ROW_HEIGHT,
    MIN_FONT_PT,
    MIN_ROW_HEIGHT,
    custom_theme_id,
    custom_theme_name,
)
from chemsuite.retro.config import CONFIG_TOOL as RETRO_CONFIG_TOOL
from chemsuite.retro.config import PYTHON_TOOL as RETRO_PYTHON_TOOL

from .commands import Portable, all_shortcuts, commands, conflicts
from .settings_hub import settings_hub

LOCK = "🔒 "
LOCKED_TIP = "Set by your IT department (machine policy)"
THEME_TITLES = {"system": "System", "light": "Light", "dark": "Dark", "groovy": "Groovy"}


class SettingsDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.hub = settings_hub()
        self.user = self.hub.user.copy()
        self._reroll = False
        self.pages = QListWidget()
        self.pages.setFixedWidth(150)
        self.stack = QStackedWidget()
        self.appearance = AppearancePage(self)
        self.keyboard = KeyboardPage(self)
        self.tools = ToolsPage(self)
        self.general = GeneralPage(self)
        for title, page in (
            ("Appearance", self.appearance),
            ("Keyboard", self.keyboard),
            ("External Tools", self.tools),
            ("General", self.general),
        ):
            self.pages.addItem(title)
            self.stack.addWidget(page)
        self.pages.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.pages.setCurrentRow(0)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.RestoreDefaults
        )
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).setText("Reset Page")
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(
            self.reset_page
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        body = QHBoxLayout()
        body.addWidget(self.pages)
        body.addWidget(self.stack, 1)
        layout = QVBoxLayout(self)
        layout.addLayout(body, 1)
        layout.addWidget(buttons)
        self.resize(820, 560)

    # -- helpers for the pages -------------------------------------------------------------------
    def lock(self, widget: QWidget, key: str) -> bool:
        locked = self.hub.locked(key)
        if locked:
            widget.setEnabled(False)
            widget.setToolTip(f"{LOCKED_TIP}: {self.hub.store.policy.path}")
        return locked

    def changed(self) -> None:
        """A page edited ``self.user``: show it everywhere now."""
        self.hub.preview(self.user, reroll=self._reroll)
        self._reroll = False

    def reroll(self) -> None:
        self._reroll = True
        self.changed()

    def reset_page(self) -> None:
        self.stack.currentWidget().reset()
        self.changed()

    # -- closing ----------------------------------------------------------------------------------
    def problem(self) -> str:
        clashes = conflicts(all_shortcuts(self.user))
        if clashes:
            sequence, ids = next(iter(clashes.items()))
            names = " and ".join(f"“{_title(i)}”" for i in ids)
            return f"{QKeySequence(sequence).toString()} is used by {names}."
        return ""

    def accept(self) -> None:
        problem = self.problem()
        if problem:
            self.pages.setCurrentRow(1)
            self.keyboard.status.setText(problem)
            return
        self.hub.save(self.user)
        super().accept()

    def reject(self) -> None:
        self.hub.revert()
        super().reject()


def _title(command_id: str) -> str:
    known = {c.id: c for c in commands()}.get(command_id)
    return known.title if known else command_id


class AppearancePage(QWidget):
    def __init__(self, dialog: SettingsDialog) -> None:
        super().__init__()
        self.dialog = dialog
        a = dialog.user.appearance
        self.theme = QComboBox()
        self.reroll = QPushButton("New Colours")
        self.reroll.setToolTip("Roll a new Groovy palette.")
        self.reroll.clicked.connect(dialog.reroll)
        edit = QPushButton("Custom Themes…")
        edit.clicked.connect(self.edit_custom)
        theme_row = QHBoxLayout()
        theme_row.addWidget(self.theme, 1)
        theme_row.addWidget(self.reroll)
        theme_row.addWidget(edit)
        self.app_font = _font_spin("System size")
        self.table_font = _font_spin("Same as application")
        self.align_h = QComboBox()
        self.align_v = QComboBox()
        for combo, choices in ((self.align_h, ALIGN_H), (self.align_v, ALIGN_V)):
            for choice in choices:
                combo.addItem(choice.title(), choice)
        self.row_height = QSpinBox()
        self.row_height.setRange(MIN_ROW_HEIGHT, MAX_ROW_HEIGHT)
        self.row_height.setSuffix(" px")
        self.status_bar = QCheckBox("Show the status bar")
        for combo in (self.align_h, self.align_v):
            combo.setToolTip("Applies to text cells; numbers stay right-aligned.")
        align = QHBoxLayout()
        align.addWidget(self.align_h)
        align.addWidget(self.align_v)
        form = QFormLayout(self)
        form.addRow("Theme:", theme_row)
        form.addRow("Application font:", self.app_font)
        form.addRow("Table font:", self.table_font)
        form.addRow("Text alignment:", align)
        form.addRow("Structure row height:", self.row_height)
        form.addRow("", self.status_bar)
        self.load(a)
        for widget, key in (
            (self.theme, "appearance.theme"),
            (self.reroll, "appearance.theme"),
            (edit, "appearance.custom_themes"),
            (self.app_font, "appearance.app_font_pt"),
            (self.table_font, "appearance.table_font_pt"),
            (self.align_h, "appearance.table_align_h"),
            (self.align_v, "appearance.table_align_v"),
            (self.row_height, "appearance.structure_row_height"),
            (self.status_bar, "appearance.status_bar"),
        ):
            dialog.lock(widget, key)
        self.theme.currentIndexChanged.connect(self._store)
        for spin in (self.app_font, self.table_font, self.row_height):
            spin.valueChanged.connect(self._store)
        for combo in (self.align_h, self.align_v):
            combo.currentIndexChanged.connect(self._store)
        self.status_bar.toggled.connect(self._store)

    def _fill_themes(self, current: str) -> None:
        self.theme.blockSignals(True)
        self.theme.clear()
        for key in BUILTIN_THEMES:
            self.theme.addItem(THEME_TITLES[key], key)
        names = sorted(self.dialog.user.appearance.custom_themes)
        if names:
            self.theme.insertSeparator(self.theme.count())
            for name in names:
                self.theme.addItem(name, custom_theme_id(name))
        self.theme.setCurrentIndex(max(0, self.theme.findData(current)))
        self.theme.blockSignals(False)

    def load(self, a) -> None:
        widgets = (self.app_font, self.table_font, self.row_height, self.align_h, self.align_v)
        for w in (*widgets, self.status_bar):
            w.blockSignals(True)
        self._fill_themes(a.theme)
        self.app_font.setValue(a.app_font_pt or MIN_FONT_PT - 1)
        self.table_font.setValue(a.table_font_pt or MIN_FONT_PT - 1)
        self.align_h.setCurrentIndex(self.align_h.findData(a.table_align_h))
        self.align_v.setCurrentIndex(self.align_v.findData(a.table_align_v))
        self.row_height.setValue(a.structure_row_height)
        self.status_bar.setChecked(a.status_bar)
        for w in (*widgets, self.status_bar):
            w.blockSignals(False)
        self.reroll.setVisible(a.theme == GROOVY)

    def _store(self, *_args) -> None:
        a = self.dialog.user.appearance
        a.theme = self.theme.currentData() or a.theme
        a.app_font_pt = _font_value(self.app_font)
        a.table_font_pt = _font_value(self.table_font)
        a.table_align_h = self.align_h.currentData()
        a.table_align_v = self.align_v.currentData()
        a.structure_row_height = self.row_height.value()
        a.status_bar = self.status_bar.isChecked()
        self.reroll.setVisible(a.theme == GROOVY)
        self.dialog.changed()

    def edit_custom(self) -> None:
        from .custom_theme_dialog import CustomThemeDialog

        a = self.dialog.user.appearance
        editor = CustomThemeDialog(a.custom_themes, custom_theme_name(a.theme), self)
        if editor.exec() != QDialog.DialogCode.Accepted:
            return
        a.custom_themes = editor.themes
        if editor.chosen:
            a.theme = custom_theme_id(editor.chosen)
        elif custom_theme_name(a.theme) not in (None, *a.custom_themes):
            a.theme = "system"  # the theme in use was deleted
        self._fill_themes(a.theme)
        self.dialog.changed()

    def reset(self) -> None:
        from chemsuite.core.settings import Appearance

        a = self.dialog.user.appearance
        fresh = Appearance(custom_themes=a.custom_themes)  # keep saved custom themes
        for key in ("theme", "app_font_pt", "table_font_pt", "table_align_h", "table_align_v",
                    "structure_row_height", "status_bar"):  # fmt: skip
            setattr(a, key, getattr(fresh, key))
        self.load(a)


def _font_spin(zero: str) -> QSpinBox:
    spin = QSpinBox()
    spin.setRange(MIN_FONT_PT - 1, MAX_FONT_PT)
    spin.setSpecialValueText(zero)  # the minimum stands for "0: default"
    spin.setSuffix(" pt")
    return spin


def _font_value(spin: QSpinBox) -> int:
    return 0 if spin.value() < MIN_FONT_PT else spin.value()


class KeyboardPage(QWidget):
    COLUMNS = ("Command", "Menu", "Shortcut")

    def __init__(self, dialog: SettingsDialog) -> None:
        super().__init__()
        self.dialog = dialog
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search commands or shortcuts")
        self.search.setClearButtonEnabled(True)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(list(self.COLUMNS))
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.verticalHeader().hide()
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.table.verticalHeader().setDefaultSectionSize(self.fontMetrics().height() + 8)
        self.editor = QKeySequenceEdit()
        self.editor.setMaximumSequenceLength(1)
        add = QPushButton("Add")
        add.setToolTip("Add the recorded shortcut to the chosen command.")
        add.clicked.connect(self.add_recorded)
        remove = QPushButton("Remove All")
        remove.setToolTip("No shortcut for the chosen command.")
        remove.clicked.connect(lambda: self.set_shortcuts(self.current(), []))
        default = QPushButton("Default")
        default.clicked.connect(lambda: self.set_shortcuts(self.current(), None))
        row = QHBoxLayout()
        row.addWidget(QLabel("Press keys:"))
        row.addWidget(self.editor, 1)
        row.addWidget(add)
        row.addWidget(remove)
        row.addWidget(default)
        self.status = QLabel()
        self.status.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.search)
        layout.addWidget(self.table, 1)
        layout.addLayout(row)
        layout.addWidget(self.status)
        self.locked = dialog.lock(self.table, "shortcuts")
        if self.locked:
            for w in (self.editor, add, remove, default):
                w.setEnabled(False)
        self.search.textChanged.connect(lambda *_: self.fill())
        self.fill()

    def fill(self) -> None:
        needle = self.search.text().strip().lower()
        shortcuts = all_shortcuts(self.dialog.user)
        clashes = {cid for ids in conflicts(shortcuts).values() for cid in ids}
        rows = []
        for c in commands():
            native = ", ".join(QKeySequence(s).toString() for s in shortcuts[c.id])
            text = f"{c.title} {c.menu} {native}".lower()
            if needle and needle not in text:
                continue
            rows.append((c, native))
        current = self.current()
        self.table.setRowCount(len(rows))
        for r, (c, native) in enumerate(rows):
            for col, value in enumerate((c.title, c.menu, native)):
                item = QTableWidgetItem(value)
                item.setData(Qt.ItemDataRole.UserRole, c.id)
                if c.id in clashes:
                    item.setForeground(Qt.GlobalColor.red)
                    item.setToolTip("This shortcut is also used by another command.")
                elif c.id in self.dialog.user.shortcuts:
                    font = item.font()
                    font.setBold(col == 2)
                    item.setFont(font)
                self.table.setItem(r, col, item)
            if c.id == current:
                self.table.selectRow(r)
        self.status.setText(self.dialog.problem())

    def current(self) -> str | None:
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return None
        return self.table.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole)

    def select(self, command_id: str) -> None:
        for r in range(self.table.rowCount()):
            if self.table.item(r, 0).data(Qt.ItemDataRole.UserRole) == command_id:
                self.table.selectRow(r)
                return

    def set_shortcuts(self, command_id: str | None, sequences: list[str] | None) -> None:
        """*sequences* None: back to the command's defaults."""
        if command_id is None or self.locked:
            return
        if sequences is None:
            self.dialog.user.shortcuts.pop(command_id, None)
        else:
            self.dialog.user.shortcuts[command_id] = sequences
        self.fill()
        self.dialog.changed()

    def add_recorded(self, sequence: QKeySequence | None = None) -> None:
        command_id = self.current()
        sequence = sequence or self.editor.keySequence()
        text = sequence.toString(Portable)
        if command_id is None or not text:
            return
        current = all_shortcuts(self.dialog.user)[command_id]
        if text not in current:
            self.set_shortcuts(command_id, [*current, text])
        self.editor.clear()

    def reset(self) -> None:
        if not self.locked:
            self.dialog.user.shortcuts.clear()
            self.fill()


RETRO_TOOLS = frozenset({RETRO_PYTHON_TOOL, RETRO_CONFIG_TOOL})
NETWORK_OFF_TIP = "Downloads are turned off by your IT department (machine policy)."
#: the AiZynthFinder data download in progress, so a reopened Settings cannot start a second
_retro_download = None


class ToolsPage(QWidget):
    def __init__(self, dialog: SettingsDialog) -> None:
        super().__init__()
        self.dialog = dialog
        self.edits: dict[str, QLineEdit] = {}
        self.results: dict[str, QLabel] = {}
        form = QFormLayout(self)
        for tool in et.visible_tools():
            edit = QLineEdit(dialog.user.tool_paths.get(tool.id, ""))
            found = et.resolve(tool.id, type(dialog.user)()) if tool.names else None
            edit.setPlaceholderText(f"Found: {found}" if found else "Not found")
            edit.setToolTip(tool.help)
            browse = QPushButton("Browse…")
            browse.clicked.connect(lambda _c=False, t=tool: self.browse(t))
            test = QPushButton("Test")
            test.clicked.connect(lambda _c=False, t=tool: self.test(t.id))
            result = QLabel()
            result.setWordWrap(True)
            result.hide()  # shown once Test has something to say
            row = QHBoxLayout()
            row.addWidget(edit, 1)
            row.addWidget(browse)
            if tool.id == RETRO_CONFIG_TOOL:
                row.addWidget(self._download_button())
            row.addWidget(test)
            box = QVBoxLayout()
            box.addLayout(row)
            box.addWidget(result)
            label = QLabel(tool.title + ":")
            label.setToolTip(tool.help)
            form.addRow(label, box)
            if dialog.lock(edit, f"tool_paths.{tool.id}"):
                browse.setEnabled(False)
                label.setText(LOCK + label.text())
            edit.textChanged.connect(lambda text, t=tool.id: self._store(t, text))
            self.edits[tool.id] = edit
            self.results[tool.id] = result

    def _store(self, tool_id: str, text: str) -> None:
        if text.strip():
            self.dialog.user.tool_paths[tool_id] = text.strip()
        else:
            self.dialog.user.tool_paths.pop(tool_id, None)

    def browse(self, tool: et.ExternalTool) -> None:
        start = self.edits[tool.id].text() or os.path.expanduser("~")
        if tool.kind == "folder":
            path = QFileDialog.getExistingDirectory(self, tool.title, start)
        else:
            path, _ = QFileDialog.getOpenFileName(self, tool.title, start)
        if path:
            self.edits[tool.id].setText(path)

    def test(self, tool_id: str) -> tuple[bool, str]:
        path = self.edits[tool_id].text().strip() or et.resolve(tool_id, self.dialog.user) or ""
        if tool_id in RETRO_TOOLS:
            ok, message = self._test_retro(tool_id)
        else:
            ok, message = et.test_tool(tool_id, path) if path else (False, "Not found.")
        self.results[tool_id].setText(("✓ " if ok else "✗ ") + message)
        self.results[tool_id].show()
        return ok, message

    def _test_retro(self, tool_id: str) -> tuple[bool, str]:
        """Test AiZynthFinder itself: the interpreter must have it, the config must be local.

        A plain ``--version`` would only prove the interpreter runs, so the worker is asked
        for AiZynthFinder's version instead; testing the config checks it alongside.
        """
        from chemsuite.retro.runner import probe

        python = self.edits[RETRO_PYTHON_TOOL].text().strip()
        config = self.edits[RETRO_CONFIG_TOOL].text().strip()
        if tool_id == RETRO_CONFIG_TOOL and not config:
            return False, "Not found."
        if not python:
            return False, "Set the AiZynthFinder Python interpreter as well."
        guard = QGuiApplication.setOverrideCursor, QGuiApplication.restoreOverrideCursor
        guard[0](Qt.CursorShape.WaitCursor)  # importing AiZynthFinder takes a few seconds
        try:
            return probe(python, config if tool_id == RETRO_CONFIG_TOOL and config else None)
        finally:
            guard[1]()

    # -- AiZynthFinder's public data -------------------------------------------------------------
    def _download_button(self) -> QPushButton:
        self.download_button = QPushButton("Download…")
        self.download_button.setToolTip(
            "Download AiZynthFinder's public USPTO policies and ZINC stock (about 1.5 GB) "
            "into a folder, and use the config.yml written beside them."
        )
        self.download_button.clicked.connect(self.download)
        if self.dialog.hub.store.policy.get("network_allowed", True) is False:
            self.download_button.setEnabled(False)
            self.download_button.setToolTip(NETWORK_OFF_TIP)
        elif self.dialog.hub.locked(f"tool_paths.{RETRO_CONFIG_TOOL}"):
            self.download_button.setEnabled(False)
            self.download_button.setToolTip(LOCKED_TIP)
        elif _retro_download is not None:
            self.download_button.setEnabled(False)
        return self.download_button

    def download(self) -> None:
        current = self.edits[RETRO_CONFIG_TOOL].text().strip()
        start = str(Path(current).expanduser().parent) if current else os.path.expanduser("~")
        folder = QFileDialog.getExistingDirectory(self, "Download AiZynthFinder Data To", start)
        if not folder:
            return
        from chemsuite.retro.download import missing

        if missing(Path(folder)):
            answer = QMessageBox.question(
                self,
                "Download AiZynthFinder Data",
                f"Download about 1.5 GB of AiZynthFinder models and stock to {folder}?",
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
        self.start_download(Path(folder))

    def start_download(self, folder: Path):
        """Fetch the data in the background; the config path is set when it is written."""
        global _retro_download
        from chemsuite.retro.download import download_public_data

        from .jobs import CANCELLED, FAILED, QUEUE

        window = self.dialog.parentWidget()
        result = self.results[RETRO_CONFIG_TOOL]
        if not hasattr(window, "_run_with_progress"):
            result.setText("✗ Open Settings from a document window to download.")
            result.show()
            return None

        def work(job):
            return download_public_data(folder, cancel=job.cancel_token, progress=job.report)

        def done(path: Path) -> None:
            def store(settings) -> None:
                settings.tool_paths[RETRO_CONFIG_TOOL] = str(path)

            settings_hub().update(store)
            job.set_summary(f"Wrote {path}")
            if shiboken6.isValid(self) and self.dialog.isVisible():
                self.edits[RETRO_CONFIG_TOOL].setText(str(path))
                result.setText(f"✓ Downloaded; config {path}")
                result.show()
            window.show_message(f"AiZynthFinder data downloaded; config set to {path}", 8000)

        def ended() -> None:
            global _retro_download
            _retro_download = None
            if not shiboken6.isValid(self):
                return
            self.download_button.setEnabled(True)
            if job.state == FAILED:
                last = (job.error_text or "").strip().splitlines()[-1:] or ["Download failed."]
                result.setText("✗ " + (last[0].partition(": ")[2] or last[0]))
            elif job.state == CANCELLED:
                result.setText("✗ Cancelled; Download again resumes where it stopped.")

        job = window._run_with_progress("Downloading AiZynthFinder data…", work, done, lane=QUEUE)
        _retro_download = job
        job.done.connect(ended)
        self.download_button.setEnabled(False)
        result.setText("Downloading…")
        result.show()
        return job

    def reset(self) -> None:
        for tool_id, edit in self.edits.items():
            if edit.isEnabled():
                edit.clear()
                self.dialog.user.tool_paths.pop(tool_id, None)


class GeneralPage(QWidget):
    def __init__(self, dialog: SettingsDialog) -> None:
        super().__init__()
        self.dialog = dialog
        self.workers = QSpinBox()
        self.workers.setRange(0, max(64, os.cpu_count() or 1))
        self.workers.setSpecialValueText(f"Automatic ({default_workers()})")
        self.workers.setToolTip(
            "Processor cores used by heavy work (descriptors, fingerprints, conformers…). "
            "Takes effect after restarting chemsuite."
        )
        self.recent_limit = QSpinBox()
        self.recent_limit.setRange(0, 50)
        clear = QPushButton("Clear Recent Files")
        clear.clicked.connect(self.clear_recent)
        self.notify = QCheckBox("Notify when a long job finishes")
        self.notify.setToolTip(
            "A desktop notification when a background job of more than 30 seconds ends "
            "while chemsuite is not in front."
        )
        form = QFormLayout(self)
        form.addRow("Processor cores:", self.workers)
        form.addRow("Recent files to keep:", self.recent_limit)
        form.addRow("", clear)
        form.addRow("Jobs:", self.notify)
        self.load()
        dialog.lock(self.workers, "performance.max_workers")
        dialog.lock(self.recent_limit, "files.recent_limit")
        dialog.lock(self.notify, "jobs.notify")
        self.workers.valueChanged.connect(self._store)
        self.recent_limit.valueChanged.connect(self._store)
        self.notify.toggled.connect(self._store)

    def load(self) -> None:
        user = self.dialog.user
        self.workers.setValue(user.performance.max_workers)
        self.recent_limit.setValue(user.files.recent_limit)
        self.notify.setChecked(user.jobs.notify)

    def _store(self, *_args) -> None:
        user = self.dialog.user
        user.performance.max_workers = self.workers.value()
        user.files.recent_limit = self.recent_limit.value()
        user.files.recent = user.files.recent[: user.files.recent_limit]
        user.jobs.notify = self.notify.isChecked()

    def clear_recent(self) -> None:
        self.dialog.user.files.recent = []

    def reset(self) -> None:
        self.workers.setValue(0)
        self.recent_limit.setValue(10)
        self.notify.setChecked(True)
