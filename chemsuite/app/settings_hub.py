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
"""The application's settings: loaded once, saved on change, broadcast to every window."""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QObject, QStandardPaths, Signal
from PySide6.QtWidgets import QApplication

from chemsuite.core.settings import Policy, Settings, SettingsStore

from .theme import apply_appearance

CONFIG_ENV = "CHEMSUITE_CONFIG_DIR"  # tests and portable installs point this elsewhere
CACHE_ENV = "CHEMSUITE_CACHE_DIR"


def config_folder() -> Path:
    override = os.environ.get(CONFIG_ENV, "").strip()
    if override:
        return Path(override)
    location = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppConfigLocation)
    return Path(location or Path.home() / ".chemsuite")


def cache_folder() -> Path:
    """Where downloaded files that can be fetched again are kept."""
    override = os.environ.get(CACHE_ENV, "").strip()
    if override:
        return Path(override)
    location = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.CacheLocation)
    return Path(location or Path.home() / ".chemsuite" / "cache")


class SettingsHub(QObject):
    """Owns the store. ``changed`` carries the *effective* settings (policy applied)."""

    changed = Signal(object)

    def __init__(self, store: SettingsStore, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.store = store
        self._effective = store.effective()

    @property
    def settings(self) -> Settings:
        """Effective settings: read these to decide how things look and behave."""
        return self._effective

    @property
    def user(self) -> Settings:
        """The user's own choices (what the Settings window edits)."""
        return self.store.user

    def locked(self, key: str) -> bool:
        return self.store.locked(key)

    def update(self, change: Callable[[Settings], None]) -> Settings:
        """Apply *change* to a copy of the user's settings, save, and broadcast."""
        user = self.store.user.copy()
        change(user)
        return self.save(user)

    def save(self, user: Settings) -> Settings:
        self.store.save(user)
        return self._broadcast(self.store.effective())

    def preview(self, user: Settings, *, reroll: bool = False) -> Settings:
        """Show *user* settings everywhere without saving (the Settings window's live preview)."""
        return self._broadcast(self.store.policy.apply(user), reroll=reroll)

    def revert(self) -> Settings:
        return self._broadcast(self.store.effective())

    def _broadcast(self, effective: Settings, *, reroll: bool = False) -> Settings:
        self._effective = effective
        app = QApplication.instance()
        if app is not None:
            apply_appearance(app, effective.appearance, reroll=reroll)
        self.changed.emit(effective)
        return effective


_hub: SettingsHub | None = None


def settings_hub() -> SettingsHub:
    """The application's hub, created (and the theme applied) on first use."""
    global _hub
    if _hub is None:
        _hub = SettingsHub(SettingsStore(config_folder(), Policy.load()))
        app = QApplication.instance()
        if app is not None:
            apply_appearance(app, _hub.settings.appearance)
    return _hub


def reset_hub() -> None:
    """Forget the hub (tests: the next use reads the config folder again)."""
    global _hub
    _hub = None
