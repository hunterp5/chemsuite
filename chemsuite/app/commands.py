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
"""Every menu command by a stable id, so its keyboard shortcut can be changed in Settings.

Windows register each command once (id, title, menu, default shortcuts); the user's overrides
live in ``Settings.shortcuts`` and replace a command's defaults entirely (an empty list removes
its shortcut). Shortcuts are stored in Qt's portable text form ("Ctrl+Shift+D"), which Qt shows
as ⌘⇧D on macOS.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from PySide6.QtGui import QAction, QKeySequence

from chemsuite.core.settings import Settings

Portable = QKeySequence.SequenceFormat.PortableText


@dataclass(frozen=True)
class Command:
    id: str
    title: str
    menu: str  # where it lives, e.g. "Data ▸ SAR"
    defaults: tuple[str, ...] = ()


_commands: dict[str, Command] = {}


def portable(sequences: Iterable[QKeySequence | QKeySequence.StandardKey | str]) -> tuple[str, ...]:
    out: list[str] = []
    for item in sequences:
        if isinstance(item, QKeySequence.StandardKey):
            keys = QKeySequence.keyBindings(item)
        else:
            keys = [QKeySequence(item)]
        for key in keys:
            text = key.toString(Portable)
            if text and text not in out:
                out.append(text)
    return tuple(out)


def register(command_id: str, title: str, menu: str, defaults=()) -> Command:
    command = Command(command_id, title.replace("&", ""), menu.replace("&", ""), portable(defaults))
    _commands[command_id] = command
    return command


def commands() -> list[Command]:
    return sorted(_commands.values(), key=lambda c: (c.menu, c.title))


def command(command_id: str) -> Command | None:
    return _commands.get(command_id)


def shortcuts_for(command_id: str, settings: Settings) -> list[str]:
    if command_id in settings.shortcuts:
        return list(settings.shortcuts[command_id])
    known = _commands.get(command_id)
    return list(known.defaults) if known else []


def apply_shortcuts(actions: Mapping[str, QAction], settings: Settings) -> None:
    for command_id, action in actions.items():
        action.setShortcuts([QKeySequence(s) for s in shortcuts_for(command_id, settings)])


def conflicts(shortcuts: Mapping[str, list[str]]) -> dict[str, list[str]]:
    """``{sequence: [command ids]}`` for every sequence bound to more than one command."""
    owners: dict[str, list[str]] = defaultdict(list)
    for command_id, sequences in shortcuts.items():
        for sequence in sequences:
            owners[QKeySequence(sequence).toString(Portable)].append(command_id)
    return {seq: ids for seq, ids in owners.items() if seq and len(ids) > 1}


def all_shortcuts(settings: Settings) -> dict[str, list[str]]:
    return {c.id: shortcuts_for(c.id, settings) for c in commands()}
