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

"""Which ligand pKa model a preparation runs with, and when chemsuite pKa is downloaded."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from PySide6.QtWidgets import QMessageBox

from chemsuite.ionisation.unipka import weights
from chemsuite.macromolecule.prepare.ligand_pka import ModelStates
from chemsuite.tools.macromolecule import pka_model


class _Ctx:
    def __init__(self, folder: Path, network: bool = True) -> None:
        self.folder, self.network, self.parent = folder, network, None

    def cache_folder(self, name: str) -> Path:
        return self.folder / name

    def network_allowed(self) -> bool:
        return self.network

    def find_panel(self, _tool_id: str):
        return None


JOB = SimpleNamespace(cancel_token=None, report=lambda done, total: None)


@pytest.fixture
def asked(monkeypatch) -> list[str]:
    monkeypatch.setattr(pka_model, "_declined", False)
    out: list[str] = []
    answer = QMessageBox.StandardButton.No

    def question(_parent, _title, text):
        out.append(text)
        return answer

    monkeypatch.setattr(QMessageBox, "question", question)
    return out


def test_declining_the_download_uses_molgpka_and_is_not_asked_again(tmp_path, asked) -> None:
    ctx = _Ctx(tmp_path)
    assert pka_model.ligand_pka(ctx)(JOB) is None
    assert pka_model.ligand_pka(ctx)(JOB) is None
    assert len(asked) == 1 and "88 MB" in asked[0]


def test_without_the_network_nothing_is_asked(tmp_path, asked) -> None:
    assert pka_model.ligand_pka(_Ctx(tmp_path, network=False))(JOB) is None
    assert not asked


def test_weights_in_the_cache_are_used_without_asking(tmp_path, asked, monkeypatch) -> None:
    found = tmp_path / "pka" / weights.NAME
    monkeypatch.setattr(weights, "cached", lambda folder: found)
    monkeypatch.setattr(weights, "predictor", lambda path: SimpleNamespace(mean=6.0, path=path))
    model = pka_model.ligand_pka(_Ctx(tmp_path))(JOB)
    assert isinstance(model, ModelStates) and model.predictor.path == found
    assert not asked
