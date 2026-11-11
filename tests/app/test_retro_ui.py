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

"""Chemistry ▸ Retrosynthesis and the Routes pane (RETRO_PLAN A2 and A3).

The search runs against the fake interpreter from ``tests/retro``, so these tests need no
AiZynthFinder install: what they check is the dialog's rules, the columns written, the routes
reaching the pane, the pane's drawing and the pane coming back with a saved document.
"""

from __future__ import annotations

import sys

import pytest

from tests.retro import recorded
from tests.retro.fake import fake_interpreter, recorded_environment, write_config

from .conftest import gui

pytestmark = [
    gui,
    pytest.mark.skipif(sys.platform == "win32", reason="the fake interpreter is a shell script"),
]

SMILES = [recorded.TARGET, "CCO", ""]


@pytest.fixture
def chem(window, tmp_path, monkeypatch):
    """A window on three structures, with the fake AiZynthFinder configured."""
    from chemsuite.app.main_window import Workspace
    from chemsuite.app.settings_hub import settings_hub
    from chemsuite.core.table import Table

    recorded_environment(
        monkeypatch,
        recorded.TREES,
        recorded.STATISTICS,
        recorded.STOCK_INFO,
        log=tmp_path / "requests.jsonl",
    )
    python = fake_interpreter(tmp_path)
    config = write_config(tmp_path)

    def paths(settings) -> None:
        settings.tool_paths["aizynthfinder_python"] = str(python)
        settings.tool_paths["aizynthfinder_config"] = str(config)

    settings_hub().update(paths)

    table = Table()
    table.add_rows({"ID": ["A", "B", "C"], "SMILES": SMILES}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


@pytest.fixture(autouse=True)
def no_warm_workers():
    """Each test gets its own worker, and none is left running afterwards."""
    from chemsuite.retro import runner

    runner.close_workers()
    yield
    runner.close_workers()


def _run(chem, qtbot, dialog):
    """Run the dialog's search and wait for the columns."""
    job = dialog.run()
    qtbot.waitUntil(lambda: "AiZynth Solved" in chem.workspace.table.column_names, timeout=60_000)
    qtbot.waitUntil(lambda: job.state != "running", timeout=60_000)
    return job


# -- the dialog -----------------------------------------------------------------------------
def test_the_dialog_offers_what_the_config_holds(chem) -> None:
    dialog = chem.open_tool("retrosynthesis")
    assert set(dialog.policies.boxes) == {"uspto", "ringbreaker"}
    assert set(dialog.stocks.boxes) == {"zinc", "bench"}
    assert set(dialog.filters.boxes) == {"uspto"}
    assert all(box.isChecked() for box in dialog.policies.boxes.values())
    assert not dialog.status.text().strip()
    assert dialog.primary.text() == "Search 2 Structures"  # the blank cell is not a target


def test_a_missing_install_blocks_the_run(window) -> None:
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"SMILES": ["CCO"]}, undoable=False)
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    dialog = window.open_tool("retrosynthesis")
    assert "not configured" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_a_config_that_would_leave_the_machine_blocks_the_run(chem, tmp_path) -> None:
    from chemsuite.app.settings_hub import settings_hub

    remote = write_config(
        tmp_path / "remote",
        "expansion:\n  uspto:\n    use_remote_models: True\n    model: m\nstock:\n  zinc: z\n",
    )
    settings_hub().update(lambda s: s.tool_paths.__setitem__("aizynthfinder_config", str(remote)))
    dialog = chem.open_tool("retrosynthesis")
    assert "use_remote_models" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_the_search_writes_the_columns(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("retrosynthesis")
    dialog.rows.set_value("all")
    dialog.time_limit.setValue(5)
    _run(chem, qtbot, dialog)
    first, second, blank = (int(oid) for oid in table.oids)
    assert table.cell("AiZynth Solved", first) == "Yes"
    assert table.cell("AiZynth Score", first) == "0.994"
    assert table.cell("AiZynth Steps", first) == "1"
    assert table.cell("AiZynth Precursors in Stock", second) != ""
    assert table.cell("AiZynth Solved", blank) == ""  # a row with no structure is untouched
    assert table.kind("AiZynth Solved") == "text"


def test_the_options_reach_the_search(chem, qtbot, tmp_path) -> None:
    from tests.retro.fake import requests_logged

    dialog = chem.open_tool("retrosynthesis")
    dialog.rows.set_value("all")
    dialog.time_limit.setValue(9)
    dialog.iterations.setValue(33)
    dialog.max_routes.setValue(4)
    dialog.exploration.setValue(2.0)
    dialog.policies.boxes["ringbreaker"].setChecked(False)
    dialog.stocks.boxes["bench"].setChecked(False)
    _run(chem, qtbot, dialog)
    entries = requests_logged(tmp_path / "requests.jsonl")
    assert entries and entries[0]["time_limit"] == 9
    assert entries[0]["iteration_limit"] == 33
    assert entries[0]["max_routes"] == 4
    assert entries[0]["algorithm_config"]["C"] == 2.0
    assert entries[0]["policies"] == ["uspto"]
    assert entries[0]["stocks"] == ["zinc"]


def test_rows_already_searched_are_skipped(chem, qtbot, tmp_path) -> None:
    from tests.retro.fake import requests_logged

    dialog = chem.open_tool("retrosynthesis")
    dialog.rows.set_value("all")
    _run(chem, qtbot, dialog)
    searched = len(requests_logged(tmp_path / "requests.jsonl"))
    again = chem.open_tool("retrosynthesis")
    again.rows.set_value("all")
    assert "searched already" in again.status.text()
    assert not again.primary.isEnabled()
    again.skip_done.setChecked(False)
    assert not again.status.text().strip()
    assert again.primary.isEnabled()
    assert len(requests_logged(tmp_path / "requests.jsonl")) == searched


def test_a_structure_that_fails_does_not_stop_the_others(chem, qtbot) -> None:
    table = chem.workspace.table
    table.set_cells("SMILES", [int(table.oids[1])], ["FAIL"])
    dialog = chem.open_tool("retrosynthesis")
    dialog.rows.set_value("all")
    _run(chem, qtbot, dialog)
    first, failed, _blank = (int(oid) for oid in table.oids)
    assert table.cell("AiZynth Solved", first) == "Yes"
    assert table.cell("AiZynth Solved", failed) == ""


# -- the Routes pane ------------------------------------------------------------------------
def _routes_panel(window):
    from chemsuite.tools.retro.routes_panel import RoutesPanel

    return next((p for p in window.panels if isinstance(p, RoutesPanel)), None)


def test_the_search_opens_the_routes_pane(chem, qtbot) -> None:
    dialog = chem.open_tool("retrosynthesis")
    dialog.rows.set_value("all")
    _run(chem, qtbot, dialog)
    panel = _routes_panel(chem)
    assert panel is not None
    assert panel.target.count() == 2
    assert panel.routes.rowCount() == 3  # the recorded target has three routes
    assert panel.routes.item(0, 1).text() == "Solved"
    assert panel.title() == "Routes (2)"


def test_the_pane_draws_the_tree_and_opens_a_step(chem, qtbot) -> None:
    from chemsuite.tools.retro.routes_panel import _MoleculeCard, _StepButton

    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    panel = _routes_panel(chem)
    items = panel.canvas.scene().items()
    cards = [item for item in items if isinstance(item, _MoleculeCard)]
    steps = [item for item in items if isinstance(item, _StepButton)]
    assert len(cards) == 3  # the target and its two precursors
    assert len(steps) == 1
    assert {card.in_stock for card in cards} == {True, False}
    # the target sits to the right of its precursors
    target = next(card for card in cards if card.smiles == recorded.TARGET)
    assert all(card.x() < target.x() for card in cards if card is not target)

    panel.show_step(steps[0])
    assert panel.steps.isVisibleTo(panel.splitter)
    fields = [panel.steps.item(row, 0).text() for row in range(panel.steps.rowCount())]
    assert "reaction" in fields
    assert "policy name" in fields
    panel.hide_step()
    assert not panel.steps.isVisibleTo(panel.splitter)


def test_a_two_step_route_draws_both_steps(chem, qtbot) -> None:
    from chemsuite.tools.retro.routes_panel import _StepButton

    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    panel = _routes_panel(chem)
    panel.routes.setCurrentCell(1, 0)  # the two-step route
    steps = [item for item in panel.canvas.scene().items() if isinstance(item, _StepButton)]
    assert sorted(step.number for step in steps) == [1, 2]


def test_a_card_copies_as_smiles(chem, qtbot) -> None:
    from PySide6.QtGui import QGuiApplication

    from chemsuite.tools.retro.routes_panel import _MoleculeCard

    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    panel = _routes_panel(chem)
    card = next(
        item
        for item in panel.canvas.scene().items()
        if isinstance(item, _MoleculeCard) and item.smiles == recorded.AMINOPHENOL
    )
    from chemsuite.app.structure_copy import copy_structures
    from chemsuite.tools.retro.routes_panel import _parse

    copy_structures([_parse(card.smiles)], "smiles")
    assert QGuiApplication.clipboard().text() == recorded.AMINOPHENOL


def test_the_pane_comes_back_with_the_document(chem, qtbot) -> None:
    from chemsuite.app import session
    from chemsuite.app.main_window import Workspace

    dialog = chem.open_tool("retrosynthesis")
    dialog.rows.set_value("all")
    _run(chem, qtbot, dialog)
    panel = _routes_panel(chem)
    panel.routes.setCurrentCell(2, 0)
    extras = session.capture(chem)
    assert any(entry["tool"] == "retro_routes" for entry in extras["panels"])

    chem._install(Workspace(chem.workspace.table, extras=extras))
    session.restore(chem, extras)
    reopened = _routes_panel(chem)
    assert reopened is not None and reopened is not panel
    assert reopened.target.count() == 2
    assert reopened.routes.rowCount() == 3
    assert reopened.routes.currentRow() == 2


# -- the canvas: sharp structures at any zoom ------------------------------------------------
def _card(panel, smiles=None):
    from chemsuite.tools.retro.routes_panel import _MoleculeCard

    return next(
        item
        for item in panel.canvas.scene().items()
        if isinstance(item, _MoleculeCard) and (smiles is None or item.smiles == smiles)
    )


def test_picture_scales_step_by_root_two() -> None:
    from chemsuite.tools.retro.routes_panel import PICTURE_MAX, picture_scale

    assert picture_scale(0.3) == 1.0
    assert picture_scale(1.0) == 1.0
    assert picture_scale(1.2) == pytest.approx(2**0.5)
    assert picture_scale(2.0) == 2.0
    assert picture_scale(2.1) == pytest.approx(2**1.5)
    assert picture_scale(1000) == PICTURE_MAX


def test_a_picture_fills_the_card_at_every_pixel_ratio(qapp) -> None:
    """At ratio 2 the picture has twice the pixels, not half the size (the old Retina bug)."""
    from chemsuite.tools.retro.routes_panel import _draw, _inner_rect

    inner = _inner_rect()
    for scale in (1.0, 2.0, 4.0):
        pixmap = _draw(recorded.TARGET, False, scale)
        assert pixmap.width() == round(inner.width() * scale)
        assert pixmap.height() == round(inner.height() * scale)


def test_zooming_in_redraws_the_structures_sharper(chem, qtbot) -> None:
    from chemsuite.tools.retro.routes_panel import picture_scale

    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    panel = _routes_panel(chem)
    canvas = panel.canvas
    canvas.resize(600, 400)
    canvas.viewport().grab()  # paint once at the fitted zoom
    ratio = canvas.viewport().devicePixelRatioF()
    canvas.zoom_by(4.0 / canvas.zoom)
    assert canvas.zoom == pytest.approx(4.0)
    target = _card(panel, recorded.TARGET)
    canvas.centerOn(target)
    canvas.viewport().grab()  # the target asks for a sharper picture
    wanted = picture_scale(4.0 * ratio)
    made = panel._pictures[(recorded.TARGET, False)]
    qtbot.waitUntil(lambda: wanted in made, timeout=5000)
    assert made[wanted].width() >= 4 * ratio * (target.boundingRect().width() - 12) - 1
    assert len(made) <= 3


def test_zoom_stays_within_its_limits_and_holds_the_pointer(chem, qtbot) -> None:
    from PySide6.QtCore import QPointF

    from chemsuite.tools.retro.routes_panel import ZOOM_MAX, ZOOM_MIN

    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    canvas = _routes_panel(chem).canvas
    canvas.resize(200, 150)
    canvas.zoom_by(1000)
    assert canvas.zoom == pytest.approx(ZOOM_MAX)
    canvas.zoom_by(1e-6)
    assert canvas.zoom == pytest.approx(ZOOM_MIN)
    canvas.zoom_by(1 / ZOOM_MIN)  # back to 1, where the route is larger than the view
    canvas.centerOn(canvas.scene().sceneRect().center())
    at = QPointF(110, 80)
    before = canvas.mapToScene(at.toPoint())
    canvas.zoom_by(1.5, at)
    after = canvas.mapToScene(at.toPoint())
    assert abs(after.x() - before.x()) < 2 and abs(after.y() - before.y()) < 2


def test_scroll_pans_and_command_scroll_zooms(chem, qtbot) -> None:
    from PySide6.QtCore import QPoint, QPointF, Qt
    from PySide6.QtGui import QWheelEvent
    from PySide6.QtWidgets import QApplication

    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    canvas = _routes_panel(chem).canvas
    canvas.resize(300, 200)
    canvas.zoom_by(2.0)

    def wheel(dy: int, modifiers) -> None:
        pos = QPointF(100, 100)
        event = QWheelEvent(
            pos,
            QPointF(canvas.viewport().mapToGlobal(pos.toPoint())),
            QPoint(0, 0),
            QPoint(0, dy),
            Qt.MouseButton.NoButton,
            modifiers,
            Qt.ScrollPhase.NoScrollPhase,
            False,
        )
        QApplication.sendEvent(canvas.viewport(), event)

    zoom = canvas.zoom
    bar = canvas.verticalScrollBar()
    start = bar.value()
    wheel(-240, Qt.KeyboardModifier.NoModifier)
    assert canvas.zoom == pytest.approx(zoom)
    assert bar.value() != start
    wheel(120, Qt.KeyboardModifier.ControlModifier)
    assert canvas.zoom > zoom


def test_the_route_is_fitted_again_on_resize_until_moved(chem, qtbot) -> None:
    dialog = chem.open_tool("retrosynthesis")
    _run(chem, qtbot, dialog)
    canvas = _routes_panel(chem).canvas
    canvas.setParent(None)  # its own window, so resizes arrive at once
    qtbot.addWidget(canvas)
    canvas.show()
    qtbot.waitExposed(canvas)
    canvas.resize(120, 80)
    qtbot.waitUntil(lambda: canvas.zoom < 1.0)
    canvas.resize(2000, 1200)
    qtbot.waitUntil(lambda: canvas.zoom == pytest.approx(1.0))  # never magnified past natural
    canvas.zoom_by(2.0)
    canvas.resize(120, 80)
    assert canvas.zoom == pytest.approx(2.0)  # the user's zoom is kept
    canvas.fit()
    assert canvas.zoom < 1.0


# -- Settings ▸ External Tools ▸ Download… ---------------------------------------------------
class _Response:
    def __init__(self, body: bytes) -> None:
        import io

        self._body = io.BytesIO(body)
        self.status = 200
        self.headers = {"Content-Length": str(len(body))}

    def read(self, size: int = -1) -> bytes:
        return self._body.read(size)

    def close(self) -> None:
        pass


def test_download_fetches_the_data_and_sets_the_config(window, qtbot, tmp_path, monkeypatch):
    from chemsuite.app.settings_hub import settings_hub
    from chemsuite.core import download
    from chemsuite.retro.config import CONFIG_TOOL

    urls = []

    def opener(request, timeout=None):
        urls.append(request.full_url)
        return _Response(b"data")

    monkeypatch.setattr(download, "_open", opener)  # every download opens through here
    dialog = window.open_settings()
    tools = dialog.tools
    assert tools.download_button.isEnabled()
    job = tools.start_download(tmp_path / "aizynth")
    assert not tools.download_button.isEnabled()  # one download at a time
    qtbot.waitUntil(lambda: not job.running, timeout=10_000)
    config = (tmp_path / "aizynth").resolve() / "config.yml"
    assert config.is_file()
    assert len(urls) == 6
    assert tools.edits[CONFIG_TOOL].text() == str(config)
    assert tools.results[CONFIG_TOOL].text().startswith("✓")
    assert tools.download_button.isEnabled()
    assert settings_hub().settings.tool_paths[CONFIG_TOOL] == str(config)
    dialog.reject()
    assert settings_hub().settings.tool_paths[CONFIG_TOOL] == str(config)  # kept on Cancel


def test_download_failure_is_reported_on_the_row(window, qtbot, tmp_path, monkeypatch):
    import urllib.error

    from chemsuite.core import download
    from chemsuite.retro.config import CONFIG_TOOL

    def opener(request, timeout=None):
        raise urllib.error.URLError("no route to host")

    monkeypatch.setattr(download, "_open", opener)  # every download opens through here
    dialog = window.open_settings()
    job = dialog.tools.start_download(tmp_path / "aizynth")
    qtbot.waitUntil(lambda: not job.running, timeout=10_000)
    text = dialog.tools.results[CONFIG_TOOL].text()
    assert text.startswith("✗") and "no route to host" in text
    assert "DownloadError" not in text
    dialog.reject()


def test_policy_can_turn_downloads_off(tmp_path, monkeypatch, qtbot) -> None:
    import json

    from chemsuite.app.settings_dialog import NETWORK_OFF_TIP, SettingsDialog
    from chemsuite.app.settings_hub import reset_hub

    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"network_allowed": False}))
    monkeypatch.setenv("CHEMSUITE_POLICY", str(policy))
    reset_hub()
    dialog = SettingsDialog()
    qtbot.addWidget(dialog)
    assert not dialog.tools.download_button.isEnabled()
    assert dialog.tools.download_button.toolTip() == NETWORK_OFF_TIP
    dialog.reject()
