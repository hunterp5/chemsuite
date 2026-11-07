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

"""The retrosynthesis core (RETRO_PLAN A1): the config rules, the parsers, the worker."""

from __future__ import annotations

import sys
import threading

import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.settings import Settings
from chemsuite.retro import routes as routes_module
from chemsuite.retro.config import (
    SearchOptions,
    catalogue,
    config_problem,
    config_sections,
    ready_problem,
    resolve_config,
    resolve_python,
)
from chemsuite.retro.routes import (
    OUTPUT_COLUMNS,
    compact_routes,
    leaf_stock_counts,
    sort_routes,
    summary_columns,
    tree_from_steps,
)
from chemsuite.retro.runner import RetroError, Worker, search, unique_targets

from . import recorded
from .fake import CONFIG_YML, fake_interpreter, recorded_environment, requests_logged, write_config

posix_only = pytest.mark.skipif(
    sys.platform == "win32", reason="the fake interpreter is a shell script"
)


# -- the config -----------------------------------------------------------------------------
def test_sections_lists_policies_filters_and_stocks() -> None:
    found = config_sections(CONFIG_YML)
    assert found["expansion"] == ("uspto", "ringbreaker")
    assert found["filter"] == ("uspto",)
    assert found["stock"] == ("zinc", "bench")


def test_sections_ignores_comments_and_nested_settings() -> None:
    body = (
        "# a comment\n"
        "expansion:\n"
        "  uspto:\n"
        "    type: template-based\n"
        "    model: uspto_model.onnx\n"
        "    cutoff_number: 50\n"
        "stock:\n"
        "  zinc: zinc_stock.hdf5\n"
        "search:\n"
        "  time_limit: 60\n"
    )
    found = config_sections(body)
    assert found["expansion"] == ("uspto",)
    assert found["stock"] == ("zinc",)
    assert found["filter"] == ()


def test_a_local_config_is_accepted() -> None:
    assert config_problem(CONFIG_YML) == ""


def test_a_remote_model_server_is_refused() -> None:
    body = CONFIG_YML.replace(
        "  uspto:\n    - uspto_model.onnx", "  uspto:\n    use_remote_models: True\n    model: x"
    )
    assert "use_remote_models" in config_problem(body)


def test_a_stock_database_on_another_host_is_refused() -> None:
    body = "stock:\n  zinc:\n    type: mongodb\n    host: chem-db.example.com\n"
    problem = config_problem(body)
    assert "chem-db.example.com" in problem


def test_a_stock_database_on_this_machine_is_allowed() -> None:
    body = (
        "expansion:\n  uspto:\n    - m\n    - t\n"
        "stock:\n  zinc:\n    type: mongodb\n    host: localhost\n"
    )
    assert config_problem(body) == ""


def test_a_host_outside_the_stock_section_does_not_make_it_remote() -> None:
    body = CONFIG_YML + "properties:\n  host: files.example.com\n"
    assert config_problem(body) == ""


def test_something_that_is_not_a_config_is_refused() -> None:
    assert "not an AiZynthFinder config" in config_problem("hello: world\n")


def test_catalogue_reads_the_file(tmp_path) -> None:
    found = catalogue(write_config(tmp_path))
    assert (found.policies, found.filters, found.stocks) == (
        ("uspto", "ringbreaker"),
        ("uspto",),
        ("zinc", "bench"),
    )
    assert found.problem == ""


def test_catalogue_without_a_file() -> None:
    assert catalogue(None).problem


def test_the_paths_come_from_settings(tmp_path) -> None:
    config = write_config(tmp_path)
    python = tmp_path / "python3"
    python.write_text("", encoding="utf-8")
    settings = Settings(
        tool_paths={"aizynthfinder_python": str(python), "aizynthfinder_config": str(config)}
    )
    assert resolve_python(settings) == str(python)
    assert resolve_config(settings) == config
    assert ready_problem(settings) == ""


def test_the_interpreter_is_not_looked_up_on_path() -> None:
    assert resolve_python(Settings()) == ""
    assert "not configured" in ready_problem(Settings())


def test_chemsuites_own_interpreter_is_refused(tmp_path) -> None:
    settings = Settings(
        tool_paths={
            "aizynthfinder_python": sys.executable,
            "aizynthfinder_config": str(write_config(tmp_path)),
        }
    )
    assert "its own environment" in ready_problem(settings)


def test_a_missing_config_is_reported(tmp_path) -> None:
    python = tmp_path / "python3"
    python.write_text("", encoding="utf-8")
    settings = Settings(tool_paths={"aizynthfinder_python": str(python)})
    assert "config.yml" in ready_problem(settings)


def test_options_survive_a_round_trip() -> None:
    options = SearchOptions(
        policies=("uspto",),
        stocks=("zinc",),
        algorithm="aizynthfinder.search.retrostar.search_tree.SearchTree",
        exploration=2.0,
        time_limit=30,
        max_routes=10,
        cluster=True,
    )
    assert SearchOptions.from_json(options.to_json()) == options


def test_damaged_options_fall_back_to_the_defaults() -> None:
    options = SearchOptions.from_json(
        {"time_limit": "soon", "algorithm": "made up", "policies": "uspto", "exploration": None}
    )
    assert options == SearchOptions()


def test_the_request_carries_every_setting() -> None:
    request = SearchOptions(policies=("uspto",), exploration=2.5, max_routes=7).request()
    assert request["policies"] == ["uspto"]
    assert request["search"]["algorithm_config"] == {"C": 2.5}
    assert request["post_processing"] == {"max_routes": 7}


# -- the parsers ----------------------------------------------------------------------------
def test_statistics_become_table_cells() -> None:
    cells = summary_columns(recorded.statistics())
    assert cells["AiZynth Solved"] == "Yes"
    assert cells["AiZynth Score"] == "0.994"
    assert cells["AiZynth Steps"] == "1"
    assert cells["AiZynth Time (s)"] == "18.4"
    assert cells["AiZynth First Solution (s)"] == "2.1"
    assert cells["AiZynth Precursors in Stock"] == f"{recorded.ACETIC}, {recorded.AMINOPHENOL}"
    assert cells["AiZynth Precursors not in Stock"] == ""
    assert cells["AiZynth Policies"] == "uspto: 137"
    assert cells["AiZynth Stock"] == "zinc"
    assert set(OUTPUT_COLUMNS) <= set(cells)


def test_an_unreported_statistic_still_gets_a_column() -> None:
    cells = summary_columns({"is_solved": False, "my_own_score": 1.5})
    assert cells["AiZynth My own score"] == "1.5"
    assert cells["AiZynth Solved"] == "No"


def test_a_missing_statistic_leaves_a_blank_cell() -> None:
    cells = summary_columns({"is_solved": True})
    assert cells["AiZynth Score"] == ""
    assert cells["AiZynth Steps"] == ""


def test_routes_are_compacted_best_first() -> None:
    found = compact_routes(recorded.TREES, recorded.STOCK_INFO)
    assert [route["score"] for route in found] == ["0.994", "0.970", "0.552"]
    assert [route["solved"] for route in found] == [True, True, False]
    best = found[0]
    assert best["tree"]["smiles"] == recorded.TARGET
    assert best["tree"]["in_stock"] is False
    reactants = best["tree"]["reactions"][0]["reactants"]
    assert [item["smiles"] for item in reactants] == [recorded.ACETIC, recorded.AMINOPHENOL]
    assert all(item["in_stock"] for item in reactants)
    assert reactants[0]["stock"] == ["zinc"]
    assert best["scores"]["number of reactions"] == 1
    assert best["metadata"]["created_at_iteration"] == 3


def test_a_step_keeps_the_reaction_and_its_template_metadata() -> None:
    best = compact_routes(recorded.TREES, recorded.STOCK_INFO)[0]
    (step,) = best["steps"]
    assert step["product"] == recorded.TARGET
    assert step["reaction"].startswith(recorded.TARGET + ">>")
    assert step["metadata"]["policy_name"] == "uspto"
    assert step["metadata"]["classification"] == "2.1.1 Amide formation"


def test_steps_are_numbered_from_the_earliest_reaction() -> None:
    two_step = compact_routes([recorded.TWO_STEP_ROUTE], recorded.STOCK_INFO)[0]
    products = [step["product"] for step in two_step["steps"]]
    assert products == [recorded.AMINOPHENOL, recorded.TARGET]  # step 1 makes the amine


def test_an_unsolved_route_counts_its_leaves() -> None:
    unsolved = compact_routes([recorded.UNSOLVED_ROUTE], recorded.STOCK_INFO)[0]
    assert unsolved["solved"] is False
    assert leaf_stock_counts(unsolved) == (1, 2)


def test_a_tree_rebuilds_from_its_steps() -> None:
    original = compact_routes([recorded.TWO_STEP_ROUTE], recorded.STOCK_INFO)[0]
    rebuilt = tree_from_steps(original["steps"])
    assert rebuilt is not None
    assert rebuilt["smiles"] == recorded.TARGET
    reactants = rebuilt["reactions"][0]["reactants"]
    assert sorted(item["smiles"] for item in reactants) == sorted(
        [recorded.ACETIC, recorded.AMINOPHENOL]
    )
    amine = next(item for item in reactants if item["smiles"] == recorded.AMINOPHENOL)
    assert amine["reactions"][0]["reactants"][0]["smiles"] == recorded.NITROPHENOL


def test_nothing_to_parse_is_not_an_error() -> None:
    assert compact_routes(None) == []
    assert compact_routes([{"smiles": ""}]) == []
    assert sort_routes([]) == []
    assert tree_from_steps([]) is None
    assert leaf_stock_counts({}) == (0, 0)


def test_sorting_keeps_the_order_within_a_tie() -> None:
    routes = [
        {"score": "0.5", "solved": False, "steps": [], "tree": {"smiles": "A"}},
        {"score": "0.5", "solved": False, "steps": [], "tree": {"smiles": "B"}},
    ]
    assert [route["tree"]["smiles"] for route in sort_routes(routes)] == ["A", "B"]


def test_a_route_without_scores_still_parses() -> None:
    bare = {"type": "mol", "smiles": "CCO", "in_stock": True}
    (route,) = compact_routes([bare])
    assert route["score"] == ""
    assert route["steps"] == []
    assert route["solved"] is True


def test_route_values_are_json_safe() -> None:
    assert routes_module._plain({"a": (1, 2), "b": {3}}) in (
        {"a": [1, 2], "b": "{3}"},
        {"a": [1, 2], "b": str({3})},
    )


def test_targets_are_searched_once_each() -> None:
    assert unique_targets(["CCO", " CCO ", "", None, "CCN"]) == ["CCO", "CCN"]


# -- the worker -----------------------------------------------------------------------------
@pytest.fixture
def fake(tmp_path, monkeypatch):
    """A fake interpreter and config, loaded with the recorded output."""
    recorded_environment(
        monkeypatch,
        recorded.TREES,
        recorded.STATISTICS,
        recorded.STOCK_INFO,
        log=tmp_path / "requests.jsonl",
    )
    python = fake_interpreter(tmp_path)
    config = write_config(tmp_path)
    worker = Worker(str(python), config)
    yield worker, tmp_path / "requests.jsonl"
    worker.close()


@posix_only
def test_the_worker_reports_its_version(fake) -> None:
    worker, _log = fake
    assert worker.hello()["version"] == "4.4.1"
    assert worker.alive


@posix_only
def test_a_search_comes_back_as_columns_and_routes(fake) -> None:
    worker, _log = fake
    reply = worker.search_one(recorded.TARGET, SearchOptions(time_limit=5))
    assert reply.ok
    assert reply.target == recorded.TARGET
    assert reply.columns["AiZynth Solved"] == "Yes"
    assert len(reply.routes) == 3
    assert reply.routes[0]["tree"]["smiles"] == recorded.TARGET


@posix_only
def test_the_options_reach_aizynthfinders_own_settings(fake) -> None:
    worker, log = fake
    options = SearchOptions(
        policies=("ringbreaker",),
        filters=("uspto",),
        stocks=("zinc",),
        exploration=2.5,
        time_limit=7,
        iteration_limit=42,
        max_transforms=4,
        return_first=True,
        max_routes=9,
    )
    assert worker.search_one(recorded.TARGET, options).ok
    (entry,) = requests_logged(log)
    assert entry["policies"] == ["ringbreaker"]
    assert entry["stocks"] == ["zinc"]
    assert entry["time_limit"] == 7
    assert entry["iteration_limit"] == 42
    assert entry["max_transforms"] == 4
    assert entry["return_first"] is True
    assert entry["max_routes"] == 9
    assert entry["algorithm_config"]["C"] == 2.5
    assert entry["algorithm_config"]["default_prior"] == 0.5  # the config's own value is kept


@posix_only
def test_an_empty_choice_means_everything_in_the_config(fake) -> None:
    worker, log = fake
    assert worker.search_one(recorded.TARGET, SearchOptions(time_limit=5)).ok
    (entry,) = requests_logged(log)
    assert entry["policies"] == ["uspto", "ringbreaker"]
    assert entry["filters"] == ["uspto"]
    assert entry["stocks"] == ["zinc", "bench"]


@posix_only
def test_one_target_that_fails_does_not_end_the_job(fake) -> None:
    worker, _log = fake
    replies = [
        worker.search_one("FAIL", SearchOptions(time_limit=5)),
        worker.search_one(recorded.TARGET, SearchOptions(time_limit=5)),
    ]
    assert not replies[0].ok
    assert "unsanitizable" in replies[0].error
    assert replies[1].ok  # the worker is still up


@posix_only
def test_a_library_writing_to_stdout_does_not_corrupt_the_replies(fake) -> None:
    worker, _log = fake
    assert worker.search_one("NOISY", SearchOptions(time_limit=5)).ok
    assert worker.search_one(recorded.TARGET, SearchOptions(time_limit=5)).ok


@posix_only
def test_the_worker_dying_is_reported_with_its_output(fake) -> None:
    worker, _log = fake
    with pytest.raises(RetroError) as raised:
        worker.search_one("CRASH", SearchOptions(time_limit=5))
    assert "stopped" in str(raised.value)
    assert not worker.alive


@posix_only
def test_a_warm_worker_searches_several_targets_in_order(tmp_path, monkeypatch) -> None:
    recorded_environment(
        monkeypatch, recorded.TREES, recorded.STATISTICS, recorded.STOCK_INFO, tmp_path / "log"
    )
    python, config = fake_interpreter(tmp_path), write_config(tmp_path)
    seen: list[tuple[int, int]] = []
    replies = search(
        str(python),
        config,
        ["CCO", "CCN", "CCC"],
        SearchOptions(time_limit=5),
        progress=lambda done, total: seen.append((done, total)),
        reuse=False,
    )
    assert [reply.target for reply in replies] == ["CCO", "CCN", "CCC"]
    assert seen[0] == (0, 3) and seen[-1] == (3, 3)
    assert [entry["target"] for entry in requests_logged(tmp_path / "log")] == ["CCO", "CCN", "CCC"]


@posix_only
def test_cancelling_keeps_the_targets_already_searched(tmp_path, monkeypatch) -> None:
    recorded_environment(monkeypatch, recorded.TREES, recorded.STATISTICS, recorded.STOCK_INFO)
    python, config = fake_interpreter(tmp_path), write_config(tmp_path)
    cancel = CancelToken()
    started = threading.Event()

    def progress(done: int, _total: int) -> None:
        if done:
            started.set()

    def cancel_once_started() -> None:
        started.wait(60)
        cancel.cancel()

    threading.Thread(target=cancel_once_started, daemon=True).start()
    with pytest.raises(Cancelled) as raised:
        search(
            str(python),
            config,
            [recorded.TARGET, "SLOW", "CCO"],
            SearchOptions(time_limit=5),
            progress=progress,
            cancel=cancel,
            reuse=False,
        )
    partial = raised.value.partial
    assert [reply.target for reply in partial] == [recorded.TARGET]
    assert partial[0].columns["AiZynth Solved"] == "Yes"


@posix_only
def test_a_worker_that_cannot_start_is_reported(tmp_path) -> None:
    worker = Worker(str(tmp_path / "no-such-python"), write_config(tmp_path))
    with pytest.raises(RetroError, match="Could not start"):
        worker.hello()
