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

"""Turning AiZynthFinder's output into table cells and into routes the panel can draw.

AiZynthFinder reports a route as alternating molecule and reaction nodes under ``children``.
That shape is compacted here into one dict per route::

    {"score": "0.993", "solved": True,
     "steps": [{"product": …, "reaction": …, "reactants": [{"smiles", "in_stock"}…],
                "metadata": {…}}…],
     "tree": {"smiles": …, "in_stock": False, "stock": ["zinc"],
              "reactions": [{"smiles": …, "metadata": {…}, "reactants": [ …the same shape… ]}]},
     "scores": {…}, "metadata": {…}}

``tree`` is what the Routes panel draws (precursors to the left of their product); ``steps``
is the same reactions flat and numbered, for the step table and the sketcher. Both are plain
JSON, so a panel's routes save with the document.
"""

from __future__ import annotations

COLUMN_SOLVED = "AiZynth Solved"
COLUMN_SCORE = "AiZynth Score"
COLUMN_STEPS = "AiZynth Steps"
COLUMN_ROUTES = "AiZynth Routes"
COLUMN_SOLVED_ROUTES = "AiZynth Solved Routes"
COLUMN_TIME = "AiZynth Time (s)"
COLUMN_IN_STOCK = "AiZynth Precursors in Stock"
COLUMN_NOT_IN_STOCK = "AiZynth Precursors not in Stock"
COLUMN_FIRST_TIME = "AiZynth First Solution (s)"
COLUMN_FIRST_ITERATION = "AiZynth First Solution Iteration"
COLUMN_NODES = "AiZynth Nodes"
COLUMN_PRECURSORS = "AiZynth Precursors"
COLUMN_PRECURSORS_IN_STOCK = "AiZynth Precursors in Stock Count"
COLUMN_POLICIES = "AiZynth Policies"
COLUMN_STOCK = "AiZynth Stock"

#: the columns a search writes, in order
OUTPUT_COLUMNS: tuple[str, ...] = (
    COLUMN_SOLVED,
    COLUMN_SCORE,
    COLUMN_STEPS,
    COLUMN_ROUTES,
    COLUMN_SOLVED_ROUTES,
    COLUMN_TIME,
    COLUMN_IN_STOCK,
    COLUMN_NOT_IN_STOCK,
    COLUMN_FIRST_TIME,
    COLUMN_FIRST_ITERATION,
    COLUMN_NODES,
    COLUMN_PRECURSORS,
    COLUMN_PRECURSORS_IN_STOCK,
    COLUMN_POLICIES,
    COLUMN_STOCK,
)

#: statistic -> (column, how to format it)
_MAPPED: tuple[tuple[str, str, str], ...] = (
    ("is_solved", COLUMN_SOLVED, "yesno"),
    ("top_score", COLUMN_SCORE, "3dp"),
    ("number_of_steps", COLUMN_STEPS, "count"),
    ("number_of_routes", COLUMN_ROUTES, "count"),
    ("number_of_solved_routes", COLUMN_SOLVED_ROUTES, "count"),
    ("search_time", COLUMN_TIME, "1dp"),
    ("precursors_in_stock", COLUMN_IN_STOCK, "list"),
    ("precursors_not_in_stock", COLUMN_NOT_IN_STOCK, "list"),
    ("first_solution_time", COLUMN_FIRST_TIME, "1dp"),
    ("first_solution_iteration", COLUMN_FIRST_ITERATION, "count"),
    ("number_of_nodes", COLUMN_NODES, "count"),
    ("number_of_precursors", COLUMN_PRECURSORS, "count"),
    ("number_of_precursors_in_stock", COLUMN_PRECURSORS_IN_STOCK, "count"),
    ("policy_used_counts", COLUMN_POLICIES, "plain"),
    ("stock_info", COLUMN_STOCK, "stocks"),
)

#: statistics that are kept elsewhere (the row itself, the routes, the panel)
_SKIP = frozenset({"target", "trees", "stock_info", "distance_matrix", "cluster_labels"})


# -- table cells --------------------------------------------------------------------------------
def summary_columns(stats: dict) -> dict[str, str]:
    """One target's statistics as table cells.

    The named columns cover what AiZynthFinder always reports; anything else it reports is
    kept too, under its own column, so a config with extra scorers or a post-processing job
    is not silently dropped.
    """
    columns = {column: _format(stats.get(key), how) for key, column, how in _MAPPED}
    known = {key for key, _column, _how in _MAPPED}
    for key, value in stats.items():
        if key in known or key in _SKIP:
            continue
        name = f"AiZynth {str(key).replace('_', ' ').strip().capitalize()}"
        if name not in columns:
            columns[name] = _format(value, "plain")
    return columns


# -- routes -------------------------------------------------------------------------------------
def compact_routes(trees: object, stock_info: object = None) -> list[dict]:
    """The routes of one target, best first (see this module's docstring for the shape)."""
    if not isinstance(trees, list):
        return []
    stocks = stock_info if isinstance(stock_info, dict) else {}
    routes: list[dict] = []
    for tree in trees:
        if not isinstance(tree, dict) or not str(tree.get("smiles") or "").strip():
            continue
        scores = tree.get("scores") if isinstance(tree.get("scores"), dict) else {}
        route = {
            "score": _score_text(scores),
            "solved": _is_solved(tree),
            "steps": reaction_steps(tree),
            "tree": _compact_molecule(tree, stocks),
        }
        if scores:
            route["scores"] = {str(key): _plain(value) for key, value in scores.items()}
        metadata = tree.get("metadata")
        if isinstance(metadata, dict) and metadata:
            route["metadata"] = {str(key): _plain(value) for key, value in metadata.items()}
        routes.append(route)
    return sort_routes(routes)


def sort_routes(routes: list[dict]) -> list[dict]:
    """Solved routes first, then by score, keeping AiZynthFinder's order within a tie."""
    ordered = list(enumerate(routes))
    ordered.sort(key=lambda pair: (_sort_key(pair[1]), -pair[0]), reverse=True)
    return [route for _index, route in ordered]


def route_tree(route: dict) -> dict | None:
    """The molecule tree of *route*, rebuilt from its steps if it has none stored."""
    tree = route.get("tree")
    if isinstance(tree, dict) and str(tree.get("smiles") or "").strip():
        return tree
    return tree_from_steps(list(route.get("steps") or []))


def tree_from_steps(steps: list[dict]) -> dict | None:
    """Rebuild a molecule tree from flat steps; each step is used once."""
    unused = [dict(step) for step in steps if isinstance(step, dict)]
    if not unused:
        return None

    def take(product: str) -> dict | None:
        for index, step in enumerate(unused):
            if str(step.get("product") or "").strip() == product:
                return unused.pop(index)
        return None

    def build(smiles: str, in_stock: bool) -> dict:
        step = take(smiles)
        node = {"smiles": smiles, "in_stock": in_stock, "reactions": []}
        if step is None:
            return node
        reactants = [
            build(str(item.get("smiles") or "").strip(), bool(item.get("in_stock")))
            for item in step.get("reactants") or []
            if isinstance(item, dict)
        ]
        reaction: dict = {"reactants": reactants}
        if str(step.get("reaction") or "").strip():
            reaction["smiles"] = str(step["reaction"]).strip()
        if isinstance(step.get("metadata"), dict) and step["metadata"]:
            reaction["metadata"] = dict(step["metadata"])
        node["reactions"] = [reaction]
        return node

    return build(root, False) if (root := _target_of(unused)) else None


def _target_of(steps: list[dict]) -> str:
    """The product that nothing else consumes: the target the steps lead to.

    The steps are numbered from the earliest reaction, so the target is the last step's
    product; a damaged list without one falls back to the first product.
    """
    consumed = {
        str(item.get("smiles") or "").strip()
        for step in steps
        for item in step.get("reactants") or []
        if isinstance(item, dict)
    }
    products = [str(step.get("product") or "").strip() for step in steps]
    for product in reversed(products):
        if product and product not in consumed:
            return product
    return next((product for product in products if product), "")


def leaf_stock_counts(route: dict) -> tuple[int, int]:
    """``(in stock, leaves)`` for *route*: how much of it can simply be bought."""
    tree = route_tree(route)
    if tree is None:
        return 0, 0
    flags = _leaf_flags(tree)
    return sum(1 for flag in flags if flag), len(flags)


def reaction_steps(molecule: dict) -> list[dict]:
    """Every reaction under *molecule* (an AiZynthFinder node), flat and in route order.

    Step 1 is the earliest reaction and the last step makes the target, which is how a chemist
    reads a scheme; AiZynthFinder's tree runs the other way, so the walk is reversed here.
    """
    steps: list[dict] = []
    _collect_steps(molecule, steps)
    steps.reverse()
    return steps


# -- the walks ----------------------------------------------------------------------------------
def _compact_molecule(node: dict, stocks: dict) -> dict:
    reactions = []
    for child in node.get("children") or []:
        if not isinstance(child, dict) or not _is_reaction(child):
            continue
        reaction: dict = {
            "reactants": [
                _compact_molecule(item, stocks)
                for item in child.get("children") or []
                if isinstance(item, dict)
            ]
        }
        if str(child.get("smiles") or "").strip():
            reaction["smiles"] = str(child["smiles"]).strip()
        metadata = child.get("metadata")
        if isinstance(metadata, dict) and metadata:
            reaction["metadata"] = {str(key): _plain(value) for key, value in metadata.items()}
        reactions.append(reaction)
    smiles = str(node.get("smiles") or "").strip()
    compact = {"smiles": smiles, "in_stock": bool(node.get("in_stock")), "reactions": reactions}
    names = _stock_names(smiles, stocks)
    if names:
        compact["stock"] = names
    return compact


def _collect_steps(node: dict, steps: list[dict]) -> None:
    for child in node.get("children") or []:
        if not isinstance(child, dict):
            continue
        if not _is_reaction(child):
            _collect_steps(child, steps)
            continue
        reactants = [item for item in child.get("children") or [] if isinstance(item, dict)]
        step = {
            "product": str(node.get("smiles") or "").strip(),
            "reaction": str(child.get("smiles") or "").strip(),
            "reactants": [
                {
                    "smiles": str(item.get("smiles") or "").strip(),
                    "in_stock": bool(item.get("in_stock")),
                }
                for item in reactants
            ],
        }
        metadata = child.get("metadata")
        if isinstance(metadata, dict) and metadata:
            step["metadata"] = {str(key): _plain(value) for key, value in metadata.items()}
        steps.append(step)
        for item in reactants:
            _collect_steps(item, steps)


def _is_reaction(node: dict) -> bool:
    return node.get("type") == "reaction" or bool(node.get("is_reaction"))


def _is_solved(node: dict) -> bool:
    """A route is solved when every leaf of it is in stock."""
    leaves: list[bool] = []

    def walk(molecule: dict) -> None:
        reactions = [
            child
            for child in molecule.get("children") or []
            if isinstance(child, dict) and _is_reaction(child)
        ]
        if not reactions:
            leaves.append(bool(molecule.get("in_stock")))
            return
        for reaction in reactions:
            for item in reaction.get("children") or []:
                if isinstance(item, dict):
                    walk(item)

    walk(node)
    return bool(leaves) and all(leaves)


def _leaf_flags(node: dict) -> list[bool]:
    reactions = [item for item in node.get("reactions") or [] if isinstance(item, dict)]
    if not reactions:
        return [bool(node.get("in_stock"))]
    flags: list[bool] = []
    for reaction in reactions:
        for item in reaction.get("reactants") or []:
            if isinstance(item, dict):
                flags.extend(_leaf_flags(item))
    return flags


def _stock_names(smiles: str, stocks: dict) -> list[str]:
    raw = stocks.get(smiles) if smiles else None
    if isinstance(raw, str):
        return [part.strip() for part in raw.split(",") if part.strip()]
    if isinstance(raw, (list, tuple)):
        return [str(item).strip() for item in raw if str(item).strip()]
    return []


# -- formatting ---------------------------------------------------------------------------------
def _sort_key(route: dict) -> tuple[int, float]:
    try:
        score = float(route.get("score") or 0.0)
    except (TypeError, ValueError):
        score = 0.0
    return (1 if route.get("solved") else 0), score


def _score_text(scores: dict) -> str:
    for key, value in scores.items():
        if str(key).lower() == "state score":
            return _format(value, "3dp")
    for value in scores.values():
        text = _format(value, "3dp")
        if text:
            return text
    return ""


def _format(value: object, how: str) -> str:
    """One cell. A missing statistic becomes a blank cell, not a misleading zero."""
    if value is None or value == "" or (isinstance(value, (list, tuple, dict)) and not value):
        return ""
    if how == "yesno":
        return _yesno(value)
    if how == "count":
        return _count(value)
    if how in ("1dp", "3dp"):
        return _number(value, 1 if how == "1dp" else 3)
    if how == "list":
        return _list(value)
    if how == "stocks":
        return _stocks(value)
    return _plain_text(value)


def _yesno(value: object) -> str:
    if isinstance(value, str):
        text = value.strip().lower()
        if text in ("true", "yes", "1"):
            return "Yes"
        if text in ("false", "no", "0"):
            return "No"
        return value.strip()
    if isinstance(value, (list, tuple)):  # a multi-objective search reports one flag per reward
        return ", ".join(_yesno(item) for item in value)
    return "Yes" if value else "No"


def _count(value: object) -> str:
    if isinstance(value, (list, tuple)):
        return ", ".join(_count(item) for item in value)
    try:
        return f"{int(float(value)):d}"  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return _plain_text(value)


def _number(value: object, digits: int) -> str:
    if isinstance(value, (list, tuple)):
        return ", ".join(_number(item, digits) for item in value)
    try:
        return f"{float(value):.{digits}f}"  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return _plain_text(value)


def _list(value: object) -> str:
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item).strip() for item in value if str(item).strip())
    return str(value).strip()


def _stocks(value: object) -> str:
    """Which stocks held the precursors, once each: the per-molecule detail is on the cards."""
    if not isinstance(value, dict):
        return _plain_text(value)
    names: list[str] = []
    for raw in value.values():
        for name in _stock_names_of(raw):
            if name not in names:
                names.append(name)
    return ", ".join(names)


def _stock_names_of(raw: object) -> list[str]:
    if isinstance(raw, str):
        return [part.strip() for part in raw.split(",") if part.strip()]
    if isinstance(raw, (list, tuple)):
        return [str(item).strip() for item in raw if str(item).strip()]
    return []


def _plain_text(value: object) -> str:
    if isinstance(value, dict):
        return ", ".join(f"{key}: {_plain_text(item)}" for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return ", ".join(_plain_text(item) for item in value)
    if isinstance(value, float):
        return f"{value:.3f}".rstrip("0").rstrip(".")
    if isinstance(value, bool):
        return "Yes" if value else "No"
    return str(value).strip()


def _plain(value: object) -> object:
    """JSON-safe copy of a metadata or score value."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, (str, bool, int, float)) or value is None:
        return value
    return str(value)
