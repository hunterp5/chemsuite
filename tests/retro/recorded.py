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

"""AiZynthFinder output in its own shape, for tests and for the fake interpreter.

The nesting, the key names and the metadata fields are AiZynthFinder 4.4.1's: a route is
alternating ``mol``/``reaction`` nodes under ``children``, the root carries ``scores`` and
``metadata`` (``dict_with_extra``), and ``extract_statistics`` reports the flat statistics.

The chemistry is paracetamol from 4-aminophenol: two routes, one solved (both precursors in
the stock) and one not (the acid chloride is not), plus a two-step route for the tree walk.
"""

from __future__ import annotations

TARGET = "CC(=O)Nc1ccc(O)cc1"
AMINOPHENOL = "Nc1ccc(O)cc1"
ACETIC = "CC(=O)O"
ACETYL_CHLORIDE = "CC(=O)Cl"
NITROPHENOL = "Oc1ccc([N+](=O)[O-])cc1"


def _mol(smiles: str, in_stock: bool, children: list | None = None) -> dict:
    node: dict = {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": in_stock,
        "hide": False,
    }
    if children:
        node["children"] = children
    return node


def _reaction(smiles: str, metadata: dict, children: list) -> dict:
    return {
        "type": "reaction",
        "smiles": smiles,
        "is_reaction": True,
        "hide": False,
        "metadata": metadata,
        "children": children,
    }


_AMIDE_METADATA = {
    "policy_probability": 0.9349,
    "policy_probability_rank": 0,
    "policy_name": "uspto",
    "template_code": 12345,
    "template": "([#8:3]-[C:2](=[O:1]))>>Cl-[C:2](=[O:1]).[#8:3]",
    "classification": "2.1.1 Amide formation",
    "library_occurrence": 6185,
}
_NITRO_METADATA = {
    "policy_probability": 0.4117,
    "policy_probability_rank": 1,
    "policy_name": "uspto",
    "template_code": 2244,
    "template": "([N:1])>>[N+:1](=O)[O-]",
    "classification": "7.1.1 Nitro to amine reduction",
    "library_occurrence": 2051,
}

#: one solved route: paracetamol from the amine and acetic acid, both in the ZINC stock
SOLVED_ROUTE = _mol(
    TARGET,
    False,
    [
        _reaction(
            f"{TARGET}>>{ACETIC}.{AMINOPHENOL}",
            _AMIDE_METADATA,
            [_mol(ACETIC, True), _mol(AMINOPHENOL, True)],
        )
    ],
)
SOLVED_ROUTE["scores"] = {
    "state score": 0.9938,
    "number of reactions": 1,
    "number of pre-cursors": 2,
    "number of pre-cursors in stock": 2,
}
SOLVED_ROUTE["metadata"] = {"created_at_iteration": 3, "is_solved": True}

#: an unsolved route: the acid chloride is not in the stock
UNSOLVED_ROUTE = _mol(
    TARGET,
    False,
    [
        _reaction(
            f"{TARGET}>>{ACETYL_CHLORIDE}.{AMINOPHENOL}",
            _AMIDE_METADATA | {"policy_probability": 0.0551, "policy_probability_rank": 4},
            [_mol(ACETYL_CHLORIDE, False), _mol(AMINOPHENOL, True)],
        )
    ],
)
UNSOLVED_ROUTE["scores"] = {"state score": 0.5520, "number of reactions": 1}
UNSOLVED_ROUTE["metadata"] = {"created_at_iteration": 11, "is_solved": False}

#: two steps: the amine comes from nitrophenol, which is in the stock
TWO_STEP_ROUTE = _mol(
    TARGET,
    False,
    [
        _reaction(
            f"{TARGET}>>{ACETIC}.{AMINOPHENOL}",
            _AMIDE_METADATA,
            [
                _mol(ACETIC, True),
                _mol(
                    AMINOPHENOL,
                    False,
                    [
                        _reaction(
                            f"{AMINOPHENOL}>>{NITROPHENOL}",
                            _NITRO_METADATA,
                            [_mol(NITROPHENOL, True)],
                        )
                    ],
                ),
            ],
        )
    ],
)
TWO_STEP_ROUTE["scores"] = {"state score": 0.9701, "number of reactions": 2}

TREES = [SOLVED_ROUTE, TWO_STEP_ROUTE, UNSOLVED_ROUTE]

STOCK_INFO = {
    ACETIC: ["zinc"],
    AMINOPHENOL: ["zinc"],
    NITROPHENOL: ["zinc"],
    ACETYL_CHLORIDE: [],
}

STATISTICS = {
    "target": TARGET,
    "search_time": 18.42,
    "first_solution_time": 2.1077,
    "first_solution_iteration": 3,
    "number_of_nodes": 137,
    "max_transforms": 4,
    "max_children": 25,
    "number_of_routes": 3,
    "number_of_solved_routes": 2,
    "top_score": 0.99375,
    "is_solved": True,
    "number_of_steps": 1,
    "number_of_precursors": 2,
    "number_of_precursors_in_stock": 2,
    "precursors_in_stock": [ACETIC, AMINOPHENOL],
    "precursors_not_in_stock": [],
    "precursors_availability": "zinc;zinc",
    "policy_used_counts": {"uspto": 137},
    "profiling": {"expansion_calls": 137, "reactants_generations": 3402},
}


def statistics(target: str = TARGET) -> dict:
    """The statistics record for *target*, with the stock information AiZynthFinder adds."""
    return dict(STATISTICS) | {"target": target, "stock_info": dict(STOCK_INFO)}
