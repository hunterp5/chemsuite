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

"""Where AiZynthFinder is, what its ``config.yml`` offers, and the settings for one search.

Two paths come from ``Settings ▸ External Tools``: the interpreter that has AiZynthFinder
installed (``aizynthfinder_python``) and its ``config.yml`` (``aizynthfinder_config``). The
config is read here only to list the policy, filter and stock names for the dialog and to
refuse a config that would leave this machine; the worker reads the file itself with
AiZynthFinder's own loader.

The config is deliberately parsed by hand, without PyYAML: chemsuite does not depend on it,
and a key-and-indent scan is enough for names and for the local-only rules.
"""

from __future__ import annotations

import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

from chemsuite.core.settings import Settings

PYTHON_TOOL = "aizynthfinder_python"
CONFIG_TOOL = "aizynthfinder_config"

#: sections of ``config.yml`` whose keys name something the dialog can choose
SECTIONS = ("expansion", "filter", "stock")

_REMOTE_MODELS = re.compile(r"use_remote_models\s*:\s*(true|yes|on)\b", re.IGNORECASE)
_LOCAL_HOSTS = frozenset({"localhost", "127.0.0.1", "::1", "0.0.0.0"})

#: label -> ``search.algorithm``; anything but Monte Carlo is an AND/OR search class
ALGORITHMS: tuple[tuple[str, str], ...] = (
    ("Monte Carlo", "mcts"),
    ("Retro*", "aizynthfinder.search.retrostar.search_tree.SearchTree"),
    ("Breadth-first", "aizynthfinder.search.breadth_first.search_tree.SearchTree"),
    ("Proof number", "aizynthfinder.search.dfpn.search_tree.SearchTree"),
)
ALGORITHM_BY_VALUE = {value: label for label, value in ALGORITHMS}


@dataclass(frozen=True)
class SearchOptions:
    """One retrosynthesis run. Policy, filter and stock names come from the user's config.

    Empty *policies*, *filters* or *stocks* mean "everything in the config": AiZynthFinder
    itself would otherwise expand with the first policy only.
    """

    policies: tuple[str, ...] = ()
    filters: tuple[str, ...] = ()
    stocks: tuple[str, ...] = ()
    algorithm: str = "mcts"
    exploration: float = 1.4
    max_transforms: int = 6
    iteration_limit: int = 100
    time_limit: int = 120
    return_first: bool = False
    exclude_target_from_stock: bool = True
    max_routes: int = 25
    cluster: bool = False

    def request(self) -> dict:
        """The ``search`` request body the worker understands (without the target)."""
        return {
            "policies": list(self.policies),
            "filters": list(self.filters),
            "stocks": list(self.stocks),
            "search": {
                "algorithm": self.algorithm,
                "max_transforms": int(self.max_transforms),
                "iteration_limit": int(self.iteration_limit),
                "time_limit": int(self.time_limit),
                "return_first": bool(self.return_first),
                "exclude_target_from_stock": bool(self.exclude_target_from_stock),
                "algorithm_config": {"C": float(self.exploration)},
            },
            "post_processing": {"max_routes": int(self.max_routes)},
            "cluster": bool(self.cluster),
        }

    def to_json(self) -> dict:
        """JSON-friendly form, for the dialog's saved settings."""
        return {
            "policies": list(self.policies),
            "filters": list(self.filters),
            "stocks": list(self.stocks),
            "algorithm": self.algorithm,
            "exploration": float(self.exploration),
            "max_transforms": int(self.max_transforms),
            "iteration_limit": int(self.iteration_limit),
            "time_limit": int(self.time_limit),
            "return_first": bool(self.return_first),
            "exclude_target_from_stock": bool(self.exclude_target_from_stock),
            "max_routes": int(self.max_routes),
            "cluster": bool(self.cluster),
        }

    @classmethod
    def from_json(cls, raw: object) -> SearchOptions:
        """Read :meth:`to_json` back, falling back to the defaults for anything damaged."""
        saved = raw if isinstance(raw, dict) else {}
        blank = cls()
        algorithm = str(saved.get("algorithm") or blank.algorithm)
        return cls(
            policies=_names(saved.get("policies")),
            filters=_names(saved.get("filters")),
            stocks=_names(saved.get("stocks")),
            algorithm=algorithm if algorithm in ALGORITHM_BY_VALUE else blank.algorithm,
            exploration=_float(saved.get("exploration"), blank.exploration),
            max_transforms=_int(saved.get("max_transforms"), blank.max_transforms),
            iteration_limit=_int(saved.get("iteration_limit"), blank.iteration_limit),
            time_limit=_int(saved.get("time_limit"), blank.time_limit),
            return_first=bool(saved.get("return_first", blank.return_first)),
            exclude_target_from_stock=bool(
                saved.get("exclude_target_from_stock", blank.exclude_target_from_stock)
            ),
            max_routes=_int(saved.get("max_routes"), blank.max_routes),
            cluster=bool(saved.get("cluster", blank.cluster)),
        )


@dataclass(frozen=True)
class Catalogue:
    """The names a ``config.yml`` offers, in file order."""

    policies: tuple[str, ...] = ()
    filters: tuple[str, ...] = ()
    stocks: tuple[str, ...] = ()
    #: why the file cannot be used, if it cannot ("" when it is fine)
    problem: str = ""
    sections: dict[str, tuple[str, ...]] = field(default_factory=dict)


# -- the two paths ------------------------------------------------------------------------------
def resolve_python(settings: Settings, text: str = "") -> str:
    """The configured interpreter path (*text* wins while a dialog is open), or ``""``.

    There is no ``PATH`` fallback: the ``python3`` on ``PATH`` is almost never the environment
    AiZynthFinder is installed in, and silently searching with the wrong one wastes minutes.
    """
    raw = (text or "").strip() or settings.tool_paths.get(PYTHON_TOOL, "").strip()
    if not raw:
        return ""
    path = Path(raw).expanduser()
    if path.is_file():
        return str(path)
    found = shutil.which(raw)
    return found or ""


def resolve_config(settings: Settings, text: str = "") -> Path | None:
    """The configured ``config.yml``, or None when it is unset or missing."""
    raw = (text or "").strip() or settings.tool_paths.get(CONFIG_TOOL, "").strip()
    if not raw:
        return None
    path = Path(raw).expanduser()
    return path if path.is_file() else None


def is_app_python(python: str) -> bool:
    """True when *python* is the interpreter running chemsuite itself."""
    try:
        return Path(python).resolve() == Path(sys.executable).resolve()
    except OSError:
        return False


# -- what may run -------------------------------------------------------------------------------
def config_problem(text: str) -> str:
    """Why *text* (a ``config.yml`` body) must not be used, or ``""`` when it may.

    The search has to stay on this machine, so a remote policy server or a stock database on
    another host is refused rather than quietly accepted.
    """
    if _REMOTE_MODELS.search(text):
        return (
            "This config sets use_remote_models: true, which would send every policy "
            "prediction to a server. Leave it false so the ONNX models run on this machine."
        )
    host = _remote_stock_host(text)
    if host:
        return (
            f"This config uses a stock database on {host}. Use a local stock file "
            "(.hdf5, .csv, .bloom or a text file) so precursor lookups stay on this machine."
        )
    if not any(f"{section}:" in text for section in SECTIONS):
        return (
            "This file has no expansion, filter or stock section, so it is not an "
            "AiZynthFinder config.yml."
        )
    return ""


def ready_problem(settings: Settings, python_text: str = "", config_text: str = "") -> str:
    """Why a search cannot start, or ``""`` when it can. Starts no process."""
    python = resolve_python(settings, python_text)
    if not python:
        return (
            "AiZynthFinder is not configured. Set its Python interpreter and config.yml "
            "under Preferences ▸ External Tools."
        )
    if is_app_python(python):
        return (
            "That interpreter is the one running chemsuite. AiZynthFinder needs its own "
            "environment (Python 3.10–3.12 with RDKit 2023)."
        )
    config = resolve_config(settings, config_text)
    if config is None:
        return "Choose the AiZynthFinder config.yml under Preferences ▸ External Tools."
    try:
        body = config.read_text(encoding="utf-8")
    except OSError as err:
        return f"Could not read {config.name}: {err}"
    return config_problem(body)


# -- reading the config -------------------------------------------------------------------------
def catalogue(config: Path | None) -> Catalogue:
    """The policy, filter and stock names in *config*, and why it cannot be used."""
    if config is None:
        return Catalogue(problem="Choose the AiZynthFinder config.yml.")
    try:
        body = config.read_text(encoding="utf-8")
    except OSError as err:
        return Catalogue(problem=f"Could not read {config.name}: {err}")
    sections = config_sections(body)
    return Catalogue(
        policies=sections["expansion"],
        filters=sections["filter"],
        stocks=sections["stock"],
        problem=config_problem(body),
        sections=sections,
    )


def config_sections(text: str) -> dict[str, tuple[str, ...]]:
    """Top-level keys under each of :data:`SECTIONS`, in file order.

    A section's entries are its keys one level in, whether the value is a path, a list of
    paths or a block of settings::

        expansion:
          uspto:
            - uspto_model.onnx
            - uspto_templates.csv.gz
        stock:
          zinc: zinc_stock.hdf5
    """
    found: dict[str, list[str]] = {name: [] for name in SECTIONS}
    section: str | None = None
    indent: int | None = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0] if raw.lstrip().startswith("#") else raw
        stripped = line.strip()
        if not stripped:
            continue
        depth = len(line) - len(line.lstrip(" "))
        if depth == 0:
            key = stripped[:-1].strip() if stripped.endswith(":") else ""
            section = key if key in found else None
            indent = None
            continue
        if section is None or stripped.startswith("- ") or ":" not in stripped:
            continue
        if indent is None:
            indent = depth  # the first entry sets the section's own indent
        if depth != indent:
            continue
        name = stripped.split(":", 1)[0].strip().strip("\"'")
        if name and name not in found[section]:
            found[section].append(name)
    return {name: tuple(values) for name, values in found.items()}


def _remote_stock_host(text: str) -> str:
    """The host of a stock database that is not on this machine, or ``""``.

    Only a stock whose ``type`` names a database is checked: a ``host:`` line elsewhere (a
    model server that is not switched on, say) does not make the config remote.
    """
    sections = _stock_blocks(text)
    for block in sections:
        kind = _value_of(block, "type").lower()
        if "mongo" not in kind:
            continue
        host = _value_of(block, "host")
        bare = host.strip("\"'").split("://")[-1].rsplit("@", 1)[-1]
        bare = bare.strip("[]").split(":", 1)[0].strip("]")
        if bare and bare.lower() not in _LOCAL_HOSTS:
            return host
        if not bare:
            return "another host"
    return ""


def _stock_blocks(text: str) -> list[list[str]]:
    """The lines of each entry under ``stock:``."""
    blocks: list[list[str]] = []
    in_stock = False
    entry_indent: int | None = None
    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        depth = len(raw) - len(raw.lstrip(" "))
        if depth == 0:
            in_stock = stripped.rstrip(":").strip() == "stock" and stripped.endswith(":")
            entry_indent = None
            continue
        if not in_stock:
            continue
        if entry_indent is None:
            entry_indent = depth
        if depth == entry_indent:
            blocks.append([stripped])
        elif blocks:
            blocks[-1].append(stripped)
    return blocks


def _value_of(lines: list[str], key: str) -> str:
    for line in lines:
        head, _, tail = line.partition(":")
        if head.strip().strip("\"'").lower() == key:
            return tail.strip()
    return ""


def _names(value: object) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(str(item).strip() for item in value if str(item).strip())


def _int(value: object, default: int) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default


def _float(value: object, default: float) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default
