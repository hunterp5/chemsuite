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

"""Move modules and packages with ``git mv`` and rewrite every reference to them.

    python scripts/move_modules.py MOVES [--dry-run]

*MOVES* is a text file of ``old new`` lines (repo-relative files or folders; ``#`` starts a
comment). Then, in ``chemsuite``, ``tests`` and ``scripts``:

- ``from ... import`` statements, absolute or relative, are rewritten to the new places. A name
  that is a module (``from chemsuite.app import theme``) follows its own module, so one import
  may become several. Relative imports stay relative while the target is in the same top-level
  package; a renamed module keeps its old name in the importing file (``import x as old``).
- Dotted names (``"chemsuite.app.jobs.Job"``, ``import chemsuite.app.jobs``) and paths
  (``chemsuite/app/main_window.py``) are rewritten anywhere in the code, in ``pyproject.toml``
  and in the Markdown docs.

Run ``ruff check --fix`` and ``ruff format`` on the touched files afterwards: rewritten imports
are written one per line and left for ruff to sort and wrap.
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve().relative_to(ROOT)  # its examples are not references
CODE_ROOTS = ("chemsuite", "tests", "scripts")
TEXT_FILES = ("pyproject.toml", "CLAUDE.md", "HANDOFF.md", "README.md")
DOTTED = re.compile(r"(?<![\w.])chemsuite(?:\.[A-Za-z_]\w*)+")
SLASHED = re.compile(r"(?<![\w/])chemsuite(?:/[A-Za-z_]\w*)+(?:\.py)?")


def module_of(path: Path) -> str:
    parts = list(path.with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def read_moves(path: Path) -> list[tuple[Path, Path]]:
    moves = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            old, new = line.split()
            moves.append((Path(old), Path(new)))
    return moves


def python_files() -> list[Path]:
    files = []
    for root in CODE_ROOTS:
        files += [
            p.relative_to(ROOT) for p in (ROOT / root).rglob("*.py") if "__pycache__" not in p.parts
        ]
    return sorted(f for f in files if f != SELF)


def file_moves(moves: list[tuple[Path, Path]], files: list[Path]) -> dict[Path, Path]:
    """Each moved Python file, old path to new path."""
    out: dict[Path, Path] = {}
    for old, new in moves:
        if (ROOT / old).is_dir():
            for f in files:
                if f.is_relative_to(old):
                    out[f] = new / f.relative_to(old)
        elif (ROOT / old).is_file():
            if old.suffix == ".py":
                out[old] = new
        else:
            sys.exit(f"nothing to move at {old}")
    return out


class Remapper:
    def __init__(self, moved: dict[Path, Path], files: list[Path]) -> None:
        self.known = {module_of(f) for f in files}
        self.table = {module_of(old): module_of(new) for old, new in moved.items()}

    def __call__(self, name: str) -> str:
        parts = name.split(".")
        for i in range(len(parts), 0, -1):
            new = self.table.get(".".join(parts[:i]))
            if new is not None:
                return ".".join([new, *parts[i:]])
        return name


def resolve(node: ast.ImportFrom, package: str) -> str:
    if not node.level:
        return node.module or ""
    base = package.split(".")
    base = base[: len(base) - node.level + 1]
    return ".".join([*base, *([node.module] if node.module else [])])


def relative(target: str, package: str) -> str | None:
    """*target* written relative to *package*, or None if they share no top-level package."""
    t, p = target.split("."), package.split(".")
    common = 0
    while common < min(len(t), len(p)) and t[common] == p[common]:
        common += 1
    if common < 2:
        return None
    return "." * (len(p) - common + 1) + ".".join(t[common:])


def rewrite_imports(
    source: str, old_path: Path, new_path: Path, remap: Remapper
) -> dict[int, tuple[int, list[str]]]:
    """Replacement lines for each import that changes: first line index -> (last, lines)."""
    old_module, new_module = module_of(old_path), module_of(new_path)
    is_package = old_path.name == "__init__.py"
    old_package = old_module if is_package else old_module.rpartition(".")[0]
    new_package = new_module if is_package else new_module.rpartition(".")[0]
    lines = source.splitlines()
    out: dict[int, tuple[int, list[str]]] = {}
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.ImportFrom):
            continue
        base = resolve(node, old_package)
        groups: dict[str, list[str]] = {}
        for alias in node.names:
            full = f"{base}.{alias.name}"
            if alias.name != "*" and full in remap.known:
                moved = remap(full)
                frm, _, name = moved.rpartition(".")
                asname = alias.asname or (alias.name if name != alias.name else None)
            else:
                frm, name, asname = remap(base), alias.name, alias.asname
            groups.setdefault(frm, []).append(f"{name} as {asname}" if asname else name)
        statements = []
        for frm, names in groups.items():
            where = relative(frm, new_package) if node.level else None
            statements.append(f"from {where or frm} import {', '.join(names)}")
        original_names = [f"{a.name} as {a.asname}" if a.asname else a.name for a in node.names]
        original_where = "." * node.level + (node.module or "")
        if statements == [f"from {original_where} import {', '.join(original_names)}"]:
            continue
        first, last = node.lineno - 1, (node.end_lineno or node.lineno) - 1
        indent = lines[first][: node.col_offset]
        tail = lines[last][node.end_col_offset :]
        if indent.strip() or (tail.strip() and not tail.strip().startswith("#")):
            raise SystemExit(f"{old_path}:{node.lineno}: import shares its line; move it by hand")
        rendered = [indent + s for s in statements]
        rendered[-1] += tail
        out[first] = (last, rendered)
    return out


def rewrite_text(text: str, remap: Remapper) -> str:
    text = DOTTED.sub(lambda m: remap(m.group(0)), text)

    def slashed(m: re.Match[str]) -> str:
        path = m.group(0)
        py = path.endswith(".py")
        dotted = remap(path.removesuffix(".py").replace("/", "."))
        return dotted.replace(".", "/") + (".py" if py else "")

    return SLASHED.sub(slashed, text)


def rewrite(source: str, old_path: Path, new_path: Path, remap: Remapper) -> str:
    imports = rewrite_imports(source, old_path, new_path, remap)
    lines = source.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        if i in imports:
            last, rendered = imports[i]
            out += rendered
            i = last + 1
        else:
            out.append(rewrite_text(lines[i], remap))
            i += 1
    return "\n".join(out) + ("\n" if source.endswith("\n") else "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("moves", type=Path)
    parser.add_argument("--dry-run", action="store_true", help="list what would change")
    args = parser.parse_args()

    moves = read_moves(args.moves)
    for _old, new in moves:
        if (ROOT / new).exists():  # git mv would move a folder *into* it
            sys.exit(f"{new} already exists")
    files = python_files()
    moved = file_moves(moves, files)
    remap = Remapper(moved, files)

    new_text: dict[Path, str] = {}
    for f in files:
        source = (ROOT / f).read_text(encoding="utf-8")
        result = rewrite(source, f, moved.get(f, f), remap)
        if result != source or f in moved:
            new_text[moved.get(f, f)] = result
    docs = [p.relative_to(ROOT) for p in (ROOT / "docs").rglob("*.md")]
    for f in [*docs, *(Path(t) for t in TEXT_FILES if (ROOT / t).exists())]:
        source = (ROOT / f).read_text(encoding="utf-8")
        result = rewrite_text(source, remap)
        if result != source:
            new_text[f] = result

    for old, new in moves:
        print(f"move  {old} -> {new}")
    for f in sorted(new_text):
        print(f"edit  {f}")
    if args.dry_run:
        return
    for old, new in moves:
        (ROOT / new).parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "mv", str(old), str(new)], cwd=ROOT, check=True)
    for f, text in new_text.items():
        (ROOT / f).write_text(text, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
