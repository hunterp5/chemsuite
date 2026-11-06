# chemsuite

A desktop medicinal chemistry workbench for working with compound tables: plotting, descriptors,
conformers and SAR analysis. Built with Python, PySide6 and RDKit. Licensed under GPLv3.

How to use it: [docs/USER_GUIDE.md](docs/USER_GUIDE.md). The design and roadmap are in
[docs/PLAN.md](docs/PLAN.md) and [docs/ROADMAP.md](docs/ROADMAP.md).

When you add a tool or change how something works, update the user guide in the same change.
Build tool dialogs from `chemsuite.app.dialogs` (header, Input / Method / Output sections, footer),
so every dialog is laid out the same way.

## Development setup

```bash
python3.12 -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

## Checks (the same ones CI runs)

```bash
ruff check chemsuite tests scripts
ruff format --check chemsuite tests scripts
python scripts/check_gpl_headers.py      # add --fix to insert missing headers
pytest                                   # unit, property-based and architecture tests
pytest -m bench -s                       # performance targets at 200k rows
```

## Layout

| Path | Contents |
|---|---|
| `chemsuite/core/` | Qt-free core: table, structures, files, jobs. Must never import Qt (enforced by tests). |
| `chemsuite/app/` | Qt shell: main window, table view, plots, panels. |
| `chemsuite/tools/` | One folder per menu feature. |
| `tests/` | pytest suite; `tests/core/test_bench.py` holds the performance gates. |
