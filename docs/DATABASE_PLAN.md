# Plan: File ▸ Import from Database (read-only SQL queries) — approved 2026-11-12, not started

## Context

Round 3 of the menu work asked for SQL database support, and the user chose "Import a query":
run a SQL query against a database and open the rows as a new table. SQLite works out of the
box; PostgreSQL, MySQL/MariaDB, Oracle and SQL Server work through optional drivers. Nothing is
ever written to the database. ROADMAP §5 already says: "SQL through the standard library or an
optional driver", and network access is off under the machine policy. PLAN.md drops sqlalchemy,
so we talk DB-API 2.0 directly. CLAUDE.md requires plan mode for a big feature.

## Design

### Qt-free core: `chemsuite/io/database.py`

- `Connection` (frozen dataclass): `kind` (`sqlite`, `postgresql`, `mysql`, `oracle`,
  `sqlserver`), `path` (SQLite), `host`, `port`, `database`, `user`.
  - The password is passed separately and **never stored**, logged or audited.
  - `label()` gives the text shown in the dialog and the status bar, e.g.
    `chembl@db.example:5432`.
- `DRIVERS`: kind → (import name, pip extra, connect function).
  - `sqlite`: `sqlite3`, standard library.
  - `postgresql`: `psycopg` 3 (LGPL-3.0).
  - `mysql`: `PyMySQL` (MIT; MariaDB too).
  - `oracle`: `oracledb` thin mode (Apache-2.0/UPL; no Oracle client needed).
  - `sqlserver`: `pyodbc` (MIT; uses the installed "ODBC Driver 17/18 for SQL Server").

  `available(kind)` uses `importlib.util.find_spec`, so nothing is imported until used.
- **Read-only, three layers:**
  1. **Statement guard.** `check_query(sql)` strips comments, allows exactly one statement,
     and requires it to start with `SELECT`, `WITH` or `VALUES`. A `WITH` that contains
     `INSERT/UPDATE/DELETE/MERGE` is refused too. Each failure is a `QueryRefused` with a
     plain message.
  2. **Read-only session**, where the server supports one:
     - SQLite opens as `file:…?mode=ro` (URI);
     - PostgreSQL sets `conn.read_only = True`;
     - MySQL runs `START TRANSACTION READ ONLY`;
     - Oracle runs `SET TRANSACTION READ ONLY`;
     - SQL Server sets `ApplicationIntent=ReadOnly`. That only routes the connection, so the
       guard and layer 3 carry the weight here.
  3. **Never commit**: always `rollback()`, then `close()`, in a `finally`.
- `read_query(conn, sql, password=None, *, structure_column=AUTO, limit=None,
  chunk_rows=CHUNK_ROWS, cancel=None, progress=None) -> Iterator[Batch]`.
  - Headers come from `cursor.description`, made unique with `delimited.unique_headers`.
  - Rows are read with `fetchmany(chunk_rows)`, with `check(cancel)` between chunks.
  - Values become text the way CSV cells do:
    - `None` → blank;
    - `bytes` → the molblock or SMILES when they decode as UTF-8, else hex;
    - `Decimal` and dates → `str`.
  - Number and text columns are then typed by `column_from_values` exactly as for CSV.
  - Structure and reaction columns reuse `delimited.guess_structure_column`,
    `reaction_columns` and `_batch`. Molblock text in a column is detected like SDF text
    (a structure column of molblocks).
  - `limit` wraps nothing: the query is fetched up to `limit` rows and then stops.
- `list_tables(conn, password)`: the catalogue query for each kind (`sqlite_master`,
  `information_schema.tables`, `all_tables`), to fill the dialog's table picker.
- **Cancel on the server.** A cancelled job calls the driver's cancel from the job's cancel
  hook: `conn.cancel()` for psycopg and oracledb, `cursor.cancel()` for pyodbc,
  `sqlite3.Connection.interrupt()`. PyMySQL has none, so its query is stopped between chunks
  only.
- **Network policy.** Every kind but SQLite requires `core.download.network_allowed()`;
  otherwise `NetworkBlocked`.
- A `statement_timeout` (default 300 s), set where the driver supports it.

### Dependencies (`pyproject.toml`)

- New extras: `postgresql = ["psycopg[binary]>=3.2"]`, `mysql = ["PyMySQL>=1.1"]`,
  `oracle = ["oracledb>=2.0"]`, `sqlserver = ["pyodbc>=5.1"]`, and a `databases` extra that
  pulls all four.
- No base-install change. Licences are all GPL-compatible; recorded in the USER_GUIDE
  licence table.

### The dialog: `chemsuite/app/document/database_dialog.py` (kit `ToolDialog`, non-modal)

- **Input** section:
  - Database type combo. A type whose driver is missing stays listable; choosing it shows
    `problem("Install the PostgreSQL driver: pip install chemsuite[postgresql]")`.
  - SQLite: a file picker.
  - Servers: host, port (the kind's default), database, user and password
    (`QLineEdit.Password`).
  - A **Recent** combo for the last 10 connections, saved in user settings
    (`database.recent`) without passwords.
- **Query** section:
  - a monospace `QPlainTextEdit` for the SQL;
  - a **Table** picker (filled by `list_tables` on demand, from a ↻ button) that inserts
    `SELECT * FROM <table>`;
  - **Preview**: runs with `limit=100` in a background job and shows the rows in a small
    read-only `QTableWidget` under the editor.
- **More Options** (folded): Structure column (Auto, a named column, None) and Row limit (none
  by default).
- Primary button **Run Query**. The SQL guard's message shows as a ⚠ problem while typing.
  No helper text; tooltips carry explanations.
- **Running.** `jobs.run("Querying <label>…", work, opened)` in the queue lane. The rows open
  like `open_path`: in this window if it is empty and unmodified, else a new one; the
  document is untitled.
- The query and connection (no password) are kept as document extras (`"source"`), so a later
  "Refresh from Database" is possible. Not built now.
- `log.info` the label and row count, and `audit("import_database", kind=…, database=…,
  rows=…)`, with no SQL and no credentials.

### Menu

- File ▸ **Import from Database…**, after Open Recent. Command id `file.import_database`.

### Docs

- USER_GUIDE:
  - a "Databases" subsection under opening files;
  - the menu table row;
  - the read-only guarantee and its limits (SQL Server's routing caveat);
  - the extras and the network policy.
- ROADMAP §5: mark the SQL row as done.

## Files

- New: `chemsuite/io/database.py`, `chemsuite/app/document/database_dialog.py`,
  `tests/io/test_database.py`, `tests/app/document/test_database_dialog.py`.
- Changed: `chemsuite/app/document/main_window.py` (one menu action; it is at its line budget,
  so the dialog opens through a small function in the dialog module), `pyproject.toml`,
  `docs/USER_GUIDE.md`, `docs/ROADMAP.md`.
- Reused:
  - `delimited.guess_structure_column`, `reaction_columns`, `unique_headers` and `_batch`;
  - `load.TableBuilder` and `parse_pending_structures`;
  - `core.download.network_allowed` and `NetworkBlocked`;
  - the kit `ToolDialog`, `Section` and `problem()`;
  - `JobHost.run`;
  - the window-reuse logic in `open_path`, factored into one helper both use.

## Verification

- `tests/io/test_database.py` (SQLite in `tmp_path`, plus a fake DB-API module for the
  servers):
  - types, NULLs and blobs;
  - unique headers;
  - SMILES and molblock columns detected;
  - chunking and `limit`;
  - cancel between chunks and through `interrupt()`;
  - the guard refuses `INSERT`, `;`-joined statements, `WITH … DELETE` and `ATTACH`;
  - SQLite `mode=ro` refuses a write even past the guard;
  - the fake server sees `read_only`/`SET TRANSACTION READ ONLY` and `rollback`, never
    `commit`;
  - the network policy blocks server kinds;
  - a missing driver gives a clear message;
  - the password is never in a log record.
- `tests/app/document/test_database_dialog.py`:
  - a SQLite file → Preview → Run Query → a new table with the right rows and a structure
    column;
  - a missing driver shows a ⚠ problem and disables Run;
  - Recent connections remembered without the password.
- Also run `tests/test_architecture.py` and `tests/test_import_boundaries.py` (`io` stays
  Qt-free). Then `ruff check`/`ruff format` on the touched paths.
- By hand: query a ChEMBL SQLite release (`SELECT canonical_smiles, chembl_id FROM … LIMIT
  1000`) and see the structures draw.
