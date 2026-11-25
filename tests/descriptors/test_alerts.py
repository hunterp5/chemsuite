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

"""Structural alerts: the Qt-free matching core."""

from __future__ import annotations

import pytest

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.descriptors.alerts import (
    ALERT_SETS,
    COMMON,
    alert_name,
    alerts_of,
    compute_alerts,
    count_column,
    names_column,
)

CATECHOL = "Oc1ccccc1O"
RHODANINE = "S=C1SC(=Cc2ccccc2)C(=O)N1"


def _mol(smiles: str):
    from rdkit import Chem

    return Chem.MolFromSmiles(smiles)


def test_every_set_has_a_catalog_and_alerts() -> None:
    from chemsuite.descriptors.alerts import catalog

    assert set(COMMON) <= set(ALERT_SETS)
    for key, spec in ALERT_SETS.items():
        assert catalog(key).GetNumEntries() > 0
        assert spec.title and spec.help
    titles = [spec.title for spec in ALERT_SETS.values()]
    assert len(set(titles)) == len(titles)  # a column name per set


def test_unknown_set_is_refused_before_any_work() -> None:
    with pytest.raises(KeyError):
        compute_alerts([CATECHOL], ["PAINS", "NOPE"])


def test_alert_names_drop_the_pains_count() -> None:
    assert alert_name("catechol_A(92)") == "catechol A"
    assert alerts_of(_mol(CATECHOL), "PAINS") == ["catechol A"]


def test_explicit_hydrogens_do_not_hide_an_alert() -> None:
    """Some sets' patterns count hydrogens, so a structure carrying them as atoms (one read
    from a file or made 3D) would miss: the core matches it without them."""
    from rdkit import Chem

    plain = _mol(CATECHOL)
    with_h = Chem.AddHs(plain)
    assert alerts_of(with_h, "NIH") == alerts_of(plain, "NIH") == ["ortho hydroquinone"]


def test_columns_count_and_name_the_alerts() -> None:
    cells = [CATECHOL, RHODANINE, "CCO", "", None, "not a structure"]
    columns = compute_alerts(cells, ["PAINS", "BRENK"])
    assert set(columns) == {
        count_column("PAINS"),
        names_column("PAINS"),
        count_column("Brenk"),
        names_column("Brenk"),
    }
    counts = columns[count_column("Brenk")]
    names = columns[names_column("Brenk")]
    shown = [counts.display(i) for i in range(len(cells))]
    assert shown == ["1", "2", "0", "", "", ""]  # a blank or unreadable cell: no answer
    assert names.display(0) == "catechol"
    assert names.display(1) == "Michael acceptor 1; Thiocarbonyl group"
    assert names.display(2) == "" and names.display(3) == ""


def test_a_cancel_keeps_the_rows_already_matched() -> None:
    token = CancelToken()
    cells = [CATECHOL] * 4

    def progress(done: int, _total: int) -> None:
        token.cancel()

    with pytest.raises(Cancelled) as stopped:
        compute_alerts(cells, ["PAINS"], progress=progress, cancel=token, chunk_rows=2)
    partial = stopped.value.partial
    assert len(partial[count_column("PAINS")]) == 2
