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

"""Permeability prediction: three models chemsuite ships, applied in numpy (Qt-free).

Each endpoint is a model fitted by **chemsuite's own Auto QSAR engine**
(:mod:`chemsuite.qsar`) and shipped as a ``.csqm`` archive in :mod:`.models`. There is no new
modelling here and no new runtime dependency: :mod:`.predict` builds the engine's feature pool
over a column of structures and asks the model, which runs on numpy alone — the app never
loads scikit-learn or a pickle to predict. docs/PERMEABILITY_PLAN.md holds the design, the
licence chain and the measured accuracy.

Two endpoints, and the honest reading of each:

* **MDR1-MDCK efflux ratio** — how much more a compound is pumped out basolateral-to-apical
  than the other way. High means a P-glycoprotein substrate. The cleaner of the two
  (Biogen's MIT-licensed set) and the one that answers an efflux question Caco-2 Papp alone
  cannot.
* **Caco-2 permeability** — apparent permeability across a Caco-2 monolayer, log10 cm/s. What
  chemists ask for by name, and the noisiest: inter-laboratory scatter on this assay is
  0.3-0.5 log units, which is the floor any model is measured against.

**PAMPA was measured and held back** (2026-11-25): on the public NCATS set, which is binary and
85% permeable, the fitted model reached ROC-AUC 0.740 but balanced accuracy 0.570 and MCC
0.252 — it ranks a series but barely separates the impermeable minority, so per-compound calls
would need more caveat than they are worth. The curation and training scripts still cover it,
and the classification path here is still exercised by the tests, so it is one line to restore
when a continuous PAMPA set is available. docs/PERMEABILITY_PLAN.md records the numbers.

Every prediction travels with an interval (or a class set) and a distance to the training set,
because a permeability number with no error bar is the thing a chemist is right to distrust —
and because these models, like every ADMET model, are only worth reading inside the chemistry
they were fitted on.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Every regression endpoint's value is a base-10 logarithm, as the training data reports it.
LOG10 = "log10"


@dataclass(frozen=True)
class Endpoint:
    """One shipped model: what it predicts, what to call it, and which archive holds it."""

    id: str
    #: what the dialog calls it
    title: str
    #: the prefix of every column written for it ("Caco-2 Papp Predicted", "… Low", …)
    column: str
    #: "regression" or "classification"
    task: str
    #: the units a chemist reads the number in; "" for a classification
    unit: str
    #: the ``.csqm`` file in ``data/``
    file: str
    #: the tooltip: what the assay measures and how to read it. Dialogs show no helper text, so
    #: this is where an explanation is allowed to live.
    tip: str


ENDPOINTS: tuple[Endpoint, ...] = (
    Endpoint(
        "mdck_er",
        "MDR1-MDCK efflux ratio",
        "MDCK ER",
        "regression",
        "log10 (B-A/A-B)",
        "mdck_er.csqm",
        "How strongly P-glycoprotein pumps the compound back out. Above about 0.5 (a ratio of "
        "3) is usually read as a P-gp substrate.",
    ),
    Endpoint(
        "caco2",
        "Caco-2 permeability",
        "Caco-2 Papp",
        "regression",
        "log10 cm/s",
        "caco2.csqm",
        "Apparent permeability across a Caco-2 monolayer. Above about -5.0 is usually read as "
        "good absorption; the assay itself scatters by 0.3 to 0.5 log units between "
        "laboratories.",
    ),
)

#: The endpoints a dialog checks first. Both of them today; the tuple stays because a held-back
#: or weaker endpoint added later should not be on by default.
COMMON = ("mdck_er", "caco2")


def by_id(endpoint_id: str) -> Endpoint:
    """The endpoint with this id."""
    for endpoint in ENDPOINTS:
        if endpoint.id == endpoint_id:
            return endpoint
    raise KeyError(f"unknown permeability endpoint {endpoint_id!r}")
