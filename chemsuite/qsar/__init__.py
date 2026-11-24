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
"""Auto QSAR/QSPR: chemsuite's own model engine, on numpy, scipy and scikit-learn.

Qt-free, and it never imports the app: a tool hands it structure cells, numbers and a
:class:`~.recipe.Recipe`, and gets a fitted model back.

The design (``docs/QSAR_PLAN.md``):

- a large feature pool — descriptors, fingerprints and the user's own columns;
- feature selection and hyperparameter search **inside** the cross-validation folds, so the
  outer folds' score is an honest estimate rather than the fit's own optimism;
- conformal prediction intervals and an applicability-domain distance for every prediction;
- **fitted with scikit-learn, applied in numpy.** A fitted estimator is exported to flat
  arrays, and the archive the app reads carries no pickle (``allow_pickle=False``), the
  pattern :mod:`chemsuite.qsar.trees` and ``chemsuite/ionisation/unipka`` already follow.
"""
