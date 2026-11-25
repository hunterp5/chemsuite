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

"""MMP memory: main process + pool workers, during and after find_pairs/transforms.

usage: python mmp_mem.py SDF N [workers]   (workers 0 = inline)
"""

import gc
import sys
import time

sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from memlib import Sampler, counters


def main():
    print("start                ", counters())
    import numpy as np
    from rdkit import Chem, rdBase

    from chemsuite.core.jobs.parallel import default_workers, process_pool, stop_pool
    from chemsuite.core.mols import mol_to_blob
    from chemsuite.sar import mmp

    path, n = sys.argv[1], int(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else default_workers()
    cells, act = [], []
    with rdBase.BlockLogs():
        for m in Chem.SDMolSupplier(path):
            if m is None:
                continue
            cells.append(mol_to_blob(m))
            act.append(float(m.GetProp("LogP")) if m.HasProp("LogP") else 0.0)
            if len(cells) >= n:
                break
    act = np.array(act)
    print(
        f"loaded {len(cells)} blobs        ", counters(), "blob MB", sum(map(len, cells)) // 2**20
    )

    pool = process_pool(workers) if workers else None
    pids = (lambda: list((pool._processes or {}).keys())) if pool else (lambda: [])
    t = time.perf_counter()
    with Sampler(pids) as s:
        pairs = mmp.find_pairs(cells, act, executor=pool)
        t1 = time.perf_counter() - t
        main, w, per = s.now()
        print(f"after find_pairs {t1:.0f}s pairs={len(pairs)}  main={main} workers={w} {per}")
        tr = mmp.transforms(pairs)
        print(f"after transforms  n={len(tr)}  main={counters()[3]}")
    print(f"PEAK main={s.peak_main}  workers={s.peak_workers}  total={s.peak_total}")
    time.sleep(2)
    main, w, per = s.now()
    print(f"idle, results kept: main={main} workers={w} per-worker={per}")
    del pairs, tr
    gc.collect()
    time.sleep(1)
    print(
        f"idle, results freed: main={counters()[3]} workers={sum(counters(p)[3] for p in pids())}"
    )
    if pool:
        stop_pool(pool)


if __name__ == "__main__":
    main()
