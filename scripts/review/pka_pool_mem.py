"""pKa memory with the app's pool: main + every worker, peak and idle afterwards.

usage: python pka_pool_mem.py SDF N [workers]
"""

import sys
import time

sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from memlib import Sampler, counters


def main():
    from rdkit import Chem, rdBase

    from chemsuite.core.jobs.parallel import default_workers, process_pool, stop_pool
    from chemsuite.core.mols import mol_to_blob
    from chemsuite.ionisation.predict import predict_pka

    path, n = sys.argv[1], int(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else default_workers()
    cells = []
    with rdBase.BlockLogs():
        for m in Chem.SDMolSupplier(path):
            if m is not None:
                cells.append(mol_to_blob(m))
            if len(cells) >= n:
                break
    pool = process_pool(workers)
    pids = lambda: list((pool._processes or {}).keys())  # noqa: E731
    t = time.perf_counter()
    last = [0.0]

    def progress(done, total):
        if time.perf_counter() - last[0] > 5:
            last[0] = time.perf_counter()
            m, w, per = s.now()
            print(
                f"  {done}/{total} {last[0] - t:.0f}s main={m} workers={w} max1={max(per or [0])}"
            )

    with Sampler(pids) as s:
        res = predict_pka(cells, executor=pool, progress=progress)
    print(f"done {len(res)} rows in {time.perf_counter() - t:.0f}s, failed {res.failed}")
    print(f"PEAK main={s.peak_main} workers={s.peak_workers} total={s.peak_total}")
    peaks = [counters(p)[1] for p in pids()]
    print("per-worker peak working set", peaks)
    time.sleep(3)
    m, w, per = s.now()
    print(f"idle: main={m} workers={w} per-worker private={per}")
    stop_pool(pool)


if __name__ == "__main__":
    main()
