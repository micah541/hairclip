"""Build + run all integrations (no set reorganization) for generators p, C, S, w and tip bumps e1, e3.
Outputs out8Z<g>.txt and out8Q<gj><gi>.txt."""
import sys, time
from concurrent.futures import ThreadPoolExecutor
import gen_hs4  # noqa: patches gen_hs3 with the adjoint generators wa, wb
import gen_hs3 as G3
BASE = ["wa", "wb", "w"]; EB = ["e1", "e3"]
JOBS = ([("Z", g) for g in BASE + EB] + [("Q", gj, gi) for gj in BASE for gi in BASE]
        + [("Q", ea, gi) for ea in EB for gi in BASE + EB])
def name(job): return job[0] + "".join(job[1:])
def one(job, prefix="8"):
    nm = name(job); path = f"in{prefix}{nm}.txt"
    if job[0] == "Z": G3.build_Z(path, job[1])
    else: G3.build_Q(path, job[1], job[2])
    t0 = time.time(); import os; ev = {"REORG": "1e300"}; ev.update({k: os.environ[k] for k in ("BSUB", "HC") if k in os.environ}); out = G3.run(path, env=ev, exe=os.environ.get("EXE", "./hsrun3"))
    open(f"out{prefix}{nm}.txt", "w").write(out)
    return nm, "RESULT OK" in out, time.time() - t0
if __name__ == "__main__":
    prefix = sys.argv[1] if len(sys.argv) > 1 else "10"
    jobs = JOBS if len(sys.argv) <= 2 else [j for j in JOBS if name(j) in sys.argv[2:]]
    with ThreadPoolExecutor(2) as ex:
        res = list(ex.map(lambda j: one(j, prefix), jobs))
    print(" ".join(f"{t}:{'OK' if ok else 'FAIL'}({dt:.0f}s)" for t, ok, dt in res))
    if not all(ok for _, ok, _ in res): sys.exit(1)
