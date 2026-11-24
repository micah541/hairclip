"""Build + run all 20 integrations (no set reorganization), outputs out7<tag>.txt."""
import sys, time, os
from concurrent.futures import ThreadPoolExecutor
import gen_hs3 as G3
TAGS = [f"Z{g}" for g in "pCSw"] + [f"Q{j}{i}" for j in "pCSw" for i in "pCSw"]
def one(tag, prefix="7", env_extra=None):
    path = f"in{prefix}{tag}.txt"
    if tag.startswith("Z"): G3.build_Z(path, tag[1:])
    else: G3.build_Q(path, tag[1], tag[2])
    env = {"REORG": "1e300"}; env.update(env_extra or {})
    t0 = time.time(); out = G3.run(path, env=env)
    open(f"out{prefix}{tag}.txt", "w").write(out)
    return tag, "RESULT OK" in out, time.time() - t0
if __name__ == "__main__":
    prefix = sys.argv[1] if len(sys.argv) > 1 else "7"
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda t: one(t, prefix), TAGS))
    print(" ".join(f"{t}:{'OK' if ok else 'FAIL'}({dt:.0f}s)" for t, ok, dt in res))
