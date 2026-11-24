import sys, os, time
from concurrent.futures import ThreadPoolExecutor
from flint import arb
import gen_hs3 as G3
from gram_post6 import DSet, ratio, na
def runcfg(cfg):
    tag, hstep, hc, order = cfg
    env = {}
    if hc: env["HC"] = hc
    os.environ["HSTEP"] = hstep
    nm = f"tune_{tag}_{hstep}_{hc}_{order}"
    if tag.startswith("Z"): G3.build_Z(nm + ".in", tag[1:])
    else: G3.build_Q(nm + ".in", tag[1], tag[2])
    t0 = time.time(); out = G3.run(nm + ".in", order=order, env=env); dt = time.time() - t0
    open(nm + ".out", "w").write(out)
    return cfg, dt, nm + ".out", "RESULT OK" in out
cfgs = []
for hstep, hc, order in [("0.05", "", 20), ("0.05", "0.02", 20), ("0.025", "0.02", 20), ("0.05", "", 30), ("0.025", "0.01", 25)]:
    for tag in ["Zw", "Qww"]:
        cfgs.append((tag, hstep, hc, order))
res = {}
with ThreadPoolExecutor(10) as ex:
    for cfg, dt, path, ok in ex.map(runcfg, cfgs):
        res[cfg] = (dt, path, ok)
for hstep, hc, order in [("0.05", "", 20), ("0.05", "0.02", 20), ("0.025", "0.02", 20), ("0.05", "", 30), ("0.025", "0.01", 25)]:
    dz, pz, okz = res[("Zw", hstep, hc, order)]; dq, pq, okq = res[("Qww", hstep, hc, order)]
    if not (okz and okq): print(hstep, hc, order, "FAILED", okz, okq); continue
    S = DSet(pz); o = 5 + na(["w"]); r1 = ratio(S, o + 7, o + 4); r3 = ratio(S, o + 9, o + 4)
    Q = DSet(pq); o2 = 5 + na(["w"])
    g = Q.affine({o2 + 1: r3, o2 + 3: -r1, o2 + 4: arb(1)}) * arb(2) ** G3.QSHIFT
    print(f"HSTEP={hstep} HC={hc or 'auto'} order={order}: rho1 rad {float(r1.rad()):.2e} rho3 rad {float(r3.rad()):.2e}  Gww = {float(g.mid()):+.6e} +- {float(g.rad()):.2e}   ({dz:.0f}s, {dq:.0f}s)")
