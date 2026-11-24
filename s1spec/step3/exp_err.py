"""Diagnostic (non-rigorous variants): which uncertainty drives the W enclosure width?"""
import sys, os, time
import gen_hs2 as G
from flint import arb
variant = sys.argv[1]
orig_params = G.params
orig_pp = G.profile_part
def params_pt():
    ps = orig_params()
    out = []
    for nm, lo, hi in ps:
        if nm == "stip": m = (lo + hi) / 2; out.append((nm, m, m))
        else: out.append((nm, lo, hi))
    return out
def pp_pt():
    tau1, xlo, xhi, Ct, r0lo, r0hi = orig_pp()
    return tau1, xlo, xhi, Ct, [0.0]*4, [0.0]*4
if "stip" in variant: G.params = params_pt
if "prof" in variant: G.profile_part = pp_pt
env = {}
if "hc" in variant: env["HC"] = "0.02"
env["SETDEBUG"] = "1"
path = f"inX_{variant}.txt"
t0 = time.time()
G.build_P(path, dsub=10.0)
out = G.run(path, env=env)
open(f"outX_{variant}.txt", "w").write(out)
for l in out.splitlines():
    if l.startswith("SET sigma=8 ") or l.startswith("SET sigma=20 "): print(variant, l)
tip, mon = G.parse(out)
if tip is None: print(variant, out[-500:]); sys.exit()
W = tip[4:10]
print(variant, "W13 tip:", W[4], "rel", (W[4][1]-W[4][0])/(W[4][1]+W[4][0]), "W01:", W[0], "rel", (W[0][1]-W[0][0])/abs(W[0][1]+W[0][0]), f"({time.time()-t0:.0f}s)")
