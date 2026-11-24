"""R=60 design: T3 coefficients (d_tjl, A0, A1) with jet box |w'| <= 0.02 and |w|-box 0.02 on [16.25,26], 0.5 beyond 26; alpha <= 1e-3, C^4 jets.
Coefficients do not depend on the weight W (only on jets/profile), so they apply to the rho-blend W of the R=60 design as is.
d_tjl are the coefficients of the EXACT polynomial D(w'',w''',w'''') = sum d_tjl w''^t w'''^j w''''^l (monomials (0,0,1),(0,1,0),(1,0,0),(1,1,0),(2,0,0),(3,0,0) only, finite
Poly in sympy, coefficients functions of (U0,U1,u0'',u0''',u0'''',profile) only) -> no box on w'',w''',w''''."""
import sys, json
sys.argv = ['x', '1e-3', '1e-2', '1e-2', '1']
import pb_nlin as M
from pb_common import *
from c4fun import PSI, H2
SUBS = [(16.25, 20.0), (20.0, 26.0), (26.0, 30.25), (30.25, 40.0), (40.0, 60.0)]
bl = make_segments([('x', a) for a in (16.25, 20.0, 26.0, 30.25, 40.0, 60.0)], [0.025, 0.025, 0.025, 0.025, 0.05])
jets = [(make_box(sa, sb),) for sa, sb in bl]
J = []
for (b,) in jets: J.append((b, PSI.jets_any(b.X, b.xc, 5), H2.jets_any(b.X, b.xc, 5)))
KEYS = {"d001": (0, 0, 1), "d010": (0, 1, 0), "d100": (1, 0, 0), "d200": (2, 0, 0), "d110": (1, 1, 0), "d300": (3, 0, 0), "A0": "A0", "A1": "A1"}
def run(e0f, e1):
    out = {}
    for b, y, h in J:
        M.A = arb(0, 1e-3); M.E0b = arb(0, e0f(b.xa)); M.E1b = arb(0, e1)
        d = M.box_vals(y, h, b.p, b.q, b.k, b.m)
        for ri, (lo, hi) in enumerate(SUBS):
            if b.xb > lo + 1e-6 and b.xa < hi - 1e-6:
                s = out.setdefault(f"[{lo:g},{hi:g}]", {})
                for n, k in KEYS.items(): s[n] = max(s.get(n, 0.0), up(d[k]))
                s["sup|u0'|"] = max(s.get("sup|u0'|", 0.0), up(d["|u0'|"])); s["sup|u0''|"] = max(s.get("sup|u0''|", 0.0), up(d["|u0''|"]))
    return out
new = run(lambda x: 0.02 if x < 26 else 0.5, 0.02)
cfg = {"new(w0 0.02|0.5, w1 0.02)": new, "old box 1e-2/1e-2": run(lambda x: 1e-2, 1e-2), "w0=0, w1=0.02": run(lambda x: 0.0, 0.02), "w0=0.5 everywhere, w1=0.02": run(lambda x: 0.5, 0.02)}
def merge(o):   # [20,30.25] = max of [20,26], [26,30.25]
    m = {k: max(o["[20,26]"][k], o["[26,30.25]"][k]) for k in o["[20,26]"]}
    r = {k: v for k, v in o.items() if k not in ("[20,26]", "[26,30.25]")}; r["[20,30.25]"] = m
    return dict(sorted(r.items(), key=lambda kv: float(kv[0][1:].split(',')[0])))
res = {k: {"T3_per_region": merge(v)} for k, v in cfg.items()}
res["sub_regions_new"] = new
json.dump(res, open('r60_T3_boxw1_0.02.json', 'w'), indent=1)
for k in cfg:
    print(k)
    for r, v in res[k]["T3_per_region"].items(): print("  ", r, {a: float('%.4g' % b) for a, b in v.items()})
