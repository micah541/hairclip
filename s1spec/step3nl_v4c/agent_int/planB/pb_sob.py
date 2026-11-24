"""(d) Plan B scale-invariant Sobolev / pointwise constants (template agent_sob/sob_constants_P2.py), W = exp(-psi) on [10.25, 260].
 (I1_W)  ||v'||^2 <= ||v|| ||v''|| + c1W ||v||^2,     c1W = (1/2) sup_[a,R] (psi'^2 - psi'')_+      (v compactly supported in (a, R))
 (I2_W)  W v^2(x) <= Psi1(x) ||v||^2_[x,R] + 2 ||v||_[x,R] ||v'||_[x,R],   Psi1(x) = sup_[x,R] psi'
         => sup_{[l,r]} v^2 <= KW (Psi1 ||v||^2 + 2 ||v|| ||v'||),  KW = exp(psi(r))  (psi increasing; psi' > 0 certified)
 same for the level-2 weight W2 = W*(1+(x/60)^4) =: exp(-psi2):  c1W2, Psi1_2, KW2 (psi2 = psi - log m2, psi2' > 0 certified).
 Cutoff derivative sups:  zeta_out = S((260-x)/10): sup|zeta_out^(j)| = supS_j / 10^j ; zeta_in = S((x-10.25)/6): supS_j / 6^j  (supS from the exact critical points).
 Cells of width 2^-6 (all breakpoints 10.25, 16.25, 20, 40, 60, 250, 260 are dyadic)."""
import sys, json, time
sys.path.insert(0, '.')
from pb_common import *
from fractions import Fraction
H = 2.0 ** -6; A_W = 10.25; R_MAX = 260.0
def cells(lo, hi):
    n = int(round((hi - lo) / H)); assert abs(lo + n * H - hi) < 1e-12
    for i in range(n): yield lo + i * H, lo + (i + 1) * H
def up4(v, sig=4):
    v = up(v) if isinstance(v, arb) else v
    if v == 0: return 0.0
    e = math.floor(math.log10(abs(v))) - (sig - 1); return float(f"{math.ceil(v / 10**e) * 10**e:.{sig}g}")
def cell_data(level):
    out = []
    for lo, hi in cells(A_W, R_MAX):
        X = I(lo, hi); P = Pj_box(X, level)
        out.append((lo, hi, P[0], P[1]))
    return out
def run():
    t0 = time.time(); R = {}
    for level, tag in ((0, "W"), (1, "W2")):
        cd = cell_data(level)
        if level == 0: assert all(c[2].lower() > 0 for c in cd), "psi' > 0 fails"
        up_d1 = [max(up(c[2]), 0.0) if level else up(c[2]) for c in cd]; up_c1 = [up(c[2] * c[2] - c[3]) for c in cd]
        lo_d1 = min(dn(c[2]) for c in cd)
        W2lo = [(W_lo(c[1]) * NW(c[0])) if level else None for c in cd]   # lower bound of W2 = W*m2 on each cell (W decreasing, m2 increasing)
        idx = lambda x: int(round((x - A_W) / H))
        o = {"psi1_min_lower_bound": lo_d1, "c1": {}, "Psi1": {}, "KW": {}}
        for Rr in (40.0, 60.0, 120.0, 260.0):
            m = max(up_c1[: idx(Rr)]); o["c1"][f"R={Rr:g}"] = up4(max(m, 0.0) / 2) if m > 0 else 0.0
            o["Psi1"][f"R={Rr:g}"] = {f"from_{l:g}": up4(max(up_d1[idx(l): idx(Rr)])) for l in (10.25, 16.25, 20.0, 40.0, 60.0, 120.0, 250.0) if l < Rr}
        regs = [(10.25, 16.25), (16.25, 20.0), (20.0, 40.0), (40.0, 80.0), (80.0, 160.0), (160.0, 260.0), (250.0, 260.0), (10.25, 260.0)]
        for l, r in regs:
            if level == 0:
                pv = psi_value(r); o["KW"][f"[{l:g},{r:g}]"] = up4(pv.exp()); o["KW"][f"[{l:g},{r:g}]_psi(right)"] = up4(pv)
            else:
                mn = min(W2lo[idx(l): idx(r)], key=lambda t: t.lower()); o["KW"][f"[{l:g},{r:g}]"] = up4(1 / mn)
        # sup of the region-wise psi' (for Psi1 at region start): done above.  Also max over [l, r] of psi' directly:
        o["note"] = "level 1: psi2 = psi - log(1+(x/60)^4); psi2' < 0 for x > ~60 (W2 INCREASES): Psi1 uses (psi2')_+, KW = 1/min W2 on the region" if level else ""
        o["sup_psi1_region"] = {f"[{l:g},{r:g}]": up4(max(up_d1[idx(l): idx(r)])) for l, r in regs}
        R[tag] = o
    # cutoff derivative sups
    import sympy as sp
    t = sp.symbols('t'); S = t**4 * (35 - 84*t + 70*t**2 - 20*t**3)
    sup_S = {}
    for j in range(1, 5):
        d = sp.diff(S, t, j); crit = [sp.Rational(0), sp.Rational(1)] + [r for r in sp.real_roots(sp.Poly(sp.diff(d, t), t)) if 0 <= r <= 1] if j < 5 else []
        sup_S[j] = max(abs(sp.N(d.subs(t, r), 40)) for r in crit)
    R["cutoff"] = {"supS_j": {j: float(v) for j, v in sup_S.items()},
                   "zeta_out_sup_j(w=10)": {j: up4(float(sup_S[j]) / 10**j * (1 + 1e-12)) for j in sup_S},
                   "zeta_in_sup_j(w=6)": {j: up4(float(sup_S[j]) / 6**j * (1 + 1e-12)) for j in sup_S}}
    R["W_values"] = {"W(10.25)": up4(W_up(10.25)), "W(16.25)": up4(W_up(16.25)), "W(250)": up4(W_up(250.0)), "W(260)": up4(W_up(260.0)),
                     "W(260)/W(10.25)": up4(W_up(260.0) / W_lo(10.25)), "W(260)/W(16.25)": up4(W_up(260.0) / W_lo(16.25)), "psi(260)": up4(psi_value(260.0)),
                     "integral_W_250_260_upper": None}
    n = 2000; tot = arb(0)
    for i in range(n):
        a = 250.0 + 10.0 * i / n; tot += arb(10.0 / n) * W_up(a)
    R["W_values"]["integral_W_250_260_upper"] = up4(tot)
    print(f"[{time.time()-t0:.0f}s]"); return R
if __name__ == "__main__":
    R = run(); json.dump(R, open('planB_sob.json', 'w'), indent=1); print(json.dumps(R, indent=1))
