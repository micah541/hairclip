"""RIGOROUS (ball arithmetic) positivity for the intermediate derivative energy (lemma L-int).
Weight W_I = exp(-phi_I) with phi_I' = (1-s) fp + s fr, fp = 2P x/(1+x^2), fr = (1/3) x (1+x^2)^(-1/3),
s = S((x - a2)/(b2 - a2)) clamped to [0,1], S(t) = t^4(35 - 84t + 70t^2 - 20t^3) (C^3).
Operator on v = w':  L1 = -d^4 - k^2 d^2 - (7km + q/4) d - 8m^2 - (5/2) k p   [(Lw)' = L1 w' + a0' w].
For v in C_c^oo((x_a, oo)):  Re<L1 v, v>_W <= -1/2 ||v''||_W^2 - int W mu1 v^2, with mu1 = -A^2/2 - (W A)''/(2W) - V,
A = 2(P1^2 - P2) + k^2,  V = c0 - (c1' - P1 c1)/2 + (W c2)''/(2W) - W''''/(2W),  P_j = phi^(j).
mu1 is a polynomial in (P1..P4, k, m, p, q); P_j are enclosed piecewise (Leibniz for the blend)."""
import sys, json, time, argparse, numpy as np, sympy as sp
from flint import arb, ctx
_ap = argparse.ArgumentParser()
_ap.add_argument('--PX', default='3'); _ap.add_argument('--A2', default='30.25'); _ap.add_argument('--B2', default='38.25')
_ap.add_argument('--XA', default='12'); _ap.add_argument('--out', default=None)
_ARGS = _ap.parse_args(); _T0 = time.time()
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
ctx.prec = 100
PX, A2, B2, XA = sp.Rational(_ARGS.PX), sp.Rational(_ARGS.A2), sp.Rational(_ARGS.B2), float(_ARGS.XA)
P1, P2, P3, P4, P5, k, m, p, q = sp.symbols('P1 P2 P3 P4 P5 k m p q')
rules = {P1: P2, P2: P3, P3: P4, P4: P5, p: -k * q, q: 1 + k * p, k: m, m: p / 4}
def D(e): return sp.expand(sum(sp.diff(e, v) * r for v, r in rules.items()))
def rr2(f): return D(D(f)) - 2 * P1 * D(f) + (P1**2 - P2) * f           # (W f)''/W,  W = e^{-phi}
c2, c1, c0 = -k**2, -(7 * k * m + q / 4), -8 * m**2 - sp.Rational(5, 2) * k * p
A = 2 * (P1**2 - P2) - c2
V = c0 - (D(c1) - P1 * c1) / 2 + rr2(c2) / 2 - rr2(P1**2 - P2) / 2
mu1 = sp.expand(-A**2 / 2 - rr2(A) / 2 - V)
assert not mu1.has(P5)
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
MU = eval("lambda P1, P2, P3, P4, k, m, p, q: " + _Pr().doprint(mu1), {"arb": arb})
# derivatives of fp, fr (orders 0..3) and of the blend s (orders 0..3)
x = sp.symbols('x')
fpd = [sp.diff(2 * PX * x / (1 + x**2), x, j) for j in range(4)]
frd = [sp.diff(x * (1 + x**2) ** sp.Rational(-1, 3) / 3, x, j) for j in range(4)]
t = (x - A2) / (B2 - A2); Ssym = t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3)
sd = [sp.diff(Ssym, x, j) for j in range(4)]
class _Pr2(_Pr):
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp()
        if ex.is_Integer: return _Pr._print_Pow(self, e)
        return f"(({self._print(b)})**({self._print(sp.Rational(ex))}))"
    def _print_Float(self, e): return f"arb('{e}')"
lam = lambda e: eval("lambda x: " + _Pr2().doprint(e), {"arb": arb})
FP = [lam(e) for e in fpd]; FR = [lam(e) for e in frd]; SD = [lam(e) for e in sd]
from math import comb
def Pj(X, region):
    """(P1..P4) enclosures on the ball X"""
    if region == "poly": return [FP[j](X) for j in range(4)]
    if region == "rho": return [FR[j](X) for j in range(4)]
    s = [SD[j](X) for j in range(4)]; d = [FR[j](X) - FP[j](X) for j in range(4)]
    return [FP[j](X) + sum(comb(j, i) * s[i] * d[j - i] for i in range(j + 1)) for j in range(4)]
STIP = H.STIP
A2f, B2f = float(A2), float(B2)
# --- profile boxes for sigma < 8 (needed when x + STIP < 8, i.e. x < 10.2491): CAPD boxes mc_prof.txt
# columns: tau_lo tau_hi | p | q | k | m (each lo hi) | STIP lo hi ; x = 8 - tau - STIP (checked: q ~ 8 - tau, p' = -kq, m' = p/4)
_BX = np.loadtxt('/home/claude/s1spec/step3nl/mc_prof.txt')
_XC = float((arb(8) - STIP).upper()) + 1e-9
def _hull(lo, hi): return arb(lo).union(arb(hi))
def box_bounds(xlo, xhi):
    """hull of (k, m, p, q) over all CAPD boxes whose x-range meets [xlo, xhi] (x <= x_c ~ 10.2491)"""
    xl = 8 - _BX[:, 1] - float(STIP.upper()); xh = 8 - _BX[:, 0] - float(STIP.lower())
    sel = (xh >= xlo - 1e-7) & (xl <= xhi + 1e-7)
    assert sel.any()
    B = _BX[sel]
    g = lambda c: _hull(float(B[:, c].min()), float(B[:, c + 1].max()))
    return g(6), g(8), g(2), g(4)          # k, m, p, q
NANFLAG = [False]
def piece(xlo, xhi):
    X = arb((xlo + xhi) / 2, (xhi - xlo) / 2)
    region = "poly" if xhi <= A2f else ("rho" if xlo >= B2f else "blend")
    if xlo + float(STIP.lower()) < 8.0:
        K_, M_, Pp, Q_ = box_bounds(xlo, min(xhi, _XC))
        if xhi > _XC:                                   # straddling: union with the sigma >= 8 tail bound
            kb, mb, pb, qs = H.profile_bounds(arb(8))
            K_ = K_.union(arb(0, float(kb.upper()))); M_ = M_.union(arb(0, float(mb.upper()))); Pp = Pp.union(arb(0, float(pb.upper())))
            Q_ = Q_.union(X + STIP + arb(0, float(qs.upper())))
    else:
        sig_lo = max(float((arb(xlo) + STIP).lower()), 8.0)
        kb, mb, pb, qs = H.profile_bounds(arb(sig_lo))
        K_, M_, Pp = arb(0, float(kb.upper())), arb(0, float(mb.upper())), arb(0, float(pb.upper()))
        Q_ = X + STIP + arb(0, float(qs.upper()))
    P = Pj(X, region)
    r = MU(*P, K_, M_, Pp, Q_)
    lo = float(r.lower())
    if lo != lo: NANFLAG[0] = True
    return r, lo, region
grid = [XA + 0.005 * i for i in range(int(round((A2f - XA) / 0.005)) + 1)]
grid += [A2f + 0.005 * i for i in range(1, int(round((B2f - A2f) / 0.005)) + 1)]
grid += [B2f + 0.02 * i for i in range(1, 1501)]
reg = {}; WTOL = 0.01; DEPTH = 14; NPIECE = [0]
def cert(a, b, d=0):
    """adaptive bisection (dependency-problem control only; every final piece is a rigorous enclosure)"""
    r, lo, rg = piece(a, b)
    if lo != lo: NANFLAG[0] = True
    if lo == lo and float(r.upper()) - lo > WTOL and d < DEPTH:
        mid = (a + b) / 2; cert(a, mid, d + 1); cert(mid, b, d + 1); return
    NPIECE[0] += 1
    if lo != lo: reg[rg] = (float('nan'), a)
    elif rg not in reg or lo < reg[rg][0]: reg[rg] = (lo, a)
for a, b in zip(grid[:-1], grid[1:]): cert(a, b)
worst = min(reg.values(), key=lambda t: t[0])
print(f"PX={PX} A2={A2f} B2={B2f}: " + "  ".join(f"{k}: {v[0]:.4f}@{v[1]:.3f}" for k, v in reg.items()) + f"  NaN={NANFLAG[0]}", flush=True)
TAILNAN = False
# x >= grid[-1]: rho region crude u-chart check (as in the original; lint_tail.py is the actual tail certificate)
u0 = arb(grid[-1]) ** (arb(-1) / 3)
N = 400; lows = []
for j in range(N):
    ulo, uhi = float(u0.upper()) * j / N, float(u0.upper()) * (j + 1) / N
    if j == 0: ulo = 1e-12
    xl, xh = uhi ** -3, ulo ** -3
    if xh > 1e12: xh = 1e12
    X = arb((xl + xh) / 2, (xh - xl) / 2)
    kb, mb, pb, qs = H.profile_bounds(arb(xl) + STIP.lower())
    P = Pj(X, "rho")
    v = MU(*P, arb(0, float(kb.upper())), arb(0, float(mb.upper())), arb(0, float(pb.upper())), X + STIP + arb(0, float(qs.upper()))) / X ** (arb(4) / 3)
    f = float(v.lower())
    if f != f: TAILNAN = True
    lows.append(f)
tail_min = float('nan') if TAILNAN else min(lows)
print("crude tail check (NOT rigorous/NaN in original; superseded by lint_tail.py):", tail_min, " NaN=", TAILNAN)
res = {"PX": float(PX), "A2": A2f, "B2": B2f, "XA": XA, "regions": {k: {"min_mu1": v[0], "at_x": v[1]} for k, v in reg.items()},
       "mu1_min": worst[0], "x_worst": worst[1], "range": [XA, grid[-1]], "crude_tail_nan": TAILNAN,
       "nan": NANFLAG[0], "n_pieces": NPIECE[0], "runtime_s": time.time() - _T0}
json.dump(res, open(_ARGS.out or f"runs/cert_PX{_ARGS.PX}_XA{_ARGS.XA}_A2{_ARGS.A2}.json", 'w'), indent=1)
