"""RIGOROUS (ball arithmetic) positivity for the intermediate derivative energy (lemma L-int).
Weight W_I = exp(-phi_I) with phi_I' = (1-s) fp + s fr, fp = 2P x/(1+x^2), fr = (1/3) x (1+x^2)^(-1/3),
s = S((x - a2)/(b2 - a2)) clamped to [0,1], S(t) = t^4(35 - 84t + 70t^2 - 20t^3) (C^3).
Operator on v = w':  L1 = -d^4 - k^2 d^2 - (7km + q/4) d - 8m^2 - (5/2) k p   [(Lw)' = L1 w' + a0' w].
For v in C_c^oo((x_a, oo)):  Re<L1 v, v>_W <= -1/2 ||v''||_W^2 - int W mu1 v^2, with mu1 = -A^2/2 - (W A)''/(2W) - V,
A = 2(P1^2 - P2) + k^2,  V = c0 - (c1' - P1 c1)/2 + (W c2)''/(2W) - W''''/(2W),  P_j = phi^(j).
mu1 is a polynomial in (P1..P4, k, m, p, q); P_j are enclosed piecewise (Leibniz for the blend)."""
import sys, json, sympy as sp
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
ctx.prec = 100
PX, A2, B2, XA = 3, 30.25, 38.25, 12.0
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
def piece(xlo, xhi):
    X = arb((xlo + xhi) / 2, (xhi - xlo) / 2)
    sig_lo = max(float((arb(xlo) + STIP).lower()), 8.0)
    kb, mb, pb, qs = H.profile_bounds(arb(sig_lo))
    K_, M_, Pp = arb(0, float(kb.upper())), arb(0, float(mb.upper())), arb(0, float(pb.upper()))
    Q_ = X + STIP + arb(0, float(qs.upper()))
    region = "poly" if xhi <= A2 else ("rho" if xlo >= B2 else "blend")
    P = Pj(X, region)
    return MU(*P, K_, M_, Pp, Q_)
grid = [XA + 0.005 * i for i in range(int(round((A2 - XA) / 0.005)) + 1)]
grid += [A2 + 0.005 * i for i in range(1, int(round((B2 - A2) / 0.005)) + 1)]
grid += [B2 + 0.02 * i for i in range(1, 1501)]
worst = (1e9, None)
for a, b in zip(grid[:-1], grid[1:]):
    lo = float(piece(a, b).lower())
    if lo < worst[0]: worst = (lo, a)
print(f"mu1 >= {worst[0]:.4f} on [{XA}, {grid[-1]:.2f}]  (worst piece at x = {worst[1]:.3f})", flush=True)
# x >= 68.25: rho region; mu1 / x^{4/3} >= c via the u-chart bound: use monotone crude bound on [68.25, oo) by pieces in u
c_tail = None
u0 = arb(grid[-1]) ** (arb(-1) / 3)
N = 400; lows = []
for j in range(N):
    ulo, uhi = float(u0.upper()) * j / N, float(u0.upper()) * (j + 1) / N
    if j == 0: ulo = 1e-12
    # x in [uhi^-3, ulo^-3]; evaluate mu1/x^(4/3) with X as a ball (crude but fine for large x)
    xl, xh = uhi ** -3, ulo ** -3
    if xh > 1e12: xh = 1e12
    X = arb((xl + xh) / 2, (xh - xl) / 2)
    kb, mb, pb, qs = H.profile_bounds(arb(xl) + STIP.lower())
    P = Pj(X, "rho")
    v = MU(*P, arb(0, float(kb.upper())), arb(0, float(mb.upper())), arb(0, float(pb.upper())), X + STIP + arb(0, float(qs.upper()))) / X ** (arb(4) / 3)
    lows.append(float(v.lower()))
print("tail pieces (crude, x up to 1e12): min mu1/x^(4/3) =", min(lows))
json.dump({"mu1_min": worst[0], "x_worst": worst[1], "range": [XA, grid[-1]], "P": PX, "a2": A2, "b2": B2,
           "tail_min_mu1_over_x43": min(lows)}, open('lint_cert.json', 'w'), indent=1)
