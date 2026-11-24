"""Plan B common helpers (R = 260).  Weight W = exp(-psi), psi' = 2 P(x) x/(1+x^2), P = 2 (x<=40), 2 - S((x-40)/20) on [40,60], 1 (x>=60);
psi(0) = 0, so W = (1+x^2)^-2 for x <= 40 (same normalisation as the P=2 design on [10.25,30.25]).  S = t^4(35-84t+70t^2-20t^3).
zeta: inner ramp S((x-10.25)/6) on [10.25,16.25]; outer ramp S((260-x)/10) on [250,260].
Level-2 weight factor  m2(x) = 1 + (x/60)^4  (effective phi2 = psi - log m2).
Boxes: sigma-boxes [sa,sb] (floats) with x = sigma - ST (tip enclosure), as in the earlier int_common.py."""
import sys, math, json, time
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
from flint import arb, acb, ctx, arb_poly
from math import comb
import sympy as sp
import h2res as HR
F, Hf, ST, STIP, dST = HR.F, HR.Hf, HR.ST, HR.STIP, HR.dST
RMAX = 260.0
UP, DN = math.inf, -math.inf
def up(v):
    assert not v.is_nan(), 'NaN ball'
    return math.nextafter(float(v.upper().mid()), UP) if v.upper().rad() == 0 else math.nextafter(float(v.upper().mid() + v.upper().rad()), UP)
def dn(v):
    assert not v.is_nan(), 'NaN ball'
    return math.nextafter(float(v.lower().mid()), DN)
def chk(v):
    assert not v.is_nan(), 'NaN ball'
    return v
def I(lo, hi): return arb(lo).union(arb(hi))
# ---------------------------------------------------------------- P(x) and its derivatives (Taylor-shifted S)
_tt = sp.symbols('tt'); _St = _tt**4 * (35 - 84*_tt + 70*_tt**2 - 20*_tt**3)
_SDP = [arb_poly([arb(int(c.p)) / int(c.q) for c in reversed(sp.Poly(sp.diff(_St, _tt, j), _tt).all_coeffs())]) for j in range(5)]
def S_derivs(X, c0, w, n=4):
    """S^(j)((x - c0)/w) as d^j/dx^j, j=0..n-1, on ball X (Taylor shift at the double midpoint, no cancellation blow-up)"""
    xm = arb(float(X.mid())); tm = (xm - arb(c0)) / w; h = (X - xm) / w
    with ctx.workprec(160):
        return [_SDP[j](arb_poly([tm, 1]))(h) / arb(w) ** j for j in range(n)]
def Pf_derivs(X, n=4):
    lo, hi = dn(X), up(X)
    if hi <= 40: return [arb(2)] + [arb(0)] * (n - 1)
    if lo >= 60: return [arb(1)] + [arb(0)] * (n - 1)
    sd = S_derivs(X, 40, 20, n); r = [-s for s in sd]; r[0] = 2 + r[0]; return r
x_ = sp.symbols('x')
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        n = int(ex); prod = "*".join([bb] * abs(n)); return f"({prod})" if n > 0 else f"(1/({prod}))"
_lam = lambda e: eval("lambda x: " + _Pr().doprint(e), {"arb": arb})
GD = [_lam(sp.diff(2 * x_ / (1 + x_**2), x_, j)) for j in range(5)]      # derivatives of g = 2x/(1+x^2)
LM2 = [_lam(sp.diff(sp.log(1 + (x_ / 60)**4), x_, j)) if j else None for j in range(5)]
_g = sp.log(1 + (x_ / 60)**4)
LM2 = [None] + [_lam(sp.diff(_g, x_, j)) for j in range(1, 5)]
def Pj_box(X, level=0):
    """[psi', psi'', psi''', psi''''] on ball X (level 1: for the level-2 weight  W*m2, i.e. psi - log m2)"""
    lo, hi = dn(X), up(X)
    cuts = [lo] + [c for c in (40.0, 60.0) if lo < c < hi] + [hi]
    out = None
    for a, b in zip(cuts[:-1], cuts[1:]):
        Xi = I(a, b); Pf = Pf_derivs(Xi, 4); G = [GD[j](Xi) for j in range(4)]
        P = [sum(comb(j, i) * Pf[i] * G[j - i] for i in range(j + 1)) for j in range(4)]
        out = P if out is None else [o.union(v) for o, v in zip(out, P)]
    if level == 1:
        L = [LM2[j](X) for j in range(1, 5)]; out = [o - l for o, l in zip(out, L)]
    for v in out: chk(v)
    return out
# ---------------------------------------------------------------- psi values and W bounds
def _integ(z, _an):
    zz = z * z; t = (z - 40) / 20; t2 = t * t; t3 = t2 * t; t4 = t2 * t2
    S = t4 * (35 - 84 * t + 70 * t2 - 20 * t3)
    return 2 * (2 - S) * z / (1 + zz)
_cache = {}
def psi_value(xv):
    if xv in _cache: return _cache[xv]
    x = arb(xv)
    with ctx.workprec(130):
        if xv <= 40: v = 2 * (1 + x * x).log()
        else:
            b40 = 2 * (arb(1) + arb(1600)).log()
            if xv <= 60: v = b40 + acb.integral(_integ, acb(40), acb(xv)).real
            else:
                b60 = b40 + acb.integral(_integ, acb(40), acb(60)).real
                v = b60 + (1 + x * x).log() - (arb(1) + arb(3600)).log()
    chk(v); _cache[xv] = v; return v
def W_up(xlo):
    x0 = math.nextafter(xlo, DN) if xlo > 0.5 else xlo
    return arb(up((-psi_value(x0)).exp()))
def W_lo(xhi):
    return arb(dn((-psi_value(math.nextafter(xhi, UP))).exp()))
# ---------------------------------------------------------------- boxes
def sig_of_x(x): return math.nextafter(float((arb(x) + ST).upper().mid()), UP)
class Box: pass
def make_box(sa, sb):
    b = Box(); b.sa, b.sb = sa, sb; sig = arb(sa).union(arb(sb)); b.sig = sig
    b.xa, b.xb = dn(arb(sa) - ST), up(arb(sb) - ST); b.xc = 0.5 * (b.xa + b.xb); b.X = I(b.xa, b.xb)
    kb, mb, pb, qs = HR.HS.profile_bounds(arb(sa)); b.kb, b.mb, b.pb, b.qs = kb, mb, pb, qs
    b.p, b.k, b.m = arb(0, pb.upper()), arb(0, kb.upper()), arb(0, mb.upper()); b.q = sig + arb(0, qs.upper())
    return b
def make_segments(bps, widths):
    pts = [sig_of_x(v) if t == 'x' else float(v) for t, v in bps]; out = []
    for (a, b), w in zip(zip(pts[:-1], pts[1:]), widths):
        n = max(1, math.ceil((b - a) / w - 1e-12)); edges = [a + (b - a) * i / n for i in range(n)] + [b]
        out += list(zip(edges[:-1], edges[1:]))
    return out
REG = [(10.25, 16.25), (16.25, 20.0), (20.0, 40.0), (40.0, 60.0), (60.0, 120.0), (120.0, 200.0), (200.0, 250.0), (250.0, 260.0)]
BPS = [('x', a) for a in (10.25, 16.25, 20.0, 40.0, 60.0, 120.0, 200.0, 250.0, 260.0)]
def region_of(xa, xb): return [i for i, (a, b) in enumerate(REG) if xb > a + 1e-6 and xa < b - 1e-6]
def reg_name(i): return f"[{REG[i][0]:g},{REG[i][1]:g}]"
def Sf(t):
    t = min(max(t, 0.0), 1.0); a = arb(t); return a * a * a * a * (35 - 84 * a + 70 * a * a - 20 * a * a * a)
def zeta_sup(xa, xb, R=260.0):
    zi = up(Sf((xb - 10.25) / 6)) if xb < 16.25 else 1.0
    zo = up(Sf((R - xa) / 10)) if xa > R - 10 else 1.0
    return min(1.0, zi) * min(1.0, zo)
def NW(x): return (1 + (x / 60) ** 4)
