"""[v2] Common helpers; weight parameter PX in {2, 3/2} (psi' = (1-s) 2 PX x/(1+x^2) + s (1/3) x (1+x^2)^(-1/3)).
Old docstring: Common rigorous helpers for the intermediate-region (weight W = exp(-psi)) certificates.
psi(0)=0, psi' = (1-s) 6x/(1+x^2) + s (1/3) x (1+x^2)^(-1/3), s = S((x-30.25)/8) clamped, S(t)=t^4(35-84t+70t^2-20t^3).
 => W = (1+x^2)^(-3) for x <= 30.25.   Boxes are sigma-boxes [sa, sb] (python floats, exact), x = sigma - ST (ST = tip enclosure)."""
import sys, math, json, time
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
from flint import arb, acb, ctx
from math import comb
import sympy as sp
import h2res as HR
F, Hf, ST, STIP, dST = HR.F, HR.Hf, HR.ST, HR.STIP, HR.dST
jets_any = HR.jets_any
A2, B2 = 30.25, 38.25
UP, DN = math.inf, -math.inf
def up(v):   # float upper bound of an arb (rigorous)
    assert not v.is_nan(), 'NaN ball'
    return math.nextafter(float(v.upper().mid()), UP) if v.upper().rad() == 0 else float(v.upper().mid() + v.upper().rad())
def dn(v):
    assert not v.is_nan(), 'NaN ball'
    return math.nextafter(float(v.lower().mid()), DN)
def chk(v):
    assert not v.is_nan(), 'NaN ball'
    return v
def I(lo, hi):
    """arb ball containing [lo, hi] (floats, exact endpoints)"""
    return arb(lo).union(arb(hi))
# ------------------------------------------------------------------ psi = integral of psi'
from fractions import Fraction
PXS = {"2": Fraction(2), "1.5": Fraction(3, 2)}
def _mkint(px):
    pn, pd = px.numerator, px.denominator
    def f(z, _an):
        qq = acb(1) + z * z
        q13 = (qq.log() * (-acb(1) / 3)).exp()
        A = 2 * acb(pn) / pd * z / qq; B = z * q13 / 3
        t = (z - acb("30.25")) / 8; t2 = t * t; t3 = t2 * t; t4 = t2 * t2
        s = t4 * (35 - 84 * t + 70 * t2 - 20 * t3)
        return (1 - s) * A + s * B
    return f
_INTS = {k: _mkint(v) for k, v in PXS.items()}
def phi_tail(x): return (arb(1) + x * x) ** (arb(2) / 3) / 4
_cache = {}
def psi_value(pk, xv):
    """rigorous ball for psi(xv) (python float xv), W = (1+x^2)^-PX for x <= 30.25; psi increasing"""
    key = (pk, xv)
    if key in _cache: return _cache[key]
    px = PXS[pk]; pxa = arb(px.numerator) / px.denominator; x = arb(xv)
    with ctx.workprec(130):
        if xv <= A2: v = pxa * (arb(1) + x * x).log()
        else:
            base = pxa * (arb(1) + arb(A2) * arb(A2)).log()
            if xv <= B2: v = base + acb.integral(_INTS[pk], acb(A2), acb(xv)).real
            else:
                p38 = base + acb.integral(_INTS[pk], acb(A2), acb(B2)).real
                v = p38 + phi_tail(x) - phi_tail(arb(B2))
    chk(v); _cache[key] = v
    return v
def W_up(pk, xlo):
    """upper bound (arb) of W = exp(-psi) on [xlo, oo): W decreasing (psi' > 0 on x > 0)"""
    assert xlo > 0
    x0 = math.nextafter(xlo, DN) if xlo > 0.5 else xlo
    return arb(up((-psi_value(pk, x0)).exp()))
def W_lo(pk, xhi):
    """lower bound of W on (0, xhi]"""
    x0 = math.nextafter(xhi, UP)
    return arb(dn((-psi_value(pk, x0)).exp()))
# ------------------------------------------------------------------ psi derivatives P1..P4 on balls (Leibniz for the blend)
x_ = sp.symbols('x')
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        return f"(({bb})**(arb({ex.p})/arb({ex.q})))"
    def _print_Float(self, e): return f"arb('{e}')"
_lam = lambda e: eval("lambda x: " + _Pr().doprint(e), {"arb": arb})
FPS = {k: [_lam(sp.diff(2 * sp.Rational(v.numerator, v.denominator) * x_ / (1 + x_**2), x_, j)) for j in range(4)] for k, v in PXS.items()}
FR = [_lam(sp.diff(x_ * (1 + x_**2) ** sp.Rational(-1, 3) / 3, x_, j)) for j in range(4)]
_t = (x_ - sp.Rational(121, 4)) / 8; _S = _t**4 * (35 - 84 * _t + 70 * _t**2 - 20 * _t**3)
# S^(j) in x: exact rational polynomials in t = (x - 30.25)/8, evaluated by Taylor shift at the (exact) centre of X (avoids the
# cancellation blow-up of the expanded form, which gave +-1.4 on a width-0.02 box near t = 1)
from flint import arb_poly
_tt = sp.symbols('tt'); _St = _tt**4 * (35 - 84 * _tt + 70 * _tt**2 - 20 * _tt**3)
_SDP = [arb_poly([arb(int(c.p)) / int(c.q) for c in reversed(sp.Poly(sp.diff(_St, _tt, j), _tt).all_coeffs())]) for j in range(4)]
def _sd(j, X):
    xm = float(X.mid()); xm_a = arb(xm)                     # exact point (a double)
    tm = (xm_a - arb(30.25)) / 8; h = (X - xm_a) / 8
    with ctx.workprec(160):
        v = _SDP[j](arb_poly([tm, 1]))(h)
    return v / arb(8) ** j
SD = [(lambda X, j=j: _sd(j, X)) for j in range(4)]
def _Pj(X, region, pk):
    FP = FPS[pk]
    if region == "poly": return [FP[j](X) for j in range(4)]
    if region == "rho": return [FR[j](X) for j in range(4)]
    s = [SD[j](X) for j in range(4)]; d = [FR[j](X) - FP[j](X) for j in range(4)]
    return [FP[j](X) + sum(comb(j, i) * s[i] * d[j - i] for i in range(j + 1)) for j in range(4)]
def Pj_box(X, pk):
    """[psi', psi'', psi''', psi''''] enclosure on the x-ball X (x >= 0), clipping at the exact breakpoints 30.25, 38.25"""
    lo, hi = dn(X), up(X)
    cuts = [lo] + [c for c in (A2, B2) if lo < c < hi] + [hi]
    out = None
    for a, b in zip(cuts[:-1], cuts[1:]):
        Xi = I(a, b)
        reg = "poly" if b <= A2 else ("blend" if b <= B2 else "rho")
        P = _Pj(Xi, reg, pk)
        out = P if out is None else [o.union(v) for o, v in zip(out, P)]
    for v in out: chk(v)
    return out
# ------------------------------------------------------------------ boxes
def sig_of_x(x):
    """smallest double >= x + ST.upper()  (so that x(sigma) >= x for all tips in the enclosure)"""
    return math.nextafter(float((arb(x) + ST).upper().mid()), UP)
class Box:
    __slots__ = ("sa", "sb", "X", "xa", "xb", "xc", "sig", "kb", "mb", "pb", "qs", "p", "q", "k", "m")
def make_box(sa, sb):
    b = Box(); b.sa, b.sb = sa, sb
    sig = arb(sa).union(arb(sb)); b.sig = sig
    b.xa, b.xb = dn(arb(sa) - ST), up(arb(sb) - ST); b.xc = 0.5 * (b.xa + b.xb)
    b.X = I(b.xa, b.xb)
    kb, mb, pb, qs = HR.HS.profile_bounds(arb(sa))
    b.kb, b.mb, b.pb, b.qs = kb, mb, pb, qs
    b.p, b.k, b.m = arb(0, pb.upper()), arb(0, kb.upper()), arb(0, mb.upper())
    b.q = sig + arb(0, qs.upper())
    return b
def make_segments(bps, widths):
    """bps: list of ('x', value) or ('s', value) breakpoints (sigma increasing); widths: sigma width per segment.
    returns list of (sa, sb) floats, covering [sigma(first), sigma(last)] with shared exact endpoints."""
    pts = [sig_of_x(v) if t == 'x' else float(v) for t, v in bps]
    out = []
    for (a, b), w in zip(zip(pts[:-1], pts[1:]), widths):
        n = max(1, math.ceil((b - a) / w - 1e-12))
        edges = [a + (b - a) * i / n for i in range(n)] + [b]
        for e0, e1 in zip(edges[:-1], edges[1:]): out.append((e0, e1))
    return out
REG = [(10.25, 16.25), (16.25, 20.0), (20.0, 30.25), (30.25, 40.0), (40.0, 60.0)]
def region_of(xa, xb):
    """regions (indices) overlapped by the x-range of the box"""
    return [i for i, (a, b) in enumerate(REG) if xb > a + 1e-6 and xa < b - 1e-6]
