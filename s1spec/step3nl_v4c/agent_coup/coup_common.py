"""Rigorous (arb) helpers for the far-field coupling constants.  See COUPLING_NOTES.md.
Conventions: x = sigma - tip (tip an interval), rho = exp(-phi), phi = (1+x^2)^(2/3)/4,
W = exp(-phi_I), phi_I' = (1-s) 2P x/(1+x^2) + s (1/3) x (1+x^2)^(-1/3), s = S((x-30.25)/8) clamped, phi_I(0) = 0.
Cutoffs are C^3 smoothsteps built from S(t) = t^4 (35-84t+70t^2-20t^3).
Pieces are float intervals [a,b] (exact endpoints); every piece enclosure is a rigorous ball."""
import sys, math, time
from math import comb
import sympy as sp
from flint import arb, acb, arb_poly, ctx
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
ctx.prec = 128
UP, DN = math.inf, -math.inf
TIP = arb('-2.2490776515195').union(arb('-2.2490776505754'))
A2, B2 = 30.25, 38.25

# ----------------------------------------------------------------- ball utilities
def chk(z):
    assert z.is_finite(), 'non-finite ball'
    return z
def upf(z):
    """rigorous float upper bound of a ball"""
    chk(z); return math.nextafter(float(z.upper()), UP)
def lof(z):
    chk(z); return math.nextafter(float(z.lower()), DN)
def Iv(a, b):
    """ball containing the float interval [a,b]"""
    mid = 0.5 * (a + b); rad = max(mid - a, b - mid)
    return arb(mid, math.nextafter(rad, UP))
def pos(z):
    """enclosure of max(z,0) (pointwise) for a ball z"""
    lo, hi = z.lower(), z.upper()
    plo = lo if lo > 0 else arb(0)
    phi_ = hi if hi > 0 else arb(0)
    return chk(plo.union(phi_))
def absb(z):
    """enclosure of |z| (pointwise)"""
    lo, hi = z.lower(), z.upper()
    alo, ahi = abs(lo), abs(hi)
    mx = ahi if ahi >= alo else alo
    if lo > 0 or hi < 0: mn = alo if alo <= ahi else ahi
    else: mn = arb(0)
    return chk(mn.union(mx))

# ----------------------------------------------------------------- smoothstep derivatives (Taylor-shifted, stable)
_tt = sp.symbols('tt'); _St = _tt**4 * (35 - 84 * _tt + 70 * _tt**2 - 20 * _tt**3)
_SDP = [arb_poly([arb(int(c.p)) / int(c.q) for c in reversed(sp.Poly(sp.diff(_St, _tt, j), _tt).all_coeffs())]) for j in range(6)]
def sderiv(j, T):
    """S^(j)(t) on ball T (T must lie in [0,1] up to rounding; polynomial extension is exact anyway)"""
    tm = arb(float(T.mid())); h = T - tm
    with ctx.workprec(200):
        v = _SDP[j](arb_poly([tm, 1]))(h)
    return chk(v)
def S_all(T, n=5): return [sderiv(j, T) for j in range(n)]

# ----------------------------------------------------------------- cutoffs: lists [c0..c4] of derivatives wrt x (or sigma)
def cut_down(Z, z0, z1):
    """C^3 smoothstep 1 at z<=z0, 0 at z>=z1: S((z1 - z)/(z1 - z0)); derivatives wrt z"""
    L = arb(z1) - arb(z0); u = (arb(z1) - Z) / L
    Sj = S_all(u); return [Sj[j] * (-1 / L) ** j for j in range(5)]
def cut_chi(SIG): return cut_down(SIG, 29, 32)            # chi_32(sigma) = S((32-sigma)/3)
def cut_zeta(X):
    """zeta = 1 - S((16.25 - x)/6)"""
    L = arb(6); u = (arb('16.25') - X) / L
    Sj = S_all(u); c = [-Sj[j] * (-1 / L) ** j for j in range(5)]; c[0] = 1 + c[0]; return c
def cut_zo(X, R): return cut_down(X, R - 3, R)             # zeta_out = S((R-x)/3)

# ----------------------------------------------------------------- symbolic coefficient polynomials
cs = sp.symbols('c0:6'); Ps = sp.symbols('P1:6'); k_, m_, p_, q_ = sp.symbols('k m p q')
_rules = {**{cs[i]: cs[i + 1] for i in range(5)}, **{Ps[i]: Ps[i + 1] for i in range(4)},
          k_: m_, m_: p_ / 4, p_: -k_ * q_, q_: 1 + k_ * p_}
def Dsym(e): return sp.expand(sum(sp.diff(e, v) * r for v, r in _rules.items() if e.has(v)))
def Dw(e): return sp.expand(Dsym(e) - Ps[0] * e)          # (Wt e)'/Wt
def AB_poly(a):
    """A/Wt, B/Wt for  -<[L,cut]f, cut f>_Wt = int A f^2 + B f'^2, drift coefficient (a k m + q/4)"""
    c0, c1, c2, c3, c4 = cs[:5]
    g1, g2, g3, g4 = c0 * c1, c0 * c2, c0 * c3, c0 * c4
    A = (-2 * Dw(Dw(Dw(g1))) + 3 * Dw(Dw(g2)) - 2 * Dw(g3) + g4 - Dw(k_**2 * g1)
         + c0 * (k_**2 * c2 + (a * k_ * m_ + q_ / 4) * c1))
    B = 6 * Dw(g1) - 6 * g2
    return sp.expand(A), sp.expand(B)
def coef_poly(a):
    c0, c1, c2, c3, c4 = cs[:5]
    return [c4 + k_**2 * c2 + (a * k_ * m_ + q_ / 4) * c1, 4 * c3 + 2 * k_**2 * c1, 6 * c2, 4 * c1]   # coef_0..coef_3
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n)); return f"({prod})" if n > 0 else f"(1/({prod}))"
        return f"(({bb})**(arb({ex.p})/arb({ex.q})))"
    def _print_Float(self, e): return f"arb('{e}')"
def compile_sym(exprs, args):
    src = "lambda " + ", ".join(map(str, args)) + ": [" + ", ".join(_Pr().doprint(e) for e in exprs) + "]"
    return eval(src, {"arb": arb})
_ARGS = list(cs[:5]) + list(Ps[:4]) + [k_, m_, p_, q_]
_AB = {}
def get_AB(a):
    if a not in _AB:
        A, B = AB_poly(a)
        for e in (A, B): assert not (e.has(cs[5]) or e.has(Ps[3]) or e.has(Ps[4])), 'unexpected high derivative'
        _AB[a] = compile_sym([A, B], _ARGS)
    return _AB[a]
_CF = {}
def get_coef(a):
    if a not in _CF: _CF[a] = compile_sym(coef_poly(a), _ARGS)
    return _CF[a]

# ----------------------------------------------------------------- profile (k, m, p, q) on sigma >= 8
_pb = {}
def prof_bounds(sig_lo):
    assert sig_lo >= 8.0, 'profile bounds only for sigma >= 8'
    key = math.floor(sig_lo * 16) / 16
    if key not in _pb:
        kb, mb, pb, qs = H.profile_bounds(arb(key))
        _pb[key] = tuple(upf(v) for v in (kb, mb, pb, qs))
    return _pb[key]
def kmpq(SIG):
    sl = lof(SIG)
    kb, mb, pb, qs = prof_bounds(sl)
    return arb(0, kb), arb(0, mb), arb(0, pb), SIG + arb(0, qs)

# ----------------------------------------------------------------- phi, rho
_x = sp.symbols('x')
def _lam(e):
    return eval("lambda x: " + _Pr().doprint(e), {"arb": arb})
_FR = [_lam(sp.diff(_x * (1 + _x**2) ** sp.Rational(-1, 3) / 3, _x, j)) for j in range(4)]    # d^j phi'
def phi_ball(X): return chk(((1 + X * X) ** (arb(2) / 3)) / 4)
def rho_data(X):
    """(phi, [P1..P4]) for rho"""
    return phi_ball(X), [chk(f(X)) for f in _FR[:4]]
_t = (_x - sp.Rational(121, 4)) / 8
class Weight:
    """W = exp(-phi_I) for given P (sympy Rational)"""
    H_ = 512                                   # table spacing 2^-9 in the blend region
    def __init__(self, PX):
        self.PX = sp.Rational(PX); self.Pa = arb(self.PX.p) / arb(self.PX.q)
        self.FP = [_lam(sp.diff(2 * self.PX * _x / (1 + _x**2), _x, j)) for j in range(4)]
        self.twoP = int(2 * self.PX); assert 2 * self.PX == self.twoP
        self.build_table()
    # exact phi_I at the blend grid nodes
    def _integrand(self, z, _an):
        qq = acb(1) + z * z
        q13 = (qq.log() * (-acb(1) / 3)).exp()
        Ap = self.twoP * z / qq; Bp = z * q13 / 3
        t = (z - acb('30.25')) / 8; t2 = t * t; t3 = t2 * t; t4 = t2 * t2
        s = t4 * (35 - 84 * t + 70 * t2 - 20 * t3)
        return (1 - s) * Ap + s * Bp
    def build_table(self):
        h = self.H_; j0, j1 = int(A2 * h), int(B2 * h)
        self.j0 = j0
        v = self.Pa * (arb(1) + arb('30.25') * arb('30.25')).log()
        tab = [v]
        with ctx.workprec(130):
            for j in range(j0, j1):
                seg = acb.integral(self._integrand, acb(j / h), acb((j + 1) / h)).real
                v = v + seg; tab.append(v)
        self.tab = [chk(t) for t in tab]
        self.phiI_B2 = self.tab[-1]
    def sblend(self, X):
        T = (X - arb('30.25')) / 8
        return [sderiv(j, T) / arb(8) ** j for j in range(4)]
    def Pj_region(self, a, b, X):
        """[P1..P4] on X=[a,b] lying in one region"""
        FP = [chk(f(X)) for f in self.FP]
        if b <= A2: return FP
        FRj = [chk(f(X)) for f in _FR[:4]]
        if a >= B2: return FRj
        s = self.sblend(X); d = [FRj[j] - FP[j] for j in range(4)]
        return [chk(FP[j] + sum(comb(j, i) * s[i] * d[j - i] for i in range(j + 1))) for j in range(4)]
    def cuts(self, a, b):
        pts = [a] + [c for c in (A2, B2) if a < c < b] + [b]
        return list(zip(pts[:-1], pts[1:]))
    def Pj(self, a, b):
        out = None
        for u, v in self.cuts(a, b):
            P = self.Pj_region(u, v, Iv(u, v))
            out = P if out is None else [o.union(w) for o, w in zip(out, P)]
        return out
    def phiI_point(self, xa):
        """phi_I at the (float or exact-ball) point xa >= 0"""
        xa = arb(xa) if not isinstance(xa, arb) else xa
        xf = float(xa.mid())
        if float(xa.upper()) <= A2: return self.Pa * (1 + xa * xa).log()
        if float(xa.lower()) >= B2: return chk(self.phiI_B2 + phi_ball(xa) - phi_ball(arb(B2)))
        assert float(xa.lower()) >= A2 - 1e-12 and float(xa.upper()) <= B2 + 1e-12, 'point straddles region boundary'
        j = round(xf * self.H_); j = min(max(j, self.j0), self.j0 + len(self.tab) - 1)
        g = arb(j) / self.H_
        a_, b_ = min(float(g), xf), max(float(g), xf)
        Pp = self.Pj_region(max(a_, A2), min(b_, B2), Iv(max(a_, A2), min(b_, B2)))[0]
        return chk(self.tab[j - self.j0] + (xa - g) * Pp)
    def phiI(self, a, b):
        """ball for phi_I over float interval [a,b] (phi_I increasing)"""
        lo = self.phiI_point(arb(a)); hi = self.phiI_point(arb(b))
        return chk(lo.lower().union(hi.upper()))
    def phiI_ball(self, X):
        """ball for phi_I over an arb ball X (outer bounds taken as exact endpoints)"""
        lo, hi = X.lower(), X.upper()
        return chk(self.phiI_point(lo).lower().union(self.phiI_point(hi).upper()))
    def W_ball(self, a, b): return chk((-self.phiI(a, b)).exp())
    def W_point(self, x): return chk((-self.phiI_point(arb(x))).exp())

# ----------------------------------------------------------------- Phi(x) = int_{xh}^x 1/W (upper bounds)
class PhiHardy:
    def __init__(self, Wt, xh=16.25, xmax=35.0, h=1.0 / 512):
        self.Wt, self.xh, self.h = Wt, xh, h
        n = int(round((xmax - xh) / h)); self.n = n
        cum = [arb(0)]; acc = arb(0)
        for i in range(n):
            g = xh + i * h
            acc = acc + self.seg_upper(g, g + h)
            cum.append(acc)
        self.cum = cum
    def seg_upper(self, g, g1):
        M = self.Wt.Pj(g, g1)[0].upper()                    # sup phi_I' on segment (>0)
        assert M > 0
        e = self.Wt.phiI_point(arb(g)).exp()
        L = arb(g1) - arb(g)
        return chk(e * (M * L).expm1() / M)
    def upper(self, xf):
        """upper bound of Phi(xf) for float xf in [xh, xmax]"""
        i = int(math.floor((xf - self.xh) / self.h + 1e-12)); i = min(max(i, 0), self.n)
        g = self.xh + i * self.h
        r = self.cum[i]
        if xf > g: r = r + self.seg_upper(g, xf)
        return chk(r.upper())

# ----------------------------------------------------------------- piece context
def piece_xs(a, b, mode):
    """mode 'x': variable is x in [a,b] (sigma = x + tip); mode 'sigma': variable is sigma in [a,b] (x = sigma - tip)"""
    V = Iv(a, b)
    return (V, V + TIP) if mode == 'x' else (V - TIP, V)

# ----------------------------------------------------------------- adaptive bisection
def _split(a, b): m = 0.5 * (a + b); return [(a, m), (m, b)]
def grid(ranges, n0):
    out = []
    for a, b in ranges:
        n = n0 if not isinstance(n0, dict) else n0.get((a, b), 64)
        pts = [a + (b - a) * i / n for i in range(n)] + [b]
        out += list(zip(pts[:-1], pts[1:]))
    return out
def bnb_sup(fun, ranges, n0=256, rel=1e-3, atol=0.0, maxpieces=400000):
    """rigorous upper bound of sup of a (>=0 or signed) quantity; fun(a,b) -> enclosure ball for the piece.
    returns (upper_float, lower_float [valid lower bound of sup], npieces)"""
    pcs = [(a, b, fun(a, b)) for a, b in grid(ranges, n0)]
    L = max(lof(e) for _, _, e in pcs)
    stack, done = pcs, []
    while stack:
        a, b, e = stack.pop()
        u = upf(e)
        if u <= L * (1 + rel) + atol or len(done) + len(stack) > maxpieces or (b - a) < 1e-14 * max(1.0, abs(b)):
            done.append((a, b, e)); continue
        for (a2, b2) in _split(a, b):
            e2 = fun(a2, b2); L = max(L, lof(e2)); stack.append((a2, b2, e2))
    return max(upf(e) for _, _, e in done), L, len(done)
def bnb_int(fun, ranges, n0=256, rel=2e-3, rounds=40, maxpieces=400000):
    """rigorous upper bound of int fun (fun >= 0 enclosures); returns (upper_float, lower_float, npieces)"""
    pcs = [(a, b, fun(a, b)) for a, b in grid(ranges, n0)]
    for _ in range(rounds):
        T = sum(upf(e) * (b - a) for a, b, e in pcs)
        wsum = sum((upf(e) - lof(e)) * (b - a) for a, b, e in pcs)
        if wsum <= rel * T or len(pcs) > maxpieces: break
        thr = rel * T / (2 * len(pcs)); new = []
        for a, b, e in pcs:
            if (upf(e) - lof(e)) * (b - a) > thr:
                for (a2, b2) in _split(a, b): new.append((a2, b2, fun(a2, b2)))
            else: new.append((a, b, e))
        pcs = new
    U = sum(upf(e) * (b - a) for a, b, e in pcs); Lw = sum(max(lof(e), 0) * (b - a) for a, b, e in pcs)
    return U * (1 + 1e-9), Lw * (1 - 1e-9), len(pcs)
