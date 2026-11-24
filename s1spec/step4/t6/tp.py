"""T6 layer 1: rigorous Taylor models in the slab time and bivariate (s, tc) polynomials with arb coefficients.

Variables: tc in [-1, 1] is the centred slab time (tau = tau0 + dtau (1 + tc)/2); s in [-1, 1] is the centred
element coordinate (x = x_e + h_e (1 + s)/2).

TP  = arb_poly in tc.  Semantics ("pointwise in tc"): the true function f(tc) lies in the ball evaluation of the
      polynomial at tc, for every tc in [-1, 1].  Remainders (truncated degrees, series tails) are folded into the
      radius of the constant coefficient.  TPs are never differentiated in tc unless they are exact polynomials
      (the Lagrange interpolants of the stored stages).
BP  = list of arb_poly in s indexed by the power of tc.  Two kinds:
      ex=True  (exact in s): for every tc the true function is a polynomial in s whose s-coefficients lie in the ball
               evaluations of the tc-polynomials; s-derivatives and s-integrals are exact operations.
      ex=False (pointwise): for every (s, tc) the true value lies in the ball evaluation.  s-integrals are allowed
               (|int_{-1}^s R| <= 2R is covered by the folding), s-derivatives are not.
Truncations fold the dropped part per s-coefficient into the tc^0 coefficient, which is valid for both kinds."""
from flint import arb, arb_poly, ctx

ZERO = arb(0)
ONE = arb(1)


def usum(vals):
    """exact upper bound of a sum of nonnegative arbs"""
    s = ZERO
    for v in vals: s = s + v
    return s.upper()


def ub(x):
    return x.abs_upper()


def addrad(c, R):
    """ball c enlarged by R >= 0 (exact arb)"""
    return c + arb(0, R) if R != 0 else c


# ----------------------------------------------------------------------------------------------- TPs (univariate)
def tconst(c):
    return arb_poly([c])


def ttrunc(p, D):
    c = p.coeffs()
    if len(c) <= D + 1: return p
    R = usum(ci.abs_upper() for ci in c[D + 1:])
    return arb_poly([addrad(c[0], R)] + c[1:D + 1])


def tmul(a, b, D):
    return ttrunc(a * b, D)


def tsup(p):
    """exact upper bound of sup_{|tc|<=1} |p(tc)|"""
    c = p.coeffs()
    return usum(v.abs_upper() for v in c) if c else ZERO


def trange(p):
    """ball containing p([-1, 1]) (constant term +- sum of the others)"""
    c = p.coeffs()
    if not c: return ZERO
    return addrad(c[0], usum(v.abs_upper() for v in c[1:]))


def tmid0(p):
    """exact midpoint of the constant coefficient"""
    c = p.coeffs()
    return arb(c[0].mid()) if c else ZERO


def tpow(g, alpha, D, tol_exp=-150, Dint=None):
    """g(tc)^alpha as a TP: g = g0 (1 + v), binomial series in v with tail <= |g0^alpha| eps^(K+1)/(1-eps), |binom| <= 1
    (valid for -1 <= alpha <= 1).  Products are kept to degree Dint (default D + 24) and truncated once at the end."""
    Dint = D + 24 if Dint is None else Dint
    alpha = arb(alpha) if not isinstance(alpha, arb) else alpha
    g0 = tmid0(g)
    if not (g0 > 0): raise ValueError('tpow: g0 must be > 0')
    v = (g - tconst(g0)) * (ONE / g0)
    eps = tsup(v)
    if not (eps < 0.5): raise ValueError(f'tpow: variation too large {eps}')
    # number of terms: eps^(K+1) < 2^tol_exp
    import math
    e = float(eps.upper())
    K = 1 if e == 0 else min(200, max(1, int(math.ceil(tol_exp * math.log(2) / math.log(e))) + 1))
    acc = tconst(ONE); vk = tconst(ONE); bk = ONE
    for k in range(1, K + 1):
        bk = bk * (alpha - (k - 1)) / k
        vk = tmul(vk, v, Dint)
        acc = acc + vk * bk
    tail = (eps ** (K + 1) / (ONE - eps)).upper() if e > 0 else ZERO
    c0 = g0 ** alpha
    out = acc * c0
    cc = out.coeffs()
    cc[0] = addrad(cc[0], (c0.abs_upper() * tail).upper())
    return ttrunc(arb_poly(cc), D)


def tcis(th, D, tol_exp=-150, Dint=None):
    """(cos th(tc), sin th(tc)) as TPs: th = th0 + d, Taylor in d with Lagrange tail eps^(K+1)/(K+1)!"""
    th0 = tmid0(th)
    d = th - tconst(th0)
    eps = tsup(d)
    import math
    e = float(eps.upper())
    if e >= 1.5: raise ValueError(f'tcis: variation too large {eps}')
    # cos d, sin d series
    K = 2
    if e > 0:
        lf = 0.0
        while True:
            K += 1; lf += math.log(K)
            if (K) * math.log(e) - lf < tol_exp * math.log(2): break
    else: K = 1
    Dint = D + 24 if Dint is None else Dint
    cs, sn = tconst(ONE), tconst(ZERO)
    dk = tconst(ONE); fk = ONE
    for k in range(1, K + 1):
        dk = tmul(dk, d, Dint); fk = fk * k
        term = dk * (ONE / fk)
        m = k % 4
        if m == 1: sn = sn + term
        elif m == 2: cs = cs - term
        elif m == 3: sn = sn - term
        else: cs = cs + term
    tail = (eps ** (K + 1) / arb.fac_ui(K + 1)).upper() if e > 0 else ZERO
    c0, s0 = th0.cos(), th0.sin()
    C = cs * c0 - sn * s0; S = cs * s0 + sn * c0
    cc = C.coeffs(); sc = S.coeffs()
    cc[0] = addrad(cc[0], tail); sc[0] = addrad(sc[0], tail)
    return ttrunc(arb_poly(cc), D), ttrunc(arb_poly(sc), D)


# ----------------------------------------------------------------------------------------------- BPs (bivariate)
class BP:
    __slots__ = ('P', 'ex')

    def __init__(self, P, ex=True):
        self.P = P if P else [arb_poly([0])]
        self.ex = ex

    # constructors
    @staticmethod
    def from_tp(tp):
        c = tp.coeffs()
        return BP([arb_poly([v]) for v in c] if c else [arb_poly([0])])

    @staticmethod
    def from_spoly(sp):
        return BP([sp])

    @staticmethod
    def zero():
        return BP([arb_poly([0])])

    def copy(self):
        return BP(list(self.P), self.ex)

    def tdeg(self): return len(self.P) - 1

    def sdeg(self): return max(p.degree() for p in self.P)

    # arithmetic
    def __add__(self, o):
        a, b = self.P, o.P
        if len(a) < len(b): a, b = b, a
        out = list(a)
        for j, q in enumerate(b): out[j] = out[j] + q
        return BP(out, self.ex and o.ex)

    def __neg__(self):
        return BP([-p for p in self.P], self.ex)

    def __sub__(self, o):
        return self + (-o)

    def scale(self, c):
        """times an arb scalar"""
        return BP([p * c for p in self.P], self.ex)

    def mul(self, o, D):
        """truncated product (tc-degree <= D); one Kronecker-packed arb_poly product (s-blocks of length L)"""
        a, b = self.P, o.P
        Sa = max(p.degree() for p in a); Sb = max(p.degree() for p in b)
        if Sa < 0 or Sb < 0: return BP([arb_poly([0])], self.ex and o.ex)
        L = Sa + Sb + 1
        z = arb(0)
        ca = []
        for p in a:
            c = p.coeffs(); ca.extend(c); ca.extend([z] * (L - len(c)))
        cb = []
        for p in b:
            c = p.coeffs(); cb.extend(c); cb.extend([z] * (L - len(c)))
        C = (arb_poly(ca) * arb_poly(cb)).coeffs()
        nt = (len(C) + L - 1) // L
        out = [arb_poly(C[m * L:(m + 1) * L]) for m in range(min(nt, D + 1))]
        if nt > D + 1:
            out[0] = _fold(out[0], [arb_poly(C[m * L:(m + 1) * L]) for m in range(D + 1, nt)])
        return BP(out, self.ex and o.ex)

    def tmul(self, tp, D):
        """times a TP (constant in s): direct convolution in tc; dropped part bounded per s-coefficient by
        R_k = sum_j |P_j[k]| * sum_{l > D-j} |t_l|"""
        t = tp.coeffs()
        if not t: return BP([arb_poly([0])], self.ex)
        P = self.P
        nt = len(t)
        out = [None] * min(len(P) + nt - 1, D + 1)
        for j, pj in enumerate(P):
            if j > D: break
            for l in range(min(nt, D + 1 - j)):
                v = pj * t[l]
                out[j + l] = v if out[j + l] is None else out[j + l] + v
        out = [p if p is not None else arb_poly([0]) for p in out]
        if len(P) + nt - 1 > D + 1:
            ta = [v.abs_upper() for v in t]
            tails = []
            for j in range(len(P)):
                lo = max(0, D + 1 - j)
                tails.append(usum(ta[lo:]) if lo < nt else ZERO)
            R = []
            for j, pj in enumerate(P):
                Tj = tails[j]
                if Tj == 0: continue
                for k, v in enumerate(pj.coeffs()):
                    u = v.abs_upper() * Tj
                    if k < len(R): R[k] = R[k] + u
                    else: R.append(u)
            if R:
                c = out[0].coeffs()
                if len(c) < len(R): c = c + [arb(0)] * (len(R) - len(c))
                for k, rk in enumerate(R): c[k] = addrad(c[k], rk.upper())
                out[0] = arb_poly(c)
        return BP(out, self.ex)

    def ttrunc(self, D):
        if len(self.P) <= D + 1: return self
        out = list(self.P[:D + 1]); out[0] = _fold(out[0], self.P[D + 1:])
        return BP(out, self.ex)

    # calculus in s
    def ds(self):
        assert self.ex, 'ds() of a pointwise BP'
        return BP([p.derivative() for p in self.P], True)

    def int_m1(self):
        """int_{-1}^s"""
        out = []
        for p in self.P:
            q = p.integral(); out.append(q - q(-ONE))
        return BP(out, self.ex)

    def at(self, s0):
        """value at s = s0 (exact arb) -> TP"""
        return arb_poly([p(s0) for p in self.P])

    def shift_s(self, a, w):
        """f(a + w sigma) as a BP in sigma (a, w exact arbs)"""
        lin = arb_poly([a, w])
        return BP([p(lin) for p in self.P], self.ex)

    def fold_pointwise(self, R):
        """add a pointwise remainder R (sup over the box) -> pointwise BP"""
        out = list(self.P)
        c = out[0].coeffs() or [arb(0)]
        c[0] = addrad(c[0], R)
        out[0] = arb_poly(c)
        return BP(out, False)

    def coef_sum_t(self, j):
        """exact upper bound of sup_{|s|<=1} |coefficient of tc^j|"""
        if j >= len(self.P): return ZERO
        return usum(v.abs_upper() for v in self.P[j].coeffs())

    def sup(self):
        """crude exact upper bound of sup |f| over [-1,1]^2"""
        return usum(usum(v.abs_upper() for v in p.coeffs()) for p in self.P)


def _fold(p0, polys):
    """add to the s-coefficients of p0 the radii sum_m |c_mk| of the polys"""
    R = []
    for q in polys:
        for k, v in enumerate(q.coeffs()):
            u = v.abs_upper()
            if k < len(R): R[k] = R[k] + u
            else: R.append(u)
    c = p0.coeffs()
    if len(c) < len(R): c = c + [arb(0)] * (len(R) - len(c))
    for k, rk in enumerate(R): c[k] = addrad(c[k], rk.upper())
    return arb_poly(c)


class CBP:
    """complex BP (re, im)"""
    __slots__ = ('re', 'im')

    def __init__(self, re, im): self.re, self.im = re, im

    def __add__(self, o): return CBP(self.re + o.re, self.im + o.im)

    def __sub__(self, o): return CBP(self.re - o.re, self.im - o.im)

    def conj(self): return CBP(self.re, -self.im)

    def mul(self, o, D):
        ac = self.re.mul(o.re, D); bd = self.im.mul(o.im, D)
        ad = self.re.mul(o.im, D); bc = self.im.mul(o.re, D)
        return CBP(ac - bd, ad + bc)

    def ttrunc(self, D):
        return CBP(self.re.ttrunc(D), self.im.ttrunc(D))

    def rmul(self, b, D):
        """times a real BP"""
        return CBP(self.re.mul(b, D), self.im.mul(b, D))

    def tmul(self, tp, D):
        return CBP(self.re.tmul(tp, D), self.im.tmul(tp, D))

    def ctmul(self, tre, tim, D):
        """times a complex TP (tre + i tim)"""
        A, B = BP.from_tp(tre), BP.from_tp(tim)
        return self.mul(CBP(A, B), D)

    def int_m1(self): return CBP(self.re.int_m1(), self.im.int_m1())

    def at(self, s0): return self.re.at(s0), self.im.at(s0)


def strunc(bp, rtol=1e-30, smin=4):
    """pointwise s-truncation: drop s-powers k > K, K the smallest with sum_{k>K} sum_j |c_jk| <= rtol * total;
    the dropped part (|s|, |tc| <= 1) is folded into the s^0 tc^0 radius.  Returns a pointwise BP."""
    deg = bp.sdeg()
    if deg <= smin: return BP(list(bp.P), False)
    col = [ZERO] * (deg + 1)
    for p in bp.P:
        for k, v in enumerate(p.coeffs()): col[k] = col[k] + v.abs_upper()
    tot = usum(col)
    if tot == 0: return BP(list(bp.P), False)
    lim = tot * rtol
    acc = ZERO; K = deg
    for k in range(deg, smin, -1):
        acc2 = acc + col[k]
        if acc2.upper() > lim: break
        acc = acc2; K = k - 1
    if K >= deg: return BP(list(bp.P), False)
    R = acc.upper()
    out = []
    for p in bp.P:
        c = p.coeffs()
        out.append(arb_poly(c[:K + 1]))
    c0 = out[0].coeffs() or [arb(0)]
    c0[0] = addrad(c0[0], R)
    out[0] = arb_poly(c0)
    return BP(out, False)


def tsup_sub(p, a=-1, b=1, K=16):
    """exact upper bound of sup |p(tc)| for tc in [a, b] (subdivision into K pieces, Taylor-shifted coefficient sums)"""
    a, b = arb(a), arb(b)
    best = ZERO
    w = (b - a) / (2 * K)
    for k in range(K):
        m = a + (2 * k + 1) * w
        q = p(arb_poly([m, w]))
        v = tsup(q)
        if v > best: best = v
        elif not (v <= best): best = max(best, v.upper())
    return best.upper() if isinstance(best, arb) else best


def texp(p, D, tol_exp=-150, Dint=None):
    """exp(p(tc)) as a TP (Taylor in p - p0 with Lagrange tail)"""
    import math
    p0 = tmid0(p); d = p - tconst(p0)
    eps = tsup(d); e = float(eps.upper())
    Dint = D + 24 if Dint is None else Dint
    K = 1; lf = 0.0
    if e > 0:
        while True:
            K += 1; lf += math.log(K)
            if K * math.log(e) - lf < tol_exp * math.log(2) and K > 2: break
    acc = tconst(ONE); dk = tconst(ONE); fk = ONE
    for k in range(1, K + 1):
        dk = tmul(dk, d, Dint); fk = fk * k
        acc = acc + dk * (ONE / fk)
    tail = (eps ** (K + 1) / arb.fac_ui(K + 1) * eps.exp()).upper() if e > 0 else ZERO
    out = acc * p0.exp()
    c = out.coeffs(); c[0] = addrad(c[0], (tail * p0.exp().abs_upper()).upper())
    return ttrunc(arb_poly(c), D)
