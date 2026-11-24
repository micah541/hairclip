"""Rigorous Taylor models on one element, local variable s in [0, 1] (x = a + h s):  f(s) in p(s) + [-r, r]  for all s in [0,1],
p an arb_poly (ball coefficients), r >= 0 an upper bound (arb, radius-free use of .upper()).
Constructors:
  TM.poly(p)                         exact polynomial
  from_series(F, a, h, D, nder)      analytic f = F(X) given as a function of an arb_series X = x(t), t = s - 1/2:
                                     coefficients c_0..c_{D-1} of the series at the midpoint; Lagrange remainder c_D(xi) t^D with
                                     c_D enclosed by the series evaluated at the element ball; derivatives f^(m), m <= nder, from the
                                     same data (f^(m) remainder D!/(D-m)! |c_D(T)| 2^-(D-m) / h^m).
Products truncate at degree DMAX (tail moved into r)."""
from flint import arb, arb_poly, arb_series
from math import factorial

DMAX = 34
from flint import ctx as _ctx
_ctx.cap = max(_ctx.cap, 48)
HALF = arb(1) / 2
T2S = arb_poly([-HALF, 1])          # t = s - 1/2

def coeffs(p, n=None):
    c = p.coeffs() if hasattr(p, 'coeffs') else list(p)
    if n is not None: c = c + [arb(0)] * (n - len(c))
    return c

def abs_sum(p):
    tot = arb(0)
    for c in coeffs(p): tot += abs(c)
    return tot

def sup01(p):
    """upper bound (arb) of |p(s)| on [0,1]"""
    if p.degree() < 0: return arb(0)
    a = abs_sum(p); b = abs(p(arb(HALF, HALF)))
    return a if a.upper() < b.upper() else b

def up(v):
    """arb -> exact nonnegative upper bound"""
    return arb(abs(v).upper())

class TM:
    __slots__ = ('p', 'r')
    def __init__(self, p, r=None):
        self.p = p if isinstance(p, arb_poly) else arb_poly([p]); self.r = arb(0) if r is None else up(r)
    @staticmethod
    def poly(p): return TM(p)
    def __add__(self, o):
        o = _tm(o); return TM(self.p + o.p, self.r + o.r)
    __radd__ = __add__
    def __neg__(self): return TM(-self.p, self.r)
    def __sub__(self, o): return self + (-_tm(o))
    def __rsub__(self, o): return _tm(o) + (-self)
    def __mul__(self, o):
        if isinstance(o, (arb, int, float)):
            o = arb(o); return TM(self.p * o, self.r * abs(o))
        o = _tm(o); pr = self.p * o.p
        r = sup01(self.p) * o.r + sup01(o.p) * self.r + self.r * o.r
        c = coeffs(pr)
        if len(c) > DMAX + 1:
            tail = arb(0)
            for v in c[DMAX + 1:]: tail += abs(v)
            pr = arb_poly(c[:DMAX + 1]); r += tail
        return TM(pr, r)
    __rmul__ = __mul__
    def sup(self): return sup01(self.p) + self.r
    def integral01(self):
        """ball containing int_0^1 f ds"""
        return self.p.integral()(arb(1)) + arb(0, self.r.upper())
    def int_against(self, q):
        """ball containing int_0^1 f q ds for an exact polynomial q"""
        I = (self.p * q).integral()(arb(1))
        return I + arb(0, (self.r * abs_int(q)).upper())
    def antideriv(self, h):
        """TM of h * int_0^s f (remainder h r s <= h r)"""
        return TM(self.p.integral() * h, self.r * abs(h))

def abs_int(q):
    """upper bound of int_0^1 |q|"""
    tot = arb(0)
    for j, c in enumerate(coeffs(q)): tot += abs(c) / (j + 1)
    return tot

def _tm(o):
    if isinstance(o, TM): return o
    if isinstance(o, arb_poly): return TM(o)
    return TM(arb_poly([arb(o)]))

def ser_coeffs(S, n):
    c = S.coeffs(); return c + [arb(0)] * (n - len(c))

def from_series(F, a, h, D=24, nder=0):
    """TMs [f, f', ..., f^(nder)] (x-derivatives) of f = F(X), X = a + h (1/2 + t), on the element [a, a+h]"""
    xm = a + h * HALF
    Sm = F(arb_series([xm, h], prec=D + 1)); Sb = F(arb_series([xm + h * arb(0, HALF), h], prec=D + 1))
    cm = ser_coeffs(Sm, D + 1); cD = ser_coeffs(Sb, D + 1)[D]
    out = []
    for m in range(nder + 1):
        # f^(m)(x) = h^-m sum_j c_{j+m} (j+m)!/j! t^j  + remainder  D!/(D-m)! cD t^(D-m) h^-m
        pt = arb_poly([cm[j + m] * (factorial(j + m) // factorial(j)) for j in range(D - m)]) * (arb(1) / h ** m)
        r = abs(cD) * (factorial(D) // factorial(D - m)) * HALF ** (D - m) / abs(h) ** m
        out.append(TM(pt(T2S), r))
    return out

def poly_eval_series(p, X, prec):
    """Horner of arb_poly p at an arb_series X (series length prec)"""
    c = coeffs(p)
    if not c: return arb_series([arb(0)], prec=prec)
    acc = arb_series([c[-1]], prec=prec)
    for v in reversed(c[:-1]): acc = acc * X + v
    return acc
