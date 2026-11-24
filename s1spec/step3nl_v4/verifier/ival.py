"""Minimal vectorised OUTWARD-ROUNDED interval arithmetic for NONNEGATIVE quantities (class P) with numpy float64.
ROUNDING (final pass): every elementary operation (+, *, /, sqrt) is IEEE-754 binary64 round-to-nearest, i.e. correctly rounded:
|fl(x) - x| <= 1/2 ulp(x) (also at binade boundaries: the double below 2^k is 2^k-2^(k-53), the rounding region below 2^k has half-width 2^(k-54)).
The result bounds are then moved one representable number outward with np.nextafter (fd toward -inf, floored at 0 since all quantities are >= 0; fu toward +inf),
so the true value lies in [fd(fl x), fu(fl x)].
An exact zero hi is kept as 0 (sums/products of exact zeros are exact zeros).  Underflow: any quantity that underflows has absolute error
<= 2.5e-324; all intermediate magnitudes are >> 1e-250 except exact zeros, and a final absolute pad (1e-150) is added in sub_signed.
Decimal constants: float() of a decimal literal is correctly rounded -> widened by one ulp each way (cc).  Constants that are computed by float
arithmetic are NOT used: they are computed in 200-bit arb (rhs_iv.Consts) and converted outward (arb2P).  No libm transcendental is used.
Signed quantities only occur at the very end (tot = POS - NEG) and are handled by sub_signed()."""
import numpy as np
from flint import arb
_INF = np.inf
def fd(x): return np.maximum(np.nextafter(x, -_INF), 0.0)
def fu(x):
    y = np.nextafter(x, _INF)
    return np.where(x == 0, 0.0, y) if np.ndim(x) else (0.0 if x == 0 else y)
FD = 1.0; FU = 1.0     # legacy names: no multiplicative inflation any more (kept so that old imports do not break)
class P:
    __slots__ = ("lo", "hi")
    def __init__(s, lo, hi): s.lo = lo; s.hi = hi
    def __add__(a, b):
        if not isinstance(b, P): b = cc(b)
        return P(fd(a.lo + b.lo), fu(a.hi + b.hi))
    __radd__ = __add__
    def __mul__(a, b):
        if not isinstance(b, P): b = cc(b)
        return P(fd(a.lo * b.lo), fu(a.hi * b.hi))
    __rmul__ = __mul__
    def __truediv__(a, b):
        if not isinstance(b, P): b = cc(b)
        assert np.all(b.lo > 0)
        return P(fd(a.lo / b.hi), fu(a.hi / b.lo))
def cc(x):
    """nonnegative decimal constant -> interval"""
    x = float(x); assert x >= 0
    return P(np.float64(fd(x)), np.float64(fu(x)))
def pt(lo, hi):
    return P(np.asarray(lo, float), np.asarray(hi, float))
def sqrt(a): return P(fd(np.sqrt(a.lo)), fu(np.sqrt(a.hi)))
def mx(a, b): return P(np.maximum(a.lo, b.lo), np.maximum(a.hi, b.hi))
def arb2P(x):
    """arb ball (nonneg) -> P, outward (float() of an endpoint is correctly rounded to nearest -> one more ulp outward)"""
    lo = float(x.lower()); hi = float(x.upper())
    return P(np.float64(fd(max(lo, 0.0))), np.float64(fu(hi)))
def psum(lst):
    r = lst[0]
    for t in lst[1:]: r = r + t
    return r
ABSPAD = 1e-150
def sub_signed(Pp, N):
    """returns (lo, hi) float arrays of Pp - N, outward"""
    return np.nextafter(Pp.lo - N.hi, -_INF) - ABSPAD, np.nextafter(Pp.hi - N.lo, _INF) + ABSPAD
def div_signed(lo, hi, L):
    """[lo,hi]/L for positive interval L=[L.lo,L.hi]; returns outward (lo,hi)"""
    rl = np.where(lo >= 0, lo / L.hi, lo / L.lo); rh = np.where(hi >= 0, hi / L.lo, hi / L.hi)
    return np.nextafter(rl, -_INF), np.nextafter(rh, _INF)

# ---------------------------------------------------------------------------------------------------------------------
# Forward-mode gradient bounds: PG carries, besides [lo,hi], an array g of shape (nv,N) of UPPER bounds of the partial
# derivatives wrt the nv physical variables (all >= 0 because every quantity is built from +,*,sqrt,max of nonneg monotone
# pieces).  inf is used when a sqrt argument can be 0.
# ---------------------------------------------------------------------------------------------------------------------
def _mg(g, h):
    with np.errstate(invalid='ignore'):
        return np.where(h == 0, 0.0, g * h)
class PG(P):
    __slots__ = ("g",)
    def __init__(s, lo, hi, g=None): s.lo = lo; s.hi = hi; s.g = g
    def __add__(a, b):
        if not isinstance(b, P): b = cc(b)
        gb = getattr(b, "g", None)
        g = a.g if gb is None else (gb if a.g is None else fu(a.g + gb))
        return PG(fd(a.lo + b.lo), fu(a.hi + b.hi), g)
    __radd__ = __add__
    def __mul__(a, b):
        if not isinstance(b, P): b = cc(b)
        gb = getattr(b, "g", None)
        t1 = None if a.g is None else fu(_mg(a.g, b.hi))
        t2 = None if gb is None else fu(_mg(gb, a.hi))
        g = t1 if t2 is None else (t2 if t1 is None else fu(t1 + t2))
        return PG(fd(a.lo * b.lo), fu(a.hi * b.hi), g)
    __rmul__ = __mul__
    def __truediv__(a, b):
        assert np.all(b.lo > 0)
        g = None if a.g is None else fu(a.g / b.lo)
        return PG(fd(a.lo / b.hi), fu(a.hi / b.lo), g)
def pg_sqrt(a):
    lo = fd(np.sqrt(a.lo)); hi = fu(np.sqrt(a.hi))
    g = None
    if a.g is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            d = fd(2.0 * fd(np.sqrt(a.lo)))
            g = fu(np.where(a.g > 0, np.where(d > 0, a.g / np.where(d > 0, d, 1.0), np.inf), 0.0))
    return PG(lo, hi, g)
def pg_mx(a, b):
    ga = getattr(a, "g", None); gb = getattr(b, "g", None)
    g = ga if gb is None else (gb if ga is None else np.maximum(ga, gb))
    return PG(np.maximum(a.lo, b.lo), np.maximum(a.hi, b.hi), g)
_sqrt_plain = sqrt; _mx_plain = mx
def sqrt(a):
    if hasattr(a, "_dsqrt"): return a._dsqrt()
    return pg_sqrt(a) if isinstance(a, PG) else _sqrt_plain(a)
def mx(a, b):
    if hasattr(a, "_dmx"): return a._dmx(b)
    if hasattr(b, "_dmx"): return b._dmx(a)
    return pg_mx(a, b) if (isinstance(a, PG) or isinstance(b, PG)) else _mx_plain(a, b)
