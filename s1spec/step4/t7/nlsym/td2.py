"""Second-order forward-mode automatic differentiation in two variables over arb balls (for nlterms7.py).

TD2(v, g0, g1, h00, h01, h11) carries a value and the first and second partial derivatives in (e0, e1).  Evaluating an
expression with interval inputs e0 = TD2.var(E0, 0), e1 = TD2.var(E1, 1) (E0, E1 balls) gives, in every component, a ball
that contains that derivative at every point of the box E0 x E1 (interval chain rule).  Supported: +, -, *, /, integer
powers, .sqrt(), mixing with arb and int (arb returns NotImplemented for foreign types, so the reflected operators run)."""
from flint import arb

Z = arb(0)


class TD2:
    __slots__ = ('v', 'g0', 'g1', 'h00', 'h01', 'h11')

    def __init__(s, v, g0=Z, g1=Z, h00=Z, h01=Z, h11=Z):
        s.v, s.g0, s.g1, s.h00, s.h01, s.h11 = v, g0, g1, h00, h01, h11

    @staticmethod
    def var(val, i):
        return TD2(val, arb(1), Z) if i == 0 else TD2(val, Z, arb(1))

    def _c(s, f):
        return TD2(f(s.v), f(s.g0), f(s.g1), f(s.h00), f(s.h01), f(s.h11))

    def __add__(a, b):
        if isinstance(b, TD2): return TD2(a.v + b.v, a.g0 + b.g0, a.g1 + b.g1, a.h00 + b.h00, a.h01 + b.h01, a.h11 + b.h11)
        return TD2(a.v + b, a.g0, a.g1, a.h00, a.h01, a.h11)
    __radd__ = __add__

    def __neg__(a): return TD2(-a.v, -a.g0, -a.g1, -a.h00, -a.h01, -a.h11)

    def __sub__(a, b): return a + (-b)

    def __rsub__(a, b): return (-a) + b

    def __mul__(a, b):
        if isinstance(b, TD2):
            return TD2(a.v * b.v, a.g0 * b.v + a.v * b.g0, a.g1 * b.v + a.v * b.g1,
                       a.h00 * b.v + 2 * (a.g0 * b.g0) + a.v * b.h00,
                       a.h01 * b.v + a.g0 * b.g1 + a.g1 * b.g0 + a.v * b.h01,
                       a.h11 * b.v + 2 * (a.g1 * b.g1) + a.v * b.h11)
        return TD2(a.v * b, a.g0 * b, a.g1 * b, a.h00 * b, a.h01 * b, a.h11 * b)
    __rmul__ = __mul__

    def inv(a):
        r = 1 / a.v; r2 = r * r; r3 = r2 * r
        return TD2(r, -a.g0 * r2, -a.g1 * r2, -a.h00 * r2 + 2 * (a.g0 * a.g0) * r3,
                   -a.h01 * r2 + 2 * (a.g0 * a.g1) * r3, -a.h11 * r2 + 2 * (a.g1 * a.g1) * r3)

    def __truediv__(a, b):
        if isinstance(b, TD2): return a * b.inv()
        return a * (1 / arb(b))

    def __rtruediv__(a, b): return a.inv() * b

    def __pow__(a, n):
        assert isinstance(n, int)
        if n < 0: return (a ** (-n)).inv()
        out = TD2(arb(1)); base = a
        while n:
            if n & 1: out = out * base
            base = base * base; n >>= 1
        return out

    def sqrt(a):
        s = a.v.sqrt(); i2 = 1 / (2 * s); i4 = 1 / (4 * s * s * s)
        return TD2(s, a.g0 * i2, a.g1 * i2, a.h00 * i2 - (a.g0 * a.g0) * i4, a.h01 * i2 - (a.g0 * a.g1) * i4,
                   a.h11 * i2 - (a.g1 * a.g1) * i4)


def comp(x, key):
    """component of a TD2 (or the value of a plain arb, whose derivatives are 0)"""
    if isinstance(x, TD2): return getattr(x, key)
    return x if key == 'v' else Z
