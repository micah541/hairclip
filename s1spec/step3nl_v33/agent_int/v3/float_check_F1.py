"""float check of F1a: d/dx of c3(x) (alpha^3 coefficient of g, k=m=p=0, q=sigma) by central differences in x of the mp finite-difference c3, vs the dual-number series coefficient."""
import sys; sys.argv = ['x', '1e-3', '1']
import mpmath as mp
from int_common import *
import int_n3d as M
from c4fun import PSI, H2
mp.mp.dps = 50
def tomp(v): return mp.mpf(v.mid().str(45, radius=False))
def c3_at(x):
    X = arb(x); y = [tomp(v) for v in PSI.jets_any(X, x, 6)]; h = [tomp(v) for v in H2.jets_any(X, x, 6)]
    q = mp.mpf(x) + tomp(ST)
    def N(u):
        Ji = 1 / mp.sqrt(1 + u[1] ** 2); Fc = 0
        for (i, j, l), f in M.CF.items():
            Fc += f(Ji, 1, 0, 0, 0, q, u[0], u[1]) * u[2] ** i * u[3] ** j * u[4] ** l
        return Fc + (u[0] - u[1] * q) / 4 + u[4] + q / 4 * u[1] - u[0] / 4
    g = lambda a: N([a * y[j] + a * a * h[j] for j in range(5)])
    d = mp.mpf('1e-6'); return (g(d) - g(-d)) / (2 * d ** 3)
for x in (12.0, 15.0, 25.5, 35.0):
    dl = 1e-4
    fd = (c3_at(x + dl) - c3_at(x - dl)) / (2 * dl)
    Xb = arb(x); y = PSI.jets_any(Xb, x, 6); h = H2.jets_any(Xb, x, 6)
    z = arb(0); tot, d3, totv, cv, cd = M.gd_bound(y, h, arb(0), arb(x) + ST.mid(), arb(0), arb(0))
    print(f"x={x}: float d/dx c3 = {mp.nstr(fd, 9)};  dual series c3' = {float(cd[3].mid()):.9g}")
