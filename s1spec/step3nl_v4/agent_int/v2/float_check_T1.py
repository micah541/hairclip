"""float (mpmath, 40 digits) check of T1: c3 = g'''(0)/6 at points x, with k=m=p=0, q=sigma (ray), jets of psi3^, h2^ from the arb functions (midpoints)."""
import sys; sys.argv = ['x', '1e-3', '1']
import mpmath as mp
from int_common import *
import int_n3 as M
mp.mp.dps = 40
def Nmp(u, p, q, k, m):
    al = 1 - k * u[0]; Ji = 1 / mp.sqrt(al * al + u[1] * u[1]); Fc = 0
    for (i, j, l), f in M.CF.items():
        Fc += f(Ji, al, k, m, p, q, u[0], u[1]) * u[2] ** i * u[3] ** j * u[4] ** l
    drift = (p + u[0] - u[1] * q / al) / 4
    lin = -u[4] - k * k * u[2] - (5 * k * m + q / 4) * u[1] + (mp.mpf(1) / 4 - 3 * m * m - k * p) * u[0]
    return Fc + drift - lin
for x in (10.3, 11.0, 16.0):
    X = arb(x); xc = x
    y = [mp.mpf(str(v.mid().str(30, radius=False))) for v in jets_any(F, X, xc)]
    h = [mp.mpf(str(v.mid().str(30, radius=False))) for v in jets_any(Hf, X, xc)]
    sg = x + float(ST.mid()); q = mp.mpf(sg)
    g = lambda a: Nmp([a * y[j] + a * a * h[j] for j in range(5)], 0, q, 0, 0)
    d = mp.mpf('1e-5'); c3 = (g(d) - g(-d)) / (2 * d ** 3)
    # rigorous series coefficient at the same point (k=m=p=0 exactly)
    ys = jets_any(F, X, xc); hs = jets_any(Hf, X, xc)
    z = arb(0)
    tot, c3a, c0 = M.g3_bound(ys, hs, arb(0), arb(sg), arb(0), arb(0))
    print(f"x={x}: mp finite-diff c3 = {mp.nstr(c3, 10)};  arb series |c3| = {float(c3a.mid()):.10g}; c3 sign-agnostic ratio {abs(float(c3))/float(c3a.mid()):.8f}")
