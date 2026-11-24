"""float (mpmath) check of T3 at x = 45 (box region [40,60]) with k=m=p=0, q=sigma: Taylor coefficients of D(w2,w3,w4) = N0(u0+w) - N0(u0) - L w
 (w0 = w1 = 0) vs the rigorous region sups;  and sup|u0'|, |u0''| at x = 59.99."""
import sys; sys.argv = ['x', '1', '1']
import mpmath as mp, json
from int_common import *
import int_n3 as M
mp.mp.dps = 40
def tomp(v): return mp.mpf(v.mid().str(40, radius=False))
al = mp.mpf('1e-3')
def Dfun(x):
    X = arb(x); y = [tomp(v) for v in jets_any(F, X, x)]; h = [tomp(v) for v in jets_any(Hf, X, x)]
    u = [al * y[j] + al * al * h[j] for j in range(5)]; q = mp.mpf(x + float(ST.mid()))
    def N(w2, w3, w4):
        uu = list(u); uu[2] += w2; uu[3] += w3; uu[4] += w4
        a = 1; Ji = 1 / mp.sqrt(1 + uu[1] ** 2); Fc = 0
        for (i, j, l), f in M.CF.items(): Fc += f(Ji, 1, 0, 0, 0, q, uu[0], uu[1]) * uu[2] ** i * uu[3] ** j * uu[4] ** l
        return Fc + (uu[0] - uu[1] * q) / 4 - (-uu[4] - q / 4 * uu[1] + uu[0] / 4)   # N0 = F + drift - L u (k=m=p=0)
    N00 = N(0, 0, 0)
    return lambda w2, w3, w4: N(w2, w3, w4) - N00 - (-w4), u
x = 45.0
f, u = Dfun(x)
print("u0' =", mp.nstr(u[1], 8), " u0'' =", mp.nstr(u[2], 8))
c = {}
c['(0,0,1)'] = mp.diff(lambda t: f(0, 0, t), 0, 1)
c['(0,1,0)'] = mp.diff(lambda t: f(0, t, 0), 0, 1)
c['(1,0,0)'] = mp.diff(lambda t: f(t, 0, 0), 0, 1)
c['(2,0,0)'] = mp.diff(lambda t: f(t, 0, 0), 0, 2) / 2
c['(3,0,0)'] = mp.diff(lambda t: f(t, 0, 0), 0, 3) / 6
c['(1,1,0)'] = mp.diff(lambda s, t: f(s, t, 0), (0, 0), (1, 1))
for kk, v in c.items(): print(kk, mp.nstr(v, 8))
print("1 - J^-4 at u0' :", mp.nstr(1 - (1 + u[1] ** 2) ** -2, 8))
