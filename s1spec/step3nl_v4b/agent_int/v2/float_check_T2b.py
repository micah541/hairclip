"""float check of the h2 residual with the (float) tight-tail profile: R(x) = -h4 - k^2 h2'' - (5km+q/4) h' + (1/4 - 2lamt - 3m^2 - kp) h + Q2(psi3^; p,q,k,m)"""
import sys; sys.argv = ['x', '1']
import numpy as np
from int_common import *
import int_res as M
from flint import ctx
def full(x):
    with ctx.workprec(300):
        X = arb(x); xc = x
        h = jets_any(Hf, X, xc); y = jets_any(F, X, xc)
        sig = arb(x) + ST.mid()
        p, q, k, m = [t.mid() for t in M.tail_state(sig)]
        a2, a1, a0 = M.coef0(p, q, k, m)
        R0 = -h[4] - a2 * h[2] - a1 * h[1] + a0 * h[0] + M.q2('Q2', y[:5], p, q, k, m)
    return float(R0.mid()), float(k)
bl = M.boxlist()
for x in (10.3, 11.0, 14.0):
    for sa, sb in bl:
        b = make_box(sa, sb)
        if b.xa <= x <= b.xb: break
    r, kk = full(x); bnd = M.h2_res_Bt(b)[0]
    print(f"x={x}: float full residual {r:+.3e} (k={kk:.2e}); rigorous box bound {float(bnd.upper().mid()):.3e}")
xs = np.arange(10.3, 22.2, 0.01); vals = []
for x in xs:
    r, _ = full(float(x)); vals.append(r * r * float(W_up("2", float(x)).mid()))
print("float ||R||_W on [12, 22.2]:", float(np.trapezoid(vals, xs)) ** .5, "(rigorous region [12,20] norm 1.8e-10 .. see json)")
xs20 = xs[xs <= 20.0]; print("float ||R||_W on [12,20]:", float(np.trapezoid(np.array(vals)[:len(xs20)], xs20)) ** .5)
