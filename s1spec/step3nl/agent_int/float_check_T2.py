"""float check of T2: high-precision pointwise ray residuals (k=m=p=0) of psi3^ and h2^ vs the rigorous box bounds, and a float W-norm of the ray parts."""
import sys; sys.argv = ['x', '1', 'W']
import mpmath as mp, json
from int_common import *
import int_res as M
from flint import ctx
mp.mp.dps = 40
def tomp(v): return mp.mpf(v.mid().str(40, radius=False))
def ray(MOD, lam, x):
    with ctx.workprec(300):
        X = arb(x); u = MOD_jets(MOD, X)
    s = x + float(STIP.mid()); u = [tomp(v) for v in u]
    return -u[4] - mp.mpf(x + tomp(STIP)) / 4 * u[1] + (mp.mpf(1) / 4 - tomp(lam)) * u[0]
def MOD_jets(MOD, X):
    xc = float(X.mid()); return jets_any(MOD, X, xc)
# pointwise vs box bound
bl = M.boxlist()
for x in (12.5, 15.0, 19.99, 24.0, 27.99, 28.01, 33.3, 45.0, 59.0):
    for sa, sb in bl:
        b = make_box(sa, sb)
        if b.xa <= x <= b.xb: break
    rp = ray(F, F.LAM, x); rh = ray(Hf, Hf.LAM, x)
    bp, _ = M.psi_res(b, b.X, b.xc); ph_ = (M.h2_res_Bt(b)[0] if sb <= M.SB_END else M.h2_res_B(b, b.X, b.xc)[0])
    print(f"x={x}: float ray residual psi3 {float(rp):+.3e} (rigorous box bound {float(bp.upper().mid()):.3e}); h2 {float(rh):+.3e} (box bound {float(ph_.upper().mid()):.3e})")
# float W-norm of the ray parts on [12,40],[12,60] (trapezoid, step 0.01)
import numpy as np
def Wf(x): 
    return float(mp.e ** (-mp.mpf(str(psi_value(float(x)).mid().str(30, radius=False)))))
for MOD, lam, name in ((F, F.LAM, 'psi3'), (Hf, Hf.LAM, 'h2')):
    xs = np.arange(12.0, 60.0 + 1e-9, 0.02); acc = {40: 0.0, 60: 0.0}
    vals = []
    for x in xs:
        r = float(ray(MOD, lam, float(x))); vals.append(r * r * Wf(x))
    vals = np.array(vals)
    for R in (40, 60):
        msk = xs <= R + 1e-9; acc[R] = float(np.trapezoid(vals[msk], xs[msk]))
    print(name, "float ray-only ||r||_W: R=40", acc[40] ** .5, " R=60", acc[60] ** .5)
