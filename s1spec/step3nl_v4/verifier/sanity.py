"""Independent sanity checks of the rigorous cell bounds (do not replace the proof; they would catch coding errors):
 (1) U >= lb and every random/float-grid shape rate (computed with the float code ../inner_closure_float2.py when ADD="") is <= U;
 (2) scale-invariance lemma used by the tail:  rate_s(shape; s alpha, s^2 Ly) <= rate_s(shape; alpha, Ly) (floors=False), s<=1."""
import os, sys, time
import numpy as np
import constants as K, bb, rhs_iv as R
from ival import *
rng = np.random.default_rng(7)
print("ADD =", sorted(K.ADD))
cells = [(0.0, 1e-12, 1e-12), (1e-3, 1e-12, 1.02e-12), (2e-3, 2e-12, 2.1e-12), (2e-4, 1.5e-11, 1.55e-11), (5e-5, 1e-14, 1.1e-14)]
fl = None
if not K.ADD:
    os.environ.update(CUBFINE="1", ALIN="spec", MU=str(K.MU), EPSI=str(K.EPSI))
    here = os.getcwd(); src = os.path.dirname(here); os.chdir(src); sys.path.insert(0, src); sys.argv = ["x", "1e-5", "2e-3", "1e-6", "3e-6"]
    import inner_closure_float2 as fl; os.chdir(here)
for a, L0, L1 in cells:
    t = time.time(); r = bb.bound_cell(0.0, a, L0, L1, tol=5e-4, rtol=0.005)
    assert r["U"] >= r["lb"]
    # random shapes at random Ly in the cell, alpha in [0,a]
    n = 60000
    y = rng.uniform([-6] * 4 + [-3], [0] * 4 + [np.log10(r["qstar"])], (n, 5)); p = 10 ** y; p[:, :4] = np.minimum(p[:, :4], 1)
    p[:20000] = np.clip(np.array(r["best"]) * 10 ** rng.normal(0, 0.15, (20000, 5)), 0, None); p[:20000, :4] = np.minimum(p[:20000, :4], 1)
    worst = -np.inf
    for Ly in (L0, L1, np.sqrt(L0 * L1)):
        for al in (a, 0.5 * a):
            sP = cc(np.sqrt(Ly)); rl, rh = R.rate_ub_lo(cc(al) if al > 0 else pt(np.float64(0), np.float64(0)), sP, sP, p, p)
            worst = max(worst, rl.max())
    msg = f"cell alpha<={a:g} Ly in [{L0:g},{L1:g}]: U={r['U']:+.5f} lb={r['lb']:+.5f} random-shape max={worst:+.5f}"
    assert worst <= r["U"] + 1e-12, msg
    if fl is not None:
        gm = -np.inf
        for Ly in (L0, L1):
            b, pb, i = fl.worst(a, Ly); gm = max(gm, b / Ly)
        msg += f" float-grid max={gm:+.5f}"
        assert gm <= r["U"] + 1e-12, msg
    print(msg, f"({time.time()-t:.1f}s)  OK", flush=True)
# scale lemma
bad = 0
for trial in range(4000):
    a = 10 ** rng.uniform(-6, -3); L = 10 ** rng.uniform(-14, -11); s = rng.uniform(0.05, 1.0)
    x = rng.uniform(0, 1, 5); x[1:3] **= 3; x[4] = 10 ** rng.uniform(-3, 1.5)
    def rate(al, Lv):
        sP = cc(np.sqrt(Lv)); return R.rate_ub_lo(cc(al), sP, sP, x[None, :], x[None, :], floors=False)[0][0]
    if rate(s * a, s * s * L) > rate(a, L) + 1e-12: bad += 1
print("scale-invariance violations (floors=False):", bad, "/ 4000")
