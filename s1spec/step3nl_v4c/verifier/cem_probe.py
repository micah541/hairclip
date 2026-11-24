"""Diagnostic (NOT a proof of an upper bound): rigorous LOWER bounds of sup rate via cross-entropy search + certified point evaluation.
Used to show where the closure fails (a positive lower bound refutes closure rigorously)."""
import sys, numpy as np, os
import bb, rhs_iv as R
from ival import *
def lower(a, L, floors=True, iters=40, n=6000, seed=1):
    al = cc(a) if a > 0 else pt(np.float64(0), np.float64(0)); sl, sh = bb._s_pts(L, L)
    qs = bb.find_qstar(al, sl, sh, floors)
    rng = np.random.default_rng(seed)
    lg_lo = np.array([-6, -6, -6, -6, -4.0]); lg_hi = np.array([0, 0, 0, 0, np.log10(qs)])
    mu = 0.5 * (lg_lo + lg_hi); sd = 0.5 * (lg_hi - lg_lo); best = None; lb = -np.inf
    for it in range(iters):
        y = np.clip(mu + sd * rng.standard_normal((n, 5)), lg_lo, lg_hi)
        if best is not None: y[:200] = np.clip(np.log10(np.maximum(best, 1e-12)) + 0.05 * sd * rng.standard_normal((200, 5)), lg_lo, lg_hi)
        p = 10 ** y; p[:, :4] = np.minimum(p[:, :4], 1.0); p[:, 4] = np.minimum(p[:, 4], qs)
        rl, rh = bb.R.rate_ub_lo(al, sl, sh, p, p, floors)
        o = np.argsort(-rl)[:150]
        if rl[o[0]] > lb: lb = rl[o[0]]; best = p[o[0]].copy()
        mu = y[o].mean(axis=0); sd = np.maximum(y[o].std(axis=0), 0.02)
    return lb, best
if __name__ == "__main__":
    for L in (1e-16, 1e-14, 1e-12):
        for a in (0.0, 1e-4):
            lb, b = lower(a, L, floors=False)
            print(f"ADD={os.environ.get('ADD','abcd')!r} alpha={a:g} Ly={L:g}: rigorous lower bound of max rate = {lb:+.5f} at {np.array2string(b, precision=3)}")
