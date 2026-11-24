"""Independent sampling check (NOT part of the proof) of the three V1 inequalities, using the float values of the interval objects' midpoints:
 (a) L(alpha,x) <= sum_k l_k(alpha) x_k            (linear majorant coded in unc.lk)
 (b) L - delta*N <= Phi(alpha)                      for random x of arbitrary scale
 (c) H - (1-delta) N <= r1 s^2                      for s <= s_up   (r1 from bb.bound_cell with R.set_mode(delta,True))
Usage: ADD=abcd FIX=1 THETA=0.3 python3 sanity_unc.py"""
import os, numpy as np
import constants as K, rhs_iv as R, unc, bb
from ival import cc, pt
theta = float(os.environ.get("THETA", 0.3)); delta = 0.2
rng = np.random.default_rng(3)
def parts(al, s, shape):
    # returns H, L, N (floats, midpoints of interval enclosures) at scale s, normalised shape xh (6,) -> physical
    xh = [pt(np.array([v]), np.array([v])) for v in shape]
    xs = [cc(s) * x for x in xh]
    H, L, N = R.rhs_pn(cc(al), *xs, floors=False)
    mid = lambda p: float(0.5 * (p.lo[0] + p.hi[0])) if np.ndim(p.lo) else float(p.lo)
    return mid(H), mid(L), mid(N), [float(s * v) for v in shape]
def shape_from(angles):
    w, t, tI, u, q = angles
    cir = lambda t: ((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))
    cw, sw = cir(w); ct, st = cir(t); cI, sI = cir(tI); cu, su = cir(u)
    return [cw * ct, cw * st / np.sqrt(K.EPS), q * cu / np.sqrt(K.EPS), sw * cI / np.sqrt(K.MU), sw * sI / np.sqrt(K.MU * K.EPSI), q * su / np.sqrt(K.MU * K.EPSI)]
R.set_theta(theta)
bad = [0, 0, 0]; n = 3000
a_hi, Lup = float(os.environ.get("A", 5e-4)), float(os.environ.get("L", 1.7e-13))
R.set_mode(0.0, False)
for i in range(n):
    al = a_hi * rng.uniform(0.05, 1); s = np.sqrt(Lup) * 10 ** rng.uniform(-3, 0)
    ang = [rng.uniform(0, 1), rng.uniform(0, 1) ** 3, rng.uniform(0, 1) ** 3, rng.uniform(0, 1), 10 ** rng.uniform(-3, 1.2)]
    sh = shape_from(ang)
    H, L, N, x = parts(al, s, sh)
    l = unc.lk(unc.A_(al), unc.Phi(a_hi, delta, theta)[1])
    x = dict(zip(("A", "B", "D4", "V0", "V2", "V4"), x))
    lin = sum(float(l[k].mid()) * x[k] for k in l)
    # NEG here is the theta-reduced N only if THETA set (it is): N = NEG - theta NEGV
    if L > lin * (1 + 1e-9): bad[0] += 1
    if L - delta * N > unc.Phi(a_hi, delta, theta)[0] * (1 + 1e-9): bad[1] += 1
print("violations of (a) L<=sum l_k x_k:", bad[0], "/", n, "; (b) L-delta N<=Phi:", bad[1], "/", n)
R.set_mode(delta, True)
r = bb.bound_cell(0.0, a_hi, Lup, Lup, tol=5e-4, rtol=0.005)
r1 = r["U"]; print("r1 =", r1, "(lb %.5f)" % r["lb"])
R.set_mode(0.0, False)
cnt = 0
for i in range(n):
    al = a_hi * rng.uniform(0.05, 1); s = np.sqrt(Lup) * 10 ** rng.uniform(-3, 0)
    ang = [rng.uniform(0, 1), rng.uniform(0, 1) ** 3, rng.uniform(0, 1) ** 3, rng.uniform(0, 1), 10 ** rng.uniform(-3, 1.2)]
    H, L, N, x = parts(al, s, shape_from(ang))
    if H - (1 - delta) * N > r1 * s * s * (1 + 1e-9) + 1e-300: cnt += 1
print("violations of (c) H-(1-delta)N<=r1 s^2:", cnt, "/", n)
