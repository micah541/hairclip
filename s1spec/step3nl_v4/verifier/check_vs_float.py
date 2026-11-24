"""Compare interval rhs (ADDED=0) with ../inner_closure_float2.py rhs at random points (degenerate boxes).  Also sanity-test
the structural claims: (i) monotonicity of the float rhs in |al|, (ii) POS(r)/r^2 nonincreasing along the radial D4,V4 direction."""
import os, sys, numpy as np
os.environ["ADD"] = ""
os.environ.update(CUBFINE="1", ALIN="spec", MU="0.08", EPSI="3e-5")
here = os.path.dirname(os.path.abspath(__file__)); src = os.path.dirname(here)
sys.path.insert(0, here)
import constants as K, ival, rhs_iv as R
K.TQ.update(d001=1.85e-2, d010=1.5e-3, d200=1.5e-2, d110=0.98, d300=3.26)   # float2 uses the pre-recertification T3 maxima
R._CC = None
cwd = os.getcwd(); os.chdir(src); sys.path.insert(0, src); sys.argv = ["x", "1e-5", "2e-3", "1e-6", "3e-6"]
import inner_closure_float2 as F
os.chdir(cwd)
rng = np.random.default_rng(1)
worst = 0
for trial in range(300):
    al = 10 ** rng.uniform(-6, -2.6); Ly = 10 ** rng.uniform(-14, -10.5)
    x = np.array([rng.uniform(0, 1), rng.uniform(0, 1) ** 3, rng.uniform(0, 1) ** 3, rng.uniform(0, 1), 10 ** rng.uniform(-3, 0.5)])
    cir = lambda t: ((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))
    cw, sw = cir(x[0]); ct, st = cir(x[1]); cI, sI = cir(x[2]); cu, su = cir(x[3]); s = np.sqrt(Ly)
    A = s * cw * ct; B = s * cw * st / np.sqrt(K.EPS); V0 = s * sw * cI / np.sqrt(K.MU); V2 = s * sw * sI / np.sqrt(K.MU * K.EPSI)
    r = s * x[4]; D4 = r * cu / np.sqrt(K.EPS); V4 = r * su / np.sqrt(K.MU * K.EPSI)
    f, _ = F.rhs(al, A, B, D4, V0, V2, V4)
    lo = x[None, :]; s_ = np.sqrt(Ly); sP = ival.cc(s_); rl, rh = R.rate_ub_lo(ival.cc(al), sP, sP, lo, lo)
    tl, th = rl * Ly, rh * Ly
    assert tl[0] <= f <= th[0] or abs(f - tl[0]) < 1e-14 * abs(f), (f, tl, th)
    worst = max(worst, (th[0] - tl[0]) / max(abs(f), 1e-300))
print("interval encloses float rhs at 300 random points; max relative width of point enclosure:", worst)
# monotone in al
for trial in range(200):
    Ly = 10 ** rng.uniform(-13, -11); x = rng.uniform(0, 1, 5); x[4] = 10 ** rng.uniform(-2, 0.5)
    cir = lambda t: ((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))
    cw, sw = cir(x[0]); ct, st = cir(x[1]); cI, sI = cir(x[2]); cu, su = cir(x[3]); s = np.sqrt(Ly)
    A = s * cw * ct; B = s * cw * st / np.sqrt(K.EPS); V0 = s * sw * cI / np.sqrt(K.MU); V2 = s * sw * sI / np.sqrt(K.MU * K.EPSI)
    r = s * x[4]; D4 = r * cu / np.sqrt(K.EPS); V4 = r * su / np.sqrt(K.MU * K.EPSI)
    vals = [F.rhs(a, A, B, D4, V0, V2, V4)[0] for a in (0, 1e-4, 5e-4, 1e-3, 2e-3)]
    assert all(np.diff(vals) >= -1e-18 * abs(vals[0]) - 1e-30), vals
print("float rhs monotone nondecreasing in |al| at 200 random points")
# radial POS/r^2 nonincreasing
bad = 0
for trial in range(200):
    al = 10 ** rng.uniform(-5, -2.6); Ly = 10 ** rng.uniform(-14, -10.5); x = rng.uniform(0, 1, 5); x[1:3] = x[1:3] ** 3
    out = []
    for q in (0.01, 0.1, 1, 10, 100, 1000):
        lo = np.array([[*x[:4], q]]); xh = R.xhat_ranges(lo, lo); sP = ival.cc(np.sqrt(Ly)); H, L_, N_ = R.rate_parts(ival.cc(al), sP, xh)
        out.append((H + L_).hi[0] / (q * q))
    if not all(np.diff(out) <= 1e-12 * np.array(out[:-1])): bad += 1
print("radial POS(r)/r^2 nonincreasing violations:", bad, "/ 200")
