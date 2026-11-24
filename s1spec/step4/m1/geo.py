"""S1 geometry on an extended grid + C4 psi3/h2 jets (cached)."""
import sys, os, numpy as np, time
sys.path[:0] = ['/home/claude/s1spec/step4/feas', '/home/claude/s1spec/step4/agent_handoff2', '/home/claude/s1spec/step4', '/home/claude/s1spec',
                '/home/claude/s1spec/step3nl/agent_int/v3', '/home/claude/s1spec/step3nl']
import chart as CH
H = CH.H
from measure import STIP
HX = 0.01; XMAX = 64.0
xg = np.arange(0, XMAX + HX / 2, HX); NX = len(xg)
n0 = len(H.xs)
tauG = np.empty(NX, complex); nuG = np.empty(NX, complex); G1G = np.empty(NX, complex); kG = np.zeros(NX)
tauG[:n0] = H.tau1; nuG[:n0] = H.nu1; G1G[:n0] = H.G1; kG[:n0] = H.M.k1
tinf = H.tau1[-1]
tauG[n0:] = tinf; nuG[n0:] = 1j * tinf; G1G[n0:] = (xg[n0:] + STIP) * tinf
RAYERR = abs(H.G1[-1] - (H.xs[-1] + STIP) * tinf)       # p(x=40) check: S1 is a ray beyond sigma~22
THINF = float(np.angle(tinf))
CACHE = '/home/claude/s1spec/step4/m1/c4_cache.npz'

def _c4_arrays():
    import c4fun as C
    from flint import arb
    out = np.zeros((4, NX))
    for i, x in enumerate(xg):
        a = C.PSI.jets_any(arb(float(x)), float(x), 2); b = C.H2.jets_any(arb(float(x)), float(x), 2)
        out[:, i] = [float(a[0].mid()), float(a[1].mid()), float(b[0].mid()), float(b[1].mid())]
    return out
if os.path.exists(CACHE): _c = np.load(CACHE)['a']
else:
    t = time.time(); _c = _c4_arrays(); np.savez(CACHE, a=_c); print("c4 cache built", time.time() - t)
PSI0, PSI1, H20, H21 = _c        # psi3, psi3', h2, h2' on xg

def c4_first(y):
    """psi3'(x), h2'(x) at y = sigma (x = y - STIP), any y (arb jets, extension beyond x=28 = power-law asymptotics of the certified object)."""
    import c4fun as C
    from flint import arb
    y = np.atleast_1d(y); o1 = np.empty(len(y)); o2 = np.empty(len(y))
    for i, yy in enumerate(y):
        x = float(yy - STIP)
        o1[i] = float(C.PSI.jets_any(arb(x), x, 2)[1].mid()); o2[i] = float(C.H2.jets_any(arb(x), x, 2)[1].mid())
    return o1, o2

def smooth(u):
    uu = np.clip(u, 0, 1); return uu**4 * (35 - 84 * uu + 70 * uu**2 - 20 * uu**3)
sig = xg + STIP
rho = np.exp(-0.25 * (1 + xg**2) ** (2 / 3))
chi32 = smooth((32 - sig) / 3)
zeta = smooth((xg - 10.25) / 6) * smooth((60 - xg) / 3)
# W = exp(-psi), psi' = (1-s) 4x/(1+x^2) + s phi', s = S((x-30.25)/8)
_s = smooth((xg - 30.25) / 8); _dpsi = (1 - _s) * 4 * xg / (1 + xg**2) + _s * xg / (3 * (1 + xg**2) ** (1 / 3))
_psi = np.concatenate([[0], np.cumsum(0.5 * (_dpsi[1:] + _dpsi[:-1])) * HX]); Wt = np.exp(-_psi); Wt = np.where(xg <= 30.25, (1 + xg**2) ** -2.0, Wt)
if __name__ == '__main__':
    print("ray err at x=40", RAYERR, "THINF", THINF, "W(30)", Wt[int(30 / HX)], (1 + 30.0**2) ** -2, "zeta(60)", zeta[-1], len(xg))
    # check U ~ 1.3056*y^1.0246 (psi3' ~ 1.4 y^1.0246)
    for y in (50., 100., 250.):
        a, b = c4_first(y); print(y, a[0], a[0] / y**1.0246, b[0])
