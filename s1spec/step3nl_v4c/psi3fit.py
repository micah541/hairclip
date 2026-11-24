"""Chebyshev pieces for psi3^ on x = sigma - sigma_tip in [0, 34] and its float eigen-residual (untrusted)."""
import sys, json, numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.integrate import solve_ivp
sys.path.insert(0, '/home/claude/s1spec')
import evans as E
import psi3hat2 as H                                  # reruns the inward construction (quick)
c = H.c; ST = E.STIP
out = solve_ivp(H.rhs, [22.0, 32.0], sum(c[i] * H.P[i] for i in range(2)), method='DOP853', rtol=1e-12, atol=1e-30, dense_output=True)
SM = 10.0
jm = sum(c[i] * H.sols[i].sol(SM) for i in range(2))
inner = solve_ivp(H.rhs, [SM, ST], jm, method='DOP853', rtol=1e-13, atol=1e-16, dense_output=True)
print("combined jet at the tip:", inner.y[:, -1])
def U(s, d=0):
    s = np.atleast_1d(s)
    if np.all(s <= SM): return inner.sol(s)[d]
    s = np.atleast_1d(s); r = np.empty(len(s)); inn = s <= 22.0
    r[inn] = sum(c[i] * H.sols[i].sol(s[inn])[d] for i in range(2)) if inn.any() else r[inn]
    if (~inn).any():
        A = 1 - 4 * H.LT; ss = s[~inn]; f = 1.0
        for j in range(d): f = f * (A - j)
        r[~inn] = c[0] * f * ss ** (A - d)          # polynomial branch (decaying part negligible beyond sigma = 22)
    return r
pieces = []; worst = 0.0; res2 = 0.0
for lo in np.arange(0.0, 34.0, 0.5):
    t = 0.5 * (1 - np.cos(np.pi * (np.arange(64) + 0.5) / 64)) * 0.5 + lo   # Chebyshev points on [lo, lo+0.5]
    cf = Ch.chebfit(4 * (t - lo) - 1, U(t + ST), 30)
    pieces.append({"lo": float(lo), "h": 0.5, "coef": cf.tolist()})
    tt = np.linspace(lo, lo + 0.5, 41); z = 4 * (tt - lo) - 1
    d = [Ch.chebval(z, Ch.chebder(cf, n)) * 4.0**n for n in range(5)]
    pr = np.array([E.PROF(s) if s <= 22 else [0.0, s, 0.0, 0.0] for s in tt + ST]).T
    p, q, k, m = pr
    res = -d[4] - k * k * d[2] - (5 * k * m + q / 4) * d[1] + (0.25 - 3 * m * m - k * p - H.LT) * d[0]
    rho = np.exp(-0.25 * (1 + tt**2) ** (2 / 3))
    res2 += 2 * np.trapezoid(rho * res**2, tt); worst = max(worst, np.max(np.abs(res)) * np.sqrt(rho.max()))
json.dump({"stip": float(ST), "lamt": H.LT, "pieces": pieces}, open('psi3hat_cheb.json', 'w'))
print(f"float eigen-residual of the Chebyshev psi3^:  ||(L - lamt) psi3^||_rho ~ {np.sqrt(res2):.2e}   (max weighted {worst:.2e});  psi3^(10) = {U(10 + ST)[0]:.4f}")
