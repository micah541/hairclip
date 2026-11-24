"""psi3^ with Newton-refined (beta, lambda): jet at sigma=10 is J1 + beta*J2 (outer solutions from sigma=22),
inward integration to the tip, even conditions u'(tip) = u'''(tip) = 0.  Then Chebyshev pieces + float residual."""
import sys, json, numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.integrate import solve_ivp
sys.path.insert(0, '/home/claude/s1spec')
import evans as E
ST, S0, SM = E.STIP, 22.0, 10.0
def mk_rhs(lam):
    def rhs(s, U):
        p, q, k, m = E.PROF(s); u, u1, u2, u3 = U
        return [u1, u2, u3, -k * k * u2 - (5 * k * m + q / 4) * u1 + (0.25 - 3 * m * m - k * p - lam) * u]
    return rhs
def outer(lam):
    A = 1 - 4 * lam; z = -(S0 / 4) ** (1 / 3)
    P1 = [S0**A, A * S0**(A - 1), A * (A - 1) * S0**(A - 2), A * (A - 1) * (A - 2) * S0**(A - 3)]
    P2 = [1.0, z, z * z, z**3]
    o = [solve_ivp(mk_rhs(lam), [S0, SM], P, method='DOP853', rtol=1e-13, atol=1e-30, dense_output=True) for P in (P1, P2)]
    n2 = np.linalg.norm(o[1].y[:, -1]); return o, n2
def inner(beta, lam):
    o, n2 = outer(lam)
    J = o[0].y[:, -1] + beta * o[1].y[:, -1] / n2
    return solve_ivp(mk_rhs(lam), [SM, ST], J, method='DOP853', rtol=1e-13, atol=1e-30, dense_output=True), o, n2
def F(v):
    sol, _, _ = inner(*v); y = sol.y[:, -1]; return np.array([y[1], y[3]]) / abs(y[0])
import psi3hat2 as H0
_o, _n2 = outer(-0.2561532224577463)
v = np.array([H0.c[1] / H0.c[0] * _n2, -0.2561532224577463])
print("initial beta from the SVD solution:", v[0])
for it in range(10):
    f0 = F(v); Jm = np.zeros((2, 2))
    for j, h in enumerate([1e-5 * max(1.0, abs(v[0])), 1e-6]):
        e = np.zeros(2); e[j] = h; Jm[:, j] = (F(v + e) - f0) / h
    v = v - np.linalg.lstsq(Jm, f0, rcond=None)[0]; print(it, v, np.abs(f0).max(), flush=True)
    if np.max(np.abs(f0)) < 1e-13: break
beta, lam = v; sol, o, n2 = inner(beta, lam); u0 = sol.y[0, -1]
print(f"beta = {beta:.15e}, lambda = {lam:.15f} (certified lamt -0.2561532224577463); tip jet / u0:", sol.y[:, -1] / u0)
A = 1 - 4 * lam
def U(s, d=0):
    s = np.atleast_1d(s).astype(float); r = np.empty(len(s))
    a = s <= SM; b = (s > SM) & (s <= S0); c = s > S0
    if a.any(): r[a] = sol.sol(s[a])[d]
    if b.any(): r[b] = o[0].sol(s[b])[d] + beta * o[1].sol(s[b])[d] / n2
    if c.any():
        f = 1.0
        for j in range(d): f *= (A - j)
        r[c] = f * s[c] ** (A - d)
    return r / u0
pieces = []; res2 = 0.0; worst = (0, 0)
for lo in np.arange(0.0, 34.0, 0.5):
    t = lo + 0.25 * (1 - np.cos(np.pi * (np.arange(64) + 0.5) / 64))
    cf = Ch.chebfit(4 * (t - lo) - 1, U(t + ST), 30); pieces.append({"lo": float(lo), "h": 0.5, "coef": cf.tolist()})
    tt = np.linspace(lo, lo + 0.5, 41); z = 4 * (tt - lo) - 1
    d = [Ch.chebval(z, Ch.chebder(cf, n)) * 4.0**n for n in range(5)]
    p, q, k, m = np.array([E.PROF(s) if s <= 22 else [0.0, s, 0.0, 0.0] for s in tt + ST]).T
    res = -d[4] - k * k * d[2] - (5 * k * m + q / 4) * d[1] + (0.25 - 3 * m * m - k * p - lam) * d[0]
    rho = np.exp(-0.25 * (1 + tt**2) ** (2 / 3)); res2 += 2 * np.trapezoid(rho * res**2, tt)
    w = np.max(np.abs(res) * np.sqrt(rho)); worst = max(worst, (w, lo))
json.dump({"stip": float(ST), "lam": float(lam), "pieces": pieces}, open('psi3hat_cheb2.json', 'w'))
print(f"float residual ||(L - lam) psi3^||_rho ~ {np.sqrt(res2):.2e}, worst weighted {worst[0]:.2e} at x = {worst[1]}; psi3^(2,5,10) =", np.round(U(np.array([2.0, 5, 10]) + ST), 4))
