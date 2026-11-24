"""psi3^ by Chebyshev collocation on x in [0, X] (even BCs u'=u'''=0 at 0, free end u''=u'''=0 at X): no shooting."""
import sys, json, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
X, N = 36.0, int(sys.argv[1]) if len(sys.argv) > 1 else 160
j = np.arange(N + 1); t = np.cos(np.pi * j / N); x = (1 - t) * X / 2           # x in [0, X]
c = np.ones(N + 1); c[0] = c[-1] = 2; c *= (-1) ** j
Tm = np.tile(t, (N + 1, 1)).T; dT = Tm - Tm.T
D = np.outer(c, 1 / c) / (dT + np.eye(N + 1)); D -= np.diag(D.sum(1)); D *= -2 / X   # d/dx
D2, D3, D4 = D @ D, D @ D @ D, D @ D @ D @ D
P = np.array([E.PROF(s + E.STIP) if s + E.STIP <= 22 else [0.0, s + E.STIP, 0.0, 0.0] for s in x]).T
p, q, k, m = P
L = -D4 - np.diag(k * k) @ D2 - np.diag(5 * k * m + q / 4) @ D + np.diag(0.25 - 3 * m * m - k * p)
B = np.eye(N + 1); A = L.copy()
i0, iX = 0, N                                    # x = 0 is node 0 (t = 1), x = X is node N
for r, row in zip([i0, i0 + 1, iX, iX - 1], [D[i0], D3[i0], D2[iX], D3[iX]]):
    A[r] = row; B[r] = 0
sh = -0.2561532224577463; v = np.ones(N + 1); lu = sl.lu_factor(A - sh * B)
for _ in range(30):
    w = sl.lu_solve(lu, B @ v); v = w / np.linalg.norm(w)
lam = sh + (v @ (B @ v)) / (v @ sl.lu_solve(lu, B @ v)) if False else None
Lv = L @ v; lam = (v[2:-2] @ Lv[2:-2]) / (v[2:-2] @ v[2:-2]); u = v / v[i0]
res = L @ u - lam * u; rho = np.exp(-0.25 * (1 + x**2) ** (2 / 3)); inner = slice(2, N - 1)
print(f"N={N}: lambda = {lam:.13f} (certified -0.2561532224577463); u(2,5,10) =", np.round(np.interp([2, 5, 10], x, u), 4),
      f"; max weighted residual at nodes = {np.max(np.abs(res[inner]) * np.sqrt(rho[inner])):.1e}")
json.dump({"X": X, "N": N, "lam": lam, "u": u.tolist()}, open(f'psi3coll_{N}.json', 'w'))
from numpy.polynomial import chebyshev as Ch
cf = Ch.chebfit(t, u, N)                                   # exact interpolant in t
xf = np.linspace(0, 34.25, 6000); tf = 1 - 2 * xf / X
dd = [Ch.chebval(tf, Ch.chebder(cf, n)) * (-2 / X) ** n for n in range(5)]
Pf = np.array([E.PROF(s + E.STIP) if s + E.STIP <= 22 else [0.0, s + E.STIP, 0.0, 0.0] for s in xf]).T
pf, qf, kf, mf = Pf
rf = -dd[4] - kf * kf * dd[2] - (5 * kf * mf + qf / 4) * dd[1] + (0.25 - 3 * mf * mf - kf * pf - lam) * dd[0]
rhof = np.exp(-0.25 * (1 + xf**2) ** (2 / 3))
print(f"fine grid on [0, 34.25]: ||(L - lam) psi3^||_rho = {np.sqrt(2 * np.trapezoid(rhof * rf**2, xf)):.2e}, max weighted {np.max(np.abs(rf) * np.sqrt(rhof)):.2e}; u'(0) = {dd[1][0]:.1e}, u'''(0) = {dd[3][0]:.1e}")
