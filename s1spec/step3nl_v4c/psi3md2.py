"""Multi-domain Chebyshev collocation for psi3^ (C^3 matching), even BCs at 0, free end at X."""
import sys, json, numpy as np, scipy.linalg as sl
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
br = np.array([0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 14, 17, 20, 23, 26, 28.0])
LT = -0.2561532224577463
def cheb(N):
    j = np.arange(N + 1); t = np.cos(np.pi * j / N); c = np.ones(N + 1); c[0] = c[-1] = 2; c *= (-1) ** j
    T = np.tile(t, (N + 1, 1)).T; D = np.outer(c, 1 / c) / (T - T.T + np.eye(N + 1)); D -= np.diag(D.sum(1)); return t, D
t, Dt = cheb(n); Ne = len(br) - 1; M = Ne * (n + 1)
A = np.zeros((M, M)); B = np.zeros((M, M)); X = []; Ds = []
for e in range(Ne):
    a, b = br[e], br[e + 1]; x = a + (1 - t) * (b - a) / 2; D = Dt * (-2 / (b - a))   # node 0 at x = a
    D2 = D @ D; D3 = D2 @ D; D4 = D3 @ D; X.append(x); Ds.append((np.eye(n + 1), D, D2, D3))
    p, q, k, m = np.array([E.PROF(s + E.STIP) if s + E.STIP <= 22 else [0.0, s + E.STIP, 0.0, 0.0] for s in x]).T
    L = -D4 - np.diag(k * k) @ D2 - np.diag(5 * k * m + q / 4) @ D + np.diag(0.25 - 3 * m * m - k * p)
    o = e * (n + 1); A[o:o + n + 1, o:o + n + 1] = L; B[o:o + n + 1, o:o + n + 1] = np.eye(n + 1)
for e in range(Ne):
    o = e * (n + 1); I, D, D2, D3 = Ds[e]
    rows_l = [o, o + 1]; rows_r = [o + n, o + n - 1]
    for r in rows_l + rows_r: A[r] = 0; B[r] = 0
    if e == 0:
        A[o, o:o + n + 1] = D[0]; A[o + 1, o:o + n + 1] = D3[0]
    else:                                          # continuity of u, u' with element e-1 (its last node)
        op = (e - 1) * (n + 1); Ip, Dp, D2p, D3p = Ds[e - 1]
        A[o, o:o + n + 1] = I[0]; A[o, op:op + n + 1] = -Ip[n]
        A[o + 1, o:o + n + 1] = D[0]; A[o + 1, op:op + n + 1] = -Dp[n]
    if e == Ne - 1:                                # admissibility at x = X: jet in span{J_poly, J_dec}
        Xs = br[-1] + E.STIP; Aex = 1 - 4 * LT
        aj = [1.0]
        for jj in range(1, 4): aj.append(aj[-1] * (Aex - 4*jj + 4) * (Aex - 4*jj + 3) * (Aex - 4*jj + 2) * (Aex - 4*jj + 1) / jj)
        Jp = np.zeros(4)
        for jj in range(4):
            ex = Aex - 4 * jj; f = 1.0
            for d in range(4): Jp[d] += aj[jj] * f * Xs ** (ex - d); f *= (ex - d)
        z = -(Xs / 4) ** (1 / 3); z1 = -(1 / 12) * (Xs / 4) ** (-2 / 3); z2 = (1 / 72) * (Xs / 4) ** (-5 / 3)
        Jd = np.array([1.0, z, z1 + z * z, z2 + 3 * z * z1 + z**3])
        Nl = sl.null_space(np.vstack([Jp, Jd]))      # 4x2: covectors annihilating the admissible plane
        R4 = np.vstack([I[n], D[n], D2[n], D3[n]])
        A[o + n, o:o + n + 1] = Nl[:, 0] @ R4; A[o + n - 1, o:o + n + 1] = Nl[:, 1] @ R4
        JP, JD = Jp, Jd
    else:                                          # continuity of u'', u''' with element e+1 (its first node)
        on = (e + 1) * (n + 1); In, Dn, D2n, D3n = Ds[e + 1]
        A[o + n, o:o + n + 1] = D2[n]; A[o + n, on:on + n + 1] = -D2n[0]
        A[o + n - 1, o:o + n + 1] = D3[n]; A[o + n - 1, on:on + n + 1] = -D3n[0]
lu = sl.lu_factor(A - LT * B); v = np.ones(M)
for _ in range(40):
    w = sl.lu_solve(lu, B @ v); v = w / np.linalg.norm(w)
Av, Bv = A @ v, B @ v; mask = np.abs(Bv) > 0; lam = (Bv[mask] @ Av[mask]) / (Bv[mask] @ Bv[mask])
u = v / v[0]; pieces = []; res2 = 0.0; worst = 0.0
for e in range(Ne):
    a, b = br[e], br[e + 1]; ue = u[e * (n + 1):(e + 1) * (n + 1)]; cf = Ch.chebfit(t, ue, n)
    pieces.append({"a": float(a), "b": float(b), "coef": cf.tolist()})
    xf = np.linspace(a, b, 400); tf = 1 - 2 * (xf - a) / (b - a)
    dd = [Ch.chebval(tf, Ch.chebder(cf, j)) * (-2 / (b - a)) ** j for j in range(5)]
    p, q, k, m = np.array([E.PROF(s + E.STIP) if s + E.STIP <= 22 else [0.0, s + E.STIP, 0.0, 0.0] for s in xf]).T
    r = -dd[4] - k * k * dd[2] - (5 * k * m + q / 4) * dd[1] + (0.25 - 3 * m * m - k * p - lam) * dd[0]
    rho = np.exp(-0.25 * (1 + xf**2) ** (2 / 3)); sel = xf <= 28.0
    res2 += 2 * np.trapezoid((rho * r**2)[sel], xf[sel]) if sel.sum() > 1 else 0; worst = max(worst, np.max((np.abs(r) * np.sqrt(rho))[sel]) if sel.any() else 0)
val = lambda x0: [Ch.chebval(1 - 2 * (x0 - pc["a"]) / (pc["b"] - pc["a"]), pc["coef"]) for pc in pieces if pc["a"] <= x0 <= pc["b"]][0]
print(f"n={n}: lambda = {lam:.13f} (cert {LT}); u(2,5,10) = {val(2.0):.5f} {val(5.0):.5f} {val(10.0):.4f}; ||res||_rho = {np.sqrt(res2):.2e}, max weighted {worst:.2e}")
jet = np.array([Ch.chebval(-1.0, Ch.chebder(pieces[-1]["coef"], d)) * (-2 / (br[-1] - br[-2])) ** d for d in range(4)])
cc = np.linalg.lstsq(np.column_stack([JP, JD]), jet, rcond=None)[0]
print("far-end jet = c1*J_poly + c2*J_dec with c =", cc, " fit residual", np.linalg.norm(np.column_stack([JP, JD]) @ cc - jet))
json.dump({"stip": float(E.STIP), "lam": float(lam), "pieces": pieces, "X": float(br[-1]), "c_poly": float(cc[0]), "c_dec": float(cc[1])}, open(f'psi3md2_{n}.json', 'w'))
