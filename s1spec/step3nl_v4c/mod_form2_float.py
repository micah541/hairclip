"""Float: exact spectral radii (on {f1,f2}^perp, half-line rho) of the two modulation forms
   alpha-level:  Q1(w) = -<w,ha><gA,w> + <w,hb><gB,w>       (current bound 268 = Minv*Lfn*||(gA,gB)||)
   E2-level:     Q2(w) = -<w,ha><D f1,w> - <w,hb><D f2,w>     (current bound 458.4 = Minv*Lfn*||(Df1,Df2)||)
D g = g'''' - 2 phi' g''' + (phi'^2 - phi'') g''.  Uses ODE-exact derivatives of the profile (no FD of p, sin th)."""
import numpy as np, json, sys
from numpy.polynomial import chebyshev as Ch
from scipy.linalg import eigh
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
h = 0.002; x = np.arange(0, 31 + h / 2, h); sig = x + E.STIP
P = np.array([E.PROF(s) if s <= 22 else [0, s, 0, 0] for s in sig]).T; p, q, k, m = P
th = -np.pi / 2 + np.concatenate([[0], np.cumsum((k[1:] + k[:-1]) / 2) * h])
phi = (1 + x**2) ** (2 / 3) / 4; rho = np.exp(-phi)
w8 = np.full(len(x), h); w8[0] = w8[-1] = h / 2; w8 *= rho
ip = lambda f, g: np.sum(w8 * f * g)
def D(f, n=1):
    for _ in range(n): f = np.gradient(f, h, edge_order=2)
    return f
uu = np.clip((32 - sig) / 3, 0, 1); chi = uu**4 * (35 - 84 * uu + 70 * uu**2 - 20 * uu**3)
J = json.load(open('psi3md2_14.json'))
def ev(xv, d):
    out = np.zeros_like(xv)
    for pc in J['pieces']:
        a, b = pc['a'], pc['b']; msk = (xv >= a) & (xv <= b); c = np.array(pc['coef'])
        for _ in range(d): c = Ch.chebder(c)
        out[msk] = Ch.chebval(1 - 2 * (xv[msk] - a) / (b - a), c) * (-2 / (b - a)) ** d
    return out
y0, y1 = ev(x, 0), ev(x, 1)
def Lstar(g):
    rg = rho * g; c0 = 0.25 - 3 * m**2 - k * p
    return (-D(rg, 4) - D(k**2 * rg, 2) + D((5 * k * m + q / 4) * rg) + c0 * rg) / rho
p1, p2 = D(phi), D(phi, 2)
def Drho(g): return D(g, 4) - 2 * p1 * D(g, 3) + (p1**2 - p2) * D(g, 2)
f1 = chi * p; f2 = -chi * np.sin(th)
M = np.array([[ip(chi * [p, -np.sin(th)][i], [f1, f2][j]) for j in range(2)] for i in range(2)]); Minv = np.linalg.inv(M)
mask = ((x > 0.3) & (x < 30.5)).astype(float)
Lf = [mask * Lstar(f1), mask * Lstar(f2)]
ha = Minv[0, 0] * Lf[0] + Minv[1, 0] * Lf[1]; hb = Minv[0, 1] * Lf[0] + Minv[1, 1] * Lf[1]
gA = chi * (y0 - q * y1); gB = chi * y1 * np.cos(th)
Df1, Df2 = mask * Drho(f1), mask * Drho(f2)
print("||Lf1||, ||Lf2||:", np.sqrt(ip(Lf[0], Lf[0])), np.sqrt(ip(Lf[1], Lf[1])), "  ||D f1||, ||D f2||:", np.sqrt(ip(Df1, Df1)), np.sqrt(ip(Df2, Df2)))
def proj(v, basis):
    G = np.array([[ip(bi, bj) for bj in basis] for bi in basis]); c = np.linalg.solve(G, [ip(v, bi) for bi in basis])
    return v - sum(ci * bi for ci, bi in zip(c, basis))
def rad(pairs, basis):
    H = []; idx = []
    for (hh, gg, s) in pairs: H += [hh, gg]
    if basis: H = [proj(v, basis) for v in H]
    n = len(H); G = np.array([[ip(a_, b_) for b_ in H] for a_ in H]); C = np.zeros((n, n))
    for i, (_, _, s) in enumerate(pairs): C[2 * i, 2 * i + 1] = C[2 * i + 1, 2 * i] = s / 2
    lam = eigh(G @ C @ G, G + 1e-14 * np.eye(n), eigvals_only=True); return lam.min(), lam.max()
for lab, basis in (("no constraint", None), ("perp f1,f2", [f1, f2])):
    print(f"{lab:14s}: alpha-form eig range {rad([(ha, gA, -1), (hb, gB, 1)], basis)}   E2-form eig range {rad([(ha, Df1, -1), (hb, Df2, -1)], basis)}")
