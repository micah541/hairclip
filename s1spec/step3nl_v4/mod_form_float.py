"""Float: exact size of the modulation alpha-term  Q(w)/alpha = -a(w) <g_A, w> + b(w) <g_B, w>  on {f1,f2,f3}^perp,
(a, b) = (<w, L*f1>, <w, L*f2>) M^-1, M_ij = <chi psi_i, f_j>, g_A = chi(psi3 - q psi3'), g_B = chi psi3' cos(theta).
Compares with the Cauchy-Schwarz bound Minv*Lfn*||(g_A, g_B)|| = 268.  f3 approximated by the float adjoint mode."""
import numpy as np, json, sys
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
h = 0.002; x = np.arange(0, 30 + h / 2, h); sig = x + E.STIP
P = np.array([E.PROF(s) if s <= 22 else [0, s, 0, 0] for s in sig]).T; p, q, k, m = P
th = -np.pi / 2 + np.concatenate([[0], np.cumsum((k[1:] + k[:-1]) / 2) * h])
rho = np.exp(-(1 + x**2) ** (2 / 3) / 4); w8 = np.full(len(x), h); w8[0] = w8[-1] = h / 2; w8 *= rho
ip = lambda f, g: np.sum(w8 * f * g)
u = np.clip((32 - sig) / 3, 0, 1); chi = u**4 * (35 - 84 * u + 70 * u**2 - 20 * u**3)
J = json.load(open('psi3md2_14.json'))
def ev(xv, d):
    out = np.zeros_like(xv)
    for pc in J['pieces']:
        a, b = pc['a'], pc['b']; msk = (xv >= a) & (xv <= b); c = np.array(pc['coef'])
        for _ in range(d): c = Ch.chebder(c)
        out[msk] = Ch.chebval(1 - 2 * (xv[msk] - a) / (b - a), c) * (-2 / (b - a)) ** d
    return out
y0, y1 = ev(x, 0), ev(x, 1)
def D(f, n=1):
    for _ in range(n): f = np.gradient(f, h, edge_order=2)
    return f
def Lstar(g):                       # rho^-1 [ -(rho g)'''' - (k^2 rho g)'' + ((5km + q/4) rho g)' + c0 rho g ]
    rg = rho * g; c0 = 0.25 - 3 * m**2 - k * p
    return (-D(rg, 4) - D(k**2 * rg, 2) + D((5 * k * m + q / 4) * rg) + c0 * rg) / rho
f1 = chi * p; f2 = -chi * np.sin(th)
# f3: float adjoint eigenfunction ~ rho^-1 * (L^T)-eigenvector; crude proxy: solve (L* - lamt) F = 0 via the
# certified relation <psi3, f3> ~ 0.032 is not available here, so use f3 := chi * (psi3 projected) as a stand-in
# (the constraint f3 only matters through the projection; we report both with and without it).
psi = [p, -np.sin(th)]
M = np.array([[ip(chi * psi[i], [f1, f2][j]) for j in range(2)] for i in range(2)])
Minv = np.linalg.inv(M)
Lf = [Lstar(f1), Lstar(f2)]
sl = slice(200, len(x) - 200)        # avoid finite-difference edge effects (x in [0.4, 29.6])
mask = np.zeros_like(x); mask[sl] = 1
gA = chi * (y0 - q * y1); gB = chi * y1 * np.cos(th)
# Q(w)/alpha = -a <gA, w> + b <gB, w>,  (a, b) = (<w,Lf1>, <w,Lf2>) Minv  =>  a = <w, ha>, b = <w, hb>
ha = mask * (Minv[0, 0] * Lf[0] + Minv[1, 0] * Lf[1]); hb = mask * (Minv[0, 1] * Lf[0] + Minv[1, 1] * Lf[1])
def proj(v, basis):
    G = np.array([[ip(bi, bj) for bj in basis] for bi in basis]); c = np.linalg.solve(G, [ip(v, bi) for bi in basis])
    return v - sum(ci * bi for ci, bi in zip(c, basis))
for label, basis in (("no constraint", None), ("perp f1,f2", [f1, f2]), ("perp f1,f2,psi3-proxy", [f1, f2, chi * y0])):
    H = [ha, hb, gA, gB] if basis is None else [proj(v, basis) for v in (ha, hb, gA, gB)]
    G = np.array([[ip(a_, b_) for b_ in H] for a_ in H])
    # symmetric form S = 1/2 sum_i (h_i g_i^T + g_i h_i^T) with pairs (ha, -gA), (hb, +gB); in the span of H:
    C = np.zeros((4, 4)); C[0, 2] = C[2, 0] = -0.5; C[1, 3] = C[3, 1] = 0.5
    # generalized eigenproblem: max |w^T S w| / |w|^2 with w = sum c_j H_j  ->  (G C G) c = lam G c
    from scipy.linalg import eigh
    lam = eigh(G @ C @ G, G, eigvals_only=True)
    print(f"{label:24s}: spectral radius of the modulation form / alpha = {np.max(np.abs(lam)):.2f}   "
          f"(||ha|| {np.sqrt(ip(H[0],H[0])):.2f}, ||hb|| {np.sqrt(ip(H[1],H[1])):.2f}, ||gA|| {np.sqrt(ip(H[2],H[2])):.2f}, ||gB|| {np.sqrt(ip(H[3],H[3])):.2f})")
print("Cauchy-Schwarz bound used so far: 268;  ||(Lf1,Lf2)|| =", np.sqrt(ip(Lf[0]*mask, Lf[0]*mask) + ip(Lf[1]*mask, Lf[1]*mask)), " ||Minv|| =", np.linalg.norm(Minv, 2))
