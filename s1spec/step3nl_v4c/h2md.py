"""Second-order slow-manifold correction h2:  (L - 2 lamt) h2 = -Q2(psi3^)  on [0, oo), even (h' = h''' = 0 at 0),
admissible at the conical end (jet at X in span{power law s^(A2-4j) series, decaying WKB}, A2 = 1 - 8 lamt).
Q2 = quadratic part of the graph equation (graph_eq.pkl), evaluated on the certified psi3^ (psi3fun) and the float
profile.  Multi-domain Chebyshev collocation on the breakpoints of psi3^ (so Q2(psi3^) is smooth on every element),
C^3 matching, same scheme as psi3md2.py.  Output h2md_<n>.json (pieces, far-end coefficients)."""
import sys, json, pickle, numpy as np, scipy.linalg as sl, sympy as sp
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec'); sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import evans as E
import psi3fun as F
from flint import arb
n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
LT = -0.2561532224577463; L2 = 2 * LT
br = np.array([pc.af for pc in F.PIECES] + [F.PIECES[-1].bf])
# Q2 as a numpy function of the jets
G = pickle.load(open('/home/claude/s1spec/step3nl/graph_eq.pkl', 'rb')); quad = G['quad']
s_ = sp.symbols('sigma'); U = sp.Function('U')(s_)
pf, qf, kf, mf = [sp.Function(nm)(s_) for nm in ('p', 'q', 'k', 'm')]
Ud = [U] + [sp.Derivative(U, (s_, j)) for j in range(1, 5)]
Js = sp.symbols('u0:5'); Pp = sp.symbols('pp qq kk mm')
ex = quad
for j in range(4, -1, -1): ex = ex.subs(Ud[j], Js[j])
ex = ex.subs({pf: Pp[0], qf: Pp[1], kf: Pp[2], mf: Pp[3]})
Q2f = sp.lambdify(Js + Pp, ex, 'numpy')
def prof(x):
    return np.array([E.PROF(s + E.STIP) if s + E.STIP <= 22 else [0.0, s + E.STIP, 0.0, 0.0] for s in x]).T
def psi_jets(x, e):
    P = F.PIECES[e]
    out = []
    for xv in x:
        X = arb(float(xv)); jp = F.jet_poly(P, X, 5); hc = F.Hcorr(P, X, 5)
        out.append([float((a + b).mid()) for a, b in zip(jp, hc)])
    return np.array(out).T
def cheb(N):
    j = np.arange(N + 1); t = np.cos(np.pi * j / N); c = np.ones(N + 1); c[0] = c[-1] = 2; c *= (-1) ** j
    T = np.tile(t, (N + 1, 1)).T; D = np.outer(c, 1 / c) / (T - T.T + np.eye(N + 1)); D -= np.diag(D.sum(1)); return t, D
t, Dt = cheb(n); Ne = len(br) - 1; M = Ne * (n + 1)
A = np.zeros((M, M)); rhs = np.zeros(M); Ds = []; Xn = []
for e in range(Ne):
    a, b = br[e], br[e + 1]; x = a + (1 - t) * (b - a) / 2; D = Dt * (-2 / (b - a))
    D2 = D @ D; D3 = D2 @ D; D4 = D3 @ D; Ds.append((np.eye(n + 1), D, D2, D3)); Xn.append(x)
    p, q, k, m = prof(x)
    Lm = -D4 - np.diag(k * k) @ D2 - np.diag(5 * k * m + q / 4) @ D + np.diag(0.25 - 3 * m * m - k * p - L2)
    o = e * (n + 1); A[o:o + n + 1, o:o + n + 1] = Lm
    J5 = psi_jets(x, e)
    rhs[o:o + n + 1] = -Q2f(*J5, p, q, k, m)
Aex = 1 - 4 * L2
for e in range(Ne):
    o = e * (n + 1); I, D, D2, D3 = Ds[e]
    for r in (o, o + 1, o + n, o + n - 1): A[r] = 0; rhs[r] = 0
    if e == 0:
        A[o, o:o + n + 1] = D[0]; A[o + 1, o:o + n + 1] = D3[0]
    else:
        op = (e - 1) * (n + 1); Ip, Dp, D2p, D3p = Ds[e - 1]
        A[o, o:o + n + 1] = I[0]; A[o, op:op + n + 1] = -Ip[n]
        A[o + 1, o:o + n + 1] = D[0]; A[o + 1, op:op + n + 1] = -Dp[n]
    if e == Ne - 1:
        Xs = br[-1] + E.STIP
        aj = [1.0]
        for jj in range(1, 4): aj.append(aj[-1] * (Aex - 4*jj + 4) * (Aex - 4*jj + 3) * (Aex - 4*jj + 2) * (Aex - 4*jj + 1) / jj)
        Jp = np.zeros(4)
        for jj in range(4):
            exn = Aex - 4 * jj; f = 1.0
            for d in range(4): Jp[d] += aj[jj] * f * Xs ** (exn - d); f *= (exn - d)
        z = -(Xs / 4) ** (1 / 3); z1 = -(1 / 12) * (Xs / 4) ** (-2 / 3); z2 = (1 / 72) * (Xs / 4) ** (-5 / 3)
        Jd = np.array([1.0, z, z1 + z * z, z2 + 3 * z * z1 + z**3])
        Nl = sl.null_space(np.vstack([Jp, Jd]))
        R4 = np.vstack([I[n], D[n], D2[n], D3[n]])
        A[o + n, o:o + n + 1] = Nl[:, 0] @ R4; A[o + n - 1, o:o + n + 1] = Nl[:, 1] @ R4
        JP, JD = Jp, Jd
    else:
        on = (e + 1) * (n + 1); In, Dn, D2n, D3n = Ds[e + 1]
        A[o + n, o:o + n + 1] = D2[n]; A[o + n, on:on + n + 1] = -D2n[0]
        A[o + n - 1, o:o + n + 1] = D3[n]; A[o + n - 1, on:on + n + 1] = -D3n[0]
h = np.linalg.solve(A, rhs)
pieces = []; res2 = 0.0; worst = 0.0; nrm2 = 0.0
for e in range(Ne):
    a, b = br[e], br[e + 1]; he = h[e * (n + 1):(e + 1) * (n + 1)]; cf = Ch.chebfit(t, he, n)
    pieces.append({"a": float(a), "b": float(b), "coef": cf.tolist()})
    xf = np.linspace(a, b, 300); tf = 1 - 2 * (xf - a) / (b - a)
    dd = [Ch.chebval(tf, Ch.chebder(cf, j)) * (-2 / (b - a)) ** j for j in range(5)]
    p, q, k, m = prof(xf)
    J5 = psi_jets(xf, e)
    r = -dd[4] - k * k * dd[2] - (5 * k * m + q / 4) * dd[1] + (0.25 - 3 * m * m - k * p - L2) * dd[0] + Q2f(*J5, p, q, k, m)
    rho = np.exp(-0.25 * (1 + xf**2) ** (2 / 3))
    res2 += np.trapezoid(rho * r**2, xf); worst = max(worst, np.max(np.abs(r) * np.sqrt(rho))); nrm2 += np.trapezoid(rho * dd[0]**2, xf)
val = lambda x0: [Ch.chebval(1 - 2 * (x0 - pc["a"]) / (pc["b"] - pc["a"]), pc["coef"]) for pc in pieces if pc["a"] <= x0 <= pc["b"]][0]
jet = np.array([Ch.chebval(-1.0, Ch.chebder(pieces[-1]["coef"], d)) * (-2 / (br[-1] - br[-2])) ** d for d in range(4)])
cc = np.linalg.lstsq(np.column_stack([JP, JD]), jet, rcond=None)[0]
print(f"n={n}: h2(0) = {val(0.0):.6f}, h2(2,5,10,20) = {val(2.0):.5f} {val(5.0):.5f} {val(10.0):.4f} {val(20.0):.3f}; "
      f"||h2||_rho(half) = {np.sqrt(nrm2):.4f}; ||res||_rho(half) = {np.sqrt(res2):.2e}, max weighted {worst:.2e}")
print("far-end jet = c1*J_poly + c2*J_dec with c =", cc, " fit residual", np.linalg.norm(np.column_stack([JP, JD]) @ cc - jet))
json.dump({"stip": float(E.STIP), "lam2": L2, "A2": Aex, "pieces": pieces, "X": float(br[-1]), "c_poly": float(cc[0]),
           "c_dec": float(cc[1])}, open(f'/home/claude/s1spec/step3nl/h2md_{n}.json', 'w'))
