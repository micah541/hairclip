"""Float exploration of the E0-level nonlinear terms (divergence form):
 alpha-linear quadratic form  int rho sum a_ij w^(i) w^(j),  a_ij = t_ij(u0) - t_ij(0), u0 = alpha psi3^;
 cubic coefficients from T2(u0)[w,w; w] (at t = 0) per monomial w^(i) w^(j) w^(l)."""
import numpy as np, pickle, sympy as sp, sys
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
from flint import arb
import psi3fun as F
P = pickle.load(open('nlin_E0.pkl', 'rb'))
f1, f2, k, m, p = sp.symbols('phi1 phi2 k m p')
U = sp.symbols('U0 U1 U2'); W = sp.symbols('W0 W1 W2'); H = sp.symbols('H0 H1 H2')
args = (f1, f2, k, m, p) + U
T1, T2 = P['T1'], P['T2']
t1 = {(i, j): sp.lambdify(args, sp.diff(T1, H[i], W[j]), 'numpy') for i in range(3) for j in range(3)}
t2 = {}
for i in range(3):
    for j in range(i, 3):
        for l in range(3):
            e = sp.diff(T2, H[i], H[j], W[l])
            if e != 0: t2[(i, j, l)] = (sp.lambdify(args, e, 'numpy'), 1 if i == j else 2)
print("nonzero cubic monomials:", sorted(t2))
# grid and profile
A = np.loadtxt('/home/claude/s1_tight/prof/prof_tight.txt', comments='END')
TT = 10.2490776510; xg = (TT - A[:, 2])[::-1]
col = lambda j: (0.5 * (A[:, 3 + 2 * j] + A[:, 4 + 2 * j]))[::-1]
x = np.concatenate([np.linspace(0.0, 10.2, 1021), np.linspace(10.25, 31.25, 841)])
kk = np.interp(x, xg, col(2), right=0.0); mm = np.interp(x, xg, col(3), right=0.0); pp = np.interp(x, xg, col(0), right=0.0)
ph1 = (x / 3) * (1 + x ** 2) ** (-1 / 3); ph2 = (1 + x ** 2) ** (-1 / 3) / 3 - 2 * x ** 2 * (1 + x ** 2) ** (-4 / 3) / 9
rho = np.exp(-(1 + x ** 2) ** (2 / 3) / 4)
def psi_jets(xv):
    X = arb(float(xv))
    if xv >= 28.0: return [float(v.mid()) for v in F.ext_jet(X, 3)]
    for Pc in F.PIECES:
        if Pc.af <= xv <= Pc.bf:
            return [float((a + b).mid()) for a, b in zip(F.jet_poly(Pc, X, 3), F.Hcorr(Pc, X, 3))]
J3 = np.array([psi_jets(v) for v in x])
Z = np.zeros_like(x)
for alpha in (1e-3, 3e-3):
    u0 = [alpha * J3[:, d] for d in range(3)]
    a = {ij: t1[ij](ph1, ph2, kk, mm, pp, *u0) - t1[ij](ph1, ph2, kk, mm, pp, Z, Z, Z) for ij in t1}
    print(f"alpha = {alpha}: sup_x |a_ij|/alpha (alpha-linear form, weight rho already factored out)")
    for reg, msk in (("x<=12", x <= 12), ("12<x<=31.25", x > 12)):
        print("   ", reg, {f"{i}{j}": float(f"{np.max(np.abs(a[(i,j)][msk]))/alpha:.3g}") for (i, j) in sorted(a)})
    print("   cubic coefficient sup (t = 0), per monomial h_i h_j w_l:")
    for reg, msk in (("x<=12", x <= 12), ("12<x<=31.25", x > 12)):
        print("   ", reg, {f"{i}{j}{l}": float(f"{np.max(np.abs(fn(ph1[msk], ph2[msk], kk[msk], mm[msk], pp[msk], *[v[msk] for v in u0]) * mult)):.3g}") for (i, j, l), (fn, mult) in sorted(t2.items())})
