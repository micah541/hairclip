"""Float: alpha-linear E0 form, sum-of-sups vs pointwise spectral bound of the symmetrised (w, w') block."""
import numpy as np, pickle, sympy as sp, json
from numpy.polynomial import chebyshev as Ch
P = pickle.load(open('nlin_E0.pkl', 'rb'))
f1, f2, k, m, p = sp.symbols('phi1 phi2 k m p'); U = sp.symbols('U0 U1 U2'); W = sp.symbols('W0 W1 W2'); H = sp.symbols('H0 H1 H2')
t1 = {(i, j): sp.lambdify((f1, f2, k, m, p) + U, sp.diff(P['T1'], H[i], W[j]), 'numpy') for i in range(3) for j in range(3)}
A = np.loadtxt('/home/claude/s1_tight/prof/prof_tight.txt', comments='END'); TT = 10.2490776510
xg = (TT - A[:, 2])[::-1]; col = lambda j: (0.5 * (A[:, 3 + 2 * j] + A[:, 4 + 2 * j]))[::-1]
x = np.linspace(0, 10.2, 2041)
kk, mm, pp = (np.interp(x, xg, col(j)) for j in (2, 3, 0))
ph1 = (x / 3) * (1 + x**2) ** (-1 / 3); ph2 = (1 + x**2) ** (-1 / 3) / 3 - 2 * x**2 * (1 + x**2) ** (-4 / 3) / 9
J = json.load(open('psi3md2_14.json'))
def ev(x, d):
    out = np.zeros_like(x)
    for pc in J['pieces']:
        a, b = pc['a'], pc['b']; msk = (x >= a) & (x <= b); c = np.array(pc['coef'])
        for _ in range(d): c = Ch.chebder(c)
        out[msk] = Ch.chebval(1 - 2 * (x[msk] - a) / (b - a), c) * (-2 / (b - a)) ** d
    return out
y = [ev(x, d) for d in range(3)]; Z = np.zeros_like(x); da = 1e-7
c = {ij: (f(ph1, ph2, kk, mm, pp, *[da * t for t in y]) - f(ph1, ph2, kk, mm, pp, Z, Z, Z)) / da for ij, f in t1.items()}
S = sum(np.abs(c[(i, j)]) for i in range(2) for j in range(2))
lam = np.zeros_like(x)
for n in range(len(x)):
    M = np.array([[c[(0, 0)][n], (c[(0, 1)][n] + c[(1, 0)][n]) / 2], [(c[(0, 1)][n] + c[(1, 0)][n]) / 2, c[(1, 1)][n]]])
    ev_ = np.linalg.eigvalsh(M); lam[n] = max(abs(ev_))
print("x<=10.2: sup of sum|c_ij| (w,w' block) =", S.max(), " at x =", x[S.argmax()])
print("          sup of spectral radius of sym block =", lam.max(), " at x =", x[lam.argmax()])
print("  c00, c01+c10, c11 at the worst point:", c[(0,0)][lam.argmax()], c[(0,1)][lam.argmax()] + c[(1,0)][lam.argmax()], c[(1,1)][lam.argmax()])
for xv in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0):
    n = int(np.argmin(abs(x - xv))); print(f"  x={xv}: c00={c[(0,0)][n]:+.2f} c01+c10={c[(0,1)][n]+c[(1,0)][n]:+.2f} c11={c[(1,1)][n]:+.2f} c02+c20={c[(0,2)][n]+c[(2,0)][n]:+.2f} c22={c[(2,2)][n]:+.2f}")
