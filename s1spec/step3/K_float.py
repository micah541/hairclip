"""How far is f3 from being an eigenfunction of L (not L*)?  K = ||(L - lamt) F|| / ||F|| on sigma <= 29,
and the same after removing the psi1, psi2 components (what enters g with phi3 = projected f3)."""
import numpy as np, sympy as sp
from math import comb
from scipy.integrate import solve_ivp
import gram_float as GF
from gram_float import coef, SIG0
import gen_hs3 as G3
STIP = GF.STIP; LAMT = G3.LAMT
def f(tau, Y):
    c = coef(SIG0 - tau, -GF.OMEGA); k, m, q, P1 = c['k'], c['m'], c['q'], c['P1']
    M = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [(0.5 - LAMT), (k * m + q / 4), -k * k, 0]])
    return -(M @ Y) - P1 * Y
sol = solve_ivp(f, [0, SIG0 - STIP], np.array(G3.wstar_init()), method='DOP853', rtol=1e-12, atol=1e-18, dense_output=True)
Xg = np.linspace(0, 29 - STIP, 30001); sig = Xg + STIP
W = np.array([sol.sol(SIG0 - s) for s in sig]); prof = np.array([GF.prof(s) for s in sig])
p, q, k, m = prof.T
w4 = -k**2 * W[:, 2] + (k * m + q / 4) * W[:, 1] + (0.5 - LAMT) * W[:, 0]
Ws = [W[:, 0], W[:, 1], W[:, 2], W[:, 3], w4]
Xs = sp.symbols('Xs'); ph = (1 + Xs**2) ** sp.Rational(2, 3) / 4
B = [sp.lambdify(Xs, sp.simplify(sp.diff(sp.exp(ph), Xs, i) / sp.exp(ph)), 'numpy') for i in range(5)]
Bv = [B[i](Xg) * np.ones_like(Xg) for i in range(5)]
Fd = [sum(comb(n, j) * Ws[j] * Bv[n - j] for j in range(n + 1)) for n in range(5)]
rho = np.exp(-(1 + Xg**2) ** (2 / 3) / 4)
LF = -Fd[4] - k**2 * Fd[2] - (5 * k * m + q / 4) * Fd[1] + (0.25 - 3 * m**2 - k * p) * Fd[0]
R = LF - LAMT * Fd[0]
ip = lambda a, b: np.trapezoid(rho * a * b, Xg)
nF = np.sqrt(ip(Fd[0], Fd[0]))
print(f"||F|| = {nF:.4e}   K = ||(L-lamt)F||/||F|| = {np.sqrt(ip(R, R))/nF:.3f}")
# remove components along psi1 = p, psi2 = <tau, tau_tip> (orthogonal projection, numerical)
th = np.array([0.0])
# tangent angle relative to the tip: theta - theta_tip = int_tip^sigma k
from scipy.integrate import cumulative_trapezoid
dth = cumulative_trapezoid(k, sig, initial=0.0)
psi2 = np.cos(dth)
Psi = np.vstack([p, psi2])
Gm = np.array([[ip(a, b) for b in Psi] for a in Psi]); bv = np.array([ip(a, R) for a in Psi])
Rp = R - Psi.T @ np.linalg.solve(Gm, bv)
print(f"after removing psi1, psi2:  ||.||/||F|| = {np.sqrt(ip(Rp, Rp))/nF:.3f}")
# how much of F lies along psi1, psi2 (cos of angle)
Fp = Fd[0] - Psi.T @ np.linalg.solve(Gm, np.array([ip(a, Fd[0]) for a in Psi]))
print(f"||(I-P)F||/||F|| = {np.sqrt(ip(Fp, Fp))/nF:.4f}")
