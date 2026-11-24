"""Numerical check of the divergence-form identity: int T1(u=0; h; w) dx = int rho w L_curv h dx  (compact supports),
and of the full nonlinear identity: int T_c(u; w) dx = int rho w F_curv(u) dx for a finite u (F_curv from graph_eq)."""
import numpy as np, pickle, sympy as sp
P = pickle.load(open('nlin_E0.pkl', 'rb'))
syms = sp.symbols('phi1 phi2 k m p U0 U1 U2 W0 W1 W2 H0 H1 H2')
fT1 = sp.lambdify(syms, P['T1'], 'numpy'); fTc = sp.lambdify(syms[:11], P['Tc'], 'numpy')
A = np.loadtxt('/home/claude/s1_tight/prof/prof_tight.txt', comments='END')
TT = 10.2490776510
tc = A[:, 2]; xg = (TT - tc)[::-1]
def col(j): return (0.5 * (A[:, 3 + 2 * j] + A[:, 4 + 2 * j]))[::-1]
pg, qg, kg, mg = col(0), col(1), col(2), col(3)
x = np.linspace(1.0, 9.0, 40001); dx = x[1] - x[0]
k = np.interp(x, xg, kg); m = np.interp(x, xg, mg); p = np.interp(x, xg, pg); q = np.interp(x, xg, qg)
ph1 = (x / 3) * (1 + x ** 2) ** (-1 / 3)
ph2 = np.gradient(ph1, dx)
rho = np.exp(-(1 + x ** 2) ** (2 / 3) / 4)
def bump(c, w):
    t = (x - c) / w; f = np.where(np.abs(t) < 1, np.exp(-1 / np.maximum(1 - t ** 2, 1e-300)), 0.0)
    return f
def ders(f, n=5):
    out = [f]
    for _ in range(n - 1): out.append(np.gradient(out[-1], dx))
    return out
h = ders(bump(4.0, 2.5) + 0.3 * bump(5.5, 1.5)); w = ders(bump(4.5, 2.8))
Z = np.zeros_like(x)
lhs = np.sum(rho * fT1(ph1, ph2, k, m, p, Z, Z, Z, w[0], w[1], w[2], h[0], h[1], h[2])) * dx
Lc = -h[4] - k ** 2 * h[2] - 5 * k * m * h[1] - (3 * m ** 2 + k * p) * h[0]
rhs = np.sum(rho * w[0] * Lc) * dx
print(f"linear identity:  int T1 = {lhs:.10e}   int rho w L_curv h = {rhs:.10e}   rel diff {abs(lhs-rhs)/abs(rhs):.2e}")
# full nonlinear identity with u = 0.05 h
G = pickle.load(open('graph_eq.pkl', 'rb'))
s_, eps = sp.symbols('sigma epsilon'); U = sp.Function('U')(s_)
pf, qf, kf, mf = [sp.Function(n)(s_) for n in ('p', 'q', 'k', 'm')]
Ud = [U] + [sp.Derivative(U, (s_, n)) for n in range(1, 5)]
J = sp.symbols('u0:5'); Pp = sp.symbols('pp qq kk mm')
expr = G['rhs']
for n in range(4, -1, -1): expr = expr.subs(Ud[n], J[n])
expr = expr.subs({pf: Pp[0], qf: Pp[1], kf: Pp[2], mf: Pp[3]})
# curvature part only: remove the drift (alpha(p+u) - u'q)/(4J) * J/alpha = (p + u - u' q/alpha)/4
al = 1 - Pp[2] * eps * J[0]
drift = (Pp[0] + eps * J[0] - eps * J[1] * Pp[1] / al) / 4
fcurv = sp.lambdify((eps,) + J + Pp, expr - drift, 'numpy')
a = 0.05
u = [a * t for t in h]
Fc = fcurv(1.0, *u, p, q, k, m)
rhs2 = np.sum(rho * w[0] * Fc) * dx
lhs2 = np.sum(rho * fTc(ph1, ph2, k, m, p, u[0], u[1], u[2], w[0], w[1], w[2])) * dx
F0 = fcurv(1.0, Z, Z, Z, Z, Z, p, q, k, m)
lhs0 = np.sum(rho * fTc(ph1, ph2, k, m, p, Z, Z, Z, w[0], w[1], w[2])) * dx; rhs0 = np.sum(rho * w[0] * F0) * dx
print(f"nonlinear identity (u = {a} h): int T_c = {lhs2:.10e}  int rho w F_curv = {rhs2:.10e}  rel diff {abs(lhs2-rhs2)/abs(rhs2):.2e}")
print(f"  (u = 0: {lhs0:.6e} vs {rhs0:.6e})")
