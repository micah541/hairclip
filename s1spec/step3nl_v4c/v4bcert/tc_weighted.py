"""Symbolic part of the v4b E0 strip / drift certificate (sympy, exact).
(1) Curvature: the divergence form of nlin_E0.py with a GENERAL weight Psi = rho * X (X = chi^2 on the strip, X = 1 on the core):
      int Psi w F_curv(u) = int T_c^X(u; w),  T_c^X = -[Phi J^-1 w'' + B1 w' + B0 w] kappa,  Phi = Psi/alpha_g,
      B1 = (Phi J^-1)_x + Phi_x J^-1, B0 = (Phi_x J^-1)_x   (two integrations by parts; boundary terms vanish: parity at x = 0, X = 0 beyond the strip).
    T_c^X / rho is linear in W = (W0, W1, W2) with coefficients C_j(U; X, X1, X2); the alpha-linear form is
      sum_ij [d C_j/d U_i](u0) w^(i) w^(j) - (same at u = 0),   i, j in {0, 1, 2}.
    Check: with X = 1, X1 = X2 = 0 the coefficients equal those of the closed v3.1 nlin_E0.pkl (T1 -> d^2/dH_i dW_j), symbolically.
(2) Euler-operator check of (1): E_w[ rho X w F_curv(u) - T_c^X ] == 0 for generic u, profile and X (so the difference is a total derivative).
(3) Drift: D(u) = -(q/4) u' G(u), G(u) = k u/(1 - k u) (nonlinear part of (1/4)(p + u - u' q/alpha_g)); its alpha-linear part
      Dlin(w) = -(q/4)[G(u0) w' + u0' G_u(u0) w];  int rho X w Dlin(w) = int rho w^2 c_dr  (after one integration by parts), with
      c_dr = (1/8) X [q' G - phi' q G + q (m u0 - k u0')/(1 - k u0)^2] + (1/8) X1 q G.      Checked: Euler operator of
      rho X w Dlin(w) - rho c_dr w^2 vanishes (difference is a total derivative).
Output: out/tc_weighted.pkl (sympy expressions), out/tc_weighted.json (check results, args)."""
import sys, os, time, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4bcommon import OUT, args_block, dump
import sympy as sp
t0 = time.time()
f1, f2, f3 = sp.symbols('phi1 phi2 phi3')
k, m, p, q = sp.symbols('k m p q')
U0, U1, U2, U3, U4 = sp.symbols('U0 U1 U2 U3 U4')
W0, W1, W2, W3, W4 = sp.symbols('W0 W1 W2 W3 W4')
X, X1, X2, X3 = sp.symbols('X X1 X2 X3')
R = sp.symbols('R')
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4, U0: U1, U1: U2, U2: U3, U3: U4, W0: W1, W1: W2, W2: W3, W3: W4,
         R: -f1 * R, f1: f2, f2: f3, X: X1, X1: X2, X2: X3}
def D(e): return sum(sp.diff(e, v) * r for v, r in rules.items())
al = 1 - k * U0
Jv = sp.sqrt(al ** 2 + U1 ** 2)
kap = (al * U2 + k * Jv ** 2 + U1 * (m * U0 + k * U1)) / Jv ** 3
Phi = R * X / al
Phx = D(Phi)
B1 = D(Phi / Jv) + Phx / Jv
B0 = D(Phx / Jv)
Tc = -(Phi / Jv * W2 + B1 * W1 + B0 * W0) * kap
TcR = Tc.subs(R, 1)
assert not TcR.has(U3) and not TcR.has(W3) and not TcR.has(f3) and not TcR.has(X3)
C = [sp.diff(TcR, Wj) for Wj in (W0, W1, W2)]
assert sp.simplify(TcR - sum(Cj * Wj for Cj, Wj in zip(C, (W0, W1, W2)))) == 0, 'T_c not linear in W'
def exact_zero_residual(expr):
    """exact check that a rational expression in the jets and sqrt(alpha_g^2 + U1^2) vanishes identically: every power of the base
    alpha_g^2 + U1^2 with half-integer exponent n/2 is replaced by S^n, denominators are cleared, and S^2 is reduced to the base."""
    Ss = sp.symbols('S', positive=True)
    base = sp.expand((1 - k * U0) ** 2 + U1 ** 2)
    rep = {}
    for a in expr.atoms(sp.Pow):
        if a.exp.is_Rational and a.exp.q == 2 and sp.expand(a.base - base) == 0: rep[a] = Ss ** int(2 * a.exp)
    e2 = expr.xreplace(rep)
    assert not any(isinstance(a, sp.Pow) and a.exp.is_Rational and a.exp.q != 1 for a in e2.atoms(sp.Pow)), 'unreplaced radical'
    nr, _d = sp.fraction(sp.together(e2))
    Pq = sp.Poly(sp.expand(nr), Ss); rr = 0
    for (ee,), cf in Pq.terms(): rr += cf * base ** (ee // 2) * Ss ** (ee % 2)
    return sp.expand(rr)
TF = {(i, j): sp.diff(C[j], (U0, U1, U2)[i]) for i in range(3) for j in range(3)}
res = {}
# (1) X = 1 reduction against the closed nlin_E0.pkl
P = pickle.load(open('/home/claude/s1spec/step3nl/nlin_E0.pkl', 'rb'))
H = sp.symbols('H0 H1 H2')
ok1 = True
for i in range(3):
    for j in range(3):
        old = sp.diff(P['T1'], H[i], (W0, W1, W2)[j])
        new = TF[(i, j)].subs({X: 1, X1: 0, X2: 0})
        rr = exact_zero_residual(old - new)
        ok1 = ok1 and (rr == 0)
        res[f'X1_reduction_{i}{j}_exact_residual'] = str(rr)
res['X1_reduction_ok'] = bool(ok1)
print(f'[{time.time()-t0:.0f}s] (1) X=1 reduction vs nlin_E0.pkl: {ok1}', flush=True)
# (2) Euler-operator check (numerical substitution of the jets of generic functions at a point is not enough for a total-derivative check,
#     so check symbolically after substituting polynomial test functions in x for all coefficient functions -> exact rational arithmetic)
xs = sp.symbols('x')
def euler_check(L_expr, wsyms, order):
    """Euler-Lagrange expression sum_j (-d)^j dL/dw^(j) with total derivative D; returns the expression"""
    E = 0
    for j, wj in enumerate(wsyms[:order + 1]):
        t = sp.diff(L_expr, wj)
        for _ in range(j): t = D(t)
        E += (-1) ** j * t
    return E
# F_curv(u) = -alpha_g^-1 d/dx(J^-1 kappa_x)
kx = D(kap); Fcurv = -(1 / al) * D(kx / Jv)
Lag = R * X * W0 * Fcurv - Tc
Eu = euler_check(Lag, (W0, W1, W2, W3, W4), 4)
# evaluate at random rational point (all jets independent symbols -> the identity must hold identically; a random-point test with
# 3 points of exact rationals is a probabilistic check, plus the polynomial-identity degree argument is not claimed)
import random
random.seed(7)
syms = sorted(Eu.free_symbols, key=str)
mx2 = 0.0
for trial in range(3):
    num = {s: sp.Rational(random.randint(1, 30), random.randint(31, 60)) for s in syms}
    num[R] = sp.Rational(1); v = sp.nsimplify(Eu.subs(num))
    mx2 = max(mx2, float(abs(sp.N(v, 60))))
res['euler_curv_max_abs_at_3_random_rational_points'] = mx2
# exact: replace sqrt(alpha_g^2 + U1^2) by a symbol S, clear denominators, reduce S^2 -> alpha_g^2 + U1^2
red = exact_zero_residual(Eu)
res['euler_curv_exact_residual'] = str(red)
res['euler_curv_ok'] = bool(mx2 == 0.0 and red == 0)
print(f'[{time.time()-t0:.0f}s] (2) Euler check curvature: {mx2}', flush=True)
# (3) drift
G = lambda u: k * u / (1 - k * u)
Gu = lambda u: k / (1 - k * u) ** 2
Dlin = -(q / 4) * (G(U0) * W1 + U1 * Gu(U0) * W0)
c_dr = sp.Rational(1, 8) * X * (D(q) * G(U0) - f1 * q * G(U0) + q * (m * U0 - k * U1) / (1 - k * U0) ** 2) + sp.Rational(1, 8) * X1 * q * G(U0)
Lag_d = R * X * W0 * Dlin - R * c_dr * W0 ** 2
Ed = euler_check(Lag_d, (W0, W1, W2), 2)
Ed_s = sp.simplify(Ed)
res['euler_drift_simplified'] = str(Ed_s)
res['euler_drift_ok'] = bool(Ed_s == 0)
# the nonlinear drift really is D(u) = drift(u) - drift_lin(u): (1/4)(p + u - u' q/alpha_g) - (1/4)(p + u - u' q) = -(q/4) u' G(u)
chk = sp.simplify((sp.Rational(1, 4) * (p + U0 - U1 * q / (1 - k * U0)) - sp.Rational(1, 4) * (p + U0 - U1 * q)) - (-(q / 4) * U1 * G(U0)))
res['drift_split_ok'] = bool(chk == 0)
# alpha-linear part of D(u0 + w) - D(u0) in w
e = sp.symbols('e')
lin = sp.diff((-(q / 4) * (U1 + e * W1) * G(U0 + e * W0)), e).subs(e, 0)
res['drift_lin_ok'] = bool(sp.simplify(lin - Dlin) == 0)
print(f'[{time.time()-t0:.0f}s] (3) drift: euler {Ed_s}, split {res["drift_split_ok"]}, lin {res["drift_lin_ok"]}', flush=True)
assert res['X1_reduction_ok'] and res['euler_curv_ok'] and res['euler_drift_ok'] and res['drift_split_ok'] and res['drift_lin_ok']
pickle.dump(dict(TF=TF, C=C, c_dr=c_dr, Dlin=Dlin), open(os.path.join(OUT, 'tc_weighted.pkl'), 'wb'))
res['runtime_s'] = round(time.time() - t0, 1)
res['args'] = args_block(['/home/claude/s1spec/step3nl/nlin_E0.pkl'])
print('sha256', dump(os.path.join(OUT, 'tc_weighted.json'), res))
