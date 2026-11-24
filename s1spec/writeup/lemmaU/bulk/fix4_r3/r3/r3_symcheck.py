#!/usr/bin/env python3
"""fix4_r3: symbolic checks of the R3 derivation (r3_forcing.py docstring), exact (sympy), residuals must be 0.
 (1) y-form vs s-form of the eps-part on U = c alpha y^p:  -d_y^2 (A(U) U'' + B(U) U'^2) == -y^-4 (d_s^2 - 5 d_s + 6) G(U),  G = A(U)(U_ss - U_s) + B(U) U_s^2.
 (2) d_s^j [ -y^-4 (d^2 - 5d + 6) G(U) ] == -y^-4 Phi_j(U) with Phi_j of r3_forcing (j = 0..3), U = c alpha e^{p s}.
 (3) the non-eps part: -U_tau - c U_s with alpha_tau = lam alpha + beta, c = 1/4 - a - b cos(th) e^-s equals
     -c_U y^p (beta + (lam + p/4) alpha - p a alpha) + c_U p b cos(th) alpha y^(p-1).
 (4) the profile-free check of the slope equation is lemmaS/lemmaS_symbolic.py (not rerun here; A, B as there).
Writes r3_symcheck.json (args block)."""
import os, sys, json, hashlib, platform
sys.dont_write_bytecode = True
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import r3_forcing as R3
y, s, al, cU, p, lam, beta, a, b, th = sp.symbols('y s alpha c_U p lam beta a b theta', positive=True)
Af = lambda x: (1 + x**2)**-2; Bf = lambda x: -3*x*(1 + x**2)**-3
res = {}
# (1) with symbolic p
U = cU*al*y**p
lhs = -sp.diff(Af(U)*sp.diff(U, y, 2) + Bf(U)*sp.diff(U, y)**2, y, 2)
Us = cU*al*sp.exp(p*s)
G = Af(Us)*(sp.diff(Us, s, 2) - sp.diff(Us, s)) + Bf(Us)*sp.diff(Us, s)**2
rhs = -sp.exp(-4*s)*(sp.diff(G, s, 2) - 5*sp.diff(G, s) + 6*G)
r1 = sp.simplify((lhs - rhs.subs(s, sp.log(y))).subs(y, sp.exp(s)))
res['1_yform_vs_sform'] = str(r1)
# (2) with p = 1.0246 exact
pq = R3.P_Q; Uq = cU*al*sp.exp(pq*s)
Gq = Af(Uq)*(sp.diff(Uq, s, 2) - sp.diff(Uq, s)) + Bf(Uq)*sp.diff(Uq, s)**2
Reps = -sp.exp(-4*s)*(sp.diff(Gq, s, 2) - 5*sp.diff(Gq, s) + 6*Gq)
for j in range(4):
    phi = R3._phi(j).subs(R3._U, Uq)
    rj = sp.simplify(sp.diff(Reps, s, j) + sp.exp(-4*s)*phi)
    res['2_level_%d' % j] = str(rj)
# (3)
alpha = sp.Function('alpha')(sp.Symbol('tau'))
Ut = cU*alpha*sp.exp(p*s); c = sp.Rational(1, 4) - a - b*sp.cos(th)*sp.exp(-s)
expr = -sp.diff(Ut, sp.Symbol('tau')) - c*sp.diff(Ut, s)
expr = expr.subs(sp.Derivative(alpha, sp.Symbol('tau')), lam*alpha + beta)
target = -cU*sp.exp(p*s)*(beta + (lam + p/4)*alpha - p*a*alpha) + cU*p*b*sp.cos(th)*alpha*sp.exp((p - 1)*s)
res['3_noneps'] = str(sp.simplify(sp.expand(expr - target)))
ok = all(v == '0' for v in res.values()); res['all_zero'] = ok
ins = [os.path.abspath(__file__), os.path.join(HERE, 'r3_forcing.py')]
res['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version.split()[0], platform=platform.platform(), sympy=sp.__version__,
                   env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, inputs={q: hashlib.sha256(open(q, 'rb').read()).hexdigest() for q in ins})
json.dump(res, open(os.path.join(HERE, 'r3_symcheck.json'), 'w'), indent=1, sort_keys=True)
print({k: v for k, v in res.items() if k != 'args'}); assert ok
