"""Residual of the frozen profile U_1 = c y^p in the far-field slope equation V_tau = -d_y^2(a V'' + b V'^2) - (1/4 - a_m) y V' + b_m cos(theta) V',
a = (1+V^2)^-2, b = -3 V (1+V^2)^-3  (curvature derivative k_sigma = a V'' + b V'^2).   With D = y d/dy = p U d/dU on functions of U:
  a(U) U'' + b(U) U'^2 = y^-2 Phi(U),  Phi = p(p-1) U (1+U^2)^-2 - 3 p^2 U^3 (1+U^2)^-3,
  residual = -y^-4 (D-2)(D-3) Phi(U),   D^j of it = -y^-4 (D-4)^j (D-2)(D-3) Phi(U).
Output: c_j >= sup_{0 < U <= 0.1} |(D-4)^j (D-2)(D-3) Phi(U)| / U, j = 0..3  (arb, 200 cells; the quotient is a polynomial in U^2 over (1+U^2)^n)."""
import json, math, sympy as sp
from flint import arb
U = sp.symbols('U', positive=True); p = sp.Rational(10246, 10000)
Phi = p*(p - 1)*U/(1 + U**2)**2 - 3*p**2*U**3/(1 + U**2)**3
D = lambda f: p*U*sp.diff(f, U)
g = sp.simplify(D(D(Phi)) - 5*D(Phi) + 6*Phi); out = []
for j in range(4):
    q = sp.together(sp.simplify(g/U)); num, den = sp.fraction(q)
    fn, fd = sp.lambdify(U, sp.expand(num), 'math'), sp.lambdify(U, sp.expand(den), 'math')
    pn, pd = sp.Poly(sp.expand(num), U), sp.Poly(sp.expand(den), U)
    ev = lambda P, x: sum((arb(sp.Integer(c.p).__int__())/arb(sp.Integer(c.q).__int__()))*x**int(m[0]) for m, c in P.terms())
    best = 0.0
    for i in range(200):
        x = arb(i).union(arb(i + 1))/2000
        v = abs(ev(pn, x)/ev(pd, x)); best = max(best, math.nextafter(float(v.upper()), math.inf))
    out.append(best); g = sp.simplify(D(g) - 4*g)
json.dump(dict(c=out, Umax=0.1, p='1.0246'), open('cres.json', 'w')); print(out)
