"""Symbolic checks for the symmetric part H_sym of L in L^2(rho), rho = exp(-phi).

(1) rho u Lu - rho(-u''^2 + A u'^2 + V u^2) is an exact derivative (variational derivative vanishes),
    using the profile relations k' = m, m' = p/4, p' = -k q, q' = 1 + k p.
(2) The scaled Jacobi system for (H_sym - lam) u = f-type problems:
        y = (u, u', u'', y4),  y4 = -u''' + phi' u'' - A u'
        y1' = y2, y2' = y3, y3' = -A y2 + phi' y3 - y4, y4' = -(V - lam) y1 + phi' y4  [ - g for forcing ]
    reproduces (rho u'')'' + (rho A u')' - rho (V - lam) u = rho g.
"""
import sympy as sp

x = sp.symbols('x')
u = sp.Function('u')(x)
phi = sp.Function('phi')(x)
k, m, p, q = [sp.Function(n)(x) for n in 'kmpq']
prof = {sp.Derivative(k, x): m, sp.Derivative(m, x): p/4, sp.Derivative(p, x): -k*q, sp.Derivative(q, x): 1 + k*p}

def D(e, n=1):
    for _ in range(n):
        e = sp.diff(e, x).subs(prof)
    return e

rho = sp.exp(-phi)
b = 5*k*m + q/4
c = sp.Rational(1, 4) - 3*m**2 - k*p
Lu = -sp.diff(u, x, 4) - k**2*sp.diff(u, x, 2) - b*sp.diff(u, x) + c*u
A = 2*(sp.diff(phi, x)**2 - sp.diff(phi, x, 2)) + k**2
r4 = sp.simplify(sp.diff(rho, x, 4)/rho)
k2 = k**2
V = (-sp.Rational(1, 2)*r4 - sp.Rational(1, 2)*(D(k2, 2) - 2*D(k2)*sp.diff(phi, x) + k2*(sp.diff(phi, x)**2 - sp.diff(phi, x, 2)))
     + sp.Rational(1, 2)*(D(b) - b*sp.diff(phi, x)) + c)
E = rho*u*Lu - rho*(-sp.diff(u, x, 2)**2 + A*sp.diff(u, x)**2 + V*u**2)
# variational derivative w.r.t. u
EL = 0
for j in range(0, 5):
    dE = sp.diff(E, sp.diff(u, x, j)) if j > 0 else sp.diff(E, u)
    EL += (-1)**j * D(dE, j) if j > 0 else dE
EL = sp.simplify(sp.expand(EL.subs(prof)))
print("(1) variational derivative of the difference:", EL)

# (2) Jacobi system check
lam, g = sp.symbols('lam'), sp.Function('g')(x)
ph1 = sp.diff(phi, x)
y1, y2, y3 = u, sp.diff(u, x), sp.diff(u, x, 2)
y4 = -sp.diff(u, x, 3) + ph1*y3 - A*y2
# claim: y4' = -(V-lam) y1 + ph1 y4 - g  <=>  (rho u'')'' + (rho A u')' - rho(V-lam)u = rho g
lhs = D(y4) - (-(V - lam)*y1 + ph1*y4)
target = -(D(D(rho*sp.diff(u, x, 2))) + D(rho*A*sp.diff(u, x)) - rho*(V - lam)*u)/rho
print("(2) Jacobi: y4' - (-(V-lam) y1 + phi' y4) + [EL operator]/rho =",
      sp.simplify(sp.expand((lhs - target).subs(prof))))
# also: EL of the form int rho(-u''^2 + A u'^2 + V u^2) is -2[(rho u'')'' + (rho A u')' - rho V u]
F = rho*(-sp.diff(u, x, 2)**2 + A*sp.diff(u, x)**2 + V*u**2)
EL2 = 0
for j in range(0, 3):
    dE = sp.diff(F, sp.diff(u, x, j)) if j > 0 else sp.diff(F, u)
    EL2 += (-1)**j * D(dE, j) if j > 0 else dE
chk = sp.simplify(sp.expand((EL2 + 2*(D(D(rho*sp.diff(u, x, 2))) + D(rho*A*sp.diff(u, x)) - rho*V*u)).subs(prof)))
print("(3) EL of the local form = -2[(rho u'')'' + (rho A u')' - rho V u]:", chk)
# explicit V in terms of phi-derivatives, for code generation
P1, P2, P3, P4 = sp.symbols('P1 P2 P3 P4')
Vexp = sp.expand(V.subs({sp.diff(phi, x, 4): P4, sp.diff(phi, x, 3): P3, sp.diff(phi, x, 2): P2, sp.diff(phi, x): P1}))
print("V =", Vexp)
