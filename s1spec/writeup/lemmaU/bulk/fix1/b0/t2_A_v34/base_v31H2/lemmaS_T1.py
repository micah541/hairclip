"""T1: symbolic check of the order-j ramp-commutator identity with variable a (Euler-operator test: LHS-RHS must be a total derivative)."""
import sympy as sp
from sympy.calculus.euler import euler_equations
s = sp.Symbol('s'); w = sp.Function('w')(s); a = sp.Function('a')(s); Th = sp.Function('Th')(s)
D = lambda f, n=1: sp.diff(f, s, n)
lhs = Th*w*D(a*D(w, 2), 2)
def rhs(sg):
    return a*Th*D(w, 2)**2 + 2*a*D(Th, 2)*w*D(w, 2) + sg*D(a*D(Th, 2))*w*D(w) - D(a)*D(Th)*D(w)**2
for sg in (-1, +1):
    d = sp.expand(lhs - rhs(sg))
    # total derivative iff Euler-Lagrange operator vanishes identically
    eqs = euler_equations(d, [w], s)
    res = sp.simplify(eqs[0].lhs) if eqs else 0
    print('sign of (a Th\'\')\' w w\' term =', sg, ' total derivative:', res == 0)
