"""G0_INDEP: symbolic check (Euler operator = 0 <=> total derivative) of the level-0 top-order identity used in G0_INDEP.md, item Q0-1:
   W v (c v''')'  ==  W c v''^2 + W (c' - 2 psi' c) v' v'' - (W psi' c)' v v''   (mod total derivatives),  psi' = -W'/W."""
import sympy as sp
x = sp.symbols('x'); W, c, v = [sp.Function(n)(x) for n in ('W', 'c', 'v')]
psi1 = -sp.diff(W, x) / W
lhs = W * v * sp.diff(c * sp.diff(v, x, 3), x)
rhs = W * c * sp.diff(v, x, 2)**2 + W * (sp.diff(c, x) - 2 * psi1 * c) * sp.diff(v, x) * sp.diff(v, x, 2) - sp.diff(W * psi1 * c, x) * v * sp.diff(v, x, 2)
Lg = sp.expand(lhs - rhs)
E = sum((-1)**j * sp.diff(sp.diff(Lg, sp.diff(v, x, j) if j else v), x, j) for j in range(0, 5))
print("Euler operator of (lhs - rhs):", sp.simplify(sp.expand(E)))
