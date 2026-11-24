"""float coarse scan of mu1 (level 0) and mu2 (level 2, weight (1+(x/60)^4)W) for Plan B weight on [10.25, 260]."""
import sympy as sp, numpy as np, sys
x = sp.symbols('x'); t = (x - 40) / 20
S = t**4 * (35 - 84*t + 70*t**2 - 20*t**3)
Pm = {'lo': sp.Integer(2), 'mid': 2 - S, 'hi': sp.Integer(1)}
P1s, P2s, P3s, P4s, P5s, k, m, p, q = sp.symbols('P1 P2 P3 P4 P5 k m p q')
rules = {P1s: P2s, P2s: P3s, P3s: P4s, P4s: P5s, p: -k*q, q: 1 + k*p, k: m, m: p/4}
def D(e): return sp.expand(sum(sp.diff(e, v)*r for v, r in rules.items()))
def rr2(f): return D(D(f)) - 2*P1s*D(f) + (P1s**2 - P2s)*f
c2, c1, c0 = -k**2, -(7*k*m + q/4), -8*m**2 - sp.Rational(5, 2)*k*p
A = 2*(P1s**2 - P2s) - c2
V = c0 - (D(c1) - P1s*c1)/2 + rr2(c2)/2 - rr2(P1s**2 - P2s)/2
mu1 = sp.expand(-A**2/2 - rr2(A)/2 - V)
f1 = sp.lambdify((P1s, P2s, P3s, P4s, k, m, p, q), mu1, 'numpy')
dphi = {r: sp.lambdify(x, [sp.diff(2*Pm[r]*x/(1 + x**2), x, j) for j in range(4)], 'numpy') for r in Pm}
def mu1_at(xv, kk=0, mm=0, pp=0):
    r = 'lo' if xv <= 40 else ('mid' if xv <= 60 else 'hi')
    P = dphi[r](xv); qq = xv - 2.2490776511
    return f1(*P, kk, mm, pp, qq)
if __name__ == "__main__":
    xs = np.concatenate([np.arange(10.25, 60, 0.05), np.arange(60, 262, 0.5)])
    v = np.array([mu1_at(a) for a in xs])
    i = v.argmin(); print("min mu1", v[i], "at", xs[i])
    for a in (10.25, 12, 16, 20, 30, 40, 45, 50, 55, 60, 80, 100, 150, 200, 260): print(a, mu1_at(a), (2*(2 if a < 40 else 1) - 1)/8)
    # coarse local minima by region
    for lo, hi in ((10.25, 16), (16, 40), (40, 60), (60, 120), (120, 260)):
        mk = (xs >= lo) & (xs <= hi); j = v[mk].argmin(); print(f"[{lo},{hi}] min {v[mk][j]:.4f} at {xs[mk][j]:.2f}")
