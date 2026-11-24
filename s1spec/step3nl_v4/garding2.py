"""Gaarding functions for the energies E0 = ||w||^2 and E2 = ||w''||^2 in L^2(rho) (half line, even w).
For M = -d^4 + c2 d^2 + c1 d + c0:  Re<Mv,v> = int rho(-v''^2 + A v'^2 + V v^2),
   A = 2(phi'^2 - phi'') - c2,  V = c0 - (c1' - phi' c1)/2 + (rho c2)''/(2 rho) - rho''''/(2 rho),
and  Re<Mv,v> <= -1/2 ||v''||^2 - int rho mu v^2,  mu = -A^2/2 - (rho A)''/(2 rho) - V.
Level 0: M = L.  Level 2: (Lw)'' = Mt w'' + e1 w' + e0 w with Mt = -d^4 - k^2 d^2 - (9km + q/4) d + c0t."""
import numpy as np, sympy as sp
X, k, m, p, q = sp.symbols('X k m p q')
phi = (1 + X**2) ** sp.Rational(2, 3) / 4
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4}
def D(e):  # total x-derivative (profile ODE + explicit X)
    return sp.diff(e, X) + sum(sp.diff(e, v) * r for v, r in rules.items())
P1 = sp.diff(phi, X)
def forms(c2, c1, c0):
    A = 2 * (P1**2 - D(P1)) - c2
    rho_rel = lambda f, n: None
    # (rho f)''/rho = f'' - 2 phi' f' + (phi'^2 - phi'') f ;  rho''''/rho via repeated use
    def rr2(f): return D(D(f)) - 2 * P1 * D(f) + (P1**2 - D(P1)) * f
    r4 = rr2(P1**2 - D(P1))   # rho''''/rho = ((rho'')/rho ...)  computed as (rho*g)''/rho with g = rho''/rho
    V = c0 - (D(c1) - P1 * c1) / 2 + rr2(c2) / 2 - r4 / 2
    mu = -A**2 / 2 - rr2(A) / 2 - V
    return A, V, mu
a2, a1, a0 = -k**2, -(5 * k * m + q / 4), sp.Rational(1, 4) - 3 * m**2 - k * p
A0, V0, mu0 = forms(a2, a1, a0)
c1t = 2 * D(a2) + a1; c0t = D(D(a2)) + 2 * D(a1) + a0
e1 = D(D(a1)) + 2 * D(a0); e0 = D(D(a0))
print("Mt first-order coefficient:", sp.simplify(c1t), "  zeroth:", sp.simplify(c0t))
A2, V2, mu2 = forms(a2, c1t, c0t)
f = sp.lambdify((X, k, m, p, q), [mu0, mu2, A0, V0, V2, e1, e0], 'numpy')
T0 = 10.249073; st = -2.249073
Dm = np.loadtxt('mc_prof.txt')[::10]
x = (8 - st) * (1 - (Dm[:, 0] + Dm[:, 1]) / 2 / T0)
P = lambda i: (Dm[:, 2 + 2 * i] + Dm[:, 3 + 2 * i]) / 2
xs = np.concatenate([x[::-1], np.linspace(10.3, 34.25, 400)])
z = np.zeros(400)
kk = np.concatenate([P(2)[::-1], z]); mm = np.concatenate([P(3)[::-1], z]); pp = np.concatenate([P(0)[::-1], z])
qq = np.concatenate([P(1)[::-1], np.linspace(10.3, 34.25, 400) + st])
M0, M2, AA, VV0, VV2, E1, E0 = f(xs, kk, mm, pp, qq)
print(f"level 0: max(-mu0) = {(-M0).max():.4f} at x={xs[np.argmax(-M0)]:.2f}   (without omega)")
print(f"level 2: max(-mu2) = {(-M2).max():.4f} at x={xs[np.argmax(-M2)]:.2f};  max|e1| = {np.abs(E1).max():.4f}, max|e0| = {np.abs(E0).max():.4f}")
for xx in [0, 1, 2, 4, 8, 12, 20, 30, 34]:
    j = np.argmin(abs(xs - xx)); print(f"  x={xs[j]:6.2f}  -mu0={-M0[j]:+8.4f}  -mu2={-M2[j]:+8.4f}  e1={E1[j]:+.4f} e0={E0[j]:+.4f}")
