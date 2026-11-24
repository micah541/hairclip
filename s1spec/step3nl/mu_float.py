"""Float profile of the positivity function mu(x) (Gaarding constant of the rho-form) on [0, 34]."""
import numpy as np, sympy as sp
X, k, m, p, q, om = sp.symbols('X k m p q om')
phi = (1 + X**2) ** sp.Rational(2, 3) / 4
P1, P2, P3, P4 = [sp.diff(phi, X, i) for i in range(1, 5)]
A = 2 * (P1**2 - P2) + k**2
Ad = 2 * (2 * P1 * P2 - P3) + 2 * k * m
Add = 2 * (2 * P2**2 + 2 * P1 * P3 - P4) + 2 * (m**2 + k * p / 4)
V = (-P1**4 / 2 + 3 * P1**2 * P2 - P1**2 * k**2 / 2 - 2 * P1 * P3 - P1 * k * m / 2 - P1 * q / 8 - 3 * P2**2 / 2
     + P2 * k**2 / 2 + P4 / 2 - k * p / 2 - 3 * m**2 / 2 + sp.Rational(3, 8))
rhoA = Add - 2 * P1 * Ad + (P1**2 - P2) * A
mu = -A**2 / 2 - rhoA / 2 - V - om
f = sp.lambdify((X, k, m, p, q, om), [mu, A, V], 'numpy')
T0 = 10.249073; st = -2.249073
D = np.loadtxt('mc_prof.txt')[::20]
x = (8 - st) * (1 - (D[:, 0] + D[:, 1]) / 2 / T0)
P = lambda i: (D[:, 2 + 2 * i] + D[:, 3 + 2 * i]) / 2
pv, qv, kv, mv = P(0), P(1), P(2), P(3)
xs = np.concatenate([x[::-1], np.linspace(10.3, 34, 300)])
kk = np.concatenate([kv[::-1], 0 * np.linspace(10.3, 34, 300)]); mm = np.concatenate([mv[::-1], 0 * np.linspace(10.3, 34, 300)])
pp = np.concatenate([pv[::-1], 0 * np.linspace(10.3, 34, 300)]); qq = np.concatenate([qv[::-1], np.linspace(10.3, 34, 300) + st])
MU, AA, VV = f(xs, kk, mm, pp, qq, 3 / 64)
i = np.argmin(MU); print(f"min mu = {MU[i]:.4f} at x = {xs[i]:.3f};  max |A| = {np.abs(AA).max():.3f}; max V = {VV.max():.3f}")
for xx in [0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 16, 20, 25, 30]:
    j = np.argmin(abs(xs - xx)); print(f"  x={xs[j]:6.2f} mu={MU[j]:+8.4f} A={AA[j]:+7.4f} V={VV[j]:+7.4f}")
