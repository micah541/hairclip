"""Untrusted part of the psi3 certificate: shoot L u = lamt u from the tip (u=1, u'=u'''=0, u''=s), pick s so the
growing Airy pair is absent (match the finite-difference eigenvector), then fit piecewise Chebyshev on [0,24]."""
import sys, json, numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP
LT = -0.2561532224577463
def prof(x):
    s = x + STIP
    return PROF(s) if s <= 22 else np.array([0.0, s, 0.0, 0.0])
def rhs(x, U):
    p, q, k, m = prof(x); u, u1, u2, u3 = U
    return [u1, u2, u3, -k * k * u2 - (5 * k * m + q / 4) * u1 + (0.25 - 3 * m * m - k * p - LT) * u]
def shoot(s, X):
    return solve_ivp(rhs, [0, X], [1.0, 0.0, s, 0.0], method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
S = np.load('slowman.npz'); x = S['x']; ps = S['psi3']
s0 = np.gradient(np.gradient(ps, x), x)[0]
# the growing modes ~ exp(0.236 x^(4/3)): bisection on s for the sign of u(X) - FD value at X = 16
X = 16.0; target = np.interp(X, x, ps)
f = lambda s: shoot(s, X).y[0, -1] - target
a, b = s0 - 0.5, s0 + 0.5
fa = f(a)
for _ in range(80):
    c = 0.5 * (a + b); fc = f(c)
    if np.sign(fc) == np.sign(fa): a, fa = c, fc
    else: b = c
sol = shoot(0.5 * (a + b), X)
xx = np.linspace(0, 12, 1201); uu = sol.sol(xx)[0]
dev = np.max(np.abs(uu - np.interp(xx, x, ps)) / (1 + np.abs(np.interp(xx, x, ps))))
pieces = []
for lo in np.arange(0, 12, 1.0):
    t = np.linspace(lo, lo + 1, 200); c = Ch.chebfit(2 * (t - lo) - 1, sol.sol(t)[0], 24)
    pieces.append({"lo": float(lo), "coef": c.tolist()})
json.dump({"s": 0.5 * (a + b), "lamt": LT, "pieces": pieces}, open('psi3hat.json', 'w'))
print(f"u''(0) = {0.5*(a+b):.12f}  (FD estimate {s0:.6f});  relative deviation from FD eigenvector on [0,12]: {dev:.2e}")
print("u at x = 2, 5, 10:", np.round(sol.sol([2, 5, 10])[0], 4), " (FD: -0.452, ~5.0, 40.8)")
