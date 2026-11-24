"""Float: Gaarding/positivity function of the derivative-level operator L1 (acting on v = w') in a weight W_I:
   (Lw)' = L1 w' + a0' w,  L1 = -d^4 - k^2 d^2 - (7km + q/4) d + (-8 m^2 - 5/2 k p).
W_I: phi' = 2p x/(1+x^2) for x <= a2, blended (smoothstep on [a2, b2]) to the rho value (1/3)x(1+x^2)^(-1/3).
Positivity: int W(v''^2 - A v'^2 - V v^2) >= 1/2 int W v''^2 + int W mu v^2."""
import numpy as np, sympy as sp, sys
X, k, m, p, q, P = sp.symbols('X k m p q P')
a2s, b2s = sp.symbols('a2 b2')
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4}
def D(e): return sp.diff(e, X) + sum(sp.diff(e, v) * r for v, r in rules.items())
def make(phip):
    def rr2(f): return D(D(f)) - 2 * phip * D(f) + (phip**2 - D(phip)) * f      # (W f)''/W
    c2, c1, c0 = -k**2, -(7 * k * m + q / 4), -8 * m**2 - sp.Rational(5, 2) * k * p
    A = 2 * (phip**2 - D(phip)) - c2
    V = c0 - (D(c1) - phip * c1) / 2 + rr2(c2) / 2 - rr2(phip**2 - D(phip)) / 2
    mu = -A**2 / 2 - rr2(A) / 2 - V
    return mu, A, V
u = sp.symbols('u')
fp = 2 * P * X / (1 + X**2); fr = X * (1 + X**2) ** sp.Rational(-1, 3) / 3
mu_p, A_p, V_p = make(fp)                       # pure polynomial region
mu_r, A_r, V_r = make(fr)                       # rho-like tail
# blend region: s = S((X - a2)/(b2 - a2)), phi' = (1-s) fp + s fr
S_ = lambda t: t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3)
sb = S_((X - a2s) / (b2s - a2s)); mu_b, A_b, V_b = make((1 - sb) * fp + sb * fr)
fpol = sp.lambdify((X, k, m, p, q, P), [mu_p, A_p, V_p], 'numpy')
frho = sp.lambdify((X, k, m, p, q), [mu_r, A_r, V_r], 'numpy')
fbl = sp.lambdify((X, k, m, p, q, P, a2s, b2s), [mu_b, A_b, V_b], 'numpy')
T0 = 10.249073; st = -2.249073
Dm = np.loadtxt('mc_prof.txt')[::10]
x = (8 - st) * (1 - (Dm[:, 0] + Dm[:, 1]) / 2 / T0)
Pp = lambda i: (Dm[:, 2 + 2 * i] + Dm[:, 3 + 2 * i]) / 2
for Pv, a2, b2 in [(3.5, 30.25, 38.25), (3.0, 30.25, 38.25), (4.0, 30.25, 38.25)]:
    xs = np.concatenate([x[::-1], np.linspace(10.3, 60, 1000)]); z = np.zeros(1000)
    kk = np.concatenate([Pp(2)[::-1], z]); mm = np.concatenate([Pp(3)[::-1], z]); pp = np.concatenate([Pp(0)[::-1], z])
    qq = np.concatenate([Pp(1)[::-1], np.linspace(10.3, 60, 1000) + st])
    MU = np.where(xs <= a2, fpol(xs, kk, mm, pp, qq, Pv)[0],
         np.where(xs >= b2, frho(xs, kk, mm, pp, qq)[0], fbl(xs, kk, mm, pp, qq, Pv, a2, b2)[0]))
    sel = xs >= 3.0
    print(f"p={Pv}: min mu1 on x>=3: {MU[sel].min():+.4f} at x={xs[sel][np.argmin(MU[sel])]:.2f}; "
          + " ".join(f"{xx}:{MU[np.argmin(abs(xs-xx))]:+.3f}" for xx in [3, 4, 5, 6, 8, 10, 15, 20, 30, 34, 38, 45, 55]))
