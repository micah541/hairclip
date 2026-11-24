"""Rescaled CDF for normal graphs over S1:  Gamma(sigma,tau) = Gamma1(sigma) + u(sigma,tau) nu1(sigma).
<Gamma_tau, N> = -k_ss + 1/4 <Gamma, N>,  N = i T.   Here alpha = 1 - k1 u, beta = u', J = |Gamma_sigma|.
  u_tau = (J/alpha) [ -J^{-1} d/dsigma (J^{-1} d/dsigma k) + (alpha (p+u) - u' q)/(4 J) ],
  k = (k1 J^2 + alpha u'' - u' alpha')/J^3.
Checks: (i) u = 0 is stationary (shrinker equation); (ii) the linearization is
  L u = -u'''' - k^2 u'' - (5 k m + q/4) u' + (1/4 - 3 m^2 - k p) u."""
import sympy as sp
s, eps = sp.symbols('sigma epsilon')
p, q, k, m = [sp.Function(n)(s) for n in ('p', 'q', 'k', 'm')]
U = sp.Function('U')(s)
rules = {sp.Derivative(p, s): -k * q, sp.Derivative(q, s): 1 + k * p, sp.Derivative(k, s): m, sp.Derivative(m, s): p / 4}
def D(e):
    return sp.expand(sp.diff(e, s).subs(rules))
u = eps * U
alpha = 1 - k * u
up = D(u)
alphap = D(alpha)
J = sp.sqrt(alpha**2 + up**2)
kap = (k * J**2 + alpha * D(up) - up * alphap) / J**3
ks = D(kap) / J
kss = D(ks) / J
rhs = (J / alpha) * (-kss + (alpha * (p + u) - up * q) / (4 * J))
r0 = sp.simplify(rhs.subs(eps, 0))
print("u = 0 residual (should be 0):", r0)
lin = sp.expand(sp.diff(rhs, eps).subs(eps, 0))
Ud = [U] + [sp.Derivative(U, (s, n)) for n in range(1, 5)]
Lref = -Ud[4] - k**2 * Ud[2] - (5 * k * m + q / 4) * Ud[1] + (sp.Rational(1, 4) - 3 * m**2 - k * p) * Ud[0]
print("linearization - L (should be 0):", sp.simplify(lin - Lref))
quad = sp.expand(sp.diff(rhs, eps, 2).subs(eps, 0) / 2)
# structure of the quadratic part: coefficient of U'''' (quasilinear) and the rest
c4 = sp.simplify(quad.coeff(Ud[4]))
print("quadratic part, coefficient of U'''':", c4)
import pickle
pickle.dump({"rhs": rhs, "lin": lin, "quad": quad}, open("graph_eq.pkl", "wb"))
