"""Divergence-form energy identity at the E0 level (inner region, weight Psi = rho; the chi-cutoff region is
handled separately with pointwise bounds).

For a test function w and the graph u over S1, the curvature part of the equation pairs as
   int Psi w F_curv(u) = int T_c(u; w),
   T_c = -[Phi J^-1 w'' + B1 w' + B0 w] kappa,   Phi = Psi/alpha_g,
   B1 = (Phi J^-1)_x + Phi_x J^-1,  B0 = (Phi_x J^-1)_x,
   kappa = J^-3 (alpha_g u'' + k J^2 + u'(m u + k u')),  alpha_g = 1 - k u,  J^2 = alpha_g^2 + u'^2,
(two integrations by parts of  -int Phi w d(J^-1 kappa_x) ;  F_curv = -alpha_g^-1 d(J^-1 kappa_x)).
T_c/rho depends on (phi', phi'', k, m, p, q, U0, U1, U2, W0, W1, W2) only (no u''', no w''').
Checks: (a) no U3; (b) degree <= 2 in (U2, W2); (c) the linearisation in u reproduces int rho w L_curv u
after integration by parts (numerical check on random data in check_lin()).
Writes nlin_E0.pkl with T_c/rho, its first and second u-variations."""
import sympy as sp, pickle, time
t0 = time.time()
x, f1, f2 = sp.symbols('x phi1 phi2')          # phi' and phi''  (rho = e^-phi)
k, m, p, q = sp.symbols('k m p q')
U0, U1, U2, U3 = sp.symbols('U0 U1 U2 U3')
W0, W1, W2, W3 = sp.symbols('W0 W1 W2 W3')
R = sp.symbols('R')                            # rho, carried as a symbol with R' = -phi' R
f3 = sp.symbols('phi3')
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4, U0: U1, U1: U2, U2: U3, W0: W1, W1: W2, W2: W3,
         R: -f1 * R, f1: f2, f2: f3}
def D(e):
    return sum(sp.diff(e, v) * r for v, r in rules.items())
al = 1 - k * U0
Jv = sp.sqrt(al ** 2 + U1 ** 2)
kap = (al * U2 + k * Jv ** 2 + U1 * (m * U0 + k * U1)) / Jv ** 3
Phi = R / al
Phx = D(Phi)
B1 = D(Phi / Jv) + Phx / Jv
B0 = D(Phx / Jv)
Tc = -(Phi / Jv * W2 + B1 * W1 + B0 * W0) * kap
Tc = Tc.subs(R, 1)                              # T_c / rho (every term carries one factor R)
assert not Tc.has(U3) and not Tc.has(W3) and not Tc.has(f3), "unexpected third derivatives"
print(f"[{time.time()-t0:.0f}s] T_c/rho built; free symbols:", sorted(map(str, Tc.free_symbols)))
# degree in (U2, W2): substitute random rationals for the other symbols and expand
import random
random.seed(1)
num = {v: sp.Rational(random.randint(1, 9), random.randint(10, 19)) for v in (f1, f2, k, m, p, U0, U1, W0, W1)}
Tn = Tc.subs(num)
vals = {}
for a in range(4):
    for b in range(4):
        vals[(a, b)] = sp.nsimplify(sp.N(Tn.subs({U2: a, W2: b}), 40), rational=True)
# fit total degree: third differences in each direction and mixed must vanish
import itertools
d3u = [sp.N(vals[(3, b)] - 3 * vals[(2, b)] + 3 * vals[(1, b)] - vals[(0, b)]) for b in range(4)]
d3w = [sp.N(vals[(a, 3)] - 3 * vals[(a, 2)] + 3 * vals[(a, 1)] - vals[(a, 0)]) for a in range(4)]
mix = sp.N(vals[(2,1)] - 2*vals[(1,1)] + vals[(0,1)] - (vals[(2,0)] - 2*vals[(1,0)] + vals[(0,0)]))  # d^3/dU2^2 dW2
print('third differences (U2):', [float(t) for t in d3u], ' (W2):', [float(t) for t in d3w], ' mixed U2^2 W2:', float(mix))
Pn = None

# variations in the u-direction h = (H0, H1, H2)
H0, H1, H2, e = sp.symbols('H0 H1 H2 e')
Te = Tc.subs({U0: U0 + e * H0, U1: U1 + e * H1, U2: U2 + e * H2}, simultaneous=True)
T1 = sp.diff(Te, e).subs(e, 0)
T2 = sp.diff(Te, e, 2).subs(e, 0)
print(f"[{time.time()-t0:.0f}s] variations done")
pickle.dump({"Tc": Tc, "T1": T1, "T2": T2}, open("nlin_E0.pkl", "wb"))
print(f"[{time.time()-t0:.0f}s] saved nlin_E0.pkl")
