"""Structure of the rescaled graph equation (curvature part)  F = -(J/alpha) J^{-1} d/ds(J^{-1} d/dsigma kappa).
J, 1/J, alpha, 1/alpha are carried as symbols with their derivative rules, so F is a polynomial in
(u2, u3, u4) with coefficients polynomial in (Ji = 1/J, ai = 1/alpha, alpha, k, m, p, q, u0, u1)."""
import sympy as sp, pickle, time
t0 = time.time()
k, m, p, q = sp.symbols('k m p q')
u0, u1, u2, u3, u4, u5 = sp.symbols('u0:6')
Jv, Ji, al, ai = sp.symbols('Jv Ji al ai')
alp = -m * u0 - k * u1                       # alpha'
Jp = (al * alp + u1 * u2) * Ji               # J'
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4, u0: u1, u1: u2, u2: u3, u3: u4, u4: u5,
         al: alp, ai: -ai**2 * alp, Jv: Jp, Ji: -Ji**2 * Jp}
def D(e):
    return sp.expand(sum(sp.diff(e, v) * r for v, r in rules.items()))
def red(e):   # use Jv*Ji = 1, al*ai = 1, Jv^2 = al^2 + u1^2 (eliminate Jv)
    e = sp.expand(e).subs(Jv, 1 / Ji)
    return sp.expand(e)
kap = (k * (al**2 + u1**2) + al * u2 - u1 * alp) * Ji**3
ks = red(D(kap) * Ji)
kss = red(D(ks) * Ji)
F = red(-(1 / Ji) * ai * kss)                 # (J/alpha)(-kss)
F = sp.expand(F.subs(ai, 1 / al))
assert not F.has(u5)
P = sp.Poly(F, u2, u3, u4)
print(f"[{time.time()-t0:.0f}s] monomials (u2^i u3^j u4^l):")
out = {}
for mon, cf in sorted(P.terms(), key=lambda t: t[0]):
    cf = sp.factor(cf); out[mon] = cf
    print("  ", mon, ":", cf)
pickle.dump(out, open("graph_struct.pkl", "wb"))
