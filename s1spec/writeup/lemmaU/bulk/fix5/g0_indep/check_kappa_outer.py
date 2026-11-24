"""G0_INDEP: arb check of kappa_out = sup_[57,60] (|zeta'| + psi' zeta), zeta = S((60-x)/3), S(t) = t^4(35-84t+70t^2-20t^3), psi' = x(1+x^2)^(-1/3)/3 (rho-type weight, x >= 38.25)."""
from flint import arb, ctx
ctx.prec = 100
best = arb(0); n = 30000
for i in range(n):
    x = arb(57 + (i + 0.5) * 3 / n, 1.5 / n); t = (60 - x) / 3
    S = t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3); dS = 140 * t**3 * (1 - t)**3
    v = dS / 3 + x * (1 + x * x) ** (arb(-1) / 3) / 3 * S
    if v.upper() > best.upper(): best = v.upper()
print("sup_[57,60] (|zeta'| + psi' zeta) <=", float(best.upper()))
