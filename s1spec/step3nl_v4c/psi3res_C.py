"""Rigorous bound of ||(L - lamt) psi3^||_{L^2(rho)} over the stage-C range (sigma in [tip, 8]) from CAPD profile
hulls (MONALL boxes) and interval evaluation of the explicit piecewise-Chebyshev psi3^ (psi3md_14.json)."""
import sys, json, numpy as np
from flint import arb
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec/step3')
from gen_hs3 import T0, LAMT
J = json.load(open('psi3md_14.json'))
P = [(arb(pc["a"]), arb(pc["b"]), [np.array(pc["coef"])] + [Ch.chebder(np.array(pc["coef"]), n) for n in range(1, 7)]) for pc in J["pieces"]]
def clen(cf, y):
    b1, b2 = arb(0), arb(0)
    for ck in cf[:0:-1]:
        b1, b2 = 2 * y * b1 - b2 + arb(float(ck)), b1
    return y * b1 - b2 + arb(float(cf[0]))
def derivs(X, xc):
    """enclosures of u^(n)(X), n = 0..4, on the pieces meeting X (mean-value form around xc)"""
    out = None
    for a, b, cfs in P:
        if X.upper() < a.lower() or X.lower() > b.upper(): continue
        sc = -2 / (b - a)
        yc = 1 + sc * (arb(xc) - a); yX = 1 + sc * (X - a)
        vals = []
        for n in range(5):
            hh = X - arb(xc)
            vc = clen(cfs[n], yc) * sc ** n; v1 = clen(cfs[n + 1], yc) * sc ** (n + 1); d2 = clen(cfs[n + 2], yX) * sc ** (n + 2)
            vals.append(vc + v1 * hh + d2 * hh * hh / 2)
        out = vals if out is None else [o.union(v) for o, v in zip(out, vals)]
    return out
tot = arb(0); worst = arb(0); lam = arb(LAMT); T0a = arb(T0)
for line in open('/tmp/mc_prof.txt'):
    ta, tb, p0, p1, q0, q1, k0, k1, m0, m1, s0, s1 = map(float, line.split())
    xl = (8 - arb(s1)) * (1 - arb(tb) / T0a); xh = (8 - arb(s0)) * (1 - arb(ta) / T0a)
    xlo, xhi = max(float(xl.lower()), 0.0), float(xh.upper())
    X = arb((xlo + xhi) / 2, (xhi - xlo) / 2); xc = (xlo + xhi) / 2
    u = derivs(X, xc)
    I = lambda lo, hi: arb((lo + hi) / 2, (hi - lo) / 2)
    p, q, k, m = I(p0, p1), I(q0, q1), I(k0, k1), I(m0, m1)
    r = -u[4] - k * k * u[2] - (5 * k * m + q / 4) * u[1] + (arb(1) / 4 - 3 * m * m - k * p - lam) * u[0]
    rho = (-(1 + arb(xlo) ** 2) ** (arb(2) / 3) / 4).exp()
    ra = abs(r).upper() if hasattr(abs(r), 'upper') else abs(r)
    contrib = 2 * arb(xhi - xlo) * rho * arb(float(abs(r).upper())) ** 2
    tot += contrib
    w = arb(float(abs(r).upper())) * rho.sqrt()
    if float(w.upper()) > float(worst.upper()): worst = w; wx = xlo
print(f"stage C (x in [0, 10.249]): ||r||_rho^2 <= {float(tot.upper()):.3e}  =>  ||r||_rho <= {float(tot.upper())**0.5:.3e};  max weighted |r| <= {float(worst.upper()):.2e} at x ~ {wx:.3f}")
