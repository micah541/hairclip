"""Rigorous residual bound for psi3^ (psi3md2_14.json + explicit far extension) on sigma in [8, 32] (stage B),
with the tail-theorem chunk bounds |k|<=Bk, |m|<=Bm, |p|<=Bp, |q - sigma|<=Bq (gen_hs2.chunks) and sigma_tip in
[-2.249127, -2.249019].  r = r_ray + delta, r_ray = -u'''' - ((x + st)/4) u' + (1/4 - lamt) u."""
import sys, json, numpy as np
from flint import arb
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec/step3')
from gen_hs3 import LAMT
from gen_hs2 import chunks
J = json.load(open('psi3md2_14.json')); Xc = J["X"]; c1, c2 = arb(J["c_poly"]), arb(J["c_dec"])
P = [(arb(pc["a"]), arb(pc["b"]), [np.array(pc["coef"])] + [Ch.chebder(np.array(pc["coef"]), n) for n in range(1, 7)]) for pc in J["pieces"]]
ST = arb((-2.249127 - 2.249019) / 2, (2.249127 - 2.249019) / 2); lam = arb(LAMT); A = 1 - 4 * lam
aj = [arb(1)]
for j in range(1, 4): aj.append(aj[-1] * (A - 4 * j + 4) * (A - 4 * j + 3) * (A - 4 * j + 2) * (A - 4 * j + 1) / j)
def clen(cf, y):
    b1, b2 = arb(0), arb(0)
    for ck in cf[:0:-1]: b1, b2 = 2 * y * b1 - b2 + arb(float(ck)), b1
    return y * b1 - b2 + arb(float(cf[0]))
def derivs_poly(X, xc):
    out = None
    for a, b, cfs in P:
        if X.upper() < a.lower() or X.lower() > b.upper(): continue
        sc = -2 / (b - a); yc = 1 + sc * (arb(xc) - a); yX = 1 + sc * (X - a); hh = X - arb(xc)
        vals = [clen(cfs[n], yc) * sc**n + clen(cfs[n + 1], yc) * sc**(n + 1) * hh + clen(cfs[n + 2], yX) * sc**(n + 2) * hh * hh / 2 for n in range(5)]
        out = vals if out is None else [o.union(v) for o, v in zip(out, vals)]
    return out
def derivs_ext(X):
    """explicit extension c1*g_poly + c2*e^{S - S(Xc)}, S' = -(s/4)^(1/3), s = x + st (Xs fixed at the junction)"""
    s = X + ST; out = []
    for d in range(5):
        v = arb(0)
        for j in range(4):
            e = A - 4 * j; f = arb(1)
            for i in range(d): f = f * (e - i)
            v += aj[j] * f * s ** (e - d)
        out.append(c1 * v)
    # decaying part: |c2 e^S| <= |c2| and its derivatives <= |c2| * (|S'|+1)^d (crude, S decreasing beyond Xc)
    z = (s / 4) ** (arb(1) / 3)
    return [o + arb(0, float((abs(c2) * (abs(z) + 1) ** d).upper())) for d, o in enumerate(out)]
ch = chunks(); tot = arb(0); worst = 0.0
for ci, (ln, hstep, bnds) in enumerate(ch):
    s_hi = 32.0 - ci * ln; s_lo = s_hi - ln
    Bk, Bm, Bp, Bq = [arb(float(max(abs(lo), abs(hi)))) for lo, hi in bnds]
    N = 100
    for i in range(N):
        sa, sb = s_lo + i * ln / N, s_lo + (i + 1) * ln / N
        xa = float((arb(sa) - ST).lower()); xb = float((arb(sb) - ST).upper())
        X = arb((xa + xb) / 2, (xb - xa) / 2)
        u = derivs_poly(X, (xa + xb) / 2) if xb <= Xc else derivs_ext(X)
        if u is None: u = derivs_ext(X)
        rr = -u[4] - ((X + ST) / 4) * u[1] + (arb(1) / 4 - lam) * u[0]
        dl = Bk * Bk * abs(u[2]) + (5 * Bk * Bm + Bq / 4) * abs(u[1]) + (3 * Bm * Bm + Bk * Bp) * abs(u[0])
        rb = float((abs(rr) + dl).upper())
        rho = (-(1 + arb(xa) ** 2) ** (arb(2) / 3) / 4).exp()
        tot += 2 * arb(xb - xa) * rho * arb(rb) ** 2
        worst = max(worst, rb * float(rho.sqrt().upper()))
print(f"stage B (sigma in [8,32]): ||r||_rho <= {float(tot.upper())**0.5:.3e}; max weighted |r| <= {worst:.2e}")
