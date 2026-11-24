"""Positivity at infinity (fourth supplement, Lemma 'positivity at infinity'), made rigorous.
Exact:  Q(u) = int rho (u''^2 - A u'^2 - (V+omega) u^2) >= 1/2 int rho u''^2 + int rho mu u^2,
        mu = -A^2/2 - (rho A)''/(2 rho) - V - omega,   (rho A)''/rho = A'' - 2 phi' A' + (phi'^2 - phi'') A.
Certify mu(x) >= c0 x^{4/3} for all |x| >= R0 in the chart u = x^{-1/3}: M(u) = u^4 mu is polynomial in (u, w, profile),
w = (1+u^6)^{-1/3}; profile from the step-1 tail bounds (valid for sigma >= 8)."""
import sys, json
import numpy as np, sympy as sp
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
ctx.prec = 120
u, w, k, m, p, q, om = sp.symbols('u w k m p q om')
P1 = w / (3 * u); P2 = u**2 * (w / 3 - 2 * w**4 / 9); P3 = u**5 * (-2 * w**4 / 3 + 16 * w**7 / 27)
P4 = u**8 * (-2 * w**4 / 3 + 32 * w**7 / 9 - 224 * w**10 / 81)
A = 2 * (P1**2 - P2) + k**2
Ad = 2 * (2 * P1 * P2 - P3) + 2 * k * m                       # A' (profile ODE: k' = m)
Add = 2 * (2 * P2**2 + 2 * P1 * P3 - P4) + 2 * (m**2 + k * p / 4)   # A''  (m' = p/4)
V = (-P1**4 / 2 + 3 * P1**2 * P2 - P1**2 * k**2 / 2 - 2 * P1 * P3 - P1 * k * m / 2 - P1 * q / 8 - 3 * P2**2 / 2
     + P2 * k**2 / 2 + P4 / 2 - k * p / 2 - 3 * m**2 / 2 + sp.Rational(3, 8))
rhoA = Add - 2 * P1 * Ad + (P1**2 - P2) * A
import os
_FM = os.environ.get('FORM_MULT', '1/1').split('/'); fm = sp.Rational(int(_FM[0]), int(_FM[1]))
_OQ = os.environ.get('OMEGA_Q', '3/64').split('/')
mu = -(fm * A)**2 / 2 - fm * rhoA / 2 - fm * (V + om)
M = sp.expand(u**4 * mu)
# sanity: chain-rule identities for the u-chart P_i (P_{i+1} = d/dx P_i, dx = -3 u^4 du ... check numerically)
Xs = sp.symbols('X', positive=True); phi = (1 + Xs**2) ** sp.Rational(2, 3) / 4
Pref = [sp.diff(phi, Xs, i) for i in range(1, 5)]
for i, Pi in enumerate([P1, P2, P3, P4]):
    val = Pi.subs({w: (1 + u**6) ** sp.Rational(-1, 3)}).subs(u, sp.Rational(1, 2))
    ref = Pref[i].subs(Xs, 8)
    assert abs(float(val - ref)) < 1e-12, (i, float(val), float(ref))
print("u-chart P_i verified at x = 8")
M0 = sp.lambdify((u, w, k, m, p, q, om), M, 'numpy')
OM = int(_OQ[0]) / int(_OQ[1])
for x in [10.25, 12, 14, 16, 18, 20, 25, 30, 40, 60, 100, 1e3, 1e5]:
    uu = x ** (-1 / 3); ww = (1 + uu**6) ** (-1 / 3)
    print(f"x = {x:9.2f}:  mu/x^(4/3) = {M0(uu, ww, 0, 0, 0, x - 2.249073, OM):+.6f}   mu = {M0(uu, ww, 0, 0, 0, x - 2.249073, OM) * x**(4/3):+.4f}")

# ------------------------------------------------------------------ certification: M(u) >= c0 on 0 < u <= R0^{-1/3}
st_ = sp.symbols('st_')
Mq = sp.expand(M.subs(q, u**-3 + st_))            # q = x + sigma_tip + dq = u^-3 + st_  (st_ = sigma_tip + dq)
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp()
        if ex.is_Integer:
            n = int(ex); bb = f"({self._print(b)})"
            prod = "*".join([bb] * abs(n)) if n != 0 else "arb(1)"
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        return super()._print_Pow(e, rational)
src = _Pr().doprint(Mq)
assert "u**" not in src
Mf = eval("lambda u, w, k, m, p, st_, om: " + src, {"arb": arb})
R0, c0 = 30, arb(int(os.environ.get("POS_C0_MILLI", "4"))) / 1000
u0 = (arb(R0) ** (arb(-1) / 3)).upper()
N = 400
OMr = arb(int(_OQ[0])) / int(_OQ[1])
ok = True; worst = None
for j in range(N):
    ulo, uhi = u0 * j / N, u0 * (j + 1) / N
    U = (ulo + uhi) / 2 + arb(0, 1) * (uhi - ulo) / 2
    W = (1 + U * U * U * U * U * U) ** (arb(-1) / 3)
    sig_lo = (1 / (uhi * uhi * uhi)) + H.STIP.lower()          # smallest sigma on the piece
    kb, mb, pb, qs = H.profile_bounds(sig_lo)
    K_, Mm, P_ = arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper())
    ST = H.STIP + arb(0, qs.upper())
    val = Mf(U, W, K_, Mm, P_, ST, OMr)
    if not bool(val > c0): ok = False
    if worst is None or bool(val.lower() < worst[0].lower()): worst = (val, j)
print(f"certified: mu(x) >= {float(c0)} x^(4/3) for all x >= {R0}:", ok, "   smallest enclosure", worst[0], "on piece", worst[1])
json.dump({"R0": R0, "c0": float(c0), "ok": ok, "pieces": N, "form_mult": str(fm), "omega": "/".join(_OQ)}, open("positivity_cert" + os.environ.get("CERT_TAG", "") + ".json", "w"))
print("worst piece lower/upper bound:", float(worst[0].lower()), float(worst[0].upper()), " (c0 =", float(c0), ")")
