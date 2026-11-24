"""Tail of lemma L-int: mu1(x) >= c x^{4/3} for x >= 68.25 (rho-type weight there), in the chart u = x^{-1/3},
w = (1+u^6)^{-1/3}, exactly as positivity.py (fourth supplement) but for the derivative-level operator L1."""
import sys, json, sympy as sp
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
ctx.prec = 120
u, w, k, m, p, q, st_ = sp.symbols('u w k m p q st_')
P1 = w / (3 * u); P2 = u**2 * (w / 3 - 2 * w**4 / 9); P3 = u**5 * (-2 * w**4 / 3 + 16 * w**7 / 27)
P4 = u**8 * (-2 * w**4 / 3 + 32 * w**7 / 9 - 224 * w**10 / 81)
Ps = sp.symbols('Q1 Q2 Q3 Q4 Q5'); kk, mm, pp, qq = sp.symbols('kk mm pp qq')
rules = {Ps[0]: Ps[1], Ps[1]: Ps[2], Ps[2]: Ps[3], Ps[3]: Ps[4], pp: -kk * qq, qq: 1 + kk * pp, kk: mm, mm: pp / 4}
def D(e): return sp.expand(sum(sp.diff(e, v) * r for v, r in rules.items()))
def rr2(f): return D(D(f)) - 2 * Ps[0] * D(f) + (Ps[0]**2 - Ps[1]) * f
c2, c1, c0 = -kk**2, -(7 * kk * mm + qq / 4), -8 * mm**2 - sp.Rational(5, 2) * kk * pp
A = 2 * (Ps[0]**2 - Ps[1]) - c2
V = c0 - (D(c1) - Ps[0] * c1) / 2 + rr2(c2) / 2 - rr2(Ps[0]**2 - Ps[1]) / 2
mu1 = sp.expand(-A**2 / 2 - rr2(A) / 2 - V)
M = sp.expand((u**4 * mu1).subs({Ps[0]: P1, Ps[1]: P2, Ps[2]: P3, Ps[3]: P4, kk: k, mm: m, pp: p, qq: u**-3 + st_}))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
src = _Pr().doprint(M); assert "u**" not in src
Mf = eval("lambda u, w, k, m, p, st_: " + src, {"arb": arb})
X0 = 68.25; c0_ = arb(3) / 1000
u0 = (arb(X0) ** (arb(-1) / 3)).upper(); N = 400; worst = None; ok = True
for j in range(N):
    ulo, uhi = u0 * j / N, u0 * (j + 1) / N
    U = (ulo + uhi) / 2 + arb(0, 1) * (uhi - ulo) / 2
    W = (1 + U**6) ** (arb(-1) / 3)
    sig_lo = 1 / (uhi * uhi * uhi) + H.STIP.lower()
    kb, mb, pb, qs = H.profile_bounds(sig_lo)
    val = Mf(U, W, arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper()), H.STIP + arb(0, qs.upper()))
    if not bool(val > c0_): ok = False
    if worst is None or val.lower() < worst.lower(): worst = val
print(f"mu1(x) >= {float(c0_)} x^(4/3) for all x >= {X0}: {ok};  smallest enclosure {worst}")
json.dump({"X0": X0, "c0": float(c0_), "ok": ok}, open('lint_tail.json', 'w'))
