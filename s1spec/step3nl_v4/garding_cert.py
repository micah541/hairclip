"""Rigorous Gaarding constants for the E2 level (inner region, weight rho, x in [0, 34.25]):
  (L w)'' = Mt w'' + e1 w' + e0 w,  Mt = -d^4 - k^2 d^2 - (9km + q/4) d + c0t   (garding2.py, exact),
  Re<Mt v, v>_rho <= -1/2 ||v''||^2 - int rho mu2 v^2   (mu2 = -A^2/2 - (rho A)''/(2 rho) - V),
and the level-0 analogue mu0 for reference.  Interval evaluation of mu2, e1, e0 on profile hull boxes (stage C:
prof_tight.txt grouped; stage B: sigma in [8, 32] with the step-1 chunk bounds).  Output garding_cert.json."""
import sys, json, time
import sympy as sp
from flint import arb
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
GROUP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
X, k, m, p, q = sp.symbols('X k m p q')
phi = (1 + X**2) ** sp.Rational(2, 3) / 4
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4}
def D(e): return sp.diff(e, X) + sum(sp.diff(e, v) * r for v, r in rules.items())
P1 = sp.diff(phi, X)
def rr2(f): return D(D(f)) - 2 * P1 * D(f) + (P1**2 - D(P1)) * f
def forms(c2, c1, c0):
    A = 2 * (P1**2 - D(P1)) - c2
    r4 = rr2(P1**2 - D(P1))
    V = c0 - (D(c1) - P1 * c1) / 2 + rr2(c2) / 2 - r4 / 2
    return -A**2 / 2 - rr2(A) / 2 - V
a2, a1, a0 = -k**2, -(5 * k * m + q / 4), sp.Rational(1, 4) - 3 * m**2 - k * p
c1t = 2 * D(a2) + a1; c0t = D(D(a2)) + 2 * D(a1) + a0
e1 = sp.expand(D(D(a1)) + 2 * D(a0)); e0 = sp.expand(D(D(a0)))
mu2 = forms(a2, c1t, c0t); mu0 = forms(a2, a1, a0)
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        return f"({bb}**(arb({ex.p})/arb({ex.q})))"
def mk(e): return eval("lambda X, k, m, p, q: " + _Pr().doprint(e), {"arb": arb})
fmu2, fmu0, fe1, fe0 = mk(mu2), mk(mu0), mk(e1), mk(e0)
def ev(Xb, pp, qq, kk, mm, reg, sup):
    vals = {"-mu2": -fmu2(Xb, kk, mm, pp, qq), "-mu0": -fmu0(Xb, kk, mm, pp, qq), "|e1|": abs(fe1(Xb, kk, mm, pp, qq)),
            "|e0|": abs(fe0(Xb, kk, mm, pp, qq))}
    d = sup.setdefault(reg, {})
    for key, v in vals.items():
        vu = float(v.upper()); assert vu == vu, 'NaN'
        d[key] = max(d.get(key, -1e300), vu)
def run():
    t0 = time.time(); sup = {}
    lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
    for g0 in range(0, len(lines), GROUP):
        grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
        ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
        hl = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        if Xb.lower() < 0: Xb = arb(0).union(arb(Xb.upper()))
        xm = float(Xb.mid())
        reg = "x<=1" if xm <= 1 else ("1<x<=3" if xm <= 3 else ("3<x<=12" if xm <= 12 else "12<x<=34.25"))
        ev(Xb, *hl, reg, sup)
    N = 50
    for ci in range(8, 32):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            Xb = (sa - HR.ST).union(sb - HR.ST)
            reg = "3<x<=12" if Xb.upper() <= 12 else "12<x<=34.25"
            ev(Xb, arb(0, pb.upper()), sa.union(sb) + arb(0, qs.upper()), arb(0, kb.upper()), arb(0, mb.upper()), reg, sup)
    print(f"[{time.time()-t0:.0f}s]")
    return sup
if __name__ == "__main__":
    sup = run()
    for reg in sorted(sup): print(f"  {reg:12s}", {kk: float(f"{v:.5g}") for kk, v in sup[reg].items()})
    json.dump({"group": GROUP, "sup": sup}, open('garding_cert.json', 'w'), indent=1)
