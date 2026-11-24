"""T4: H^2-level Gaarding constants for the derivative-level operator L1 in the intermediate weight W = exp(-psi).
 L1 v = -v'''' - k^2 v'' - (7km + q/4) v' - (8m^2 + (5/2) k p) v,   v = w'.
 (L1 v)'' = M1 v'' + e1 v' + e0 v,  M1 = -d^4 + a2 d^2 + c1t d + c0t  (acting on V = v''),
   a2 = -k^2, a1 = -(7km+q/4), a0 = -8m^2-(5/2)kp, c1t = 2D(a2)+a1, c0t = D(D(a2))+2D(a1)+a0, e1 = D(D(a1))+2D(a0), e0 = D(D(a0)),
   D = total derivative along the profile ODE p' = -kq, q' = 1+kp, k' = m, m' = p/4.
 Re<M1 V, V>_W <= -1/2 ||V''||_W^2 - int W mu2 V^2,  mu2 = -A^2/2 - (W A)''/(2W) - Vpot  (garding_cert.py forms(), phi -> psi);
 W-derivatives: rr2(f) = (W f)''/W = D(D f) - 2 P1 D f + (P1^2 - P2) f with P_j = psi^(j), D P_j = P_{j+1}.
 mu1 (level-0 form for L1, lint_cert.py) is evaluated as a check (certified there: mu1 >= 0.2435 on [12, 30.25]).
Usage: python3 int_garding.py [width_scale]"""
import sys
from int_common import *
WS = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
P1, P2, P3, P4, P5, k, m, p, q = sp.symbols('P1 P2 P3 P4 P5 k m p q')
rules = {P1: P2, P2: P3, P3: P4, P4: P5, p: -k * q, q: 1 + k * p, k: m, m: p / 4}
def D(e): return sp.expand(sum(sp.diff(e, v) * r for v, r in rules.items()))
def rr2(f): return D(D(f)) - 2 * P1 * D(f) + (P1**2 - P2) * f
def forms(c2, c1, c0):
    A = 2 * (P1**2 - P2) - c2
    V = c0 - (D(c1) - P1 * c1) / 2 + rr2(c2) / 2 - rr2(P1**2 - P2) / 2
    return sp.expand(-A**2 / 2 - rr2(A) / 2 - V)
a2, a1, a0 = -k**2, -(7 * k * m + q / 4), -8 * m**2 - sp.Rational(5, 2) * k * p
c1t = 2 * D(a2) + a1; c0t = D(D(a2)) + 2 * D(a1) + a0
e1 = sp.expand(D(D(a1)) + 2 * D(a0)); e0 = sp.expand(D(D(a0)))
mu1 = forms(a2, a1, a0); mu2 = forms(a2, c1t, c0t)
assert not mu1.has(P5) and not mu2.has(P5), 'P5 appears'
print("mu2 free symbols:", mu2.free_symbols, " e1:", e1, " e0:", e0)
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
mk = lambda e: eval("lambda P1, P2, P3, P4, k, m, p, q: " + _Pr().doprint(e), {"arb": arb})
fmu1, fmu2, fe1, fe0 = mk(mu1), mk(mu2), mk(e1), mk(e0)
def run():
    t0 = time.time(); sup = {}
    bps = [('x', 12.0), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)]
    bl = make_segments(bps, [0.01 * WS, 0.01 * WS, 0.01 * WS, 0.02 * WS])
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        P = Pj_box(b.X)
        vals = {"-mu2": -fmu2(*P, b.k, b.m, b.p, b.q), "|e1|": abs(fe1(*P, b.k, b.m, b.p, b.q)), "|e0|": abs(fe0(*P, b.k, b.m, b.p, b.q)),
                "-mu1": -fmu1(*P, b.k, b.m, b.p, b.q)}
        for ri in region_of(b.xa, b.xb):
            s = sup.setdefault(f"[{REG[ri][0]:g},{REG[ri][1]:g}]", {})
            for key, v in vals.items():
                vu = up(v); s[key] = max(s.get(key, -1e300), vu)
        if n % 500 == 0: print(n, len(bl), b.xa, f"{time.time()-t0:.0f}s", flush=True)
    print(f"[{time.time()-t0:.0f}s]")
    return sup
if __name__ == "__main__":
    sup = run()
    for reg in sup: print(reg, {kk: float(f"{v:.5g}") for kk, v in sup[reg].items()})
    json.dump({"width_scale": WS, "sup": sup}, open(f'int_T4_{WS:g}.json', 'w'), indent=1)
