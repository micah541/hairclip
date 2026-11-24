"""T3: nonlinear coefficient bounds on [12,20],[20,30.25],[30.25,40],[40,60] (adapted from nlin_E2.py).
 D := F(u0+w) - F(u0) - L w = sum d_tjl w''^t w'''^j w''''^l + d_000,  u0 = alpha psi3^ + alpha^2 h2^, |alpha| <= AMAX,
 jet box |w| <= ETA0, |w'| <= ETA1;  d_000 bounded by A0 |w| + A1 |w'| (mean value, sup over the box of the U-derivatives).
 Also sup|u0'|, sup|u0''| (alpha-box), and sup|psi3^'|, |psi3^''|, |h2^'|, |h2^''|.   Usage: python3 int_nlin.py AMAX ETA0 ETA1 [width_scale]"""
import sys, pickle
from int_common import *
AMAX, ETA0, ETA1 = (float(a) for a in sys.argv[1:4])
WS = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
GS = pickle.load(open('/home/claude/s1spec/step3nl/graph_struct.pkl', 'rb'))
Ji, al, k, m, p, q, u0, u1 = sp.symbols('Ji al k m p q u0 u1')
U0, U1, a2, a3, a4, w2, w3, w4 = sp.symbols('U0 U1 a2 a3 a4 w2 w3 w4')
J = sp.sqrt((1 - k * U0) ** 2 + U1 ** 2)
sub = {Ji: 1 / J, al: 1 - k * U0, u0: U0, u1: U1}
Fpert = sum(cf.subs(sub) * (a2 + w2) ** i * (a3 + w3) ** j * (a4 + w4) ** l for (i, j, l), cf in GS.items())
Fpert += (p + U0 - U1 * q / (1 - k * U0)) / 4
Lw = -w4 - k ** 2 * w2
Pw = sp.Poly(sp.expand(Fpert - Lw - Fpert.subs({w2: 0, w3: 0, w4: 0})), w2, w3, w4)
coef = {mon: c for mon, c in zip(Pw.monoms(), Pw.coeffs())}
print("w-monomials:", sorted(coef))
Flow = Fpert.subs({w2: 0, w3: 0, w4: 0})
dU0 = sp.diff(Flow, U0) - (sp.Rational(1, 4) - 3 * m ** 2 - k * p)
dU1 = sp.diff(Flow, U1) + (5 * k * m + q / 4)
# dU1 is affine in q and the q-part cancels identically on the ray (coefficient O(k)): split it to avoid the width/4 dependency artefact
assert sp.simplify(sp.diff(dU1, q, 2)) == 0
dU1_q = sp.factor(sp.simplify(sp.diff(dU1, q))); print('dU1_q =', dU1_q)
dU1_0 = dU1.subs(q, 0)
_chk = {U0: 0.37, U1: 0.21, k: 0.013, m: 0.02, p: -0.3, q: 7.3, a2: 0.1, a3: 0.2, a4: 0.3}
assert abs(float((dU1 - (dU1_0 + q * dU1_q)).subs(_chk))) < 1e-12
dU1 = dU1_0 + q * dU1_q
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        if ex.is_Rational and ex.q == 2:
            n = int(ex.p); base = f"({bb}.sqrt())"; prod = "*".join([base] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        raise ValueError(ex)
ARG = "U0, U1, a2, a3, a4, p, q, k, m"
mk = lambda e: eval(f"lambda {ARG}: " + _Pr().doprint(e), {"arb": arb})
CFN = {mon: mk(c) for mon, c in coef.items()}
FU0, FU1 = mk(dU0), mk(dU1)
A = arb(0, AMAX); E0b, E1b = arb(0, ETA0), arb(0, ETA1)
def absup(v): return chk(abs(v)).upper()
def box_vals(y, h, pp, qq, kk, mm):
    uj = [A * y[j] + A * A * h[j] for j in range(5)]
    Ub0, Ub1 = uj[0] + E0b, uj[1] + E1b
    # J^2 = (1 - k U0)^2 + U1^2 > 0 : explicit guard (sqrt of a ball must not contain 0)
    al_ = 1 - kk * Ub0; J2 = al_ * al_ + Ub1 * Ub1; assert J2.lower() > 0, 'J^2 ball touches 0'
    out = {mon: abs(f(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm)) for mon, f in CFN.items()}
    out["A0"] = abs(FU0(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm))
    out["A1"] = abs(FU1(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm))
    out["|u0'|"] = abs(uj[1]); out["|u0''|"] = abs(uj[2]); out["|u0|"] = abs(uj[0])
    out["|psi3'|"] = abs(y[1]); out["|psi3''|"] = abs(y[2]); out["|h2'|"] = abs(h[1]); out["|h2''|"] = abs(h[2])
    out["|u0'''|"] = abs(uj[3]); out["|u0''''|"] = abs(uj[4])
    return out
def run():
    t0 = time.time(); sup = {}
    bps = [('x', 12.0), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)]
    bl = make_segments(bps, [0.02 * WS, 0.025 * WS, 0.025 * WS, 0.05 * WS])
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        y = jets_any(F, b.X, b.xc); h = jets_any(Hf, b.X, b.xc)
        for v in y + h: chk(v)
        d = box_vals(y, h, b.p, b.q, b.k, b.m)
        for ri in region_of(b.xa, b.xb):
            s = sup.setdefault(f"[{REG[ri][0]:g},{REG[ri][1]:g}]", {})
            for key, v in d.items():
                vu = up(v); s[str(key)] = max(s.get(str(key), 0.0), vu)
        if n % 200 == 0: print(n, len(bl), b.xa, f"{time.time()-t0:.0f}s", flush=True)
    print(f"[{time.time()-t0:.0f}s]")
    return sup
if __name__ == "__main__":
    sup = run()
    for reg in sup: print(reg, {kk: float(f"{v:.4g}") for kk, v in sup[reg].items()})
    json.dump({"alpha_max": AMAX, "eta0": ETA0, "eta1": ETA1, "width_scale": WS, "sup": sup}, open(f'int_T3_{AMAX:g}_{ETA0:g}_{ETA1:g}_{WS:g}.json', 'w'), indent=1)
