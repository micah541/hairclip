"""E2-level nonlinear coefficients (expanded form, paired with w'''' after two integrations by parts on the test side):
   D := F(u0 + w) - F(u0) - L w = sum_{t,j,l} d_tjl(x) w''^t w'''^j w''''^l  +  d_000,
   F = sum_mon c_mon(U0, U1) u2^i u3^j u4^l + (p + u - u'q/alpha_g)/4      (graph_struct.pkl monomials),
   u0 = alpha psi3^ + alpha^2 h2^,  U0 = u0 + w0, U1 = u0' + w1  (w0, w1 in the jet box [-eta0, eta0] x [-eta1, eta1]).
 d_000 = [F_lower(U0, U1; u0'', u0''', u0'''') - F_lower(u0, u0'; ...)] - (L3 part on w0, w1), bounded by its mean-value form
   |d_000| <= A0 |w| + A1 |w'|,  A_i = sup over the box of |d/dU_i (F - L)(U0, U1; u0'', u0''', u0'''')|.
Output (per region): A0, A1 and sup|d_tjl| for (t,j,l) in {(1,0,0), (0,1,0), (0,0,1), (2,0,0), (1,1,0), (3,0,0)}.
Usage: python3 nlin_E2.py AMAX ETA0 ETA1 [GROUP]"""
import sys, json, time, pickle
from flint import arb
import sympy as sp
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX, ETA0, ETA1 = (float(a) for a in sys.argv[1:4])
GROUP = int(sys.argv[4]) if len(sys.argv) > 4 else 50
GS = pickle.load(open('/home/claude/s1spec/step3nl/graph_struct.pkl', 'rb'))
Ji, al, k, m, p, q, u0, u1 = sp.symbols('Ji al k m p q u0 u1')
U0, U1, a2, a3, a4, w2, w3, w4 = sp.symbols('U0 U1 a2 a3 a4 w2 w3 w4')   # a_i: u0-jets (alpha psi3^ + alpha^2 h2^)^(i)
J = sp.sqrt((1 - k * U0) ** 2 + U1 ** 2)
sub = {Ji: 1 / J, al: 1 - k * U0, u0: U0, u1: U1}
Fpert = sum(cf.subs(sub) * (a2 + w2) ** i * (a3 + w3) ** j * (a4 + w4) ** l for (i, j, l), cf in GS.items())
Fpert += (p + U0 - U1 * q / (1 - k * U0)) / 4
Lw = -w4 - k ** 2 * w2                      # parts of L acting on w2..w4 (L has no w3 term); w0, w1 parts go into d_000
# polynomial coefficients in (w2, w3, w4) at fixed (U0, U1)
Pw = sp.Poly(sp.expand(Fpert - Lw - Fpert.subs({w2: 0, w3: 0, w4: 0})), w2, w3, w4)
coef = {mon: c for mon, c in zip(Pw.monoms(), Pw.coeffs())}
print("w-monomials:", sorted(coef))
# d_000 = Flow(U0,U1) - Flow(u00,u01) - L01(w0, w1); its U-derivatives (mean value form):
Flow = Fpert.subs({w2: 0, w3: 0, w4: 0})
L01 = lambda W0, W1: -(5 * k * m + q / 4) * W1 + (sp.Rational(1, 4) - 3 * m ** 2 - k * p) * W0
dU0 = sp.diff(Flow, U0) - (sp.Rational(1, 4) - 3 * m ** 2 - k * p)
dU1 = sp.diff(Flow, U1) + (5 * k * m + q / 4)
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
def box_vals(y, h, pp, qq, kk, mm):
    uj = [A * y[j] + A * A * h[j] for j in range(5)]
    Ub0, Ub1 = uj[0] + E0b, uj[1] + E1b
    out = {mon: abs(f(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm)) for mon, f in CFN.items()}
    out["A0"] = abs(FU0(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm))
    out["A1"] = abs(FU1(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm))
    return out
def run():
    t0 = time.time(); sup = {}
    def upd(reg, d):
        s = sup.setdefault(reg, {})
        for key, v in d.items():
            vu = float(v.upper()); assert vu == vu, 'NaN'
            s[str(key)] = max(s.get(str(key), 0.0), vu)
    lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
    for g0 in range(0, len(lines), GROUP):
        grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
        ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
        hl = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        if Xb.lower() < 0: Xb = arb(0).union(arb(Xb.upper()))
        for Pf, Ph in zip(F.PIECES, Hf.PIECES):
            if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
            y = HR.jets(F, Pf, Xb, 5); h = HR.jets(Hf, Ph, Xb, 5)
            upd("x<=12", box_vals(y, h, *hl))
    N = 50
    for ci in range(8, 32):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2)
            y = HR.jets_any(F, X, xc); h = HR.jets_any(Hf, X, xc)
            reg = "x<=12" if xb <= 12 else "12<x<=34.25"
            upd(reg, box_vals(y, h, arb(0, pb.upper()), sa.union(sb) + arb(0, qs.upper()), arb(0, kb.upper()), arb(0, mb.upper())))
    print(f"[{time.time()-t0:.0f}s]")
    return sup
if __name__ == "__main__":
    sup = run()
    for reg in sorted(sup): print(reg, {kk: float(f"{v:.4g}") for kk, v in sup[reg].items()})
    json.dump({"alpha_max": AMAX, "eta0": ETA0, "eta1": ETA1, "sup": sup}, open(f'nlin_E2_{AMAX:g}_{ETA0:g}_{ETA1:g}.json', 'w'), indent=1)
