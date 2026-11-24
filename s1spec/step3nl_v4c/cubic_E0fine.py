"""Rigorous cubic(+) coefficients of the E0 identity (divergence form, core chi = 1):
   R3 = int_0^1 (1-t) T2(u0 + t w)[w, w; w] dt,  T2 = D_u^2 (T_c/rho) (nlin_E0.pkl), u0 = alpha psi3^ + alpha^2 h2^.
T2 = sum s_{ijl,d}(phi', phi'', k, m, p, U0, U1) H_i H_j W_l U2^d (polynomial in H, W, U2; d <= 2).
With H = W = w-jets and U2 = u0'' + t w'':  R3 <= (1/2) sum_mono C_mono |w|^a |w'|^b |w''|^c  pointwise, where C_mono
collects |s| C(d,e) |u0''|^(d-e) over the jet box U0 in u0 + [-eta0, eta0], U1 in u0' + [-eta1, eta1] (t in [0,1]).
Output: per region (x<=2, 2<x<=12, 12<x<=31.25) the sup of C_mono for every monomial (a, b, c), a + b + c >= 3."""
import sys, json, time, pickle
from math import comb
from flint import arb
import sympy as sp
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX, ETA0, ETA1 = (float(a) for a in sys.argv[1:4])
GROUP = int(sys.argv[4]) if len(sys.argv) > 4 else 100
P = pickle.load(open('/home/claude/s1spec/step3nl/nlin_E0.pkl', 'rb'))
f1, f2, k, m, p = sp.symbols('phi1 phi2 k m p')
U0, U1, U2 = sp.symbols('U0 U1 U2'); W = sp.symbols('W0 W1 W2'); H = sp.symbols('H0 H1 H2')
T2 = P['T2']
zero = {**{h: 0 for h in H}, **{w_: 0 for w_ in W}, U2: 0}
terms = []      # (i, j, l, d, sympy coefficient)
for i in range(3):
    for j in range(i, 3):
        for l in range(3):
            e = sp.diff(T2, H[i], H[j], W[l])
            if e == 0: continue
            mult = 1 if i == j else 2          # T2 = sum_{i,j} (full double sum); the (i,j),(j,i) pair counted once here
            # NOTE: diff wrt H_i H_j of sum_{a,b} c_ab H_a H_b gives c_ij + c_ji; for i == j it gives 2 c_ii.
            # The quadratic form value is sum_{i<=j} (d^2/dH_i dH_j) H_i H_j / (2 if i == j else 1).
            fac = sp.Rational(1, 2) if i == j else 1
            for d in range(3):
                cd = sp.diff(e, U2, d).subs(zero) / sp.factorial(d)
                cd = cd * fac
                if cd != 0: terms.append((i, j, l, d, cd))
print("nonzero (i,j,l,d):", [(t[0], t[1], t[2], t[3]) for t in terms])
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
FN = [(i, j, l, d, eval("lambda phi1, phi2, k, m, p, U0, U1: " + _Pr().doprint(cd), {"arb": arb})) for (i, j, l, d, cd) in terms]
A = arb(0, AMAX); E0b, E1b = arb(0, ETA0), arb(0, ETA1)
def phis(Xv):
    E = (1 + Xv * Xv) ** (arb(-1) / 3)
    return Xv * E / 3, E / 3 - 2 * Xv * Xv * E * E * E * E / 9
def box(y, h, ph1, ph2, kk, mm, pp, reg, sup):
    uj = [A * y[j] + A * A * h[j] for j in range(3)]
    Ub0, Ub1 = uj[0] + E0b, uj[1] + E1b
    acc = {}
    for (i, j, l, d, f) in FN:
        s = abs(f(ph1, ph2, kk, mm, pp, Ub0, Ub1))
        for e in range(d + 1):
            mono = [0, 0, 0]
            for o in (i, j, l): mono[o] += 1
            mono[2] += e
            c = s * comb(d, e)
            for _ in range(d - e): c = c * abs(uj[2])      # no ball**n: NaN for balls containing 0
            key = tuple(mono); acc[key] = acc.get(key, arb(0)) + c
    dd = sup.setdefault(reg, {})
    for key, v in acc.items():
        vu = float(v.upper()); assert vu == vu, 'NaN coefficient bound'
        dd[str(key)] = max(dd.get(str(key), 0.0), vu)
BRK = [0, 1, 2, 4, 8, 12, 16.25, 31.25]
def regions(lo, hi):
    return [f"{a}<x<={b}" for a, b in zip(BRK[:-1], BRK[1:]) if hi > a and lo <= b]
def run():
    t0 = time.time(); sup = {}
    lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
    for g0 in range(0, len(lines), GROUP):
        grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
        ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
        pp, qq, kk, mm = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        if Xb.lower() < 0: Xb = arb(0).union(arb(Xb.upper()))
        regs = regions(float(Xb.lower()), float(Xb.upper()))
        ph1, ph2 = phis(Xb)
        for Pf, Ph in zip(F.PIECES, Hf.PIECES):
            if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
            y = HR.jets(F, Pf, Xb, 3); h = HR.jets(Hf, Ph, Xb, 3)
            for reg in regs: box(y, h, ph1, ph2, kk, mm, pp, reg, sup)
    N = 50
    for ci in range(8, 29):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2); ph1, ph2 = phis(X)
            y = HR.jets_any(F, X, xc)[:3]; h = HR.jets_any(Hf, X, xc)[:3]
            for reg in regions(xa, xb):
                box(y, h, ph1, ph2, arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper()), reg, sup)
    print(f"[{time.time()-t0:.0f}s]")
    return sup
if __name__ == "__main__":
    sup = run()
    for reg in sup:
        print(reg, {kk: float(f"{v:.3g}") for kk, v in sorted(sup.get(reg, {}).items())})
    json.dump({"alpha_max": AMAX, "eta0": ETA0, "eta1": ETA1, "note": "R3 <= (1/2) sum C_mono |w|^a|w'|^b|w''|^c", "sup": sup},
              open(f'cubic_E0fine_{AMAX:g}_{ETA0:g}_{ETA1:g}.json', 'w'), indent=1)
