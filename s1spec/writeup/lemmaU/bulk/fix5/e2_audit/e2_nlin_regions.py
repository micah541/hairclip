"""E2 audit: E2-level nonlinear coefficients d_tjl, A0, A1 (as step3nl/nlin_E2.py) on REGIONS WITH THEIR OWN JET BOXES.
   D := F(u0 + w) - F(u0) - L w = sum d_tjl w''^t w'''^j w''''^l + d_000,  |d_000| <= A0 |w| + A1 |w'|,
   u0 = alpha psi3 + alpha^2 h2, alpha in [-AMAX, AMAX], (w, w') in [-eta0, eta0] x [-eta1, eta1] on the region.
Differences to nlin_E2.py (same machinery, read-only imports of h2res):
  (1) regions and boxes are arguments: REG = name:xlo:xhi:eta0:eta1 ; a profile box is evaluated for every region it intersects;
  (2) the drift part of A1 is differentiated analytically:  d/dU1[-U1 q/(4(1-kU0))] + q/4 = -(q/4) k U0/(1-k U0)
      (nlin_E2.py evaluates q/4 - q/(4(1-kU0)) in balls: dependency artefact 5e-3 in the far field); identity checked in sympy below;
  (3) PAD=1 widens the u0-jets by the certified sup|psi3~ - psi3^|^(j), |h2~ - h2^|^(j) (c4switch_S1.json, sup_abs), so the result
      covers the C^4 functions psi3~, h2~ of the v4 decomposition as well as the C^3 ones.
usage: python3 -B e2_nlin_regions.py AMAX PAD OUT.json REG [REG ...]"""
import sys, json, time, pickle, hashlib, os
from flint import arb
import sympy as sp
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX = float(sys.argv[1]); PAD = int(sys.argv[2]); OUT = sys.argv[3]
REGS = []
for a in sys.argv[4:]:
    n, lo, hi, e0, e1 = a.split(':'); REGS.append((n, float(lo), float(hi), float(e0), float(e1)))
GROUP = 50
P_GS = '/home/claude/s1spec/step3nl/graph_struct.pkl'; P_PROF = '/home/claude/s1_tight/prof/prof_tight.txt'
P_C4 = '/home/claude/s1spec/step3nl_v4/agent_int/c4switch/c4switch_S1.json'
GS = pickle.load(open(P_GS, 'rb'))
Ji, al, k, m, p, q, u0, u1 = sp.symbols('Ji al k m p q u0 u1')
U0, U1, a2, a3, a4, w2, w3, w4 = sp.symbols('U0 U1 a2 a3 a4 w2 w3 w4')
J = sp.sqrt((1 - k * U0) ** 2 + U1 ** 2)
sub = {Ji: 1 / J, al: 1 - k * U0, u0: U0, u1: U1}
Fnd = sum(cf.subs(sub) * (a2 + w2) ** i * (a3 + w3) ** j * (a4 + w4) ** l for (i, j, l), cf in GS.items())      # F without the drift
drift = (p + U0 - U1 * q / (1 - k * U0)) / 4
Lw = -w4 - k ** 2 * w2
Pw = sp.Poly(sp.expand(Fnd - Lw - Fnd.subs({w2: 0, w3: 0, w4: 0})), w2, w3, w4)       # the drift has no w2..w4
coef = {mon: c for mon, c in zip(Pw.monoms(), Pw.coeffs())}
Flow = Fnd.subs({w2: 0, w3: 0, w4: 0})
# analytic drift derivatives minus the drift part of L01 (= -(q/4) W1 + (1/4) W0):
dr0 = -U1 * q * k / (4 * (1 - k * U0) ** 2)                # d/dU0 drift - 1/4
dr1 = -(q / 4) * k * U0 / (1 - k * U0)                     # d/dU1 drift + q/4
assert sp.simplify(sp.diff(drift, U0) - sp.Rational(1, 4) - dr0) == 0 and sp.simplify(sp.diff(drift, U1) + q / 4 - dr1) == 0
dU0 = sp.diff(Flow, U0) + (3 * m ** 2 + k * p) + dr0       # d/dU0 (F - L)
dU1 = sp.diff(Flow, U1) + 5 * k * m + dr1                  # d/dU1 (F - L)
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
A = arb(0, AMAX)
if PAD:
    c4 = json.load(open(P_C4)); PY = [arb(0, v * (1 + 1e-12)) for v in c4['psi3']['sup_abs']]; PH = [arb(0, v * (1 + 1e-12)) for v in c4['h2']['sup_abs']]
else:
    PY = PH = [arb(0)] * 5
def box_vals(y, h, pp, qq, kk, mm, e0, e1):
    uj = [A * (y[j] + PY[j]) + A * A * (h[j] + PH[j]) for j in range(5)]
    Ub0, Ub1 = uj[0] + arb(0, e0), uj[1] + arb(0, e1)
    out = {str(mon): abs(f(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm)) for mon, f in CFN.items()}
    out["A0"] = abs(FU0(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm))
    out["A1"] = abs(FU1(Ub0, Ub1, uj[2], uj[3], uj[4], pp, qq, kk, mm))
    out["sup|kU0|"] = abs(kk * Ub0); out["sup|U1|"] = abs(Ub1)
    return out
def run():
    t0 = time.time(); sup = {r[0]: {} for r in REGS}; nb = {r[0]: 0 for r in REGS}
    def upd(reg, d):
        s = sup[reg]; nb[reg] += 1
        for key, v in d.items():
            vu = float(v.upper()); assert vu == vu, 'NaN'
            s[key] = max(s.get(key, 0.0), vu)
    def regs_of(xlo, xhi):
        return [r for r in REGS if xhi >= r[1] and xlo <= r[2]]
    lines = [l for l in open(P_PROF) if not l.startswith('END')]
    for g0 in range(0, len(lines), GROUP):
        grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
        ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
        hl = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        if Xb.lower() < 0: Xb = arb(0).union(arb(Xb.upper()))
        rr = regs_of(float(Xb.lower()), float(Xb.upper()))
        for Pf, Ph in zip(F.PIECES, Hf.PIECES):
            if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
            y = HR.jets(F, Pf, Xb, 5); h = HR.jets(Hf, Ph, Xb, 5)
            for r in rr: upd(r[0], box_vals(y, h, *hl, r[3], r[4]))
    N = 50
    for ci in range(8, 32):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2)
            y = HR.jets_any(F, X, xc); h = HR.jets_any(Hf, X, xc)
            for r in regs_of(xa, xb):
                upd(r[0], box_vals(y, h, arb(0, pb.upper()), sa.union(sb) + arb(0, qs.upper()), arb(0, kb.upper()), arb(0, mb.upper()), r[3], r[4]))
    return sup, nb, time.time() - t0
if __name__ == "__main__":
    sup, nb, secs = run()
    for r in REGS:
        print(r, 'boxes', nb[r[0]]); print('   ', {kk: float(f"{v:.6g}") for kk, v in sup[r[0]].items()})
    sha = lambda pth: hashlib.sha256(open(pth, 'rb').read()).hexdigest()
    ins = [P_GS, P_PROF, P_C4, '/home/claude/s1spec/step3nl/h2res.py', '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py',
           '/home/claude/s1spec/step3nl/psi3md2_14.json', '/home/claude/s1spec/step3nl/h2md_14.json', '/home/claude/s1spec/step3/hstail.py', os.path.abspath(__file__)]
    json.dump({"alpha_max": AMAX, "pad_c4": PAD, "regions": [dict(name=r[0], xlo=r[1], xhi=r[2], eta0=r[3], eta1=r[4], nboxes=nb[r[0]]) for r in REGS], "sup": sup,
               "label": "cert (arb interval sups over hull boxes; inputs taken as certified)",
               "args": dict(argv=sys.argv, env={k_: os.environ.get(k_) for k_ in ("PYTHONHASHSEED",)}, inputs={pth: sha(pth) for pth in ins}), "secs": round(secs, 1)},
              open(OUT, 'w'), indent=1)
