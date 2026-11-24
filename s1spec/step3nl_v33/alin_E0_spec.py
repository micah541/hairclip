"""Rigorous alpha-linear coefficients of the E0 identity (divergence form, core region chi = 1):
   D N(u0)[w; w] = int rho sum_ij a_ij(x; alpha) w^(i) w^(j),  a_ij = t_ij(u0) - t_ij(0),  u0 = alpha psi3^ + alpha^2 h2^,
where t_ij = d^2 T1/dH_i dW_j (T1 = first u-variation of T_c/rho, bilinear in (h; w)).
For |alpha| <= AMAX:  |a_ij|/|alpha| <= |c1(0)| + AMAX |c2(A)|  (arb power series in alpha: c1 at the exact point 0,
c2 enclosed over A = [-AMAX, AMAX]).  Boxes: CAPD hull boxes (grouped by GROUP) on x <= tau_tip, chunk-bound boxes
on sigma in [8, 29] (x up to 31.25).  Output alin_E0.json: per region the sup of |a_ij|/|alpha|."""
import sys, json, time, pickle
from flint import arb, arb_series
import sympy as sp
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 3e-3
GROUP = int(sys.argv[2]) if len(sys.argv) > 2 else 10
P = pickle.load(open('/home/claude/s1spec/step3nl/nlin_E0.pkl', 'rb'))
f1, f2, k, m, p = sp.symbols('phi1 phi2 k m p')
U = sp.symbols('U0 U1 U2'); W = sp.symbols('W0 W1 W2'); H = sp.symbols('H0 H1 H2')
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        if ex == sp.Rational(1, 2): return f"({bb}.sqrt())"
        if ex == -sp.Rational(1, 2): return f"(1/({bb}.sqrt()))"
        if ex.is_Rational and ex.q == 2:   # half-integer powers
            n = int(ex.p); base = f"({bb}.sqrt())"; prod = "*".join([base] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        raise ValueError(f"power {ex}")
T1 = P['T1']
TF = {}
for i in range(3):
    for j in range(3):
        e = sp.diff(T1, H[i], W[j])
        if e != 0:
            TF[(i, j)] = eval("lambda phi1, phi2, k, m, p, U0, U1, U2: " + _Pr().doprint(e), {"arb": arb})
A = arb(0, AMAX)
def coeff_bounds(y, h, ph1, ph2, kk, mm, pp):
    """SIGNED enclosures of a_ij(alpha)/alpha for |alpha| <= AMAX, then spectral data of the symmetrised form."""
    z = arb(0); c = {}
    u0s = [arb_series([z, y[j], h[j]], prec=3) for j in range(3)]
    uAs = [arb_series([A * y[j] + A * A * h[j], y[j] + 2 * A * h[j], h[j]], prec=3) for j in range(3)]
    for ij, f in TF.items():
        c0 = f(ph1, ph2, kk, mm, pp, *u0s).coeffs() + [arb(0)] * 3
        cA = f(ph1, ph2, kk, mm, pp, *uAs).coeffs() + [arb(0)] * 3
        c[ij] = c0[1] + arb(0, AMAX) * cA[2]          # a_ij/alpha in this ball for all |alpha| <= AMAX
    g = lambda i, j: c.get((i, j), arb(0))
    m00, m11 = g(0, 0), g(1, 1); m01 = (g(0, 1) + g(1, 0)) / 2
    tr = (m00 + m11) / 2; df = (m00 - m11) / 2
    up = lambda z: arb(abs(z).upper())                  # exact nonnegative upper bound (abs of a 0-straddling ball is not >= 0)
    adf, am = up(df), up(m01)
    lam = up(tr) + (adf * adf + am * am).sqrt()       # >= spectral radius of the symmetric 2x2 (w, w') block
    return {"spec01": lam, "|a02+a20|": abs(g(0, 2) + g(2, 0)), "|a12+a21|": abs(g(1, 2) + g(2, 1)), "|a22|": abs(g(2, 2)),
            "sum01": abs(g(0, 0)) + abs(g(0, 1) + g(1, 0)) + abs(g(1, 1))}
def phis(X):
    E = (1 + X * X) ** (arb(-1) / 3)
    return X * E / 3, E / 3 - 2 * X * X * E * E * E * E / 9
def run():
    t0 = time.time(); sup = {"x<=2": {}, "2-12": {}, "12-31.25": {}}
    def upd(reg, d):
        for ij, v in d.items():
            vu = float(v.upper()); assert vu == vu, 'NaN'
            sup[reg][ij] = max(sup[reg].get(ij, 0.0), vu)
    lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
    for g0 in range(0, len(lines), GROUP):
        grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
        ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
        hl = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        Xb = arb(0).union(Xb) if Xb.lower() < 0 else Xb
        pp_, qq_, kk_, mm_ = hl
        ph1, ph2 = phis(Xb)
        for Pf, Ph in zip(F.PIECES, Hf.PIECES):
            if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
            y = HR.jets(F, Pf, Xb, 3); h = HR.jets(Hf, Ph, Xb, 3)
            upd("x<=2" if Xb.upper() <= 2 else "2-12", coeff_bounds(y, h, ph1, ph2, kk_, mm_, pp_))
    print(f"  stage C done [{time.time()-t0:.0f}s]", flush=True)
    N = 50
    for ci in range(8, 29):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2); ph1, ph2 = phis(X)
            y = HR.jets_any(F, X, xc)[:3]; h = HR.jets_any(Hf, X, xc)[:3]
            reg = "2-12" if xb <= 12 else "12-31.25"
            upd(reg, coeff_bounds(y, h, ph1, ph2, arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper())))
    print(f"  stage B done [{time.time()-t0:.0f}s]", flush=True)
    return sup
if __name__ == "__main__":
    sup = run()
    for reg, d in sup.items():
        print(reg, {kk: float(f"{v:.4g}") for kk, v in sorted(d.items())})
    json.dump({"alpha_max": AMAX, "note": "sup over region of |alpha|^-1 x (spectral radius of sym (w,w') block | coefficient sums)", "sup": sup},
              open(f'alin_E0_spec_{AMAX:g}.json', 'w'), indent=1)
