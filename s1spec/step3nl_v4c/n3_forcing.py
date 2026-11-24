"""Rigorous bound of the third-order forcing on the second-order slow manifold u0 = alpha psi3^ + alpha^2 h2^:
   g(alpha) = N0(alpha psi3^ + alpha^2 h2^),  g(0) = g'(0) = 0,  g''(0)/2 = Q2(psi3^),
   |g(alpha) - alpha^2 Q2(psi3^)| <= |alpha|^3 (|c3(0)| + amax |c4(0)| + amax^2 sup_A |c5|), c_j = g^(j)/j! (pointwise in x).
g'''(xi)/6 is enclosed as the t^3 coefficient of the power series of N0 at alpha = A + t, A = [-alpha_max, alpha_max]
(python-flint arb_series), with the jets of psi3^, h2^ and the profile enclosed on each box.
N0 = F_curv + drift - L u, F_curv from the exact monomial structure (graph_struct.pkl), drift = (p + u - u'q/alpha_g)/4.
Output: C3(alpha_max) = || sup |g'''|/6 ||_{L^2(rho; [0, 34.25])}  (the inner cutoff chi_32 vanishes beyond 34.25)."""
import sys, json, time, pickle
from flint import arb, arb_series
import sympy as sp
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 3e-3
STRIDE = int(sys.argv[2]) if len(sys.argv) > 2 else 1
PREC = 6
# ---- F_curv coefficient functions (Ji = 1/J, al = alpha_g)
GS = pickle.load(open('/home/claude/s1spec/step3nl/graph_struct.pkl', 'rb'))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"({e.p}/{e.q})"
    def _print_Integer(self, e): return f"({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
CF = {mon: eval("lambda Ji, al, k, m, p, q, u0, u1: " + _Pr().doprint(cf)) for mon, cf in GS.items()}
def N0(u, p, q, k, m):
    """u = [u0..u4] (arb_series or arb); profile arbs"""
    al = 1 - k * u[0]
    J2 = al * al + u[1] * u[1]
    Ji = J2.rsqrt() if hasattr(J2, 'rsqrt') else 1 / J2.sqrt()
    Fc = 0
    for (i, j, l), f in CF.items():
        term = f(Ji, al, k, m, p, q, u[0], u[1])
        for _ in range(i): term = term * u[2]
        for _ in range(j): term = term * u[3]
        for _ in range(l): term = term * u[4]
        Fc = Fc + term
    drift = (p + u[0] - u[1] * q / al) / 4
    lin = -u[4] - k * k * u[2] - (5 * k * m + q / 4) * u[1] + (arb(1) / 4 - 3 * m * m - k * p) * u[0]
    return Fc + drift - lin
A = arb(0, AMAX)
def g3_bound(y, h, p, q, k, m):
    """sup_{|alpha| <= AMAX} |g(alpha) - alpha^2 Q2| / |alpha|^3 <= |c3(0)| + AMAX |c4(0)| + AMAX^2 sup_{xi in A} |c5(xi)|
    (c_j = g^(j)/j!; the coefficients at the exact point 0 avoid the dependency blow-up of an interval expansion point)"""
    z = arb(0)
    u0s = [arb_series([z, y[j], h[j]], prec=PREC) for j in range(5)]
    c0 = N0(u0s, p, q, k, m).coeffs() + [arb(0)] * PREC
    uAs = [arb_series([A * y[j] + A * A * h[j], y[j] + 2 * A * h[j], h[j]], prec=PREC) for j in range(5)]
    cA = N0(uAs, p, q, k, m).coeffs() + [arb(0)] * PREC
    return abs(c0[3]) + AMAX * abs(c0[4]) + AMAX * AMAX * abs(cA[5])
def check_Q2():
    """float sanity: series coefficient of alpha^2 at alpha = 0 equals Q2(y)"""
    y = [arb(0.9), arb(0.3), arb(1.2), arb(-0.4), arb(0.7)]; h = [arb(0.1)] * 5
    p, q, k, m = arb(-2.4), arb(0.3), arb(1.5), arb(0.05)
    us = [arb_series([arb(0), y[j], h[j]], prec=PREC) for j in range(5)]
    c = N0(us, p, q, k, m).coeffs()
    q2 = HR.q2('Q2', y, p, q, k, m)
    print("sanity: N0 series at alpha=0:", [float(t.mid()) for t in c], "  Q2(y) =", float(q2.mid()))
def stage_C(path='/home/claude/s1_tight/prof/prof_tight.txt'):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0); nb = 0
    for li, line in enumerate(open(path)):
        if line.startswith('END'): break
        if li % STRIDE: continue
        v = list(map(float, line.split())); ta, tb, tc = v[0], v[1], v[2]
        hl = [HR.I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        xlo = max(float(Xb.lower()), 0.0); xhi = float(Xb.upper())
        p, q, k, m = hl
        best = arb(0)
        for Pf, Ph in zip(F.PIECES, Hf.PIECES):
            if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
            y = HR.jets(F, Pf, Xb, 5); h = HR.jets(Hf, Ph, Xb, 5)
            g = g3_bound(y, h, p, q, k, m)
            best = arb(max(best.upper(), g.upper()))
        rho = HR.rho_of(xlo); ln = arb(xhi - xlo) * STRIDE
        tot += ln * rho * best * best; nb += 1
        w = float(best.upper()) * float(rho.sqrt().upper())
        if w > worst[0]: worst = (w, xlo)
    print(f"C: {nb} boxes: ||g'''/6||^2 <= {float(tot.upper()):.4e}; max weighted {worst[0]:.3e} at x={worst[1]:.3f} [{time.time()-t0:.0f}s]", flush=True)
    return tot, worst
def stage_B(N=50):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0)
    for ci in range(8, 32):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            p, k, m = arb(0, pb.upper()), arb(0, kb.upper()), arb(0, mb.upper())
            q = sa.union(sb) + arb(0, qs.upper())
            xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2)
            y = HR.jets_any(F, X, xc); h = HR.jets_any(Hf, X, xc)
            g = g3_bound(y, h, p, q, k, m)
            rho = HR.rho_of(xa)
            tot += arb(xb - xa) * rho * g * g
            w = float(g.upper()) * float(rho.sqrt().upper())
            if w > worst[0]: worst = (w, xa)
    print(f"B: ||g'''/6||^2 <= {float(tot.upper()):.4e}; max weighted {worst[0]:.3e} at x={worst[1]:.3f} [{time.time()-t0:.0f}s]", flush=True)
    return tot, worst
if __name__ == "__main__":
    check_Q2()
    totC, wC = stage_C(); totB, wB = stage_B()
    C3 = (totC + totB).sqrt()
    out = {"alpha_max": AMAX, "C3": float(C3.upper()), "C2part": float(totC.upper()), "B2part": float(totB.upper()),
           "worstC": wC, "worstB": wB, "stride": STRIDE}
    if STRIDE == 1: json.dump(out, open(f'n3_forcing_{AMAX:g}.json', 'w'), indent=1)
    print(f"||chi (N0(alpha psi3^ + alpha^2 h2^) - alpha^2 Q2(psi3^))||_rho <= {float(C3.upper()):.4e} |alpha|^3  for |alpha| <= {AMAX}"
          + ("" if STRIDE == 1 else f"   (stride {STRIDE}: stage C is an estimate only)"))
