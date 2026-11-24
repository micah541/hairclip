"""Rigorous 8x8 Gram matrix (half-line L^2(rho)) of  V = (L*f1, L*f2, gA, gB, D f1, D f2)  and  F = (f1, f2).
Same boxes (hulls of the Lyapunov CAPD run outA5.txt: 2400 'MB' rows in sigma in [8,32] with chunk profile bounds, 1025 'MC'
rows near the tip) and the same symbolic operators as /home/claude/s1spec/step3/lyap_consts2.py (read-only reuse by copying).
G_ij = sum_boxes  (dX rho)(box) * u_i(box) * u_j(box)   (ball product of the pointwise enclosures).
Conventions: C' = k S, S' = -k C;  psi2 = -sin(theta) = Ct C + St S;  cos(theta) = CSIGN (Ct S - St C) with CSIGN fixed below
(checked numerically against the float theta profile in check_frame.py)."""
import sys, json, time
import sympy as sp
from math import comb
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import gen_hs3 as G3
from gen_hs2 import chunks
from gram_post7 import DSet
import h2res as HR
F = HR.F
ctx.prec = 100
U = arb(0, 1)
CSIGN = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] in ('1', '-1') else None
def A(lo, hi): return arb((lo + hi) / 2) + U * ((hi - lo) / 2 * (1 + 1e-15) + 1e-300)
def hull(a, b):
    lo = min(a.lower(), b.lower()); hi = max(a.upper(), b.upper()); return (lo + hi) / 2 + U * ((hi - lo) / 2)
R_ = lambda p, q: arb(p) / arb(q)
X, P1, P2, P3, P4, k, m, q, p, lam = sp.symbols('X P1 P2 P3 P4 k m q p lam')
d = sp.symbols('d0:5')
Bp = [sp.Integer(1), P1, P1**2 + P2, P1**3 + 3*P1*P2 + P3, P1**4 + 6*P1**2*P2 + 3*P2**2 + 4*P1*P3 + P4]
Em = [sp.Integer(1), -P1, P1**2 - P2, -P1**3 + 3*P1*P2 - P3, P1**4 - 6*P1**2*P2 + 3*P2**2 + 4*P1*P3 - P4]
def der(n): return sum(comb(n, j) * Em[n - j] * d[j] for j in range(n + 1))   # e^phi (e^{-phi} g)^(n)
Ls = sp.expand(-der(4) - k**2 * der(2) + (k*m + q/4) * der(1) + d[0] / 2)       # = rho^-1 [...] of the task statement
aS = {n: sp.expand(Ls.coeff(d[n])) for n in range(5)}
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"R_({e.p},{e.q})"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp()
        if ex.is_Integer:
            n = int(ex); bb = f"({self._print(b)})"
            prod = "*".join([bb] * abs(n)) if n != 0 else "arb(1)"
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        return super()._print_Pow(e, rational)
def fn(expr, args):
    return eval("lambda " + ",".join(str(a) for a in args) + ": " + _Pr().doprint(sp.expand(expr)), {"arb": arb, "R_": R_})
VARS = (P1, P2, P3, P4, k, m, q, p)
cS = {n: fn(aS[n], VARS) for n in range(5)}
Cs, Ss = sp.symbols('Cs Ss')
rules = {p: -k*q, q: 1 + k*p, k: m, m: p/4, Cs: k*Ss, Ss: -k*Cs}
def D(e): return sp.expand(sum(sp.diff(e, v) * rv for v, rv in rules.items()))
def ders(e):
    out = [e]
    for _ in range(4): out.append(D(out[-1]))
    return out
PV = (k, m, q, p, Cs, Ss)
dp = [fn(e, PV) for e in ders(p)]; dC = [fn(e, PV) for e in ders(Cs)]; dS = [fn(e, PV) for e in ders(Ss)]
def weight(Xv):
    E = (1 + Xv * Xv) ** (arb(-1) / 3)
    X2 = Xv * Xv; E2 = E * E; E4 = E2 * E2; E7 = E4 * E2 * E; E10 = E7 * E2 * E
    P = [Xv * E / 3, E / 3 - 2 * X2 * E4 / 9, -2 * Xv * E4 / 3 + 16 * X2 * Xv * E7 / 27,
         -2 * E4 / 3 + 32 * X2 * E7 / 9 - 224 * X2 * X2 * E10 / 81]
    return P, (-(1 + Xv * Xv) ** (arb(2) / 3) / 4).exp()
Su = sp.symbols('u'); Sm = Su**4 * (35 - 84*Su + 70*Su**2 - 20*Su**3)
dSm_pts = [fn(sp.diff(Sm, Su, j), (Su,)) for j in range(8)]
def centered(ders_pts, I, j):
    m_ = arb(I.mid()); dl = I - m_; tot_ = arb(0); pw = arb(1); fac = 1
    for i in range(0, len(ders_pts) - j):
        if i > 0: pw = pw * dl; fac *= i
        tot_ += ders_pts[j + i](m_) * pw / fac
    return tot_
def _pw(x, j):
    r = arb(1)
    for _ in range(j): r = r * x
    return r
def leib(chi, ps, n): return sum((comb(n, j) * chi[j] * ps[n - j] for j in range(n + 1)), arb(0))
NAMES = ["Lf1", "Lf2", "gA", "gB", "Df1", "Df2", "f1", "f2"]
def frame_data():
    fin = DSet('/home/claude/s1spec/step3/outA5.txt', [0, 1, 2, 3, 4])
    C8h, S8h, Ct2h, St2h = fin.hull[5], fin.hull[6], fin.hull[7], fin.hull[8]
    CtB = C8h * Ct2h - S8h * St2h; StB = S8h * Ct2h + C8h * St2h
    return CtB, StB, Ct2h, St2h
def boxes():
    """yield per-box dict with ball enclosures of the 8 functions, weight dX*rho, and diagnostics"""
    CtB, StB, Ct2h, St2h = frame_data()
    T0 = arb(G3.T0); CH = chunks()
    for line in open('/home/claude/s1spec/step3/outA5.txt'):
        tag = line[:2]
        if tag not in ("MB", "MC"): continue
        v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
        h = [A(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
        if tag == "MB":
            st, C_, S_ = h[0], h[1], h[2]; Ct, St = CtB, StB
            c = min(int(ta + 1e-9), len(CH) - 1)
            kv, mv, pv_, dq = [A(lo, hi) for lo, hi in CH[c][2]]
            sig = A(32 - tb, 32 - ta); Xv = sig - st; qv = sig + dq; dX = arb(tb) - arb(ta)
            if tb <= 3.0 + 1e-12:
                uu = A(ta / 3, tb / 3); chi = [centered(dSm_pts, uu, j) * _pw(arb(-1) / 3, j) for j in range(5)]
            else: chi = [arb(1)] + [arb(0)] * 4
        else:
            pv_, qv, kv, mv, st, C_, S_ = h[0], h[1], h[2], h[3], h[4], h[7], h[8]; Ct, St = Ct2h, St2h
            cc = (8 - st) / T0
            Xv = hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta))); dX = cc * (arb(tb) - arb(ta))
            chi = [arb(1)] + [arb(0)] * 4
        lo = max(float(Xv.lower()), 0.0); hi = float(Xv.upper())
        Xv = arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15) + 1e-300)       # x >= 0 (domain of the half line)
        P, rho = weight(Xv)
        pvs = (kv, mv, qv, pv_, C_, S_)
        ps1 = [dp[n](*pvs) for n in range(5)]
        ps2 = [Ct * dC[n](*pvs) + St * dS[n](*pvs) for n in range(5)]
        a = {n: cS[n](P[0], P[1], P[2], P[3], kv, mv, qv, pv_) for n in range(5)}
        out = {}
        for key, ps in (("1", ps1), ("2", ps2)):
            g = [leib(chi, ps, n) for n in range(5)]
            out["Lf" + key] = sum((a[n] * g[n] for n in range(5)), arb(0))
            out["Df" + key] = g[4] - 2 * P[0] * g[3] + (P[0] * P[0] - P[1]) * g[2]
            out["f" + key] = g[0]
        # psi3^ jets (value, first derivative) on the ball Xv
        xc = float(Xv.mid()); y = HR.jets_any(F, Xv, xc)[:2]
        cth = (St * C_ - Ct * S_) if CSIGN is None else CSIGN * (Ct * S_ - St * C_)
        out["gA"] = chi[0] * (y[0] - qv * y[1]); out["gB"] = chi[0] * y[1] * cth
        yield dict(tag=tag, ta=ta, tb=tb, X=Xv, w=dX * rho, u=out, psi2=ps2[0], cth=cth, p=pv_, chi=chi[0])
if __name__ == "__main__":
    t0 = time.time(); n = len(NAMES)
    G = [[arb(0) for _ in range(n)] for _ in range(n)]; nb = 0
    for b in boxes():
        u = [b['u'][nm] for nm in NAMES]
        for i in range(n):
            for j in range(i, n): G[i][j] += b['w'] * u[i] * u[j]
        nb += 1
    for i in range(n):
        for j in range(i):  G[i][j] = G[j][i]
    for i in range(n):
        for j in range(n): assert G[i][j].is_finite(), (i, j)
    print("boxes", nb, f"[{time.time()-t0:.0f}s]")
    json.dump({"names": NAMES, "G": [[[float(G[i][j].lower()), float(G[i][j].upper())] for j in range(n)] for i in range(n)], "boxes": nb},
              open('gram8_raw.json', 'w'), indent=1)
    for i in range(n): print(NAMES[i], "  ".join(f"{float(G[i][j].mid()):+.5g}(r{float(G[i][j].rad()):.1e})" for j in range(n)))
