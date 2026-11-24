"""Rigorous ||(L - 2 lamt) h2^ + Q2(psi3^)||_{L^2(rho; [0,oo))} for the C^3 functions h2^ (h2fun.py) and psi3^ (psi3fun.py),
with the TIGHT step-1 profile.  Q2 = quadratic part of the graph equation (12-term polynomial, q2_derivs.pkl):
  Q2 = -4k^3 u u'' - 3k^3 u'^2 - 26k^2 m u u' - (5/2)k^2 p u^2 - 15 k m^2 u^2 - 4k u u'''' - 6k u'u''' - 3k u''^2
       - 6m u u''' - 8m u'u'' - p u u'' - (1/2) p u'^2        (vanishes on the ray k = m = p = 0).
Stage C: second-order Taylor form per CAPD box (as psi3res3.py), Q2 and its total x-derivatives from the ODE.
Stage B (sigma in [8, 32]): ray + chunk bounds for the h2 part (as psi3res2/3), |Q2| by interval evaluation with the chunk
bounds of (k, m, p).  Stage T (sigma >= 32): explicit extension bounds + |Q2| tail bound.
h2^ and psi3^ share their breakpoints, so the residual is smooth on every element."""
import sys, json, time, pickle
from flint import arb, arb_poly
import sympy as sp
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import psi3fun as F
import h2fun as Hf
from gen_hs2 import chunks
import hstail as HS
lam2 = Hf.LAM; STIP = Hf.STIP
assert abs(float((F.STIP - Hf.STIP).mid())) < 1e-15
TT = arb((10.249077650575472 + 10.249077651519443) / 2, (10.249077651519443 - 10.249077650575472) / 2 * 1.000001)
TT = TT.union(arb(10.249077650575472)).union(arb(10.249077651519443))
ST = 8 - TT; dST = abs(ST - STIP)
STRIDE = int(sys.argv[1]) if (len(sys.argv) > 1 and sys.argv[1].isdigit() and __name__ == '__main__') else 1
def rho_of(xlo): return (-(1 + arb(xlo) ** 2) ** (arb(2) / 3) / 4).exp()
# ---- Hermite corrections as arb polynomials (both modules)
NJ = 7
for MOD in (F, Hf):
    HBP = {(s, d): arb_poly([arb(c.p) / arb(c.q) for c in MOD.HB[(s, d)]]) for s in ('L', 'R') for d in range(4)}
    for P in MOD.PIECES:
        L = P.b - P.a; Hp = arb_poly([arb(0)])
        for d in range(4):
            if P.Hl[d] != 0: Hp = Hp + HBP[('L', d)] * (P.Hl[d] * L ** d)
            if P.Hr[d] != 0: Hp = Hp + HBP[('R', d)] * (P.Hr[d] * L ** d)
        P.Hd = [Hp]
        for kk_ in range(1, NJ + 1): P.Hd.append(P.Hd[-1].derivative())
        P.Lx = L
def jets(MOD, P, X, n=NJ):
    yy = P.y(X); hv = (X - P.a) / P.Lx
    return [MOD.clen(P.cfs[d], yy) * P.sc ** d + P.Hd[d](hv) / P.Lx ** d for d in range(n)]
# ---- Q2 and its total derivatives in arb
Q = pickle.load(open('/home/claude/s1spec/step3nl/q2_derivs.pkl', 'rb'))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
ARGS = "u0,u1,u2,u3,u4,u5,u6,pp,qq,kk,mm"
QF = {key: eval(f"lambda {ARGS}: " + _Pr().doprint(Q[key]), {"arb": arb}) for key in ('Q2', 'DQ2', 'D2Q2')}
def q2(key, y, p, q, k, m):
    yy = list(y) + [arb(0)] * (7 - len(y))
    return QF[key](*yy[:7], p, q, k, m)
# ------------------------------------------------------------------ stage C
def coef0(p, q, k, m): return k * k, 5 * k * m + q / 4, arb(1) / 4 - lam2 - 3 * m * m - k * p
def coef1(p, q, k, m): return 2 * k * m, 5 * m * m + arb(3) / 2 * k * p + arb(1) / 4, -arb(5) / 2 * m * p + k * k * q
def coef2(p, q, k, m):
    return (2 * m * m + k * p / 2, 4 * m * p - arb(3) / 2 * k * k * q,
            -arb(5) / 8 * p * p + arb(9) / 2 * k * m * q + k * k + k * k * k * p)
def I(lo, hi): return arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15)).union(arb(lo)).union(arb(hi))
def stage_C(path='/home/claude/s1_tight/prof/prof_tight.txt'):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0); parts = [0.0, 0.0, 0.0]; nb = 0
    for li, line in enumerate(open(path)):
        if line.startswith('END'): break
        if li % STRIDE: continue
        v = list(map(float, line.split())); ta, tb, tc = v[0], v[1], v[2]
        c = [I(v[3 + 2 * j], v[4 + 2 * j]) for j in range(4)]; hl = [I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(4)]
        Xc = TT - arb(tc); Hh = arb(tc - tb).union(arb(tc - ta)); Xb = Xc + Hh
        if Xb.upper() < 0: continue
        xlo = max(float(Xb.lower()), 0.0); xhi = float(Xb.upper())
        p, q, k, m = c; P_, Q_, K_, M_ = hl
        a2, a1, a0 = coef0(p, q, k, m); a2p, a1p, a0p = coef1(p, q, k, m)
        b2, b1, b0 = coef0(P_, Q_, K_, M_); b2p, b1p, b0p = coef1(P_, Q_, K_, M_); b2pp, b1pp, b0pp = coef2(P_, Q_, K_, M_)
        best = None
        for e, (Pf, Ph) in enumerate(zip(F.PIECES, Hf.PIECES)):
            if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
            h = jets(Hf, Ph, Xc); Hb = jets(Hf, Ph, Xb); y = jets(F, Pf, Xc); Yb = jets(F, Pf, Xb)
            R0 = -h[4] - a2 * h[2] - a1 * h[1] + a0 * h[0] + q2('Q2', y[:5], p, q, k, m)
            R1 = (-h[5] - a2p * h[2] - a2 * h[3] - a1p * h[1] - a1 * h[2] + a0p * h[0] + a0 * h[1]
                  + q2('DQ2', y[:6], p, q, k, m))
            R2 = (-Hb[6] - b2pp * Hb[2] - 2 * b2p * Hb[3] - b2 * Hb[4] - b1pp * Hb[1] - 2 * b1p * Hb[2] - b1 * Hb[3]
                  + b0pp * Hb[0] + 2 * b0p * Hb[1] + b0 * Hb[2] + q2('D2Q2', Yb[:7], P_, Q_, K_, M_))
            r = abs(R0 + R1 * Hh + R2 * Hh * Hh / 2)
            best = r if best is None else arb(max(best.upper(), r.upper()))
            hm = float(abs(Hh).upper())
            parts[0] = max(parts[0], float(abs(R0).upper())); parts[1] = max(parts[1], float(abs(R1).upper()) * hm)
            parts[2] = max(parts[2], float(abs(R2).upper()) * hm * hm / 2)
        ru = arb(best.upper()); rho = rho_of(xlo)
        tot += arb(xhi - xlo) * rho * ru * ru; nb += 1
        w = float(ru.upper()) * float(rho.sqrt().upper())
        if w > worst[0]: worst = (w, xlo)
    print(f"C: {nb} boxes, ||r||^2 <= {float(tot.upper()):.4e} (||r|| <= {float(tot.upper())**.5:.3e}); max weighted |r| "
          f"{worst[0]:.3e} at x = {worst[1]:.4f}; max |R(xc)| {parts[0]:.2e}, |R'|h {parts[1]:.2e}, |R''|h^2/2 {parts[2]:.2e} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return tot, worst, parts

# ------------------------------------------------------------------ stage Bt: sigma in [8, SB_END], tight tail profile
SB_END = 20
sys.path.insert(0, '/home/claude/s1_tight')
import tail as T1mod
from tail2 import nl_radii, yhat
from laplace import K as Klap
_BR = json.load(open('/home/claude/s1_tight/bracket2.json'))
EPS = arb(_BR['eps_a']).union(arb(_BR['eps_b']))
_R = nl_radii(arb(8), arb(1.62e-4 * (1 + 1e-12))); assert _R['verified']
RN = _R['rn0'] + 2 * _R['rn1']; RTOT = _R['R']
_Y, _Z = yhat(arb(8)); CN = (_Y[0] / Klap(arb(8), 0)).real
def Efun(s): return T1mod.A_E * (arb(s) ** (arb(4) / 3) - arb(8) ** (arb(4) / 3))
def tail_state(sig):
    """enclosure of (p, q, k, m) of S1 at the point sigma = sig (arb point), from the tight tail theorem"""
    mu = (sig / 4) ** (arb(1) / 3); w = Efun(sig).exp() * (sig / 8) ** (arb(1) / 3)
    Kj = [Klap(sig, j) for j in range(3)]
    k = EPS * CN * Kj[0] + arb(0, (RN / w).upper())
    m = EPS * CN * Kj[1] + arb(0, (mu * RN / w).upper())
    p = 4 * EPS * CN * Kj[2] + arb(0, (4 * mu * mu * RN / w).upper())
    J = RTOT ** 2 * (arb(8) / 4) ** (arb(2) / 3) * (-2 * Efun(sig)).exp() / (2 * mu)
    q = sig + 4 * k * m + 4 * J * arb(0.5, 0.5)       # J in [0, Jb]
    return p, q, k, m
def stage_Bt(N=200):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0); parts = [0.0, 0.0, 0.0]
    for ci in range(8, SB_END):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N; sc = (sa + sb) / 2
            p, q, k, m = tail_state(sc)
            kb, mb, pb, qs = HS.profile_bounds(sa)
            P_, K_, M_ = arb(0, pb.upper()), arb(0, kb.upper()), arb(0, mb.upper())
            Q_ = sa.union(sb) + arb(0, qs.upper())
            Xc = sc - ST; Hh = (sa - sc).union(sb - sc); Xb = Xc + Hh
            xlo = float(Xb.lower()); xhi = float(Xb.upper())
            a2, a1, a0 = coef0(p, q, k, m); a2p, a1p, a0p = coef1(p, q, k, m)
            b2, b1, b0 = coef0(P_, Q_, K_, M_); b2p, b1p, b0p = coef1(P_, Q_, K_, M_); b2pp, b1pp, b0pp = coef2(P_, Q_, K_, M_)
            best = None
            for e, (Pf, Ph) in enumerate(zip(F.PIECES, Hf.PIECES)):
                if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
                h = jets(Hf, Ph, Xc); Hb = jets(Hf, Ph, Xb); y = jets(F, Pf, Xc); Yb = jets(F, Pf, Xb)
                R0 = -h[4] - a2 * h[2] - a1 * h[1] + a0 * h[0] + q2('Q2', y[:5], p, q, k, m)
                R1 = (-h[5] - a2p * h[2] - a2 * h[3] - a1p * h[1] - a1 * h[2] + a0p * h[0] + a0 * h[1]
                      + q2('DQ2', y[:6], p, q, k, m))
                R2 = (-Hb[6] - b2pp * Hb[2] - 2 * b2p * Hb[3] - b2 * Hb[4] - b1pp * Hb[1] - 2 * b1p * Hb[2] - b1 * Hb[3]
                      + b0pp * Hb[0] + 2 * b0p * Hb[1] + b0 * Hb[2] + q2('D2Q2', Yb[:7], P_, Q_, K_, M_))
                r = abs(R0 + R1 * Hh + R2 * Hh * Hh / 2)
                best = r if best is None else arb(max(best.upper(), r.upper()))
                hm = float(abs(Hh).upper())
                parts[0] = max(parts[0], float(abs(R0).upper()) * float(rho_of(xlo).sqrt().upper()))
                parts[1] = max(parts[1], float(abs(R1).upper()) * hm * float(rho_of(xlo).sqrt().upper()))
                parts[2] = max(parts[2], float(abs(R2).upper()) * hm * hm / 2 * float(rho_of(xlo).sqrt().upper()))
            ru = arb(best.upper()); rho = rho_of(xlo)
            tot += arb(xhi - xlo) * rho * ru * ru
            wv = float(ru.upper()) * float(rho.sqrt().upper())
            if wv > worst[0]: worst = (wv, xlo)
    print(f"Bt (sigma in [8,{SB_END}], tight tail): ||r||^2 <= {float(tot.upper()):.4e} (||r|| <= {float(tot.upper())**.5:.3e}); "
          f"max weighted |r| {worst[0]:.3e} at x = {worst[1]:.3f}; max weighted R0 {parts[0]:.2e}, R'h {parts[1]:.2e}, R''h^2/2 {parts[2]:.2e} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return tot, worst, parts
# ------------------------------------------------------------------ stage B (ray + chunk bounds) on [SB_END, 32]
def u_box(MOD, X, xc, n=5):
    out = None
    for P in MOD.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        jc = MOD.jet_poly(P, arb(xc), n + 2); jX = MOD.jet_poly(P, X, n + 2); hh = X - arb(xc); Hc = MOD.Hcorr(P, X, n)
        vals = [jc[d] + jc[d + 1] * hh + jX[d + 2] * hh * hh / 2 + Hc[d] for d in range(n)]
        out = vals if out is None else [o.union(v) for o, v in zip(out, vals)]
    return out
def r_ray_poly(MOD, lam, X, xc):
    best = None
    for P in MOD.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        c = arb(xc); jc = MOD.jet_poly(P, c, 7); jX = MOD.jet_poly(P, X, 7); hh = X - c
        r0 = -jc[4] - (c + STIP) / 4 * jc[1] + (arb(1) / 4 - lam) * jc[0]
        r1 = -jc[5] - jc[1] / 4 - (c + STIP) / 4 * jc[2] + (arb(1) / 4 - lam) * jc[1]
        r2 = -jX[6] - jX[2] / 2 - (X + STIP) / 4 * jX[3] + (arb(1) / 4 - lam) * jX[2]
        Hc = MOD.Hcorr(P, X, 5)
        rH = -Hc[4] - (X + STIP) / 4 * Hc[1] + (arb(1) / 4 - lam) * Hc[0]
        v = abs(r0 + r1 * hh + r2 * hh * hh / 2) + abs(rH)
        best = v if best is None else best.union(v)
    return best
Xs = Hf.XC + STIP
def r_ray_ext(MOD, lam, X):
    s = X + STIP
    z = -(s / 4) ** (arb(1) / 3); z1 = -(s / 4) ** (arb(-2) / 3) / 12; z2 = (s / 4) ** (arb(-5) / 3) / 72
    z3 = -5 * (s / 4) ** (arb(-8) / 3) / 864
    E = (-3 * (s / 4) ** (arb(4) / 3) + 3 * (Xs / 4) ** (arb(4) / 3)).exp()
    return abs(MOD.C1 * 4 * MOD.AJ[4] * s ** (MOD.A - 16)) + abs(MOD.C2 * E * (-(z3 + 4 * z * z2 + 3 * z1 * z1 + 6 * z * z * z1) + arb(1) / 4 - lam))
def jets_any(MOD, X, xc):
    if X.upper() <= float(MOD.XC.lower()): return u_box(MOD, X, xc)
    if X.lower() >= float(MOD.XC.upper()): return MOD.ext_jet(X, 5)
    return [a.union(b) for a, b in zip(u_box(MOD, X, xc), MOD.ext_jet(X, 5))]
def stage_B():
    t0 = time.time(); totB = arb(0); worstB = (0, 0); parts = {"ray": 0.0, "st": 0.0, "prof": 0.0, "Q2": 0.0}
    for ci, (ln, hstep, bnds) in enumerate(chunks()):
        s_hi = 32.0 - ci; s_lo = s_hi - 1
        if s_lo < SB_END: continue
        Bk, Bm, Bp, Bq = [arb(float(max(abs(lo), abs(hi)))) for lo, hi in bnds]
        N = 100
        for i in range(N):
            sa, sb = s_lo + i / N, s_lo + (i + 1) / N
            xa = float((arb(sa) - ST).lower()); xb = float((arb(sb) - ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2)
            u = jets_any(Hf, X, xc)
            if xb <= float(Hf.XC.lower()): rr = r_ray_poly(Hf, lam2, X, xc)
            elif xa >= float(Hf.XC.upper()): rr = r_ray_ext(Hf, lam2, X)
            else: rr = r_ray_poly(Hf, lam2, X, xc).union(r_ray_ext(Hf, lam2, X))
            y = jets_any(F, X, xc)
            q2v = abs(q2('Q2', y, arb(0, Bp.upper()), X + STIP, arb(0, Bk.upper()), arb(0, Bm.upper())))
            rs = dST / 4 * abs(u[1])
            dl = Bk * Bk * abs(u[2]) + (5 * Bk * Bm + Bq / 4) * abs(u[1]) + (3 * Bm * Bm + Bk * Bp) * abs(u[0])
            rb = float((rr + rs + dl + q2v).upper()); rho = rho_of(xa)
            totB += arb(xb - xa) * rho * arb(rb) ** 2
            sr = float(rho.sqrt().upper())
            for key, val in (("ray", rr), ("st", rs), ("prof", dl), ("Q2", q2v)): parts[key] = max(parts[key], float(val.upper()) * sr)
            if rb * sr > worstB[0]: worstB = (rb * sr, xa)
    print(f"B: ||r||^2 <= {float(totB.upper()):.4e} (||r|| <= {float(totB.upper())**.5:.3e}), max weighted {worstB[0]:.2e} at "
          f"x={worstB[1]:.2f}; max weighted parts {parts} [{time.time()-t0:.0f}s]", flush=True)
    return totB, worstB, parts
# ------------------------------------------------------------------ stage T
def Mbound(MOD, d, s1):
    v = arb(0)
    for j in range(4):
        e = MOD.A - 4 * j; f = arb(1)
        for i in range(d): f = f * (e - i)
        v += abs(MOD.AJ[j] * f) * s1 ** (-4 * j - d)
    zmax = 1 + (s1 / 4) ** (arb(1) / 3)
    return abs(MOD.C1) * v + abs(MOD.C2) * zmax ** d * 4
def stage_T():
    x1 = arb(32) - ST.lower(); s1 = x1 + STIP
    kb, mb, pb, qs = HS.profile_bounds(32)
    M = [Mbound(Hf, d, s1) for d in range(3)]
    rext = abs(Hf.C1 * 4 * Hf.AJ[4]) * s1 ** (Hf.A - 16) / s1 ** Hf.A + abs(Hf.C2) * (arb(1) / 2 + 1) / s1 ** Hf.A
    Mr = rext + dST / 4 * M[1] + kb * kb * M[2] + (5 * kb * mb + qs / 4) * M[1] + (3 * mb * mb + kb * pb) * M[0]
    # |Q2(psi3^)| <= sum|coef| * (max_d Mpsi_d)^2 * s^(2 A_psi), coefficients bounded with kb, mb, pb (k <= 1)
    Mp = max(float(Mbound(F, d, s1).upper()) for d in range(5)) * (1 + 1e-12)
    Mp = arb(Mp)   # crude: |psi^(d)| <= Mp s^A_psi for d <= 4 (derivatives only decrease the power)
    cq = 4 * kb ** 3 + 3 * kb ** 3 + 26 * kb * kb * mb + arb(5) / 2 * kb * kb * pb + 15 * kb * mb * mb + 4 * kb + 6 * kb + 3 * kb + 6 * mb + 8 * mb + pb + pb / 2
    Mq = cq * Mp * Mp                                  # |Q2| <= Mq s^(2 A_psi)
    B = 2 * F.A if 2 * F.A > Hf.A else Hf.A            # common exponent (s >= 1)
    Mtot = Mr + Mq
    totT = Mtot * Mtot * x1 ** (2 * B) * (-(x1 ** (arb(4) / 3)) / 4).exp() / (x1 ** (arb(1) / 3) / 3 - 2 * B / x1)
    print(f"T: ||r||^2 <= {float(totT.upper()):.3e}")
    return totT
def h2_norm():
    """||h2^||_rho on the half line (stage-C/B grids with Taylor-form values; tail by the M-bound)."""
    tot = arb(0); N = 4000
    for i in range(N):
        xa, xb = 34.25 * i / N, 34.25 * (i + 1) / N; xc = (xa + xb) / 2; X = arb(xc, (xb - xa) / 2)
        u0 = jets_any(Hf, X, xc)[0]
        tot += arb(xb - xa) * rho_of(xa) * abs(u0) * abs(u0)
    x1 = arb(34.25); s1 = x1 + STIP; M0 = Mbound(Hf, 0, s1)
    tot += M0 * M0 * x1 ** (2 * Hf.A) * (-(x1 ** (arb(4) / 3)) / 4).exp() / (x1 ** (arb(1) / 3) / 3 - 2 * Hf.A / x1)
    return tot
if __name__ == "__main__":
    print("tip interval ST =", ST, " |ST - STIP| <=", float(dST.upper()))
    totC, worstC, partsC = stage_C()
    if STRIDE > 1:
        print("(stride > 1: timing/diagnostic run only)"); sys.exit(0)
    totBt, worstBt, partsBt = stage_Bt()
    totB, worstB, partsB = stage_B(); totT = stage_T()
    tot = totC + totBt + totB + totT
    n2 = h2_norm()
    res = {"C2": float(totC.upper()), "Bt2": float(totBt.upper()), "worstBt": worstBt, "B2": float(totB.upper()), "T2": float(totT.upper()), "r_half": float(tot.upper()) ** 0.5,
           "worstC": worstC, "partsC": partsC, "worstB": worstB, "partsB": partsB, "h2_norm2_half": [float(n2.lower()), float(n2.upper())],
           "h2_0": float(jets(Hf, Hf.PIECES[0], arb(0))[0].mid())}
    json.dump(res, open('h2res.json', 'w'), indent=1)
    print(f"TOTAL (half line): ||(L - 2 lamt) h2^ + Q2(psi3^)||_rho <= {res['r_half']:.4e};  ||h2^||_rho^2 in {res['h2_norm2_half']}")
