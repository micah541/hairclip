"""Rigorous ||(L - lamt) psi3^||_{L^2(rho; [0,oo))} with the TIGHT step-1 profile (s1_tight).
Stage C (x in [0, tau_tip], i.e. sigma in [tip, 8]): second-order Taylor form of the residual itself on each of the
  51,246 CAPD boxes of s1_tight/prof/prof_tight.txt (centre enclosure + hull per box):
     R(x) in R(xc) + R'(xc) h + R''(box) h^2/2,   h = x - xc in [tc - tb, tc - ta],  xc = tau_tip - tc,
  with the profile derivatives from the ODE (k' = m, m' = p/4, p' = -k q, q' = 1 + k p, ' = d/dx = d/dsigma).
  R = -u4 - a2 u2 - a1 u1 + a0 u0,  a2 = k^2, a1 = 5 k m + q/4, a0 = 1/4 - lamt - 3 m^2 - k p.
Stage B (sigma in [8, 32]) and T (sigma >= 32): as psi3res2.py, with the tight tip interval.
Usage: python3 psi3res3.py [stageC_stride]   (stride > 1 only for timing tests: then C is NOT a bound)"""
import sys, json, time
from flint import arb, arb_poly
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import psi3fun as F
from gen_hs3 import T0
from gen_hs2 import chunks
import hstail as H
lam = F.LAM; STIP = F.STIP
TT = arb((10.249077650575472 + 10.249077651519443) / 2, (10.249077651519443 - 10.249077650575472) / 2 * 1.000001)
TT = TT.union(arb(10.249077650575472)).union(arb(10.249077651519443))     # tau_tip enclosure (s1curve TIP, tight bracket)
ST = 8 - TT                                      # sigma_tip enclosure
dST = abs(ST - STIP)
STRIDE = int(sys.argv[1]) if len(sys.argv) > 1 else 1
def rho_of(xlo): return (-(1 + arb(xlo) ** 2) ** (arb(2) / 3) / 4).exp()
# ---- Hermite corrections as arb polynomials in h (per piece), derivatives precomputed
def hb_poly(side, d):
    cf = F.HB[(side, d)]
    return arb_poly([arb(c.p) / arb(c.q) for c in cf])
HBP = {(s, d): hb_poly(s, d) for s in ('L', 'R') for d in range(4)}
NJ = 7
for P in F.PIECES:
    L = P.b - P.a
    Hp = arb_poly([arb(0)])
    for d in range(4):
        if P.Hl[d] != 0: Hp = Hp + HBP[('L', d)] * (P.Hl[d] * L ** d)
        if P.Hr[d] != 0: Hp = Hp + HBP[('R', d)] * (P.Hr[d] * L ** d)
    P.Hd = [Hp]
    for k in range(1, NJ + 1): P.Hd.append(P.Hd[-1].derivative())
    P.Lx = L
def jets(P, X, n=NJ):
    """u^(d)(X), d < n, for the C^3 piece P (polynomial + Hermite correction), X arb"""
    yy = P.y(X); hv = (X - P.a) / P.Lx
    return [F.clen(P.cfs[d], yy) * P.sc ** d + P.Hd[d](hv) / P.Lx ** d for d in range(n)]
# ------------------------------------------------------------------ stage C
def coef0(p, q, k, m):
    a2 = k * k; a1 = 5 * k * m + q / 4; a0 = arb(1) / 4 - lam - 3 * m * m - k * p
    return a2, a1, a0
def coef1(p, q, k, m):
    a2p = 2 * k * m; a1p = 5 * m * m + arb(3) / 2 * k * p + arb(1) / 4; a0p = -arb(5) / 2 * m * p + k * k * q
    return a2p, a1p, a0p
def coef2(p, q, k, m):
    a2pp = 2 * m * m + k * p / 2; a1pp = 4 * m * p - arb(3) / 2 * k * k * q
    a0pp = -arb(5) / 8 * p * p + arb(9) / 2 * k * m * q + k * k + k * k * k * p
    return a2pp, a1pp, a0pp
def I(lo, hi): return arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15)) .union(arb(lo)).union(arb(hi))
def stage_C(path='/home/claude/s1_tight/prof/prof_tight.txt'):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0); parts = [0.0, 0.0, 0.0]; nb = 0
    for li, line in enumerate(open(path)):
        if line.startswith('END'): break
        if li % STRIDE: continue
        v = list(map(float, line.split()))
        ta, tb, tc = v[0], v[1], v[2]
        c = [I(v[3 + 2 * j], v[4 + 2 * j]) for j in range(4)]       # p q k m at tc
        hl = [I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(4)]    # hulls
        Xc = TT - arb(tc)
        Hh = arb(tc - tb).union(arb(tc - ta))                         # h = x - xc (exact doubles)
        Xb = Xc + Hh
        if Xb.upper() < 0: continue
        xlo = max(float(Xb.lower()), 0.0); xhi = float(Xb.upper())
        p, q, k, m = c
        a2, a1, a0 = coef0(p, q, k, m); a2p, a1p, a0p = coef1(p, q, k, m)
        P_, Q_, K_, M_ = hl
        b2, b1, b0 = coef0(P_, Q_, K_, M_); b2p, b1p, b0p = coef1(P_, Q_, K_, M_); b2pp, b1pp, b0pp = coef2(P_, Q_, K_, M_)
        best = None
        for Pc in F.PIECES:
            if Xb.upper() < Pc.a.lower() or Xb.lower() > Pc.b.upper(): continue
            u = jets(Pc, Xc); U = jets(Pc, Xb)
            R0 = -u[4] - a2 * u[2] - a1 * u[1] + a0 * u[0]
            R1 = -u[5] - a2p * u[2] - a2 * u[3] - a1p * u[1] - a1 * u[2] + a0p * u[0] + a0 * u[1]
            R2 = (-U[6] - b2pp * U[2] - 2 * b2p * U[3] - b2 * U[4] - b1pp * U[1] - 2 * b1p * U[2] - b1 * U[3]
                  + b0pp * U[0] + 2 * b0p * U[1] + b0 * U[2])
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
# ------------------------------------------------------------------ stage B (as psi3res2, tight ST)
def u_box(X, xc, n=5):
    out = None
    for P in F.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        jc = F.jet_poly(P, arb(xc), n + 2); jX = F.jet_poly(P, X, n + 2); hh = X - arb(xc); Hc = F.Hcorr(P, X, n)
        vals = [jc[d] + jc[d + 1] * hh + jX[d + 2] * hh * hh / 2 + Hc[d] for d in range(n)]
        out = vals if out is None else [o.union(v) for o, v in zip(out, vals)]
    return out
def r_ray_poly(X, xc):
    best = None
    for P in F.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        c = arb(xc); jc = F.jet_poly(P, c, 7); jX = F.jet_poly(P, X, 7); hh = X - c
        r0 = -jc[4] - (c + STIP) / 4 * jc[1] + (arb(1) / 4 - lam) * jc[0]
        r1 = -jc[5] - jc[1] / 4 - (c + STIP) / 4 * jc[2] + (arb(1) / 4 - lam) * jc[1]
        r2 = -jX[6] - jX[2] / 2 - (X + STIP) / 4 * jX[3] + (arb(1) / 4 - lam) * jX[2]
        Hc = F.Hcorr(P, X, 5)
        rH = -Hc[4] - (X + STIP) / 4 * Hc[1] + (arb(1) / 4 - lam) * Hc[0]
        v = abs(r0 + r1 * hh + r2 * hh * hh / 2) + abs(rH)
        best = v if best is None else best.union(v)
    return best
Xs = F.XC + STIP
def r_ray_ext(X):
    s = X + STIP
    z = -(s / 4) ** (arb(1) / 3); z1 = -(s / 4) ** (arb(-2) / 3) / 12; z2 = (s / 4) ** (arb(-5) / 3) / 72
    z3 = -5 * (s / 4) ** (arb(-8) / 3) / 864
    E = (-3 * (s / 4) ** (arb(4) / 3) + 3 * (Xs / 4) ** (arb(4) / 3)).exp()
    return abs(F.C1 * 4 * F.AJ[4] * s ** (F.A - 16)) + abs(F.C2 * E * (-(z3 + 4 * z * z2 + 3 * z1 * z1 + 6 * z * z * z1) + arb(1) / 4 - lam))
def stage_B():
    t0 = time.time(); totB = arb(0); worstB = (0, 0); parts = {"ray": 0.0, "st": 0.0, "prof": 0.0}
    for ci, (ln, hstep, bnds) in enumerate(chunks()):
        s_hi = 32.0 - ci; s_lo = s_hi - 1
        Bk, Bm, Bp, Bq = [arb(float(max(abs(lo), abs(hi)))) for lo, hi in bnds]
        N = 100
        for i in range(N):
            sa, sb = s_lo + i / N, s_lo + (i + 1) / N
            xa = float((arb(sa) - ST).lower()); xb = float((arb(sb) - ST).upper()); xc = (xa + xb) / 2
            X = arb(xc, (xb - xa) / 2)
            if xb <= float(F.XC.lower()):
                u = u_box(X, xc); rr = r_ray_poly(X, xc)
            elif xa >= float(F.XC.upper()):
                u = F.ext_jet(X, 5); rr = r_ray_ext(X)
            else:
                u = [a.union(b) for a, b in zip(u_box(X, xc), F.ext_jet(X, 5))]; rr = r_ray_poly(X, xc).union(r_ray_ext(X))
            rs = dST / 4 * abs(u[1])
            dl = Bk * Bk * abs(u[2]) + (5 * Bk * Bm + Bq / 4) * abs(u[1]) + (3 * Bm * Bm + Bk * Bp) * abs(u[0])
            rb = float((rr + rs + dl).upper()); rho = rho_of(xa)
            totB += arb(xb - xa) * rho * arb(rb) ** 2
            sr = float(rho.sqrt().upper())
            parts["ray"] = max(parts["ray"], float(rr.upper()) * sr); parts["st"] = max(parts["st"], float(rs.upper()) * sr)
            parts["prof"] = max(parts["prof"], float(dl.upper()) * sr)
            if rb * sr > worstB[0]: worstB = (rb * sr, xa)
    print(f"B: ||r||^2 <= {float(totB.upper()):.4e} (||r|| <= {float(totB.upper())**.5:.3e}), max weighted {worstB[0]:.2e} at "
          f"x={worstB[1]:.2f}; max weighted parts {parts} [{time.time()-t0:.0f}s]", flush=True)
    return totB, worstB, parts
def stage_T():
    x1 = arb(32) - ST.lower(); s1 = x1 + STIP
    kb, mb, pb, qs = H.profile_bounds(32)
    def Mbound(d):
        v = arb(0)
        for j in range(4):
            e = F.A - 4 * j; f = arb(1)
            for i in range(d): f = f * (e - i)
            v += abs(F.AJ[j] * f) * s1 ** (-4 * j - d)
        zmax = 1 + (s1 / 4) ** (arb(1) / 3)
        return abs(F.C1) * v + abs(F.C2) * zmax ** d * 4
    M0, M1, M2 = Mbound(0), Mbound(1), Mbound(2)
    rext = abs(F.C1 * 4 * F.AJ[4]) * s1 ** (F.A - 16) / s1 ** F.A + abs(F.C2) * (arb(1) / 2 + 1) / s1 ** F.A
    Mr = rext + dST / 4 * M1 + kb * kb * M2 + (5 * kb * mb + qs / 4) * M1 + (3 * mb * mb + kb * pb) * M0
    totT = Mr * Mr * x1 ** (2 * F.A) * (-(x1 ** (arb(4) / 3)) / 4).exp() / (x1 ** (arb(1) / 3) / 3 - 2 * F.A / x1)
    print(f"T: ||r||^2 <= {float(totT.upper()):.3e}")
    return totT
if __name__ == "__main__":
    print("tip interval ST =", ST, " |ST - STIP| <=", float(dST.upper()))
    totC, worstC, partsC = stage_C()
    if STRIDE > 1:
        print("(stride > 1: timing/diagnostic run only)"); sys.exit(0)
    totB, worstB, partsB = stage_B()
    totT = stage_T()
    tot = totC + totB + totT
    res = {"C2": float(totC.upper()), "B2": float(totB.upper()), "T2": float(totT.upper()), "r_half": float(tot.upper()) ** 0.5,
           "worstC": worstC, "partsC": partsC, "worstB": worstB, "partsB": partsB, "ST": [float(ST.lower()), float(ST.upper())]}
    json.dump(res, open('psi3res3.json', 'w'), indent=1)
    print(f"TOTAL (half line): ||(L - lamt) psi3^||_rho <= {res['r_half']:.4e}   (full line x sqrt2: {res['r_half']*2**.5:.4e})")
