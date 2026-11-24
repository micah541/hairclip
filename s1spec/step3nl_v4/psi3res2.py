"""Rigorous ||(L - lamt) psi3^||_{L^2(rho; [0,oo))} (half line, as in the linear certificate) for the C^3 function
psi3^ of psi3fun.py.  Three ranges:
 C: x in [0, 8 - st]   (CAPD profile hulls on 51,246 boxes, /tmp/mc_prof.txt);
 B: sigma in [8, 32]    (profile = ray + chunk bounds; ray part by a second-order Taylor form of the residual itself);
 T: sigma >= 32         (explicit extension; closed-form residual; integral to infinity by a convexity bound)."""
import sys, json, time
from flint import arb
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import psi3fun as F
from gen_hs3 import T0
from gen_hs2 import chunks
import hstail as H
lam = F.LAM; STIP = F.STIP
ST = arb((-2.249127 - 2.249019) / 2, (2.249127 - 2.249019) / 2)
dST = abs(ST - STIP)                       # |st_true - stip|
def rho_of(xlo): return (-(1 + arb(xlo) ** 2) ** (arb(2) / 3) / 4).exp()
def u_box(X, xc, n=5):
    """enclosure of u^(d)(X), d < n, union over the pieces meeting X (Taylor form for P, direct for H)"""
    out = None
    for P in F.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        jc = F.jet_poly(P, arb(xc), n + 2); jX = F.jet_poly(P, X, n + 2); hh = X - arb(xc); Hc = F.Hcorr(P, X, n)
        vals = [jc[d] + jc[d + 1] * hh + jX[d + 2] * hh * hh / 2 + Hc[d] for d in range(n)]
        out = vals if out is None else [o.union(v) for o, v in zip(out, vals)]
    return out
# ------------------------------------------------------------------ stage C
t0 = time.time(); totC = arb(0); worstC = (0, 0)
for line in open('/tmp/mc_prof.txt'):
    ta, tb, p0, p1, q0, q1, k0, k1, m0, m1, s0, s1 = map(float, line.split())
    xl = (8 - arb(s1)) * (1 - arb(tb) / arb(T0)); xh = (8 - arb(s0)) * (1 - arb(ta) / arb(T0))
    xlo, xhi = max(float(xl.lower()), 0.0), float(xh.upper())
    X = arb((xlo + xhi) / 2, (xhi - xlo) / 2); u = u_box(X, (xlo + xhi) / 2)
    I = lambda lo, hi: arb((lo + hi) / 2, (hi - lo) / 2)
    p, q, k, m = I(p0, p1), I(q0, q1), I(k0, k1), I(m0, m1)
    r = -u[4] - k * k * u[2] - (5 * k * m + q / 4) * u[1] + (arb(1) / 4 - 3 * m * m - k * p - lam) * u[0]
    ru = float(abs(r).upper()); rho = rho_of(xlo)
    totC += arb(xhi - xlo) * rho * arb(ru) ** 2
    w = ru * float(rho.sqrt().upper())
    if w > worstC[0]: worstC = (w, xlo)
print(f"C: ||r||^2 <= {float(totC.upper()):.4e}  (||r|| <= {float(totC.upper())**.5:.3e}), max weighted {worstC[0]:.2e} at x={worstC[1]:.3f}  [{time.time()-t0:.0f}s]", flush=True)
# ------------------------------------------------------------------ stage B
def r_ray_poly(X, xc):
    """bound of |R_stip[P + H](X)| on the pieces meeting X, R_s = -d^4 - ((x + stip)/4) d + 1/4 - lamt"""
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
print(f"B: ||r||^2 <= {float(totB.upper()):.4e}  (||r|| <= {float(totB.upper())**.5:.3e}), max weighted {worstB[0]:.2e} at x={worstB[1]:.2f}; "
      f"max weighted parts {parts}  [{time.time()-t0:.0f}s]", flush=True)
# ------------------------------------------------------------------ tail sigma >= 32, i.e. x >= x1 = 32 - st_lo
x1 = arb(32) - ST.lower(); s1 = x1 + STIP           # s1 >= 32
kb, mb, pb, qs = H.profile_bounds(32)
# |u^(d)| <= M_d s^A for s >= s1, d <= 2 (sum of |a_j| prod|A-4j-i| s^(-4j-d) <= value at s1, decreasing; |c2 E stuff| <= |c2|(|z|+1)^d with |z| <= (s/4)^(1/3) <= s^A/... )
def Mbound(d):
    v = arb(0)
    for j in range(4):
        e = F.A - 4 * j; f = arb(1)
        for i in range(d): f = f * (e - i)
        v += abs(F.AJ[j] * f) * s1 ** (-4 * j - d)
    zmax = 1 + (s1 / 4) ** (arb(1) / 3)              # (s/4)^(1/3)/s^A is decreasing; bound |c2 E (e^S)^(d)/e^S| <= |c2| zmax^d (s/s1)^{d/3}
    return abs(F.C1) * v + abs(F.C2) * zmax ** d * 4    # crude: zmax^d (s/s1)^(d/3) <= zmax^d s^A/s1^A * 4 for d<=2
M0, M1, M2 = Mbound(0), Mbound(1), Mbound(2)
rext = abs(F.C1 * 4 * F.AJ[4]) * s1 ** (F.A - 16) / s1 ** F.A + abs(F.C2) * (arb(1) / 2 + 1) / s1 ** F.A   # |R g|, |R E| over s^A, s >= s1
Mr = rext + dST / 4 * M1 + kb * kb * M2 + (5 * kb * mb + qs / 4) * M1 + (3 * mb * mb + kb * pb) * M0
# int_{x1}^oo rho r^2 <= Mr^2 int x^{2A} e^{-x^{4/3}/4} <= Mr^2 x1^{2A} e^{-x1^{4/3}/4} / (x1^{1/3}/3 - 2A/x1)
totT = Mr * Mr * x1 ** (2 * F.A) * (-(x1 ** (arb(4) / 3)) / 4).exp() / (x1 ** (arb(1) / 3) / 3 - 2 * F.A / x1)
print(f"T: ||r||^2 <= {float(totT.upper()):.3e}")
tot = totC + totB + totT
res = {"C2": float(totC.upper()), "B2": float(totB.upper()), "T2": float(totT.upper()), "r_half": float(tot.upper()) ** 0.5,
       "worstC": worstC, "worstB": worstB, "jumps": F.JUMPS}
json.dump(res, open('psi3res2.json', 'w'), indent=1)
print(f"TOTAL (half line): ||(L - lamt) psi3^||_rho <= {res['r_half']:.4e}   (full line x sqrt2: {res['r_half']*2**.5:.4e})")
