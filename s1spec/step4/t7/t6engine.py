"""T6 layer 3: per-element bivariate (s, tc) enclosures of the base and of its residual.

For element e (x = x_e + h (1+s)/2):
  exact in s (BP ex=True):  kappa, kappa_x .. kappa_x^(5), Q = int_0^x kappa_x^2, Phi, rho_kappa = kappa_tau|_x - Phi_x,
                            rho_theta = int_0^x rho_kappa (= theta_tau - Phi), theta = -pi/2 + int_0^x kappa
  pointwise (BP ex=False):  E = exp(i phi), phi = theta - theta(x_e);
                            z = r + i A = E [ z_e + (h/2) int_{-1}^s conj(E) rho_theta ]   (r: normal, A: tangential residual)
                            r' = rho_theta - kappa A, r'' = rho_kappa - kappa' A - kappa^2 r,
                            r''' = rho_kappa' - kappa'' A - 3 kappa kappa' r - kappa^2 r', S1 = int_0^x r, S2 = int_0^x S1
                            G = conj(E) [G_e + (h/2) int E] = q + i p (position in the tangent frame), e^{i theta} = P_e E
The identities (T1 resid.py / signs.py; re-derived in T6_REPORT): with V* := -k'' + c p + b n and w* := -k k' + Q + c(q-x) + b t1,
the true normal / tangential velocities of the base are V~ = V* + r and varpi = w* + A, and (r, A)' = (rho_theta - k A, k r)."""
import math
from flint import arb, arb_poly
import tp6 as T
from tp6 import BP, CBP
from t6base import expand

HALF = arb(1) / 2
TWO = arb(2)
MONE = -arb(1)


def tp_scal(c):
    return T.tconst(c)


def xpoly(xe, h):
    """x(s, tc) = x_e + (h/2)(1+s) as a BP"""
    a = (xe + h * HALF).coeffs(); b = (h * HALF).coeffs()
    n = max(len(a), len(b)); a = a + [arb(0)] * (n - len(a)); b = b + [arb(0)] * (n - len(b))
    return BP([arb_poly([a[j], b[j]]) for j in range(n)], True)


def lin_s(c0, c1):
    """c0(tc) + c1(tc) s"""
    a = c0.coeffs() or [arb(0)]; b = c1.coeffs() or [arb(0)]
    n = max(len(a), len(b)); a = a + [arb(0)] * (n - len(a)); b = b + [arb(0)] * (n - len(b))
    return BP([arb_poly([a[j], b[j]]) for j in range(n)], True)


def element_a(B, e, Qe, rte, the):
    """type-(a) quantities of element e; Qe, rte, the: TPs (Q, rho_theta, theta at x_e)"""
    D = B.D; Ns = B.Ns
    d, dd = B.local_coeffs(e)
    h, ih = B.h[e], B.ih[e]
    f2 = ih * TWO                                       # ds/dx
    f2_2 = T.tmul(f2, f2, D); f2_3 = T.tmul(f2_2, f2, D)
    kap = expand(d, Ns, D)
    ks = kap.ds(); kss = ks.ds(); ksss = kss.ds(); k4s = ksss.ds(); k5s = k4s.ds()
    kx = ks.tmul(f2, D); kxx = kss.tmul(f2_2, D); kxxx = ksss.tmul(f2_3, D)
    # Q = Qe + (2/h) int_{-1}^s kappa_s^2
    Q = BP.from_tp(Qe) + ks.mul(ks, D).int_m1().tmul(f2, D)
    x = xpoly(B.x[e], h)
    c = B.c
    k2 = kap.mul(kap, D)
    Phi = (-kxxx) - k2.mul(kx, D) + kap.mul(Q, D) - x.mul(kap, D).tmul(c, D)
    Phix = Phi.ds().tmul(f2, D)
    # kappa_tau at fixed x:  sum dd_q N_q + elld * ALE
    gl = B.gl[e]
    fac = []
    for q in range(len(d)):
        if q < 5: fac.append(T.tmul(d[q], (gl - B.dHH[e]) * arb(q), D))
        elif q < 10: fac.append(T.tmul(d[q], (gl - B.dHH[e + 1]) * arb(q - 5), D))
        else: fac.append(T.tconst(arb(0)))
    ALE = expand(fac, Ns, D) - kx.mul(lin_s(B.dx[e] + B.dh[e] * HALF, B.dh[e] * HALF), D)
    kt = expand(dd, Ns, D) + ALE.tmul(B.elld, D)
    rk = kt - Phix
    hh = h * HALF
    rt = BP.from_tp(rte) + rk.int_m1().tmul(hh, D)
    th = BP.from_tp(the) + kap.int_m1().tmul(hh, D)
    return dict(kap=kap, ks=ks, kss=kss, ksss=ksss, k4s=k4s, k5s=k5s, kx=kx, kxx=kxx, kxxx=kxxx, Q=Q, x=x,
                Phi=Phi, Phix=Phix, kt=kt, rk=rk, rt=rt, th=th, f2=f2, f2_2=f2_2, f2_3=f2_3, k2=k2, hh=hh)


# ------------------------------------------------------------------------------------------- rotation E = e^{i phi}
def exp_i(phi, DS, DT, tol=1e-26):
    """E = exp(i phi(s, tc)) as a pointwise CBP (s-degree DS, tc-degree DT), phi: exact BP.
    phi = phi0(tc) + psi, psi with zero s-constant;  e^{i psi} = sum F_k s^k by (k+1) F_{k+1} = i sum (m+1) psi_{m+1} F_{k-m};
    tail for |s| <= 1 by Cauchy: |F_k| <= M_R R^-k, M_R <= exp(sum_k sup|psi_k| R^k)."""
    # s-coefficients of phi as TPs
    sdeg = phi.sdeg()
    cols = []
    for k in range(sdeg + 1):
        cols.append(T.ttrunc(arb_poly([p.coeffs()[k] if k < len(p.coeffs()) else arb(0) for p in phi.P]), DT))
    phi0 = cols[0]
    psi = cols[1:]                      # psi[m] = coefficient of s^(m+1)
    # adaptive s-degree: smallest DS (>= 8, <= DSmax) whose Cauchy tail is below tol
    sups0 = [float(T.tsup(p).upper()) for p in psi]
    def tail_est(ds):
        b = math.inf
        for R in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
            S = sum(sv * R ** (m + 1) for m, sv in enumerate(sups0))
            if S > 600: continue
            b = min(b, S - (ds + 1) * math.log(R) - math.log(1 - 1 / R))
        return b
    DSmax = DS
    DS = 8
    while DS < DSmax and tail_est(DS) > math.log(tol): DS += 2
    Fr = [T.tconst(arb(1))]; Fi = [T.tconst(arb(0))]
    for k in range(DS):
        ar = arb_poly([0]); ai = arb_poly([0])
        for m in range(min(k + 1, len(psi))):
            w = psi[m] * arb(m + 1)
            ar = ar + T.tmul(w, Fi[k - m], DT)
            ai = ai + T.tmul(w, Fr[k - m], DT)
        inv = arb(1) / (k + 1)
        Fr.append(ar * (-inv)); Fi.append(ai * inv)
    # Cauchy tail
    sups = [float(T.tsup(p).upper()) for p in psi]
    best = math.inf; bestR = None
    for R in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
        S = sum(sv * R ** (m + 1) for m, sv in enumerate(sups))
        if S > 600: continue
        lb = S - (DS + 1) * math.log(R) - math.log(1 - 1 / R)
        if lb < best: best, bestR = lb, R
    if bestR is None: raise ValueError('exp_i: phi too large')
    R = arb(bestR)
    S = sum((T.tsup(p) * R ** (m + 1) for m, p in enumerate(psi)), arb(0))
    tail = (S.exp() * R ** (-(DS + 1)) / (1 - 1 / R)).upper()
    # assemble as BPs (rows: tc power, poly in s)
    def tobp(F):
        nt = max(len(f.coeffs()) for f in F)
        P = []
        for j in range(nt):
            P.append(arb_poly([f.coeffs()[j] if j < len(f.coeffs()) else arb(0) for f in F]))
        return BP(P, False)
    Er, Ei = tobp(Fr), tobp(Fi)
    Er = Er.fold_pointwise(tail); Ei = Ei.fold_pointwise(tail)
    c0, s0 = T.tcis(phi0, DT)
    E = CBP(Er, Ei).ctmul(c0, s0, DT)
    return E, dict(tail=float(tail), R=bestR, phi0=phi0, DS=DS)


def element_rot(B, e, A_, ze, Ge, Pe, DS=30, DTE=8, DJ=8, DTZ=4, DW=16):
    """rotation part; A_: dict from element_a; ze, Ge, Pe: complex TPs (pairs) at x_e.
    DTE: tc-degree of E (its tc-variation is small; the tail is folded); DJ: tc-degree of the geometry (G, Pg)."""
    D = B.D
    hh = A_['hh']
    phi = A_['th'] - BP.from_tp(A_['th'].at(MONE))
    E, info = exp_i(phi, DS, DTE)
    Ez = E.ttrunc(DTZ)                       # the residual products need E only to ~1e-15 relative
    W = Ez.conj().rmul(A_['rt'].ttrunc(DW), DW).int_m1().tmul(hh, DW)
    W = CBP(W.re + BP.from_tp(T.ttrunc(ze[0], DW)), W.im + BP.from_tp(T.ttrunc(ze[1], DW)))
    z = Ez.mul(W, DW)
    # position in the local tangent frame (geometry: tc-degree DJ suffices)
    EJ = CBP(E.re.ttrunc(DJ), E.im.ttrunc(DJ))
    I = EJ.int_m1().tmul(T.ttrunc(hh, DJ), DJ)
    Gs = CBP(I.re + BP.from_tp(T.ttrunc(Ge[0], DJ)), I.im + BP.from_tp(T.ttrunc(Ge[1], DJ)))
    G = EJ.conj().mul(Gs, DJ)
    # global e^{i theta} = Pe E
    Pg = EJ.ctmul(T.ttrunc(Pe[0], DJ), T.ttrunc(Pe[1], DJ), DJ)
    return dict(E=E, z=z, G=G, Pg=Pg, info=info)
