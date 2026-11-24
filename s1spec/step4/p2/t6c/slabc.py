"""Direct T6 residual of the (SE)-corrected base (review B1): run_slab_c(i) -> (summary dict, jets dict), as T6's run_slab.

The corrected base is the one of p2/se_corr.py:
    kap^ = kap~ + dk,   dk = (a1/ell) B(u) + (a2/ell^2) B'(u),   u = x/ell(tau),
    B(u) = psi((u - 0.6)/0.2) / (0.2 I_psi),  psi = t^6 (1-t)^6,  I_psi = 1/12012   (int B du = 1, int B' du = 0),
so that Dtheta := theta^ - theta~ = a1 Bbar(u) + (a2/ell) B(u), Bbar(u) = int_0.6^u B (Dtheta = a1 for u >= 0.8).
T6's position is W(x) = -DS1 + int_0^x e^{i theta} (tip fixed, theta(0) = -pi/2), so the axis is Im W = 0 and (SE) holds when
    theta^(ell) = -3 pi/2   <=>  a1 = -delta,  delta := theta~(ell) + 3 pi/2  (exact: Dtheta(ell) = a1),
    Im W^(ell) = 0          <=>  F(a2) := h~ + a1 ell J1 + a2 J2 + N(a2) = 0,
with h~ = Im W~(ell), J1 = int_0.6^1 cos theta~(ell u) Bbar du, J2 = int_0.6^0.8 cos theta~ B du and
N = Im int_0^ell e^{i theta~} (e^{i Dtheta} - 1 - i Dtheta) dx, |N| <= int Dtheta^2/2.  For every tc, a2 -> -(h~ + a1 ell J1 + N)/J2
maps the ball a2 in m + [-eps, eps] (m = -(h~ + a1 ell J1)/J2) into itself and contracts (|dN/da2| <= sup|Dtheta| << |J2|),
so the root exists and is unique there (Banach); it is C^1 in tau by the implicit function theorem, with
    a2' = -(h~' + (a1 ell J1)' + a2 J2' + d_tau N) / (J2 + d_a2 N),
    h~' = Re int_0^ell theta~_tau e^{i theta~} dx + ell' sin theta~(ell),   theta~_tau|x = rho_theta + Phi - Phi(0)   (T6),
    a1' = -(theta~_tau|x + ell' kap~)(ell),  J' from d/dtau of the u-integrals (no boundary terms).
All of these are TPs (T6's tp.py) from pass 1, the T6 run of the base.  Pass 2 re-runs T6 on the corrected base from the
element that contains u = 0.6 to the far end: the elements that contain u = 0.6 and u = 0.8 are split there (the pieces are
compositions s = a(tc) + w(tc) sigma of the parent; x-derivatives are composed from the parent, so no 1/h of a piece is used),
and on the window pieces kap^ = kap~ + dk, kap^_tau|x = kap~_tau|x + d_tau dk|x (d_tau dk from a1, a1', a2, a2', ell, ell').
Every norm, scalar and jet then comes from T6's own code (slab.py), now for the corrected base.
Slab 11: the knot 92 crosses u = 0.6 inside the slab (unique crossing: g' = x_92' - 0.6 ell' has a fixed sign), so the slab is
cut in time at the crossing t* (an arb ball) into two sub-slabs, each with a fixed element structure; norms are the max over the
two, norms_int the root of the sum of squares, jets the union."""
import os, sys, math, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S4 = os.path.normpath(os.path.join(HERE, '..', '..'))
T6DIR = os.path.join(S4, 't6')
if T6DIR not in sys.path: sys.path.insert(0, T6DIR)
from flint import arb, arb_poly, fmpq, fmpq_poly, ctx
ctx.prec = 128                                      # before the module constants (T6 sets 128 only when a slab is loaded)
import tp as T
from tp import BP, CBP, strunc
from base import SlabBase
from engine import element_a, element_rot, xpoly, lin_s, exp_i
import weights as WT
import slab as S6

ONE, ZERO = arb(1), arb(0)
HALF = arb(1) / 2
MONE = -ONE
U1, U2 = arb(3) / 5, arb(4) / 5                    # window u in [0.6, 0.8]
_t = fmpq_poly([0, 1]); _PSI = _t ** 6 * (1 - _t) ** 6
CB = 60060                                          # 1/(0.2 I_psi), I_psi = 1/12012
BMAX = arb(CB) / 4096                               # sup B = CB psi(1/2)
B1MAX = arb(CB) * 30 / 1024                         # sup |B'| <= CB 5 * 6 sup t^5(1-t)^5
BPOLS = []
_p = _PSI
for _k in range(8):                                 # B^(k)(u) = CB 5^k psi^(k)(t), t = 5u - 3
    BPOLS.append([arb(c.p) * (CB * 5 ** _k) / arb(c.q) for c in _p.coeffs()]); _p = _p.derivative()
_I = _PSI.integral()
BBAR = [arb(c.p) * CB / (5 * arb(c.q)) for c in _I.coeffs()]          # Bbar = (CB/5) int_0^t psi (exact rationals -> balls)
PI32 = 3 * arb.pi() / 2


# ------------------------------------------------------------------------------------------------------------ BP tools
def cpoly(c):
    return BP.from_tp(T.tconst(c))


def poly_bp(coeffs, t, D):
    """sum_k coeffs[k] t^k for an exact BP t (Horner)"""
    acc = cpoly(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        acc = acc.mul(t, D) + cpoly(c)
    return acc


def compose(F, a, w, D):
    """F(a(tc) + w(tc) sigma, tc) as a BP in sigma (a, w TPs; the piece lies inside the parent, |a| + |w| <= 1)"""
    sd = F.sdeg()
    if sd < 0: return BP.zero()
    cols = []
    for k in range(sd + 1):
        cols.append(arb_poly([p.coeffs()[k] if k < len(p.coeffs()) else ZERO for p in F.P]))
    L = lin_s(a, w)
    acc = BP.from_tp(cols[sd])
    for k in range(sd - 1, -1, -1):
        acc = acc.mul(L, D) + BP.from_tp(cols[k])
    return BP(acc.P, F.ex)


def ccompose(Fc, a, w, D):
    return CBP(compose(Fc.re, a, w, D), compose(Fc.im, a, w, D))


def tball(p, r):
    """TP p with r >= 0 added to the radius of its constant coefficient"""
    c = p.coeffs() or [ZERO]
    c[0] = T.addrad(c[0], arb(r).abs_upper() if not isinstance(r, float) else arb(r))
    return arb_poly(c)


def tinv(p, D):
    """1/p for a TP of fixed sign"""
    m = T.tmid0(p)
    return T.tpow(p, -1, D) if m > 0 else -T.tpow(-p, -1, D)


def trng(p):
    return T.trange(p)


def parent_derivs(B, e, st):
    """on parent element e, with the BASE state st at x_e: the x-derivatives kap..kap^(5), the base residual rho_kappa~ and
    its x-derivative (all BPs in the parent s, computed by engine.element_a exactly as T6 does)"""
    A = std_element(B, e, st)
    return [A['kap'], A['kx'], A['kxx'], A['kxxx'], A['k4'], A['k5']], A['rk'], A['rkx']


def element_c(Kb, dK, dKT, rkb, rkxb, xe, h, c, st, Qbe, D):
    """T6 type-(a) quantities of a piece of the corrected base.
    Kb: base x-derivatives kap~..kap~^(5) on the piece; dK, dKT: the bump's x-derivatives and (d_tau dk, d_x d_tau dk), or None;
    rkb, rkxb: the base residual rho_kappa~ (with the base Q~) and its x-derivative on the piece; st: corrected state at the
    piece start; Qbe: the base Q~ there.  rho_kappa^ = rho_kappa~ + d_tau dk - DPhi_x, with DPhi_x = Phi_x[kap~ + dk, Q^] -
    Phi_x[kap~, Q~] expanded in (dk, DQ) so that no O(1) terms cancel; Phi_x = -k4 - k k1^2 - k^2 k2 + k1 Q - c (k + x k1) and
    Phi_xx = -k5 - 4 k k1 k2 - k^2 k3 + k2 Q - c (2 k1 + x k2) (k_m: m-th x-derivative; Q_x = k1^2)."""
    hh = h * HALF
    x = xpoly(xe, h)
    kb, kb1, kb2, kb3 = Kb[0], Kb[1], Kb[2], Kb[3]
    dQ = BP.from_tp(st['Qe'] - Qbe)
    if dK is None: K = Kb
    else:
        K = [p + q for p, q in zip(Kb, dK)]
        d0, d1, d2, d3, d4, d5 = dK
        dQ = dQ + d1.mul(kb1.scale(arb(2)) + d1, D).int_m1().tmul(hh, D)
    kap, kx, kxx, kxxx, k4, k5 = K
    k2 = kap.mul(kap, D); kx2 = kx.mul(kx, D)
    Q = BP.from_tp(st['Qe']) + kx2.int_m1().tmul(hh, D)
    Phi = (-kxxx) - k2.mul(kx, D) + kap.mul(Q, D) - x.mul(kap, D).tmul(c, D)
    DPx = kb1.mul(dQ, D); DPxx = kb2.mul(dQ, D)
    rk, rkx = rkb, rkxb
    if dK is not None:
        kb_sq = kb.mul(kb, D)
        DPx = (DPx - d4 - d0.mul(kx2, D) - kb.mul(d1.mul(kb1.scale(arb(2)) + d1, D), D)
               - d0.mul(kb.scale(arb(2)) + d0, D).mul(kxx, D) - kb_sq.mul(d2, D) + d1.mul(Q, D)
               - (d0 + x.mul(d1, D)).tmul(c, D))
        DPxx = (DPxx - d5 - (d0.mul(kx, D).mul(kxx, D) + kb.mul(d1.mul(kxx, D) + kb1.mul(d2, D), D)).scale(arb(4))
                - d0.mul(kb.scale(arb(2)) + d0, D).mul(kxxx, D) - kb_sq.mul(d3, D) + d2.mul(Q, D)
                - (d1.scale(arb(2)) + x.mul(d2, D)).tmul(c, D))
        rk = rk + dKT[0]; rkx = rkx + dKT[1]
    rk = rk - DPx; rkx = rkx - DPxx
    rt = BP.from_tp(st['rte']) + rk.int_m1().tmul(hh, D)
    th = BP.from_tp(st['the']) + kap.int_m1().tmul(hh, D)
    return dict(kap=kap, kx=kx, kxx=kxx, kxxx=kxxx, k4=k4, k5=k5, Q=Q, x=x, Phi=Phi, rk=rk, rkx=rkx, rt=rt, th=th, k2=k2, hh=hh)


def std_element(B, e, st):
    """engine.element_a for a whole element, with rkx, k4, k5 added (as slab.py computes them)"""
    D = B.D
    A = element_a(B, e, st['Qe'], st['rte'], st['the'])
    f4 = T.tmul(A['f2_2'], A['f2_2'], D); f5 = T.tmul(A['f2_3'], A['f2_2'], D)
    A['rkx'] = A['rk'].ds().tmul(A['f2'], D); A['k4'] = A['k4s'].tmul(f4, D); A['k5'] = A['k5s'].tmul(f5, D)
    return A


# ------------------------------------------------------------------------------------------------------------ the bump
class Bump:
    """dk and d_tau dk on a piece with x = xe + h(1+sigma)/2, given the TPs a1, a1', a2, a2' (and ell, ell')"""
    def __init__(self, B, a1, a1p, a2, a2p):
        D = B.D; self.D = D
        L, Lp = B.ell, B.elld
        il = T.tpow(L, -1, D); self.il = il
        self.ilp = [T.tconst(ONE)]
        for _ in range(7): self.ilp.append(T.tmul(self.ilp[-1], il, D))
        il2, il3 = self.ilp[2], self.ilp[3]
        self.A1 = T.tmul(a1, il, D); self.A2 = T.tmul(a2, il2, D)
        self.A1p = T.tmul(a1p, il, D) - T.tmul(T.tmul(a1, Lp, D), il2, D)
        self.A2p = T.tmul(a2p, il2, D) - T.tmul(T.tmul(a2, Lp, D), il3, D) * 2
        self.lr = T.tmul(Lp, il, D)                  # ell'/ell
        self.lr2 = T.tmul(Lp, il2, D)                # ell'/ell^2

    def on(self, xe, h):
        D = self.D
        u = xpoly(xe, h).tmul(self.il, D)
        t = u.scale(arb(5)) - cpoly(arb(3))
        Bk = [poly_bp(BPOLS[k], t, D) for k in range(8)]
        ut = u.tmul(-self.lr, D)                     # u_tau = -u ell'/ell
        dK = [Bk[m].tmul(T.tmul(self.A1, self.ilp[m], D), D) + Bk[m + 1].tmul(T.tmul(self.A2, self.ilp[m], D), D) for m in range(6)]
        dkt = (Bk[0].tmul(self.A1p, D) + Bk[1].mul(ut, D).tmul(self.A1, D) + Bk[1].tmul(self.A2p, D)
               + Bk[2].mul(ut, D).tmul(self.A2, D))
        dktx = (Bk[1].tmul(T.tmul(self.A1p, self.il, D), D)
                + (Bk[2].mul(ut, D).tmul(self.il, D) - Bk[1].tmul(self.lr2, D)).tmul(self.A1, D)
                + Bk[2].tmul(T.tmul(self.A2p, self.il, D), D)
                + (Bk[3].mul(ut, D).tmul(self.il, D) - Bk[2].tmul(self.lr2, D)).tmul(self.A2, D))
        return dK, [dkt, dktx]


def window_weights(B, xe, h):
    """B(u), Bbar(u), u on a piece (pass 1)"""
    D = B.D
    il = T.tpow(B.ell, -1, D)
    u = xpoly(xe, h).tmul(il, D)
    t = u.scale(arb(5)) - cpoly(arb(3))
    return poly_bp(BPOLS[0], t, D), poly_bp(BBAR, t, D), u


# ------------------------------------------------------------------------------------------------------------ geometry
def split_points(B):
    D = B.D
    return T.tmul(B.ell, T.tconst(U1), D), T.tmul(B.ell, T.tconst(U2), D)


def find_parent(B, X, forced=None):
    """index e with x_e < X < x_{e+1} for every tc (T6 knots), or the forced index"""
    if forced is not None: return forced
    D = B.D
    hits = [e for e in range(B.mesh.JT, B.N) if trng(B.x[e] - X).upper() < 0 and trng(B.x[e + 1] - X).lower() > 0]
    if len(hits) != 1: raise ValueError(f'slab {B.i}: split point not inside one element over the slab ({hits})')
    return hits[0]


def pieces(B, e, X, D):
    """the two pieces of parent e cut at X: (xe, h, a, w) with s = a + w sigma; and the TP s* of the cut"""
    ss = T.tmul(X - B.x[e], B.ih[e], D) * 2 - T.tconst(ONE)          # s* = 2 (X - x_e)/h_e - 1
    a1 = (ss - T.tconst(ONE)) * HALF; w1 = (ss + T.tconst(ONE)) * HALF
    a2 = (ss + T.tconst(ONE)) * HALF; w2 = (T.tconst(ONE) - ss) * HALF
    return [(B.x[e], X - B.x[e], a1, w1), (X, B.x[e + 1] - X, a2, w2)], ss


# ------------------------------------------------------------------------------------------------------------ per piece: T6's slab.py block
class Ctx:
    def __init__(self, B, o):
        self.B, self.o = B, o
        D, J, MT = o['D'], o['J'], o['MT']
        self.tsub = [(-1 + 2 * m / MT, -1 + 2 * (m + 1) / MT) for m in range(MT)]
        self.ellr = [S6.tp_at_ball(B.ell, a, b) for a, b in self.tsub]
        self.llr = [(WT.lo_(e - WT.X0), WT.up(e - WT.X0)) for e in self.ellr]
        ell_all = S6.tp_at_ball(B.ell, -1, 1); self.ll_all = (WT.lo_(ell_all - WT.X0), WT.up(ell_all - WT.X0))
        self.acc = S6.Acc(MT)
        self.crude = dict(supr_tip=0.0, supA_tip=0.0, supr1_tip=0.0, supA1_tip=0.0, supArho_tip=0.0,
                          supA_far=0.0, supAx_far=0.0, supA1_far=0.0, supr_far=0.0)
        self.maxrad = dict(rk=0.0, r=0.0, S2=0.0, S1=0.0)
        self.s_deg = dict(r=0, r3=0, E=0)

    def snapshot(self):
        return (dict(self.acc.d), dict(self.crude), dict(self.maxrad), dict(self.s_deg))

    def restore(self, snap):
        self.acc.d = dict(snap[0]); self.crude = dict(snap[1]); self.maxrad = dict(snap[2]); self.s_deg = dict(snap[3])


def post(C, A_, R_, xe, h, tipel, st, norms=True):
    """slab.py's per-element block (r derivatives, S1, S2, weighted norm polynomials, crude sups) for one piece; returns the new
    state and the jet quantities (BPs)"""
    B, o = C.B, C.o
    D, J, MT = o['D'], o['J'], o['MT']
    hh = A_['hh']
    z = R_['z']
    r = strunc(z.re, o['strtol']); Aim = strunc(z.im, o['strtol'])
    rt = strunc(A_['rt'], o['strtol'])
    kap, kx, kxx, k2 = A_['kap'], A_['kx'], A_['kxx'], A_['k2']
    rk = A_['rk']; rkx = A_['rkx']
    rkp = strunc(rk.ttrunc(J), o['strtol']); rkxp = strunc(rkx.ttrunc(J), o['strtol'])
    kapJ, kxJ, kxxJ, k2J = kap.ttrunc(J), kx.ttrunc(J), kxx.ttrunc(J), k2.ttrunc(J)
    rJ, AJ_, rtJ = r.ttrunc(J), Aim.ttrunc(J), rt.ttrunc(J)
    kkxJ = kapJ.mul(kxJ, J)
    r1 = strunc(rtJ - kapJ.mul(AJ_, J), o['strtol'])
    r2 = strunc(rkp - kxJ.mul(AJ_, J) - k2J.mul(rJ, J), o['strtol'])
    r3 = strunc(rkxp - kxxJ.mul(AJ_, J) - kkxJ.mul(rJ, J).scale(arb(3)) - k2J.mul(r1, J), o['strtol'])
    A1 = kapJ.mul(rJ, J)
    S1 = strunc(BP.from_tp(st['S1e']) + r.int_m1().tmul(hh, D), o['strtol'])
    S2 = strunc(BP.from_tp(st['S2e']) + S1.int_m1().tmul(hh, D), o['strtol'])
    x_end = xe + h
    xe0, xe1 = S6.tp_at_ball(xe, -1, 1), S6.tp_at_ball(x_end, -1, 1)
    xlo_el, xhi_el = WT.lo_(xe0), WT.up(xe1)
    if norms:
        if tipel:
            xa, xb = xe(ZERO), x_end(ZERO)
            pts = [xa] + [c for c in S6.CUTS if c > xa and c < xb] + [xb]
            sb = [MONE]
            for k in range(1, len(pts)):
                seg = pts[k] - pts[k - 1]
                n = max(1, int(math.ceil(float(seg.mid()) / o['DXT'])))
                for q in range(1, n + 1):
                    xq = pts[k - 1] + seg * q / n
                    sq_ = arb(((2 * (xq - xa) / (xb - xa)) - 1).mid()) if (k < len(pts) - 1 or q < n) else ONE
                    sb.append(sq_)
            sb[-1] = ONE
        else:
            sb = [MONE + 2 * arb(q) / o['NCO'] for q in range(o['NCO'] + 1)]
        ncell = len(sb) - 1
        need = []
        if xlo_el < S6.XRHO: need += ['S2', 'r']
        if xhi_el > S6.XW: need += ['S1', 'r1', 'r3']
        if xhi_el > S6.XF: need += ['r', 'S1', 'r1', 'r2', 'r3']
        need = [f for f in S6.FUN_ALL if f in need]
        F = dict(S2=S2, r=r, S1=S1, r1=r1, r2=r2, r3=r3)
        sq = {}
        hJ = T.ttrunc(hh, 2 * J)
        for f in need:
            g = F[f].ttrunc(J)
            g2 = g.mul(g, 2 * J).int_m1()
            vals = [g2.at(sv) for sv in sb]
            sq[f] = [T.tmul(vals[c + 1] - vals[c], hJ, 2 * J) for c in range(ncell)]
        if 'S1' in need:
            g1 = F['S1'].ttrunc(2 * J).int_m1()
            v1 = [g1.at(sv) for sv in sb]
            M1 = [T.tmul(v1[c + 1] - v1[c], hJ, 2 * J) for c in range(ncell)]
            Lc = [hJ * (sb[c + 1] - sb[c]) for c in range(ncell)]
        msubs = range(MT) if not tipel else [None]
        for c in range(ncell):
            for m in msubs:
                if tipel:
                    xa_c = xe(ZERO) + h(ZERO) * HALF * (1 + sb[c]); xb_c = xe(ZERO) + h(ZERO) * HALF * (1 + sb[c + 1])
                    xa_c, xb_c = WT.lo_(xa_c), WT.up(xb_c); lla, llb = C.ll_all
                else:
                    ta, tb = C.tsub[m]
                    xb_ = S6.tp_at_ball(xe, ta, tb); hb_ = S6.tp_at_ball(h, ta, tb)
                    xa_c = WT.lo_(xb_ + hb_ * HALF * (1 + sb[c])); xb_c = WT.up(xb_ + hb_ * HALF * (1 + sb[c + 1]))
                    lla, llb = C.llr[m]
                ms = range(MT) if tipel else [m]
                ya, yb = xa_c - WT.X0, xb_c - WT.X0
                items = []
                if 'S2' in need and xa_c < S6.XRHO:
                    wr = WT.rho_sup(xa_c, xb_c)
                    items += [('R2', 'S2', wr), ('Rs', 'r', wr), ('Rs_chi', 'r', WT.rho_chi2_sup(xa_c, xb_c))]
                if xb_c > S6.XW and 'S1' in need:
                    ww = WT.Wzeta2_sup(xa_c, xb_c)
                    if ww > 0: items += [('R1', 'S1', ww), ('RIs', 'r1', ww), ('RIs3', 'r3', ww)]
                if xb_c > S6.XF and 'r2' in need:
                    for j, f in enumerate(('S1', 'r', 'r1', 'r2', 'r3')):
                        wl = WT.Tlit_sup(j, xa_c, xb_c)
                        if wl > 0: items.append((f'RFlit{j}', f, wl))
                        for dn, ramps in WT.DESIGNS.items():
                            wt = WT.T_sup(j, ya, yb, lla, llb, ramps[j])
                            if wt > 0:
                                fj = 'S1sh' if j == 0 else f
                                items.append((f'RF{j}_{dn}', fj, wt))
                                hA, hB = WT.hats_sup(j, ya, yb, lla, llb, ramps[j])
                                items.append((f'RFa{j}_{dn}', fj, hA)); items.append((f'RFb{j}_{dn}', fj, hB))
                for name, f, w in items:
                    if f == 'S1sh':
                        for mm in ms:
                            C.acc.add(name + '|P', mm, sq['S1'][c] * w); C.acc.add(name + '|M', mm, M1[c] * w)
                            C.acc.add(name + '|L', mm, Lc[c] * w)
                    else:
                        for mm in ms: C.acc.add(name, mm, sq[f][c] * w)
        up = S6.up
        sr, sA, sr1 = up(r.sup()), up(Aim.sup()), up(r1.sup())
        sA1 = up(A1.sup())
        cr = C.crude
        if xlo_el < S6.XRHO:
            cr['supr_tip'] = max(cr['supr_tip'], float(sr)); cr['supA_tip'] = max(cr['supA_tip'], float(sA))
            cr['supr1_tip'] = max(cr['supr1_tip'], float(sr1)); cr['supA1_tip'] = max(cr['supA1_tip'], float(sA1))
            rr = up(xhi_el / 3 * (1 + xlo_el ** 2) ** (-arb(1) / 3))
            cr['supArho_tip'] = max(cr['supArho_tip'], float(up(sA * rr)))
        if xhi_el > S6.XW:
            cr['supA_far'] = max(cr['supA_far'], float(sA)); cr['supA1_far'] = max(cr['supA1_far'], float(sA1))
            cr['supr_far'] = max(cr['supr_far'], float(sr))
            cr['supAx_far'] = max(cr['supAx_far'], float(up(sA / max(xlo_el, S6.XW))))
        mr = C.maxrad
        mr['rk'] = max(mr['rk'], max(float(cf.rad()) for cf in rk.P[0].coeffs()))
        mr['r'] = max(mr['r'], max(float(cf.rad()) for cf in r.P[0].coeffs()))
        mr['S1'] = max(mr['S1'], max(float(cf.rad()) for cf in S1.P[0].coeffs()))
        mr['S2'] = max(mr['S2'], max(float(cf.rad()) for cf in S2.P[0].coeffs()))
        C.s_deg['r'] = max(C.s_deg['r'], r.sdeg()); C.s_deg['r3'] = max(C.s_deg['r3'], r3.sdeg())
        C.s_deg['E'] = max(C.s_deg['E'], R_['info']['DS'])
    # jet quantities (slab.py)
    JTd = o['JT']
    G, Pg = R_['G'], R_['Pg']
    c_, b_ = B.c, B.b
    x = A_['x']
    p_, q_ = G.im.ttrunc(JTd), G.re.ttrunc(JTd)
    t1_, n_ = Pg.re.ttrunc(JTd), (-Pg.im).ttrunc(JTd)
    cJ, bJ = T.ttrunc(c_, JTd), T.ttrunc(b_, JTd)
    kxJ4, kxxJ4 = kx.ttrunc(JTd), kxx.ttrunc(JTd)
    kkx4 = kkxJ.ttrunc(JTd)
    Q_J = A_['Q'].ttrunc(JTd)
    Vs = (-kxxJ4) + p_.tmul(cJ, JTd) + n_.tmul(bJ, JTd)
    ws = (-kkx4) + Q_J + (q_ - x.ttrunc(JTd)).tmul(cJ, JTd) + t1_.tmul(bJ, JTd)
    rJ4, AJ4 = r.ttrunc(JTd), Aim.ttrunc(JTd)
    b1 = (-(kkx4.scale(arb(6)))) + Q_J - x.ttrunc(JTd).tmul(cJ, JTd) + AJ4
    quants = [('x', x), ('k0', kap), ('k1', kx), ('k2', kxx), ('k3', A_['kxxx']), ('k4', A_['k4']), ('k5', A_['k5']), ('Q', A_['Q']),
              ('theta', A_['th']), ('p', p_), ('q', q_), ('n', n_), ('t1', t1_), ('V', Vs + rJ4), ('varpi', ws + AJ4),
              ('b1', b1), ('Vstar', Vs), ('varpistar', ws), ('r', r), ('A', Aim), ('r1', r1), ('r2', r2), ('r3', r3), ('A1', A1)]
    new = dict(Qe=A_['Q'].at(ONE), rte=A_['rt'].at(ONE), the=A_['th'].at(ONE), ze=z.at(ONE), Ge=G.at(ONE), Pe=Pg.at(ONE),
               S1e=S1.at(ONE), S2e=S2.at(ONE))
    return new, quants


def jet_range(f, JTd, slo, shi):
    """ball containing f over sigma in [slo, shi] and tc in [-1, 1] (slab.py's jet hull)"""
    fJ = f.ttrunc(JTd)
    sball = WT.iv(slo, shi)
    vals = [pp(sball) for pp in fJ.P]
    return vals[0] + arb(0, T.usum(v.abs_upper() for v in vals[1:])) if len(vals) > 1 else vals[0]


def jets_whole(quants, NC, JTd):
    out = {}
    for nm, f in quants:
        out[nm] = [jet_range(f, JTd, MONE + 2 * arb(cc) / NC, MONE + 2 * arb(cc + 1) / NC) for cc in range(NC)]
    return out


def jets_split(qs, ss, NC, JTd):
    """parent rows from the two pieces (qs[0]: s in [-1, s*], qs[1]: [s*, 1]); union over the pieces that meet each sub-cell"""
    r_ = trng(ss); slo_, shi_ = r_.lower(), r_.upper()
    out = {nm: [None] * NC for nm, _ in qs[0]}
    for cc in range(NC):
        sc, sd = MONE + 2 * arb(cc) / NC, MONE + 2 * arb(cc + 1) / NC
        rng = []
        # piece 1: sigma = -1 + 2 (s + 1)/(s* + 1), decreasing in s*
        if sc < shi_:
            if slo_ + 1 > 0:
                lo = max(MONE, MONE + 2 * (sc + 1) / (shi_ + 1)); hi = min(ONE, MONE + 2 * (sd + 1) / (slo_ + 1))
            else: lo, hi = MONE, ONE
            if lo <= hi: rng.append((0, lo, hi))
        # piece 2: sigma = -1 + 2 (s - s*)/(1 - s*), decreasing in s*
        if sd > slo_:
            if 1 - shi_ > 0:
                lo = max(MONE, MONE + 2 * (sc - shi_) / (1 - shi_)); hi = min(ONE, MONE + 2 * (sd - slo_) / (1 - slo_))
            else: lo, hi = MONE, ONE
            if lo <= hi: rng.append((1, lo, hi))
        assert rng, 'jets_split: empty sub-cell'
        for k, (nm, _) in enumerate(qs[0]):
            v = None
            for (pi, lo, hi) in rng:
                w = jet_range(qs[pi][k][1], JTd, lo, hi)
                v = w if v is None else v.union(w)
            out[nm][cc] = v
    return out


# ------------------------------------------------------------------------------------------------------------ one (sub)slab
def geo_step(A_, g, Phi0, o, D):
    """accurate geometry along the base (tc-degree D, not T6's 8): W += int e^{i theta}, IW += int theta_tau|x e^{i theta},
    P = e^{i theta} at the element end.  g = (W, P, IW) as pairs of TPs.  (T6's DJ = 8 gives e^{i theta} to ~1e-15 only.)"""
    W, P, IW = g                                          # P = e^{i theta(x_e)}: tcis of theta at the previous end (radius ~1e-34)
    phi = A_['th'] - BP.from_tp(A_['th'].at(MONE))
    E, info = exp_i(phi, o['DSmax'], D)
    hh = A_['hh']
    I = E.int_m1().tmul(hh, D)
    tht = A_['rt'] + A_['Phi'] - BP.from_tp(Phi0)
    IT = E.rmul(tht, D).int_m1().tmul(hh, D)
    cm = lambda z, w: (T.tmul(z[0], w[0], D) - T.tmul(z[1], w[1], D), T.tmul(z[0], w[1], D) + T.tmul(z[1], w[0], D))
    dW = cm(P, I.at(ONE)); dI = cm(P, IT.at(ONE))
    return (W[0] + dW[0], W[1] + dW[1]), T.tcis(A_['th'].at(ONE), D), (IW[0] + dI[0], IW[1] + dI[1])


def run_part(B, o, forced=None, verbose=False):
    """the corrected base on one (sub)slab B; forced = (ea, eb) parent indices (slab 11)"""
    D, J, MT = o['D'], o['J'], o['MT']
    N = B.N; JTn = B.mesh.JT
    Z = T.tconst(ZERO)
    st = dict(Qe=Z, rte=Z, the=T.tconst(-arb.pi() / 2), ze=(Z, Z), Ge=(Z, T.tconst(-arb(2.500806276619802))),
              Pe=(Z, T.tconst(-ONE)), S1e=Z, S2e=Z)
    C = Ctx(B, o)
    rot = lambda A_, s_: element_rot(B, 0, A_, s_['ze'], s_['Ge'], s_['Pe'], DS=o['DSmax'], DTE=o['DTE'], DJ=o['DJ'])
    X1, X2 = split_points(B)
    ea = find_parent(B, X1, forced[0] if forced else None); eb = find_parent(B, X2, forced[1] if forced else None)
    assert JTn <= ea < eb < N, (ea, eb)
    NC, JTd = o['NC'], o['JT']
    rows = {}                                    # parent element -> {name: [ball per sub-cell]}
    g = ((T.tconst(-arb(2.500806276619802)), Z), (Z, T.tconst(-ONE)), (Z, Z))   # accurate (W, e^{i theta}, IW), W(0) = -DS1
    Phi0 = None
    t0 = time.process_time()
    # ---------------------------------------------------------------- pass 1a: T6 on the base before the window (norms, jets)
    for e in range(ea):
        A_ = std_element(B, e, st)
        if e == 0: Phi0 = A_['Phi'].at(MONE)
        R_ = rot(A_, st)
        g = geo_step(A_, g, Phi0, o, D)
        st, qs = post(C, A_, R_, B.x[e], B.h[e], e < JTn, st)
        rows[e] = jets_whole(qs, NC, JTd)
    snap, st0, g0 = C.snapshot(), dict(st), g
    t1 = time.process_time()
    # ---------------------------------------------------------------- pass 1b: the base from x_ea to ell (no norms): J's, far end
    J1, J2, J1p, J2p = Z, Z, Z, Z
    il = T.tpow(B.ell, -1, D); Lp = B.elld
    tht_far = ZERO
    def win_add(A_, R_, xe, h, inwin):
        nonlocal J1, J2, J1p, J2p, g, tht_far
        g = geo_step(A_, g, Phi0, o, D)
        tht = (A_['rt'] + A_['Phi'] - BP.from_tp(Phi0)).ttrunc(o['DJ'])
        tht_far = max(tht_far, S6.up(tht.sup()))
        if not inwin: return
        DJ = o['DJ']
        Bu, Bb, u = window_weights(B, xe, h)
        Bu, Bb, u = Bu.ttrunc(DJ), Bb.ttrunc(DJ), u.ttrunc(DJ)
        hl = T.tmul(T.ttrunc(A_['hh'], DJ), T.ttrunc(il, DJ), DJ)
        cs, sn = R_['Pg'].re, R_['Pg'].im
        Lam = tht + u.mul(A_['kap'].ttrunc(DJ), DJ).tmul(T.ttrunc(Lp, DJ), DJ)
        I = lambda f: f.int_m1().tmul(hl, DJ).at(ONE)
        J2 = J2 + I(cs.mul(Bu, DJ)); J1 = J1 + I(cs.mul(Bb, DJ))
        sL = sn.mul(Lam, DJ)
        J2p = J2p - I(sL.mul(Bu, DJ)); J1p = J1p - I(sL.mul(Bb, DJ))
    s1 = dict(st)
    bst = {}                                     # base state at the start of every piece / element from x_ea on (for pass 2)
    def split_base(e, X, inwin_k, tag):
        """parent e cut at X: the base residual and derivatives are composed from the parent (T6's rho_kappa, no cancellation)"""
        nonlocal s1
        K, rk, rkx = parent_derivs(B, e, s1)
        pcs, ss = pieces(B, e, X, D)
        comp = []
        for k, (xe, h, a, w) in enumerate(pcs):
            Kp = [compose(f, a, w, D) for f in K]; rkp = compose(rk, a, w, D); rkxp = compose(rkx, a, w, D)
            comp.append((Kp, rkp, rkxp))
            bst[(tag, k)] = dict(s1)
            A_ = element_c(Kp, None, None, rkp, rkxp, xe, h, B.c, s1, s1['Qe'], D); R_ = rot(A_, s1)
            win_add(A_, R_, xe, h, k == inwin_k)
            s1, _ = post(C, A_, R_, xe, h, False, s1, norms=False)
            if tag == 'eb' and k == 0: yield_X2.append(g)
        return ss, comp
    yield_X2 = []
    ssa, compa = split_base(ea, X1, 1, 'ea')
    for e in range(ea + 1, eb):
        bst[e] = dict(s1)
        A_ = std_element(B, e, s1); R_ = rot(A_, s1)
        win_add(A_, R_, B.x[e], B.h[e], True)
        s1, _ = post(C, A_, R_, B.x[e], B.h[e], False, s1, norms=False)
    ssb, compb = split_base(eb, X2, 0, 'eb')
    gX2 = yield_X2[0]
    for e in range(eb + 1, N):
        A_ = std_element(B, e, s1); R_ = rot(A_, s1)
        win_add(A_, R_, B.x[e], B.h[e], False)
        if e == N - 1: PhiN, kapN, kxN, kxxN = A_['Phi'].at(ONE), A_['kap'].at(ONE), A_['kx'].at(ONE), A_['kxx'].at(ONE)
        s1, _ = post(C, A_, R_, B.x[e], B.h[e], False, s1, norms=False)
    t2 = time.process_time()
    # ---------------------------------------------------------------- far end of the base; a1, a2 and their tau-derivatives
    cm = lambda z, w: (T.tmul(z[0], w[0], D) - T.tmul(z[1], w[1], D), T.tmul(z[0], w[1], D) + T.tmul(z[1], w[0], D))
    (WN, PN, IWN), (WX2, PX2, IWX2) = g, gX2                              # accurate geometry (geo_step)
    L = B.ell
    J1 = J1 + T.tmul(WN[0] - WX2[0], il, D)
    # tail u in [0.8, 1] (Bbar = 1): d/dtau int_0.8^1 cos theta~(ell u) du = -(1/ell) Im int_X2^ell theta~_tau e^{i theta~} dx
    #   - (ell'/ell^2) int_X2^ell x kap~ sin theta~ dx,  int x kap sin theta dx = [-x cos theta] + int cos theta dx  (kap = theta_x)
    J1p = (J1p - T.tmul(IWN[1] - IWX2[1], il, D)
           - T.tmul(T.tmul(Lp, T.tmul(il, il, D), D), -T.tmul(L, PN[0], D) + T.tmul(X2, PX2[0], D) + (WN[0] - WX2[0]), D))
    ht = WN[1]
    htp_kin = IWN[0] + T.tmul(Lp, PN[1], D)                              # kinematic: Re int theta_tau e^{i theta} + ell' sin theta(ell)
    # the same rate from T6's velocity identity W_tau|x = e^{i theta}(varpi + i V~), varpi = w* + A, V~ = V* + r (T6_REPORT 1):
    # h~' = Im[e^{i theta(ell)} (varpi + ell' + i V~)(ell)] with varpi + ell' = -k k_x + eps_ell + c q + b t1 + A at ell (no O(ell) terms)
    GN = (T.tmul(PN[0], WN[0], D) + T.tmul(PN[1], WN[1], D), T.tmul(PN[0], WN[1], D) - T.tmul(PN[1], WN[0], D))   # conj(P) W = q + i p
    epsl = Lp - (-s1['Qe'] + T.tmul(B.c, L, D))
    Dp = -T.tmul(kapN, kxN, D) + epsl + T.tmul(B.c, GN[0], D) + T.tmul(B.b, PN[0], D) + s1['ze'][1]
    Vt = -kxxN + T.tmul(B.c, GN[1], D) - T.tmul(B.b, PN[1], D) + s1['ze'][0]          # n = -Im P
    htp = T.tmul(PN[1], Dp, D) + T.tmul(PN[0], Vt, D)
    idgap = float(T.tsup_sub(htp - htp_kin, -1, 1, 16))                  # identity check (both enclose the same function)
    delta = s1['the'] + T.tconst(PI32)
    a1 = -delta
    a1p = -(s1['rte'] + PhiN - Phi0 + T.tmul(Lp, kapN, D))
    iJ2 = tinv(J2, D)
    rL, rLp = trng(L), trng(Lp)
    Lmin, Lmax, Lpmax = rL.lower(), rL.upper(), rLp.abs_upper()
    J2min = trng(J2).abs_lower()
    a1s, a1ps = T.tsup(a1), T.tsup(a1p)
    m = -T.tmul(ht + T.tmul(T.tmul(a1, L, D), J1, D), iJ2, D)
    ms = T.tsup(m)
    win = Lmax * 2 / 5                                                      # ell - X1 <= 0.4 ell
    eps = ZERO
    for _ in range(8):                                                      # self-map of the ball m + [-eps, eps]
        dth = a1s + BMAX * (ms + eps) / Lmin
        Nbar = (dth * dth * win / 2).upper()
        new = (Nbar / J2min).upper()
        if new <= eps: break
        eps = (new * 2).upper()
    else: raise RuntimeError('a2: no self-map')
    assert new <= eps
    lip = (dth / J2min).upper()
    assert lip < 1, f'a2: no contraction ({lip})'
    a2 = tball(m, eps)
    a2s = T.tsup(a2)
    dJ = T.tmul(T.tmul(a1p, L, D), J1, D) + T.tmul(T.tmul(a1, Lp, D), J1, D) + T.tmul(T.tmul(a1, L, D), J1p, D)
    num = htp + dJ + T.tmul(a2, J2p, D)
    dtht = a1ps + a1s * BMAX * Lpmax / Lmin + a2s * Lpmax * (BMAX + B1MAX) / Lmin ** 2       # |d_tau Dtheta| at fixed a2
    Nt = (Lpmax * a1s * a1s / 2 + win * (tht_far * dth * dth / 2 + dth * dtht)).upper()
    a2p = -T.tmul(tball(num, Nt), tinv(tball(J2, dth.upper()), D), D)
    t3 = time.process_time()
    # ---------------------------------------------------------------- pass 2: the corrected base from x_ea to ell (norms, jets)
    C.restore(snap); st = dict(st0); g2 = g0                 # g2: accurate geometry of the corrected base (diagnostic)
    bmp = Bump(B, a1, a1p, a2, a2p) if not o.get('zero_bump') else Bump(B, Z, Z, Z, Z)
    def split_corr(e, X, comp, ss, bump_k, tag):
        nonlocal st, g2
        pcs, _ = pieces(B, e, X, D)
        qq = []
        for k, (xe, h, a, w) in enumerate(pcs):
            Kp, rkp, rkxp = comp[k]
            dK, dKT = bmp.on(xe, h) if k == bump_k else (None, None)
            A_ = element_c(Kp, dK, dKT, rkp, rkxp, xe, h, B.c, st, bst[(tag, k)]['Qe'], D); R_ = rot(A_, st)
            g2 = geo_step(A_, g2, Phi0, o, D)
            st, qs = post(C, A_, R_, xe, h, False, st); qq.append(qs)
        rows[e] = jets_split(qq, ss, NC, JTd)
    split_corr(ea, X1, compa, ssa, 1, 'ea')
    for e in range(ea + 1, eb):
        K, rk, rkx = parent_derivs(B, e, bst[e])
        dK, dKT = bmp.on(B.x[e], B.h[e])
        A_ = element_c(K, dK, dKT, rk, rkx, B.x[e], B.h[e], B.c, st, bst[e]['Qe'], D); R_ = rot(A_, st)
        g2 = geo_step(A_, g2, Phi0, o, D)
        st, qs = post(C, A_, R_, B.x[e], B.h[e], False, st)
        rows[e] = jets_whole(qs, NC, JTd)
    split_corr(eb, X2, compb, ssb, 0, 'eb')
    for e in range(eb + 1, N):
        A_ = std_element(B, e, st); R_ = rot(A_, st)
        g2 = geo_step(A_, g2, Phi0, o, D)
        st, qs = post(C, A_, R_, B.x[e], B.h[e], False, st)
        rows[e] = jets_whole(qs, NC, JTd)
    t4 = time.process_time()
    sub = lambda p: float(T.tsup_sub(p, -1, 1, 16))
    se = dict(exact_end=dict(imW=sub(g2[0][1]), turn=sub(st['the'] + T.tconst(PI32)),
                             note='corrected base: sup |Im W^(ell)| (accurate geometry) and sup |theta^(ell) + 3 pi/2| (enclosures of 0)'),
              sharp=dict(a1=sub(a1), a1p=sub(a1p), a2=sub(a2), a2p=sub(a2p), h=sub(ht), hp=sub(htp)),
              rate_check=dict(identity_minus_kinematic=idgap, rad_identity=float(htp.coeffs()[0].rad()), rad_kinematic=float(htp_kin.coeffs()[0].rad())),
              ea=ea, eb=eb, a1=float(a1s), a1p=float(a1ps), a2=float(a2s), a2p=float(T.tsup(a2p)), eps_a2=float(eps),
              lip=float(lip), Nbar=float(Nbar), Nt=float(Nt), J2=[float(trng(J2).lower()), float(trng(J2).upper())],
              J1=[float(trng(J1).lower()), float(trng(J1).upper())], h=float(T.tsup(ht)), hp=float(T.tsup(htp)),
              cpu=dict(pass1a=t1 - t0, pass1b=t2 - t1, a2=t3 - t2, pass2=t4 - t3))
    return C, st, rows, se


def finish(B, o, parts):
    """norms and scalars of the slab from its parts [(B_k, C_k, st_k)] (one, or two sub-slabs)"""
    D, J, MT = o['D'], o['J'], o['MT']
    norms, norms_int = {}, {}
    sc = {}
    rng_ = {}
    for (Bk, C, st) in parts:
        acc = C.acc; S1N = st['S1e']
        names = sorted(set(k[0] for k in acc.d))
        for base_nm in sorted(set(nm.split('|')[0] for nm in names if '|' in nm)):
            for m in range(MT):
                P_ = acc.d.get((base_nm + '|P', m)); M_ = acc.d.get((base_nm + '|M', m)); L_ = acc.d.get((base_nm + '|L', m))
                if P_ is None: continue
                tot = P_ - T.tmul(M_, S1N, 2 * J) * 2 + T.tmul(L_, T.tmul(S1N, S1N, 2 * J), 2 * J)
                acc.add(base_nm, m, tot)
        for nm in sorted(set(k[0] for k in acc.d if '|' not in k[0])):
            best = ZERO; integ = ZERO
            for m in range(MT):
                p = acc.d.get((nm, m))
                if p is None: continue
                a, b = C.tsub[m]
                v = T.tsup_sub(p, a, b, o['KSUB'])
                if v > best: best = v
                integ = integ + S6.int_upper(p, a, b)
            v = float(S6.up(best.sqrt())) if best > 0 else 0.0
            norms[nm] = max(norms.get(nm, 0.0), v)
            iv_ = S6.up(integ * Bk.dt / 2)
            norms_int[nm] = norms_int.get(nm, ZERO) + iv_
        n = Bk.n
        a_, c_, b_ = Bk.a, Bk.c, Bk.b
        ell, elld = Bk.ell, Bk.elld
        QN = st['Qe']
        lnl = Bk.Y[n + 1]
        lam = T.texp(lnl, D)
        eps_ell = elld - (-QN + T.tmul(c_, ell, D))
        eps_lam = Bk.Yd[n + 1] - (a_ - T.tconst(arb(1) / 4))
        eps_c = Bk.Yd[n + 2] + T.tmul(b_, lam, D)
        l2 = T.tmul(lam, lam, D)
        eps_t = Bk.Yd[n + 3] - T.tmul(l2, l2, D)
        qN, t1N = st['Ge'][0], st['Pe'][0]
        rN, AN = st['ze'][0], st['ze'][1]
        thN = st['the']
        Dend = T.tmul(c_, qN, D) + T.tmul(b_, t1N, D) + AN + eps_ell
        for nm, p in [('eps_ell', eps_ell), ('eps_lnlam', eps_lam), ('eps_c', eps_c), ('eps_t', eps_t), ('closure_q', qN),
                      ('turning_t1', t1N), ('turn_def', thN + PI32), ('r_end', rN), ('A_end', AN), ('D_end', Dend),
                      ('S1_end', S1N)]:
            sc['sup_' + nm] = max(sc.get('sup_' + nm, 0.0), float(T.tsup_sub(p, -1, 1, 16)))
        for nm, p in [('a', a_), ('b', b_), ('c', c_), ('ell', ell), ('ell_tau', elld), ('lnlam', lnl), ('Q_end', QN)]:
            v = S6.tp_at_ball(p, -1, 1); lo, hi = float(v.lower()), float(v.upper())
            if nm in rng_: rng_[nm] = [min(rng_[nm][0], lo), max(rng_[nm][1], hi)]
            else: rng_[nm] = [lo, hi]
    for nm, v in rng_.items(): sc['range_' + nm] = v
    norms_int = {k: (float(S6.up(v.sqrt())) if v > 0 else 0.0) for k, v in norms_int.items()}
    return norms, norms_int, sc


def crossing(B, k):
    """the unique tc* with x_k(tc*) = 0.6 ell(tc*) (T6 knot k; g' = x_k' - 0.6 ell' of fixed sign), as a ball"""
    D = B.D
    g = B.x[k] - T.tmul(B.ell, T.tconst(U1), D)
    gp = B.dx[k] - T.tmul(B.elld, T.tconst(U1), D)
    rg = trng(gp)
    assert rg.lower() > 0 or rg.upper() < 0, 'crossing: g\' changes sign'
    lo, hi = MONE, ONE
    glo = g(lo)
    for _ in range(200):
        mid = (lo + hi) / 2
        gm = g(mid)
        if gm.contains(0): break
        if (gm > 0) == (glo > 0): lo, glo = mid, gm
        else: hi = mid
    # rigorous: sign change on [lo, hi]
    a_, b_ = g(lo), g(hi)
    assert (a_ > 0 and b_ < 0) or (a_ < 0 and b_ > 0), 'crossing: no certified sign change'
    return lo, hi, (1 if rg.lower() > 0 else -1)


class SubSlabBase(SlabBase):
    """SlabBase on the time piece tc in [t0, t1] (arb balls) of slab i, re-centred: tc = al + be tc'"""
    def __init__(self, i, t0, t1, D=20, poldir=None):
        super().__init__(i, D=D, poldir=poldir)
        al = (t0 + t1) / 2; be = (t1 - t0) / 2
        q = arb_poly([al, be])
        self.Y = [p(q) for p in self.Y]; self.Yd = [p(q) for p in self.Yd]
        self.dt = self.dt * be
        self._mesh_tms(); self._scalars()


CROSS = {11: 92}                                    # slab -> T6 knot that crosses u = 0.6 inside the slab


def run_slab_c(i, opts=None, poldir=None, verbose=False):
    o = S6.opts_default()
    if opts: o.update(opts)
    t_start = time.process_time()
    B = SlabBase(i, D=o['D'], poldir=poldir)
    parts, seinfo, rows_all = [], [], []
    if i in CROSS:
        k = CROSS[i]
        lo, hi, sgn = crossing(B, k)
        tstar = WT.iv(lo, hi)                       # ball containing the crossing
        X1 = T.tmul(B.ell, T.tconst(U1), B.D)
        for (t0, t1, before) in ((MONE, tstar, True), (tstar, ONE, False)):
            Bk = SubSlabBase(i, t0, t1, D=o['D'], poldir=poldir)
            # before the crossing x_k - 0.6 ell has the sign of g(-1); the cut lies in element k if x_k < X1, else in k - 1
            gneg = (sgn > 0) == before                 # g increasing: g < 0 before the crossing
            ea = k if gneg else k - 1
            eb = find_parent(Bk, T.tmul(Bk.ell, T.tconst(U2), Bk.D))
            C, st, rows, se = run_part(Bk, o, forced=(ea, eb), verbose=verbose)
            parts.append((Bk, C, st)); seinfo.append(dict(se, tc=[float(t0.lower()) if isinstance(t0, arb) else -1.0,
                                                                 float(t1.upper()) if isinstance(t1, arb) else 1.0]))
            rows_all.append(rows)
    else:
        C, st, rows, se = run_part(B, o, verbose=verbose)
        parts.append((B, C, st)); seinfo.append(se); rows_all.append(rows)
    norms, norms_int, sc = finish(B, o, parts)
    # jets: T6's element rows, union over the parts
    NC = o['NC']; N = B.N
    names = [nm for nm, _ in rows_all[0][0].items()] if False else list(rows_all[0][0].keys())
    jets = {nm: np.zeros((N, NC, 2)) for nm in names}
    for e in range(N):
        for nm in names:
            for cc in range(NC):
                v = None
                for rows in rows_all:
                    w = rows[e][nm][cc]; v = w if v is None else v.union(w)
                jets[nm][e, cc, 0] = np.nextafter(float(v.lower()), -np.inf)
                jets[nm][e, cc, 1] = np.nextafter(float(v.upper()), np.inf)
    crude = {}
    for (_, C, _) in parts:
        for k, v in C.crude.items(): crude[k] = max(crude.get(k, 0.0), v)
    out = dict(slab=i, tau0=B.tau0, k=B.k, n=B.n, N=N, opts=o, norms=norms, norms_int=norms_int, scalars=sc, t7=crude,
               se_direct=seinfo, diag=dict(maxrad=parts[0][1].maxrad, sdeg=parts[0][1].s_deg, cpu_total=time.process_time() - t_start))
    jets['_names'] = np.array(names)
    jets['_meta'] = np.array([i, B.tau0, float(B.dt.mid()), B.k, N, NC], float)
    return out, jets
