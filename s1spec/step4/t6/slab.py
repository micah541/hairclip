"""T6: rigorous residual enclosure of one slab of the T1 base.  run_slab(i) -> (summary dict, jets dict of numpy arrays).

Norms (all certified upper bounds of the sup over the slab; squares are exact polynomials in tc, see README in T6_REPORT):
  N^2(tc) <= sum_cells wbar_c * (h/2) int_c f(s, tc)^2 ds,   wbar_c = sup of the weight over (cell x time-subinterval),
  the right side is a polynomial in tc (coefficient balls), bounded on each time-subinterval by subdivision."""
import os, math, time
import numpy as np
from flint import arb, arb_poly, ctx
import tp as T
from tp import BP, CBP, strunc
from base import SlabBase, expand
from engine import element_a, element_rot
import weights as WT

ONE, ZERO = arb(1), arb(0)
MONE = -ONE
HALF = arb(1) / 2
CUTS = [arb('10.25'), arb('16.25'), arb('34.25'), WT.X0 + 14, WT.X0 + 22, WT.X0 + 28, WT.X0 + 30, WT.X0 + 29, WT.X0 + 32,
        14 - WT.STIP, 28 - WT.STIP]
XRHO = arb('34.25')
XW = arb('10.25')
XF = arb('16.24')          # below every ramp foot (14 + X0 = 16.249, 14 - STIP = 16.249)
FUN_ALL = ('S2', 'r', 'S1', 'r1', 'r2', 'r3')


def opts_default():
    return dict(D=20, J=12, DTE=8, DJ=8, DSmax=30, DXT=0.25, NCO=2, MT=4, KSUB=8, NC=4, JT=4, strtol=1e-30)


def up(x): return x.upper()
def fl(x): return float(x.mid())


def tp_at_ball(p, lo, hi):
    """ball enclosing p over tc in [lo, hi]"""
    return p(WT.iv(arb(lo), arb(hi)))


def tp_range_mid(p):
    c = p.coeffs()
    return [float(v.mid()) for v in c[:3]]


def int_upper(p, a, b):
    """exact upper bound of int_a^b (mid_p(t) + sum_k rad_k |t|^k) dt, i.e. of the integral of every function enclosed
    pointwise by the ball polynomial p on [a, b]"""
    a, b = arb(a), arb(b)
    cs = p.coeffs()
    P = arb_poly([arb(c.mid()) for c in cs]).integral()
    tot = P(b) - P(a)
    for k, c in enumerate(cs):
        r = c.rad()
        if r == 0: continue
        if a >= 0: Jk = (b ** (k + 1) - a ** (k + 1)) / (k + 1)
        elif b <= 0: Jk = ((-a) ** (k + 1) - (-b) ** (k + 1)) / (k + 1)
        else: Jk = ((-a) ** (k + 1) + b ** (k + 1)) / (k + 1)
        tot = tot + r * Jk.upper()
    return tot.upper()


class Acc:
    """accumulators of weighted squared-norm polynomials per (name, time-subinterval)"""
    def __init__(self, MT):
        self.MT = MT; self.d = {}

    def add(self, name, m, tp):
        key = (name, m)
        self.d[key] = tp if key not in self.d else self.d[key] + tp


def run_slab(i, opts=None, poldir=None, verbose=False):
    o = opts_default()
    if opts: o.update(opts)
    t_start = time.process_time()
    B = SlabBase(i, D=o['D'], poldir=poldir)
    D, J, MT = o['D'], o['J'], o['MT']
    N = B.N; JT = B.mesh.JT
    Z = T.tconst(ZERO)
    st = dict(Qe=Z, rte=Z, the=T.tconst(-arb.pi() / 2), ze=(Z, Z), Ge=(Z, T.tconst(-arb(2.500806276619802))),
              Pe=(Z, T.tconst(-ONE)), S1e=Z, S2e=Z)
    acc = Acc(MT)
    tsub = [(-1 + 2 * m / MT, -1 + 2 * (m + 1) / MT) for m in range(MT)]
    ellr = [tp_at_ball(B.ell, a, b) for a, b in tsub]
    llr = [(WT.lo_(e - WT.X0), WT.up(e - WT.X0)) for e in ellr]
    ell_all = tp_at_ball(B.ell, -1, 1); ll_all = (WT.lo_(ell_all - WT.X0), WT.up(ell_all - WT.X0))
    jets = {}; JQ = None
    crude = dict(supr_tip=0.0, supA_tip=0.0, supr1_tip=0.0, supA1_tip=0.0, supArho_tip=0.0,
                 supA_far=0.0, supAx_far=0.0, supA1_far=0.0, supr_far=0.0)
    maxrad = dict(rk=0.0, r=0.0, S2=0.0, S1=0.0)
    s_deg = dict(r=0, r3=0, E=0)
    tim = dict(a=0.0, rot=0.0, der=0.0, norm=0.0, jet=0.0)
    for e in range(N):
        t0 = time.process_time()
        A_ = element_a(B, e, st['Qe'], st['rte'], st['the'])
        t1 = time.process_time()
        R_ = element_rot(B, e, A_, st['ze'], st['Ge'], st['Pe'], DS=o['DSmax'], DTE=o['DTE'], DJ=o['DJ'])
        t2 = time.process_time()
        hh = A_['hh']; f2 = A_['f2']
        z = R_['z']
        r = strunc(z.re, o['strtol']); Aim = strunc(z.im, o['strtol'])
        rt = strunc(A_['rt'], o['strtol'])
        kap, kx, kxx, k2 = A_['kap'], A_['kx'], A_['kxx'], A_['k2']
        rk = A_['rk']
        rkx = rk.ds().tmul(f2, D)
        # derivatives of r at tc-degree J (the norms use degree J; higher degrees are folded)
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
        t3 = time.process_time()
        # ---------------------------------------------------------------- geometry of the element
        tipel = e < JT              # fixed knots
        xe, h = B.x[e], B.h[e]
        xe0, xe1 = tp_at_ball(xe, -1, 1), tp_at_ball(B.x[e + 1], -1, 1)
        xlo_el, xhi_el = WT.lo_(xe0), WT.up(xe1)
        # cells in s
        if tipel:
            xa, xb = xe(ZERO), B.x[e + 1](ZERO)
            pts = [xa] + [c for c in CUTS if c > xa and c < xb] + [xb]
            sb = [MONE]
            for k in range(1, len(pts)):
                seg = pts[k] - pts[k - 1]
                n = max(1, int(math.ceil(float(seg.mid()) / o['DXT'])))
                for q in range(1, n + 1):
                    xq = pts[k - 1] + seg * q / n
                    sq = arb(((2 * (xq - xa) / (xb - xa)) - 1).mid()) if (k < len(pts) - 1 or q < n) else ONE
                    sb.append(sq)
            sb[-1] = ONE
        else:
            sb = [MONE + 2 * arb(q) / o['NCO'] for q in range(o['NCO'] + 1)]
        ncell = len(sb) - 1
        # functions needed on this element
        need = []
        if xlo_el < XRHO: need += ['S2', 'r']
        if xhi_el > XW: need += ['S1', 'r1', 'r3']
        if xhi_el > XF: need += ['r', 'S1', 'r1', 'r2', 'r3']
        need = [f for f in FUN_ALL if f in need]
        F = dict(S2=S2, r=r, S1=S1, r1=r1, r2=r2, r3=r3)
        # per-cell s-integrals of f^2 (and of S1 for the j=0 shift), times h/2
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
        # ---------------------------------------------------------------- weights per cell and time-subinterval
        msubs = range(MT) if not tipel else [None]
        for c in range(ncell):
            for m in msubs:
                if tipel:
                    xa_c = xe(ZERO) + h(ZERO) * HALF * (1 + sb[c]); xb_c = xe(ZERO) + h(ZERO) * HALF * (1 + sb[c + 1])
                    xa_c, xb_c = WT.lo_(xa_c), WT.up(xb_c); lla, llb = ll_all
                else:
                    ta, tb = tsub[m]
                    xb_ = tp_at_ball(xe, ta, tb); hb_ = tp_at_ball(h, ta, tb)
                    xa_c = WT.lo_(xb_ + hb_ * HALF * (1 + sb[c])); xb_c = WT.up(xb_ + hb_ * HALF * (1 + sb[c + 1]))
                    lla, llb = llr[m]
                ms = range(MT) if tipel else [m]
                ya, yb = xa_c - WT.X0, xb_c - WT.X0
                items = []
                if 'S2' in need and xa_c < XRHO:
                    wr = WT.rho_sup(xa_c, xb_c)
                    items += [('R2', 'S2', wr), ('Rs', 'r', wr), ('Rs_chi', 'r', WT.rho_chi2_sup(xa_c, xb_c))]
                if xb_c > XW and 'S1' in need:
                    ww = WT.Wzeta2_sup(xa_c, xb_c)
                    if ww > 0: items += [('R1', 'S1', ww), ('RIs', 'r1', ww), ('RIs3', 'r3', ww)]
                if xb_c > XF and 'r2' in need:
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
                            acc.add(name + '|P', mm, sq['S1'][c] * w); acc.add(name + '|M', mm, M1[c] * w)
                            acc.add(name + '|L', mm, Lc[c] * w)
                    else:
                        for mm in ms: acc.add(name, mm, sq[f][c] * w)
        t4 = time.process_time()
        # ---------------------------------------------------------------- crude pointwise sups (T7 corrections)
        sr, sA, sr1 = up(r.sup()), up(Aim.sup()), up(r1.sup())
        sA1 = up(A1.sup())
        if xlo_el < XRHO:
            crude['supr_tip'] = max(crude['supr_tip'], float(sr)); crude['supA_tip'] = max(crude['supA_tip'], float(sA))
            crude['supr1_tip'] = max(crude['supr1_tip'], float(sr1)); crude['supA1_tip'] = max(crude['supA1_tip'], float(sA1))
            rr = up(xhi_el / 3 * (1 + xlo_el ** 2) ** (-arb(1) / 3))          # |rho'/rho| <= x/3 (1+x^2)^(-1/3)
            crude['supArho_tip'] = max(crude['supArho_tip'], float(up(sA * rr)))
        if xhi_el > XW:
            crude['supA_far'] = max(crude['supA_far'], float(sA)); crude['supA1_far'] = max(crude['supA1_far'], float(sA1))
            crude['supr_far'] = max(crude['supr_far'], float(sr))
            crude['supAx_far'] = max(crude['supAx_far'], float(up(sA / max(xlo_el, XW))))
        maxrad['rk'] = max(maxrad['rk'], max(float(cf.rad()) for cf in rk.P[0].coeffs()))
        maxrad['r'] = max(maxrad['r'], max(float(cf.rad()) for cf in r.P[0].coeffs()))
        maxrad['S1'] = max(maxrad['S1'], max(float(cf.rad()) for cf in S1.P[0].coeffs()))
        maxrad['S2'] = max(maxrad['S2'], max(float(cf.rad()) for cf in S2.P[0].coeffs()))
        s_deg['r'] = max(s_deg['r'], r.sdeg()); s_deg['r3'] = max(s_deg['r3'], r3.sdeg()); s_deg['E'] = max(s_deg['E'], R_['info']['DS'])
        # ---------------------------------------------------------------- jets (hulls on NC sub-cells x whole slab)
        JTd = o['JT']
        G, Pg = R_['G'], R_['Pg']
        c_, b_ = B.c, B.b
        x = A_['x']
        k3 = A_['kxxx']; k4 = A_['k4s'].tmul(T.tmul(A_['f2_2'], A_['f2_2'], D), D)
        k5 = A_['k5s'].tmul(T.tmul(A_['f2_3'], A_['f2_2'], D), D)
        p_, q_ = G.im.ttrunc(JTd), G.re.ttrunc(JTd)
        t1_, n_ = Pg.re.ttrunc(JTd), (-Pg.im).ttrunc(JTd)
        cJ, bJ = T.ttrunc(c_, JTd), T.ttrunc(b_, JTd)
        kxJ4, kxxJ4 = kx.ttrunc(JTd), kxx.ttrunc(JTd)
        kkx4 = kkxJ.ttrunc(JTd)
        Q_J = A_['Q'].ttrunc(JTd)
        Vs = (-kxxJ4) + p_.tmul(cJ, JTd) + n_.tmul(bJ, JTd)
        ws = (-kkx4) + Q_J + (q_ - x.ttrunc(JTd)).tmul(cJ, JTd) + t1_.tmul(bJ, JTd)
        rJ, AJ = r.ttrunc(JTd), Aim.ttrunc(JTd)
        b1 = (-(kkx4.scale(arb(6)))) + Q_J - x.ttrunc(JTd).tmul(cJ, JTd) + AJ
        quants = [('x', x), ('k0', kap), ('k1', kx), ('k2', kxx), ('k3', k3), ('k4', k4), ('k5', k5), ('Q', A_['Q']),
                  ('theta', A_['th']), ('p', p_), ('q', q_), ('n', n_), ('t1', t1_), ('V', Vs + rJ), ('varpi', ws + AJ),
                  ('b1', b1), ('Vstar', Vs), ('varpistar', ws), ('r', r), ('A', Aim), ('r1', r1), ('r2', r2), ('r3', r3),
                  ('A1', A1)]
        if JQ is None:
            JQ = [nm for nm, _ in quants]
            for nm in JQ: jets[nm] = np.zeros((N, o['NC'], 2))
        NC = o['NC']
        for nm, f in quants:
            fJ = f.ttrunc(JTd)
            for cc in range(NC):
                sball = WT.iv(MONE + 2 * arb(cc) / NC, MONE + 2 * arb(cc + 1) / NC)
                vals = [pp(sball) for pp in fJ.P]
                rng = vals[0] + arb(0, T.usum(v.abs_upper() for v in vals[1:])) if len(vals) > 1 else vals[0]
                jets[nm][e, cc, 0] = np.nextafter(float(rng.lower()), -np.inf)
                jets[nm][e, cc, 1] = np.nextafter(float(rng.upper()), np.inf)
        t5 = time.process_time()
        tim['a'] += t1 - t0; tim['rot'] += t2 - t1; tim['der'] += t3 - t2; tim['norm'] += t4 - t3; tim['jet'] += t5 - t4
        # ---------------------------------------------------------------- accumulate states
        st['Qe'], st['rte'], st['the'] = A_['Q'].at(ONE), A_['rt'].at(ONE), A_['th'].at(ONE)
        st['ze'] = z.at(ONE); st['Ge'] = G.at(ONE); st['Pe'] = Pg.at(ONE)
        st['S1e'] = S1.at(ONE); st['S2e'] = S2.at(ONE)
        if verbose and e % 20 == 0:
            print(f'  e={e} x={fl(xe0):.1f} cells={ncell} t={time.process_time() - t_start:.1f}s', flush=True)
    # ==================================================================== end of curve: scalars, far-end defects
    out = dict(slab=i, tau0=B.tau0, k=B.k, n=B.n, N=N, opts=o)
    S1N = st['S1e']
    # finish the shifted j=0 norms
    names = sorted(set(k[0] for k in acc.d))
    shifted = sorted(set(nm.split('|')[0] for nm in names if '|' in nm))
    for base_nm in shifted:
        for m in range(MT):
            P_ = acc.d.get((base_nm + '|P', m)); M_ = acc.d.get((base_nm + '|M', m)); L_ = acc.d.get((base_nm + '|L', m))
            if P_ is None: continue
            tot = P_ - T.tmul(M_, S1N, 2 * J) * 2 + T.tmul(L_, T.tmul(S1N, S1N, 2 * J), 2 * J)
            acc.add(base_nm, m, tot)
    norms, norms_int = {}, {}
    for nm in sorted(set(k[0] for k in acc.d if '|' not in k[0])):
        best = ZERO; integ = ZERO
        for m in range(MT):
            p = acc.d.get((nm, m))
            if p is None: continue
            a, b = tsub[m]
            v = T.tsup_sub(p, a, b, o['KSUB'])
            if v > best: best = v
            integ = integ + int_upper(p, a, b)
        norms[nm] = float(up(best.sqrt())) if best > 0 else 0.0
        # (int over the slab of N^2 dtau)^(1/2) <= ((dtau/2) int_{-1}^{1} bound(tc) dtc)^(1/2)
        iv_ = up(integ * B.dt / 2)
        norms_int[nm] = float(up(iv_.sqrt())) if iv_ > 0 else 0.0
    out['norms'] = norms
    out['norms_int'] = norms_int
    # scalars (TPs)
    n = B.n
    a_, c_, b_ = B.a, B.c, B.b
    ell, elld = B.ell, B.elld
    QN = st['Qe']
    lnl = B.Y[n + 1]
    lam = T.texp(lnl, D)
    eps_ell = elld - (-QN + T.tmul(c_, ell, D))
    eps_lam = B.Yd[n + 1] - (a_ - T.tconst(arb(1) / 4))
    eps_c = B.Yd[n + 2] + T.tmul(b_, lam, D)
    l2 = T.tmul(lam, lam, D)
    eps_t = B.Yd[n + 3] - T.tmul(l2, l2, D)
    qN, t1N = st['Ge'][0], st['Pe'][0]
    rN, AN = st['ze'][0], st['ze'][1]
    thN = st['the']
    Dend = T.tmul(c_, qN, D) + T.tmul(b_, t1N, D) + AN + eps_ell
    sc = {}
    for nm, p in [('eps_ell', eps_ell), ('eps_lnlam', eps_lam), ('eps_c', eps_c), ('eps_t', eps_t), ('closure_q', qN),
                  ('turning_t1', t1N), ('turn_def', thN + 3 * arb.pi() / 2), ('r_end', rN), ('A_end', AN), ('D_end', Dend),
                  ('S1_end', S1N)]:
        sc['sup_' + nm] = float(T.tsup_sub(p, -1, 1, 16))
    def rng(p):
        v = tp_at_ball(p, -1, 1); return [float(v.lower()), float(v.upper())]
    for nm, p in [('a', a_), ('b', b_), ('c', c_), ('ell', ell), ('ell_tau', elld), ('lnlam', lnl), ('Q_end', QN)]:
        sc['range_' + nm] = rng(p)
    out['scalars'] = sc
    out['t7'] = crude
    out['diag'] = dict(maxrad=maxrad, sdeg=s_deg, cpu=tim, cpu_total=time.process_time() - t_start)
    jets['_names'] = np.array(JQ)
    jets['_meta'] = np.array([i, B.tau0, float(B.dt.mid()), B.k, N, o['NC']], float)
    return out, jets
