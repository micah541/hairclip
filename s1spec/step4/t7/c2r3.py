"""T7 (C2), rank-3 part, and (C3), per slab.  Rigorous bound of the variation of the rank-3 term of the projected form,
     V_s >= int_slab || d Pi_tau / d tau ||_Z dtau ,   Pi_tau(w) = a_tau(w)^T M_tau^-T l_tau(w)  (w = rho^1/2 u),
Z = {w : <w, e^{-phi/2} phi_j> = 0} (the constraint space, tau-independent).  For a certificate interval with tau_c a slab
boundary, c_R = max(sum of V_s over the slabs right of tau_c, sum over the slabs left of tau_c) bounds ||Pi_tau - Pi_c||_Z.

Method.  d Pi/d tau = sum_i delta_i m_i + sum_i ahat_i l_tau,i, delta = a_tau - M_tau M^-1 a (the part of the change of the
directions that is not a recombination of them; Pi is invariant under g -> T g), m = M^-T l, ahat = M^-1 a, with
M = <chi g, phi>_rho, M_tau = <chi g_tau, phi>_rho.  Its exact operator norm on Z at a time point is the spectral bound of
the 6x6 forms against the projected Gram matrices (cert7.rate_bound): sharp (it reproduces the FEM value).  Between grid
points tau_k (spacing h): ||Pi'(tau)|| <= ||Pi'(tau_k)|| + |tau - tau_k| C2, with C2 >= sup_slab ||d/dtau(Pi')||_Z from a crude
product-of-norms bound on one tc-ball covering the slab (second tau-derivatives of the base; p_tautau needs V~_tau, bounded
through r_tau from T6's residual identity).  Hence int over [tau_k, tau_k+1] <= h max(F_k, F_k+1) + h^2 C2 / 4.
Base: T6's bivariate (s, tc) enclosures (t6engine, copied), all quantities formed polynomial-first in tc.
(C3) at the grid points (exact) and propagated with crude Lipschitz constants from the same ball.
usage: python3 c2r3.py --slabs a-b --out DIR"""
import os, sys, json, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, arb_poly, arb_mat, ctx
import numpy as np
import base7
import r3core as RC
from tm import TM, up
import tp6 as T
from tp6 import BP
import t6engine as EN
ctx.prec = 128
D_T = 20


def bp_dt(bp, k=1):
    """exact tc-derivative (k times) of an exact-in-tc BP"""
    P = list(bp.P)
    for _ in range(k):
        P = [P[m] * m for m in range(1, len(P))] or [arb_poly([0])]
    return BP(P, bp.ex)


class SlabBP:
    """per tip element BPs (s, tc) of slab S (T6 engine chaining); derived BPs formed polynomial-first"""
    def __init__(self, S, xmax):
        B = S.SB; self.B = B; self.S = S; Z = T.tconst(arb(0)); ONE = arb(1); D = B.D
        st = dict(Qe=Z, rte=Z, the=T.tconst(-arb.pi() / 2), ze=(Z, Z), Ge=(Z, T.tconst(-arb(base7.DS1))), Pe=(Z, T.tconst(-ONE)))
        two_dt = 2 / B.dt; self.two_dt = two_dt
        c, b = B.c, B.b; ct = T.tmul(B.c.derivative(), T.tconst(two_dt), D); bt = T.tmul(B.b.derivative(), T.tconst(two_dt), D)
        self.sc = dict(c=c, b=b, ct=ct, bt=bt, ctt=T.tmul(ct.derivative(), T.tconst(two_dt), D), btt=T.tmul(bt.derivative(), T.tconst(two_dt), D))
        self.els = []
        for e in range(B.N):
            if float(B.x[e](arb(0)).lower()) > float(xmax.upper()): break
            A_ = EN.element_a(B, e, st['Qe'], st['rte'], st['the'])
            R_ = EN.element_rot(B, e, A_, st['ze'], st['Ge'], st['Pe'], DS=30, DTE=8, DJ=8)
            st['Qe'], st['rte'], st['the'] = A_['Q'].at(ONE), A_['rt'].at(ONE), A_['th'].at(ONE)
            st['ze'] = R_['z'].at(ONE); st['Ge'] = R_['G'].at(ONE); st['Pe'] = R_['Pg'].at(ONE)
            d = {}
            d['k0'], d['k1'], d['k2'] = A_['kap'], A_['kx'], A_['kxx']
            for k in ('k0', 'k1', 'k2'):
                d[k + 't'] = bp_dt(d[k]).tmul(T.tconst(two_dt), D); d[k + 'tt'] = bp_dt(d[k], 2).tmul(T.tconst(two_dt * two_dt), D)
            d['Q'] = A_['Q']; d['Qt'] = bp_dt(A_['Q']).tmul(T.tconst(two_dt), D); d['Qtt'] = bp_dt(A_['Q'], 2).tmul(T.tconst(two_dt * two_dt), D)
            d['tht'] = bp_dt(A_['th']).tmul(T.tconst(two_dt), D); d['thtt'] = bp_dt(A_['th'], 2).tmul(T.tconst(two_dt * two_dt), D)
            d['rt'] = A_['rt']; d['rtt'] = bp_dt(A_['rt']).tmul(T.tconst(two_dt), D)
            d['thx'] = A_['Phi'] + A_['rt']                 # theta_tau at fixed x (all elements; T6 identity)
            p, q, t1 = R_['G'].im, R_['G'].re, R_['Pg'].re; n = -R_['Pg'].im
            d['p'], d['q'], d['n'], d['t1'], d['r'], d['A'] = p, q, n, t1, R_['z'].re, R_['z'].im
            Vs = (-d['k2']) + p.tmul(c, D) + n.tmul(b, D)
            pt = Vs + d['r'] - d['tht'].mul(q, D)
            nt = -(d['tht'].mul(t1, D))
            Vst = (-d['k2t']) + p.tmul(ct, D) + pt.tmul(c, D) + n.tmul(bt, D) + nt.tmul(b, D)
            d['Vs'], d['pt'], d['nt'], d['Vst'] = Vs, pt, nt, Vst
            if e < B.mesh.JT:                               # second tau-derivatives, polynomial-first (fixed knots: th_tau exact)
                xb = EN.xpoly(B.x[e], B.h[e])
                w_ = (-(d['k0'].mul(d['k1'], D))) + d['Q'] + (q - xb).tmul(c, D) + t1.tmul(b, D) + d['A']
                qt = w_ + d['tht'].mul(p, D)
                ptt0 = Vst - d['thtt'].mul(q, D) - d['tht'].mul(qt, D)              # + r_tau (bounded separately)
                ntt = (-(d['thtt'].mul(t1, D))) - d['tht'].mul(d['tht'], D).mul(n, D)
                Vstt0 = (-d['k2tt']) + p.tmul(self.sc['ctt'], D) + pt.tmul(ct, D).scale(arb(2)) + ptt0.tmul(c, D) \
                        + n.tmul(self.sc['btt'], D) + nt.tmul(bt, D).scale(arb(2)) + ntt.tmul(b, D)
                d['ptt0'], d['ntt'], d['Vstt0'] = ptt0, ntt, Vstt0
            d['x'] = B.x[e](arb(0)); d['h'] = B.h[e](arb(0))
            self.els.append(d)
        self.xk = [float(el['x'].mid()) for el in self.els]

    def at(self, t, keys):
        out = []
        for el in self.els:
            out.append({k: base7.bp_at(el[k], t) for k in keys} | dict(x=el['x'], h=el['h']))
        return out, {k: v(t) for k, v in self.sc.items()}


KEYS1 = ('k0', 'k1', 'k2', 'k0t', 'k1t', 'k2t', 'Q', 'Qt', 'p', 'n', 'Vs', 'pt', 'nt', 'Vst')
KEYS2 = KEYS1 + ('k0tt', 'k1tt', 'k2tt', 'Qtt', 'tht', 'thtt', 'q', 't1', 'r', 'A', 'rt', 'rtt', 'ptt0', 'ntt', 'Vstt0')


def bs_of(k0, k1, k2, Q, x, c):
    return (-(k0 * k0), -2 * k0 * k1, -2 * k1 * k1 - 2 * k0 * k2, -6 * k0 * k1 + Q - x * c,
            -5 * k1 * k1 - 6 * k0 * k2 - arb_poly([c]), arb_poly([c]) - 3 * k1 * k1 - 4 * k0 * k2)


def bt_of(k0, k1, k2, t0, t1, t2, Qt, x, ct):
    return (-2 * k0 * t0, -2 * (k1 * t0 + k0 * t1), -2 * (k2 * t0 + 2 * k1 * t1 + k0 * t2),
            -6 * (t0 * k1 + k0 * t1) + Qt - x * ct, -10 * k1 * t1 - 6 * (t0 * k2 + k0 * t2) - arb_poly([ct]),
            arb_poly([ct]) - 6 * k1 * t1 - 4 * (t0 * k2 + k0 * t2))


def btt_of(k0, k1, k2, t0, t1, t2, s0, s1, s2, Qtt, x, ctt):
    return (-2 * (t0 * t0 + k0 * s0), -2 * (s0 * k1 + 2 * t0 * t1 + k0 * s1),
            -4 * (t1 * t1 + k1 * s1) - 2 * (s0 * k2 + 2 * t0 * t2 + k0 * s2),
            -6 * (s0 * k1 + 2 * t0 * t1 + k0 * s1) + Qtt - x * ctt,
            -10 * (t1 * t1 + k1 * s1) - 6 * (s0 * k2 + 2 * t0 * t2 + k0 * s2) - arb_poly([ctt]),
            arb_poly([ctt]) - 6 * (t1 * t1 + k1 * s1) - 4 * (s0 * k2 + 2 * t0 * t2 + k0 * s2))


def funcs(SB, t, FX, order=1):
    """per coarse element TMs: order 1 -> [a (3), a_tau (3), l (3), l_tau (3)]; order 2 -> + [a_tautau (3), l_tautau (3)]"""
    els, sc = SB.at(t, KEYS1 if order == 1 else KEYS2)
    c, b, ct, bt = sc['c'], sc['b'], sc['ct'], sc['bt']
    out = []; xs = SB.xk; Iacc = arb(0)
    for fx in FX.data:
        a, h = fx['a'], fx['h']
        E = max(0, min(len(xs) - 1, int(np.searchsorted(xs, float(a.mid()) + 1e-12, side='right') - 1)))
        el = els[E]; S_ = arb_poly([2 * (a - el['x']) / el['h'] - 1, 2 * h / el['h']])
        P = {k: el[k](S_) for k in el if k not in ('x', 'h')}
        tmk = lambda k: TM(P[k])
        x = arb_poly([a, h])
        bs = bs_of(P['k0'], P['k1'], P['k2'], P['Q'], x, c)
        bt_ = bt_of(P['k0'], P['k1'], P['k2'], P['k0t'], P['k1t'], P['k2t'], P['Qt'], x, ct)
        phis = fx['phis']; chi0 = fx['chi'][0]
        F = [chi0 * tmk(k) for k in ('p', 'n', 'Vs', 'pt', 'nt', 'Vst')]
        F += [RC.ell_strong(phis[j], bs, fx['d']) for j in range(3)] + [RC.ell_b(phis[j], bt_, fx['d']) for j in range(3)]
        if order == 2:
            ctt, btt = sc['ctt'], sc['btt']
            p, q, n, t1 = tmk('p'), tmk('q'), tmk('n'), tmk('t1'); tht, thtt = tmk('tht'), tmk('thtt')
            pt, nt = tmk('pt'), tmk('nt')
            # r_tau: z = r + iA = e^{i th} int_0^x e^{-i th} rho_theta  =>  |r_tau| <= |th_tau| |z| + int_0^x (|th_tau| |rho_th| + |rho_th,tau|)
            zab = up((tmk('r') * tmk('r') + tmk('A') * tmk('A')).sup()).sqrt()
            loc = (up(tht.sup()) * up(TM(P['rt']).sup()) + up(TM(P['rtt']).sup())) * up(h)
            rtb = up(tht.sup()) * zab + Iacc + loc; Iacc = Iacc + loc
            rT = TM(arb_poly([0]), rtb)
            ptt = tmk('ptt0') + rT; ntt = tmk('ntt'); Vstt = tmk('Vstt0') + c * rT
            btt_ = btt_of(P['k0'], P['k1'], P['k2'], P['k0t'], P['k1t'], P['k2t'], P['k0tt'], P['k1tt'], P['k2tt'], P['Qtt'], x, ctt)
            F = F[0:6] + [chi0 * v for v in (ptt, ntt, Vstt)] + F[6:12] + [RC.ell_b(phis[j], btt_, fx['d']) for j in range(3)]
        out.append(F)
    return out


def gram(FX, F):
    m = len(F[0]) + 3; G = [[arb(0)] * m for _ in range(m)]
    for fx, fs in zip(FX.data, F):
        allf = fs + [fx['phis'][j][0] for j in range(3)]
        w = [fx['e1'] * f for f in allf]
        for i in range(m):
            for j in range(i, m):
                v = (w[i] * allf[j]).integral01() * fx['h']; G[i][j] += v
                if j != i: G[j][i] += v
    return RC.amat(G)


def tr_up(C, G):
    """sqrt(trace(C G C^T)) upper bound (Frobenius-type norm of the functionals C h)"""
    S = C * G * C.transpose()
    return arb(max(float(sum((S[i, i] for i in range(S.nrows())), arb(0)).upper()), 0.0)).sqrt()


def crude(G21):
    """min of the crude bounds with and without the direction rescaling (both valid)"""
    a = crude_one(G21, True); b = crude_one(G21, False)
    if not math.isfinite(float(b['C2'].upper())) or float(a['C2'].upper()) <= float(b['C2'].upper()): out = dict(a)
    else: out = dict(b)
    for k in ('Lip_mbar', 'Lip_Cx', 'Lip_lbar', 'Lip_nV'): out[k] = a[k]
    out['C2_scaled'], out['C2_unscaled'] = a['C2'], b['C2']
    return out


def crude_one(G21, scale):
    """crude sup over the ball of ||d/dtau (Pi')||_Z and Lipschitz constants for (C3); blocks a0 a1 a2 l0 l1 l2 phi
    (direction blocks rescaled by 1/||chi g_i||, which leaves Pi and its derivatives invariant)"""
    G0 = G21; G21 = RC.gscale(G21, 3) if scale else G21
    blk = lambda i, j: arb_mat(3, 3, [G21[3 * i + a, 3 * j + b] for a in range(3) for b in range(3)])
    M, Mt, Mtt = blk(0, 6), blk(1, 6), blk(2, 6); Mi = M.inv()
    Tm = Mt * Mi; Tp = Mtt * Mi - Mt * Mi * Mt * Mi
    I3 = arb_mat(3, 3, [arb(1) if i == j else arb(0) for i in range(3) for j in range(3)])
    def C(*parts):
        out = arb_mat(3, 21)
        for blkid, Mat in parts:
            for i in range(3):
                for j in range(3): out[i, 3 * blkid + j] += Mat[i, j]
        return out
    MiT = Mi.transpose()
    nd = tr_up(C((1, I3), (0, -1 * Tm)), G21); ndp = tr_up(C((2, I3), (0, -1 * Tp), (1, -1 * Tm)), G21)
    nm = tr_up(C((3, MiT)), G21); nmp = tr_up(C((3, -1 * (MiT * Mt.transpose() * MiT)), (4, MiT)), G21)
    na = tr_up(C((0, Mi)), G21); nap = tr_up(C((0, -1 * (Mi * Mt * Mi)), (1, Mi)), G21)
    nl = tr_up(C((4, I3)), G21); nlp = tr_up(C((5, I3)), G21)
    C2 = ndp * nm + nd * nmp + nap * nl + na * nlp
    # Lipschitz data for (C3) (unscaled Grams): normalised M_n = Dg^-1 M Dphi^-1
    G21 = G0; blk = lambda i, j: arb_mat(3, 3, [G21[3 * i + a, 3 * j + b] for a in range(3) for b in range(3)])
    M, Mt = blk(0, 6), blk(1, 6)
    sq_lo = lambda v: arb(max(float(v.lower()), 0.0)).sqrt(); sq_up = lambda v: arb(max(float(v.upper()), 0.0)).sqrt()
    ng = [sq_lo(G21[i, i]) for i in range(3)]; ngt = [sq_up(G21[3 + i, 3 + i]) for i in range(3)]; nphi = [sq_lo(G21[18 + j, 18 + j]) for j in range(3)]
    Mn = arb_mat(3, 3, [M[i, j] / (ng[i] * nphi[j]) for i in range(3) for j in range(3)])
    Mnt = arb_mat(3, 3, [(Mt[i, j] / ng[i] + M[i, j] * ngt[i] / (ng[i] * ng[i])) / nphi[j] for i in range(3) for j in range(3)])  # |entries| bound
    Mni = Mn.inv()
    fro = lambda A: arb(float(sum((abs(A[i, j]) ** 2 for i in range(3) for j in range(3)), arb(0)).upper())).sqrt()
    nphimin = min(nphi, key=lambda v: float(v.lower()))
    Lam = tr_up(C((3, I3)), G21) / nphimin; Lamt = tr_up(C((4, I3)), G21) / nphimin
    LipM = fro(Mni) ** 2 * fro(Mnt)
    return dict(C2=C2, nd=nd, nm=nm, na=na, nl=nl, Lip_mbar=LipM, Lip_Cx=LipM * Lam + fro(Mni) * Lamt, Lip_lbar=Lamt,
                Lip_nV=sq_up(G21[5, 5]))


def run_slab(s, FX, nmin=2, nmax=12, overhead=0.5, verbose=True):
    t0 = time.time(); S = base7.Slab(s); SB = SlabBP(S, RC.XE + 1); tb = time.time() - t0
    G21 = gram(FX, funcs(SB, arb(0, 1), FX, order=2)); cr = crude(G21); tcr = time.time() - t0 - tb
    C2 = cr['C2']; L = S.dtf
    def point(tc):
        G = gram(FX, funcs(SB, arb(tc), FX, order=1)); rt, info = RC.rate_bound(G); c3 = RC.c3_from_gram(G)
        return float(rt.upper()), {k: float(v.upper()) if k != 'nV' else float(v.lower()) for k, v in c3.items()}
    Fm, c3m = point(-1); Fp, c3p = point(1)
    # grid: second-order part L^2 C2 / (4 n) <= overhead * L min(F)
    C2f = float(C2.upper())
    n = int(min(nmax, max(nmin, math.ceil(L * C2f / (4 * overhead * max(min(Fm, Fp), 1e-300))))))
    tcs = [-1 + 2 * k / n for k in range(n + 1)]; Fs = [Fm]; C3s = [c3m]
    for tc in tcs[1:-1]:
        f, c3 = point(tc); Fs.append(f); C3s.append(c3)
    Fs.append(Fp); C3s.append(c3p)
    h = L / n
    V = sum(h * max(Fs[k], Fs[k + 1]) + h * h * C2f / 4 for k in range(n))
    hh = h / 2
    c3 = dict(mbar=max(c['mbar'] for c in C3s) + hh * float(cr['Lip_mbar'].upper()) * 1.0,
              lbar=max(c['lbar'] for c in C3s) + hh * float(cr['Lip_lbar'].upper()),
              Cx=max(c['Cx'] for c in C3s) + hh * float(cr['Lip_Cx'].upper()),
              Cx_full=max(c['Cx_full'] for c in C3s) + hh * float(cr['Lip_Cx'].upper()),
              Cx_i=[max(c[f'Cx{i}'] for c in C3s) + hh * float(cr['Lip_Cx'].upper()) for i in (1, 2, 3)],
              nV_min=min(c['nV'] for c in C3s) - hh * float(cr['Lip_nV'].upper()))
    out = dict(slab=s, tau=[S.tau0, S.tau0 + L], k=S.k0, V_rank3=V, n_grid=n, rate_pts=Fs, C2=C2f,
               crude={k: float(v.upper()) for k, v in cr.items()}, C3=c3, C3_pts=C3s, t=time.time() - t0, t_bp=tb, t_crude=tcr)
    if verbose:
        print(f"slab {s} (k {S.k0:.1f}): V_rank3 <= {V:.3e} (rate {min(Fs):.3e}..{max(Fs):.3e}, C2 {C2f:.2e}, n {n}); "
              f"mbar {c3['mbar']:.3f} lbar {c3['lbar']:.3f} Cx {c3['Cx']:.3f} nV >= {c3['nV_min']:.4e} ({out['t']:.0f}s)", flush=True)
    return out


def fixed_coarse(hc=0.5, h=0.1):
    """the fixed functions on the coarse mesh: from the cache fixed_coarse.json if present (portable), else built with cert7
    (needs step 3's c4fun) and cached"""
    S5 = base7.Slab(base7.slab_range()[0]); xtip = [float(v) for v in S5.SB.mesh.xtip]      # tip knots: the same on every slab
    if os.path.exists(RC.CACHE):
        FC = RC.FixedCache()
        if FC.xtip is None or FC.xtip == xtip: return FC
        print('fixed_coarse.json was built for other tip knots: rebuilding (needs step 3 c4fun)', flush=True)
    import cert7 as C7, certdata as CD
    BL = base7.late_base(CD.XE + 1); Bc = S5.base(arb(-1), CD.XE + 1)
    FX = C7.Fixed(C7.coarse_nodes(Bc, CD.mesh(Bc, h), hc), BL); RC.save_fixed(FX, xtip=xtip); return RC.FixedCache()



def parse_slabs(spec):
    """'a-b' or comma-separated slabs / ranges, e.g. '5-20' or '5,13,21-23'"""
    out = []
    for part in spec.split(','):
        if '-' in part: a_, b_ = [int(v) for v in part.split('-')]; out += list(range(a_, b_ + 1))
        elif part.strip(): out.append(int(part))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--slabs', default='5-5'); ap.add_argument('--out', default=os.path.join(HERE, 'out_c2r3'))
    base7.base_arg(ap)
    A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); base7.set_base(A.base)
    FX = fixed_coarse()
    for s in parse_slabs(A.slabs):
        fn = os.path.join(A.out, 'c2r3_%03d.json' % s)
        if os.path.exists(fn): continue
        o = run_slab(s, FX)
        json.dump(o, open(fn, 'w'), indent=1)
