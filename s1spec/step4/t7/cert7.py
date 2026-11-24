"""T7 production of (C1) + (C2): one certificate per interval I = [tau_lo, tau_hi] (a group of polished slabs), tau_c = a slab
boundary inside I (exact stored row 0 of slab s_c).  theta_I = 0.025 (T5's theta_0 of the closure).

Claim certified for each interval (all tau in I, all even u in H^4_rho supported in [0, XE], XE = 32 - STIP, <u, phi_j>_rho = 0):
     B_tau(u) <= -c_I ||u||_rho^2 - theta_I ||u''||_rho^2 .
Proof = the prototype certificate (cert.py: Hermite-cubic block + Haynsworth inertia + tail + leak) applied at tau_c to
     Q_c(w) - sum_e (mu0_e ||w||_e^2 + mu1_e ||w'||_e^2) - c_R ||w||^2      (w = rho^(1/2) u),
since for tau in I, with Delta = (tau) - (tau_c) (C22, C21, C11 are tau-independent; boundary terms vanish by symmetry):
     Q_tau(w) - Q_c(w) = int [-2 dC20 w'^2 + (dC20'' - dC10' + dC00) w^2] + (Pi_tau - Pi_c)(w),
     dC20 = -db2/2,  dC10 = -(db2 phi' + db1)/2,  dC00 = -db2 c - db1 phi'/2 - db0,   c = phi''/2 + phi'^2/4,
mu1_e >= sup_{e x I} |db2|, mu0_e >= sup_{e x I} |dC20'' - dC10' + dC00| (ball evaluation of the base on 2 tc-balls per slab), and
c_R >= sup_I sup_{w in Z} |Pi_tau(w) - Pi_c(w)| / ||w||^2, Z = the constraint space (Gram matrices of the representers projected
orthogonally to e^{-phi/2} phi_j; 4-term expansion with the ball matrix M_tau).
Also returned per slab (C3): ||M^-1||, lbar, C_x (= sup_Z |M^-T l(u)| / ||u||, normalised directions), n_s = ||chi V_s||_rho.
Residual corrections (T6): V = V_s + r in g_3 and b1 += int_0^x kappa r are NOT included; the certified c_I holds for the
polished base with V_s; the stated sensitivity (report) absorbs T6's bounds.
usage: python3 cert7.py --times i-j --out DIR [--base STORE] [--slabs a-b] [--h 0.1] [--nbis 7]   (index into schedule(); one JSON each)"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, arb_poly, arb_mat, ctx
import base7
import certdata as CD
import cert as C0
from tm import TM, sup01, up
ctx.prec = 128
THETA_I = arb(25) / 1000
C0.THETA = THETA_I
PAIRS = C0.PAIRS


def schedule(slabs=None):
    """the certificate intervals on the selected store (base7.schedule: size 2/4/10/20 by k at the group start); slabs = 'a-b'
    restricts the plan to slabs a..b (default: the whole store; the v1 store used 5-187)"""
    a_, b_ = ([int(v) for v in slabs.split('-')] if slabs else (None, None))
    return base7.schedule(None, a_, b_)


# --------------------------------------------------------------------------------------------- fixed functions per element
class Fixed:
    """caches of the tau-independent functions on a mesh: chi, phi_j jets, phi' jets, e^{-phi/2}, e^{-phi}"""
    def __init__(self, nodes, BL):
        self.nodes = nodes; self.data = []
        geoL = CD.Geo(BL)
        for e in range(len(nodes) - 1):
            a = nodes[e]; h = nodes[e + 1] - a
            GL = geoL.element(CD.t1_index(BL, a), a, h, nder=4)
            d, em, e1 = CD.phi_tms(a, h)
            chi = [TM(c) for c in CD.chi_jets(a, h)]
            pe = CD.psi_index(a)
            ps = CD.psi3_piece(pe, a, h) if pe is not None else CD.psi3_ext(a, h)
            phis = [[sum((CD.comb(m, kk) * chi[kk] * bs[m - kk] for kk in range(m + 1)), TM(arb_poly([0]))) for m in range(5)]
                    for bs in (GL['P'], GL['N'], ps)]
            cphi = 0.5 * d[1] + 0.25 * (d[0] * d[0])
            self.data.append(dict(a=a, h=h, d=d, em=em, e1=e1, chi=chi, phis=phis, cphi=cphi))


def ell_b(f, b, d):
    """b-linear part of rho^-1 L^dag(rho f): (b2 f)''-type terms - (b1 f)'-type terms + b0 f (rho ratios r1, r2)"""
    p1, p2, p3, p4 = d
    r1 = -1 * p1; r2 = p1 * p1 - p2
    f0, f1, f2 = f[0], f[1], f[2]; b2, b2p, b2pp, b1, b1p, b0 = [TM(v) for v in b]
    h0 = b2 * f0; h1 = b2p * f0 + b2 * f1; h2 = b2pp * f0 + 2 * (b2p * f1) + b2 * f2
    t2 = h2 + 2 * (r1 * h1) + r2 * h0
    t1 = (b1p * f0 + b1 * f1) + r1 * (b1 * f0)
    return t2 - t1 + b0 * f0


def r3_elem(B, geo, fx, e):
    """per element: directions g = (p, n, V_s) and the b-polys (for ell) of base B"""
    a, h = fx['a'], fx['h']; E = CD.t1_index(B, a)
    G = geo.element(E, a, h); L = G['L']
    bs = CD.coeff_polys(L, B, a, h)
    g = [G['p'], G['n'], TM(-1 * L['K'][2]) + B.c * G['p'] + B.b * G['n']]
    return g, bs, G


Z3 = lambda: [[arb(0)] * 3 for _ in range(3)]


def r3_center(Bc, FX):
    """Gram data at tau_c; all integrals int e^{-phi} (.)(.) over [0, XE]"""
    gc = CD.Geo(Bc); C = dict(GG=Z3(), GL=Z3(), GP=Z3(), M=Z3(), PL=Z3())
    for e, fx in enumerate(FX.data):
        h = fx['h']; e1 = fx['e1']; chi0 = fx['chi'][0]; phis = fx['phis']
        g, bs, _ = r3_elem(Bc, gc, fx, e)
        ells = [CD.ell_strong(phis[j], bs, fx['d']) for j in range(3)]
        cg = [chi0 * gi for gi in g]
        for i in range(3):
            for j in range(3):
                C['GG'][i][j] += (e1 * cg[i] * cg[j]).integral01() * h
                C['GL'][i][j] += (e1 * ells[i] * ells[j]).integral01() * h
                C['GP'][i][j] += (e1 * phis[i][0] * phis[j][0]).integral01() * h
                C['M'][i][j] += (e1 * cg[i] * phis[j][0]).integral01() * h
                C['PL'][i][j] += (e1 * phis[i][0] * ells[j]).integral01() * h
    return C


def frame_pq(K, p0, q0, h, tol=1e-32, itmax=80):
    """p' = -k q, q' = 1 + k p on [a, a+h] (x = a + h s), Picard iteration in TMs; rigorous remainder by contraction
    (L = h sup|k| < 1): the true (p, q) lies within ||Phi(y_n) - y_n|| / (1 - L) of y_n."""
    L = h * K.sup()
    if not float(L.upper()) < 0.9: raise RuntimeError('frame_pq: element too long')
    p = TM(arb_poly([p0])); q = TM(arb_poly([q0])); one = TM(arb_poly([1]))
    for it in range(itmax):
        pn = (-1 * (K * q)).antideriv(h) + p0; qn = (one + K * p).antideriv(h) + q0
        d = max(float((pn - p).sup().upper()), float((qn - q).sup().upper()))
        p, q = pn, qn
        if d < tol: break
    pn = (-1 * (K * q)).antideriv(h) + p0; qn = (one + K * p).antideriv(h) + q0
    d = max((pn - p).sup(), (qn - q).sup(), key=lambda v: float(v.upper()))
    R = up(d / (1 - L))
    return TM(pn.p, pn.r + R), TM(qn.p, qn.r + R)


def r3_funcs(B, FX):
    """per coarse element the 12 functions [chi g (3), chi g_tau (3), ell (3), ell_tau (3)] of base B (point or ball), sharp
    geometry: p, q from the Frenet ODE (no rotation of accumulated positions), n = -sin theta, t1 = cos theta;
    kinematics: theta_tau = int k_tau (tip angle fixed), V = V_s (+ r, T6), p_tau = V_s - theta_tau q, n_tau = -theta_tau t1,
    V_s,tau = -k_tau'' + c_tau p + c p_tau + b_tau n + b n_tau; b_tau from k_tau, Q_tau = int 2 k' k_tau', c_tau, b_tau."""
    out = []; p0 = -arb(base7.DS1); q0 = arb(0); tht0 = arb(0); Qt0 = arb(0)
    for e, fx in enumerate(FX.data):
        a, h = fx['a'], fx['h']; E = CD.t1_index(B, a); L = B.local(E, a, h)
        cs, sn = CD.trig_tms(L['th'], a, h); K0t = TM(L['K'][0])
        p, q = frame_pq(K0t, p0, q0, h)
        n = -1 * sn; t1 = cs
        bs = CD.coeff_polys(L, B, a, h)
        Vs = TM(-1 * L['K'][2]) + B.c * p + B.b * n
        tht = TM(L['Kt'][0].integral() * h + tht0)
        pt = Vs - tht * q; nt = -1 * (tht * t1)
        Vst = TM(-1 * L['Kt'][2]) + B.ct * p + B.c * pt + B.bt * n + B.b * nt
        k0, k1, k2 = L['K'][0], L['K'][1], L['K'][2]; t0_, t1_, t2_ = L['Kt'][0], L['Kt'][1], L['Kt'][2]
        Qt = (2 * k1 * t1_).integral() * h + Qt0; x = arb_poly([a, h])
        bt = (-2 * k0 * t0_, -2 * (k1 * t0_ + k0 * t1_), -2 * (k2 * t0_ + 2 * k1 * t1_ + k0 * t2_),
              -6 * (t0_ * k1 + k0 * t1_) + Qt - x * B.ct, -10 * k1 * t1_ - 6 * (t0_ * k2 + k0 * t2_) - arb_poly([B.ct]),
              arb_poly([B.ct]) - 6 * k1 * t1_ - 4 * (t0_ * k2 + k0 * t2_))
        phis = fx['phis']; chi0 = fx['chi'][0]
        ells = [CD.ell_strong(phis[j], bs, fx['d']) for j in range(3)]
        lts = [ell_b(phis[j], bt, fx['d']) for j in range(3)]
        out.append([chi0 * v for v in (p, n, Vs)] + [chi0 * v for v in (pt, nt, Vst)] + ells + lts)
        p0, q0 = p.p(arb(1)) + arb(0, p.r.upper()), q.p(arb(1)) + arb(0, q.r.upper())
        tht0 = tht0 + L['Kt'][0].integral()(arb(1)) * h; Qt0 = Qt(arb(1))
    return out


class TipBP:
    """T6 bivariate (s, tc) enclosures of the base on the tip elements (x_e < XE) of slab S (t6engine element_a/element_rot,
    chained exactly as T6's run_slab)"""
    def __init__(self, S, xmax):
        import tp6 as T, t6engine as EN
        B = S.SB; self.B = B; self.T = T; Z = T.tconst(arb(0)); ONE = arb(1)
        st = dict(Qe=Z, rte=Z, the=T.tconst(-arb.pi() / 2), ze=(Z, Z), Ge=(Z, T.tconst(-arb(base7.DS1))), Pe=(Z, T.tconst(-ONE)))
        self.els = []
        for e in range(B.N):
            if float(B.x[e](arb(0)).lower()) > float(xmax.upper()): break
            A_ = EN.element_a(B, e, st['Qe'], st['rte'], st['the'])
            R_ = EN.element_rot(B, e, A_, st['ze'], st['Ge'], st['Pe'], DS=30, DTE=8, DJ=8)
            st['Qe'], st['rte'], st['the'] = A_['Q'].at(ONE), A_['rt'].at(ONE), A_['th'].at(ONE)
            st['ze'] = R_['z'].at(ONE); st['Ge'] = R_['G'].at(ONE); st['Pe'] = R_['Pg'].at(ONE)
            f2 = A_['f2']
            ktx = A_['kt'].ds().tmul(f2, B.D); ktxx = ktx.ds().tmul(f2, B.D)
            self.els.append(dict(x=B.x[e](arb(0)), h=B.h[e](arb(0)), kap=A_['kap'], kx=A_['kx'], kxx=A_['kxx'], Q=A_['Q'],
                                 Phi=A_['Phi'], rt=A_['rt'], kt=A_['kt'], ktx=ktx, ktxx=ktxx,
                                 q=R_['G'].re, p=R_['G'].im, t1=R_['Pg'].re, sn=R_['Pg'].im, r=R_['z'].re))
        self.xk = [el['x'] for el in self.els]

    def at(self, t):
        """per element s-polys at tc = t (ball): dict of arb_polys in s in [-1, 1]; scalars"""
        B = self.B; two_dt = 2 / B.dt
        sc = dict(c=B.c(t), b=B.b(t), ct=-(B.a.derivative()(t) * two_dt), bt=B.b.derivative()(t) * two_dt)
        out = []
        for el in self.els:
            out.append({k: base7.bp_at(v, t) for k, v in el.items() if k not in ('x', 'h')} | dict(x=el['x'], h=el['h']))
        return out, sc


def r3_funcs_bp(TB_, t, FX):
    """the 12 functions [chi (p, n, V_s), chi (p_tau, n_tau, V_s,tau), ell_1..3, ell_tau,1..3] at tc-ball t from T6's sharp
    bivariate base: p, q (tangent-frame position), n = -sin theta, t1 = cos theta, kappa jets exact; kinematics
    theta_tau = Phi + rho_theta, V~ = V_s + r (T6), p_tau = V~ - theta_tau q, n_tau = -theta_tau t1,
    V_s,tau = -k_tau'' + c_tau p + c p_tau + b_tau n + b n_tau; b_tau from k_tau, Q_tau = int 2 k' k_tau', c_tau, b_tau."""
    els, sc = TB_.at(t); c, b, ct, bt = sc['c'], sc['b'], sc['ct'], sc['bt']
    xs = [float(v.mid()) for v in TB_.xk]
    out = []; Qt0 = arb(0)
    for fx in FX.data:
        a, h = fx['a'], fx['h']
        E = max(0, min(len(xs) - 1, int(__import__('numpy').searchsorted(xs, float(a.mid()) + 1e-12, side='right') - 1)))
        el = els[E]; S_ = arb_poly([2 * (a - el['x']) / el['h'] - 1, 2 * h / el['h']])
        L = {k: TM(el[k](S_)) for k in ('kap', 'kx', 'kxx', 'kt', 'ktx', 'ktxx', 'p', 'q', 't1', 'sn', 'r', 'Phi', 'rt', 'Q')}
        k0, k1, k2 = el['kap'](S_), el['kx'](S_), el['kxx'](S_)
        t0_, t1_, t2_ = el['kt'](S_), el['ktx'](S_), el['ktxx'](S_)
        p, q, n, tt1 = L['p'], L['q'], -1 * L['sn'], L['t1']
        Vs = -1 * L['kxx'] + c * p + b * n
        tht = L['Phi'] + L['rt']
        pt = (Vs + L['r']) - tht * q; nt = -1 * (tht * tt1)
        Vst = -1 * L['ktxx'] + ct * p + c * pt + bt * n + b * nt
        Qt = (2 * k1 * t1_).integral() * h + Qt0; x = arb_poly([a, h])
        bs = (-(k0 * k0), -2 * k0 * k1, -2 * k1 * k1 - 2 * k0 * k2, -6 * k0 * k1 + el['Q'](S_) - x * c,
              -5 * k1 * k1 - 6 * k0 * k2 - arb_poly([c]), arb_poly([c]) - 3 * k1 * k1 - 4 * k0 * k2)
        btau = (-2 * k0 * t0_, -2 * (k1 * t0_ + k0 * t1_), -2 * (k2 * t0_ + 2 * k1 * t1_ + k0 * t2_),
                -6 * (t0_ * k1 + k0 * t1_) + Qt - x * ct, -10 * k1 * t1_ - 6 * (t0_ * k2 + k0 * t2_) - arb_poly([ct]),
                arb_poly([ct]) - 6 * k1 * t1_ - 4 * (t0_ * k2 + k0 * t2_))
        phis = fx['phis']; chi0 = fx['chi'][0]
        ells = [CD.ell_strong(phis[j], bs, fx['d']) for j in range(3)]
        lts = [ell_b(phis[j], btau, fx['d']) for j in range(3)]
        out.append([chi0 * v for v in (p, n, Vs)] + [chi0 * v for v in (pt, nt, Vst)] + ells + lts)
        Qt0 = Qt(arb(1))
    return out


def gram15(FX, funcs):
    """15x15 Gram (int e^{-phi} f g) of [12 functions, phi_1..3]"""
    m = 15; G = [[arb(0)] * m for _ in range(m)]
    for fx, fs in zip(FX.data, funcs):
        allf = fs + [fx['phis'][j][0] for j in range(3)]
        w = [fx['e1'] * f for f in allf]
        for i in range(m):
            for j in range(i, m):
                v = (w[i] * allf[j]).integral01() * fx['h']; G[i][j] += v
                if j != i: G[j][i] += v
    return C0.amat(G)


def np_mid(A):
    import numpy as np
    return np.array([[float(A[i, j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])


def form_norm(G, K):
    """rigorous upper bound of sup |h^T K h| / |w|^2 over w (in the space where G is the Gram of the representers of h):
    float eigenbasis V of mid(G); G' = V^-1 G V^-T ~ diag(s); d_i = max(s_i, 0) + delta (delta = 100 x the largest radius /
    off-diagonal of G'); certify D^-1/2 G' D^-1/2 <= (1 + eps) I (Gershgorin); then |h^T K h| <= rho(D^1/2 V^T K V D^1/2)(1 + eps)|w|^2
    (y = V^-1 h... here h = V y, |y_i|^2 weighted by D^-1 is <= (1+eps)|w|^2)."""
    import numpy as np
    n = G.nrows(); Gm = np_mid(G); Gm = 0.5 * (Gm + Gm.T); s, V = np.linalg.eigh(Gm)
    Va = C0.np2arb(V); Vi = Va.inv(); Gp = Vi * G * Vi.transpose()
    radmax = max(float(Gp[i, j].rad()) + (abs(float(Gp[i, j].mid())) if i != j else 0.0) for i in range(n) for j in range(n))
    delta = max(100 * radmax, 1e-40)
    d = [max(float(Gp[i, i].mid()), 0.0) + delta for i in range(n)]
    Dih = arb_mat(n, n, [(1 / arb(d[i]).sqrt() if i == j else arb(0)) for i in range(n) for j in range(n)])
    Gn = Dih * Gp * Dih
    eps = max((sum((abs(Gn[i, j]) for j in range(n) if j != i), arb(0)) + Gn[i, i] - 1).upper() for i in range(n))
    eps = arb(max(float(eps), 0.0))
    Dh = arb_mat(n, n, [(arb(d[i]).sqrt() if i == j else arb(0)) for i in range(n) for j in range(n)])
    S = Dh * (Va.transpose() * K * Va) * Dh; S = (S + S.transpose()) * arb(0.5)
    rho = lmax_abs_up(S)
    return rho * (1 + eps), dict(eps=float(eps.upper()), delta=delta, smax=float(max(s)), smin=float(min(s)))


def projected(G15, idx):
    """Gram of the representers idx projected orthogonally to the phi representers (indices 12..14)"""
    n = len(idx)
    G11 = arb_mat(n, n, [G15[i, j] for i in idx for j in idx])
    G12 = arb_mat(n, 3, [G15[i, 12 + j] for i in idx for j in range(3)])
    G22 = arb_mat(3, 3, [G15[12 + i, 12 + j] for i in range(3) for j in range(3)])
    return G11 - G12 * G22.inv() * G12.transpose()


def lmax_abs_up(S):
    """rigorous upper bound of the spectral radius of a symmetric ball matrix: float eigenbasis congruence + Gershgorin"""
    import numpy as np
    n = S.nrows(); Sm = np_mid(S); Sm = 0.5 * (Sm + Sm.T); w, V = np.linalg.eigh(Sm); Va = C0.np2arb(V)
    VtV = Va.transpose() * Va
    e = max((sum((abs(VtV[i, j]) for j in range(n) if j != i), arb(0)) + abs(VtV[i, i] - 1)).upper() for i in range(n))
    e = arb(float(e))
    if not float(e.upper()) < 0.5: return None
    C = Va.transpose() * S * Va
    m = max((abs(C[i, i]) + sum((abs(C[i, j]) for j in range(n) if j != i), arb(0))).upper() for i in range(n))
    return arb(float(m)) / (1 - e)        # |y^T C y| <= m |y|^2, x = V y, |x|^2 >= (1 - e) |y|^2


def rate_bound(G15):
    """||d Pi / d tau||_Z, d Pi/d tau = sum_i delta_i m_i + sum_i ahat_i l_tau,i with delta = a_tau - M_tau M^-1 a (in Z),
    m = M^-T l, ahat = M^-1 a (functionals; M = <chi g, phi>_rho, M_tau = <chi g_tau, phi>_rho); each part by form_norm on the
    projected 6x6 Gram with K = [[0, I/2], [I/2, 0]]."""
    M = arb_mat(3, 3, [G15[i, 12 + j] for i in range(3) for j in range(3)])
    Mt = arb_mat(3, 3, [G15[3 + i, 12 + j] for i in range(3) for j in range(3)])
    Mi = M.inv(); Tm = Mt * Mi
    G12 = projected(G15, list(range(12)))
    Ca = arb_mat(6, 12); Cb = arb_mat(6, 12); MiT = Mi.transpose()
    for i in range(3):
        Ca[i, 3 + i] = arb(1); Cb[3 + i, 9 + i] = arb(1)
        for j in range(3):
            Ca[i, j] = -Tm[i, j]; Ca[3 + i, 6 + j] = MiT[i, j]; Cb[i, j] = Mi[i, j]
    K6 = arb_mat(6, 6)
    for i in range(3): K6[i, 3 + i] = K6[3 + i, i] = arb(1) / 2
    na, ia = form_norm(Ca * G12 * Ca.transpose(), K6)
    nb, ib = form_norm(Cb * G12 * Cb.transpose(), K6)
    return na + nb, dict(part_dir=float(na.upper()), part_ell=float(nb.upper()), eps=max(ia['eps'], ib['eps']), smin=min(ia['smin'], ib['smin']))


def c3_from_gram(G15):
    """(C3) at this base: mbar = ||M_n^-1|| (normalised g, phi), lbar = (sum ||ell_j||^2 / ||phi_j||^2)^1/2, C_x on Z, n_V = ||chi V_s||"""
    import numpy as np
    ng = [G15[i, i].sqrt() for i in range(3)]; nphi = [G15[12 + j, 12 + j].sqrt() for j in range(3)]
    Mn = arb_mat(3, 3, [G15[i, 12 + j] / (ng[i] * nphi[j]) for i in range(3) for j in range(3)]); Mni = Mn.inv()
    GP = arb_mat(3, 3, [G15[12 + i, 12 + j] for i in range(3) for j in range(3)]); GPi = GP.inv()
    PL = arb_mat(3, 3, [G15[12 + i, 6 + j] for i in range(3) for j in range(3)])
    GL = arb_mat(3, 3, [G15[6 + i, 6 + j] for i in range(3) for j in range(3)])
    GLZ = GL - PL.transpose() * GPi * PL
    Dinv = arb_mat(3, 3, [(1 / nphi[i] if i == j else arb(0)) for i in range(3) for j in range(3)])
    Cx = lmax_abs_up(Mni.transpose() * (Dinv * GLZ * Dinv) * Mni).sqrt()
    Cx_full = lmax_abs_up(Mni.transpose() * (Dinv * GL * Dinv) * Mni).sqrt()
    mbar = lmax_abs_up(Mni * Mni.transpose()).sqrt()
    lbar = sum((GL[j, j] / GP[j, j] for j in range(3)), arb(0)).sqrt()
    return dict(mbar=mbar, lbar=lbar, Cx=Cx, Cx_full=Cx_full, nV=ng[2])


def lmax_up(G):
    """upper bound of lambda_max of a symmetric ball matrix (Gershgorin), clipped at 0"""
    n = G.nrows()
    m = max((G[i, i] + sum((abs(G[i, j]) for j in range(n) if j != i), arb(0))).upper() for i in range(n))
    return arb(max(float(m), 0.0))


def nrm_up(A):
    """Frobenius upper bound"""
    return sum((A[i, j] ** 2 for i in range(A.nrows()) for j in range(A.ncols())), arb(0)).sqrt()


def rate_norms(C, R):
    """per unit tau, sup over the ball: ||A_tau||_Z, ||Lambda_tau||_Z, |M_tau| (entrywise), ||chi g_i,tau|| (unprojected)"""
    A = C0.amat; GPi = A(C['GP']).inv()
    Mt = A(R['M']); PLt = A(R['PL'])
    GgZ = A(R['Gg']) - Mt * GPi * Mt.transpose(); GlZ = A(R['Gl']) - PLt.transpose() * GPi * PLt
    return dict(dA=lmax_up(GgZ).sqrt(), dL=lmax_up(GlZ).sqrt(), dM=[[up(Mt[i, j]) for j in range(3)] for i in range(3)],
                dg=[up(A(R['Gg'])[i, i]).sqrt() for i in range(3)])


def r3_bounds(C, acc):
    """c_R and (C3) data from the tau_c Grams C and accumulated change bounds acc = dict(dA, dL, dM (3x3 entry bounds), dg)"""
    A = C0.amat
    GP = A(C['GP']); GPi = GP.inv(); Mc = A(C['M'])
    dM = arb_mat(3, 3, [arb(0, float(acc['dM'][i][j].upper())) for i in range(3) for j in range(3)])
    Mb = Mc + dM; Mci = Mc.inv(); Mbi = Mb.inv()
    GaaZ = A(C['GG']) - Mc * GPi * Mc.transpose()
    PL = A(C['PL']); GllZc = A(C['GL']) - PL.transpose() * GPi * PL
    nA = lmax_up(Mci * GaaZ * Mci.transpose()).sqrt(); nAb = lmax_up(Mbi * GaaZ * Mbi.transpose()).sqrt()
    nL = lmax_up(Mbi.transpose() * GllZc * Mbi).sqrt()
    ndA, ndL, ndM, nMbi = acc['dA'], acc['dL'], nrm_up(dM), nrm_up(Mbi)
    cR = nA * ndM * nL + nAb * ndL + ndA * nL + ndA * nMbi * ndL
    GG = A(C['GG']); GLc = A(C['GL'])
    ng = [GG[i, i].sqrt() + arb(0, float(acc['dg'][i].upper())) for i in range(3)]; nphi = [GP[i, i].sqrt() for i in range(3)]
    Mn = arb_mat(3, 3, [Mb[i, j] / (ng[i] * nphi[j]) for i in range(3) for j in range(3)]); Mni = Mn.inv()
    nphimin = min(nphi, key=lambda v: float(v.lower()))
    lbar = sum((GLc[j, j] / GP[j, j] for j in range(3)), arb(0)).sqrt() + ndL / nphimin
    Dinv = arb_mat(3, 3, [(1 / nphi[i] if i == j else arb(0)) for i in range(3) for j in range(3)])
    Cx = lmax_up(Mni.transpose() * (Dinv * GllZc * Dinv) * Mni).sqrt() + nrm_up(Mni) * ndL / nphimin
    nV = GG[2, 2].sqrt() - acc['dg'][2]
    return dict(cR=cR, nA=nA, nL=nL, ndA=ndA, ndL=ndL, ndM=ndM, mbar=nrm_up(Mni), lbar=lbar, Cx=Cx, nV=nV)


def coarse_nodes(Bc, fine, hc, xsw=3.0):
    """fine nodes for x < xsw (keeps the phi-weight series remainders small), the hc-mesh beyond"""
    co = CD.mesh(Bc, hc)
    return [v for v in fine if float(v.mid()) < xsw - 1e-9] + [v for v in co if float(v.mid()) >= xsw - 1e-9]


def c_phi(FX):
    """(C3) c_phi >= max_j ||rho^-1 (rho phi_j)''||_rho  (= f'' + 2 r1 f' + r2 f, r1 = -phi', r2 = phi'^2 - phi'')"""
    acc = [arb(0)] * 3
    for fx in FX.data:
        p1, p2 = fx['d'][0], fx['d'][1]
        for j in range(3):
            f = fx['phis'][j]; v = f[2] - 2 * (p1 * f[1]) + (p1 * p1 - p2) * f[0]
            acc[j] += (fx['e1'] * v * v).integral01() * fx['h']
    nphi = [sum(((fx['e1'] * fx['phis'][j][0] * fx['phis'][j][0]).integral01() * fx['h'] for fx in FX.data), arb(0)).sqrt() for j in range(3)]
    return [float((a.sqrt() / n).upper()) for a, n in zip(acc, nphi)], [float(a.sqrt().upper()) for a in acc]


# --------------------------------------------------------------------------------------------- assembly at tau_c
def assemble(Bc, FXf, th):
    t0 = time.time(); nodes = FXf.nodes; ne = len(nodes) - 1; ndof = 2 * (ne + 1)
    drop = {1, ndof - 2, ndof - 1}; free = [i for i in range(ndof) if i not in drop]; fmap = {g: k for k, g in enumerate(free)}; n = len(free)
    Band = C0.Band
    KQ, MW, K1, K2 = Band(n), Band(n), Band(n), Band(n)
    Gv = [[arb(0)] * 3 for _ in range(n)]; Pv = [[arb(0)] * 3 for _ in range(n)]; Lv = [[arb(0)] * 3 for _ in range(n)]
    z = lambda: [[arb(0)] * 3 for _ in range(3)]
    M3, GG, GP, GL_ = z(), z(), z(), z()
    geo = CD.Geo(Bc); elinfo = []
    unit = {k: TM(arb_poly([1 if k == (0, 0) else 0])) for k in PAIRS}
    unit1 = {k: TM(arb_poly([1 if k == (1, 1) else 0])) for k in PAIRS}
    unit2 = {k: TM(arb_poly([1 if k == (2, 2) else 0])) for k in PAIRS}
    for e in range(ne):
        fx = FXf.data[e]; a, h = fx['a'], fx['h']; E = CD.t1_index(Bc, a)
        G = geo.element(E, a, h); L = G['L']; b2, b2p, b2pp, b1, b1p, b0 = CD.coeff_polys(L, Bc, a, h)
        d, em, e1 = fx['d'], fx['em'], fx['e1']; p1, p2, p3, p4 = d; cphi = fx['cphi']
        e2 = TM(b2) - p1 * p1 + p2
        Cs = {(2, 2): TM(arb_poly([1 - th])), (2, 1): -th * p1, (2, 0): (1 - th) * cphi - 0.5 * (p1 * p1) - 0.5 * e2,
              (1, 1): -(1 + th) * (p1 * p1), (1, 0): -th * (p1 * cphi) - 0.5 * (p1 * p1 * p1) - 0.5 * (e2 * p1) - 0.5 * TM(b1),
              (0, 0): (1 - th) * (cphi * cphi) - (p1 * p1) * cphi - e2 * cphi - 0.5 * (TM(b1) * p1) - TM(b0)}
        Ke = C0.local_matrix(Cs, h); Me = C0.local_matrix(unit, h); K1e = C0.local_matrix(unit1, h); K2e = C0.local_matrix(unit2, h)
        gdof = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        for k in range(4):
            if gdof[k] not in fmap: continue
            for l in range(k, 4):
                if gdof[l] not in fmap: continue
                i, j = fmap[gdof[k]], fmap[gdof[l]]
                KQ.add(i, j, Ke[k][l]); MW.add(i, j, Me[k][l]); K1.add(i, j, K1e[k][l]); K2.add(i, j, K2e[k][l])
        chi = fx['chi']; phis = fx['phis']
        g = [G['p'], G['n'], TM(-1 * L['K'][2]) + Bc.c * G['p'] + Bc.b * G['n']]
        ells = [CD.ell_strong(f, (b2, b2p, b2pp, b1, b1p, b0), d) for f in phis]
        for ii in range(3):
            vg = C0.local_vec(em * chi[0] * g[ii], h); vp = C0.local_vec(em * phis[ii][0], h); vl = C0.local_vec(em * ells[ii], h)
            for k in range(4):
                if gdof[k] in fmap:
                    r = fmap[gdof[k]]; Gv[r][ii] += vg[k]; Pv[r][ii] += vp[k]; Lv[r][ii] += vl[k]
            for jj in range(3):
                M3[ii][jj] += (e1 * chi[0] * g[ii] * phis[jj][0]).integral01() * h
                GG[ii][jj] += (e1 * (chi[0] * g[ii]) * (chi[0] * g[jj])).integral01() * h
                GP[ii][jj] += (e1 * phis[ii][0] * phis[jj][0]).integral01() * h
                GL_[ii][jj] += (e1 * ells[ii] * ells[jj]).integral01() * h
        C21p = -th * p2
        cphip = 0.5 * p3 + 0.5 * (p1 * p2); e2p = TM(b2p) - 2 * (p1 * p2) + p3
        C20p = (1 - th) * cphip - (p1 * p2) - 0.5 * e2p
        S = {k: Cs[k].sup() for k in PAIRS}
        elinfo.append(dict(h=h, S=S, S21p=C21p.sup(), S20p=C20p.sup(), x=a, E=E))
    return dict(n=n, ne=ne, nodes=nodes, KQ=KQ, MW=MW, K1=K1, K2=K2, Gv=Gv, Pv=Pv, Lv=Lv, M3=M3, GG=GG, GP=GP, GL=GL_,
                elinfo=elinfo, k=Bc.k, t_assemble=time.time() - t0, free=free)


def c2_local(D, FXf, Bc, Bs, mu):
    """update per-element margins mu = [mu0_e, mu1_e] with the balls Bs"""
    for e, (fx, el) in enumerate(zip(FXf.data, D['elinfo'])):
        a, h = fx['a'], fx['h']; E = el['E']; p1, p2 = fx['d'][0], fx['d'][1]; cphi = fx['cphi']
        bc = CD.coeff_polys(Bc.local(E, a, h), Bc, a, h)
        for B in Bs:
            bb = CD.coeff_polys(B.local(E, a, h), B, a, h)
            db2, db2p, db2pp, db1, db1p, db0 = [x - y for x, y in zip(bb, bc)]
            m1 = sup01(db2)
            # -db2''/2 + db1'/2 - db0 = Delta(1.5 k'^2 + 2 k k'' - 1.5 c)  (b2'' = -2k'^2 - 2kk'', b1' = -5k'^2 - 6kk'' - c, b0 = c - 3k'^2 - 4kk'')
            Lb, Lc = B.local(E, a, h), Bc.local(E, a, h)
            core = lambda L_, B_: 1.5 * L_['K'][1] * L_['K'][1] + 2 * L_['K'][0] * L_['K'][2] - arb_poly([1.5 * B_.c])
            V0 = TM(core(Lb, B) - core(Lc, Bc)) + 0.5 * (TM(db2p) * p1 + TM(db2) * p2) - TM(db2) * cphi - 0.5 * (TM(db1) * p1)
            m0 = V0.sup()
            mu[0][e] = max(mu[0][e], up(m0), key=lambda v: float(v.upper()))
            mu[1][e] = max(mu[1][e], up(m1), key=lambda v: float(v.upper()))


def stepB(D, mu, cmax):
    """cert.stepB_static with the margins folded into the coefficient sups (S00 + cmax + mu0_e, S11 + mu1_e)"""
    for el, m0, m1 in zip(D['elinfo'], mu[0], mu[1]):
        S = dict(el['S']); S[(0, 0)] = S[(0, 0)] + cmax + m0; S[(1, 1)] = S[(1, 1)] + m1; el['S'] = S
    return C0.stepB_static(D)


def margin_band(D, mu):
    dl = [[mu[0][e], mu[1][e], arb(0)] for e in range(D['ne'])]
    return C0.margin_band(D, dl)


def check(D, SB, cs, Mb):
    """C0.check with KQ replaced by KQ - band(mu)"""
    KQ0 = D['KQ']
    if Mb is not None: D['KQ'] = KQ0.lin([(arb(-1), Mb)])
    try: r = C0.check(D, SB, cs)
    finally: D['KQ'] = KQ0
    return r


def bisect(D, SB, Mb, lo, hi, nbis, label):
    hist = []; best = None
    r = check(D, SB, arb(lo), Mb); hist.append((lo, r.get('ok'), r.get('why')))
    print(f'  [{label}] c = {lo:.5f}: ok={r["ok"]} {r.get("why", "")}', flush=True)
    if not r['ok']: return None, hist
    best = lo
    for _ in range(nbis):
        mid = 0.5 * (lo + hi); r = check(D, SB, arb(mid), Mb); hist.append((mid, r.get('ok'), r.get('why')))
        print(f'  [{label}] c = {mid:.5f}: ok={r["ok"]} {r.get("why", "")}', flush=True)
        if r['ok']: lo = best = mid
        else: hi = mid
    return best, hist


def run_cert(T, h=0.1, hc=0.5, nbis=7, nsub=2, outdir='.', FXc=None, FXf=None, BL=None):
    t0 = time.time(); lo, hi = T['slabs']; sc = T['slab_c']
    Sc = base7.Slab(sc); Bc = Sc.base(arb(-1), CD.XE + 1)
    print(f"cert {T['id']}: slabs {lo}-{hi}, tau_c = {Sc.tau0} (slab {sc}), k_c = {Bc.k:.2f}", flush=True)
    # (C2) local part: element margins from nsub tc-balls per slab of I (sharp base, BP-first in tc)
    ne = len(FXf.nodes) - 1; mu = [[arb(0)] * ne, [arb(0)] * ne]; per_slab = {}
    for s in range(lo, hi + 1):
        S = Sc if s == sc else base7.Slab(s); ts = time.time()
        balls = [arb((2 * j + 1) / nsub - 1, 1 / nsub) for j in range(nsub)]
        c2_local(T['_D'], FXf, Bc, [S.base(t, CD.XE + 1) for t in balls], mu)
        per_slab[s] = dict(slab=s, tau=[S.tau0, S.tau0 + S.dtf], k=S.k0, t=time.time() - ts,
                           mu0_cum=max(float(v.upper()) for v in mu[0]), mu1_cum=max(float(v.upper()) for v in mu[1]))
    # (C2) rank-3 part: exact rate ||dPi/dtau||_Z at tau_c (c2r3); reserve c_R = 1.6 n_side (F_c e^{0.3 H} L + 1e-3); c2r3.py certifies
    # sum_{slabs on each side} V_rank3 <= c_R (merge.py checks)
    import c2r3 as R3
    tr = time.time(); SBc = R3.SlabBP(Sc, CD.XE + 1); Gc = R3.gram(FXc, R3.funcs(SBc, arb(-1), FXc, order=1))
    import r3core as _RC
    Fc, finfo = _RC.rate_bound(Gc); c3c = c3_from_gram(Gc)
    H = max(sum(per_slab[s]['tau'][1] - per_slab[s]['tau'][0] for s in range(sc, hi + 1)),
            sum(per_slab[s]['tau'][1] - per_slab[s]['tau'][0] for s in range(lo, sc)) if sc > lo else 0.0)
    import math as _m
    L_ = Sc.dtf; nside = H / L_
    cR = arb(1.6 * nside * (float(Fc.upper()) * _m.exp(0.3 * H) * L_ + 1.0e-3))   # first order + c2r3's second-order term (<= 1.03e-3 per slab)
    print(f"   rank-3: ||dPi/dtau||_Z(tau_c) <= {float(Fc.upper()):.4e}, H = {H:.3f} -> reserved c_R = {float(cR.upper()):.3e}; "
          f"C3(tau_c): mbar {float(c3c['mbar'].upper()):.3f} lbar {float(c3c['lbar'].upper()):.3f} Cx {float(c3c['Cx'].upper()):.3f} "
          f"nV {float(c3c['nV'].lower()):.4e} ({time.time() - tr:.0f}s)", flush=True)
    per_slab = [per_slab[s] for s in sorted(per_slab)]
    print(f"   local margins: mu0 <= {max(float(v.upper()) for v in mu[0]):.3e} (x<10: {max(float(v.upper()) for v, el in zip(mu[0], T['_D']['elinfo']) if float(el['x'].mid()) < 10):.3e}), mu1 <= {max(float(v.upper()) for v in mu[1]):.3e}", flush=True)
    D = T['_D']; cmax = arb(0.09)
    tB = time.time(); SB = stepB(D, mu, cmax + up(cR)); tB = time.time() - tB
    Mb = margin_band(D, mu)
    cI, hI = bisect(D, SB, Mb, T.get('lo', 0.04), T.get('hi', 0.08), nbis, 'C1+C2')
    cIc = None if cI is None else cI - float(up(cR).upper())             # the band certificate is at shift cs = c_I + c_R
    cc, hc_ = bisect(D, SB, None, cI if cI is not None else 0.04, T.get('hi', 0.08) + 0.005, max(3, nbis - 3), 'C1 at tau_c')
    out = dict(id=T['id'], slabs=[lo, hi], slab_c=sc, tau_c=Sc.tau0, tau_I=[per_slab[0]['tau'][0], per_slab[-1]['tau'][1]],
               k_c=Bc.k, theta=float(THETA_I), c_I=cIc, c_shift_band=cI, c_R=float(up(cR).upper()), c_tau_c=cc,
               mu0_max=max(float(v.upper()) for v in mu[0]), mu1_max=max(float(v.upper()) for v in mu[1]),
               mu0_tip=max(float(v.upper()) for v, el in zip(mu[0], D['elinfo']) if float(el['x'].mid()) < 10),
               n=D['n'], ne=D['ne'], h=h, per_slab=per_slab, hist_I=hI, hist_c=hc_, t_stepB=tB, t_total=time.time() - t0,
               rate_Pi_tau_c=float(Fc.upper()), H=H, c_R_reserved=float(up(cR).upper()),
               C3_tau_c={k: (float(v.upper()) if k != 'nV' else float(v.lower())) for k, v in c3c.items()},
               mu_profile=[[float(el['x'].mid()), float(m0.upper()), float(m1.upper())] for el, m0, m1 in zip(D['elinfo'], mu[0], mu[1])])
    json.dump(out, open(os.path.join(outdir, 'cert7_%02d.json' % T['id']), 'w'), indent=1, default=str)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--times', default='0-0'); ap.add_argument('--out', default=os.path.join(HERE, 'out'))
    ap.add_argument('--h', type=float, default=0.1); ap.add_argument('--hc', type=float, default=0.5); ap.add_argument('--nbis', type=int, default=7)
    ap.add_argument('--nsub', type=int, default=2); ap.add_argument('--lo', type=float, default=0.04); ap.add_argument('--hi', type=float, default=0.08)
    ap.add_argument('--slabs', default=None, help='a-b: plan the certificates on slabs a..b only (default: whole store)')
    base7.base_arg(ap)
    A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); base7.set_base(A.base)
    ids = []
    for part in A.times.split(','):
        if '-' in part: a_, b_ = [int(v) for v in part.split('-')]; ids += list(range(a_, b_ + 1))
        elif part.strip(): ids.append(int(part))
    SCH = schedule(A.slabs); t0 = time.time()
    BL = base7.late_base(CD.XE + 1)
    FXf = None; FXc = None
    for T in [SCH[i] for i in ids]:
        if os.path.exists(os.path.join(A.out, 'cert7_%02d.json' % T['id'])): continue          # restartable
        Sc = base7.Slab(T['slab_c']); Bc = Sc.base(arb(-1), CD.XE + 1)
        if FXf is None:
            nodes = CD.mesh(Bc, A.h); FXf = Fixed(nodes, BL); FXc = Fixed(coarse_nodes(Bc, nodes, A.hc), BL)
            CPHI = c_phi(FXc)
            print(f'fixed functions: {len(nodes) - 1} fine / {len(FXc.nodes) - 1} coarse elements; c_phi (normalised, raw) = {CPHI} ({time.time() - t0:.0f}s)', flush=True)
        T['_D'] = assemble(Bc, FXf, THETA_I); T['lo'], T['hi'] = A.lo, A.hi
        print(f"assembled n = {T['_D']['n']} ({T['_D']['t_assemble']:.0f}s)", flush=True)
        o = run_cert(T, A.h, A.hc, A.nbis, A.nsub, A.out, FXc, FXf, BL)
        o['c_phi_normalised'], o['c_phi_raw'] = CPHI; o.pop('_D', None)
        json.dump(o, open(os.path.join(A.out, 'cert7_%02d.json' % T['id']), 'w'), indent=1, default=str)
        print(f"CERT {T['id']}: c_I = {o['c_I']}, c(tau_c) >= {o['c_tau_c']}, c_R = {o['c_R']:.2e}, mu0 <= {o['mu0_max']:.2e}, "
              f"mu1 <= {o['mu1_max']:.2e}  ({o['t_total']:.0f}s)", flush=True)
