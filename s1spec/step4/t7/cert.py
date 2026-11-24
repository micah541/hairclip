"""T7 arb prototype of the (C1) certificate at one time (T5 2.4/1.2): rigorous lower bound for the projected form.
Claim certified: for all even u in H^4_rho supported in [0, XE] (XE = 32 - STIP, end of supp chi_32) with <u, phi_j>_rho = 0,
     B(u) <= -c_s ||u||_rho^2 - theta ||u''||_rho^2,     theta = 0.03,
for the largest c_s found by bisection.  With w = rho^(1/2) u this is Q_{c_s}(w) >= 0 (wfem.py), proved by
 (1) FEM block: Hermite cubic space V_h (all of it; banded), Q restricted to V_h ∩ null(Phi^T) minus explicit margins is positive:
     interval LDL^T of the banded part (inertia p) + Haynsworth on the 9x9 Schur complement of the rank-6 term and the 3 constraints;
 (2) tail: w = Pi_h w + w_r (Hermite interpolant; w_r = w_r' = 0 at nodes, int w_h'' w_r'' = 0 per element), explicit local
     constants (clamped beam 4.7300^-2 h^2, Wirtinger h/2pi), margins delta_{a,e} absorb all cross terms;
 (3) constraint leak: Pi_h w misses the constraints by eps_j = <w_r, e^{-phi/2} phi_j>, corrected inside V_h by eta (Phi^T eta = I).
usage: python3 cert.py [h] [slab] [theta_rel] -> cert_<slab>_<h>.json"""
import sys, time, json, numpy as np
from math import comb
from flint import arb, arb_poly, arb_mat, ctx
ctx.prec = 128
import certdata as CD
import tm
from tm import TM, sup01, up

THETA = arb(3) / 100
K0C = 1 / arb(4.7300040743) ** 2      # clamped-clamped: ||f|| <= (h/4.7300)^2 ||f''||  (4.73004074... = first root of cos z cosh z = 1)
K1C = 1 / (2 * arb.pi())                # f' in H^1_0 with zero mean: ||f'|| <= (h / 2 pi) ||f''||
CSMAX = arb(7) / 100                    # c_s range searched: [0, 0.07]
H = [arb_poly([1, 0, -3, 2]), arb_poly([0, 1, -2, 1]), arb_poly([0, 0, 3, -2]), arb_poly([0, 0, -1, 1])]
HD = [[q, q.derivative(), q.derivative().derivative()] for q in H]
PAIRS = [(2, 2), (2, 1), (2, 0), (1, 1), (1, 0), (0, 0)]
def _prepare():
    MOM, ABS = {}, {}
    for (a, b) in PAIRS:
        P = [[HD[k][a] * HD[l][b] + (HD[k][b] * HD[l][a] if a != b else arb_poly([0])) for l in range(4)] for k in range(4)]
        M = arb_mat(16, tm.DMAX + 1); Ab = [[None] * 4 for _ in range(4)]
        for k in range(4):
            for l in range(4):
                c = tm.coeffs(P[k][l])
                for j in range(tm.DMAX + 1):
                    M[4 * k + l, j] = sum((ci / (i + j + 1) for i, ci in enumerate(c)), arb(0))
                Ab[k][l] = tm.abs_int(P[k][l])
        MOM[(a, b)], ABS[(a, b)] = M, Ab
    return MOM, ABS
MOM, ABS = _prepare()

def local_matrix(Cs, h):
    """4x4 element matrix of sum_ab C_ab (sym) N^(a) N^(b) for TMs Cs[(a,b)], element length h"""
    sc = [arb(1), h, arb(1), h]; out = [[arb(0)] * 4 for _ in range(4)]
    for (a, b) in PAIRS:
        C = Cs[(a, b)]; c = tm.coeffs(C.p, tm.DMAX + 1)
        v = MOM[(a, b)] * arb_mat(tm.DMAX + 1, 1, c); f = h ** (1 - a - b); rr = C.r.upper()
        for k in range(4):
            for l in range(4):
                out[k][l] += (v[4 * k + l, 0] + arb(0, (rr * ABS[(a, b)][k][l]).upper())) * f * sc[k] * sc[l]
    return out

def local_vec(f, h):
    sc = [arb(1), h, arb(1), h]
    return [f.int_against(H[k]) * h * sc[k] for k in range(4)]

class Band:
    def __init__(self, n): self.n = n; self.b = [[arb(0)] * 4 for _ in range(n)]
    def add(self, i, j, v):
        if i > j: i, j = j, i
        if j - i <= 3: self.b[i][j - i] += v
        else: raise ValueError('band')
    def lin(self, others):
        """sum of coef * Band"""
        out = Band(self.n)
        for i in range(self.n):
            for d in range(4): out.b[i][d] = self.b[i][d] + sum((c * o.b[i][d] for c, o in others), arb(0))
        return out
    def matvec(self, x):
        y = [arb(0)] * self.n
        for i in range(self.n):
            for d in range(4):
                j = i + d
                if j >= self.n: break
                y[i] += self.b[i][d] * x[j]
                if d: y[j] += self.b[i][d] * x[i]
        return y
    def quad(self, x, z=None):
        z = x if z is None else z; y = self.matvec(z); return sum((xi * yi for xi, yi in zip(x, y)), arb(0))
    def mid(self):
        A = np.zeros((self.n, self.n))
        for i in range(self.n):
            for d in range(4):
                if i + d < self.n: A[i, i + d] = A[i + d, i] = float(self.b[i][d].mid())
        return A

def ldl(Bd):
    n = Bd.n; L = [[arb(0)] * 4 for _ in range(n)]; Dg = [arb(0)] * n
    for i in range(n):
        s = Bd.b[i][0]
        for k in range(max(0, i - 3), i): s -= L[k][i - k] * L[k][i - k] * Dg[k]
        if s.contains(0): return None, None, i
        Dg[i] = s
        for j in range(i + 1, min(n, i + 4)):
            t = Bd.b[i][j - i]
            for k in range(max(0, j - 3), i): t -= L[k][j - k] * L[k][i - k] * Dg[k]
            L[i][j - i] = t / s
    return L, Dg, None

def solve(L, Dg, b):
    n = len(Dg); y = list(b)
    for i in range(n):
        for k in range(max(0, i - 3), i): y[i] -= L[k][i - k] * y[k]
    z = [y[i] / Dg[i] for i in range(n)]; x = list(z)
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, min(n, i + 4)): x[i] -= L[i][j - i] * x[j]
    return x

def inertia_small(S):
    """rigorous inertia of a small symmetric ball matrix (list of lists) via float eigenvectors + Gershgorin; None if undecided"""
    m = len(S); Sm = np.array([[float(S[i][j].mid()) for j in range(m)] for i in range(m)]); Sm = 0.5 * (Sm + Sm.T)
    w, V = np.linalg.eigh(Sm); Va = arb_mat(m, m, [arb(float(V[i, j])) for i in range(m) for j in range(m)])
    Sa = arb_mat(m, m, [S[i][j] for i in range(m) for j in range(m)])
    C = Va.transpose() * Sa * Va
    if Va.det().contains(0): return None
    pos = neg = 0
    for i in range(m):
        rad = sum((abs(C[i, j]) for j in range(m) if j != i), arb(0)).upper() + C[i, i].rad()
        c = C[i, i].mid()
        if c > 0 and (c - rad) > 0: pos += 1
        elif c < 0 and (c + rad) < 0: neg += 1
        else: return None
    return pos, neg

def assemble(slab=0, theta=0.5, htarget=0.1, late=(187, 1.0), verbose=True):
    t0 = time.time()
    ab = CD.ArbBase(slab, theta, CD.XE + 1); abL = CD.ArbBase(*late, CD.XE + 1)
    nodes = CD.mesh(ab, htarget); ne = len(nodes) - 1; ndof = 2 * (ne + 1)
    drop = {1, ndof - 2, ndof - 1}; free = [i for i in range(ndof) if i not in drop]; fmap = {g: k for k, g in enumerate(free)}; n = len(free)
    KQ, MW, K1, K2 = Band(n), Band(n), Band(n), Band(n)
    Gv = [[arb(0)] * 3 for _ in range(n)]; Pv = [[arb(0)] * 3 for _ in range(n)]; Lv = [[arb(0)] * 3 for _ in range(n)]
    M3 = [[arb(0)] * 3 for _ in range(3)]; GG = [[arb(0)] * 3 for _ in range(3)]; GP = [[arb(0)] * 3 for _ in range(3)]; GL_ = [[arb(0)] * 3 for _ in range(3)]
    geo = CD.Geo(ab); geoL = CD.Geo(abL); elinfo = []
    one = TM(arb_poly([1])); th = THETA
    for e in range(ne):
        a = nodes[e]; h = nodes[e + 1] - a
        E = CD.t1_index(ab, a); EL = CD.t1_index(abL, a)
        G = geo.element(E, a, h); GL = geoL.element(EL, a, h, nder=4)
        L = G['L']; b2, b2p, b2pp, b1, b1p, b0 = CD.coeff_polys(L, ab, a, h)
        d, em, e1 = CD.phi_tms(a, h); p1, p2, p3, p4 = d
        cphi = 0.5 * p2 + 0.25 * (p1 * p1); e2 = TM(b2) - p1 * p1 + p2
        Cs = {(2, 2): TM(arb_poly([1 - th])), (2, 1): -th * p1, (2, 0): (1 - th) * cphi - 0.5 * (p1 * p1) - 0.5 * e2,
              (1, 1): -(1 + th) * (p1 * p1), (1, 0): -th * (p1 * cphi) - 0.5 * (p1 * p1 * p1) - 0.5 * (e2 * p1) - 0.5 * TM(b1),
              (0, 0): (1 - th) * (cphi * cphi) - (p1 * p1) * cphi - e2 * cphi - 0.5 * (TM(b1) * p1) - TM(b0)}
        Ke = local_matrix(Cs, h)
        unit = {(a_, b_): (TM(arb_poly([1])) if (a_, b_) == (0, 0) else TM(arb_poly([0]))) for (a_, b_) in PAIRS}
        Me = local_matrix(unit, h)
        unit1 = {k: TM(arb_poly([1 if k == (1, 1) else 0])) for k in PAIRS}; K1e = local_matrix(unit1, h)
        unit2 = {k: TM(arb_poly([1 if k == (2, 2) else 0])) for k in PAIRS}; K2e = local_matrix(unit2, h)
        gdof = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        for k in range(4):
            if gdof[k] not in fmap: continue
            for l in range(k, 4):
                if gdof[l] not in fmap: continue
                i, j = fmap[gdof[k]], fmap[gdof[l]]
                KQ.add(i, j, Ke[k][l]); MW.add(i, j, Me[k][l]); K1.add(i, j, K1e[k][l]); K2.add(i, j, K2e[k][l])
        # directions (base at tau_c) and constraints (fixed)
        chi = [TM(c) for c in CD.chi_jets(a, h)]
        g = [G['p'], G['n'], TM(-1 * L['K'][2]) + ab.c * G['p'] + ab.b * G['n']]
        pe = CD.psi_index(a)
        ps = CD.psi3_piece(pe, a, h) if pe is not None else CD.psi3_ext(a, h)
        bases = [GL['P'], GL['N'], ps]
        phis = [[sum((comb(m, kk) * chi[kk] * bs[m - kk] for kk in range(m + 1)), TM(arb_poly([0]))) for m in range(5)] for bs in bases]
        ells = [CD.ell_strong(f, (b2, b2p, b2pp, b1, b1p, b0), d) for f in phis]
        for ii in range(3):
            fg = em * chi[0] * g[ii]; vg = local_vec(fg, h)
            vp = local_vec(em * phis[ii][0], h); vl = local_vec(em * ells[ii], h)
            for k in range(4):
                if gdof[k] in fmap:
                    r = fmap[gdof[k]]; Gv[r][ii] += vg[k]; Pv[r][ii] += vp[k]; Lv[r][ii] += vl[k]
            for jj in range(3):
                M3[ii][jj] += (e1 * chi[0] * g[ii] * phis[jj][0]).integral01() * h
                GG[ii][jj] += (e1 * (chi[0] * g[ii]) * (chi[0] * g[jj])).integral01() * h
                GP[ii][jj] += (e1 * phis[ii][0] * phis[jj][0]).integral01() * h
                GL_[ii][jj] += (e1 * ells[ii] * ells[jj]).integral01() * h
        # Step-B local sups (C' only needed for C21, C20)
        C21p = -th * p2
        cphip = 0.5 * p3 + 0.5 * (p1 * p2); e2p = TM(b2p) - 2 * (p1 * p2) + p3
        C20p = (1 - th) * cphip - (p1 * p2) - 0.5 * e2p
        S = {k: Cs[k].sup() for k in PAIRS}; S[(0, 0)] = S[(0, 0)] + CSMAX
        elinfo.append(dict(h=h, S=S, S21p=C21p.sup(), S20p=C20p.sup(), x=a))
        if verbose and e % 50 == 0: print(f'  element {e}/{ne} x={float(a.mid()):.2f} ({time.time()-t0:.0f}s)', flush=True)
    return dict(n=n, ne=ne, nodes=nodes, KQ=KQ, MW=MW, K1=K1, K2=K2, Gv=Gv, Pv=Pv, Lv=Lv, M3=M3, GG=GG, GP=GP, GL=GL_,
                elinfo=elinfo, k=ab.k, kL=abL.k, t_assemble=time.time() - t0, free=free)

def amat(A): return arb_mat(len(A), len(A[0]), [A[i][j] for i in range(len(A)) for j in range(len(A[0]))])
def tr(A): return sum((A[i, i] for i in range(A.nrows())), arb(0))

def stepB_static(D):
    """c_s-independent tail constants: local beta_{a,e}, gamma_e; ||R||_op; Phi_kappa"""
    th = THETA; hmax = max(float(el['h'].upper()) for el in D['elinfo']); hm = arb(hmax)
    Mi = amat(D['M3']).inv(); GG, GL, GP = amat(D['GG']), amat(D['GL']), amat(D['GP'])
    Rop = tr(Mi * GG * Mi.transpose() * GL).sqrt()          # >= ||Gamma_G^1/2 M^-T Gamma_L^1/2||_2 = sup |g(a)^T M^-T l(b)| / (|a||b|)
    Phik = tr(GP).sqrt() * K0C * hm * hm                     # |eps| <= Phik ||w_r''||
    betas, gmax = [], arb(0); bmax = [arb(0)] * 3
    for el in D['elinfo']:
        hh = el['h']; S = el['S']; k1h = K1C * hh; k0h = K0C * hh * hh
        b2_ = 4 * S[(2, 1)] * k1h + 2 * S[(2, 0)] * k0h
        b1_ = 2 * el['S21p'] * k1h + 2 * S[(2, 0)] * k1h + 2 * S[(1, 1)] * k1h + 2 * S[(1, 0)] * k0h
        b0_ = 2 * el['S20p'] * k1h + 2 * S[(1, 0)] * k1h + 2 * S[(0, 0)] * k0h
        g_ = el['S21p'] * k1h * k1h + 2 * el['S20p'] * k0h * k1h + 2 * S[(2, 0)] * k1h * k1h + S[(1, 1)] * k1h * k1h + 2 * S[(1, 0)] * k0h * k1h + S[(0, 0)] * k0h * k0h
        betas.append([b0_, b1_, b2_]); gmax = max(gmax, g_, key=lambda v: float(v.upper()))
        bmax = [max(bmax[i], [b0_, b1_, b2_][i], key=lambda v: float(v.upper())) for i in range(3)]
    b0G = 2 * Rop * K0C * hm * hm                            # rank-3 term on w_r:  |R(w)-R(w_h)| <= Rop (2|w_h||w_r| + |w_r|^2)
    gR = Rop * (K0C * hm * hm) ** 2
    cg = congr_inertia(to_dense(D['MW']), D['MW'].mid())
    lmin = min((arb(c) - r) for c, r in zip(cg['cen'], cg['rad']))
    if not float(lmin.lower()) > 0: raise RuntimeError('lambda_min(M_w) not certified')
    return dict(Mi=Mi, Rop=Rop, Phik=Phik, betas=betas, gmax=gmax, bmax=bmax, b0G=b0G, gR=gR, hm=hm, lmin_M=arb(lmin.lower()))

def to_dense(B):
    n = B.n; M = arb_mat(n, n)
    for i in range(n):
        for d in range(4):
            j = i + d
            if j < n: M[i, j] = B.b[i][d]; M[j, i] = B.b[i][d]
    return M

def np2arb(A):
    return arb_mat(A.shape[0], A.shape[1], [arb(float(v)) for v in np.asarray(A).ravel()])

def gersh(C):
    """per-row Gershgorin data of a square arb_mat: centers (float mid), radii (exact upper bounds as arb)"""
    n = C.nrows(); cen, rad = [], []
    for i in range(n):
        r = arb(C[i, i].rad())
        for j in range(n):
            if j != i: r += abs(C[i, j])
        cen.append(float(C[i, i].mid())); rad.append(arb(r.upper()))
    return cen, rad

def congr_inertia(Kd, Kmid):
    """rigorous inertia of the symmetric ball matrix Kd: V^T Kd V (V float eigenvectors of Kmid) + Gershgorin; V invertible checked"""
    w, V = np.linalg.eigh(0.5 * (Kmid + Kmid.T)); Va = np2arb(V)
    VtV = Va.transpose() * Va; c2, r2 = gersh(VtV)
    if not all(float((arb(c) - r).lower()) > 0 for c, r in zip(c2, r2)): return None
    C = Va.transpose() * Kd * Va; cen, rad = gersh(C)
    pos = neg = 0
    for c, r in zip(cen, rad):
        if float((arb(c) - r).lower()) > 0: pos += 1
        elif float((arb(c) + r).upper()) < 0: neg += 1
        else: return None
    return dict(Va=Va, C=C, cen=cen, rad=rad, pos=pos, neg=neg)

def schur_quad(cg, E):
    """enclosure (arb_mat m x m) of E^T K^-1 E from the congruence data cg of K (K = V^-T C V^-1, C = Lam + Delta)"""
    Va, C, cen = cg['Va'], cg['C'], cg['cen']; n = C.nrows(); m = E.ncols()
    lam = [arb(c) for c in cen]; Y = Va.transpose() * E
    LY = arb_mat(n, m, [Y[i, j] / lam[i] for i in range(n) for j in range(m)])
    Dl = C - arb_mat(n, n, [lam[i] if i == j else arb(0) for i in range(n) for j in range(n)])
    q = max((cg['rad'][i] / abs(lam[i])) for i in range(n))
    q = arb(max(float(v.upper()) for v in [q]))
    if not float(q.upper()) < 0.5: return None
    T = q * q / (1 - q)
    S1 = Y.transpose() * LY; S2 = LY.transpose() * (Dl * LY)
    out = arb_mat(m, m)
    for i in range(m):
        n1 = sum((abs(Y[r, i]) for r in range(n)), arb(0))
        for j in range(m):
            ninf = max((abs(LY[r, j]) for r in range(n)), key=lambda v: float(v.upper()))
            out[i, j] = S1[i, j] - S2[i, j] + arb(0, (n1 * T * ninf).upper())
    return out

def full_mv(D, KB, x, Mi):
    """K_f x = (KQ - cs MW) x + U S U^T x  (S = 1/2 [[0, M^-T],[M^-1, 0]])"""
    n = D['n']; y = KB.matvec(x)
    gx = [sum((D['Gv'][r][i] * x[r] for r in range(n)), arb(0)) for i in range(3)]
    lx = [sum((D['Lv'][r][i] * x[r] for r in range(n)), arb(0)) for i in range(3)]
    a = [sum((Mi[j, i] * lx[j] for j in range(3)), arb(0)) / 2 for i in range(3)]      # (M^-T l)_i / 2  -> coefficient of G_i
    bb = [sum((Mi[i, j] * gx[j] for j in range(3)), arb(0)) / 2 for i in range(3)]     # (M^-1 g)_i / 2  -> coefficient of L_i
    return [y[r] + sum((D['Gv'][r][i] * a[i] + D['Lv'][r][i] * bb[i] for i in range(3)), arb(0)) for r in range(n)]

def eta_data(D, KB, cs, SB):
    """leak corrector: eta~ = float K_f-orthogonal corrector, eta = eta~ (Phi^T eta~)^-1 exactly (Phi^T eta = I);
    |Q_h(w_0, eta_j)| <= kappa_rho ||w_0|| for w_0 in null(Phi^T) (residual of K_f eta~ against range(Phi));
    returns etaK >= ||eta^T K_f eta||_F, N^eta, kappa_rho"""
    n = D['n']; Mi = SB['Mi']
    Kf = KB.mid(); Gm = np.array([[float(D['Gv'][r][i].mid()) for i in range(3)] for r in range(n)])
    Lm = np.array([[float(D['Lv'][r][i].mid()) for i in range(3)] for r in range(n)]); Pm = np.array([[float(D['Pv'][r][i].mid()) for i in range(3)] for r in range(n)])
    M3 = np.array([[float(D['M3'][i][j].mid()) for j in range(3)] for i in range(3)]); Mf = np.linalg.inv(M3)
    Kf = Kf + 0.5 * (Gm @ Mf.T @ Lm.T + Lm @ Mf @ Gm.T)
    X = np.linalg.solve(Kf, Pm); et = X @ np.linalg.inv(Pm.T @ X)
    cols = [[arb(float(et[r, j])) for r in range(n)] for j in range(3)]
    PT = arb_mat(3, 3, [sum((D['Pv'][r][i] * cols[j][r] for r in range(n)), arb(0)) for i in range(3) for j in range(3)]).inv()
    eta = [[sum((cols[k][r] * PT[k, j] for k in range(3)), arb(0)) for r in range(n)] for j in range(3)]
    Kfe_f = [full_mv(D, KB, cols[j], Mi) for j in range(3)]            # K_f eta~_j (arb)
    Xl = np.linalg.lstsq(Pm, np.array([[float(v.mid()) for v in Kfe_f[j]] for j in range(3)]).T, rcond=None)[0]
    res = [[Kfe_f[j][r] - sum((D['Pv'][r][i] * arb(float(Xl[i, j])) for i in range(3)), arb(0)) for r in range(n)] for j in range(3)]
    r2 = sum((sum((v * v for v in res[j]), arb(0)) for j in range(3)), arb(0))
    # |w0^T K_f eta_j| = |w0^T K_f eta~ PT| ; |c^T rho| <= ||w|| |rho| / sqrt(lam_min(M_w))
    PTn = sum((PT[i, j] ** 2 for i in range(3) for j in range(3)), arb(0)).sqrt()
    kap = (r2 / SB['lmin_M']).sqrt() * PTn
    Kfe = [full_mv(D, KB, eta[j], Mi) for j in range(3)]
    etaK = sum((sum((eta[i][r] * Kfe[j][r] for r in range(n)), arb(0)) ** 2 for i in range(3) for j in range(3)), arb(0)).sqrt()
    Neta = [sum((B.quad(eta[j]) for j in range(3)), arb(0)).sqrt() for B in (D['MW'], D['K1'], D['K2'])]
    return etaK, Neta, kap

def deltas(D, SB, etaK, Neta, kap):
    th = THETA; Phik = SB['Phik']; K0h2 = K0C * SB['hm'] * SB['hm']
    gtot = SB['gmax'] + Phik * Phik * etaK + sum((SB['bmax'][i] * Phik * Neta[i] for i in range(3)), arb(0)) + SB['gR'] \
           + SB['Rop'] * 2 * K0h2 * Phik * Neta[0]
    epsL = (1 - th - gtot) / 2
    if not float(epsL.lower()) > 0: return None, None
    b = [SB['b0G'] + 2 * Phik * kap, arb(0), arb(0)]
    epsL = arb(epsL.lower())                                       # exact lower bound -> margins are exact upper bounds (conservative)
    delta = [[up(3 * bb[i] ** 2 / (4 * epsL) + 3 * b[i] ** 2 / (4 * epsL)) for i in range(3)] for bb in SB['betas']]
    info = dict(Rop=float(SB['Rop'].upper()), Phik=float(Phik.upper()), etaK=float(etaK.upper()), Neta=[float(v.upper()) for v in Neta], kappa_rho=float(kap.upper()),
                b0G=float(SB['b0G'].upper()), gtot=float(gtot.upper()),
                delta_max=[max(float(dd[i].upper()) for dd in delta) for i in range(3)],
                delta_tip=[max(float(dd[i].upper()) for dd, el in zip(delta, D['elinfo']) if float(el['x'].mid()) < 10) for i in range(3)])
    return delta, info

def margin_band(D, delta):
    """sum_e (delta0 M_e + delta1 K1_e + delta2 K2_e), assembled"""
    n = D['n']; Bm = Band(n); fmap = {g: k for k, g in enumerate(D['free'])}
    unit = {i: {k: TM(arb_poly([1 if k == kk else 0])) for k in PAIRS} for i, kk in enumerate(((0, 0), (1, 1), (2, 2)))}
    for e, el in enumerate(D['elinfo']):
        h = el['h']; loc = [[arb(0)] * 4 for _ in range(4)]
        for i in range(3):
            m = local_matrix(unit[i], h)
            for k in range(4):
                for l in range(4): loc[k][l] += delta[e][i] * m[k][l]
        gdof = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        for k in range(4):
            if gdof[k] not in fmap: continue
            for l in range(k, 4):
                if gdof[l] in fmap: Bm.add(fmap[gdof[k]], fmap[gdof[l]], loc[k][l])
    return Bm

def check(D, SB, cs):
    """certificate test at shift cs: margins from step B, In(K') = (n-p, p, 0) and In(Schur) = (3+p, 6-p, 0), all rigorous"""
    n = D['n']; KB = D['KQ'].lin([(-cs, D['MW'])])
    etaK, Neta, kap = eta_data(D, KB, cs, SB)
    delta, info = deltas(D, SB, etaK, Neta, kap)
    if delta is None: return dict(ok=False, why='step B gamma too large')
    Kp = KB.lin([(arb(-1), margin_band(D, delta))])
    cg = congr_inertia(to_dense(Kp), Kp.mid())
    if cg is None: return dict(ok=False, why='K\' inertia undecided')
    p = cg['neg']
    E = arb_mat(n, 9, [ (D['Gv'][r][j] if j < 3 else (D['Lv'][r][j - 3] if j < 6 else D['Pv'][r][j - 6])) for r in range(n) for j in range(9)])
    Q9 = schur_quad(cg, E)
    if Q9 is None: return dict(ok=False, why='Neumann q too large')
    S = [[arb(0)] * 9 for _ in range(9)]
    for i in range(3):
        for j in range(3): S[i][3 + j] = -2 * D['M3'][i][j]; S[3 + j][i] = -2 * D['M3'][i][j]
    for i in range(9):
        for j in range(9): S[i][j] = S[i][j] - Q9[i, j]
    inn = inertia_small(S)
    S_rad = max(float(S[i][j].rad()) for i in range(9) for j in range(9))
    ok = inn is not None and inn == (3 + p, 6 - p)
    gmin = min(abs(c) - float(r.upper()) for c, r in zip(cg['cen'], cg['rad']))
    return dict(ok=ok, p=p, inertia=inn, S_rad=S_rad, gersh_gap_min=gmin, stepB=info)

if __name__ == '__main__':
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1; slab = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    br = [float(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else [0.0, 0.07]; nbis = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    T0 = time.time(); D = assemble(slab, thr, h)
    print(f"assembled: n = {D['n']} ({D['ne']} elements), k = {D['k']:.2f}, late k = {D['kL']:.0f}, {D['t_assemble']:.0f}s", flush=True)
    t1 = time.time(); SB = stepB_static(D); tB = time.time() - t1
    print(f"step B static: ||R||_op <= {float(SB['Rop'].upper()):.3e}, Phi_kappa = {float(SB['Phik'].upper()):.3e}, gamma_loc = {float(SB['gmax'].upper()):.2e}, beta_max = {[f'{float(v.upper()):.2e}' for v in SB['bmax']]} ({tB:.0f}s)", flush=True)
    t2 = time.time(); lo, hi = br; hist = []; info = None
    tc = time.time(); r0 = check(D, SB, arb(lo)); t_check = time.time() - tc; hist.append((lo, r0)); print(f'  c_s = {lo}:', r0, f'({t_check:.0f}s)', flush=True)
    if not r0['ok']: lo = float('nan')
    if r0['ok']:
        info = r0['stepB']
        for it in range(nbis):
            mid = 0.5 * (lo + hi); r = check(D, SB, arb(mid)); hist.append((mid, r))
            if r['ok']: info = r['stepB']
            if r['ok']: lo = mid
            else: hi = mid
            print(f'  c_s = {mid:.6f}: {r}', flush=True)
    tC = time.time() - t2
    out = dict(h=h, n=D['n'], ne=D['ne'], k=D['k'], slab=slab, theta=thr, c_cert=lo, stepB=info, t_assemble=D['t_assemble'], t_stepB=tB, t_bisect=tC, t_check=t_check,
               t_total=time.time() - T0, hist=[(c, {k: v for k, v in r.items()}) for c, r in hist])
    json.dump(out, open(f'cert_{slab}_{h}.json', 'w'), indent=1, default=str)
    print(f"CERTIFIED c_s = {lo:.6f}  (theta = 0.03; float truth 0.06887 at k=20)  total {time.time()-T0:.0f}s", flush=True)
