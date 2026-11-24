"""T7 core functions without step-3 dependencies (portable): constraint-operator pieces, Gram/form-norm bounds, and a lossless
cache of the fixed functions (chi_32, phi_j jets, phi-weight jets, e^{-phi/2}, e^{-phi}) on the coarse mesh, so that the batch
scripts do not need step 3's c4fun (phi_3 = chi_32 psi3^) at run time.  The cache is produced once by cert7.Fixed (which uses
c4fun) and stored in fixed_coarse.json (every arb as exact (mantissa, exponent) pairs of mid and radius)."""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import numpy as np
from flint import arb, arb_poly, arb_mat, ctx
from tm import TM, from_series
import polstore
ctx.prec = 128
XE = arb(32) - arb('-2.249077651121112')            # end of supp chi_32 (x = sigma - STIP), as in certdata
CACHE = os.path.join(HERE, 'fixed_coarse.json')


# --------------------------------------------------------------------------------------- operator pieces
def ell_strong(f, b, d):
    """rho^-1 L^dag(rho f) from TMs f = (f..f''''), b = (b2, b2', b2'', b1, b1', b0) polys, d = (phi'..phi'''') (certdata.ell_strong)"""
    p1, p2, p3, p4 = d
    r1 = -1 * p1; r2 = p1 * p1 - p2; r3 = -1 * (p1 * p1 * p1) + 3 * (p1 * p2) - p3
    r4 = p1 * p1 * p1 * p1 - 6 * (p1 * p1 * p2) + 3 * (p2 * p2) + 4 * (p1 * p3) - p4
    f0, f1, f2, f3, f4 = f; b2, b2p, b2pp, b1, b1p, b0 = [TM(v) for v in b]
    t4 = f4 + 4 * (r1 * f3) + 6 * (r2 * f2) + 4 * (r3 * f1) + r4 * f0
    h0 = b2 * f0; h1 = b2p * f0 + b2 * f1; h2 = b2pp * f0 + 2 * (b2p * f1) + b2 * f2
    t2 = h2 + 2 * (r1 * h1) + r2 * h0
    t1 = (b1p * f0 + b1 * f1) + r1 * (b1 * f0)
    return -1 * t4 + t2 - t1 + b0 * f0


def ell_b(f, b, d):
    """b-linear part of rho^-1 L^dag(rho f) (cert7.ell_b)"""
    p1, p2 = d[0], d[1]
    r1 = -1 * p1; r2 = p1 * p1 - p2
    f0, f1, f2 = f[0], f[1], f[2]; b2, b2p, b2pp, b1, b1p, b0 = [TM(v) for v in b]
    h0 = b2 * f0; h1 = b2p * f0 + b2 * f1; h2 = b2pp * f0 + 2 * (b2p * f1) + b2 * f2
    t2 = h2 + 2 * (r1 * h1) + r2 * h0
    t1 = (b1p * f0 + b1 * f1) + r1 * (b1 * f0)
    return t2 - t1 + b0 * f0


def phi_tms(a, h, D=24):
    """phi' .. phi'''' (TMs), e^{-phi/2}, e^{-phi} (certdata.phi_tms)"""
    F = lambda X: (1 + X * X) ** (arb(2) / 3) / 4
    d = from_series(F, a, h, D, nder=4)[1:]
    em = from_series(lambda X: (-(1 + X * X) ** (arb(2) / 3) / 8).exp(), a, h, D)[0]
    e1 = from_series(lambda X: (-(1 + X * X) ** (arb(2) / 3) / 4).exp(), a, h, D)[0]
    return d, em, e1


def t1_index(ab, x):
    xf = float(x.mid()); xs = [float(v.mid()) for v in ab.xk]
    return max(0, min(int(np.searchsorted(xs, xf + 1e-12, side='right') - 1), len(ab.els) - 1))


# --------------------------------------------------------------------------------------- matrices
def amat(A): return arb_mat(len(A), len(A[0]), [A[i][j] for i in range(len(A)) for j in range(len(A[0]))])
def np2arb(A): return arb_mat(A.shape[0], A.shape[1], [arb(float(v)) for v in np.asarray(A).ravel()])
def np_mid(A): return np.array([[float(A[i, j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])


def lmax_abs_up(S):
    """rigorous upper bound of the spectral radius of a symmetric ball matrix: float eigenbasis congruence + Gershgorin"""
    n = S.nrows(); Sm = np_mid(S); Sm = 0.5 * (Sm + Sm.T); w, V = np.linalg.eigh(Sm); Va = np2arb(V)
    VtV = Va.transpose() * Va
    e = max((sum((abs(VtV[i, j]) for j in range(n) if j != i), arb(0)) + abs(VtV[i, i] - 1)).upper() for i in range(n))
    e = arb(float(e))
    if not float(e.upper()) < 0.5: return None
    C = Va.transpose() * S * Va
    m = max((abs(C[i, i]) + sum((abs(C[i, j]) for j in range(n) if j != i), arb(0))).upper() for i in range(n))
    return arb(float(m)) / (1 - e)


def form_norm(G, K):
    """rigorous upper bound of sup |h^T K h| / |w|^2, G = Gram of the representers of h (see cert7.form_norm)"""
    n = G.nrows(); Gm = np_mid(G); Gm = 0.5 * (Gm + Gm.T); s, V = np.linalg.eigh(Gm)
    Va = np2arb(V); Vi = Va.inv(); Gp = Vi * G * Vi.transpose()
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
    n = len(idx)
    G11 = arb_mat(n, n, [G15[i, j] for i in idx for j in idx])
    G12 = arb_mat(n, 3, [G15[i, 12 + j] for i in idx for j in range(3)])
    G22 = arb_mat(3, 3, [G15[12 + i, 12 + j] for i in range(3) for j in range(3)])
    return G11 - G12 * G22.inv() * G12.transpose()


def gscale(G, nblk):
    """rescale the direction functionals (blocks 0..nblk-1 of size 3: a, a_tau, ...) by 1/||chi g_i|| (float, tau-independent
    constants): Pi and all its tau-derivatives are invariant under g -> D g, and the scaling keeps the 3x3 systems well conditioned
    when ||chi V~|| is small (late times)"""
    n = [max(float(G[i, i].mid()), 1e-300) ** 0.5 for i in range(3)]
    m = G.nrows(); d = [arb(1 / n[i % 3]) if i < 3 * nblk else arb(1) for i in range(m)]
    D = arb_mat(m, m, [(d[i] if i == j else arb(0)) for i in range(m) for j in range(m)])
    return D * G * D


def rate_bound(G15, scale=True):
    """||d Pi/d tau||_Z (cert7.rate_bound): sum_i delta_i m_i + sum_i ahat_i l_tau,i, each part by form_norm on 6x6 projected Grams"""
    if scale: G15 = gscale(G15, 2)
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
    """(C3) at a time point: mbar = ||M_n^-1||, lbar, C_x on Z and on the full space, n_V = ||chi V_s||_rho (normalised g, phi)"""
    ng = [G15[i, i].sqrt() for i in range(3)]; nphi = [G15[12 + j, 12 + j].sqrt() for j in range(3)]
    Mn = arb_mat(3, 3, [G15[i, 12 + j] / (ng[i] * nphi[j]) for i in range(3) for j in range(3)]); Mni = Mn.inv()
    GP = arb_mat(3, 3, [G15[12 + i, 12 + j] for i in range(3) for j in range(3)]); GPi = GP.inv()
    PL = arb_mat(3, 3, [G15[12 + i, 6 + j] for i in range(3) for j in range(3)])
    GL = arb_mat(3, 3, [G15[6 + i, 6 + j] for i in range(3) for j in range(3)])
    GLZ = GL - PL.transpose() * GPi * PL
    Dinv = arb_mat(3, 3, [(1 / nphi[i] if i == j else arb(0)) for i in range(3) for j in range(3)])
    SZ = Mni.transpose() * (Dinv * GLZ * Dinv) * Mni
    Cx = lmax_abs_up(SZ).sqrt()
    Cxi = [arb(max(float(SZ[i, i].upper()), 0.0)).sqrt() for i in range(3)]          # |x_i| <= Cx_i ||u||
    Cx_full = lmax_abs_up(Mni.transpose() * (Dinv * GL * Dinv) * Mni).sqrt()
    mbar = lmax_abs_up(Mni * Mni.transpose()).sqrt()
    lbar = sum((GL[j, j] / GP[j, j] for j in range(3)), arb(0)).sqrt()
    return dict(mbar=mbar, lbar=lbar, Cx=Cx, Cx_full=Cx_full, nV=ng[2], Cx1=Cxi[0], Cx2=Cxi[1], Cx3=Cxi[2])


# --------------------------------------------------------------------------------------- fixed-function cache
def _enc(v): return polstore.enc(v)
def _dec(t): return polstore.dec(tuple(t))
def _etm(T_): return dict(p=[_enc(c) for c in T_.p.coeffs()], r=_enc(arb(T_.r)))
def _dtm(d): return TM(arb_poly([_dec(c) for c in d['p']]), _dec(d['r']))


def save_fixed(FX, path=CACHE, xtip=None):
    data = []
    for fx in FX.data:
        data.append(dict(a=_enc(fx['a']), h=_enc(fx['h']), d=[_etm(v) for v in fx['d']], em=_etm(fx['em']), e1=_etm(fx['e1']),
                         chi=[_etm(v) for v in fx['chi']], phis=[[_etm(v) for v in row] for row in fx['phis']], cphi=_etm(fx['cphi'])))
    json.dump(dict(nodes=[_enc(v) for v in FX.nodes], data=data, xtip=xtip, note='T7 fixed functions on the coarse mesh (cert7.Fixed)'), open(path, 'w'))


class FixedCache:
    def __init__(self, path=CACHE):
        J = json.load(open(path)); self.nodes = [_dec(t) for t in J['nodes']]; self.data = []; self.xtip = J.get('xtip')
        for d in J['data']:
            self.data.append(dict(a=_dec(d['a']), h=_dec(d['h']), d=[_dtm(v) for v in d['d']], em=_dtm(d['em']), e1=_dtm(d['e1']),
                                  chi=[_dtm(v) for v in d['chi']], phis=[[_dtm(v) for v in row] for row in d['phis']], cphi=_dtm(d['cphi'])))
