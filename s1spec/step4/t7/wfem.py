"""T7 float check of the w-formulation used by the arb certificate (cert.py).
u = e^{phi/2} w, phi = (1+x^2)^(2/3)/4, so ||u||_rho = ||w||, and with U2 = w'' + phi' w' + c w, U1 = w' + phi' w/2, c = phi''/2 + phi'^2/4:
  Q(w) := -B(u) - c_s ||u||_rho^2 - theta ||u''||_rho^2 = int [C22 w''^2 + 2 C21 w''w' + 2 C20 w''w + C11 w'^2 + 2 C10 w'w + C00 w^2] + R(w)
  C22 = 1-theta, C21 = -theta phi', C20 = (1-theta) c - phi'^2/2 - e2/2, C11 = -(1+theta) phi'^2,
  C10 = -theta phi' c - phi'^3/2 - e2 phi'/2 - b1/2,  C00 = (1-theta) c^2 - phi'^2 c - e2 c - b1 phi'/2 - (b0 + c_s),   e2 = b2 - phi'^2 + phi''
  R(w) = g(w)^T M^-T l(w),  g_i(w) = int e^{-phi/2} chi g_i w,  l_j(w) = int e^{-phi/2} ell_j w,  ell_j = rho^-1 L^dag(rho phi_j) (strong form).
Hermite cubic FEM in w on [0, X], X = 32 - STIP (end of supp chi_32), w'(0) = 0, w(X) = w'(X) = 0."""
import numpy as np, scipy.linalg as sla
from numpy.polynomial.legendre import leggauss
import sys; sys.path.insert(0, "/home/claude/s1spec/step4/t1")
import bgeo, fem
from rep import STIP

def phi_jets(x):
    """phi', phi'', phi''', phi'''' of phi = (1+x^2)^(2/3)/4"""
    g = 1 + x**2
    p1 = x * g ** (-1 / 3) / 3
    p2 = g ** (-1 / 3) / 3 - 2 * x**2 * g ** (-4 / 3) / 9
    p3 = -2 * x * g ** (-4 / 3) / 9 - 4 * x * g ** (-4 / 3) / 9 + 16 * x**3 * g ** (-7 / 3) / 27
    p4 = (-6 * g ** (-4 / 3) / 9 + (6 * 8 / 27) * x**2 * g ** (-7 / 3) + 48 * x**2 * g ** (-7 / 3) / 27 - (16 * 14 / 81) * x**4 * g ** (-10 / 3))
    return p1, p2, p3, p4

def rho_ratios(p1, p2, p3, p4):
    """rho^(k)/rho, rho = e^{-phi}"""
    r1 = -p1; r2 = p1**2 - p2; r3 = -p1**3 + 3 * p1 * p2 - p3; r4 = p1**4 - 6 * p1**2 * p2 + 3 * p2**2 + 4 * p1 * p3 - p4
    return r1, r2, r3, r4

def ell_strong(f, b, r):
    """rho^-1 L^dag(rho f), L^dag psi = -psi'''' + (b2 psi)'' - (b1 psi)' + b0 psi.  f = (f..f''''), b = (b2,b2',b2'', b1,b1', b0), r = rho ratios"""
    f0, f1, f2, f3, f4 = f; b2, b2p, b2pp, b1, b1p, b0 = b; r1, r2, r3, r4 = r
    t4 = f4 + 4 * r1 * f3 + 6 * r2 * f2 + 4 * r3 * f1 + r4 * f0
    h0 = b2 * f0; h1 = b2p * f0 + b2 * f1; h2 = b2pp * f0 + 2 * b2p * f1 + b2 * f2
    t2 = h2 + 2 * r1 * h1 + r2 * h0
    k0 = b1 * f0; k1 = b1p * f0 + b1 * f1
    t1 = k1 + r1 * k0
    return -t4 + t2 - t1 + b0 * f0

class WFEM:
    def __init__(self, h=0.1, nq=8):
        X = 32 - STIP; self.X = X; x1 = 29 - STIP
        n1 = int(np.ceil(x1 / h)); n2 = int(np.ceil((X - x1) / h))
        self.xn = np.r_[np.linspace(0, x1, n1 + 1), np.linspace(x1, X, n2 + 1)[1:]]
        self.ne = len(self.xn) - 1; he = np.diff(self.xn); self.he = he
        gx, gw = leggauss(nq); xi = (gx + 1) / 2
        self.xq = self.xn[:-1, None] + he[:, None] * xi[None, :]; self.wq = he[:, None] * (gw / 2)[None, :]
        self.Bl = [fem.hermite(xi, hh) for hh in he]          # per element (B0, B1, B2) arrays (4, nq)
        self.ndof = 2 * (self.ne + 1)
        self.free = np.array([i for i in range(self.ndof) if i not in (1, self.ndof - 2, self.ndof - 1)])
        self.idx = 2 * np.arange(self.ne)[:, None] + np.arange(4)[None, :]
    def scatter(self, loc):
        A = np.zeros((self.ndof, self.ndof))
        for e in range(self.ne): A[np.ix_(self.idx[e], self.idx[e])] += loc[e]
        return A[np.ix_(self.free, self.free)]
    def vec(self, loc):
        v = np.zeros(self.ndof)
        for e in range(self.ne): v[self.idx[e]] += loc[e]
        return v[self.free]
    def bilin(self, C):
        """C = dict (a, b) -> coefficient array at quadrature points, symmetric quadratic form sum_ab C_ab N^(a) N^(b)"""
        loc = np.zeros((self.ne, 4, 4))
        for e in range(self.ne):
            Bs = self.Bl[e]; w = self.wq[e]
            for (a, b), Cab in C.items():
                M = np.einsum('q,iq,jq->ij', w * Cab[e], Bs[a], Bs[b])
                loc[e] += M if a == b else (M + M.T)
        return self.scatter(loc)
    def load(self, f):
        return self.vec(np.array([np.einsum('q,iq->i', self.wq[e] * f[e], self.Bl[e][0]) for e in range(self.ne)]))

def build(geom, th, cs, h=0.1, phis=None):
    """float matrices of Q (theta = th, shift cs) for BaseGeom geom; returns dict"""
    F = WFEM(h); x = F.xq; s = geom.at(x.ravel(), mmax=3); sh = lambda a: a.reshape(x.shape)
    K = [sh(s['K'][m]) for m in range(4)]; Q = sh(s['Q'])
    b2 = -K[0]**2; b2p = -2 * K[0] * K[1]; b2pp = -2 * K[1]**2 - 2 * K[0] * K[2]
    b1 = -6 * K[0] * K[1] + Q - geom.c * x; b1p = -6 * K[1]**2 - 6 * K[0] * K[2] + K[1]**2 - geom.c
    b0 = 0.25 - geom.a - 3 * K[1]**2 - 4 * K[0] * K[2]
    p1, p2, p3, p4 = phi_jets(x); c = p2 / 2 + p1**2 / 4; e2 = b2 - p1**2 + p2
    C = {(2, 2): (1 - th) * np.ones_like(x), (2, 1): -th * p1, (2, 0): (1 - th) * c - p1**2 / 2 - e2 / 2, (1, 1): -(1 + th) * p1**2,
         (1, 0): -th * p1 * c - p1**3 / 2 - e2 * p1 / 2 - b1 / 2, (0, 0): (1 - th) * c**2 - p1**2 * c - e2 * c - b1 * p1 / 2 - (b0 + cs)}
    KQ = F.bilin(C)
    em = np.exp(-(1 + x**2) ** (2 / 3) / 8)
    chi = [bgeo.chi32(x, n) for n in range(5)]
    g = [sh(s['p']), sh(s['n']), sh(s['Vs'])]
    Gv = np.array([F.load(em * chi[0] * gi) for gi in g]).T
    if phis is None: phis = phi_funcs(x, chi)
    r = rho_ratios(p1, p2, p3, p4)
    ells = [ell_strong(f, (b2, b2p, b2pp, b1, b1p, b0), r) for f in phis]
    Lv = np.array([F.load(em * l) for l in ells]).T
    # weak form a(psi_k, phi_j) in w variables (check of the strong ell_j)
    Lw = []
    for f in phis:
        Xj = (-(p1**2 - p2) + b2) * f[0] + 2 * p1 * f[1] - f[2]
        loc = []
        for e in range(F.ne):
            B0, B1_, B2_ = F.Bl[e]; w = F.wq[e] * em[e]
            U2 = B2_ + p1[e] * B1_ + c[e] * B0; U1 = B1_ + p1[e] / 2 * B0
            loc.append(np.einsum('q,iq->i', w * Xj[e], U2) + np.einsum('q,iq->i', w * b1[e] * f[0][e], U1) + np.einsum('q,iq->i', w * b0[e] * f[0][e], B0))
        Lw.append(F.vec(np.array(loc)))
    Lw = np.array(Lw).T
    Phiv = np.array([F.load(em * f[0]) for f in phis]).T
    M3 = np.array([[np.sum(F.wq * em**2 * chi[0] * gi * f[0]) for f in phis] for gi in g])
    Mw = F.bilin({(0, 0): np.ones_like(x)}); K2 = F.bilin({(2, 2): np.ones_like(x)})
    return dict(F=F, KQ=KQ, Gv=Gv, Lv=Lv, Lw=Lw, Phiv=Phiv, M3=M3, Mw=Mw, K2=K2, ells=ells, C=C)

def phi_funcs(x, chi):
    """S1 constraint functions chi p_S1, chi n_S1, chi psi3 with 4 derivatives (float, finite-difference-free: splines of the jets)"""
    from scipy.interpolate import CubicSpline
    xf = np.linspace(0, 36, 36001)
    P, N = fem.s1_jets(xf); Ps = fem.psi3_jets(xf)
    out = []
    for base in (P[0], N[0], Ps[0]):
        sp = CubicSpline(xf, base)
        # derivatives 0..4 by spline differentiation of a 5th-order smoothing: use numpy gradient on fine grid (float check only)
        d = [base]
        for _ in range(4): d.append(np.gradient(d[-1], xf, edge_order=2))
        dj = [np.interp(x.ravel(), xf, dd).reshape(x.shape) for dd in d]
        out.append(leibniz(chi, dj))
    return out

def leibniz(chi, f):
    from math import comb
    return [sum(comb(n, k) * chi[k] * f[n - k] for k in range(n + 1)) for n in range(5)]

def certify_float(D, delta=0.0):
    """float version of the inertia certificate: In(K') + In(S) == (n+3, 6, 0)?"""
    Kp = D['KQ'] - delta * D['Mw']; n = Kp.shape[0]
    U = np.c_[D['Gv'], D['Lv']]; E = np.c_[U, D['Phiv']]
    Sh = np.zeros((9, 9)); Sh[0:3, 3:6] = -2 * D['M3']; Sh[3:6, 0:3] = -2 * D['M3'].T
    wK = np.linalg.eigvalsh(Kp); p = int(np.sum(wK < 0))
    S = Sh - E.T @ np.linalg.solve(Kp, E); wS = np.linalg.eigvalsh(0.5 * (S + S.T))
    return p, int(np.sum(wS > 0)), int(np.sum(wS < 0)), wK[:p + 2], wS

if __name__ == '__main__':
    g = bgeo.BaseGeom(0, 0.5)
    for th in (0.0, 0.03):
        D = build(g, th, 0.0, 0.1)
        Bm = D['KQ'] + 0.5 * (np.c_[D['Gv'], D['Lv']] @ np.block([[np.zeros((3, 3)), np.linalg.inv(D['M3']).T], [np.linalg.inv(D['M3']), np.zeros((3, 3))]]) @ np.c_[D['Gv'], D['Lv']].T)
        Z = sla.null_space(D['Phiv'].T); w = sla.eigh(Z.T @ Bm @ Z, Z.T @ D['Mw'] @ Z, eigvals_only=True)
        print(f"theta={th}: smallest eigenvalues of Q (= c(theta) - c_s): {np.round(w[:3], 5)}  n={Bm.shape[0]}  |Lv-Lw|/|Lw| = {np.linalg.norm(D['Lv']-D['Lw'])/np.linalg.norm(D['Lw']):.2e}")
        Bw = D['KQ'] + 0.5 * (np.c_[D['Gv'], D['Lw']] @ np.block([[np.zeros((3, 3)), np.linalg.inv(D['M3']).T], [np.linalg.inv(D['M3']), np.zeros((3, 3))]]) @ np.c_[D['Gv'], D['Lw']].T)
        w2 = sla.eigh(Z.T @ Bw @ Z, Z.T @ D['Mw'] @ Z, eigvals_only=True); print('     with weak-form l_j:', np.round(w2[:3], 5))
    # weak vs strong l_j: compare with the weak form of fem.FEM (rho-weighted) on the same functions
    for cs in (0.06, 0.065, 0.0689, 0.07, 0.075):
        D = build(g, 0.03, cs, 0.1); print('c_s', cs, 'certificate (p, n+, n-) of S:', certify_float(D)[:3])
