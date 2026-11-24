"""T7 float FEM for the projected form (C1) of T5 2.4 (rho-weighted, C^1 cubic Hermite, even class u'(0) = 0).
Domain [0, X]: X = 34.25 with u(X) = u'(X) = 0 ('clamped': the class of (C1), u supported in [0, 34.25]) or free end.
  a(u, v) = <L u, v>_rho = int u''[(-rho'' + rho b2) v - 2 rho' v' - rho v''] + rho b1 v u' + rho b0 v u
  B(u) = a(u, u) - sum_ij <chi g_i, u> (M^-T)_ij a(u, phi_j),  M_ij = <chi g_i, phi_j>,  on {<u, phi_j> = 0}.
Constraints phi = (chi p_S1, chi n_S1, chi psi3^) (T5 D2), directions g = (p~, n~, V~) of the base."""
import sys, numpy as np, scipy.linalg as sla
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, '/home/claude/s1spec/step4/feas'); sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_int/v3')
import bgeo

def rho_jets(x):
    g = (1 + x**2) ** (2 / 3); g1 = (4 / 3) * x * (1 + x**2) ** (-1 / 3)
    g2 = (4 / 3) * (1 + x**2) ** (-1 / 3) - (8 / 9) * x**2 * (1 + x**2) ** (-4 / 3)
    r = np.exp(-g / 4); return r, -g1 / 4 * r, (-g2 / 4 + g1**2 / 16) * r

def hermite(xi, h):
    t = xi; t2, t3 = t * t, t ** 3
    B0 = np.array([2*t3 - 3*t2 + 1, h * (t3 - 2*t2 + t), -2*t3 + 3*t2, h * (t3 - t2)])
    B1 = np.array([6*t2 - 6*t, h * (3*t2 - 4*t + 1), -6*t2 + 6*t, h * (3*t2 - 2*t)]) / h
    B2 = np.array([12*t - 6, h * (6*t - 4), -12*t + 6, h * (6*t - 2)]) / h**2
    return B0, B1, B2

class FEM:
    def __init__(self, X=34.25, h=0.1, nq=6, end='clamped'):
        self.ne = int(round(X / h)); self.h = X / self.ne; self.X = X
        gx, gw = leggauss(nq); xi = (gx + 1) / 2
        self.xq = np.arange(self.ne)[:, None] * self.h + self.h * xi[None, :]
        self.wq = np.tile(gw / 2 * self.h, (self.ne, 1)); self.B = hermite(xi, self.h)
        self.ndof = 2 * (self.ne + 1)
        drop = [1] + ([self.ndof - 2, self.ndof - 1] if end == 'clamped' else [])
        self.free = np.array([i for i in range(self.ndof) if i not in drop])
        self.rho, self.rho1, self.rho2 = rho_jets(self.xq)
        self.idx = 2 * np.arange(self.ne)[:, None] + np.arange(4)[None, :]
    def scatter(self, loc):
        A = np.zeros((self.ndof, self.ndof))
        for e in range(self.ne): A[np.ix_(self.idx[e], self.idx[e])] += loc[e]
        return A[np.ix_(self.free, self.free)]
    def vec(self, loc):
        v = np.zeros(self.ndof)
        for e in range(self.ne): v[self.idx[e]] += loc[e]
        return v[self.free]
    def mass(self):
        w = self.wq * self.rho; B0 = self.B[0]; return self.scatter(np.einsum('eq,iq,jq->eij', w, B0, B0))
    def stiff2(self):
        w = self.wq * self.rho; B2 = self.B[2]; return self.scatter(np.einsum('eq,iq,jq->eij', w, B2, B2))
    def form(self, b2, b1, b0):
        """K[k, l] = a(psi_l, psi_k)"""
        B0, B1, B2 = self.B; w = self.wq
        cA = w * (-self.rho2 + self.rho * b2); cB = w * (-2 * self.rho1); cC = -w * self.rho
        cD = w * self.rho * b1; cE = w * self.rho * b0
        loc = (np.einsum('eq,iq,jq->eij', cA, B0, B2) + np.einsum('eq,iq,jq->eij', cB, B1, B2) + np.einsum('eq,iq,jq->eij', cC, B2, B2)
               + np.einsum('eq,iq,jq->eij', cD, B0, B1) + np.einsum('eq,iq,jq->eij', cE, B0, B0))
        return self.scatter(loc)
    def load(self, f):
        """<f, psi_k>_rho"""
        return self.vec(np.einsum('eq,iq->ei', self.wq * self.rho * f, self.B[0]))
    def aload(self, b2, b1, b0, f, f1, f2):
        """a(psi_l, f) for every basis function psi_l (f with derivatives f1, f2)"""
        B0, B1, B2 = self.B; w = self.wq
        X = (-self.rho2 + self.rho * b2) * f - 2 * self.rho1 * f1 - self.rho * f2
        return self.vec(np.einsum('eq,iq->ei', w * X, B2) + np.einsum('eq,iq->ei', w * self.rho * b1 * f, B1) + np.einsum('eq,iq->ei', w * self.rho * b0 * f, B0))
    def ip(self, f, g):
        return float(np.sum(self.wq * self.rho * f * g))

_S1 = None
def s1_jets(x):
    """S1 constraint functions with two derivatives: p, n = cos(int kappa) (translation mode), on points x"""
    from f3geom import s1_geometry
    from scipy.interpolate import CubicSpline
    global _S1
    if _S1 is None:
        xg = np.arange(0, 60.0005, 0.0025); g = s1_geometry(xg)
        Th = CubicSpline(xg, g['kap']).antiderivative()(xg)
        _S1 = {k: CubicSpline(xg, v) for k, v in dict(p=g['p'], q=g['q'], k=g['kap'], m=g['kap1'], Th=Th).items()}
    p, q, k, m, Th = [_S1[c](x) for c in ('p', 'q', 'k', 'm', 'Th')]
    P = (p, -k * q, -m * q - k * (1 + k * p))
    N = (np.cos(Th), -k * np.sin(Th), -m * np.sin(Th) - k**2 * np.cos(Th))
    return P, N

def psi3_jets(x):
    import c4fun
    from flint import arb
    P = c4fun.PSI; out = np.zeros((3,) + x.shape); xf = x.ravel(); o = out.reshape(3, -1)
    ends = [float(p.b.mid()) for p in P.pieces]; XC = float(P.M.XC.mid())
    for j, xv in enumerate(xf):
        if xv <= XC:
            e = min(int(np.searchsorted(ends, xv)), len(ends) - 1); jj = P.piece_jets(e, arb(float(xv)), 3)
        else: jj = c4fun.ext_jets(P.M, arb(float(xv)), 3)
        o[:, j] = [float(v.mid()) for v in jj]
    return out

def cut_prod(c, f):
    """jets 0..2 of c * f from jets of c and f"""
    return (c[0] * f[0], c[1] * f[0] + c[0] * f[1], c[2] * f[0] + 2 * c[1] * f[1] + c[0] * f[2])

class Setup:
    """FEM + fixed constraint data (independent of time)"""
    def __init__(self, X=34.25, h=0.1, end='clamped', nq=6):
        self.fem = F = FEM(X, h, nq, end); x = F.xq
        ch = [bgeo.chi32(x, n) for n in range(3)]
        P, N = s1_jets(x); Ps = psi3_jets(x)
        self.phi = [cut_prod(ch, P), cut_prod(ch, N), cut_prod(ch, Ps)]
        self.chi = ch[0]
        self.M = F.mass(); self.K2 = F.stiff2()
        self.Phiv = np.array([F.load(f[0]) for f in self.phi]).T

def operator(st, geom):
    """K, Gv, Lv, M3 of the base geometry (BaseGeom) on the setup st"""
    F = st.fem; x = F.xq; s = geom.at(x.ravel(), mmax=3)
    sh = lambda a: a.reshape(x.shape)
    cf = geom.coeffs(s); b2, b1, b0 = sh(cf['b2']), sh(cf['b1']), sh(cf['b0'])
    K = F.form(b2, b1, b0)
    g = [sh(s['p']), sh(s['n']), sh(s['Vs'])]
    Gv = np.array([F.load(st.chi * gi) for gi in g]).T
    Lv = np.array([F.aload(b2, b1, b0, *f) for f in st.phi]).T
    M3 = np.array([[F.ip(st.chi * gi, f[0]) for f in st.phi] for gi in g])          # M_ij = <chi g_i, phi_j>
    return dict(K=K, Gv=Gv, Lv=Lv, M3=M3, b=(b2, b1, b0), g=g, s=s)

def projected(st, op, Phiv=None):
    """matrix of B on the FEM space, and the constraint matrix"""
    Phiv = st.Phiv if Phiv is None else Phiv
    Bm = op['K'] - op['Gv'] @ np.linalg.inv(op['M3']).T @ op['Lv'].T
    return Bm, Phiv

def top_eigs(st, Bm, Phiv, theta=0.0, nev=3, M=None):
    M = st.M if M is None else M
    Z = sla.null_space(Phiv.T)
    A = Z.T @ (Bm + theta * st.K2) @ Z; A = 0.5 * (A + A.T)
    w, V = sla.eigh(A, Z.T @ M @ Z)
    return w[::-1][:nev], Z @ V[:, ::-1][:, :nev]
