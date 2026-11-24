"""C^1 cubic-Hermite Galerkin discretization (even class: u'(0)=0 essential, u'''(0)=0 natural; free at sigma=Sigma) of the rho-weighted
form a(u,v) = <A u, v>_rho,   A u = -u'''' + b2 u'' + b1 u' + b0 u,
   b2 = -kap^2,  b1 = -5 kap kap' - q/4 + tau,  b0 = 1/4 - 3 kap'^2 - 4 kap kap''     (see FEM_NOTES in report; tau = 0 for S1)
Weak form: a(u,v) = int u'' [(-rho'' + rho b2) v - 2 rho' v' - rho v''] + rho b1 v u' + rho b0 v u."""
import numpy as np, scipy.linalg as sla
from numpy.polynomial.legendre import leggauss

def rho_fun(x):
    g = (1 + x**2) ** (2 / 3); g1 = (4 / 3) * x * (1 + x**2) ** (-1 / 3)
    g2 = (4 / 3) * (1 + x**2) ** (-1 / 3) - (8 / 9) * x**2 * (1 + x**2) ** (-4 / 3)
    r = np.exp(-g / 4); return r, -g1 / 4 * r, (-g2 / 4 + g1**2 / 16) * r

def hermite_basis(xi, h):
    t = xi; t2, t3 = t * t, t * t * t
    B0 = np.array([2*t3 - 3*t2 + 1, h * (t3 - 2*t2 + t), -2*t3 + 3*t2, h * (t3 - t2)])
    B1 = np.array([6*t2 - 6*t, h * (3*t2 - 4*t + 1), -6*t2 + 6*t, h * (3*t2 - 2*t)]) / h
    B2 = np.array([12*t - 6, h * (6*t - 4), -12*t + 6, h * (6*t - 2)]) / h**2
    return B0, B1, B2

class Mesh:
    def __init__(self, Sigma, h, nq=5):
        self.ne = int(round(Sigma / h)); self.h = Sigma / self.ne; self.Sigma = Sigma; self.nq = nq
        gx, gw = leggauss(nq); xi = (gx + 1) / 2
        self.xe = np.arange(self.ne) * self.h
        self.xq = self.xe[:, None] + self.h * xi[None, :]            # (ne,nq)
        self.wq = np.tile(gw / 2 * self.h, (self.ne, 1))
        self.B = hermite_basis(xi, self.h)
        self.ndof = 2 * (self.ne + 1)
        self.free = np.r_[0, np.arange(2, self.ndof)]              # drop slope dof at node 0
        r, r1, r2 = rho_fun(self.xq); self.rho, self.rho1, self.rho2 = r, r1, r2
    def scatter(self, loc):
        ne = self.ne; A = np.zeros((self.ndof, self.ndof))
        for e in range(ne):
            ix = np.arange(2 * e, 2 * e + 4); A[np.ix_(ix, ix)] += loc[e]
        return A[np.ix_(self.free, self.free)]
    def mass(self):
        B0 = self.B[0]; w = self.wq * self.rho
        return self.scatter(np.einsum('eq,iq,jq->eij', w, B0, B0))
    def stiff2(self):                     # int rho u'' v''
        B2 = self.B[2]; w = self.wq * self.rho
        return self.scatter(np.einsum('eq,iq,jq->eij', w, B2, B2))
    def form(self, cf):
        """cf: dict of coefficient values at the quadrature points (ne,nq): b2,b1,b0.  Returns K[i,j] = a(phi_j, phi_i)."""
        B0, B1, B2 = self.B; w = self.wq
        cA = w * (-self.rho2 + self.rho * cf['b2']); cB = w * (-2 * self.rho1); cC = -w * self.rho
        cD = w * self.rho * cf['b1']; cE = w * self.rho * cf['b0']
        loc = (np.einsum('eq,iq,jq->eij', cA, B0, B2) + np.einsum('eq,iq,jq->eij', cB, B1, B2) + np.einsum('eq,iq,jq->eij', cC, B2, B2)
               + np.einsum('eq,iq,jq->eij', cD, B0, B1) + np.einsum('eq,iq,jq->eij', cE, B0, B0))
        return self.scatter(loc)

def coeffs(mesh, sp, use_tau=True):
    x = mesh.xq; k, k1, k2, q, tau = sp['kap'](x), sp['kap1'](x), sp['kap2'](x), sp['q'](x), sp['tau'](x)
    return dict(b2=-k**2, b1=-5 * k * k1 - q / 4 + (tau if use_tau else 0.0), b0=0.25 - 3 * k1**2 - 4 * k * k2)

def eig_nonsym(K, M, nev=12):
    d = 1 / np.sqrt(np.diag(M)); Ks = K * d[:, None] * d[None, :]; Ms = M * d[:, None] * d[None, :]
    w, vr = sla.eig(Ks, Ms)
    o = np.argsort(-w.real); return w[o], (vr * d[:, None])[:, o]

def eig_sym(K, M, G=None):
    """eigenvalues/vectors of the symmetric part (K+K^T)/2 relative to the inner product matrix G (default M);
    returns M-orthonormal... vectors normalised w.r.t. G."""
    S = 0.5 * (K + K.T); G = M if G is None else G
    d = 1 / np.sqrt(np.diag(G)); Ss = S * d[:, None] * d[None, :]; Gs = G * d[:, None] * d[None, :]
    w, v = sla.eigh(Ss, Gs)
    o = np.argsort(-w); return w[o], (v * d[:, None])[:, o]
