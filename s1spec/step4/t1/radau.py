"""Radau IIA collocation per slab: Y(tau0 + c dt) is the degree-s polynomial through Y0 and the s stage values,
Y'(tau0 + c_i dt) = F(Y_i).  Exactly continuous across slabs (c_s = 1).  Simplified Newton with the Jacobian frozen
per slab, block-diagonalised through the eigenvectors of A^{-1} (one complex LU per conjugate pair)."""
import numpy as np, scipy.linalg as sla
import mpmath as mp


def radau_coeffs(s, dps=60):
    mp.mp.dps = dps
    P = lambda n, x: mp.legendre(n, x)
    f = lambda x: P(s, x) - P(s - 1, x)
    # initial guesses from numpy, then refine
    from numpy.polynomial import legendre as LG
    cs = np.zeros(s + 1); cs[s] = 1; cs[s - 1] = -1
    xr = np.sort(LG.legroots(cs).real)
    xs = [mp.findroot(f, mp.mpf(float(x))) if abs(x - 1) > 1e-12 else mp.mpf(1) for x in xr]
    c = [(1 + x) / 2 for x in xs]
    A = mp.matrix(s, s)
    for j in range(s):
        others = [c[k] for k in range(s) if k != j]
        den = mp.fprod([c[j] - o for o in others])
        # integrate prod (t - o) from 0 to c_i: expand polynomial coefficients
        coef = [mp.mpf(1)]
        for o in others:
            new = [mp.mpf(0)] * (len(coef) + 1)
            for i, a in enumerate(coef):
                new[i + 1] += a; new[i] -= a * o
            coef = new
        for i in range(s):
            A[i, j] = mp.fsum([a * c[i] ** (k + 1) / (k + 1) for k, a in enumerate(coef)]) / den
    return c, A


class Radau:
    def __init__(self, s):
        self.s = s
        cm, Am = radau_coeffs(s)
        self.cm, self.Am = cm, Am
        self.c = np.array([float(x) for x in cm]); self.A = np.array([[float(Am[i, j]) for j in range(s)] for i in range(s)])
        self.Ainv = np.linalg.inv(self.A)
        lam, T = np.linalg.eig(self.Ainv)
        o = np.argsort(lam.imag); lam, T = lam[o], T[:, o]
        self.lam, self.T, self.Ti = lam, T, np.linalg.inv(T)
        # representative per conjugate pair (imag >= 0), partner index for imag < 0
        self.rep = [k for k in range(s) if lam[k].imag >= -1e-12]
        self.partner = {}
        for k in range(s):
            if lam[k].imag < -1e-12:
                j = int(np.argmin(np.abs(lam - np.conj(lam[k])))); self.partner[k] = j
        # nodes for the full polynomial: 0 and c
        self.nodes = np.r_[0.0, self.c]

    def factor(self, J, dt):
        n = J.shape[0]; self.lus = {}
        for k in self.rep:
            Mk = (self.lam[k] / dt) * np.eye(n) - J
            if abs(self.lam[k].imag) < 1e-12: Mk = Mk.real
            self.lus[k] = sla.lu_factor(Mk)
        self.dt = dt

    def solve(self, F, Y0, dt, Z0=None, tol=1e-14, maxit=40, scale=None, verbose=False):
        s = self.s; n = len(Y0)
        Z = np.zeros((s, n)) if Z0 is None else Z0.copy()
        sc = np.ones(n) if scale is None else scale
        hist = []
        for it in range(maxit):
            Fz = np.array([F(Y0 + Z[i]) for i in range(s)])
            Rz = -(self.Ainv / dt) @ Z + Fz
            Rw = self.Ti @ Rz
            dW = np.empty((s, n), complex)
            for k in self.rep: dW[k] = sla.lu_solve(self.lus[k], Rw[k])
            for k, j in self.partner.items(): dW[k] = np.conj(dW[j])
            dZ = (self.T @ dW).real
            Z += dZ
            err = np.max(np.abs(dZ) / sc)
            hist.append(err)
            if verbose: print(f"    it {it} |dZ| = {err:.2e}", flush=True)
            if err < tol or (it > 4 and err > 0.5 * hist[-2] and err < 30 * tol): break
        return Y0 + Z, hist

    def interp_matrix(self, theta):
        """Lagrange interpolation weights on nodes (0, c) at relative times theta (array): W (len(theta), s+1)"""
        nd = self.nodes; th = np.atleast_1d(theta)
        W = np.ones((len(th), len(nd)))
        for j in range(len(nd)):
            for k in range(len(nd)):
                if k != j: W[:, j] *= (th - nd[k]) / (nd[j] - nd[k])
        return W

    def deriv_matrix(self, theta):
        """d/dtheta of the Lagrange basis at theta (per unit of relative time)"""
        nd = self.nodes; th = np.atleast_1d(theta); m = len(nd)
        Wd = np.zeros((len(th), m))
        for j in range(m):
            for l in range(m):
                if l == j: continue
                term = np.ones(len(th)) / (nd[j] - nd[l])
                for k in range(m):
                    if k != j and k != l: term *= (th - nd[k]) / (nd[j] - nd[k])
                Wd[:, j] += term
        return Wd
