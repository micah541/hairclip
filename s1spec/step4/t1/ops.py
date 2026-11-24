"""T1 float operator.  Curvature form of the rescaled CDF on the half curve (arclength from the tip, tip fixed at -DS1):
   kappa_tau|_x = Phi_x,  Phi = theta_tau = -k''' - k^2 k' + k Q - (1/4 - a) x k,   Q(x) = int_0^x k'^2,
   a = 1/4 + (k''''(0) + K0^2 k''(0))/K0   (keeps kappa(0) = K0),
   ell_tau = -Q(ell) + (1/4 - a) ell,   (ln lam)_tau = a - 1/4,   b = k''(0) + (1/4 - a) DS1 (tip fixed), c_tau = -b lam, t_tau = lam^4.
Test space (Petrov-Galerkin): Leg_0..Leg_{P-6} on every element (moments; flux form, exact quadrature) and, at every node
j >= 1, the trial basis function of the nodal 4th derivative (Galerkin row; removes the zigzag near-null mode of pure moments).
With Leg_0 in the test space, int kappa is conserved exactly; I^1..I^(P-5) of the kappa-residual vanish at element ends.
ALE: knots x_j(ell), d/dtau|_x = d/dtau|_u - kappa_x (x_e' + u h_e') ell_tau  (+ nodal rescaling term)."""
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from rep import Rep, Mesh, K0, DS1

NS = 4          # extra scalars: ell, ln lam, c, t


class Dyn:
    def __init__(self, R):
        self.R = R; B = R.B; self.B = B
        N = R.N; self.n = R.n; self.ntot = R.n + NS
        self.nm = nm = B.nmom - 1                                   # Legendre moments per element
        self.MV = B.Mom[:nm] @ B.V[0]                                # (nm, nb)
        self.G0 = (B.V[0].T * B.wq[None, :]) @ B.V[0]                # local Gram (nb, nb)
        # row index: element e moments -> e*(nm+1) + m ; Galerkin row of node j>=1 -> (j-1)*(nm+1) + nm
        r, c, kind, e_, a_, q_ = [], [], [], [], [], []
        for e in range(N):
            for m in range(nm):
                for q in range(B.nb):
                    g = R.L[e, q]
                    if g >= 0: r.append(e * (nm + 1) + m); c.append(g); kind.append(0); e_.append(e); a_.append(m); q_.append(q)
        for j in range(1, N + 1):
            for e, qt in ((j - 1, 9), (j, 4)):
                if e >= N: continue
                for q in range(B.nb):
                    g = R.L[e, q]
                    if g >= 0: r.append((j - 1) * (nm + 1) + nm); c.append(g); kind.append(1); e_.append(e); a_.append(qt); q_.append(q)
        self.mr, self.mc = np.array(r), np.array(c)
        self.mk, self.me, self.ma, self.mq = np.array(kind), np.array(e_), np.array(a_), np.array(q_)
        self.mj = self.mr // (nm + 1)                                # node index - 1 for Galerkin rows
        self.i02 = R.nidx[0, 2]; self.i04 = R.nidx[0, 4]

    def hn(self, h):
        """row normaliser of the Galerkin row of node j (j = 1..N): mean length of the adjacent elements"""
        return np.r_[0.5 * (h[:-1] + h[1:]), h[-1]]

    def mass(self, ell):
        x, dx = self.R.mesh.knots(ell); h = np.diff(x)
        sc = self.R.lscale(h); hn = self.hn(h)
        e, a, q = self.me, self.ma, self.mq
        v0 = self.MV[np.minimum(a, self.nm - 1), q] * sc[e, q]
        v1 = h[e] * sc[e, a] * self.G0[a, q] * sc[e, q] / hn[self.mj]
        vals = np.where(self.mk == 0, v0, v1)
        return sps.csc_matrix((vals, (self.mr, self.mc)), shape=(self.n, self.n))

    def parts(self, Y):
        """everything at quadrature points; Y = [y (n), ell, lnlam, c, t]"""
        R, B = self.R, self.B
        y = Y[:self.n]; ell = Y[self.n]
        x, dx = R.mesh.knots(ell); h = np.diff(x); dh = np.diff(dx)
        d = R.local(y, h)
        H = h[:, None]
        k0 = d @ B.V[0].T; k1 = d @ B.V[1].T / H; k2 = d @ B.V[2].T / H**2; k3 = d @ B.V[3].T / H**3
        f = k1 * k1
        Qloc = H * (f @ B.Icum.T)
        Qe = np.r_[0.0, np.cumsum(h * (f @ B.wq))]
        Q = Qe[:-1, None] + Qloc
        H0 = h[0]
        k2t, k4t = y[self.i02] / H0**2 / B.nrm[2], y[self.i04] / H0**4 / B.nrm[4]
        a = 0.25 + (k4t + K0**2 * k2t) / K0
        xq = x[:-1, None] + H * B.uq[None, :]
        Phi = -k3 - k0 * k0 * k1 + k0 * Q - (0.25 - a) * xq * k0
        Dn = R.nodal(y, h)
        PhiN = -Dn[:, 3] - Dn[:, 0]**2 * Dn[:, 1] + Dn[:, 0] * Qe - (0.25 - a) * x * Dn[:, 0]
        elld = -Qe[-1] + (0.25 - a) * ell
        # ALE generator per unit ell_dot: sum_m m (h'/h - H'/H) d_m N_m - k1 (x_e' + u h_e')
        Hn, dHn = R.nodal_scale(h, dh)
        mm = np.arange(5)[None, :]
        dA = np.zeros_like(d)
        dA[:, :5] = d[:, :5] * mm * (dh / h - dHn[:-1] / Hn[:-1])[:, None]
        dA[:, 5:10] = d[:, 5:10] * mm * (dh / h - dHn[1:] / Hn[1:])[:, None]
        ALE = dA @ B.V[0].T - k1 * (dx[:-1, None] + B.uq[None, :] * dh[:, None])
        return dict(k0=k0, k1=k1, k2=k2, k3=k3, k2t=k2t, k4t=k4t, x=x, dx=dx, h=h, dh=dh, d=d, Q=Q, Qe=Qe, a=a, Phi=Phi, PhiN=PhiN,
                    elld=elld, ALE=ALE, xq=xq, Dn=Dn)

    def rhs(self, p, elld):
        """test-function integrals of Phi_x - elld * ALE (rows as in mass())"""
        R, B = self.R, self.B; N, nm = R.N, self.nm; h = p['h']; sc = R.lscale(h); hn = self.hn(h)
        out = np.zeros((N, nm + 1))
        GM = (B.mfac[None, :nm] * (p['PhiN'][1:, None] - B.sgn[None, :nm] * p['PhiN'][:-1, None]) - p['Phi'] @ B.MomD[:nm].T) / h[:, None]
        out[:, :nm] = GM - elld * (p['ALE'] @ B.Mom[:nm].T)
        # Galerkin rows: node j = e+1 gets the right-end function (local 9) of element e and the left-end function (local 4) of element e+1
        Wq = B.wq[None, :]
        PV1 = (p['Phi'] * Wq) @ B.V[1]            # (N, nb): int Phi dN_q/du du
        AV0 = (p['ALE'] * Wq) @ B.V[0]
        gR = sc[:, 9] * (-PV1[:, 9] - elld * h * AV0[:, 9])
        gL = sc[:, 4] * (-PV1[:, 4] - elld * h * AV0[:, 4])
        g = gR.copy(); g[:-1] += gL[1:]
        out[:, nm] = g / hn
        return out.ravel()

    def F(self, Y, ret_parts=False):
        p = self.parts(Y)
        Mm = self.mass(Y[self.n])
        yd = spla.spsolve(Mm, self.rhs(p, p['elld']))
        lam = np.exp(Y[self.n + 1]); a = p['a']
        b = p['k2t'] + (0.25 - a) * DS1
        out = np.r_[yd, p['elld'], a - 0.25, -b * lam, lam**4]
        if ret_parts: p.update(b=b, Mm=Mm); return out, p
        return out

    def jac(self, Y, F0=None, eps=1e-7):
        F0 = self.F(Y) if F0 is None else F0
        n = self.ntot; J = np.empty((n, n))
        sc = np.maximum(np.abs(Y), 1e-3)
        for j in range(n):
            Yp = Y.copy(); dj = eps * sc[j]; Yp[j] += dj
            J[:, j] = (self.F(Yp) - F0) / dj
        return J
