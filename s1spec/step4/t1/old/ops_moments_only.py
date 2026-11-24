"""T1 float operator.  Curvature form of the rescaled CDF on the half curve (arclength from the tip, tip fixed at -DS1):
   kappa_tau|_x = -k'''' - k^2 k'' - k k'^2 + k' Q - (1/4 - a)(k + x k'),   Q(x) = int_0^x k'^2,
   a = 1/4 + (k''''(0) + K0^2 k''(0))/K0   (keeps kappa(0) = K0),
   ell_tau = -Q(ell) + (1/4 - a) ell,   (ln lam)_tau = a - 1/4,   b = k''(0) + (1/4 - a) DS1 (tip fixed), c_tau = -b lam, t_tau = lam^4.
Petrov-Galerkin: residual tested against Leg_0..Leg_{P-5} on every element (moments exact by Gauss quadrature).
ALE: knots x_j(ell), d/dtau|_x = d/dtau|_u - kappa_x (x_e' + u h_e') ell_tau  (+ nodal rescaling term)."""
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from rep import Rep, Mesh, K0, DS1

NS = 4          # extra scalars: ell, ln lam, c, t


class Dyn:
    def __init__(self, R):
        self.R = R; B = R.B; self.B = B
        self.MV = [B.Mom @ B.V[k] for k in range(7)]                 # (nmom, nb): moments of basis derivatives
        self.n = R.n; self.ntot = R.n + NS
        N, nm = R.N, B.nmom
        # sparsity pattern of the mass matrix
        rows, cols, loc = [], [], []
        for e in range(N):
            for m in range(nm):
                for q in range(B.nb):
                    g = R.L[e, q]
                    if g >= 0: rows.append(e * nm + m); cols.append(g); loc.append((e, m, q))
        self.mrows = np.array(rows); self.mcols = np.array(cols)
        loc = np.array(loc); self.me, self.mm, self.mq = loc[:, 0], loc[:, 1], loc[:, 2]
        self.i02 = R.nidx[0, 2]; self.i04 = R.nidx[0, 4]

    def mass(self, ell):
        x, dx = self.R.mesh.knots(ell); h = np.diff(x)
        sc = self.R.lscale(h)
        vals = self.MV[0][self.mm, self.mq] * sc[self.me, self.mq]
        n = self.n
        return sps.csc_matrix((vals, (self.mrows, self.mcols)), shape=(n, n))

    def parts(self, Y):
        """returns dict with everything at quadrature points; Y = [y (n), ell, lnlam, c, t]"""
        R, B = self.R, self.B
        y = Y[:self.n]; ell = Y[self.n]
        x, dx = R.mesh.knots(ell); h = np.diff(x); dh = np.diff(dx)
        d = R.local(y, h)
        H = h[:, None]
        k0 = d @ B.V[0].T; k1 = d @ B.V[1].T / H; k2 = d @ B.V[2].T / H**2; k3 = d @ B.V[3].T / H**3; k4 = d @ B.V[4].T / H**4
        f = k1 * k1
        Qloc = H * (f @ B.Icum.T)
        Qe = np.r_[0.0, np.cumsum(h * (f @ B.wq))]
        Q = Qe[:-1, None] + Qloc
        H0 = h[0]
        k2t, k4t = y[self.i02] / H0**2 / B.nrm[2], y[self.i04] / H0**4 / B.nrm[4]
        a = 0.25 + (k4t + K0**2 * k2t) / K0
        xq = x[:-1, None] + H * B.uq[None, :]
        G = -k4 - k0 * k0 * k2 - k0 * f + k1 * Q - (0.25 - a) * (k0 + xq * k1)
        # flux Phi = theta_tau = -k3 - k^2 k1 + k Q - (1/4-a) x k ;  G = d Phi / dx
        Phi = -k3 - k0 * k0 * k1 + k0 * Q - (0.25 - a) * xq * k0
        Dn = R.nodal(y, h)
        PhiN = -Dn[:, 3] - Dn[:, 0]**2 * Dn[:, 1] + Dn[:, 0] * Qe - (0.25 - a) * x * Dn[:, 0]
        GM = (B.mfac[None, :] * (PhiN[1:, None] - B.sgn[None, :] * PhiN[:-1, None]) - Phi @ B.MomD.T) / H   # moments of G (flux form)
        elld = -Qe[-1] + (0.25 - a) * ell
        # ALE generator per unit ell_dot: sum_m m (h'/h) d_m H_m - k1 (x_e' + u h_e')
        Hn, dHn = R.nodal_scale(h, dh)
        mm = np.arange(5)[None, :]
        dA = np.zeros_like(d)
        dA[:, :5] = d[:, :5] * mm * (dh / h - dHn[:-1] / Hn[:-1])[:, None]
        dA[:, 5:10] = d[:, 5:10] * mm * (dh / h - dHn[1:] / Hn[1:])[:, None]
        ALE = dA @ B.V[0].T - k1 * (dx[:-1, None] + B.uq[None, :] * dh[:, None])
        return dict(GM=GM, Phi=Phi, PhiN=PhiN, k3=k3, k2t=k2t, x=x, dx=dx, h=h, d=d, k0=k0, k1=k1, k2=k2, k4=k4, Q=Q, Qe=Qe, a=a, G=G, elld=elld, ALE=ALE, xq=xq)

    def F(self, Y, ret_parts=False):
        p = self.parts(Y)
        Mm = self.mass(Y[self.n])
        rhs = (p['GM'] - p['elld'] * (p['ALE'] @ self.B.Mom.T)).ravel()      # (N*nmom)
        yd = spla.spsolve(Mm, rhs)
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

    def residual_moments(self, Y, Yd):
        """moment residual of the PDE for a given state and time derivative (diagnostic)"""
        p = self.parts(Y); Mm = self.mass(Y[self.n])
        rhs = (p['GM'] - Yd[self.n] * (p['ALE'] @ self.B.Mom.T)).ravel()
        return Mm @ Yd[:self.n] - rhs
