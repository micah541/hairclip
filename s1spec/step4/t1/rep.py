"""T1 representation: kappa(x, tau) on the half curve x in [0, ell(tau)] (arclength from the tip, rescaled frame).
Per element: polynomial of degree P in u in [0,1] (x = x_e + h_e u) = Hermite-C^4 part (nodal kappa^(m)(x_j), m = 0..4)
+ (P-9) bubbles u^5 (1-u)^5 Leg_i(2u-1).  So kappa is a C^4 piecewise polynomial of degree P.
Mesh: tip zone [0, XT] fixed knots (density 1/(h0 + cT x)); outer zone [XT, ell] knots x_j(ell) from
   (1/c) ln(x/XT) + (x - XT)/(beta ell) = K(ell) s_j   (Lambert W; geometric, then uniform in physical arclength).
Constraints built in: kappa(0) = K0, kappa'(0) = kappa'''(0) = kappa'(ell) = kappa'''(ell) = 0 (exact symmetry at both ends)."""
import numpy as np
from math import comb, factorial
from numpy.polynomial import legendre as LG
from numpy.polynomial import chebyshev as CH
from scipy.special import lambertw

K0 = 1.5398206378719341          # S1 tip curvature (evans.PROF at the tip)
DS1 = 2.500806276619802          # S1 tip at -DS1 (p(0) = -DS1, n(0) = 1, theta(0) = -pi/2)
STIP = -2.249077651121112        # y = x + STIP


def mono_mul(a, b): return np.convolve(a, b)


def hermite_local(P):
    """monomial coefficient matrix C (P+1, nb): columns = H0m (m=0..4), H1m (m=0..4), bubbles (P-9)."""
    cols = []
    om5 = np.array([comb(5, k) * (-1) ** k for k in range(6)], float)        # (1-u)^5
    for m in range(5):
        s = np.zeros(5 - m); s[:] = [comb(4 + k, k) for k in range(5 - m)]
        c = mono_mul(np.r_[np.zeros(m), 1.0 / factorial(m)], mono_mul(om5, s))
        cols.append(np.pad(c, (0, P + 1 - len(c))))
    for m in range(5):     # H1m(u) = (-1)^m H0m(1-u)
        c0 = cols[m]; c = np.zeros(P + 1)
        # expand sum_k c0[k] (1-u)^k
        for k, ck in enumerate(c0):
            if ck == 0: continue
            for i in range(k + 1): c[i] += ck * comb(k, i) * (-1) ** i
        cols.append((-1) ** m * c)
    b5 = mono_mul(np.r_[np.zeros(5), 1.0], om5)                                  # u^5 (1-u)^5
    for i in range(P - 9):
        li = LG.leg2poly(np.r_[np.zeros(i), 1.0])                               # Leg_i(t), t = 2u-1
        # substitute t = 2u - 1
        c = np.zeros(i + 1)
        for k, lk in enumerate(li):
            for j in range(k + 1): c[j] += lk * comb(k, j) * 2.0 ** j * (-1.0) ** (k - j)
        bc = mono_mul(b5, c); cols.append(np.pad(bc, (0, P + 1 - len(bc))))
    return np.array(cols).T


def cheb_exact(mono):
    """exact Chebyshev coefficients (rounded to float) of sum_k mono[k] u^k on u in [0,1], t = 2u - 1"""
    from fractions import Fraction
    n = len(mono)
    # p(t) coefficients: u = (t+1)/2
    pt = [Fraction(0)] * n
    for k, ck in enumerate(mono):
        if ck == 0: continue
        for i in range(k + 1): pt[i] += ck * comb(k, i) / Fraction(2) ** k
    # power series in t -> Chebyshev: t^k = 2^(1-k) sum_j C(k, j) T_|k-2j| (halved for j = k/2)
    ch = [Fraction(0)] * n
    for k, a in enumerate(pt):
        if a == 0: continue
        for j in range(k + 1):
            m = abs(k - 2 * j); w = Fraction(comb(k, j), 2 ** k)
            ch[m] += a * w
    return np.array([float(v) for v in ch])


def mono_der(c, k):
    for _ in range(k): c = np.polynomial.polynomial.polyder(c)
    return c


class Basis:
    """local basis values at Gauss points, derivatives 0..6, moments against Leg_m(2u-1), m = 0..P-5"""
    def __init__(self, P, nq=None, kmax=6):
        self.P = P; self.nb = P + 1; self.nmom = P - 4; self.nbub = P - 9
        self.nq = nq or 2 * P
        g, w = LG.leggauss(self.nq); self.uq = (g + 1) / 2; self.wq = w / 2
        C = hermite_local(P)
        uu = np.linspace(0, 1, 2001); Vu = np.vander(uu, P + 1, increasing=True)
        nrm = np.max(np.abs(Vu @ C), axis=0)                 # unknowns are multiplied by nrm (basis normalised to max ~1)
        nrm[5:10] = nrm[:5]                                    # identical factor for left/right functions of equal order (exact C^4)
        self.nrm = 2.0 ** np.round(np.log2(nrm))               # dyadic: exact in float and arb
        self.C = C / self.nrm[None, :]
        # accurate float tables: exact rational Chebyshev series (t = 2u - 1) of the normalised basis, rounded once
        from arbops import hermite_exact
        from fractions import Fraction
        cols = hermite_exact(P); nrmF = [Fraction(float(v)) for v in self.nrm]
        self.cheb = np.array([cheb_exact([c / nrmF[j] for c in cols[j]]) for j in range(self.nb)]).T     # (P+1, nb)
        self.V = [self.evalpts(self.uq, k) for k in range(kmax + 1)]                    # (nq, nb)
        # test functions Leg_m(2u-1) at nodes, times weights -> moment operator (nmom, nq)
        self.Pm = np.array([LG.legval(2 * self.uq - 1, np.r_[np.zeros(m), 1.0]) for m in range(self.nmom)])
        self.Mom = self.Pm * self.wq[None, :] * (2 * np.arange(self.nmom) + 1)[:, None]   # Legendre coefficients of the residual
        # flux form: (2m+1) int_0^1 Phi * d/du Leg_m(2u-1) du  and end values Leg_m(1) = 1, Leg_m(-1) = (-1)^m
        self.dPm = np.array([2 * LG.legval(2 * self.uq - 1, LG.legder(np.r_[np.zeros(m), 1.0])) if m > 0 else np.zeros(self.nq) for m in range(self.nmom)])
        self.MomD = self.dPm * self.wq[None, :] * (2 * np.arange(self.nmom) + 1)[:, None]
        self.sgn = (-1.0) ** np.arange(self.nmom); self.mfac = 2 * np.arange(self.nmom) + 1.0
        # cumulative integration matrix on [0, u_q] exact for degree <= nq-1
        Vl = LG.legvander(2 * self.uq - 1, self.nq - 1); Vinv = np.linalg.inv(Vl)
        Icum = np.zeros((self.nq, self.nq))
        for j in range(self.nq):
            cl = LG.legint(np.eye(self.nq)[j], lbnd=-1) / 2.0                  # d/du = 2 d/dt
            Icum[:, j] = LG.legval(2 * self.uq - 1, cl)
        self.Icum = Icum @ Vinv
        # end values of basis functions and derivatives (u = 0, 1)
        self.E0 = [self.evalpts(np.array([0.0]), k)[0] for k in range(kmax + 1)]
        self.E1 = [self.evalpts(np.array([1.0]), k)[0] for k in range(kmax + 1)]

    def evalpts(self, u, k):
        t = 2 * np.atleast_1d(np.asarray(u, float)) - 1
        out = np.empty((len(t), self.nb))
        for j in range(self.nb):
            c = CH.chebder(self.cheb[:, j], k) * 2.0**k if k else self.cheb[:, j]
            out[:, j] = CH.chebval(t, c)
        return out


class Mesh:
    """tip zone: JT fixed elements on [0, XT] (h ~ h0 + cT x).
    outer zone: knots x_j(ell) = XT (1 + C_j ell^p)^(1/p) (j < N), x_N = ell:  dx/d ln(ell) = x (1 - (XT/x)^p),
    i.e. Lagrangian far out (x ~ ell), at rest near XT; elements are split (exact refinement) when too large.
    The initial outer knots come from (1/c) ln(x/XT) + (x - XT)/(beta ell) = K(ell) s (geometric, then uniform physical)."""
    PEXP = 3.0
    def __init__(self, JT=40, XT=40.0, h0=0.35, cT=0.08, NO=80, c=0.12, beta=0.02, ell0=None, C=None):
        self.JT, self.XT, self.h0, self.cT = JT, XT, h0, cT
        A = (1 + cT * XT / h0)
        self.xtip = (h0 / cT) * (A ** (np.arange(JT + 1) / JT) - 1); self.xtip[-1] = XT
        if C is None:
            xo = self.lambert(ell0, NO, c, beta)
            C = ((xo[1:-1] / XT) ** self.PEXP - 1) / ell0 ** self.PEXP
        self.C = np.asarray(C, float)
        self.N = JT + len(self.C) + 1

    def lambert(self, ell, NO, c, b):
        XT, s = self.XT, np.arange(NO + 1) / NO
        A, B = 1 / c, 1 / (b * ell)
        K = A * np.log(ell / XT) + B * (ell - XT)
        lnz = np.log(B / A * XT) + B * XT / A + K * s / A
        x = (A / B) * np.real(lambertw(np.exp(lnz))); x[0] = XT; x[-1] = ell
        return x

    def knots(self, ell):
        p = self.PEXP; XT = self.XT
        g = 1 + self.C * ell**p
        xo = XT * g ** (1 / p); dxo = XT * self.C * ell**(p - 1) * g ** (1 / p - 1)
        x = np.r_[self.xtip, xo, ell]; dx = np.r_[np.zeros(self.JT + 1), dxo, 1.0]
        return x, dx

    def inserted(self, xnew, ell):
        """new mesh with an outer knot at x = xnew (XT < xnew < x_{N-1}) inserted"""
        Cn = ((xnew / self.XT) ** self.PEXP - 1) / ell ** self.PEXP
        m = Mesh.__new__(Mesh); m.JT, m.XT, m.h0, m.cT, m.xtip = self.JT, self.XT, self.h0, self.cT, self.xtip
        m.C = np.sort(np.r_[self.C, Cn]); m.N = m.JT + len(m.C) + 1; return m


_BASIS = {}
def basis(P, nq=None):
    key = (P, nq)
    if key not in _BASIS: _BASIS[key] = Basis(P, nq)
    return _BASIS[key]


class Rep:
    """global DOF layout: nodal D[j, m] (j = 0..N, m = 0..4), bubbles b[e, i]; free vector y (without ell)."""
    def __init__(self, P=9, mesh=None, nq=None):
        self.B = basis(P, nq); self.P = P; self.mesh = mesh
        N = self.mesh.N; self.N = N; nbub = self.B.nbub
        fixed = {(0, 0): K0, (0, 1): 0.0, (0, 3): 0.0, (N, 1): 0.0, (N, 3): 0.0}
        self.fixed = fixed
        idx = -np.ones((N + 1, 5), int); k = 0
        bidx = np.zeros((N, nbub), int)
        for j in range(N + 1):                 # interleaved: node j, then bubbles of element j (banded mass matrix)
            for m in range(5):
                if (j, m) not in fixed: idx[j, m] = k; k += 1
            if j < N:
                bidx[j] = k + np.arange(nbub); k += nbub
        self.nidx = idx; self.bidx = bidx; self.n = k
        # local -> global map per element: 10 nodal + nbub (index -1 = fixed)
        L = np.zeros((N, 10 + nbub), int)
        for e in range(N):
            L[e, :5] = idx[e]; L[e, 5:10] = idx[e + 1]; L[e, 10:] = self.bidx[e]
        self.L = L
        self.fixval = np.zeros((N, 10 + nbub))
        for (j, m), v in fixed.items():
            if j < N: self.fixval[j, m] = v * self.B.nrm[m]
            if j > 0: self.fixval[j - 1, 5 + m] = v * self.B.nrm[5 + m]

    @staticmethod
    def nodal_scale(h, dh=None):
        """nodal length scale H_j = sqrt(h_{j-1} h_j) (one-sided at the ends) and dH_j/dell"""
        H = np.r_[h[0], np.sqrt(h[:-1] * h[1:]), h[-1]]
        if dh is None: return H
        g = dh / h
        dH = H * np.r_[g[0], 0.5 * (g[:-1] + g[1:]), g[-1]]
        return H, dH

    def lscale(self, h):
        """factor (h_e/H_j)^m for the nodal columns of every element"""
        H = self.nodal_scale(h); m = np.arange(5)
        sc = np.ones((self.N, self.B.nb))
        sc[:, :5] = (h / H[:-1])[:, None] ** m[None, :]; sc[:, 5:10] = (h / H[1:])[:, None] ** m[None, :]
        return sc

    def local(self, y, h):
        """local coefficient vectors d (N, nb); unknowns are yhat_{j,m} = H_j^m kappa^(m)(x_j) and bubble coefficients"""
        ye = np.r_[y, 0.0]
        d = np.where(self.L >= 0, ye[self.L], self.fixval)
        return d * self.lscale(h)

    def to_y(self, D, bub, h):
        """global unknown vector from physical nodal derivatives D (N+1,5) and bubbles (N, nbub)"""
        H = self.nodal_scale(h); y = np.zeros(self.n)
        nr = self.B.nrm
        for j in range(self.N + 1):
            for m in range(5):
                if self.nidx[j, m] >= 0: y[self.nidx[j, m]] = D[j, m] * H[j] ** m * nr[m]
        for e in range(self.N): y[self.bidx[e]] = bub[e] * nr[10:]
        return y

    def nodal(self, y, h):
        """physical nodal derivatives D (N+1, 5)"""
        H = self.nodal_scale(h); D = np.zeros((self.N + 1, 5))
        nr = self.B.nrm
        for (j, m), v in self.fixed.items(): D[j, m] = v
        for j in range(self.N + 1):
            for m in range(5):
                if self.nidx[j, m] >= 0: D[j, m] = y[self.nidx[j, m]] / H[j] ** m / nr[m]
        return D

    def eval_q(self, y, ell, k):
        """kappa^(k) at quadrature points (N, nq)"""
        x, _ = self.mesh.knots(ell); h = np.diff(x)
        d = self.local(y, h)
        return (d @ self.B.V[k].T) / h[:, None] ** k

    def eval_at(self, y, ell, xs, k=0):
        x, _ = self.mesh.knots(ell); h = np.diff(x)
        d = self.local(y, h)
        e = np.clip(np.searchsorted(x, xs, side='right') - 1, 0, self.N - 1)
        u = (xs - x[e]) / h[e]
        out = np.empty(len(xs))
        for ee in np.unique(e):
            m = e == ee
            out[m] = (self.B.evalpts(u[m], k) @ d[ee]) / h[ee] ** k
        return out
