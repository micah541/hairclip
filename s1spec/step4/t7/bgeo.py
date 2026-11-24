"""T7 float geometry of the T1 base (t1/slabs_P11.pkl) at (slab i, relative time theta), at arbitrary points x in [0, ell].
Per T1 element (t = 2u - 1): kappa = exact Chebyshev series (Rep basis), derivatives by chebder, Q = int kappa'^2 and
theta = -pi/2 + int kappa by chebint (exact), Gamma = -D_S1 + int e^{i theta} (48-point Chebyshev fit, ~1e-15).
Linearisation (T5 1.1, signs of T1 signs.py; b1 simplified exactly with T1's identity for varpi):
   L e = -e'''' + b2 e'' + b1 e' + b0 e,  b2 = -k^2,  b1 = -6 k k' + Q - (1/4 - a) x,  b0 = 1/4 - a - 3 k'^2 - 4 k k''.
Derivative level (v = e'):  (L e)' = L1 v + b0' e,  L1 = -d^4 + c2 d^2 + c1 d + c0,  c2 = -k^2, c1 = -8 k k' + Q - (1/4 - a) x,
   c0 = -8 k'^2 - 10 k k''.   Directions g = (p, n, Vs), Vs = -k'' + (1/4 - a) p + b n (residual part of V~ dropped)."""
import sys, numpy as np
from numpy.polynomial import chebyshev as CH
sys.path.insert(0, '/home/claude/s1spec/step4/t1')
from base2 import Base
from rep import K0, DS1, STIP

_B = None
def base():
    global _B
    if _B is None: _B = Base('P11')
    return _B

class BaseGeom:
    def __init__(self, i, theta, nfit=48):
        B = base(); R, D = B.rep(i); self.i, self.theta = i, theta
        self.tau = B.tau0[i] + theta * B.dt[i]
        Y = B.Y(self.tau, i); n = R.n; self.ell = ell = Y[n]
        x, _ = R.mesh.knots(ell); h = np.diff(x); d = R.local(Y[:n], h)
        p = D.parts(Y); self.a = p['a']; self.c = 0.25 - self.a; self.b = p['k2t'] + self.c * DS1
        self.k = K0 * np.exp(-Y[n + 1]); self.xk, self.h, self.N = x, h, R.N
        self.ck = [R.B.cheb @ d[e] for e in range(R.N)]                 # kappa(t), t in [-1,1]
        tg = np.cos(np.pi * (np.arange(nfit) + 0.5) / nfit)            # Chebyshev points
        self.cQ, self.cth, self.cG = [], [], []
        Q0, th0, G0 = 0.0, -np.pi / 2, -DS1 + 0j
        for e in range(R.N):
            s = h[e] / 2                                                 # dx = s dt
            k1 = CH.chebder(self.ck[e]) / s
            cq = CH.chebint(CH.chebmul(k1, k1), lbnd=-1) * s; cq[0] += Q0; self.cQ.append(cq)
            ct = CH.chebint(self.ck[e], lbnd=-1) * s; ct[0] += th0; self.cth.append(ct)
            E = np.exp(1j * CH.chebval(tg, ct))
            cE = CH.chebfit(tg, E.real, nfit - 1) + 1j * CH.chebfit(tg, E.imag, nfit - 1)
            cg = CH.chebint(cE, lbnd=-1) * s; cg[0] += G0; self.cG.append(cg)
            Q0, th0, G0 = CH.chebval(1.0, cq), CH.chebval(1.0, ct), CH.chebval(1.0, cg)
        self.Qend, self.closure = Q0, abs(G0.imag)
    def elem(self, xs):
        return np.clip(np.searchsorted(self.xk, xs, side='right') - 1, 0, self.N - 1)
    def at(self, xs, mmax=5):
        """dict of arrays at points xs: K[m] (m <= mmax), Q, th, G, p, q, n, t1, Vs"""
        xs = np.atleast_1d(np.asarray(xs, float)); e = self.elem(xs)
        out = dict(x=xs, K=np.zeros((mmax + 1, len(xs))), Q=np.zeros(len(xs)), th=np.zeros(len(xs)), G=np.zeros(len(xs), complex))
        for ee in np.unique(e):
            m = e == ee; t = 2 * (xs[m] - self.xk[ee]) / self.h[ee] - 1; s = self.h[ee] / 2
            c = self.ck[ee]
            for j in range(mmax + 1):
                out['K'][j, m] = CH.chebval(t, c) / s**j; c = CH.chebder(c)
            out['Q'][m] = CH.chebval(t, self.cQ[ee]); out['th'][m] = CH.chebval(t, self.cth[ee]); out['G'][m] = CH.chebval(t, self.cG[ee])
        T = np.exp(1j * out['th']); nu = 1j * T; G = out['G']
        out['p'] = (G * np.conj(nu)).real; out['q'] = (G * np.conj(T)).real; out['n'] = nu.real; out['t1'] = T.real
        out['Vs'] = -out['K'][2] + self.c * out['p'] + self.b * out['n']
        return out
    def coeffs(self, g):
        """b2, b1, b0 and derivative-level c2, c1, c0 (+ needed derivatives) from a sample dict g"""
        K = g['K']; x = g['x']; c = self.c
        k, k1, k2, k3 = K[0], K[1], K[2], K[3]
        b2 = -k**2; b1 = -6 * k * k1 + g['Q'] - c * x; b0 = 0.25 - self.a - 3 * k1**2 - 4 * k * k2
        return dict(b2=b2, b1=b1, b0=b0)

def smoothstep(u):
    u = np.clip(u, 0, 1); return u**4 * (35 - 84 * u + 70 * u**2 - 20 * u**3)
def smoothstep_der(u, n):
    """derivatives of S(u) w.r.t. u (0 outside (0,1) for n >= 1)"""
    P = np.array([0, 0, 0, 0, 35, -84, 70, -20], float)
    for _ in range(n): P = np.polynomial.polynomial.polyder(P)
    uc = np.clip(u, 0, 1); v = np.polynomial.polynomial.polyval(uc, P)
    if n >= 1: v = np.where((u > 0) & (u < 1), v, 0.0)
    return v

def chi32(x, n=0):
    """chi_32(sigma) = S((32 - sigma)/3), sigma = x + STIP; derivatives in x (0..4)"""
    u = (32 - (x + STIP)) / 3
    return smoothstep_der(u, n) * (-1 / 3) ** n if n else smoothstep(u)
