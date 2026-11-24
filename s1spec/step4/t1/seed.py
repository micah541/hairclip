"""HLS snapshot -> exact Fourier evaluation of the rescaled curvature jets along the half curve from the tip (alpha = pi -> 0)."""
import numpy as np
from rep import K0

class HLSHalf:
    def __init__(self, path):
        z = np.load(path); self.N = int(z['N']); self.L = float(z['L']); self.h = z['th_hat']; self.t = float(z['t'])
        self.X0 = z['X0']; self.kmax = float(z['kmax'])
        self.kk = np.arange(self.N // 2 + 1)
        self.w = np.ones(self.N // 2 + 1); self.w[1:] = 2.0
        if self.N % 2 == 0: self.w[-1] = 1.0
        self.c = 2 * np.pi / self.L
        # physical tip curvature (with the reversed orientation sign)
        th1 = self.theta_der(np.array([np.pi]), 1)[0]
        self.ktip = -self.c * th1
        self.lam = K0 / self.ktip
        self.ell = self.L / (2 * self.lam)

    def theta_der(self, al, m):
        """d^m/dalpha^m theta(alpha), theta = alpha + periodic"""
        al = np.atleast_1d(al); out = np.zeros(len(al))
        coef = (1j * self.kk) ** m * self.h * self.w / self.N
        for i0 in range(0, len(al), 2000):
            a = al[i0:i0 + 2000]
            out[i0:i0 + 2000] = np.real(np.exp(1j * np.outer(a, self.kk)) @ coef)
        if m == 0: out += al
        if m == 1: out += 1.0
        return out

    def kap_raw(self, x, m=0):
        """rescaled kappa^(m)(x) along the strand from the tip (x<0: the other strand)"""
        lam = self.lam; al = np.pi - self.c * lam * np.asarray(x)
        return (-1) ** (m + 1) * (lam * self.c) ** (m + 1) * self.theta_der(al, m + 1)

    def kap(self, x, m=0):
        """even part about the tip (exactly symmetric at both ends)"""
        x = np.asarray(x, float)
        return 0.5 * (self.kap_raw(x, m) + (-1) ** m * self.kap_raw(-x, m))
