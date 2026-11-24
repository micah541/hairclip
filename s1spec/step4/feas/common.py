"""Shared helpers for the feasibility study (float only)."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec')
K0S1 = 1.5398   # nominal tip curvature of S1 (k1[0]); lam = K0S1 / kmax (task convention for F1)

class Snap:
    """Fourier data of a HLS snapshot.  theta(alpha) = alpha + th(alpha),  th = (1/N)[h0 + 2 Re sum_k h_k e^{ik alpha} (+ Nyquist)],  s = L alpha/2pi."""
    def __init__(self, path):
        z = np.load(path)
        self.N = int(z['N']); self.L = float(z['L']); self.X0 = np.array(z['X0'], float); self.t = float(z['t'])
        self.h = np.array(z['th_hat']); self.kmax_stored = float(z['kmax'])
        self.kap = np.arange(len(self.h), dtype=float)
        w = np.full(len(self.h), 2.0 / self.N); w[0] = 1.0 / self.N; w[-1] = 1.0 / self.N   # Nyquist term: real, weight 1/N
        self.c = self.h * w            # th(alpha) = Re sum c_k e^{ik alpha}  (Nyquist: c_N/2 cos)
    def th_derivs(self, a, nd=0, chunk=48):
        return self.th_derivs_d(np.asarray(a, float) - np.pi, nd, chunk)
    def th_derivs_d(self, dd, nd=0, chunk=48):
        """d^nd/dalpha^nd of the periodic part th at alpha = pi + dd (dd: offsets from the tip, given directly to avoid rounding of alpha near pi)."""
        a = np.asarray(dd, float) + np.pi; out = np.empty(len(a)); dd = np.asarray(dd, float)
        ph = np.exp(1j * np.pi * self.kap)           # e^{ik pi}
        cc = self.c * ph * (1j * self.kap) ** nd      # e^{ik a} = e^{ik pi} e^{ik (a-pi)}
        for i in range(0, len(a), chunk):
            d = dd[i:i + chunk]
            E = np.exp(1j * np.outer(d, self.kap))
            out[i:i + chunk] = (E @ cc).real
        return out
    def theta_alpha(self, a):          # d theta / d alpha
        return 1.0 + self.th_derivs(a, 1)
    def kmax(self):
        a = np.linspace(np.pi - 0.02, np.pi + 0.02, 4001)
        return (2 * np.pi / self.L) * np.max(np.abs(self.theta_alpha(a)))
    def ktip(self):
        return (2 * np.pi / self.L) * abs(float(self.theta_alpha(np.array([np.pi]))[0]))
