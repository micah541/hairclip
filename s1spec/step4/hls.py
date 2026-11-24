"""Curve diffusion flow  X_t = -k_ss n  (n = i t, k = theta_s) for closed plane curves.
Hou-Lowengrub-Shelley tangent-angle / equal-arclength formulation, Fourier pseudo-spectral in alpha in [0, 2pi),
ETDRK4 (Kassam-Trefethen) in time.  theta = omega*alpha + th (th periodic), s_alpha = L/2pi.
   U = -k_ss = -(2pi/L)^3 th_aaa,     f = U theta_a,     L_t = -2pi mean(f),
   T_v = int_0^alpha (f - mean f),    theta_t = (2pi/L)(U_a + T_v theta_a).
Stiff linear part treated exactly: -(2pi/L)^4 (d^4 + omega^2 d^2).  Adaptive N (spectral tail) and dt = c/max|k|^4.
"""
import numpy as np
from numpy.fft import rfft, irfft, fft, ifft


class HLS:
    def __init__(self, th, L, X0, omega=1, cdt=1e-2, tail_tol=1e-11, Nmax=1 << 17, filt=True):
        self.N = len(th); self.th_hat = rfft(th); self.L = float(L); self.X0 = np.array(X0, float)
        self.om = omega; self.cdt = cdt; self.tail_tol = tail_tol; self.Nmax = Nmax; self.t = 0.0
        self.filt = filt
        self._setup()
        self.coef_key = None

    def _setup(self):
        N = self.N
        self.kap = np.arange(N // 2 + 1, dtype=float)
        self.ik = 1j * self.kap
        kc = self.kap / (N // 2)
        self.F = np.exp(-36.0 * kc**36) if self.filt else np.ones_like(kc)

    # ---------- spectral helpers
    def d(self, vh, n):
        return irfft((self.ik**n) * vh, self.N)

    def nonlin(self, th_hat):
        """time tau' with dtau'/dt = (2pi/L)^4:  th_tau' = -(d^4 + om^2 d^2) th + N(th),
        N = -G theta_a + om^2 th_aa,  G = int_0^alpha (g - mean g),  g = th_aaa theta_a.
        Slaved: (log L)_tau' = mean g,  t_tau' = (L/2pi)^4,  X0_tau' = -(L/2pi) th_aaa(0) n(0)."""
        N = self.N
        th_a = self.om + self.d(th_hat, 1)
        th3 = self.d(th_hat, 3)
        g = th3 * th_a
        gh = rfft(g)
        gbar = gh[0].real / N
        gh[0] = 0.0
        gh[1:] = gh[1:] / self.ik[1:]
        G = irfft(gh, N); G -= G[0]
        Nh = -rfft(G * th_a) + (self.om**2) * (-(self.kap**2)) * th_hat
        th0 = irfft(th_hat, N)[0]
        n0 = np.array([-np.sin(th0), np.cos(th0)])
        return Nh, gbar, th3[0], n0

    def coeffs(self, h):
        key = (h, self.N)
        if self.coef_key == key: return
        lin = -(self.kap**4 - self.om**2 * self.kap**2)
        M = 32
        r = np.exp(1j * np.pi * (np.arange(1, M + 1) - 0.5) / M)
        LR = h * lin[:, None] + r[None, :]
        self.E = np.exp(h * lin); self.E2 = np.exp(h * lin / 2)
        self.Q = h * np.real(np.mean((np.exp(LR / 2) - 1) / LR, axis=1))
        self.f1 = h * np.real(np.mean((-4 - LR + np.exp(LR) * (4 - 3 * LR + LR**2)) / LR**3, axis=1))
        self.f2 = h * np.real(np.mean((2 + LR + np.exp(LR) * (-2 + LR)) / LR**3, axis=1))
        self.f3 = h * np.real(np.mean((-4 - 3 * LR - LR**2 + np.exp(LR) * (4 - LR)) / LR**3, axis=1))
        self.coef_key = key

    def _slaved(self, gbar, th3, n0, logL):
        Lr = np.exp(logL) / (2 * np.pi)
        return gbar, Lr**4, -Lr * th3 * n0

    def step(self, h):
        """one ETDRK4 step of size h in tau'"""
        self.coeffs(h)
        u = self.th_hat; lL = np.log(self.L); t0 = self.t; X0 = self.X0
        Nu, g1, a1, n1 = self.nonlin(u); s1 = self._slaved(g1, a1, n1, lL)
        a = self.E2 * u + self.Q * Nu
        Na, g2, a2, n2 = self.nonlin(a); s2 = self._slaved(g2, a2, n2, lL + h / 2 * s1[0])
        b = self.E2 * u + self.Q * Na
        Nb, g3, a3, n3 = self.nonlin(b); s3 = self._slaved(g3, a3, n3, lL + h / 2 * s2[0])
        c = self.E2 * a + self.Q * (2 * Nb - Nu)
        Nc, g4, a4, n4 = self.nonlin(c); s4 = self._slaved(g4, a4, n4, lL + h * s3[0])
        self.th_hat = self.E * u + self.f1 * Nu + 2 * self.f2 * (Na + Nb) + self.f3 * Nc
        if self.filt: self.th_hat *= self.F
        comb = lambda i: (s1[i] + 2 * s2[i] + 2 * s3[i] + s4[i]) / 6
        self.L = np.exp(lL + h * comb(0)); self.t = t0 + h * comb(1); self.X0 = X0 + h * comb(2)

    # ---------- diagnostics
    def kmax(self):
        return (2 * np.pi / self.L) * np.max(np.abs(self.om + self.d(self.th_hat, 1)))

    def tail(self):
        a = np.abs(self.th_hat[1:])
        n = len(a)
        return np.max(a[int(0.8 * n):]) / max(np.max(a), 1e-300)

    def refine(self):
        N2 = 2 * self.N
        h2 = np.zeros(N2 // 2 + 1, complex)
        h2[:self.N // 2 + 1] = self.th_hat * 2.0
        h2[self.N // 2] *= 0.5
        self.th_hat = h2; self.N = N2; self._setup(); self.coef_key = None

    def coarsen(self):
        N2 = self.N // 2
        h2 = self.th_hat[:N2 // 2 + 1] / 2.0
        h2 = h2.copy(); h2[-1] = h2[-1].real
        self.th_hat = h2; self.N = N2; self._setup(); self.coef_key = None

    def tail_half(self):
        """tail measure the solution would have at N/2"""
        a = np.abs(self.th_hat[1:self.N // 4 + 1])
        n = len(a)
        return np.max(a[int(0.8 * n):]) / max(np.max(a), 1e-300)

    def theta(self, up=1):
        N = self.N * up
        h = np.zeros(N // 2 + 1, complex); h[:self.N // 2 + 1] = self.th_hat * up
        if up > 1: h[self.N // 2] *= 0.5
        al = 2 * np.pi * np.arange(N) / N
        return al, self.om * al + irfft(h, N)

    def curve(self, up=1):
        al, th = self.theta(up)
        N = len(al)
        g = np.exp(1j * th)
        gh = fft(g); mean = gh[0] / N
        kk = np.fft.fftfreq(N, 1.0 / N)
        Gh = np.zeros_like(gh); nz = kk != 0
        Gh[nz] = gh[nz] / (1j * kk[nz])
        G = ifft(Gh); G = G - G[0]
        Z = (self.X0[0] + 1j * self.X0[1]) + (self.L / (2 * np.pi)) * (G + mean * al)
        return al, th, Z, abs(mean) * self.L

    def area(self):
        al, th, Z, _ = self.curve(2)
        dZ = (self.L / (2 * np.pi)) * np.exp(1j * th)
        return 0.5 * np.sum(np.imag(np.conj(Z) * dZ)) * (2 * np.pi / len(al))

    def run_to(self, tend=None, kstop=None, every=None, cb=None, maxsteps=10**9, verbose=False):
        n = 0
        while n < maxsteps:
            km = self.kmax()
            if kstop is not None and km >= kstop: break
            tha = np.max(np.abs(self.om + self.d(self.th_hat, 1)))
            h = self.cdt / tha**4
            hprev = getattr(self, 'hprev', None)
            if hprev is None: hprev = 1e-4 * h          # ramp up from a small first step (initial layer)
            h = min(h, 2.0**0.25 * hprev)
            h = 2.0 ** (np.floor(4 * np.log2(h)) / 4)
            self.hprev = h
            if tend is not None:
                dtdtau = (self.L / (2 * np.pi))**4
                if self.t + h * dtdtau > tend:
                    h = (tend - self.t) / dtdtau
                    if h <= 1e-300: break
            if self.tail() > self.tail_tol and self.N < self.Nmax:
                self.refine()
                if verbose: print(f"  refine -> N={self.N} at t={self.t:.6e}", flush=True)
            elif n % 50 == 0 and self.N > 1024 and self.tail_half() < 1e-3 * self.tail_tol:
                self.coarsen()
                if verbose: print(f"  coarsen -> N={self.N} at t={self.t:.6e}", flush=True)
            self.step(h); n += 1
            if cb is not None and every is not None and n % every == 0: cb(self)
            if not np.all(np.isfinite(self.th_hat)): raise FloatingPointError("NaN")
        return n


def from_points(Zfun, omega, N, M=200000):
    """equal-arclength theta from a parametric closed curve Zfun(phi), phi in [0, 2pi)"""
    ph = np.linspace(0, 2 * np.pi, M + 1)
    Z = Zfun(ph)
    dZ = np.diff(Z)
    ds = np.abs(dZ)
    s = np.concatenate([[0], np.cumsum(ds)])
    L = s[-1]
    ang = np.unwrap(np.angle(dZ))
    smid = 0.5 * (s[1:] + s[:-1])
    sj = L * np.arange(N) / N
    th = np.interp(sj, smid, ang, period=None)
    # fix ends: use periodic interpolation of ang - omega*2pi*s/L
    per = ang - omega * 2 * np.pi * smid / L
    thj = np.interp(sj, np.concatenate([smid - L, smid, smid + L]), np.concatenate([per, per, per]))
    al = 2 * np.pi * np.arange(N) / N
    X0 = np.interp(0.0, s, Z.real), np.interp(0.0, s, Z.imag)
    return thj, L, X0
