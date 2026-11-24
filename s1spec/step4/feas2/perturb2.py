"""Closure-preserving even perturbation of an HLS state (generalisation of feas/perturb.py to bumps centred at +-s_c).
Shapes (all even in s = arclength from alpha=pi; pair of Gaussians at +-s_c, b = g(w) - c g(2w), c from  int kappa b = 0):
  peak : centre s_c = arclength of max|kappa| (0 for t>=0.3, where the tip is the curvature peak), w = 2*ell, amp = ell,  ell = 1.5398/kmax(t_j)
  fixed: centre 0 (tip, alpha=pi), w = 0.05, amp = 0.05   (slot/chamber scale)
  far  : centre at the farthest point of the curve in +x (far loop/bulb), w = 0.3, amp = 0.3
peak normal displacement = delta*amp.   dtheta = e' + kappa * int_0^s kappa e   (first order; keeps L and closure)."""
import sys, numpy as np
sys.dont_write_bytecode = True
from numpy.fft import rfft, irfft
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from common import *
import hls

def pairs(s, sc, w):
    g = lambda c: np.exp(-(s - c) ** 2 / (2 * w ** 2))
    gp = lambda c: -(s - c) / w ** 2 * g(c)
    return g(sc) + g(-sc), gp(sc) + gp(-sc)

def make(path_in, path_out, shape, delta):
    S = Snap(path_in); N, L = S.N, S.L
    A = hls.HLS(np.zeros(N), L, S.X0, omega=1); A.th_hat = S.h.copy(); _, th0, Z0, clo0 = A.curve(1)
    kmax = float(np.load(path_in)['kmax']); ell = K0S1 / kmax
    j = np.arange(N); s = L * (j - N // 2) / N
    thal = 1.0 + irfft(1j * S.kap * S.h, N); kap = (2 * np.pi / L) * thal
    if shape == 'peak':
        sc = abs(s[np.argmax(np.where(s > 0, np.abs(kap), 0))]); w, amp = 2 * ell, ell
        if sc < 2 * w: sc = 0.0 if sc < 0.5 * w else sc
    elif shape == 'fixed': sc, w, amp = 0.0, 0.05, 0.05
    elif shape == 'far':
        sc = float(s[np.argmax(np.where(s > 0, Z0.real, -1e9))]); w, amp = 0.3, 0.3
    g1, g1p = pairs(s, sc, w); g2, g2p = pairs(s, sc, 2 * w)
    c = np.sum(kap * g1) / np.sum(kap * g2)
    b = g1 - c * g2; bp = g1p - c * g2p
    m = np.max(np.abs(b)); b, bp = b / m, bp / m
    e, ep = delta * amp * b, delta * amp * bp
    f = kap * e; mean = f.mean(); f = f - mean
    ks = 2 * np.pi * S.kap / L; fh = rfft(f); Fh = np.zeros_like(fh); Fh[1:] = fh[1:] / (1j * ks[1:])
    F = irfft(Fh, N); F = F - F[N // 2]
    dth = ep + kap * F
    h2 = S.h + rfft(dth)
    z = dict(np.load(path_in)); z['th_hat'] = h2
    np.savez(path_out, **z)
    B = hls.HLS(np.zeros(N), L, S.X0, omega=1); B.th_hat = h2; _, th1, Z1, clo1 = B.curve(1)
    Nn = 1j * np.exp(1j * th0); ech = ((Z1 - Z0) * np.conj(Nn)).real
    sup = (np.abs(np.abs(s) - sc) < 5 * w)
    err = np.max(np.abs(ech[sup] - e[sup])) / max(np.max(np.abs(e)), 1e-300)   # includes the O(delta) tangential/second-order bits
    return dict(kmax=kmax, ell=ell, shape=shape, delta=delta, sc=float(sc), w=float(w), amp=float(amp), c=float(c), mean_kappa_e=float(mean),
                clo_before=float(clo0), clo_after=float(clo1), normal_disp_relerr=float(err), peak_e=float(np.max(np.abs(e))),
                dth_sym=float(np.max(np.abs(dth - dth[(-j) % N]))), dth_max=float(np.max(np.abs(dth))), N=N)
