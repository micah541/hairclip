"""Smooth reflection-symmetric, closure-preserving perturbation of a HLS snapshot.
Normal displacement e(s) (even in the arclength s from the tip) with  int kappa e ds = 0  (so L is unchanged):
   delta theta(s) = e'(s) + kappa(s) * int_0^s kappa e        (first order in delta; derived in FEAS_REPORT)
added to th on the equal-arclength grid; th_hat += rfft(dtheta).  Shapes (b = g(s/w) - c g(s/2w), Gaussians, c fixed by int kappa b = 0):
   tip:     w = 2 lam_j,   e = delta * lam_j * b            (lam_j = 1.5398/kmax: 'current scale', relative amplitude delta)
   chamber: w = 0.05,      e = delta * 0.05 * b             (slot/chamber scale, fixed physical width)
usage: perturb.py in.npz out.npz shape delta"""
import sys, numpy as np
from numpy.fft import rfft, irfft
from common import *
RHO = lambda x: np.exp(-0.25 * (1 + x**2) ** (2 / 3))

def profile(s, w, kap):
    g1 = np.exp(-s**2 / (2 * w**2)); g2 = np.exp(-s**2 / (8 * w**2))
    ds = s[1] - s[0]
    c = np.sum(kap * g1) / np.sum(kap * g2)
    b = g1 - c * g2
    bp = -s / w**2 * g1 + c * s / (4 * w**2) * g2
    m = np.max(np.abs(b)); b, bp = b / m, bp / m          # peak normal displacement = delta * amp
    return b, bp, c

def make(path_in, path_out, shape, delta):
    S = Snap(path_in); N, L = S.N, S.L
    km = S.kmax(); lam = K0S1 / km
    j = np.arange(N); s = L * (j - N // 2) / N
    thal = 1.0 + irfft(1j * S.kap * S.h, N)                    # theta_alpha on the grid
    kap = (2 * np.pi / L) * thal
    w, amp = (2 * lam, lam) if shape == 'tip' else (0.05, 0.05)
    b, bp, c = profile(s, w, kap)
    e, ep = delta * amp * b, delta * amp * bp
    f = kap * e; mean = f.mean()
    f = f - mean                                                 # int kappa e = 0 up to roundoff (reported)
    ks = 2 * np.pi * S.kap / L; fh = rfft(f); Fh = np.zeros_like(fh); Fh[1:] = fh[1:] / (1j * ks[1:])
    F = irfft(Fh, N); F = F - F[N // 2]                          # int_0^s kappa e  (periodic, zero at the tip)
    dth = ep + kap * F
    h2 = S.h + rfft(dth)
    z = dict(np.load(path_in)); z['th_hat'] = h2; z['kmax'] = km
    np.savez(path_out, **z)
    # diagnostics: closure and normal displacement check
    import hls
    A = hls.HLS(np.zeros(N), L, S.X0, omega=1); A.th_hat = S.h.copy()
    B = hls.HLS(np.zeros(N), L, S.X0, omega=1); B.th_hat = h2
    _, th0, Z0, clo0 = A.curve(1); _, th1, Z1, clo1 = B.curve(1)
    Nn = 1j * np.exp(1j * th0); ech = ((Z1 - Z0) * np.conj(Nn)).real
    # displaced curve point at the same arclength: Z1(S) = Z0(S) + e N + (tangential); compare normal comp. in the support
    sup = np.abs(s) < 6 * w
    err = np.max(np.abs(ech[sup] - e[sup])) / max(np.max(np.abs(e)), 1e-300)
    xr = np.abs(s) / lam
    m = xr <= 28; hx = np.diff(xr[N // 2:][:1])  # unused
    dx = (s[1] - s[0]) / lam
    size_rho = float(np.sqrt(np.sum((RHO(xr) * (e / lam) ** 2)[(xr <= 28) & (s >= 0)]) * dx))
    size_W = float(np.sqrt(np.sum(((1 + xr**2) ** -2 * (e / lam) ** 2)[(xr <= 28) & (s >= 0)]) * dx))
    info = dict(kmax=km, lam=lam, shape=shape, delta=delta, w=w, c=c, mean_kappa_e=float(mean), clo_before=float(clo0), clo_after=float(clo1),
                normal_disp_relerr=float(err), max_e_over_lam=float(np.max(np.abs(e)) / lam), size_rho=size_rho, size_W=size_W,
                dtheta_far=float(np.max(np.abs(dth[np.abs(s) > 8 * 2 * w]))) if np.any(np.abs(s) > 16 * w) else 0.0,
                dtheta_max=float(np.max(np.abs(dth))), dth_sym=float(np.max(np.abs(dth + dth[(-j) % N]))) )
    return info

if __name__ == '__main__':
    print(make(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])))
