"""Base-curve geometry for the linearization L_Gamma on the half line sigma >= 0 (rescaled arclength from the tip, orientation of the S1 chart).
Provides kappa, kappa', kappa'', q = Y.T, p = Y.N, tau = int_0^sigma kappa V_Gamma (V_Gamma = -kappa'' + p/4) as cubic splines (CubicSpline) in sigma."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/feas2'); sys.path.insert(0, '/home/claude/s1spec/step4/agent_handoff2'); sys.path.insert(0, '/home/claude/s1spec')
from scipy.interpolate import CubicSpline
from numpy.fft import rfft, irfft, fft, ifft
from common import *
import handoff2 as H
from evans import PROF, STIP

def s1_geometry(xgrid):
    """S1 profile at x = sigma - sigma_tip (x>=0).  Beyond sigma_tot=22 the curve is the straight ray (p=k=m=0, q=sigma_tot), as in measure.py."""
    sig = xgrid + STIP; P = np.zeros((4, len(xgrid))); ins = sig <= 22
    P[:, ins] = PROF(sig[ins]); P[1, ~ins] = sig[~ins]
    p, q, k, m = P
    return dict(x=xgrid, kap=k, kap1=m, kap2=p / 4, q=q, p=p, tau=np.zeros_like(xgrid))

def numeric_geometry(path, xc=None, lam=None, frame='fit', Smax=45.0, dsig=0.01):
    """Geometry of the rescaled numerical curve on sigma in [0,Smax].  xc, lam: frame (complex centre, scale).  Returns dict on a uniform grid (spacing <= dsig)."""
    S = Snap(path); N, L = S.N, S.L
    tha0 = 1.0 + irfft(1j * S.kap * S.h, N)[N // 2]
    kt = abs((2 * np.pi / L) * tha0)
    if lam is None: lam = K0S1 / kt
    # fine periodic grid with M points, spacing in sigma <= dsig
    M = N
    while L / M / lam > dsig: M *= 2
    up = M // N
    hh = np.zeros(M // 2 + 1, complex); hh[:N // 2 + 1] = S.h * up
    if up > 1: hh[N // 2] *= 0.5
    kk = np.arange(M // 2 + 1)
    al = 2 * np.pi * np.arange(M) / M
    thp = irfft(hh, M); th = al + thp
    d1 = 1.0 + irfft(1j * kk * hh, M); d2 = irfft((1j * kk) ** 2 * hh, M); d3 = irfft((1j * kk) ** 3 * hh, M)
    g = np.exp(1j * th); gh = fft(g); mean = gh[0] / M
    kf = np.fft.fftfreq(M, 1.0 / M); Gh = np.zeros_like(gh); nz = kf != 0; Gh[nz] = gh[nz] / (1j * kf[nz])
    G = ifft(Gh); G = G - G[0]
    Z = (S.X0[0] + 1j * S.X0[1]) + (L / (2 * np.pi)) * (G + mean * al)
    # strand from the tip, alpha decreasing from pi
    nS = int(Smax * lam / (L / M)) + 2; it = M // 2
    idx = it - np.arange(nS)
    sig = np.arange(nS) * (L / M) / lam
    Zs, ths = Z[idx], th[idx] + np.pi            # reversed orientation: theta_tilde = theta + pi
    kap = -lam * (2 * np.pi / L) * d1[idx]
    kap1 = lam**2 * (2 * np.pi / L) ** 2 * d2[idx]
    kap2 = -lam**3 * (2 * np.pi / L) ** 3 * d3[idx]
    tip = Zs[0]
    if xc is None: xc = tip.real + H.D_S1 * lam + 1j * tip.imag
    Y = (Zs - xc) / lam; T = np.exp(1j * ths)
    q = (Y * np.conj(T)).real; p = (Y * np.conj(1j * T)).real
    V = -kap2 + p / 4
    # tau = int_0^sigma kappa V  (spline antiderivative)
    tau = CubicSpline(sig, kap * V).antiderivative()(sig)
    return dict(x=sig, kap=kap, kap1=kap1, kap2=kap2, q=q, p=p, tau=tau, V=V, lam=lam, xc=xc, Y=Y, T=T, M=M, kt=kt, tip=tip)

def as_splines(g):
    return {k: CubicSpline(g['x'], g[k]) for k in ('kap', 'kap1', 'kap2', 'q', 'p', 'tau')}

def s1_splines(Smax=45.0, dx=0.01):
    x = np.arange(0, Smax + dx / 2, dx); return as_splines(s1_geometry(x))

if __name__ == '__main__':
    # consistency: numeric geometry at k=1600 vs S1
    g = numeric_geometry('/home/claude/s1spec/step4/feas/snaps/base_1600.npz')
    s = s1_geometry(g['x'][:int(30 / 0.0066)])
    n = len(s['x'])
    for key in ('kap', 'kap1', 'kap2', 'q', 'p'):
        d = g[key][:n] - s[key]
        print(key, 'max|num-S1| on [0,30]:', np.max(np.abs(d)), ' at x<=5:', np.max(np.abs(d[g['x'][:n] <= 5])), 'scale', np.max(np.abs(s[key])))
    print('kt', g['kt'], 'lam', g['lam'], 'xc', g['xc'], 'tau(30)', g['tau'][n-1], 'V max', np.max(np.abs(g['V'][:n])))
