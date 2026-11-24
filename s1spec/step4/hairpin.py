"""The hairpin gamma0 of the notes: boundary of {dist(p,H) <= 1/2}, H = {(x,+-0.55): 0<=x<=12} + left half-circle
of radius 0.55.  CCW traversal from the leftmost point (-1.05, 0), theta(0) = -pi/2.  Curvature jumps mollified over
arclength 0.05 (C-infinity bump of half-width 0.025)."""
import numpy as np
from numpy.fft import rfft, irfft

def curvature_fn(M=1 << 20, arm=12.0, a_moll=0.025, W=1.0, w=0.1):
    R_out, r_cap, r_slot = w / 2 + W, W / 2, w / 2
    import os
    rc = float(os.environ.get("TIPRC", "0"))          # 0: round caps; rc > 0: squared tips with corner radius rc
    cap = [(np.pi * r_cap, 1 / r_cap)] if rc <= 0 else [(np.pi / 2 * rc, 1 / rc), (W - 2 * rc, 0.0), (np.pi / 2 * rc, 1 / rc)]
    Rch = float(os.environ.get("CHAMBER", "0"))      # >0: keyhole chamber of radius Rch at the slot end
    if Rch > 0:
        rcn = 0.02; al = np.arccos((w + float(os.environ.get("CHDW", "0")) + 2 * rcn) / (2 * Rch + 2 * rcn))
        send = [(al * rcn, 1 / rcn), ((np.pi + 2 * al) * Rch, -1 / Rch), (al * rcn, 1 / rcn)]
    else:
        send = [(np.pi * r_slot, -1 / r_slot)]
    segs = [(np.pi / 2 * R_out, 1 / R_out), (arm, 0.0)] + cap + [(arm, 0.0)] + send + [(arm, 0.0)] + cap + [(arm, 0.0), (np.pi / 2 * R_out, 1 / R_out)]
    L = sum(l for l, _ in segs)
    s = L * np.arange(M) / M
    ds = L / M
    # exact cell averages on cells [s_j - ds/2, s_j + ds/2] (symmetric, integral-exact); the profile is periodic
    edges = np.cumsum([0] + [l for l, _ in segs])
    lo = s - ds / 2; hi = s + ds / 2
    k = np.zeros(M)
    for (l, kk), e0, e1 in zip(segs, edges[:-1], edges[1:]):
        for sh in (-L, 0.0, L):
            ov = np.clip(np.minimum(hi, e1 + sh) - np.maximum(lo, e0 + sh), 0, None)
            k += kk * ov / ds
    # periodic mollification by FFT convolution
    ds = L / M
    x = (np.arange(M) - M // 2) * ds
    phi = np.zeros(M); ins = np.abs(x) < a_moll
    phi[ins] = np.exp(-1.0 / (1 - (x[ins] / a_moll)**2))
    phi /= phi.sum()
    phi = np.roll(phi, -M // 2)
    km = irfft(rfft(k) * rfft(phi), M)
    return s, km, L

def hairpin_theta(N, arm=12.0, M=1 << 20, a_moll=0.025, W=1.0, w=0.1):
    s, k, L = curvature_fn(M, arm, a_moll, W, w)
    ds = L / M
    th = -np.pi / 2 + np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
    Z = -(w / 2 + W) + np.concatenate([[0.0], np.cumsum(0.5 * (np.exp(1j * th[1:]) + np.exp(1j * th[:-1])) * ds)])
    tot = th[-1] + 0.5 * (k[-1] + k[0]) * ds - th[0]
    clo = Z[-1] + 0.5 * (np.exp(1j * th[-1]) + np.exp(1j * th[0])) * ds - Z[0]
    # close the curve: theta -> theta + eps sin(2 pi s / L) (keeps the mirror symmetry and the turning)
    wv = np.sin(2 * np.pi * s / L)
    for _ in range(3):
        e = np.exp(1j * th)
        D = np.sum(e) * ds
        dD = np.sum(1j * e * wv) * ds
        eps = -D.imag / dD.imag
        th = th + eps * wv
    e = np.exp(1j * th); clo = np.sum(e) * ds
    sj = L * np.arange(N) / N
    thj = np.interp(sj, s, th)
    al = 2 * np.pi * np.arange(N) / N
    return thj - al, L, (-(w / 2 + W), 0.0), tot, clo

if __name__ == "__main__":
    th, L, X0, tot, clo = hairpin_theta(16384)
    print(f"L = {L:.6f}, total turning/2pi = {tot/2/np.pi:.12f}, closure defect = {clo:.3e}")
    from hls import HLS
    S = HLS(th, L, X0, omega=1)
    print(f"area = {S.area():.6f}, kmax = {S.kmax():.4f}, tail = {S.tail():.2e}, closure(spectral) = {S.curve()[3]:.2e}")
