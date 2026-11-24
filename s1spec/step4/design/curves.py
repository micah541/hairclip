"""Embedded mirror-symmetric closed curves (symmetric about the x-axis), built from a half-curve of (length, curvature) segments.
disk_notch(Rd, w, Rc, phi, Lw, ...): disk of radius Rd with a (flared) slot cut in from the right ending in a keyhole chamber of radius Rc.
Half curve: start at the leftmost point (-Rd,0), heading -pi/2 -> disk arc to the slot mouth -> mouth fillet -> flared wall (angle phi) ->
fillet -> straight slot (width ~w, length l0) -> fillet rcn -> half chamber (curvature -1/Rc).  The mouth angle beta is shot so the half curve ends on the axis."""
import numpy as np
from numpy.fft import rfft, irfft

def half_segments(beta, Rd, phi, Lw, w, Rc, l0, rm=0.05, rf=0.2, rcn=0.02):
    al = np.arccos((w + 2 * rcn) / (2 * Rc + 2 * rcn))
    t1 = np.pi / 2 - phi + beta
    segs = [(Rd * (np.pi - beta), 1 / Rd), (rm * t1, 1 / rm), (Lw, 0.0)]
    if phi > 1e-9: segs.append((rf * phi, 1 / rf))
    segs += [(l0, 0.0), (al * rcn, 1 / rcn), ((np.pi / 2 + al) * Rc, -1 / Rc)]
    return segs

def integrate(segs, M):
    L = 2 * sum(l for l, _ in segs)
    s = L * np.arange(M) / M; ds = L / M
    full = segs + segs[::-1]
    edges = np.cumsum([0] + [l for l, _ in full]); lo = s - ds / 2; hi = s + ds / 2
    k = np.zeros(M)
    for (l, kk), e0, e1 in zip(full, edges[:-1], edges[1:]):
        ov = np.clip(np.minimum(hi, e1) - np.maximum(lo, e0), 0, None); k += kk * ov / ds
    return s, k, L

def mid_y(segs):
    th = -np.pi / 2; z = 0j
    for l, k in segs:
        if abs(k) < 1e-14: z += l * np.exp(1j * th)
        else: z += (np.exp(1j * (th + k * l)) - np.exp(1j * th)) / (1j * k); th += k * l
    return z.imag, th

def solve_beta(Rd, phi, Lw, w, Rc, l0, **kw):
    f = lambda b: mid_y(half_segments(b, Rd, phi, Lw, w, Rc, l0, **kw))[0]
    b0, b1 = 0.02, 0.03; f0, f1 = f(b0), f(b1)
    for _ in range(40):
        b2 = b1 - f1 * (b1 - b0) / (f1 - f0); b0, f0, b1 = b1, f1, b2; f1 = f(b1)
        if abs(f1) < 1e-13: break
    return b1, f1

def disk_theta(N, Rd, w, Rc, phi=0.0, Lw=None, l0=0.3, a_moll=0.025, M=1 << 20, **kw):
    if Lw is None: Lw = 0.8 * Rd
    beta, err = solve_beta(Rd, phi, Lw, w, Rc, l0, **kw)
    segs = half_segments(beta, Rd, phi, Lw, w, Rc, l0, **kw)
    s, k, L = integrate(segs, M); ds = L / M
    x = (np.arange(M) - M // 2) * ds; ph = np.zeros(M); ins = np.abs(x) < a_moll
    ph[ins] = np.exp(-1.0 / (1 - (x[ins] / a_moll)**2)); ph /= ph.sum(); ph = np.roll(ph, -M // 2)
    k = irfft(rfft(k) * rfft(ph), M)
    th = -np.pi / 2 + np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
    wv = np.sin(2 * np.pi * s / L)
    for _ in range(3):
        e = np.exp(1j * th); D = np.sum(e) * ds; dD = np.sum(1j * e * wv) * ds; th = th - D.imag / dD.imag * wv
    sj = L * np.arange(N) / N; thj = np.interp(sj, s, th); al = 2 * np.pi * np.arange(N) / N
    return thj - al, L, (-Rd, 0.0), beta, err

if __name__ == "__main__":
    import sys; sys.path.insert(0, '/home/claude/s1spec/step4')
    from hls import HLS
    for args in [(2.5, 0.05, 0.3, 0.0), (3.0, 0.05, 0.3, np.radians(20))]:
        th, L, X0, beta, err = disk_theta(16384, *args)
        S = HLS(th, L, X0, omega=1); al, t, Z, clo = S.curve(4)
        # embeddedness: min distance between non-neighbouring points (coarse)
        print(args, 'L', L, 'beta', beta, 'ymid err', err, 'area', S.area(), 'kmax', S.kmax(), 'closure', clo)
