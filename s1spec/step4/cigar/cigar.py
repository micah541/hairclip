"""Quick numerical test (not rigorous): does a long thin convex stadium ('cigar') lose embeddedness / blow up under curve
diffusion?  Same HLS solver as the hairpin runs.  Stadium: half-width r, straight length ell, mirror symmetric about y = 0,
CCW from the leftmost point.  Reports kmax, number of extra axis crossings (self-intersections), min half-width of the strip."""
import sys, time, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4')
from numpy.fft import rfft, irfft
from hls import HLS
from run_hairpin import axis_diag
r = float(sys.argv[1]); ell = float(sys.argv[2]); N = int(sys.argv[3]); tend = float(sys.argv[4]); wall = float(sys.argv[5])
M = 1 << 20; amoll = min(0.025, r / 8)
segs = [(np.pi / 2 * r, 1 / r), (ell, 0.0), (np.pi * r, 1 / r), (ell, 0.0), (np.pi / 2 * r, 1 / r)]
L = sum(l for l, _ in segs); s = L * np.arange(M) / M; ds = L / M
edges = np.cumsum([0] + [l for l, _ in segs]); lo, hi = s - ds / 2, s + ds / 2; k = np.zeros(M)
for (l, kk), e0, e1 in zip(segs, edges[:-1], edges[1:]):
    for sh in (-L, 0.0, L):
        k += kk * np.clip(np.minimum(hi, e1 + sh) - np.maximum(lo, e0 + sh), 0, None) / ds
x = (np.arange(M) - M // 2) * ds; phi = np.zeros(M); ins = np.abs(x) < amoll
phi[ins] = np.exp(-1.0 / (1 - (x[ins] / amoll) ** 2)); phi /= phi.sum(); phi = np.roll(phi, -M // 2)
k = irfft(rfft(k) * rfft(phi), M)
th = -np.pi / 2 + np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
wv = np.sin(2 * np.pi * s / L)
for _ in range(3):
    e = np.exp(1j * th); D = np.sum(e) * ds; dD = np.sum(1j * e * wv) * ds; th = th + (-D.imag / dD.imag) * wv
sj = L * np.arange(N) / N; thj = np.interp(sj, s, th); al = 2 * np.pi * np.arange(N) / N
S = HLS(thj - al, L, (-(ell / 2 + r), 0.0), omega=1, cdt=1e-2, tail_tol=1e-10, Nmax=1 << 17)
NMIN = int(sys.argv[6]) if len(sys.argv) > 6 else 0
_co = S.coarsen
S.coarsen = (lambda: _co() if S.N > NMIN else None)
A0 = S.area(); t0 = time.time(); out = []
def minwidth(S):
    a, th_, Z, clo = S.curve(2); n = len(a); y = Z.imag[: n // 2]; xx = Z.real[: n // 2]
    xs = np.sort(xx); band = (xx > xs[int(0.02 * len(xs))]) & (xx < xs[int(0.98 * len(xs))])
    i = np.argmin(np.abs(y[band])); return np.abs(y[band])[i], xx[band][i], Z.real.min(), Z.real.max()
nxt = 0.0
while S.t < tend and time.time() - t0 < wall:
    S.run_to(maxsteps=20)
    if S.t >= nxt or S.kmax() > 50 / r:
        d = axis_diag(S); mw, xmw, xl, xr = minwidth(S)
        line = (f"t={S.t:.5f} t/r^4={S.t / r**4:9.1f} N={S.N} L={S.L:.4f} dA={S.area()/A0-1:+.1e} kmax*r={S.kmax()*r:7.3f} "
                f"extra_axis_crossings={len(d['cross'])} min_halfwidth/r={mw / r:.4f} at x={xmw:.3f} span=[{xl:.3f},{xr:.3f}]")
        print(line, flush=True); out.append(line); nxt = S.t * 1.35 + 1e-4 * r**4
        if S.kmax() * r > 200: print("STOP: large curvature"); break
        if len(d["cross"]) > 0: print("   crossings at x =", [round(c[1], 4) for c in d["cross"]])
print("wall", round(time.time() - t0), "s")
