"""Stage 1: symmetrized (even) fit of (xc, lam, alpha) with C4 psi3/h2 and extraction of the even profiles.
Outputs (npz): inner grid x (0..64): ue, ue' ; ray profile f(y), u(y) on the base grid; geometry arrays for jets."""
import sys, time, numpy as np
from geo import *
from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline
from numpy.fft import irfft
class Arr2(CH.Arr):
    pass
def prep2(path, up=8):
    S, z = CH.H.load(path)
    al, th, Z, _ = S.curve(up); N = len(al); L = S.L
    Nb = S.N
    hh = np.zeros(N // 2 + 1, complex); hh[:Nb // 2 + 1] = S.th_hat * up; hh[Nb // 2] *= 0.5
    thd = 1 + irfft(1j * np.arange(N // 2 + 1) * hh, N)
    j = (-np.arange(N)) % N; y0 = Z[N // 2].imag
    Zm = np.conj(Z[j] - 1j * y0) + 1j * y0; thm = np.pi - th[j]; thdm = thd[j]
    A1 = CH.Arr(Z, th, L, N); A1.th = th; A1.thd = thd
    A2 = CH.Arr(Zm, thm, L, N); A2.th = thm; A2.thd = thdm
    return dict(S=S, z=z, up=A1, lo=A2, tip=Z[N // 2], L=L, N=N)
XFIT = 28.0
XCR = 24.25          # crawl (normal graph over S1) for x < XCR; ray frame (S1 is a straight ray, p < 1e-13) beyond

def graph2(A, xc, lam, nmax, slope=False):
    """chart.graph with optional slope u' = (1-k1 u) tan(phi); x-grid xg."""
    da = A.da; it = A.it
    win = int(150 * lam / (A.L / A.N)); sl = slice(max(it - win, 0), it + 1)
    Zw, Tw = A.Z[sl][::-1], -A.T[sl][::-1]; Ms = len(Zw)
    thw = A.th[sl][::-1] + np.pi; tdw = -A.thd[sl][::-1]
    th1 = np.angle(tauG)
    Pp = xc + lam * G1G[:nmax]
    def evth(jj, t):
        h = da; t2 = t * t; t3 = t2 * t
        return (2*t3 - 3*t2 + 1)*thw[jj] + (t3 - 2*t2 + t)*h*tdw[jj] + (-2*t3 + 3*t2)*thw[jj+1] + (t3 - t2)*h*tdw[jj+1]
    def ev(jj, t):
        h = da; t2 = t * t; t3 = t2 * t
        return (2*t3 - 3*t2 + 1)*Zw[jj] + (t3 - 2*t2 + t)*h*Tw[jj] + (-2*t3 + 3*t2)*Zw[jj+1] + (t3 - t2)*h*Tw[jj+1]
    def evd(jj, t):
        h = da; t2 = t * t
        return ((6*t2 - 6*t)*Zw[jj] + (3*t2 - 4*t + 1)*h*Tw[jj] + (-6*t2 + 6*t)*Zw[jj+1] + (3*t2 - 2*t)*h*Tw[jj+1]) / h
    u = np.full(nmax, np.nan); up = np.full(nmax, np.nan)
    u[0] = ((Zw[0] - Pp[0]) * np.conj(nuG[0])).real / lam; up[0] = 0.0
    j = 0
    for i in range(1, nmax):
        tt, nv = tauG[i], nuG[i]
        f = lambda jj: ((Zw[jj] - Pp[i]) * np.conj(tt)).real
        jj = max(j - 5, 0); fj = f(jj); found = False
        while jj < Ms - 2:
            fn = f(jj + 1)
            if (fj <= 0 < fn) or (fj >= 0 > fn):
                t = fj / (fj - fn)
                for _ in range(6):
                    g = ((ev(jj, t) - Pp[i]) * np.conj(tt)).real; gp = (evd(jj, t) * np.conj(tt)).real * da
                    if gp == 0: break
                    dt = -g / gp; t += dt
                    if abs(dt) < 1e-15: break
                if -0.01 <= t <= 1.01:
                    s = ((ev(jj, t) - Pp[i]) * np.conj(nv)).real
                    if abs(s) < 0.5 * lam * (1 + xg[i]):
                        u[i] = s / lam; j = jj; found = True
                        if slope:
                            ph = evth(jj, t) - th1[i]; up[i] = (1 - kG[i] * u[i]) * np.tan(np.angle(np.exp(1j * ph)))
                break
            jj += 1; fj = fn
            if jj - j > 4000: break
        if not found: break
    return (u, up) if slope else u

def psi_fit(P, a0=None, p0=(0.0, 0.0), xe=XFIT, verbose=False):
    """constrained even-part fit u_even = alpha psi3 + alpha^2 h2 (+ w), minimizing rho-weighted half-line L2 over x<=xe; params (dxc/lam0, dlog lam, alpha)."""
    S = P['S']; tip = P['tip']; N = P['N']
    al, th = S.theta(1); kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[len(al) // 2])
    lam0 = H.K0 / kt; xc0 = tip.real + H.D_S1 * lam0 + 1j * tip.imag
    n = int(round(xe / HX)) + 1
    wh = H.wh[:n]; sw = np.sqrt(wh)
    def res(p, ret=False):
        xc = xc0 + p[0] * lam0; lam = lam0 * np.exp(p[1])
        a = graph2(P['up'], xc, lam, n); b = graph2(P['lo'], xc, lam, n)
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))): return np.full(n, 1e3)
        ue = 0.5 * (a + b)
        r = (ue - p[2] * PSI0[:n] - p[2]**2 * H20[:n]) * sw
        return (r, ue, a, b) if ret else r
    if a0 is None: a0 = -3.3 / kt
    sol = least_squares(res, [p0[0], p0[1], a0], x_scale=[1e-3, 1e-3, abs(a0)], diff_step=1e-6, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    r, ue, a, b = res(sol.x, True)
    return dict(xc=xc0 + sol.x[0] * lam0, lam=lam0 * np.exp(sol.x[1]), alpha=sol.x[2], p=sol.x.copy(), r2=float(np.linalg.norm(r)), lam0=lam0, kt=kt, nfev=sol.nfev,
                oddn=float(np.sqrt(np.sum(wh * (0.5 * (a - b))**2))))

def strand_arrays(A, xc, lam, ymax_arc):
    """strand from the tip outward on the BASE-resolution subsample; returns rescaled Gamma, unit outward tangent."""
    it = A.it
    win = int(ymax_arc * lam / (A.L / A.N)); sl = slice(max(it - win, 0), it + 1)
    Zw, Tw = A.Z[sl][::-1], -A.T[sl][::-1]
    return (Zw - xc) / lam, Tw / np.abs(Tw)

def ray_profile(A, xc, lam, ymax_arc=3000.0, ds=0.1):
    """f(q) = tan(phi), u(q) along the strand in the ray frame (tau_inf). Returns q, u, f (monotone-q portion, q>=12), chart-end info."""
    G, Tn = strand_arrays(A, xc, lam, ymax_arc)
    stride = max(1, int(round(ds / (A.L / A.N / lam))))
    G = G[::stride]; Tn = Tn[::stride]
    c = G * np.exp(-1j * THINF); q = c.real; u = c.imag
    d = Tn * np.exp(-1j * THINF); f = d.imag / d.real
    ok = q > 12.0
    i0 = np.argmax(ok)
    dq = np.diff(q[i0:])
    bad = np.where(dq <= 0)[0]
    iend = i0 + (bad[0] if len(bad) else len(dq))     # q increasing on [i0, iend]
    return q[i0:iend + 1], u[i0:iend + 1], f[i0:iend + 1]

def analyse(path, up=8, a0=None, p0=(0.0, 0.0), xe=XFIT, ymax_arc=3000.0, save=None):
    t0 = time.time()
    P = prep2(path, up=up)
    fit = psi_fit(P, a0=a0, p0=p0, xe=xe)
    xc, lam, alpha = fit['xc'], fit['lam'], fit['alpha']
    ncr = int(round((XCR + 2.0) / HX)) + 1
    out = {}
    uu = []; upp = []
    for k in ('up', 'lo'):
        u, ups = graph2(P[k], xc, lam, ncr, slope=True); uu.append(u); upp.append(ups)
    out['x_crawl'] = xg[:ncr]; out['u_crawl'] = 0.5 * (uu[0] + uu[1]); out['up_crawl'] = 0.5 * (upp[0] + upp[1])
    out['u_crawl_odd'] = 0.5 * (uu[0] - uu[1]); out['up_crawl_odd'] = 0.5 * (upp[0] - upp[1])
    rp = [ray_profile(P[k], xc, lam, ymax_arc) for k in ('up', 'lo')]
    # common y grid: spline each strand, average
    qlo = max(rp[0][0][0], rp[1][0][0]); qhi = min(rp[0][0][-1], rp[1][0][-1])
    yy = np.exp(np.arange(np.log(qlo + 1e-9), np.log(qhi), 0.002))
    fs = [CubicSpline(r[0], r[2])(yy) for r in rp]; us = [CubicSpline(r[0], r[1])(yy) for r in rp]
    out.update(y=yy, f_even=0.5 * (fs[0] + fs[1]), f_odd=0.5 * (fs[0] - fs[1]), uy_even=0.5 * (us[0] + us[1]), uy_odd=0.5 * (us[0] - us[1]))
    # chart extent of each strand (raw, before common range): last q with |f|<=0.5
    ext = []
    for r in rp:
        big = np.where(np.abs(r[2]) > 0.5)[0]
        ext.append((r[0][big[0]] if len(big) else np.nan, r[0][-1]))
    out['chart_ext'] = np.array(ext)       # [strand][first q with |f|>0.5, last q of monotone-q range]
    out.update(fit)
    out['W'] = 0
    out['kmax'] = float(P['z']['kmax']); out['L'] = P['L']; out['N'] = P['N']; out['wall'] = time.time() - t0
    out['tipy'] = P['tip']
    if save: np.savez(save, **{k: np.asarray(v) for k, v in out.items()})
    return out, P

if __name__ == '__main__':
    path = sys.argv[1]; up = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    o, P = analyse(path, up=up, save=sys.argv[3] if len(sys.argv) > 3 else None)
    print({k: o[k] for k in ('kmax', 'lam', 'alpha', 'r2', 'oddn', 'nfev', 'chart_ext', 'wall')}); print("a*k", o['alpha'] * o['kt'])
