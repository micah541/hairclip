"""Even-part chart analysis of a HLS snapshot (reuses agent_handoff2/handoff2.py: S1 chart, Hermite normal graph, psi3, h2, rho weights).
u_even = (u_upper + u_lower)/2 with both strands in the SAME frame (xc, lam), axis through the tip (y = tip.imag); (xc_real, lam) and (alpha, beta2)
fitted jointly by variable projection on u_even (removes the odd part left by the ~4e-9 physical asymmetry of the snapshots).
Norms: HALF line, trapezoid, weight rho (hn) or W = (1+x^2)^-2 (hnW)."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/agent_handoff2')
import handoff2 as H
from handoff2 import xs, HX, tau1, nu1, G1, K0, D_S1, PSI3, H2, wh, XT, load, strand
from scipy.optimize import least_squares
WW = (1 + xs[:len(PSI3)] ** 2) ** -2          # W weight
whW = np.concatenate([[HX / 2], np.full(len(PSI3) - 1, HX)]) * WW
def hn(f, n): return float(np.sqrt(np.sum(wh[:n] * f[:n] ** 2)))
def hnW(f, n): return float(np.sqrt(np.sum(whW[:n] * f[:n] ** 2)))

class Arr:
    def __init__(self, Z, th, L, N):
        it = N // 2; self.da = 2 * np.pi / N
        T = (L / (2 * np.pi)) * np.exp(1j * th)
        self.Z = Z; self.T = T; self.it = it; self.N = N; self.L = L

def prep(path, up=8):
    S, z = load(path)
    al, th, Z, _ = S.curve(up); N = len(al); L = S.L
    j = (-np.arange(N)) % N; y0 = Z[N // 2].imag
    Zm = np.conj(Z[j] - 1j * y0) + 1j * y0; thm = np.pi - th[j]
    return dict(S=S, z=z, up=Arr(Z, th, L, N), lo=Arr(Zm, thm, L, N), tip=Z[N // 2], L=L, N=N)

def graph(A, xc, lam, nmax):
    da = A.da; it = A.it
    win = int(150 * lam / (A.L / A.N)); sl = slice(max(it - win, 0), it + 1)
    Zw, Tw = A.Z[sl][::-1], -A.T[sl][::-1]; Ms = len(Zw)
    Pp = xc + lam * G1[:nmax]
    def ev(jj, t):
        h = da; t2 = t * t; t3 = t2 * t
        return (2*t3 - 3*t2 + 1)*Zw[jj] + (t3 - 2*t2 + t)*h*Tw[jj] + (-2*t3 + 3*t2)*Zw[jj+1] + (t3 - t2)*h*Tw[jj+1]
    def evd(jj, t):
        h = da; t2 = t * t
        return ((6*t2 - 6*t)*Zw[jj] + (3*t2 - 4*t + 1)*h*Tw[jj] + (-6*t2 + 6*t)*Zw[jj+1] + (3*t2 - 2*t)*h*Tw[jj+1]) / h
    u = np.full(nmax, np.nan)
    u[0] = ((Zw[0] - Pp[0]) * np.conj(nu1[0])).real / lam
    j = 0
    for i in range(1, nmax):
        tt, nv = tau1[i], nu1[i]
        f = lambda jj: ((Zw[jj] - Pp[i]) * np.conj(tt)).real
        jj = max(j - 5, 0); fj = f(jj); found = False
        while jj < Ms - 2:
            fn = f(jj + 1)
            if (fj <= 0 < fn) or (fj >= 0 > fn):
                t = fj / (fj - fn)
                for _ in range(6):
                    g = ((ev(jj, t) - Pp[i]) * np.conj(tt)).real; gp = (evd(jj, t) * np.conj(tt)).real
                    if gp == 0: break
                    dt = -g / gp; t += dt
                    if abs(dt) < 1e-14: break
                if -0.01 <= t <= 1.01:
                    s = ((ev(jj, t) - Pp[i]) * np.conj(nv)).real
                    if abs(s) < 0.5 * lam * (1 + xs[i]): u[i] = s / lam; j = jj; found = True
                break
            jj += 1; fj = fn
            if jj - j > 4000: break
        if not found: break
    return u

def analyse(path, xe_cap=XT, nonlin_start=None, verbose=False, P=None):
    P = P or prep(path)
    S = P['S']; tip = P['tip']; N = P['N']
    kt = abs((2 * np.pi / S.L) * np.gradient(np.angle(np.exp(1j * P['up'].T / abs(P['up'].T) * 1)), 1)[0]) if False else None
    al, th = S.theta(8); kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2 if False else len(al) // 2])
    lam0 = K0 / kt; xc0 = tip.real + D_S1 * lam0 + 1j * tip.imag
    nT = int(round(xe_cap / HX)) + 1
    ug = graph(P['up'], xc0, lam0, nT); ul = graph(P['lo'], xc0, lam0, nT)
    xend = min(xs[:nT][np.isfinite(ug)][-1], xs[:nT][np.isfinite(ul)][-1])
    xe = min(xe_cap, xend - 0.3); n = int(round(xe / HX)) + 1
    sw = np.sqrt(wh[:n]); B = np.column_stack([PSI3[:n], H2[:n]]) * sw[:, None]
    def res(p, ret=False):
        xc = xc0 + p[0] * lam0; lam = lam0 * np.exp(p[1])
        a = graph(P['up'], xc, lam, n); b = graph(P['lo'], xc, lam, n)
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))): return np.full(n, 1e3)
        ue = 0.5 * (a + b); y = ue * sw; c, *_ = np.linalg.lstsq(B, y, rcond=None); r = y - B @ c
        return (r, c, ue, a, b) if ret else r
    p0 = [0.0, 0.0] if nonlin_start is None else nonlin_start
    sol = least_squares(res, p0, x_scale=[1e-3, 1e-3], diff_step=1e-6, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    r, c, ue, a, b = res(sol.x, True)
    xc = xc0 + sol.x[0] * lam0; lam = lam0 * np.exp(sol.x[1])
    w = ue - c[0] * PSI3[:n] - c[1] * H2[:n]
    odd = 0.5 * (a - b)
    return dict(kmax=float(P['z']['kmax']), xc=xc, lam=lam, alpha=c[0], beta2=c[1], n=n, xe=xe, w=w, ue=ue, odd=odd, r2=hn(w, n), lam0=lam0,
                p=sol.x.copy(), tip=tip, nfev=sol.nfev)

if __name__ == '__main__':
    import time
    for f in sys.argv[1:]:
        t = time.time(); d = analyse(f)
        print(f"{f.split('/')[-1]} k={d['kmax']:.1f} lam={d['lam']:.4e} alpha={d['alpha']:+.4e} b2/a2={d['beta2']/d['alpha']**2:+.3f} xe={d['xe']:.2f} "
              f"r2_even={d['r2']:.3e} r2/|a|^3={d['r2']/abs(d['alpha'])**3:.1f} ||odd||={hn(d['odd'], d['n']):.2e} [{time.time()-t:.0f}s]", flush=True)
