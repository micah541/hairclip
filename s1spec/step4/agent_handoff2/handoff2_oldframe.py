"""handoff2: how small is w = u - alpha psi3 - beta2 h2 for snapshots of a CDF hairpin run?  (non-rigorous)
Usage: python3 handoff2.py out.txt snapfile1 snapfile2 ...   (norms: HALF line x>=0, weight rho, trapezoid)"""
import sys, json, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec')
import measure as M
# --- correct S1 geometry: measure.py integrates theta1 by the trapezoid rule (error up to 4e-6 -> ~1e-5 floor in u).
from scipy.interpolate import CubicSpline
M.th1 = -np.pi / 2 + CubicSpline(M.xs, M.k1).antiderivative()(M.xs)
M.tau1 = np.exp(1j * M.th1); M.nu1 = 1j * M.tau1; M.G1 = M.q1 * M.tau1 + M.p1 * M.nu1
M.F2 = -M.chi * np.sin(M.th1)
from measure import xs, HX, rho, tau1, nu1, G1, K0, D_S1, load, strand, modulate
from numpy.polynomial import chebyshev as C
from scipy.optimize import least_squares

def cheb_fn(path):
    d = json.load(open(path)); out = np.full(len(xs), np.nan)
    for p in d['pieces']:
        a, b = p['a'], p['b']; m = (xs >= a - 1e-12) & (xs <= b + 1e-12)
        out[m] = C.chebval(1 - 2 * (xs[m] - a) / (b - a), p['coef'])
    return out
PSI3 = cheb_fn('/home/claude/s1spec/step3nl/psi3md2_14.json'); H2 = cheb_fn('/home/claude/s1spec/step3nl/h2md_14.json')
XT = 28.0; nT = int(round(XT / HX)) + 1
PSI3, H2 = PSI3[:nT], H2[:nT]
wh = (M.wq / 2)[:nT].copy()            # half-line weight (trapezoid ends halved at x=0 and x=40; x=28 end negligible)
wh[0] = HX * rho[0] / 2 * 1.0           # trapezoid end at x=0 (wq/2 already has HX/2*rho at 0)
def hn(f, n=None):
    n = n or nT; return float(np.sqrt(np.sum(wh[:n] * f[:n]**2)))

def hermite_graph(S, xc, lam, nmax, up=8, win=None):
    """u(x), x on the S1 grid, intersect S1 normals with the strand, cubic-Hermite interpolation of the curve (exact tangents)."""
    al, th, Z, _ = S.curve(up); N = len(al); it = N // 2
    T = (S.L / (2 * np.pi)) * np.exp(1j * th)                  # dZ/dalpha
    da = 2 * np.pi / N
    if win is None: win = int(150 * lam / (S.L / N))            # about 150 lam of arclength each side of the tip
    sl = slice(max(it - win, 0), it + 1)
    Zw, Tw = Z[sl][::-1], -T[sl][::-1]                           # strand from the tip, parameter b = -alpha direction
    Ms = len(Zw)
    Pp = xc + lam * G1[:nmax]
    def ev(jj, t):
        h = da; t2 = t * t; t3 = t2 * t
        h00 = 2*t3 - 3*t2 + 1; h10 = t3 - 2*t2 + t; h01 = -2*t3 + 3*t2; h11 = t3 - t2
        return h00*Zw[jj] + h10*h*Tw[jj] + h01*Zw[jj+1] + h11*h*Tw[jj+1]
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
                for _ in range(6):                              # Newton on the Hermite cubic
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

def fit(S, xc0, lam0, xe, cols, nfree=None):
    """variable projection: nonlinear (xc, lam) = (xc0 + dx*lam0, lam0*exp(s)), linear coefficients on cols (list of arrays)."""
    n = int(round(xe / HX)) + 1; sw = np.sqrt(wh[:n]); B = np.column_stack([c[:n] for c in cols]) * sw[:, None]
    def res(p, ret=False):
        u = hermite_graph(S, xc0 + p[0] * lam0, lam0 * np.exp(p[1]), n)
        if not np.all(np.isfinite(u)): return np.full(n, 1e3)
        y = u * sw; c, *_ = np.linalg.lstsq(B, y, rcond=None); r = y - B @ c
        return (r, c, u) if ret else r
    sol = least_squares(res, [0.0, 0.0], x_scale=[1e-3, 1e-3], diff_step=1e-6, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    r, c, u = res(sol.x, True)
    return dict(xc=xc0 + sol.x[0] * lam0, lam=lam0 * np.exp(sol.x[1]), c=c, nrm=float(np.linalg.norm(r)), u=u, n=n)

def fit_constrained(S, xc0, lam0, xe, a0):
    n = int(round(xe / HX)) + 1; sw = np.sqrt(wh[:n])
    def res(p):
        u = hermite_graph(S, xc0 + p[0] * lam0, lam0 * np.exp(p[1]), n)
        if not np.all(np.isfinite(u)): return np.full(n, 1e3)
        return (u - p[2] * PSI3[:n] - p[2]**2 * H2[:n]) * sw
    sol = least_squares(res, [0.0, 0.0, a0], x_scale=[1e-3, 1e-3, max(abs(a0), 1e-5)], diff_step=1e-6, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return dict(alpha=sol.x[2], nrm=float(np.linalg.norm(sol.fun)))

def analyse(path):
    S, z = load(path); Zs, als = strand(S); tip = Zs[0]
    al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    lam0 = K0 / kt; xc0 = tip.real + D_S1 * lam0
    v, u0, ok, r = modulate(Zs, xc0, lam0, 34.0)              # measure.py modulation (<u,f1>=<u,f2>=0), start values
    xc1, lam1 = v
    # chart end: from the Hermite graph at the start modulation
    ug = hermite_graph(S, xc1, lam1, int(round(XT / HX)) + 1)
    fin = np.isfinite(ug); xend = xs[:len(ug)][fin][-1] if fin.any() else 0
    xe = min(XT, xend - 0.3)
    n = int(round(xe / HX)) + 1
    r0 = hn(ug, n)                                           # ||u|| at the measure.py modulation, on [0, xe]
    f1 = fit(S, xc1, lam1, xe, [PSI3])
    f2 = fit(S, xc1, lam1, xe, [PSI3, H2])
    alpha, beta2 = f2['c']
    fc = fit_constrained(S, f2['xc'], f2['lam'], xe, alpha)
    ac = fc['alpha']
    return dict(kmax=float(z['kmax']), t=S.t, lam=f2['lam'], alpha=alpha, beta2=beta2, ratio=beta2 / alpha**2,
                r0=r0, r1=f1['nrm'], r1_alpha=f1['c'][0], r2=f2['nrm'], alphaC=ac, r2c=fc['nrm'], xend=xend, xe=xe,
                lam_tip=lam0, u=f2['u'])

HDR = ("  kmax      alpha      a*kmax   b2/a^2      r0       r1       r2    r2/|a|^3   alphaC   r2C(b2=a^2)  r2C/|aC|^3  chart_end")
def row(d):
    return (f"{d['kmax']:8.1f} {d['alpha']:+.4e} {d['alpha']*d['kmax']:+8.3f} {d['ratio']:+8.3f} {d['r0']:.3e} {d['r1']:.3e} "
            f"{d['r2']:.3e} {d['r2']/abs(d['alpha'])**3:9.3f}  {d['alphaC']:+.3e} {d['r2c']:.3e} {d['r2c']/abs(d['alphaC'])**3:9.3f}  {d['xend']:6.2f}")
if __name__ == "__main__":
    out = sys.argv[1]; files = sys.argv[2:]
    if '--norms' in files:
        print('half-line norms: psi3', hn(PSI3), 'h2', hn(H2), 'psi3(0)', PSI3[0], 'h2(0)', H2[0]); sys.exit()
    with open(out, 'w') as fo:
        fo.write("norms: half-line (x>=0), weight rho, window [0, min(28, chart_end-0.3)]; alpha: psi3(0)=1 normalisation\n" + HDR + "\n")
        for f in files:
            try: d = analyse(f)
            except Exception as e:
                fo.write(f"{f} failed: {e!r}\n"); fo.flush(); continue
            fo.write(row(d) + "\n"); fo.flush(); print(row(d), flush=True)
