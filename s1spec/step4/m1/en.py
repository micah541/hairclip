"""Stage 2: hand-off quantities from the stage-1 profile file."""
import sys, numpy as np
from math import factorial
from geo import *
from scipy.interpolate import CubicSpline
DS = 0.002; MM = 200; DEG = 12
def stencils(m=MM, deg=DEG, jmax=4):
    t = np.arange(-m, m + 1) * DS; hw = m * DS
    P = np.linalg.pinv(np.vander(t / hw, deg + 1, increasing=True))
    return np.array([P[j] * factorial(j) / hw**j for j in range(jmax + 1)])
OPS = stencils()
EPS, EPSI, MU = 1e-5, 3e-5, 0.08
WIN = [(14, 60), (50, 75), (75, 112), (112, 169), (169, 253)]
THR = dict(E0=2e-7, E1=1e-4, E2=0.2)   # on [14,60]
THW = 3e-6                              # windows [50,253], levels 1..3 (and E4 at tau1); E0 <= 2e-7

def d4(f, h):          # second derivative, 4th order, interior; ends 2nd-order
    g = np.empty_like(f); g[2:-2] = (-f[4:] + 16 * f[3:-1] - 30 * f[2:-2] + 16 * f[1:-3] - f[:-4]) / (12 * h * h)
    g[:2] = g[2]; g[-2:] = g[-3]; return g
def d1(f, h):
    g = np.empty_like(f); g[2:-2] = (-f[4:] + 8 * f[3:-1] - 8 * f[1:-3] + f[:-4]) / (12 * h)
    g[:2] = g[2]; g[-2:] = g[-3]; return g

def build(d):
    """returns dict with v(s) grid etc."""
    al = float(d['alpha']); lam = float(d['lam'])
    y = d['y']; fe = d['f_even']
    xc_ = d['x_crawl']; upc = d['up_crawl']
    # f(y) function: crawl for y<22 (x < 24.249), ray beyond
    ymax = y[-1]
    s0 = np.log(8.6); s1 = np.log(ymax) - 1e-9
    sg = np.arange(s0, s1, DS); yg = np.exp(sg)
    fray = CubicSpline(np.log(y), fe)(np.minimum(sg, np.log(ymax)))
    fcr = CubicSpline(xc_, upc)(np.clip(yg - STIP, 0, xc_[-1]))
    use_cr = yg < 22.0
    f = np.where(use_cr, fcr, fray)
    a, b = c4_first(yg); U = al * a + al**2 * b
    v = f - U
    # crawl/ray consistency in the overlap
    ov = (yg > 14) & (yg < 22)
    cons = float(np.max(np.abs(fcr[ov] - fray[ov]))) if np.any(ov & (yg > y[0])) else np.nan
    ovm = ov & (yg > y[0]); cons = float(np.max(np.abs(fcr[ovm] - fray[ovm])))
    return dict(sg=sg, yg=yg, f=f, U=U, v=v, cons=cons, ymax=ymax)

def wint(f, sg, sa, sb):
    """trapezoid of f over [sa,sb] on the uniform grid sg with exact (linearly interpolated) end points"""
    m = (sg > sa) & (sg < sb)
    xs_ = np.concatenate([[sa], sg[m], [sb]]); ys_ = np.concatenate([[np.interp(sa, sg, f)], f[m], [np.interp(sb, sg, f)]])
    return float(np.trapezoid(ys_, xs_))

def shell(B, ymax_ok=None):
    sg, yg, v = B['sg'], B['yg'], B['v']
    n = len(sg); m = MM
    z = np.full((5, n), np.nan)
    for j in range(5):
        z[j, m:n - m] = np.convolve(v, OPS[j][::-1], mode='valid')
    om = 28.0 / yg
    out = {}
    for (ya, yb) in WIN:
        sel = (yg >= ya) & (yg <= yb) & np.isfinite(z[4])
        if sel.sum() < 10 or yg[sel][-1] < 0.995 * yb or yg[sel][0] > 1.005 * ya:
            out[(ya, yb)] = [np.nan] * 5; continue
        out[(ya, yb)] = [wint(om * np.nan_to_num(z[j])**2, sg, np.log(ya), np.log(yb)) for j in range(5)]
    # sliding windows [y,1.5y], y in [50,253]
    slide = {j: (0.0, None) for j in range(5)}
    for ya in np.exp(np.arange(np.log(50), np.log(253) + 1e-9, 0.02)):
        yb = 1.5 * ya
        sel = (yg >= ya) & (yg <= yb) & np.isfinite(z[4])
        if sel.sum() < 10 or yg[sel][-1] < 0.995 * yb: continue
        for j in range(5):
            e = wint(om * np.nan_to_num(z[j])**2, sg, np.log(ya), np.log(yb))
            if e > slide[j][0]: slide[j] = (e, ya)
    return out, slide, z

def inner(d):
    """A, B, V0, V2, Lambda on the x grid (half-line)."""
    al = float(d['alpha'])
    y = d['y']; xc_ = d['x_crawl']
    NXc = len(xc_)
    sg = np.log(y)
    spf = CubicSpline(sg, d['f_even']); spu = CubicSpline(sg, d['uy_even'])
    up = np.empty(NX); u = np.empty(NX)
    up[:NXc] = d['up_crawl']; u[:NXc] = d['u_crawl']
    yy = xg[NXc:] + STIP
    ok = yy <= y[-1]
    up[NXc:] = np.where(ok, spf(np.log(np.minimum(yy, y[-1]))), np.nan); u[NXc:] = np.where(ok, spu(np.log(np.minimum(yy, y[-1]))), np.nan)
    wp = up - al * PSI1 - al**2 * H21          # w'
    w = u - al * PSI0 - al**2 * H20
    # w'' (odd extension of w' about 0), w''' 
    wpe = np.concatenate([-wp[1:3][::-1], wp]); wpp = d1(wpe, HX)[2:]; wpp[0] = 0.0   # w''(0)=? from symmetry w'(0)=0, w''(0) generally nonzero
    # better at x=0: w''(0) = (-wp[2]+16wp[1]... ) one-sided via 4th order: use odd extension derivative (valid, wp odd)
    wp3 = d1(wpp, HX)
    chi = chi32
    # chi derivatives
    chip = d1(chi, HX); chipp = d2 = d4(chi, HX)
    wt = chi * w; wtpp = chipp * w + 2 * chip * wp + chi * wpp
    wh = np.full(NX, HX) * rho; wh[0] *= 0.5
    m34 = xg <= 34.25
    A = float(np.sqrt(np.sum((wh * wt**2)[m34]))); Bn = float(np.sqrt(np.sum((wh * wtpp**2)[m34])))
    V = zeta * wp
    Vp = d1(V, HX); Vpp = d4(V, HX)
    wW = np.full(NX, HX) * Wt; wW[0] *= 0.5
    fin = np.isfinite(V)
    V0 = float(np.sqrt(np.sum((wW * np.where(fin, V, 0)**2)))); V2 = float(np.sqrt(np.sum(wW * np.where(fin, Vpp, 0)**2)))
    E = A**2 + EPS * Bn**2; EI = V0**2 + EPSI * V2**2; Lam = E + MU * EI
    sel = (xg >= 16.25) & (xg <= 57) & fin
    return dict(A=A, B=Bn, V0=V0, V2=V2, E=E, EI=EI, sqrtLam=float(np.sqrt(Lam)), sqrtE=float(np.sqrt(E)), sqrtmuEI=float(np.sqrt(MU * EI)),
                supwp=float(np.max(np.abs(wp[sel]))), supw_in=float(np.max(np.abs(w[xg <= 16.25]))), wp=wp, w=w, wpp=wpp)

if __name__ == '__main__':
    d = dict(np.load(sys.argv[1]))
    I = inner(d); print({k: v for k, v in I.items() if not isinstance(v, np.ndarray)})
    B = build(d); print("cons crawl/ray", B['cons'], "ymax", B['ymax'])
    out, slide, z = shell(B)
    for k, v in out.items(): print(k, ["%.3e" % e for e in v])
    print("slide", {j: ("%.2e" % e, ya) for j, (e, ya) in slide.items()})
