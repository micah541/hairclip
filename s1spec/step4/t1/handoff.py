"""Milestone-1 hand-off quantities computed directly from the T1 base (exactly symmetric, rescaled frame), reusing step4/m1 code:
modulation fit (xc, lam, alpha) as m1/prof.psi_fit, Lambda by m1/en.inner, shells/windows by m1/en.shell (U = alpha psi3' + alpha^2 h2'
as in M1, and the v3 choice U = 1.3056 alpha y^1.0246), cone jets by m1/jets.cyl_jets (centre r = 200, L = r/4, X = 3).
usage: python3 handoff.py TAG k1 [k2 ...]   -> handoff_<TAG>.json"""
import sys, json, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/m1')
from numpy.polynomial import legendre as LG
from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline
import geo, en, jets
from geo import xg, HX, G1G, tauG, nuG, kG, PSI0, PSI1, H20, H21, THINF, CH
from base2 import Base
from rep import K0, DS1, STIP

NQ = 40
gq, wq = LG.leggauss(NQ); UQ = (gq + 1) / 2

class Curve:
    """rescaled base curve at one time: theta, Gamma, kappa jets at arbitrary x (per-element Legendre series)"""
    def __init__(self, R, Y):
        self.R = R; self.Y = Y; n = R.n; self.ell = Y[n]
        x, _ = R.mesh.knots(self.ell); self.x = x; h = np.diff(x); self.h = h
        d = R.local(Y[:n], h); self.d = d
        B = R.B
        Vq = B.evalpts(UQ, 0)
        kq = d @ Vq.T                                              # kappa at Gauss points
        # theta: exact cumulative (Legendre series of kappa per element, integrated)
        cK = np.array([LG.legfit(gq, kq[e], NQ - 1) for e in range(R.N)])
        self.cTh = []; th0 = -np.pi / 2; self.th_e = [th0]
        for e in range(R.N):
            ci = LG.legint(cK[e], lbnd=-1) * (h[e] / 2); ci[0] += th0
            self.cTh.append(ci); th0 = LG.legval(1.0, ci); self.th_e.append(th0)
        self.cE = []; G0 = -DS1 + 0j; self.G_e = [G0]
        for e in range(R.N):
            thq = LG.legval(gq, self.cTh[e]); Eq = np.exp(1j * thq)
            cE = LG.legfit(gq, Eq.real, NQ - 1) + 1j * LG.legfit(gq, Eq.imag, NQ - 1)
            ci = LG.legint(cE, lbnd=-1) * (h[e] / 2); ci[0] += G0
            self.cE.append(ci); G0 = LG.legval(1.0, ci); self.G_e.append(G0)
        self.cTh = np.array(self.cTh); self.cE = np.array(self.cE)

    def elem(self, xs):
        xs = np.asarray(xs, float); e = np.clip(np.searchsorted(self.x, xs, side='right') - 1, 0, self.R.N - 1)
        t = 2 * (xs - self.x[e]) / self.h[e] - 1; return e, t

    def theta(self, xs):
        e, t = self.elem(xs); out = np.empty(len(t))
        for ee in np.unique(e): m = e == ee; out[m] = LG.legval(t[m], self.cTh[ee])
        return out

    def Gam(self, xs):
        e, t = self.elem(xs); out = np.empty(len(t), complex)
        for ee in np.unique(e): m = e == ee; out[m] = LG.legval(t[m], self.cE[ee])
        return out

    def kap(self, xs, m=0):
        return self.R.eval_at(self.Y[:self.R.n], self.ell, np.asarray(xs, float), m)


def graph_over_S1(C, xc, lam, nmax, slope=False):
    """normal graph u(x_i) of the base (frame (xc, lam)) over S1 on the 0.01 grid, Newton on the foot point"""
    P = xc + lam * G1G[:nmax]; tt = tauG[:nmax]; nv = nuG[:nmax]
    xb = xg[:nmax].copy() * lam                                  # base arclength ~ lam * S1 arclength
    for _ in range(30):
        Gb = C.Gam(xb); Tb = np.exp(1j * C.theta(xb))
        g = ((Gb - P) * np.conj(tt)).real; gp = (Tb * np.conj(tt)).real
        dx = -g / gp; xb += dx
        if np.max(np.abs(dx)) < 1e-14: break
    Gb = C.Gam(xb); u = ((Gb - P) * np.conj(nv)).real / lam
    if not slope: return u
    ph = C.theta(xb) - np.angle(tt)
    up = (1 - kG[:nmax] * u) * np.tan(np.angle(np.exp(1j * ph)))
    return u, up


def analyse(C, ymax_arc=3000.0):
    n = int(round(28.0 / HX)) + 1
    wh = geo.H.wh[:n] if hasattr(geo.H, 'wh') else None
    import chart as CHm
    wh = CHm.H.wh[:n]; sw = np.sqrt(wh)
    def res(p, ret=False):
        xc = p[0]; lam = np.exp(p[1])
        ue = graph_over_S1(C, xc, lam, n)
        r = (ue - p[2] * PSI0[:n] - p[2]**2 * H20[:n]) * sw
        return (r, ue) if ret else r
    a0 = -3.3 / (K0 / 1.0) * 0 - 1e-4
    sol = least_squares(res, [0.0, 0.0, a0], x_scale=[1e-3, 1e-3, 1e-4], diff_step=1e-6, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    xc, lam, alpha = sol.x[0], np.exp(sol.x[1]), sol.x[2]
    ncr = int(round((24.25 + 2.0) / HX)) + 1
    u, up = graph_over_S1(C, xc, lam, ncr, slope=True)
    out = dict(xc=xc, lamfit=lam, alpha=alpha, r2=float(np.linalg.norm(res(sol.x))), x_crawl=xg[:ncr], u_crawl=u, up_crawl=up)
    # ray profile along the strand (base arclength sampled at ds = 0.1 lam), frame (xc, lam), rotated by THINF
    xs = np.arange(0.0, min(ymax_arc * lam, C.ell * 0.999), 0.1 * lam)
    G = (C.Gam(xs) - xc) / lam; T = np.exp(1j * C.theta(xs))
    c = G * np.exp(-1j * THINF); q = c.real; uu = c.imag
    dd = T * np.exp(-1j * THINF); f = dd.imag / dd.real
    ok = q > 12.0; i0 = np.argmax(ok); dq = np.diff(q[i0:]); bad = np.where(dq <= 0)[0]
    iend = i0 + (bad[0] if len(bad) else len(dq))
    q, uu, f = q[i0:iend + 1], uu[i0:iend + 1], f[i0:iend + 1]
    yy = np.exp(np.arange(np.log(q[0] + 1e-9), np.log(q[-1]), 0.002))
    out.update(y=yy, f_even=CubicSpline(q, f)(yy), uy_even=CubicSpline(q, uu)(yy))
    big = np.where(np.abs(f) > 0.5)[0]; out['chart_ext'] = float(q[big[0]]) if len(big) else float(q[-1])
    return out


def cone(C, xc, lam, rc=200.0, X=3.0):
    """cone-cylinder jets K2, K3, K4 (cylinder units) and slope at centre radius rc (frame (xc, lam))"""
    ds = 0.02
    xs = np.arange(0.0, min(3000.0 * lam, C.ell * 0.999), ds * lam)
    G = (C.Gam(xs) - xc) / lam
    A = dict(phi=C.theta(xs) + np.pi, D1=C.kap(xs, 0) * lam, D2=C.kap(xs, 1) * lam**2, D3=C.kap(xs, 2) * lam**3, G=G, ds=ds, N=len(xs))
    r = np.abs(G); jc = int(np.argmax(r >= rc)); Lc = 0.25 * r[jc]
    dj = jets.cyl_jets(A, jc, Lc, X)
    return dict(r=float(r[jc]), L=float(Lc), slope=dj['hx'], K2=Lc * dj['h2'], K3=Lc**2 * dj['h3'], K4=Lc**3 * dj['h4'], N3=Lc**1.5 * dj['n3'])


def quantities(B, kk):
    tau, i = B.tau_of_k(kk); R, D = B.rep(i); Y = B.Y(tau, i)
    C = Curve(R, Y)
    d = analyse(C)
    d['lam'] = d['lamfit']; d['kmax'] = kk
    I = en.inner(d)
    Bb = en.build(d)
    out = dict(k=kk, tau=tau, alpha=d['alpha'], alpha_k=d['alpha'] * kk, xc=d['xc'], lamfit=d['lamfit'], fit_r2=d['r2'],
               A=I['A'], sqrtLam=I['sqrtLam'], sqrtE=I['sqrtE'], sqrtmuEI=I['sqrtmuEI'], chart_ext=d['chart_ext'], ymax=float(Bb['ymax']))
    for lbl, U in (('M1', Bb['U']), ('v3', 1.3056 * d['alpha'] * Bb['yg']**1.0246)):
        Bv = dict(Bb); Bv['U'] = U; Bv['v'] = Bb['f'] - U
        sh, slide, z = en.shell(Bv)
        out[f'shell14_60_{lbl}'] = sh[(14, 60)]
        out[f'win_{lbl}'] = {f'{a}-{b}': v for (a, b), v in sh.items() if a >= 50}
        out[f'slide_{lbl}'] = [float(slide[j][0]) for j in range(5)]
        out[f'slide_at_{lbl}'] = [slide[j][1] for j in range(5)]
    out['v_over_cubic'] = {str(yy): float(-(Bb['f'] - Bb['U'])[np.argmin(abs(Bb['yg'] - yy))] / (yy / kk) ** 3) for yy in (50, 112, 253, 350)}
    out['cone'] = cone(C, d['xc'], d['lamfit'])
    return out


if __name__ == '__main__':
    tag = sys.argv[1]; B = Base(tag); res = []
    for kk in map(float, sys.argv[2:]):
        o = quantities(B, kk); res.append(o)
        print(json.dumps({k: (v if not isinstance(v, (list, dict)) else v) for k, v in o.items()}, default=float)[:1800], flush=True)
    json.dump(res, open(f'handoff_{tag}.json', 'w'), indent=1, default=float)
