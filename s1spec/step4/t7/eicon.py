"""T7 float: E_I constants along the T1 base for a family of weights W = e^{-psi} (T5 (I0), (I2)), v = zeta e'.
Weights: psi' = (1-s) 2P x/(1+x^2) + s psi_far',  s = S((x - Xb)/wb),  psi_far' = 4x/(1+x^2) ('p2') or phi' = x(1+x^2)^(-1/3)/3 ('rho');
Xb = inf: pure power (1+x^2)^(-P) to the far end (T5 D4).
Outputs per (weight, k):
  mu1  = inf mu(x) on the E_I zone (Gaarding function of L1 in W, lint form), split [10.25,16.25] / [16.25,30.25] / beyond; Hardy term of b0'e;
  mu2  = same for the level-2 operator L1^(2) = -d^4 + c2 d^2 + (c1 + 2c2') d + (c0 + 2c1' + c2'');
  a_f, b_f, c_0..c_3  (zeta ramp [10.25,16.25], base drift c1);  gam0 = max_i ||zeta g_i'||_W, gam2 = max_i ||(zeta g_i')''||_W
  (g normalised: p/||chi p||, n/||chi n||, V/||chi V||);  R1, RIs proxies are taken from T1 (see budget.py)."""
import sys, json, numpy as np, sympy as sp
import bgeo
from bgeo import smoothstep_der

X = sp.symbols('X')
def S_(t): return t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3)

class Weight:
    def __init__(self, P, Xb=np.inf, wb=8.0, far='p2'):
        self.P, self.Xb, self.wb, self.far = P, Xb, wb, far
        self.name = f"P={P:g}" + ("" if not np.isfinite(Xb) else f"|{far}@{Xb:g}")
        fp = 2 * sp.nsimplify(P) * X / (1 + X**2)
        ff = (4 * X / (1 + X**2)) if far == 'p2' else X * (1 + X**2) ** sp.Rational(-1, 3) / 3
        d = [sp.lambdify(X, sp.diff(fp, X, n), 'numpy') for n in range(4)]
        df = [sp.lambdify(X, sp.diff(ff, X, n), 'numpy') for n in range(4)]
        self.dp, self.df = d, df
        # psi itself (for W values): psi = P ln(1+x^2) (+ far part by quadrature)
    def psid(self, x):
        """psi', psi'', psi''', psi'''' at x"""
        out = np.array([self.dp[n](x) * np.ones_like(x) for n in range(4)])
        if not np.isfinite(self.Xb): return out
        u = (x - self.Xb) / self.wb; S = [smoothstep_der(u, n) / self.wb**n for n in range(4)]
        S[0] = np.where(u >= 1, 1.0, np.where(u <= 0, 0.0, S[0]))
        f = np.array([self.df[n](x) * np.ones_like(x) for n in range(4)]); dlt = f - out
        res = out.copy()
        for n in range(4):                                  # (s * dlt)^(n)
            for m in range(n + 1):
                from math import comb
                res[n] = res[n] + comb(n, m) * S[m] * dlt[n - m]
        return res
    def W(self, x):
        """W(x) (normalised W(0) = 1) by integration of psi' on a grid"""
        xs = np.concatenate([[0.0], np.atleast_1d(x)]); xg = np.unique(np.concatenate([np.linspace(0, 40, 8001), np.geomspace(40, max(xs.max(), 41), 20001), xs]))
        ps = np.concatenate([[0.0], np.cumsum(0.5 * (self.psid(xg[1:])[0] + self.psid(xg[:-1])[0]) * np.diff(xg))])
        return np.exp(-np.interp(np.atleast_1d(x), xg, ps))

def gaarding_mu(ps, c2, c2p, c2pp, c1, c1p, c0):
    """mu = -A^2/2 - (W A)''/(2W) - Vv  for L = -d^4 + c2 d^2 + c1 d + c0 in the weight W = e^{-psi}"""
    p1, p2, p3, p4 = ps
    w2 = p1**2 - p2; w2p = 2 * p1 * p2 - p3; w2pp = 2 * p2**2 + 2 * p1 * p3 - p4
    rr2 = lambda f, fp, fpp: fpp - 2 * p1 * fp + w2 * f
    A = 2 * w2 - c2; Ap = 2 * w2p - c2p; App = 2 * w2pp - c2pp
    # rr2(w2) needs w2'' (above);  rr2(c2) needs c2', c2''
    Vv = c0 - (c1p - p1 * c1) / 2 + rr2(c2, c2p, c2pp) / 2 - rr2(w2, w2p, w2pp) / 2
    return -A**2 / 2 - rr2(A, Ap, App) / 2 - Vv

def grid_far(g, x0=10.25):
    x0 = min(x0, 10.25)
    ell = g.ell
    xs = np.unique(np.concatenate([np.arange(x0, 40.0, 0.01), np.geomspace(40.0, ell, 40000), [ell]]))
    return xs[xs <= ell]

def g_derivs(s, g):
    """derivatives 0..3 of p, n, V (analytic recursions; K[m] = kappa^(m))"""
    K = s['K']; k, k1, k2, k3, k4, k5 = K[:6]; p, q, n, t1 = s['p'], s['q'], s['n'], s['t1']
    P0 = p; P1 = -k * q; P2 = -k1 * q - k * (1 + k * p); P3 = -k2 * q - 2 * k1 * (1 + k * p) - k * (k1 * p + k * P1)
    N0 = n; N1 = -k * t1; N2 = -k1 * t1 - k**2 * n; N3 = -k2 * t1 - 3 * k * k1 * n + k**3 * t1
    c, b = g.c, g.b
    V = [-k2 + c * P0 + b * N0, -k3 + c * P1 + b * N1, -k4 + c * P2 + b * N2, -k5 + c * P3 + b * N3]
    return [P0, P1, P2, P3], [N0, N1, N2, N3], V

def rho_norm_chi(g, f):
    x = np.linspace(0, 34.25, 30001); s = g.at(x, mmax=2)
    rho = np.exp(-0.25 * (1 + x**2) ** (2 / 3)); chi = bgeo.chi32(x)
    return float(np.sqrt(np.trapezoid(rho * (chi * f(s, g)) ** 2, x)))

def constants(g, wt, xr=(10.25, 16.25)):
    xs = grid_far(g, xr[0]); xs = xs[xs >= xr[0]]; s = g.at(xs, mmax=5); K = s['K']; k, k1, k2 = K[0], K[1], K[2]
    c = g.c; Q = s['Q']
    c2 = -k**2; c2p = -2 * k * k1; c2pp = -2 * k1**2 - 2 * k * k2
    c1 = -8 * k * k1 + Q - c * xs; c1p = -7 * k1**2 - 8 * k * k2 - c
    c0 = -8 * k1**2 - 10 * k * k2
    ps = wt.psid(xs)
    mu = gaarding_mu(ps, c2, c2p, c2pp, c1, c1p, c0)
    # level 2: c1^(2) = c1 + 2 c2', c0^(2) = c0 + 2 c1' + c2''  (c1^(2)' = c1' + 2 c2'')
    k3 = K[3]; c2ppp = -6 * k1 * k2 - 2 * k * k3
    mu2 = gaarding_mu(ps, c2, c2p, c2pp, c1 + 2 * c2p, c1p + 2 * c2pp, c0 + 2 * c1p + c2pp)
    b0p = -10 * k1 * k2 - 4 * k * k3
    zones = {'ramp': (xr[0], xr[1]), 'inner': (xr[1], 30.25), 'mid': (30.25, 1000.0), 'far': (1000.0, np.inf)}
    out = dict(k=g.k, ell=g.ell, weight=wt.name)
    for nm, (a, bnd) in zones.items():
        m = (xs >= a) & (xs <= bnd)
        if m.any():
            i = np.argmin(mu[m]); out[f'mu1_{nm}'] = float(mu[m][i]); out[f'mu1_{nm}_at'] = float(xs[m][i])
            out[f'mu2_{nm}'] = float(np.min(mu2[m]))
            out[f'xb0p_{nm}'] = float(np.max(np.abs(xs[m] * b0p[m])))
    m = xs >= xr[1]
    Peff = np.maximum(xs * ps[0] / 2, 0.5)        # local power: W ~ x^(-2 Peff)
    out['hardy'] = float(np.max(np.abs(xs[m] * b0p[m]) * 2 / (2 * Peff[m] + 1)))
    out['mu1'] = min(out[f'mu1_{nm}'] for nm in zones if f'mu1_{nm}' in out) - out['hardy']
    MUF = float(__import__('os').environ.get('MU1FLOOR', 'nan'))
    if np.isfinite(MUF) and out['mu1'] < MUF:
        # (I0) split: -int W mu v^2 <= -MUF V0^2 + int_bad W (MUF - mu) zeta^2 e'^2 <= -MUF V0^2 + a_mu W1^2
        Wb = wt.W(xs); zb = np.where(xs >= xr[1], 1.0, bgeo.smoothstep((xs - xr[0]) / (xr[1] - xr[0])))
        rhob = np.exp(-0.25 * (1 + xs**2) ** (2 / 3)); bad = mu - out['hardy'] < MUF
        out['a_mu'] = float(np.max(np.where(bad, Wb * (MUF - mu + out['hardy']) * zb**2 / rhob, 0.0))); out['mu1_raw'] = out['mu1']; out['mu1'] = MUF
    else: out['a_mu'] = 0.0
    out['mu2'] = min(out[f'mu2_{nm}'] for nm in zones if f'mu2_{nm}' in out)
    # ---- modulation norms gam0, gam2 (normalised directions), weight W on [10.25, ell]
    W = wt.W(xs)
    zeta = [np.where(xs >= xr[1], 1.0, bgeo.smoothstep((xs - xr[0]) / (xr[1] - xr[0])))]
    for n in range(1, 4): zeta.append(smoothstep_der((xs - xr[0]) / (xr[1] - xr[0]), n) / (xr[1] - xr[0])**n)
    P, N, V = g_derivs(s, g)
    nrm = [rho_norm_chi(g, lambda s_, g_: s_['p']), rho_norm_chi(g, lambda s_, g_: s_['n']), rho_norm_chi(g, lambda s_, g_: s_['Vs'])]
    out['nrm_p'], out['nrm_n'], out['n_V'] = nrm
    I = lambda f: float(np.sqrt(np.trapezoid(W * f**2, xs)))
    gam0, gam2 = [], []
    for D, nm in zip((P, N, V), nrm):
        d1 = [D[1] / nm, D[2] / nm, D[3] / nm]                      # g', g'', g''' (normalised)
        gam0.append(I(zeta[0] * d1[0]))
        z2 = zeta[2] * d1[0] + 2 * zeta[1] * d1[1] + zeta[0] * d1[2]  # (zeta g')''
        gam2.append(I(z2))
    out['gam0_list'], out['gam2_list'] = gam0, gam2; out['gam0'], out['gam2'] = max(gam0), max(gam2)
    # W-Gram matrices of zeta g_i' and (zeta g_i')'' for the UNnormalised g = (p, n, Vs) (basis-free modulation constants with cx.json)
    D1 = [D[1] for D in (P, N, V)]; D2 = [zeta[2] * D[1] + 2 * zeta[1] * D[2] + zeta[0] * D[3] for D in (P, N, V)]
    out['Gam0'] = [[float(np.trapezoid(W * zeta[0]**2 * a * b, xs)) for b in D1] for a in D1]
    out['Gam2'] = [[float(np.trapezoid(W * a * b, xs)) for b in D2] for a in D2]
    # where gam0 lives (V-direction): share from x > 1000
    mfar = xs > 0.5 * g.ell
    out['gam0_V_bend_share'] = float(np.sqrt(np.trapezoid((W * (zeta[0] * V[1] / nrm[2])**2)[mfar], xs[mfar])) / max(gam0[2], 1e-300))
    # ---- zeta-ramp constants (base drift c1, c2 ~ 0 on the ramp)
    mr = (xs >= xr[0]) & (xs <= xr[1]); xq = xs[mr]
    Wq = W[mr]; rho = np.exp(-0.25 * (1 + xq**2) ** (2 / 3))
    Lr = xr[1] - xr[0]; u = (xq - xr[0]) / Lr
    z = [bgeo.smoothstep(u)] + [smoothstep_der(u, n) / Lr**n for n in range(1, 5)]
    psr = ps[:, mr]
    # derivatives of gj = W zeta zeta^(j):  use dense finite differences on the ramp grid (dx = 0.01) of the exact values
    gj = [Wq * z[0] * z[j] for j in range(5)]
    def der(f, n):
        for _ in range(n): f = np.gradient(f, xq, edge_order=2)
        return f
    c1r = c1[mr]
    Az = -2 * der(gj[1], 3) + 3 * der(gj[2], 2) - 2 * der(gj[3], 1) + gj[4] - c1r * gj[1]
    Bz = 6 * der(gj[1], 1) - 6 * gj[2]
    out['a_f'] = float(np.max(np.maximum(Az, 0) / rho)) + out['a_mu']; out['b_f'] = float(np.max(np.maximum(Bz, 0) / rho))
    cf = [z[4] - c1r * z[1], 4 * z[3], 6 * z[2], 4 * z[1]]
    for j in range(4): out[f'c{j}'] = float(np.max(np.abs(cf[j]) * np.sqrt(Wq / rho)))
    out['Psi1'] = float(np.max(np.abs(ps[0][xs >= xr[0]])))
    out['Psi2'] = float(np.max(np.abs((ps[0]**2 - ps[1])[xs >= xr[0]])))
    out['xr'] = list(xr)
    return out

if __name__ == '__main__':
    B = bgeo.base()
    KS = [(20.4, 0, 0.5), (40, None, None), (200, None, None), (2e3, None, None), (2e4, None, None)]
    WTS = [Weight(2), Weight(1.5), Weight(1), Weight(1.5, 60, 8, 'p2'), Weight(1, 60, 8, 'p2'), Weight(1, 40, 8, 'rho'), Weight(2, 40, 8, 'rho')]
    XR = tuple(float(v) for v in sys.argv[1].split(',')) if len(sys.argv) > 1 else (10.25, 16.25)
    if len(sys.argv) > 2: WTS = [w for w in WTS if w.name in sys.argv[2].split(';')]
    res = []
    for kk, i, th in KS:
        if i is None: tau, i = B.tau_of_k(kk); th = (tau - B.tau0[i]) / B.dt[i]
        g = bgeo.BaseGeom(i, th)
        for wt in WTS:
            o = constants(g, wt, XR); o['slab'] = i; o['theta'] = th; res.append(o)
            print(f"k={g.k:8.1f} {wt.name:12s} mu1={o['mu1']:+.4f} (ramp {o['mu1_ramp']:+.3f} inner {o['mu1_inner']:+.3f} mid {o.get('mu1_mid', np.nan):+.3f}"
                  f" far {o.get('mu1_far', np.nan):+.3f}; hardy {o['hardy']:.1e}) mu2={o['mu2']:+.3f} a_f={o['a_f']:.3e} b_f={o['b_f']:.3e} "
                  f"c=({o['c0']:.3f},{o['c1']:.3f},{o['c2']:.3f},{o['c3']:.3f}) gam0={o['gam0']:.3e} gam2={o['gam2']:.3e} (bend share {o['gam0_V_bend_share']:.2f}) nV={o['n_V']:.3e}", flush=True)
    json.dump(res, open('eicon.json' if XR == (10.25, 16.25) else f'eicon_{XR[0]:g}_{XR[1]:g}.json', 'w'), indent=1)
