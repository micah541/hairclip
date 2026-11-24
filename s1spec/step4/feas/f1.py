"""F1: Chebyshev spectral-element DOF for theta(sigma) of the half curve, geometric partition.
usage: python3 f1.py K [K ...]   (uses design/snap_clip_24_0.4_0.05_0.85_K.npz)"""
import sys, json, time, numpy as np
from scipy.fft import dct
from numpy.polynomial import chebyshev as C
from common import *
P = 128                                  # reference degree (interpolant used as 'truth')
TOLS = [1e-10, 1e-11, 1e-12, 1e-13]

def partition(r, send):
    e = [0., 1., 2., 3., 4.]
    while e[-1] < send: e.append(e[-1] * (1 + r))
    e[-1] = send
    return np.array(e)

def cheb_coeffs(f):                      # f at Chebyshev-Lobatto x_m = cos(pi m/P), m=0..P
    c = dct(f, type=1) / (len(f) - 1); c[0] /= 2; c[-1] /= 2
    return c

def analyse(S, lam, edges, tag=''):
    sig_end = edges[-1]
    xm = np.cos(np.pi * np.arange(P + 1) / P)
    # node values theta(sigma) on the half curve, reversed orientation: alpha = pi - 2 pi lam sigma / L
    sg = np.concatenate([0.5 * (edges[i] + edges[i + 1]) + 0.5 * (edges[i + 1] - edges[i]) * xm for i in range(len(edges) - 1)])
    al = np.pi - 2 * np.pi * lam * sg / S.L
    th = np.pi + al + S.th_derivs(al, 0)             # theta + pi (reversed orientation), additive constants irrelevant
    th = th.reshape(len(edges) - 1, P + 1)
    res = []
    for i in range(len(edges) - 1):
        c = cheb_coeffs(th[i]); h = edges[i + 1] - edges[i]
        ac = np.abs(c)
        rec = dict(a=edges[i], b=edges[i + 1], c=ac)
        # fourth derivative of the degree-P interpolant (truth), in sigma units
        d4 = C.chebder(c, 4) * (2 / h) ** 4
        xf = np.cos(np.linspace(0, np.pi, 400)); rec['d4sup'] = float(np.max(np.abs(C.chebval(xf, d4))))
        res.append(rec)
    out = {}
    for tol in TOLS:
        ps = []
        for rec in res:
            idx = np.where(rec['c'] >= tol)[0]
            ps.append(int(idx.max()) if len(idx) else 0)     # smallest p with all c_j<tol for j>p
        ps = np.array(ps)
        out[tol] = dict(p=ps.tolist(), dof=int(np.sum(ps + 1)), pmax=int(ps.max()), saturated=bool(np.any(ps >= P - 2)))
    # sup error of theta'''' when truncating at p(1e-13)  and  DOF needed s.t. d4 err <= tol4
    p13 = out[1e-13]['p']; errs = []
    for rec, p, e in zip(res, p13, range(len(res))):
        h = rec['b'] - rec['a']; c = C.Chebyshev(np.concatenate([np.zeros(p + 1), np.zeros(0)]))
        ct = np.zeros(P + 1); ct[p + 1:] = rec['c'] * 0  # placeholder (sign irrelevant for sup bound below)
        errs.append(None)
    return res, out, sg, th

NOISE = 3e-15    # Chebyshev coefficients below this are rounding noise of the theta values (~1e-15)
_J = np.arange(P + 1, dtype=float)
_T4 = _J**2 * (_J**2 - 1) * (_J**2 - 4) * (_J**2 - 9) / 105.0     # sup |T_j^(4)| on [-1,1]
def d4_err_for_p(c_signed, h, p):
    """bound sum |c_j| sup|T_j^(4)| (2/h)^4 over the truncated tail j>p, ignoring coefficients below the noise floor"""
    a = np.abs(c_signed).copy(); a[:p + 1] = 0; a[a < NOISE] = 0
    return float(np.sum(a * _T4) * (2 / h) ** 4)

def p_extrap(c, h, t4, fitlo=1e-12, fithi=1e-5):
    """degree p such that the extrapolated (geometric) coefficient tail satisfies  sum_{j>p} A e^{-g j} sup|T_j^(4)| (2/h)^4 <= t4.
    A, g fitted to the coefficient envelope in the window fitlo<|c|<fithi; falls back to the noise-cut result if the fit is not possible."""
    a = np.abs(c); env = np.maximum(a[:-1], a[1:]); j = np.arange(len(env))
    m = (env > fitlo) & (env < fithi)
    jmax = j[env > NOISE].max() if np.any(env > NOISE) else 0
    if m.sum() < 4: return None
    g, lnA = np.polyfit(j[m], np.log(env[m]), 1)
    if g > -0.05: return None
    for p in range(0, 400):
        jj = np.arange(p + 1, 600.0)
        tail = np.exp(lnA + g * jj) * jj**2 * (jj**2 - 1) * (jj**2 - 4) * (jj**2 - 9) / 105.0 * (2 / h) ** 4
        if tail.sum() <= t4: return p
    return None

def run(K, rs=(0.25, 0.5, 1.0), snap=None):
    K = int(K)
    S = Snap(snap or f'/home/claude/s1spec/step4/design/snap_clip_24_0.4_0.05_0.85_{K}.npz')
    km = S.kmax_stored; lam = K0S1 / km; send = 0.5 * S.L / lam
    rep = dict(K=K, kmax=km, lam=lam, sigma_end=send, N=S.N, L=S.L, runs={})
    for r in rs:
        t0 = time.time(); edges = partition(r, send)
        xm = np.cos(np.pi * np.arange(P + 1) / P)
        sg = np.concatenate([0.5 * (edges[i] + edges[i + 1]) + 0.5 * (edges[i + 1] - edges[i]) * xm for i in range(len(edges) - 1)])
        dd = -2 * np.pi * lam * sg / S.L
        th = (2 * np.pi + dd + S.th_derivs_d(dd, 0)).reshape(len(edges) - 1, P + 1)
        cs = [cheb_coeffs(th[i]) for i in range(len(edges) - 1)]
        hs = np.diff(edges)
        ent = dict(nel=len(hs), tol={})
        for tol in TOLS:
            ps = np.array([int(np.where(np.abs(c) >= max(tol, NOISE))[0].max()) if np.any(np.abs(c) >= max(tol, NOISE)) else 0 for c in cs])
            e4 = max(d4_err_for_p(c, h, p) for c, h, p in zip(cs, hs, ps))
            # per-region breakdown
            core = float(np.sum((ps + 1)[edges[1:] <= 4 + 1e-9]))
            ent['tol'][f'{tol:.0e}'] = dict(dof=int(np.sum(ps + 1)), pmax=int(ps.max()), pmin=int(ps.min()), d4err_sup=e4, dof_core=core,
                                           sat=bool(ps.max() >= P - 2))
        # DOF needed so that sup err of theta'''' <= t4 on every element
        for t4 in (1e-8, 1e-10, 1e-12):
            tot = 0; pm = 0
            for c, h in zip(cs, hs):
                p = 0
                while p < P - 1 and d4_err_for_p(c, h, p) > t4: p += 1
                tot += p + 1; pm = max(pm, p)
            ent[f'dof_d4_{t4:.0e}'] = tot; ent[f'pmax_d4_{t4:.0e}'] = pm
        for t4 in (1e-8, 1e-10, 1e-12):
            ps = [p_extrap(c, h, t4) for c, h in zip(cs, hs)]
            if any(p is None for p in ps):   # elements without a usable geometric fit (very smooth / flat): use p13 there
                p13 = [int(np.where(np.abs(c) >= max(1e-13, NOISE))[0].max()) if np.any(np.abs(c) >= 1e-13) else 0 for c in cs]
                ps = [p if p is not None else q for p, q in zip(ps, p13)]
            ent[f'dof_d4x_{t4:.0e}'] = int(sum(p + 1 for p in ps)); ent[f'pmax_d4x_{t4:.0e}'] = int(max(ps))
        ent['p13_list'] = [int(np.where(np.abs(c) >= 1e-13)[0].max()) if np.any(np.abs(c) >= 1e-13) else 0 for c in cs]
        ent['edges'] = edges.tolist()
        ent['d4sup_ref'] = float(max(d4_err_for_p(c, h, -1) for c, h in zip(cs, hs)))
        ent['secs'] = time.time() - t0
        rep['runs'][str(r)] = ent
        print(f"K={K} r={r} nel={ent['nel']} " + ' '.join(f"tol{t}: DOF={v['dof']} pmax={v['pmax']} d4err={v['d4err_sup']:.1e}{'*' if v['sat'] else ''}" for t, v in ent['tol'].items())
              + f" | DOF(d4<=1e-8/-10/-12): noise-cut {ent['dof_d4_1e-08']}/{ent['dof_d4_1e-10']}/{ent['dof_d4_1e-12']} extrap {ent['dof_d4x_1e-08']}/{ent['dof_d4x_1e-10']}/{ent['dof_d4x_1e-12']} d4sup={ent['d4sup_ref']:.2e} [{ent['secs']:.0f}s]", flush=True)
    json.dump(rep, open(f'/home/claude/s1spec/step4/feas/logs/f1_{K}.json', 'w'), indent=1)
    return rep

if __name__ == '__main__':
    for a in sys.argv[1:]:
        if '=' in a: K, p = a.split('='); run(K, snap=p)
        else: run(a)
