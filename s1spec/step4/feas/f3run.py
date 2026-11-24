"""F3: tail dimension of sym_rho(L_Gamma + transport).  usage: python3 f3run.py [--frame tip|fit] [--notau] [--sigma 35] K ...
(i) #eig of sym_rho above -c;  (ii) smallest M s.t. max eig of sym_rho on the rho-orthogonal complement of the first M S1-sym eigenvectors <= -c."""
import sys, time, json, argparse, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/feas')
from f3geom import *; from f3fem import *
ap = argparse.ArgumentParser(); ap.add_argument('--frame', default='tip'); ap.add_argument('--notau', action='store_true'); ap.add_argument('--sigma', type=float, default=35.0)
ap.add_argument('--h', type=float, default=0.0); ap.add_argument('--eps2', type=float, default=0.0); ap.add_argument('--tag', default=''); ap.add_argument('ks', nargs='*')
a = ap.parse_args()
CS = (0.05, 0.1, 0.25, 0.5, 1.0)

def minM(U, S_G, c, Mmax=None):
    """smallest M s.t. max eig of U[:,M:]^T S U[:,M:] <= -c (monotone in M by interlacing -> bisection)"""
    SU = U.T @ S_G @ U; n = SU.shape[0]
    top = lambda M: np.linalg.eigvalsh(SU[M:, M:])[-1]
    lo, hi = 0, n - 1
    if top(0) <= -c: return 0, top(0)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if top(mid) <= -c: hi = mid
        else: lo = mid
    return hi, top(hi)

def hfor(k): return a.h if a.h > 0 else (0.1 if k >= 100 else (0.05 if k >= 30 else 0.025))

def analyse_one(label, sp, mesh, M, G, K1, U1, ws1):
    K = mesh.form(coeffs(mesh, sp, use_tau=not a.notau))
    w, V = eig_sym(K, M, G=G)
    S = 0.5 * (K + K.T)
    out = dict(label=label, top=[float(x) for x in w[:8]])
    for c in CS:
        out[f'n_{c}'] = int(np.sum(w > -c))
        Mc, tv = minM(U1, S, c); out[f'M_{c}'] = int(Mc); out[f'maxres_{c}'] = float(tv)
    return out

res = []
Sigma = a.sigma
ref = {}
for k in ['S1'] + a.ks:
    t0 = time.time(); kk = 100 if k == 'S1' else int(k); h = hfor(kk)
    mesh = Mesh(Sigma, h)
    key = h
    if key not in ref:
        sp1 = s1_splines(Smax=Sigma + 5, dx=0.01); M = mesh.mass()
        K1 = mesh.form(coeffs(mesh, sp1, use_tau=False))
        G = M if a.eps2 == 0 else M + a.eps2 * mesh.stiff2()
        dsc = 1 / np.sqrt(np.diag(M))
        if a.eps2 > 0:   # H^2-type inner product G = M + eps2*K2 ; sym part wrt G of the L^2-Galerkin operator A = M^-1 K (computed in Jacobi-scaled variables)
            def Kg(Kx):
                Ms = M * dsc[:, None] * dsc[None, :]; Ks = Kx * dsc[:, None] * dsc[None, :]; Gs = G * dsc[:, None] * dsc[None, :]
                GA = Gs @ np.linalg.solve(Ms, Ks); Ss = 0.5 * (GA + GA.T)
                return Ss / (dsc[:, None] * dsc[None, :])
        else:
            Kg = lambda Kx: Kx
        # S1 reference eigenvectors w.r.t. G
        if a.eps2 > 0:
            S1m = Kg(K1); d = 1 / np.sqrt(np.diag(G)); w1, v1 = sla.eigh(S1m * d[:, None] * d[None, :], G * d[:, None] * d[None, :]); o = np.argsort(-w1); ws1, U1 = w1[o], (v1 * d[:, None])[:, o]
        else:
            ws1, U1 = eig_sym(K1, M, G=G)
        ref[key] = (M, G, K1, U1, ws1, Kg, sp1)
    M, G, K1, U1, ws1, Kg, sp1 = ref[key]
    if k == 'S1':
        sp = sp1; label = 'S1'
    else:
        P = None; kw = {}
        path = f'/home/claude/s1spec/step4/feas/snaps/base_{kk}.npz'
        if a.frame == 'fit' and kk >= 100:
            import chart
            d = chart.analyse(path); kw = dict(xc=d['xc'], lam=d['lam']); fr = f"fit lam={d['lam']:.4e}"
        else: fr = 'tip'
        g = numeric_geometry(path, Smax=Sigma + 5, dsig=0.01 if kk >= 100 else 0.005, **kw)
        sp = as_splines(g); label = f'k={kk} ({fr}, kt={g["kt"]:.3f}, lam={g["lam"]:.4e})'
    K = mesh.form(coeffs(mesh, sp, use_tau=(not a.notau) and k != 'S1'))
    Kx = Kg(K); S = 0.5 * (Kx + Kx.T)
    d = 1 / np.sqrt(np.diag(G)); w = sla.eigvalsh(S * d[:, None] * d[None, :], G * d[:, None] * d[None, :])[::-1]
    out = dict(label=label, top=[round(float(x), 4) for x in w[:8]], ndof=len(M), h=h, Sigma=Sigma)
    for c in CS:
        out[f'n_{c}'] = int(np.sum(w > -c)); Mc, tv = minM(U1, S, c) if a.eps2 == 0 else minM(U1, S, c); out[f'M_{c}'] = int(Mc)
    res.append(out)
    print(f"{label}: n(>-.05/.1/.25/.5/1)={out['n_0.05']}/{out['n_0.1']}/{out['n_0.25']}/{out['n_0.5']}/{out['n_1.0']} | M(c=.05/.1/.25/.5/1)={out['M_0.05']}/{out['M_0.1']}/{out['M_0.25']}/{out['M_0.5']}/{out['M_1.0']} | top {out['top']} [{time.time()-t0:.0f}s]", flush=True)
json.dump(res, open(f'/home/claude/s1spec/step4/feas/logs/f3_{a.frame}{"_notau" if a.notau else ""}{"_eps" if a.eps2 else ""}{a.tag}.json', 'w'), indent=1)
