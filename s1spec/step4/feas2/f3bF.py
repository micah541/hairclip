"""F3b variant F: T-clock frame  lam = (T-t)^(1/4)  (mu = 1/4 exactly), FIXED centre x* = tip position at T, transport tau included.
usage: python3 f3bF.py statefiles..."""
import sys, json, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from f3geom import *; from f3fem import *
from scipy.interpolate import CubicSpline
CS = (0.05, 0.1, 0.25, 0.5, 1.0); TT = 4.598863926611
S_T = Snap('/home/claude/s1spec/step4/feas/snaps/base_1600.npz'); import hls
A_ = hls.HLS(np.zeros(S_T.N), S_T.L, S_T.X0, omega=1); A_.th_hat = S_T.h.copy(); _, _, ZT, _ = A_.curve(1)
XSTAR = ZT[S_T.N // 2]; print('x* = tip(T) =', XSTAR, flush=True)
def minM(U, S, c):
    SU = U.T @ S @ U; n = SU.shape[0]; top = lambda M: np.linalg.eigvalsh(SU[M:, M:])[-1]
    lo, hi = 0, min(n - 1, 40)
    if top(0) <= -c: return 0
    if top(hi) > -c: return -1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if top(mid) <= -c: hi = mid
        else: lo = mid
    return hi
def cf(mesh, sp, tau=None):
    x = mesh.xq; k, k1, k2, q = sp['kap'](x), sp['kap1'](x), sp['kap2'](x), sp['q'](x)
    return dict(b2=-k**2, b1=-5*k*k1 - q/4 + (tau(x) if tau else 0.0), b0=0.25 - 3*k1**2 - 4*k*k2)
res = {}
for f in sys.argv[1:]:
    z = np.load(f); t = float(z['t']); L = float(z['L']); lam = (TT - t) ** 0.25
    Sig = min(35.0, L / (2 * lam) - 5.5); mesh = Mesh(Sig, 0.05); M = mesh.mass()
    spS1 = s1_splines(Smax=Sig + 5, dx=0.01); wS1, U1 = eig_sym(mesh.form(cf(mesh, spS1)), M)
    g = numeric_geometry(f, xc=XSTAR, lam=lam, Smax=Sig + 5, dsig=0.01); sp = as_splines(g)
    tipY = (g['tip'] - XSTAR) / lam
    K = mesh.form(cf(mesh, sp, tau=sp['tau'])); S = 0.5 * (K + K.T); w, V = eig_sym(K, M)
    o = dict(t=t, lam=lam, Sigma=Sig, tipY=[tipY.real, tipY.imag], kt_lam=float(g['kt'] * lam), top=[round(float(x), 4) for x in w[:10]],
             n={str(c): int(np.sum(w > -c)) for c in CS}, M={str(c): minM(U1, S, c) for c in CS})
    res[f't={t:g}'] = o
    print(f"t={t:.4f} lamF={lam:.4f} Sigma={Sig:.1f} tipY={tipY.real:.2f} kt*lam={o['kt_lam']:.3f}: top {o['top'][:8]} n={[o['n'][str(c)] for c in CS]} M={[o['M'][str(c)] for c in CS]}", flush=True)
json.dump(res, open('logs/f3bF.json', 'w'), indent=1)
