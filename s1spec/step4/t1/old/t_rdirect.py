"""independent check: r = V~ + k_xx - (1/4 - a) p - b n from the reconstructed base (exact Gamma_tau) vs the kernel formula"""
import numpy as np
from base2 import Base; from arbops import SlabPolish; import resid
from rep import K0, DS1
B = Base('B300'); i = 10; R, D = B.rep(i); q = B.slabs[i]
SP = SlabPolish(R, D, B.RD, q['nodes'], q['dt'])          # not polished: float base, r ~ 1e-10
sm = resid.samples(SP, 0.37); o, f = resid.norms(sm)
X, H = sm['X'], sm['H']
th = -np.pi / 2 + resid.cumint(sm['k0'], H)
G = -DS1 + resid.cumint(np.exp(1j * th), H)
Gt = resid.cumint(1j * sm['tt'] * np.exp(1j * th), H)     # Gamma_tau at fixed x (tip fixed, theta_tau(0) = 0)
nu = 1j * np.exp(1j * th)
V = (Gt * np.conj(nu)).real; p = (G * np.conj(nu)).real; nn = nu.real
rd = V + sm['k2'] - (0.25 - sm['a']) * p - sm['b'] * nn
rk = f['r']
for xm in (5, 20, 100, sm['ell']):
    m = X <= xm
    print(f'x<={xm:8.1f}: max|r_direct| {np.abs(rd[m]).max():.3e}  max|r_kernel| {np.abs(rk[m]).max():.3e}  max|diff| {np.abs(rd[m]-rk[m]).max():.3e}')
# ODE check: r' = rho_theta - k A, A = int k r ; derivative by per-element Legendre fit
from numpy.polynomial import legendre as LG
def dxf(F):
    out = np.empty_like(F)
    for e in range(F.shape[0]):
        c = LG.legfit(resid.g, F[e], resid.NQ - 1); out[e] = LG.legval(resid.g, LG.legder(c)) * 2 / H[e]
    return out
for nm, rr in (('direct', rd), ('kernel', rk)):
    Ar = resid.cumint(sm['k0'] * rr, H); lhs = dxf(rr); rhs = sm['rt'] - sm['k0'] * Ar
    m = X <= 20
    print(nm, ': max |r\' - (rho_t - k A)| on x<=20:', np.abs(lhs - rhs)[m].max(), '  max|r\'|', np.abs(lhs[m]).max())
# check theta_tau consistency: d/dx tt should be kappa_tau|_x; and rt = tt - Phi:  compare tt - rt (= Phi) with float Phi
T = np.exp(1j * th); qq = (G * np.conj(T)).real; t1 = T.real
dd = (rd - rk); m = X <= 30
Bm = np.c_[p[m], nn[m], qq[m], t1[m], np.ones(m.sum())]
c, *_ = np.linalg.lstsq(Bm, dd[m], rcond=None)
print('fit of r_direct - r_kernel on [0,30] by (p, n, q, t1, 1):', c, ' residual', np.abs(dd[m] - Bm @ c).max(), ' max|diff|', np.abs(dd[m]).max())
Xr = X.ravel(); dr = dd.ravel(); ii = np.linspace(0, int(m.sum()) - 1, 8).astype(int)
print('diff at sample x:', [(round(float(Xr[j]), 2), float(dr[j])) for j in ii])
