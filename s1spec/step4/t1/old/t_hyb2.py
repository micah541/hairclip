import numpy as np
from rep import *; from ops import *; from seed import HLSHalf; from run import seed_closed, geom; from radau import Radau
S = HLSHalf('../feas2/snaps/p2_k20.npz'); R = Rep(9, Mesh()); D = Dyn(R)
y = seed_closed(R, D, S); Y = np.r_[y, S.ell, np.log(S.lam), 0.0, S.t]
for ell in (S.ell, 4e4, 3.7e5):
    M = D.mass(ell).toarray(); s = np.linalg.svd(M, compute_uv=False); print(f'ell={ell:.0f} cond(M)={s[0]/s[-1]:.2e}')
Fv, p = D.F(Y, True); B = R.B; h = p['h']
# exact d/dtau int kappa
dxl = p['dx']; dhh = np.diff(dxl) * Fv[R.n]
dl = np.where(R.L >= 0, np.r_[Fv[:R.n], 0.0][R.L], 0.0) * R.lscale(h)
Hn, dHn = R.nodal_scale(h, np.diff(dxl)); mm = np.arange(5)[None, :]
dA = np.zeros_like(p['d']); dA[:, :5] = p['d'][:, :5] * mm * (np.diff(dxl) / h - dHn[:-1] / Hn[:-1])[:, None]; dA[:, 5:10] = p['d'][:, 5:10] * mm * (np.diff(dxl) / h - dHn[1:] / Hn[1:])[:, None]
ddu = dl + Fv[R.n] * dA
print('exact d/dtau int kappa =', np.sum(dhh * (p['d'] @ (B.V[0].T @ B.wq)) + h * (ddu @ (B.V[0].T @ B.wq))))
J = D.jac(Y, Fv); w = np.linalg.eigvals(J[:R.n + 1, :R.n + 1]); o = np.argsort(-w.real)
print('top eigenvalues', np.round(w[o][:5], 4), ' min Re', w.real.min(), ' max |Im|', np.abs(w.imag).max())
RD = Radau(8); RD.factor(J, 0.15); sc = np.ones(len(Y)); sc[R.n:] = np.maximum(np.abs(Y[R.n:]), 1)
Ys, hist = RD.solve(D.F, Y, 0.15, scale=sc, tol=1e-15, maxit=15); print('Newton:', ' '.join(f'{x:.1e}' for x in hist))
# late state from run A (slab 125) as a test of the stretched mesh
from base import Base
BA = Base('A'); Y2 = BA.nodes[125, 0]
F2 = D.F(Y2); J2 = D.jac(Y2, F2); RD.factor(J2, 0.15); sc2 = np.ones(len(Y2)); sc2[R.n:] = np.maximum(np.abs(Y2[R.n:]), 1)
Ys2, hist2 = RD.solve(D.F, Y2, 0.15, scale=sc2, tol=1e-15, maxit=15); print('Newton late:', ' '.join(f'{x:.1e}' for x in hist2))
w2 = np.linalg.eigvals(J2[:R.n + 1, :R.n + 1]); o2 = np.argsort(-w2.real); print('late top eigenvalues', np.round(w2[o2][:5], 4), ' min Re', w2.real.min())
