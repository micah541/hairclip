import numpy as np
from rep import *; from ops import *; from seed import HLSHalf; from run import seed_closed, geom
S = HLSHalf('../feas2/snaps/p2_k20.npz'); R = Rep(9, Mesh()); D = Dyn(R)
y = seed_closed(R, D, S); Y = np.r_[y, S.ell, np.log(S.lam), 0.0, S.t]
p = D.parts(Y)
GMq = p['G'] @ R.B.Mom.T
print('flux vs quadrature moments: max diff', np.max(np.abs(p['GM'] - GMq)), ' max |GM|', np.max(np.abs(GMq)))
Fv = D.F(Y)
for e_ in (1e-6,):
    gp = geom(R, D, Y + e_ * Fv); gm = geom(R, D, Y - e_ * Fv)
    print('d/dtau turn (FD):', (gp['turn'] - gm['turn']) / (2 * e_), '  d/dtau Im G(ell):', (gp['G_end'].imag - gm['G_end'].imag) / (2 * e_))
for e_ in (1e-3, 1e-4):
    gp = geom(R, D, Y + e_ * Fv); gm = geom(R, D, Y - e_ * Fv)
    print(e_, 'd/dtau turn (FD):', (gp['turn'] - gm['turn']) / (2 * e_), '  d/dtau Im G(ell):', (gp['G_end'].imag - gm['G_end'].imag) / (2 * e_))
# exact: d/dtau int kappa = sum_e [ hdot int P du + h int Pdot du ]
p = D.parts(Y); B = R.B; h = p['h']; dh = np.diff(p['dx']) * Fv[R.n]
yd = Fv[:R.n]
dl = R.local(np.r_[yd], h) - np.where(R.L >= 0, 0, 0)  # local coefficients of ydot (fixed entries contribute fixval -> must remove)
dl = np.where(R.L >= 0, np.r_[yd, 0.0][R.L], 0.0) * R.lscale(h)
Hn, dHn = R.nodal_scale(h, np.diff(p['dx']))
mm = np.arange(5)[None, :]
dA = np.zeros_like(p['d']); dA[:, :5] = p['d'][:, :5] * mm * (np.diff(p['dx']) / h - dHn[:-1] / Hn[:-1])[:, None]; dA[:, 5:10] = p['d'][:, 5:10] * mm * (np.diff(p['dx']) / h - dHn[1:] / Hn[1:])[:, None]
ddu = dl + Fv[R.n] * dA          # d/dtau of local coefficients at fixed u
intP = p['d'] @ (B.V[0].T @ B.wq); intPd = ddu @ (B.V[0].T @ B.wq)
print('exact d/dtau int kappa =', np.sum(dh * intP + h * intPd))
