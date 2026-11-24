import numpy as np
from rep import *; from ops import *; from seed import HLSHalf; from run import seed_closed, geom; from radau import Radau
S = HLSHalf('../feas2/snaps/p2_k20.npz'); R = Rep(9, Mesh(ell0=S.ell)); D = Dyn(R)
y = seed_closed(R, D, S); Y = np.r_[y, S.ell, np.log(S.lam), 0.0, S.t]
Fv, p = D.F(Y, True); B = R.B; h = p['h']; dxl = p['dx']
dhh = np.diff(dxl) * Fv[R.n]
dl = np.where(R.L >= 0, np.r_[Fv[:R.n], 0.0][R.L], 0.0) * R.lscale(h)
Hn, dHn = R.nodal_scale(h, np.diff(dxl)); mm = np.arange(5)[None, :]
dA = np.zeros_like(p['d']); dA[:, :5] = p['d'][:, :5] * mm * (np.diff(dxl) / h - dHn[:-1] / Hn[:-1])[:, None]; dA[:, 5:10] = p['d'][:, 5:10] * mm * (np.diff(dxl) / h - dHn[1:] / Hn[1:])[:, None]
ddu = dl + Fv[R.n] * dA
print('exact d/dtau int kappa =', np.sum(dhh * (p['d'] @ (B.V[0].T @ B.wq)) + h * (ddu @ (B.V[0].T @ B.wq))))
print('cond M', np.linalg.cond(p['Mm'].toarray()))
J = D.jac(Y, Fv); w = np.linalg.eigvals(J[:R.n + 1, :R.n + 1]); o = np.argsort(-w.real); print('top eig', np.round(w[o][:4], 3))
RD = Radau(8); RD.factor(J, 0.15); sc = np.ones(len(Y)); sc[R.n:] = np.maximum(np.abs(Y[R.n:]), 1)
Ys, hist = RD.solve(D.F, Y, 0.15, scale=sc, tol=1e-15, maxit=40); print('Newton:', ' '.join(f'{x:.1e}' for x in hist))
Ys1, h1 = RD.solve(D.F, Y, 0.15, Z0=Ys - Y, scale=sc, tol=1e-15, maxit=1)
dZ = (Ys1 - Ys) / sc[None, :]
i, j = np.unravel_index(np.argmax(np.abs(dZ)), dZ.shape)
typ = {}
for jj in range(R.N + 1):
    for m in range(5):
        if R.nidx[jj, m] >= 0: typ[R.nidx[jj, m]] = ('node', jj, m)
for e in range(R.N):
    for q, g in enumerate(R.bidx[e]): typ[g] = ('bub', e, q)
for kk, nm in enumerate(['ell', 'lnlam', 'c', 't']): typ[R.n + kk] = (nm,)
top = np.argsort(-np.abs(dZ).max(0))[:8]
x, _ = R.mesh.knots(S.ell)
print('slow components:', [(typ[k], f'{np.abs(dZ[:, k]).max():.1e}') for k in top])
# linearity / Jacobian check along the slow direction
v = np.zeros(len(Y)); v[top[0]] = 1.0
for d_ in (1e-6, 1e-8):
    Fd = (D.F(Y + d_ * v) - D.F(Y - d_ * v)) / (2 * d_); print(d_, 'J col err', np.max(np.abs(Fd - J[:, top[0]])), ' |col|', np.max(np.abs(Fd)))
