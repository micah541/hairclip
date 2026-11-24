import numpy as np
from base import Base
B = Base('A'); R, D, RD = B.R, B.D, B.RD
i = int(__import__('sys').argv[1]); Y0 = B.nodes[i, 0]; dt = B.dt[i]
sc = np.ones(len(Y0)); sc[R.n:] = np.maximum(np.abs(Y0[R.n:]), 1.0)
F0 = D.F(Y0); J = D.jac(Y0, F0); RD.factor(J, dt)
Ys, hist = RD.solve(D.F, Y0, dt, scale=sc, tol=1e-15, maxit=12)
print(' '.join(f'{h:.1e}' for h in hist))
# component of the last increment
Ys2, h2 = RD.solve(D.F, Y0, dt, Z0=Ys - Y0, scale=sc, tol=1e-15, maxit=1)
dZ = (Ys2 - Ys)[-1] / sc
j = np.argsort(-np.abs(dZ))[:10]
typ = {}
for jj in range(R.N + 1):
    for m in range(5):
        if R.nidx[jj, m] >= 0: typ[R.nidx[jj, m]] = ('node', jj, m)
for e in range(R.N):
    for q, g in enumerate(R.bidx[e]): typ[g] = ('bub', e, q)
for k, nm in enumerate(['ell', 'lnlam', 'c', 't']): typ[R.n + k] = (nm,)
x, _ = R.mesh.knots(Y0[R.n])
print('largest last-increment components:', [(typ[k], f'{dZ[k]:.1e}', f'x={x[typ[k][1]]:.0f}' if typ[k][0] in ('node', 'bub') else '') for k in j])
print('cond of Newton matrix block 0:', np.linalg.cond(RD.lam[RD.rep[0]] / dt * np.eye(len(Y0)) - J))
