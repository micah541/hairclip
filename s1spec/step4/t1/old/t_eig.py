import sys, numpy as np
from base import Base
B = Base('A'); R, D, RD = B.R, B.D, B.RD
for i in map(int, sys.argv[1:]):
    Y0 = B.nodes[i, 0]; F0 = D.F(Y0); J = D.jac(Y0, F0)
    w, V = np.linalg.eig(J[:R.n + 1, :R.n + 1])
    o = np.argsort(-w.real)
    x, _ = R.mesh.knots(Y0[R.n])
    def loc(v):
        a = np.zeros(R.N + 1)
        for j in range(R.N + 1):
            ii = R.nidx[j][R.nidx[j] >= 0]; a[j] = np.sum(np.abs(v[ii])**2)
        return x[np.argmax(a)]
    print(f"slab {i} k={B.k(B.tau0[i]):.1f}: top real parts:", [(np.round(w[j], 3), round(loc(V[:, j]))) for j in o[:6]])
    big = np.abs(w) > 1e3
    print('   |w|>1e3: count', big.sum(), ' max Re', w.real[big].max(), ' max |Im|', np.abs(w.imag[big]).max(), ' max |w|', np.abs(w).max())
    print('   Radau lam/dt:', np.round(RD.lam / B.dt[i], 1))
