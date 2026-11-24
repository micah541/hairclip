import numpy as np
from rep import *; from ops import *; from base import Base
R = Rep(9, Mesh()); D = Dyn(R); BA = Base('A')
for i in (40, 70, 100, 125):
    Y2 = BA.nodes[i, 0]; F2 = D.F(Y2); J2 = D.jac(Y2, F2)
    w, V = np.linalg.eig(J2[:R.n + 1, :R.n + 1]); o = np.argsort(-w.real)
    x, _ = R.mesh.knots(Y2[R.n]); h = np.diff(x)
    def loc(v):
        a = np.zeros(R.N + 1)
        for j in range(R.N + 1):
            ii = R.nidx[j][R.nidx[j] >= 0]; a[j] = np.sum(np.abs(v[ii])**2)
        return int(np.argmax(a))
    print(f"slab {i} ell={Y2[R.n]:.0f}: ", [(np.round(w[j], 2), loc(V[:, j]), round(x[loc(V[:, j])]), round(h[min(loc(V[:, j]), R.N-1)], 1)) for j in o[:5]])
