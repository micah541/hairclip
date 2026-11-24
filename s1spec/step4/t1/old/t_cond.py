import numpy as np
from rep import *; from ops import *
for P in (9, 11):
    R = Rep(P, Mesh()); D = Dyn(R); M = D.mass(314.5).toarray()
    cn = np.linalg.norm(M, axis=0); Ms = M / cn[None, :]
    rn = np.linalg.norm(Ms, axis=1); Ms2 = Ms / rn[:, None]
    s = np.linalg.svd(Ms2, compute_uv=False)
    print(P, 'cond raw %.2e  col-equil %.2e  row+col %.2e' % (np.linalg.cond(M), np.linalg.cond(Ms), s[0] / s[-1]), ' smallest sv:', s[-6:])
    U, S, Vt = np.linalg.svd(Ms2); v = Vt[-1] / cn
    # which DOF types dominate the weakest mode
    typ = np.zeros(R.n, int) - 1; node = np.zeros(R.n, int)
    for j in range(R.N + 1):
        for m in range(5):
            if R.nidx[j, m] >= 0: typ[R.nidx[j, m]] = m; node[R.nidx[j, m]] = j
    for e in range(R.N):
        typ[R.bidx[e]] = 9; node[R.bidx[e]] = e
    i = np.argsort(-np.abs(Vt[-1]))[:8]
    print('   weakest mode components (type, node, val):', [(typ[k], node[k], round(Vt[-1][k], 3)) for k in i])
