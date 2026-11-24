import numpy as np
from base2 import Base
B = Base('P11')
for i in (60, 120, 160, 187):
    R, D = B.rep(i); Y = B.slabs[i]['nodes'][0]; J = D.jac(Y)
    w = np.linalg.eigvals(J[:R.n + 1, :R.n + 1]); o = np.argsort(-w.real)
    pos = w[(w.real > 0) & (np.abs(w - (0.25 - D.parts(Y)['a'])) > 1e-3)]
    print(f"slab {i} k={B.k(B.tau0[i]):.0f} n={R.n}: top {np.round(w[o][:4], 4)}; other Re>0: {len(pos)} (max {pos.real.max() if len(pos) else 0:.2e}); min Re {w.real.min():.2e}", flush=True)
