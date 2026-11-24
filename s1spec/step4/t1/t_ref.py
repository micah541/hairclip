import numpy as np
from rep import *; from ops import *; from seed import HLSHalf; from run import seed_closed
from run2 import refine
S = HLSHalf('../feas2/snaps/p2_k20.npz')
R = Rep(9, Mesh(ell0=S.ell)); D = Dyn(R); y = seed_closed(R, D, S)
x, _ = R.mesh.knots(S.ell); print('N', R.N, 'first outer knot', x[R.mesh.JT + 1], 'h last tip', x[R.mesh.JT] - x[R.mesh.JT - 1], 'h first outer', x[R.mesh.JT + 1] - x[R.mesh.JT])
R1, y1 = refine(R, y, S.ell, 0.5 * (x[R.mesh.JT] + x[R.mesh.JT + 1]))
xs = np.linspace(0, S.ell, 20001)
for m in range(6): print(m, 'max |kappa^(m) diff| after insertion:', np.max(np.abs(R.eval_at(y, S.ell, xs, m) - R1.eval_at(y1, S.ell, xs, m))))
