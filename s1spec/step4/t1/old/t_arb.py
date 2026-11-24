import sys, time, numpy as np
from base2 import Base; from arbops import *
B = Base('B300'); i = int(sys.argv[1]) if len(sys.argv) > 1 else 10
R, D = B.rep(i); RD = B.RD; q = B.slabs[i]; nodes = q['nodes']; dt = q['dt']
t0 = time.time(); SP = SlabPolish(R, D, RD, nodes, dt); print('setup', time.time() - t0)
# consistency: float implicit residual at stage 1 vs arb
Yd = (RD.deriv_matrix(RD.c) @ nodes) / dt
p = D.parts(nodes[1]); Mm = D.mass(nodes[1][R.n]); Ef = Mm @ Yd[0][:R.n] - D.rhs(p, Yd[0][R.n])
t0 = time.time(); res = SP.residuals(); print('arb residual time', time.time() - t0)
Ea = np.array([float(v.mid()) for v in res[0][0]])
print('stage1: |E_float| max', np.abs(Ef).max(), ' |E_arb| max', np.abs(Ea).max(), ' |diff|', np.abs(Ef - Ea).max(), ' scalars', [float(v.mid()) for v in res[0][1]])
# polish
Yj = nodes[len(nodes) // 2]; J = D.jac(Yj); RD.factor(J, dt)
for it in range(12):
    t0 = time.time(); em, dz = SP.step(J); print(f'polish it {it}: max|E| before = {em:.2e}, max|dZ| = {dz:.2e}  ({time.time()-t0:.1f}s)', flush=True)
    if em < 1e-34: break
