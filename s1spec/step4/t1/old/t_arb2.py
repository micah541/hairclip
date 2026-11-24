import sys, numpy as np
from base2 import Base; from arbops import *
B = Base('B300'); i = 10
R, D = B.rep(i); RD = B.RD; q = B.slabs[i]; nodes = q['nodes']; dt = q['dt']
SP = SlabPolish(R, D, RD, nodes, dt)
Yd = (RD.deriv_matrix(RD.c) @ nodes) / dt
p = D.parts(nodes[1]); Mm = D.mass(nodes[1][R.n]); Ef = Mm @ Yd[0][:R.n] - D.rhs(p, Yd[0][R.n])
res = SP.residuals(); Ea = np.array([float(v.mid()) for v in res[0][0]])
nm = D.nm; dif = np.abs(Ef - Ea).reshape(R.N, nm + 1)
print('max diff per row type (moments m=0..nm-1, galerkin):', dif.max(0))
e = np.argmax(dif.max(1)); print('worst element', e, 'of', R.N, ' x =', p['x'][e], ' h =', p['h'][e], ' row diffs', dif[e])
print('elements with diff > 1e-10:', np.where(dif.max(1) > 1e-10)[0])
