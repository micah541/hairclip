import numpy as np
from handoff2 import *
d = analyse('snap_extA_2240.npz'); u = d['u']; n = len(u); w = u - d['alpha']*PSI3[:n] - d['beta2']*H2[:n]
sw = np.sqrt(wh[:n])
for name, cols in [('p1,sin th1', [M.p1, np.sin(M.th1)]), ('+ q-rot (q1)', [M.p1, np.sin(M.th1), M.q1]), ('+psi3,h2 re-fit', [M.p1, np.sin(M.th1), PSI3, H2])]:
    B = np.column_stack([c[:n] for c in cols]) * sw[:, None]; c, *_ = np.linalg.lstsq(B, w*sw, rcond=None)
    print(name, 'coef', c, 'resid/||w||', np.linalg.norm(w*sw - B@c)/np.linalg.norm(w*sw))
