import numpy as np
from rep import *; from ops import *
R = Rep(9, Mesh()); D = Dyn(R)
for ell in (314.5, 3e3, 1e4, 4e4, 3.7e5):
    M = D.mass(ell).toarray(); s = np.linalg.svd(M, compute_uv=False)
    x, _ = R.mesh.knots(ell); h = np.diff(x); r = h[1:] / h[:-1]
    print(f"ell={ell:9.1f}: cond(M)={s[0]/s[-1]:.2e} smallest sv {s[-3:]}  h ratio range [{r.min():.3f},{r.max():.3f}]  h_max={h.max():.1f}")
